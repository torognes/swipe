/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2014 Torbjorn Rognes, University of Oslo,
    Oslo University Hospital and Sencel Bioinformatics AS

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

// The 7-bit search kernel with AVX2 (32 database sequences at once,
// __m256i), compiled with -mavx2 and selected at run time
// (cpu_features.avx2): the same algorithm as search7_ssse3(), twice
// as many channels.

#include "swipe.h"
#include "intrinsics_to_functions.h"  // v_load, v_store, v_shuffle_8, ...
#include "align_cells.h"  // Ops_7_avx2, align_cells(), No_mask, Mask
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstring>  // std::memset
#include <iterator>  // std::next

#ifndef __AVX2__
#error "search7_avx2.cc must be compiled with -mavx2"
#endif

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

constexpr std::size_t CHANNELS = channels_7_avx2;
constexpr std::size_t CDEPTH = 4;

static_assert(CHANNELS == sizeof(__m256i), "one byte lane per channel");

// the byte 0x80 (-128 has the same bits)
constexpr char byte_0x80 = static_cast<char>(-128);

// the shuffle indices of the database symbols of one position:
// symbols 0 to 15 select the low table (bit 7 set for 16 to 31: zero),
// symbols 16 to 31 the high table
struct ShuffleIndices
{
  __m256i low;
  __m256i high;
};

inline auto shuffle_indices(__m256i const symbols) -> ShuffleIndices
{
  auto const bit_4 = v256_dup_i8(0x10);
  auto const bit_7 = v256_dup_i8(byte_0x80);
  // bit 4 of the symbol moved to bit 7: set for symbols 16 to 31
  auto const high = v_shift_left_i16<3>(v_and(symbols, bit_4));
  return {v_or(symbols, high), v_or(symbols, v_xor(high, bit_7))};
}

inline auto profile_vector(__m256i const low_table,
                           __m256i const high_table,
                           ShuffleIndices const & indices) -> __m256i
{
  return v_or(v_shuffle_8(low_table, indices.low),
              v_shuffle_8(high_table, indices.high));
}

// The score profile of a block of CDEPTH x CHANNELS database residues:
// for each of the 32 query symbols, CDEPTH vectors (one per database
// position), each holding the score of that query symbol against the
// residue of every channel. A row of the transposed score matrix (32
// bytes: the scores of one query symbol against the 32 database
// symbols) is two 16-entry tables; vpshufb looks up within each
// 128-bit half, so each table is broadcast to both halves. The index
// vectors select the low table for symbols 0 to 15 (bit 7 set
// otherwise: zero) and the high one for symbols 16 to 31.
inline auto dprofile_shuffle7(BYTE * dprofile,
                              BYTE const * score_matrix,
                              BYTE const * dseq_byte) -> void
{
  static_assert(CDEPTH == 4, "four database positions per block");
  auto const * const dseq = reinterpret_cast<__m256i const *>(dseq_byte);
  auto const position_0 = shuffle_indices(v_load(dseq));
  auto const position_1 = shuffle_indices(v_load(std::next(dseq, 1)));
  auto const position_2 = shuffle_indices(v_load(std::next(dseq, 2)));
  auto const position_3 = shuffle_indices(v_load(std::next(dseq, 3)));

  auto const * const matrix = reinterpret_cast<__m128i const *>(score_matrix);
  auto * const profile = reinterpret_cast<__m256i *>(dprofile);
  for (std::ptrdiff_t row = 0; row < static_cast<std::ptrdiff_t>(score_matrix_width); ++row)
  {
    auto const low_table = v256_broadcast_128(std::next(matrix, 2 * row));
    auto const high_table = v256_broadcast_128(std::next(matrix, (2 * row) + 1));
    auto * const line = std::next(profile, static_cast<std::ptrdiff_t>(CDEPTH) * row);
    v_store(line, profile_vector(low_table, high_table, position_0));
    v_store(std::next(line, 1), profile_vector(low_table, high_table, position_1));
    v_store(std::next(line, 2), profile_vector(low_table, high_table, position_2));
    v_store(std::next(line, 3), profile_vector(low_table, high_table, position_3));
  }
}

}  // anonymous namespace


auto search7_avx2(BYTE * * q_start,
                  BYTE gap_open_penalty,
                  BYTE gap_extend_penalty,
                  BYTE const * score_matrix,
                  BYTE * dprofile,
                  BYTE * hearray,
                  db_thread_s & dbt,
                  long sequences,
                  long const * seqnos,
                  long * scores,
                  long qlen) -> void
{
  auto * const hep = reinterpret_cast<__m256i*>(hearray);
  auto ** const qp = reinterpret_cast<__m256i**>(q_start);
  std::array<BYTE const *, CHANNELS> d_begin;
  std::array<BYTE const *, CHANNELS> d_end;

  // the database residues of the channels, aligned for the loads
  alignas(__m256i) std::array<BYTE, CDEPTH * CHANNELS> dseq {{}};
  // the lanes of the channels that start a new sequence (0x80), and
  // the scores, read lane by lane
  alignas(__m256i) std::array<BYTE, CHANNELS> restart {{}};
  alignas(__m256i) std::array<BYTE, CHANNELS> lane_scores {{}};
  BYTE const zero = 0;

  std::array<long, CHANNELS> seq_id;
  long next_id = 0;
  long done = 0;

  std::memset(hearray, 0x80, static_cast<std::size_t>(qlen) * hearray_row_bytes_avx2);

  auto const Z = v256_dup_i8(byte_0x80);
  auto const Q = v256_dup_i8(static_cast<char>(gap_open_penalty));
  auto const R = v256_dup_i8(static_cast<char>(gap_extend_penalty));

  auto S = Z;

  for (std::size_t c = 0; c < CHANNELS; c++)
  {
    d_begin[c] = &zero;
    d_end[c] = d_begin[c];
    seq_id[c] = -1;
  }

  // the next CDEPTH residues of channel c (zeros after its end);
  // false when the channel reached the end of its sequence
  auto const fill_channel = [&](std::size_t const c) -> bool
  {
    for (std::size_t j = 0; j < CDEPTH; j++)
    {
      if (d_begin[c] < d_end[c])
      {
        dseq[(CHANNELS * j) + c] = *d_begin[c];
        d_begin[c] = std::next(d_begin[c]);
      }
      else
      {
        dseq[(CHANNELS * j) + c] = 0;
      }
    }
    return d_begin[c] != d_end[c];
  };

  bool easy = false;

  while (true)
  {
    if (easy)
    {
      // fill all channels
      for (std::size_t c = 0; c < CHANNELS; c++)
      {
        if (not fill_channel(c))
        {
          easy = false;
        }
      }

      dprofile_shuffle7(dprofile, score_matrix, dseq.data());
      align_cells<Ops_7_avx2>(S, hep, qp, Q, R, qlen, Z, No_mask{});
      continue;
    }

    // One or more sequences ended in the previous block: switch over
    // to new sequences
    easy = true;
    restart.fill(0);
    v_store(reinterpret_cast<__m256i *>(lane_scores.data()), S);

    for (std::size_t c = 0; c < CHANNELS; c++)
    {
      if (d_begin[c] < d_end[c])
      {
        // this channel has more sequence
        if (not fill_channel(c))
        {
          easy = false;
        }
        continue;
      }

      // sequence in channel c ended: change of sequence
      restart[c] = 0x80;

      long const cand_id = seq_id[c];
      if (cand_id >= 0)
      {
        // save score
        scores[cand_id] = lane_scores[c] - 0x80;
        done++;
      }

      if (next_id < sequences)
      {
        // get next sequence
        seq_id[c] = next_id;
        long const seqnosf = seqnos[next_id];
        View<char> const sequence =
          db_getsequence(dbt, entry_seqno(seqnosf), entry_where(seqnosf), c);
        d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
        d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
        next_id++;
        if (not fill_channel(c))
        {
          easy = false;
        }
      }
      else
      {
        // no more sequences, empty channel
        seq_id[c] = -1;
        d_begin[c] = &zero;
        d_end[c] = d_begin[c];
        static_cast<void>(fill_channel(c));
      }
    }

    if (done == sequences)
    {
      break;
    }

    dprofile_shuffle7(dprofile, score_matrix, dseq.data());
    auto const mask = v_load(reinterpret_cast<__m256i const *>(restart.data()));
    align_cells<Ops_7_avx2>(S, hep, qp, Q, R, qlen, Z, Mask<Ops_7_avx2>{mask});
  }
}
