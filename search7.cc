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

#include "swipe.h"
#include "search_kernels.h"  // the prototype of the kernel
#include "intrinsics_to_functions.h"  // v128_t, v_load, v_store, v_merge_*, v_dup_*, ...
#include "align_cells.h"  // Ops_7, align_cells(), No_mask, Mask
#include "view.h"  // View
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstring>  // std::memset
#include <iterator>  // std::next

constexpr std::size_t CHANNELS = channels_7;
static_assert(sizeof(v128_t) == vector_bytes, "an SSE vector");
constexpr std::size_t CDEPTH = 4;

// the byte 0x80 (the lanes of _mm_set_epi8() are char: 0x80 does not
// fit in a signed char, -128 has the same bits)
constexpr char byte_0x80 = static_cast<char>(-128);

// the zero score of a 7-bit lane (the scores are stored with a bias of
// 0x80), which is also the lane value that resets a score to zero
// (Ops_7::mask(), align_cells.h)
constexpr BYTE zero_score = 0x80;

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

#ifdef SWIPE_SSSE3

// profile_row(row) strides, in 16-byte vectors: a row of the 32 x 32 score
// matrix is two vectors, a row of the profile one vector per CDEPTH
constexpr std::ptrdiff_t matrix_row_vectors = 2;
constexpr auto profile_row_vectors = static_cast<std::ptrdiff_t>(CDEPTH);

// the shuffle indices of the database symbols of one position:
// symbols 0 to 15 select the low table (bit 7 set for 16 to 31: zero),
// symbols 16 to 31 the high table
struct ShuffleIndices
{
  v128_t low;
  v128_t high;
};

inline auto shuffle_indices(v128_t const symbols) -> ShuffleIndices
{
  auto const bit_4 = v_dup_i8(0x10);
  auto const bit_7 = v_dup_i8(byte_0x80);
  // bit 4 of the symbol moved to bit 7: set for symbols 16 to 31
  auto const high = v_shift_left_i16<3>(v_and(symbols, bit_4));
  return {v_or(symbols, high), v_or(symbols, v_xor(high, bit_7))};
}

inline auto profile_vector(v128_t const low_table,
                           v128_t const high_table,
                           ShuffleIndices const & indices) -> v128_t
{
  return v_or(v_shuffle_8(low_table, indices.low),
              v_shuffle_8(high_table, indices.high));
}

inline auto dprofile_shuffle7(BYTE * dprofile,
			      BYTE const * score_matrix,
			      BYTE const * dseq_byte) -> void
{
  static_assert(CDEPTH == 4, "four database positions per block");
  auto const * const dseq = reinterpret_cast<v128_t const *>(dseq_byte);
  
  // 16 x 4 = 64 db symbols
  // ca 458 instructions

  // make masks

  /* Note: pshufb only on modern Intel cpus (SSSE3), not AMD */
  /* SSSE3: Supplemental SSE3 */

  auto const position_0 = shuffle_indices(v_load(dseq));
  auto const position_1 = shuffle_indices(v_load(std::next(dseq, 1)));
  auto const position_2 = shuffle_indices(v_load(std::next(dseq, 2)));
  auto const position_3 = shuffle_indices(v_load(std::next(dseq, 3)));

  // one row of the score matrix, shuffled into the four vectors of a
  // row of the profile (the former profline(j) macro)
  auto const * const matrix = reinterpret_cast<v128_t const *>(score_matrix);
  auto * const profile = reinterpret_cast<v128_t *>(dprofile);
  auto const profile_row = [&](std::ptrdiff_t const row) -> void
  {
    auto const * const matrix_row = std::next(matrix, matrix_row_vectors * row);
    auto * const profile_line = std::next(profile, profile_row_vectors * row);
    auto const low_table = v_load(matrix_row);
    auto const high_table = v_load(std::next(matrix_row));
    v_store(profile_line, profile_vector(low_table, high_table, position_0));
    v_store(std::next(profile_line, 1), profile_vector(low_table, high_table, position_1));
    v_store(std::next(profile_line, 2), profile_vector(low_table, high_table, position_2));
    v_store(std::next(profile_line, 3), profile_vector(low_table, high_table, position_3));
  };

  // four rows per iteration: GCC does not unroll the 32 rows of a
  // one-row loop, whose overhead cost +0.54% instructions (callgrind,
  // sprot50k); four rows cost +0.09% against the 32 unrolled macros
  static_assert(score_matrix_width % 4 == 0, "rows in groups of four");
  for (std::ptrdiff_t row = 0; row < static_cast<std::ptrdiff_t>(score_matrix_width); row += 4)
  {
    profile_row(row);
    profile_row(row + 1);
    profile_row(row + 2);
    profile_row(row + 3);
  }
}

#else

// the score matrix row of the database symbol of each channel (byte
// offsets: a row is 32 bytes)
using ChannelRows = std::array<std::ptrdiff_t, CHANNELS>;

// the score matrix rows are transposed by groups of 8 columns (8 bytes
// per channel); a row is two vectors of 16 columns
constexpr std::ptrdiff_t group_columns = 8;
constexpr auto vector_columns = static_cast<std::ptrdiff_t>(sizeof(v128_t));

// The 8 columns (query symbols) from 'column' of the score matrix rows
// of the 16 channels (8 x 16 bytes), transposed into the 8 profile
// vectors of these query symbols, database position j: four rounds of
// unpacks, of 8-, 16-, 32- and 64-bit elements
template <std::ptrdiff_t column>
inline auto dprofile_fill7_columns(v128_t * const profile,
                                   BYTE const * const score_matrix,
                                   ChannelRows const & rows,
                                   std::size_t const j) -> void
{
  auto const row = [&](std::size_t const c) -> v128_t const *
  {
    return reinterpret_cast<v128_t const *>(std::next(score_matrix, rows[c] + column));
  };
  // the 8 bytes of an odd channel: when they start a 16-byte boundary,
  // a full load, which the unpack takes as its memory operand
  // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.
  auto const odd_row = [&](std::size_t const c) -> v128_t
  {
    return (column % vector_columns == 0) ? v_load(row(c)) : v_load_64(row(c));
  };
  auto const pair = [&](std::size_t const c) -> v128_t
  {
    return v_merge_lo_8(v_load_64(row(c)), odd_row(c + 1));
  };
  auto const b0 = pair(0);
  auto const b1 = pair(2);
  auto const b2 = pair(4);
  auto const b3 = pair(6);
  auto const b4 = pair(8);
  auto const b5 = pair(10);
  auto const b6 = pair(12);
  auto const b7 = pair(14);

  auto const w0 = v_merge_lo_16(b0, b1);
  auto const w1 = v_merge_hi_16(b0, b1);
  auto const w2 = v_merge_lo_16(b2, b3);
  auto const w3 = v_merge_hi_16(b2, b3);
  auto const w4 = v_merge_lo_16(b4, b5);
  auto const w5 = v_merge_hi_16(b4, b5);
  auto const w6 = v_merge_lo_16(b6, b7);
  auto const w7 = v_merge_hi_16(b6, b7);

  auto const l0 = v_merge_lo_32(w0, w2);
  auto const l1 = v_merge_hi_32(w0, w2);
  auto const l2 = v_merge_lo_32(w1, w3);
  auto const l3 = v_merge_hi_32(w1, w3);
  auto const l4 = v_merge_lo_32(w4, w6);
  auto const l5 = v_merge_hi_32(w4, w6);
  auto const l6 = v_merge_lo_32(w5, w7);
  auto const l7 = v_merge_hi_32(w5, w7);

  // profile vector of query symbol column + s, database position j
  auto const store = [&](std::ptrdiff_t const s, v128_t const vector) -> void
  {
    auto const symbol = column + s;
    v_store(std::next(profile, (static_cast<std::ptrdiff_t>(CDEPTH) * symbol) + static_cast<std::ptrdiff_t>(j)), vector);
  };
  store(0, v_merge_lo_64(l0, l4));
  store(1, v_merge_hi_64(l0, l4));
  store(2, v_merge_lo_64(l1, l5));
  store(3, v_merge_hi_64(l1, l5));
  store(4, v_merge_lo_64(l2, l6));
  store(5, v_merge_hi_64(l2, l6));
  store(6, v_merge_lo_64(l3, l7));
  store(7, v_merge_hi_64(l3, l7));
}

inline auto dprofile_fill7(BYTE * dprofile,
			   BYTE const * score_matrix,
			   BYTE const * dseq) -> void
{
  static_assert(score_matrix_width == 4 * group_columns, "four groups of 8 columns");
  auto * const profile = reinterpret_cast<v128_t *>(dprofile);
  
  // 4 x 16 db symbols
  // ca (60x2+68x2)x4 = 976 instructions

  for (std::size_t j = 0; j < CDEPTH; j++)
  {
    ChannelRows rows {{}};
    for (std::size_t i = 0; i < CHANNELS; i++)
    {
      rows[i] = static_cast<std::ptrdiff_t>(dseq[(j * CHANNELS) + i] * score_matrix_width);
    }

    dprofile_fill7_columns<0>(profile, score_matrix, rows, j);
    dprofile_fill7_columns<group_columns>(profile, score_matrix, rows, j);
    dprofile_fill7_columns<2 * group_columns>(profile, score_matrix, rows, j);
    dprofile_fill7_columns<3 * group_columns>(profile, score_matrix, rows, j);
  }

}

#endif

}  // anonymous namespace

void
#ifdef SWIPE_SSSE3
search7_ssse3
#else
search7
#endif
       (BYTE * * q_start,
	BYTE gap_open_penalty,
	BYTE gap_extend_penalty,
	BYTE const * score_matrix,
	BYTE * dprofile,
	BYTE * hearray,
	db_thread_s & dbt,
	long sequences,
	long const * seqnos,
	long * scores,
	long qlen)
{
  auto * const hep = reinterpret_cast<v128_t*>(hearray);
  auto const * const * const qp = reinterpret_cast<v128_t const * const *>(q_start);
  std::array<BYTE const *, CHANNELS> d_begin {{}};
  std::array<BYTE const *, CHANNELS> d_end {{}};
  
  // the database residues of the channels, 16-byte aligned for the loads
  alignas(v128_t) std::array<BYTE, CDEPTH * CHANNELS> dseq {{}};
  // the lanes of the channels that start a new sequence (0x80), and
  // the scores, read lane by lane
  alignas(v128_t) std::array<BYTE, CHANNELS> restart {{}};
  alignas(v128_t) std::array<BYTE, CHANNELS> lane_scores {{}};
  BYTE const zero = 0;

  std::array<long, CHANNELS> seq_id {{}};
  long next_id = 0;
  long done = 0;
  
  std::memset(hearray, zero_score, static_cast<std::size_t>(qlen) * hearray_row_bytes);

  auto const Z = v_dup_i8(byte_0x80);
  GapVectors<Ops_7> const gaps {v_dup_i8(static_cast<char>(gap_open_penalty)),
                                v_dup_i8(static_cast<char>(gap_extend_penalty))};

  auto S = Z;

  d_begin.fill(&zero);
  d_end.fill(&zero);
  seq_id.fill(-1);

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

  while(true)
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

#ifdef SWIPE_SSSE3
      dprofile_shuffle7(dprofile, score_matrix, dseq.data());
#else
      dprofile_fill7(dprofile, score_matrix, dseq.data());
#endif

      align_cells<Ops_7>(S, hep, qp, gaps, qlen, Z, No_mask{});
      continue;
    }

    // One or more sequences ended in the previous block 
    // We have to switch over to a new sequence

    easy = true;
    restart.fill(0);
    v_store(reinterpret_cast<v128_t *>(lane_scores.data()), S);

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

      // sequence in channel c ended
      // change of sequence

      restart[c] = zero_score;

      long const cand_id = seq_id[c];
		  
      if (cand_id >= 0)
      {
        // save score
        scores[cand_id] = lane_scores[c] - zero_score;
        done++;
      }

      if (next_id < sequences)
      {
        // get next sequence
        seq_id[c] = next_id;
        long const seqnosf = seqnos[next_id];

        View<char> const sequence =
          db_getsequence(dbt, entry_seqno(seqnosf), entry_where(seqnosf), c);
		      
        // printf("Seqno: %ld Address: %p\n", seqno, address);
        d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
        d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
        next_id++;
		      
        // fill channel
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

#ifdef SWIPE_SSSE3
    dprofile_shuffle7(dprofile, score_matrix, dseq.data());
#else
    dprofile_fill7(dprofile, score_matrix, dseq.data());
#endif

    auto const M = v_load(reinterpret_cast<v128_t const *>(restart.data()));
    align_cells<Ops_7>(S, hep, qp, gaps, qlen, Z, Mask<Ops_7>{M});
  }
}
