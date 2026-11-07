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

// What search7() (search7.cc) and search7_ssse3()
// (arch/x86_64/SSSE3/search7_ssse3.cc) share: the 16 byte lanes of an
// SSE vector, and the search loop over blocks of 4 database residues
// per channel, generic over the score profile builder (SSE2 unpacks or
// SSSE3 shuffles). The AVX2 version is in arch/x86_64/AVX2/search7_avx2.cc.
// [moved from search7.cc, formerly compiled twice, with and without
// SWIPE_SSSE3]

#ifndef SWIPE_SEARCH7_SHARED_H
#define SWIPE_SEARCH7_SHARED_H

#include "swipe.h"  // BYTE, channels_7, vector_bytes, db_getsequence, entry_seqno, ...
#include "arch/intrinsics_to_functions.h"  // v128_t, v_load, v_store, v_dup_i8
#include "align_cells.h"  // Ops_7, align_cells(), GapVectors, No_mask, Mask
#include "view.h"  // View
#include <array>
#include <cstddef>  // std::size_t
#include <cstring>  // std::memset
#include <iterator>  // std::next


constexpr std::size_t channels_7_sse = channels_7;
static_assert(sizeof(v128_t) == vector_bytes, "an SSE vector");
// the database residues of a channel in a block
constexpr std::size_t depth_7 = 4;

// the byte 0x80 (the lanes of _mm_set_epi8() are char: 0x80 does not
// fit in a signed char, -128 has the same bits)
constexpr char byte_0x80 = static_cast<char>(-128);

// the zero score of a 7-bit lane (the scores are stored with a bias of
// 0x80), which is also the lane value that resets a score to zero
// (Ops_7::mask(), align_cells.h)
constexpr BYTE zero_score_7 = 0x80;


// The search of the sequences seqnos by blocks of depth_7 residues of
// each of the 16 channels: build_profile(dprofile, score_matrix, dseq)
// fills the score profile of a block
template <typename BuildProfile>
inline auto search7_blocks(BuildProfile const & build_profile,
                           BYTE * * q_start,
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
  auto * const hep = reinterpret_cast<v128_t*>(hearray);
  auto const * const * const qp = reinterpret_cast<v128_t const * const *>(q_start);
  std::array<BYTE const *, channels_7_sse> d_begin {{}};
  std::array<BYTE const *, channels_7_sse> d_end {{}};
  
  // the database residues of the channels, 16-byte aligned for the loads
  alignas(v128_t) std::array<BYTE, depth_7 * channels_7_sse> dseq {{}};
  // the lanes of the channels that start a new sequence (0x80), and
  // the scores, read lane by lane
  alignas(v128_t) std::array<BYTE, channels_7_sse> restart {{}};
  alignas(v128_t) std::array<BYTE, channels_7_sse> lane_scores {{}};
  BYTE const zero = 0;

  std::array<long, channels_7_sse> seq_id {{}};
  long next_id = 0;
  long done = 0;
  
  std::memset(hearray, zero_score_7, static_cast<std::size_t>(qlen) * hearray_row_bytes);

  auto const Z = v_dup_i8(byte_0x80);
  GapVectors<Ops_7> const gaps {v_dup_i8(static_cast<char>(gap_open_penalty)),
                                v_dup_i8(static_cast<char>(gap_extend_penalty))};

  auto S = Z;

  d_begin.fill(&zero);
  d_end.fill(&zero);
  seq_id.fill(-1);

  // the next depth_7 residues of channel c (zeros after its end);
  // false when the channel reached the end of its sequence
  auto const fill_channel = [&](std::size_t const c) -> bool
  {
    for (std::size_t j = 0; j < depth_7; j++)
    {
      if (d_begin[c] < d_end[c])
      {
        dseq[(channels_7_sse * j) + c] = *d_begin[c];
        d_begin[c] = std::next(d_begin[c]);
      }
      else
      {
        dseq[(channels_7_sse * j) + c] = 0;
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

      for (std::size_t c = 0; c < channels_7_sse; c++)
      {
        if (not fill_channel(c))
        {
          easy = false;
        }
      }

      build_profile(dprofile, score_matrix, dseq.data());

      align_cells<Ops_7>(S, hep, qp, gaps, qlen, Z, No_mask{});
      continue;
    }

    // One or more sequences ended in the previous block 
    // We have to switch over to a new sequence

    easy = true;
    restart.fill(0);
    v_store(reinterpret_cast<v128_t *>(lane_scores.data()), S);

    for (std::size_t c = 0; c < channels_7_sse; c++)
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

      restart[c] = zero_score_7;

      long const cand_id = seq_id[c];
		  
      if (cand_id >= 0)
      {
        // save score
        scores[cand_id] = lane_scores[c] - zero_score_7;
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

    build_profile(dprofile, score_matrix, dseq.data());

    auto const M = v_load(reinterpret_cast<v128_t const *>(restart.data()));
    align_cells<Ops_7>(S, hep, qp, gaps, qlen, Z, Mask<Ops_7>{M});
  }
}

#endif  // SWIPE_SEARCH7_SHARED_H
