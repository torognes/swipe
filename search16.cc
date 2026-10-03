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
#include "intrinsics_to_functions.h"  // v_load, v_store, v_merge_*, v_dup_*, ...
#include "search16_shared.h"  // dprofile_fill16(), LaneScores16, RestartLanes16
#include "align_cells.h"  // Ops_16, align_cells(), No_mask, Mask
#include "view.h"  // View
#include <algorithm>  // std::fill_n
#include <array>
#include <cstddef>  // std::size_t
#include <iterator>  // std::next

constexpr std::size_t CHANNELS = channels_16_sse;
constexpr std::size_t CDEPTH = 4;

// the word 0x8000 (the lanes of _mm_set_epi16() are short: 0x8000
// does not fit in a signed short, -32768 has the same bits)
constexpr short word_0x8000 = static_cast<short>(-32768);
// the zero score of a 16-bit lane (the scores are stored with a bias
// of 0x8000), which is also the lane value that resets a score to zero
// (Ops_16::mask(), align_cells.h)
constexpr WORD zero_score = 0x8000;

auto search16(WORD * * q_start,
	      WORD gap_open_penalty,
	      WORD gap_extend_penalty,
	      WORD * score_matrix,
	      WORD * dprofile,
	      WORD * hearray,
	      db_thread_s & dbt,
	      long sequences,
	      long const * seqnos,
	      long * scores,
	      long * bestpos,
	      int qlen) -> void
{
  auto * const hep = reinterpret_cast<__m128i*>(hearray);
  auto const * const * const qp = reinterpret_cast<__m128i const * const *>(q_start);
  std::array<BYTE const *, CHANNELS> d_begin {{}};
  std::array<BYTE const *, CHANNELS> d_pos {{}};
  std::array<BYTE const *, CHANNELS> d_best {{}};
  std::array<BYTE const *, CHANNELS> d_end {{}};

  // the database residues of the channels, 16-byte aligned for the loads
  alignas(__m128i) std::array<BYTE, CDEPTH * CHANNELS> dseq {{}};
  BYTE const zero = 0;

  std::array<long, CHANNELS> seq_id {{}};
  long next_id = 0;
  long done = 0;

  auto const Z = v_dup_i16(word_0x8000);
  GapVectors<Ops_16> const gaps {v_dup_i16(static_cast<short>(gap_open_penalty)),
                                 v_dup_i16(static_cast<short>(gap_extend_penalty))};

  auto S = Z;
  auto SL = Z;

  // the H and E scores of each query residue
  std::fill_n(hep, 2 * static_cast<std::size_t>(qlen), Z);

  d_begin.fill(&zero);
  d_pos.fill(&zero);
  d_best.fill(&zero);
  d_end.fill(&zero);
  seq_id.fill(-1);

  // the next CDEPTH residues of channel c (zeros after its end)
  auto const fill_channel = [&](std::size_t const c) -> void
  {
    for (std::size_t j = 0; j < CDEPTH; j++)
    {
      if (d_pos[c] < d_end[c])
      {
        dseq[(CHANNELS * j) + c] = *d_pos[c];
        d_pos[c] = std::next(d_pos[c]);
      }
      else
      {
        dseq[(CHANNELS * j) + c] = 0;
      }
    }
  };

  /* save column address if new highscore */
  auto const save_best_positions = [&](unsigned int const mask) -> void
  {
    for (std::size_t c = 0; c < CHANNELS; c++)
    {
      if (lane16_is_set(mask, c))
      {
        d_best[c] = d_pos[c];
      }
    }
  };

  bool easy = false;

  while (true)
  {
    if (easy)
    {
      for (std::size_t c = 0; c < CHANNELS; c++)
      {
        fill_channel(c);
        if ((d_pos[c] == d_end[c]) && (seq_id[c] > -1))
        {
          easy = false;
        }
      }

      dprofile_fill16<CDEPTH>(dprofile, score_matrix, dseq.data());
      align_cells<Ops_16>(S, hep, qp, gaps, qlen, Z, No_mask{});
      save_best_positions(v_mask_gt_i16(S, SL));
      SL = S;
      continue;
    }

    easy = true;
    RestartLanes16 restart;
    LaneScores16 const lane_scores(S);

    for (std::size_t c = 0; c < CHANNELS; c++)
    {
      if (d_pos[c] < d_end[c])
      {
        fill_channel(c);
        if (d_pos[c] == d_end[c])
        {
          easy = false;
        }
        continue;
      }

      restart.lanes[c] = zero_score;

      long const cand_id = seq_id[c];
      if (cand_id >= 0)
      {
        long const score = lane_scores.lanes[c] ^ zero_score;
        scores[cand_id] = score;
        bestpos[cand_id] = d_best[c] - d_begin[c];
        done++;
      }

      if (next_id < sequences)
      {
        seq_id[c] = next_id;
        long const seqnosf = seqnos[next_id];
        View<char> const sequence =
          db_getsequence(dbt, entry_seqno(seqnosf), entry_where(seqnosf), c);
        d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
        d_pos[c] = d_begin[c];
        d_best[c] = d_begin[c];
        d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
        next_id++;
        fill_channel(c);
        if (d_pos[c] == d_end[c])
        {
          easy = false;
        }
      }
      else
      {
        seq_id[c] = -1;
        d_begin[c] = &zero;
        d_pos[c] = d_begin[c];
        d_best[c] = d_begin[c];
        d_end[c] = d_begin[c];
        fill_channel(c);
      }
    }

    if (done == sequences)
    {
      break;
    }

    auto const M = restart.vector();
    dprofile_fill16<CDEPTH>(dprofile, score_matrix, dseq.data());
    align_cells<Ops_16>(S, hep, qp, gaps, qlen, Z, Mask<Ops_16>{M});

    /* save column address if new highscore */
    SL = v_adds_i16(SL, M);
    SL = v_adds_i16(SL, M);
    save_best_positions(v_mask_gt_i16(S, SL));
    SL = S;
  }
}
