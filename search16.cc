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
#include "align_cells.h"  // Ops_16, align_cells(), No_mask, Mask
#include "view.h"  // View
#include <algorithm>  // std::fill_n
#include <array>
#include <cstddef>  // std::size_t
#include <iterator>  // std::next

constexpr std::size_t CHANNELS = channels_16;
constexpr std::size_t CDEPTH = 4;

// the word 0x8000 (the lanes of _mm_set_epi16() are short: 0x8000
// does not fit in a signed short, -32768 has the same bits)
constexpr short word_0x8000 = static_cast<short>(-32768);
constexpr WORD lane_restart = 0x8000;

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

inline auto dprofile_fill16(WORD * dprofile_word,
			    WORD * score_matrix_word,
			    BYTE const * dseq) -> void
{
  __m128i xmm0;
  __m128i xmm1;
  __m128i xmm2;
  __m128i xmm3;
  __m128i xmm4;
  __m128i xmm5;
  __m128i xmm6;
  __m128i xmm7;
  __m128i xmm8;
  __m128i xmm9;
  __m128i xmm10;
  __m128i xmm11;
  __m128i xmm12;
  __m128i xmm13;
  __m128i xmm14;
  __m128i xmm15;
  __m128i xmm16;
  __m128i xmm17;
  __m128i xmm18;
  __m128i xmm19;
  __m128i xmm20;
  __m128i xmm21;
  __m128i xmm22;
  __m128i xmm23;
  __m128i xmm24;
  __m128i xmm25;
  __m128i xmm26;
  __m128i xmm27;
  __m128i xmm28;
  __m128i xmm29;
  __m128i xmm30;
  __m128i xmm31;
  
  for (std::size_t j = 0; j < CDEPTH; j++)
  {
    std::array<int, CHANNELS> d;
    for (std::size_t z = 0; z < CHANNELS; z++)
    {
      d[z] = dseq[(j * CHANNELS) + z] << 5;
    }

    //      for(int i=0; i<24; i += 8)
    for(std::size_t i=0; i<32; i += 8)
    {
      xmm0  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[0] + i));
      xmm1  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[1] + i));
      xmm2  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[2] + i));
      xmm3  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[3] + i));
      xmm4  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[4] + i));
      xmm5  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[5] + i));
      xmm6  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[6] + i));
      xmm7  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[7] + i));
      
      xmm8  = v_merge_lo_16(xmm0,  xmm1);
      xmm9  = v_merge_hi_16(xmm0,  xmm1);
      xmm10 = v_merge_lo_16(xmm2,  xmm3);
      xmm11 = v_merge_hi_16(xmm2,  xmm3);
      xmm12 = v_merge_lo_16(xmm4,  xmm5);
      xmm13 = v_merge_hi_16(xmm4,  xmm5);
      xmm14 = v_merge_lo_16(xmm6,  xmm7);
      xmm15 = v_merge_hi_16(xmm6,  xmm7);
      
      xmm16 = v_merge_lo_32(xmm8,  xmm10);
      xmm17 = v_merge_hi_32(xmm8,  xmm10);
      xmm18 = v_merge_lo_32(xmm12, xmm14);
      xmm19 = v_merge_hi_32(xmm12, xmm14);
      xmm20 = v_merge_lo_32(xmm9,  xmm11);
      xmm21 = v_merge_hi_32(xmm9,  xmm11);
      xmm22 = v_merge_lo_32(xmm13, xmm15);
      xmm23 = v_merge_hi_32(xmm13, xmm15);
      
      xmm24 = v_merge_lo_64(xmm16, xmm18);
      xmm25 = v_merge_hi_64(xmm16, xmm18);
      xmm26 = v_merge_lo_64(xmm17, xmm19);
      xmm27 = v_merge_hi_64(xmm17, xmm19);
      xmm28 = v_merge_lo_64(xmm20, xmm22);
      xmm29 = v_merge_hi_64(xmm20, xmm22);
      xmm30 = v_merge_lo_64(xmm21, xmm23);
      xmm31 = v_merge_hi_64(xmm21, xmm23);
      
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+0)) + (CHANNELS*j)), xmm24);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+1)) + (CHANNELS*j)), xmm25);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+2)) + (CHANNELS*j)), xmm26);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+3)) + (CHANNELS*j)), xmm27);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+4)) + (CHANNELS*j)), xmm28);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+5)) + (CHANNELS*j)), xmm29);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+6)) + (CHANNELS*j)), xmm30);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+7)) + (CHANNELS*j)), xmm31);
    }
  }
}

// the score of each channel, read lane by lane
struct LaneScores
{
  alignas(__m128i) std::array<WORD, CHANNELS> lanes {{}};

  explicit LaneScores(__m128i const scores)
  {
    v_store(reinterpret_cast<__m128i *>(lanes.data()), scores);
  }
};

// the lanes of the channels that start a new sequence (0x8000)
struct RestartLanes
{
  alignas(__m128i) std::array<WORD, CHANNELS> lanes {{}};

  auto vector() const -> __m128i
  {
    return v_load(reinterpret_cast<__m128i const *>(lanes.data()));
  }
};

// lane c of a comparison mask (two bits per 16-bit lane)
inline auto lane_is_set(unsigned int const mask, std::size_t const c) -> bool
{
  return (mask & (3U << (2 * c))) != 0;
}

}  // anonymous namespace

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
  auto * const * const qp = reinterpret_cast<__m128i**>(q_start);
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
      if (lane_is_set(mask, c))
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

      dprofile_fill16(dprofile, score_matrix, dseq.data());
      align_cells<Ops_16>(S, hep, qp, gaps, qlen, Z, No_mask{});
      save_best_positions(v_mask_gt_i16(S, SL));
      SL = S;
      continue;
    }

    easy = true;
    RestartLanes restart;
    LaneScores const lane_scores(S);

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

      restart.lanes[c] = lane_restart;

      long const cand_id = seq_id[c];
      if (cand_id >= 0)
      {
        long const score = lane_scores.lanes[c] ^ 0x8000;
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
    dprofile_fill16(dprofile, score_matrix, dseq.data());
    align_cells<Ops_16>(S, hep, qp, gaps, qlen, Z, Mask<Ops_16>{M});

    /* save column address if new highscore */
    SL = v_adds_i16(SL, M);
    SL = v_adds_i16(SL, M);
    save_best_positions(v_mask_gt_i16(S, SL));
    SL = S;
  }
}
