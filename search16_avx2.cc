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

// The 16-bit kernels with AVX2 (16 database sequences at once,
// __m256i), compiled with -mavx2 and selected at run time
// (cpu_features.avx2): search16_avx2() and search16s_avx2(), the same
// algorithms as search16() and search16s(), twice as many channels.

#include "swipe.h"
#include "intrinsics_to_functions.h"  // v256_load, v256_store, v_merge_*, ...
#include "align_cells.h"  // Ops_16_avx2, align_cells(), align_cells_single()
#include <algorithm>  // std::fill_n
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator>  // std::next

#ifndef __AVX2__
#error "search16_avx2.cc must be compiled with -mavx2"
#endif

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

constexpr std::size_t CHANNELS = channels_16_avx2;

static_assert(CHANNELS * sizeof(WORD) == sizeof(__m256i), "one word lane per channel");

// the word 0x8000 (-32768 has the same bits)
constexpr short word_0x8000 = static_cast<short>(-32768);
constexpr WORD lane_restart = 0x8000;

// The score profile of a block of depth x CHANNELS database residues:
// for each of the 32 query symbols, depth vectors (one per database
// position) of the scores of that query symbol against the residue of
// every channel. The score matrix rows of the database residues are
// transposed by 8 x 8 blocks, as in dprofile_fill16(): the unpacks of
// AVX2 work within each 128-bit half, so a vector holding the row of
// channel c in its low half and the row of channel c + 8 in its high
// half transposes both blocks at once, into the 16 lanes of a profile
// vector.
template <std::size_t depth>
inline auto dprofile_fill16(WORD * dprofile_word,
                            WORD const * score_matrix_word,
                            BYTE const * dseq) -> void
{
  constexpr std::size_t half = CHANNELS / 2;
  constexpr std::size_t row_words = score_matrix_width;
  auto * const profile = reinterpret_cast<__m256i *>(dprofile_word);

  for (std::size_t j = 0; j < depth; j++)
  {
    std::array<std::ptrdiff_t, CHANNELS> d;
    for (std::size_t c = 0; c < CHANNELS; c++)
    {
      d[c] = static_cast<std::ptrdiff_t>(dseq[(j * CHANNELS) + c] * row_words);
    }

    for (std::size_t i = 0; i < score_matrix_width; i += half)
    {
      // eight query symbols (i to i + 7) of the row of channel c
      auto const row = [&](std::size_t const c) -> __m128i const *
      {
        return reinterpret_cast<__m128i const *>(
          std::next(score_matrix_word, d[c] + static_cast<std::ptrdiff_t>(i)));
      };
      auto const in0 = v256_load_halves(row(0), row(half + 0));
      auto const in1 = v256_load_halves(row(1), row(half + 1));
      auto const in2 = v256_load_halves(row(2), row(half + 2));
      auto const in3 = v256_load_halves(row(3), row(half + 3));
      auto const in4 = v256_load_halves(row(4), row(half + 4));
      auto const in5 = v256_load_halves(row(5), row(half + 5));
      auto const in6 = v256_load_halves(row(6), row(half + 6));
      auto const in7 = v256_load_halves(row(7), row(half + 7));

      auto const a0 = v256_merge_lo_16(in0, in1);
      auto const a1 = v256_merge_hi_16(in0, in1);
      auto const a2 = v256_merge_lo_16(in2, in3);
      auto const a3 = v256_merge_hi_16(in2, in3);
      auto const a4 = v256_merge_lo_16(in4, in5);
      auto const a5 = v256_merge_hi_16(in4, in5);
      auto const a6 = v256_merge_lo_16(in6, in7);
      auto const a7 = v256_merge_hi_16(in6, in7);

      auto const b0 = v256_merge_lo_32(a0, a2);
      auto const b1 = v256_merge_hi_32(a0, a2);
      auto const b2 = v256_merge_lo_32(a4, a6);
      auto const b3 = v256_merge_hi_32(a4, a6);
      auto const b4 = v256_merge_lo_32(a1, a3);
      auto const b5 = v256_merge_hi_32(a1, a3);
      auto const b6 = v256_merge_lo_32(a5, a7);
      auto const b7 = v256_merge_hi_32(a5, a7);

      // profile vector of query symbol s, database position j
      auto const store = [&](std::size_t const s, __m256i const vector) -> void
      {
        v256_store(std::next(profile, static_cast<std::ptrdiff_t>((depth * s) + j)), vector);
      };
      store(i + 0, v256_merge_lo_64(b0, b2));
      store(i + 1, v256_merge_hi_64(b0, b2));
      store(i + 2, v256_merge_lo_64(b1, b3));
      store(i + 3, v256_merge_hi_64(b1, b3));
      store(i + 4, v256_merge_lo_64(b4, b6));
      store(i + 5, v256_merge_hi_64(b4, b6));
      store(i + 6, v256_merge_lo_64(b5, b7));
      store(i + 7, v256_merge_hi_64(b5, b7));
    }
  }
}

// the score of each channel, read lane by lane
struct LaneScores
{
  alignas(__m256i) std::array<WORD, CHANNELS> lanes {{}};

  explicit LaneScores(__m256i const scores)
  {
    v256_store(reinterpret_cast<__m256i *>(lanes.data()), scores);
  }
};

// the lanes of the channels that start a new sequence (0x8000)
struct RestartLanes
{
  alignas(__m256i) std::array<WORD, CHANNELS> lanes {{}};

  auto vector() const -> __m256i
  {
    return v256_load(reinterpret_cast<__m256i const *>(lanes.data()));
  }
};

// lane c of a comparison mask (two bits per 16-bit lane)
inline auto lane_is_set(unsigned int const mask, std::size_t const c) -> bool
{
  return (mask & (3U << (2 * c))) != 0;
}

}  // anonymous namespace


auto search16_avx2(WORD * * q_start,
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
  constexpr std::size_t CDEPTH = 4;
  auto * const hep = reinterpret_cast<__m256i*>(hearray);
  auto ** const qp = reinterpret_cast<__m256i**>(q_start);
  std::array<BYTE const *, CHANNELS> d_begin;
  std::array<BYTE const *, CHANNELS> d_pos;
  std::array<BYTE const *, CHANNELS> d_best;
  std::array<BYTE const *, CHANNELS> d_end;

  // the database residues of the channels, aligned for the loads
  alignas(__m256i) std::array<BYTE, CDEPTH * CHANNELS> dseq {{}};
  BYTE const zero = 0;

  std::array<long, CHANNELS> seq_id;
  long next_id = 0;
  long done = 0;

  auto const Z = v256_dup_i16(word_0x8000);
  auto const Q = v256_dup_i16(static_cast<short>(gap_open_penalty));
  auto const R = v256_dup_i16(static_cast<short>(gap_extend_penalty));

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

  // save column address if new highscore
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

      dprofile_fill16<CDEPTH>(dprofile, score_matrix, dseq.data());
      align_cells<Ops_16_avx2>(S, hep, qp, Q, R, qlen, Z, No_mask{});
      save_best_positions(v256_mask_gt_i16(S, SL));
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
    dprofile_fill16<CDEPTH>(dprofile, score_matrix, dseq.data());
    align_cells<Ops_16_avx2>(S, hep, qp, Q, R, qlen, Z, Mask<Ops_16_avx2>{M});

    SL = v256_adds_i16(SL, M);
    SL = v256_adds_i16(SL, M);
    save_best_positions(v256_mask_gt_i16(S, SL));
    SL = S;
  }
}


auto search16s_avx2(WORD * * q_start,
                    WORD gap_open_penalty,
                    WORD gap_extend_penalty,
                    WORD * score_matrix,
                    WORD * dprofile,
                    WORD * hearray,
                    DbThread const * dbta,
                    long sequences,
                    long const * seqnos,
                    long * scores,
                    long * bestpos,
                    long * bestq,
                    int qlen) -> void
{
  constexpr std::size_t CDEPTH = 1;
  auto * const hep = reinterpret_cast<__m256i*>(hearray);
  auto ** const qp = reinterpret_cast<__m256i**>(q_start);
  std::array<BYTE const *, CHANNELS> d_begin;
  std::array<BYTE const *, CHANNELS> d_pos;
  std::array<BYTE const *, CHANNELS> d_end;
  std::array<BYTE const *, CHANNELS> d_best;
  std::array<long, CHANNELS> q_best;

  // the database residues of the channels, aligned for the loads
  alignas(__m256i) std::array<BYTE, CDEPTH * CHANNELS> dseq {{}};
  BYTE const zero = 0;

  std::array<long, CHANNELS> seq_id;
  long next_id = 0;
  long done = 0;

  auto const Z = v256_dup_i16(word_0x8000);
  auto const Q = v256_dup_i16(static_cast<short>(gap_open_penalty));
  auto const R = v256_dup_i16(static_cast<short>(gap_extend_penalty));

  auto S = Z;
  auto SL = Z;

  // the H and E scores of each query residue
  std::fill_n(hep, 2 * static_cast<std::size_t>(qlen), Z);

  d_begin.fill(&zero);
  d_pos.fill(&zero);
  d_best.fill(&zero);
  d_end.fill(&zero);
  q_best.fill(-1);
  seq_id.fill(-1);

  // the next residue of channel c (zero after its end)
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

  // save column address if new highscore, and the query position of
  // that score (the first one, scanning backwards)
  auto const save_best_positions = [&](unsigned int const mask) -> void
  {
    if (mask == 0)
    {
      return;
    }
    for (std::size_t c = 0; c < CHANNELS; c++)
    {
      if (lane_is_set(mask, c))
      {
        d_best[c] = std::prev(d_pos[c]);
      }
    }
    for (long i = qlen - 1; i >= 0; i--)
    {
      auto const m2 = mask & v256_mask_eq_i16(hep[2 * i], S);
      if (m2 == 0)
      {
        continue;
      }
      for (std::size_t c = 0; c < CHANNELS; c++)
      {
        if (lane_is_set(m2, c))
        {
          q_best[c] = i;
        }
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
      align_cells_single<Ops_16_avx2>(S, hep, qp, Q, R, qlen, Z, No_mask{});
      save_best_positions(v256_mask_gt_i16(S, SL));
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
        bestq[cand_id] = q_best[c];
        done++;
      }

      if (next_id < sequences)
      {
        seq_id[c] = next_id;
        long const seqnosf = seqnos[next_id];
        long const seqno = entry_seqno(seqnosf);
        db_mapsequences(*dbta[c], seqno, seqno);
        View<char> const sequence =
          db_getsequence(*dbta[c], seqno, entry_where(seqnosf), c);
        d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
        d_pos[c] = d_begin[c];
        d_best[c] = d_begin[c];
        d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
        q_best[c] = -1;
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
        d_pos[c] = &zero;
        d_end[c] = d_pos[c];
        fill_channel(c);
      }
    }

    if (done == sequences)
    {
      break;
    }

    auto const M = restart.vector();
    dprofile_fill16<CDEPTH>(dprofile, score_matrix, dseq.data());
    align_cells_single<Ops_16_avx2>(S, hep, qp, Q, R, qlen, Z, Mask<Ops_16_avx2>{M});

    SL = v256_adds_i16(SL, M);
    SL = v256_adds_i16(SL, M);
    save_best_positions(v256_mask_gt_i16(S, SL));
    SL = S;
  }
}
