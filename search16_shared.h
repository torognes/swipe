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

// What search16() and search16s() (search16.cc, search16s.cc) share:
// the score profile builder, over the number of database residues per
// channel, and the lane-by-lane access to the scores and the restart
// mask (the 8 word lanes of an SSE vector). The AVX2 versions are in
// search16_avx2.cc.

#ifndef SWIPE_SEARCH16_SHARED_H
#define SWIPE_SEARCH16_SHARED_H

#include "swipe.h"  // BYTE, WORD, channels_16, score_matrix_width
#include "arch/intrinsics_to_functions.h"  // v128_t, v_load, v_store, v_merge_*
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator>  // std::next


constexpr std::size_t channels_16_sse = channels_16;

static_assert(channels_16_sse * sizeof(WORD) == sizeof(v128_t), "one word lane per channel");

// The score profile of a block of depth x channels database residues:
// for each of the 32 query symbols, depth vectors (one per database
// position) of the scores of that query symbol against the residue of
// every channel. The score matrix rows of the database residues are
// transposed by 8 x 8 blocks.
template <std::size_t depth>
inline auto dprofile_fill16(WORD * dprofile_word,
                            WORD const * score_matrix_word,
                            BYTE const * dseq) -> void
{
  constexpr std::size_t row_words = score_matrix_width;
  constexpr std::size_t block = channels_16_sse;  // 8 x 8 words
  auto * const profile = reinterpret_cast<v128_t *>(dprofile_word);

  for (std::size_t j = 0; j < depth; j++)
  {
    std::array<std::ptrdiff_t, channels_16_sse> d {{}};
    for (std::size_t c = 0; c < channels_16_sse; c++)
    {
      d[c] = static_cast<std::ptrdiff_t>(dseq[(j * channels_16_sse) + c] * row_words);
    }

    //      for(int i=0; i<24; i += 8)
    for (std::size_t i = 0; i < score_matrix_width; i += block)
    {
      // eight query symbols (i to i + 7) of the row of channel c
      auto const row = [&](std::size_t const c) -> v128_t
      {
        return v_load(reinterpret_cast<v128_t const *>(
          std::next(score_matrix_word, d[c] + static_cast<std::ptrdiff_t>(i))));
      };
      auto const in0 = row(0);
      auto const in1 = row(1);
      auto const in2 = row(2);
      auto const in3 = row(3);
      auto const in4 = row(4);
      auto const in5 = row(5);
      auto const in6 = row(6);
      auto const in7 = row(7);

      auto const a0 = v_merge_lo_16(in0, in1);
      auto const a1 = v_merge_hi_16(in0, in1);
      auto const a2 = v_merge_lo_16(in2, in3);
      auto const a3 = v_merge_hi_16(in2, in3);
      auto const a4 = v_merge_lo_16(in4, in5);
      auto const a5 = v_merge_hi_16(in4, in5);
      auto const a6 = v_merge_lo_16(in6, in7);
      auto const a7 = v_merge_hi_16(in6, in7);

      auto const b0 = v_merge_lo_32(a0, a2);
      auto const b1 = v_merge_hi_32(a0, a2);
      auto const b2 = v_merge_lo_32(a4, a6);
      auto const b3 = v_merge_hi_32(a4, a6);
      auto const b4 = v_merge_lo_32(a1, a3);
      auto const b5 = v_merge_hi_32(a1, a3);
      auto const b6 = v_merge_lo_32(a5, a7);
      auto const b7 = v_merge_hi_32(a5, a7);

      // profile vector of query symbol s, database position j
      auto const store = [&](std::size_t const s, v128_t const vector) -> void
      {
        v_store(std::next(profile, static_cast<std::ptrdiff_t>((depth * s) + j)), vector);
      };
      store(i + 0, v_merge_lo_64(b0, b2));
      store(i + 1, v_merge_hi_64(b0, b2));
      store(i + 2, v_merge_lo_64(b1, b3));
      store(i + 3, v_merge_hi_64(b1, b3));
      store(i + 4, v_merge_lo_64(b4, b6));
      store(i + 5, v_merge_hi_64(b4, b6));
      store(i + 6, v_merge_lo_64(b5, b7));
      store(i + 7, v_merge_hi_64(b5, b7));
    }
  }
}

// the score of each channel, read lane by lane
struct LaneScores16
{
  alignas(v128_t) std::array<WORD, channels_16_sse> lanes {{}};

  explicit LaneScores16(v128_t const scores)
  {
    v_store(reinterpret_cast<v128_t *>(lanes.data()), scores);
  }
};

// the lanes of the channels that start a new sequence (0x8000)
struct RestartLanes16
{
  alignas(v128_t) std::array<WORD, channels_16_sse> lanes {{}};

  auto vector() const -> v128_t
  {
    return v_load(reinterpret_cast<v128_t const *>(lanes.data()));
  }
};

// lane c of a comparison mask (two bits per 16-bit lane)
inline auto lane16_is_set(unsigned int const mask, std::size_t const c) -> bool
{
  return (mask & (3U << (2 * c))) != 0;
}

#endif  // SWIPE_SEARCH16_SHARED_H
