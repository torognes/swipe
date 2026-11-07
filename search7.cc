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

// The 7-bit kernel with SSE2 only: the score profile built by
// transposing the score matrix rows (unpacks). The search loop is
// shared with the SSSE3 version (search7_shared.h).

#include "swipe.h"
#include "search_kernels.h"  // the prototype of the kernel
#include "arch/intrinsics_to_functions.h"  // v128_t, v_load, v_load_64, v_store, v_merge_*
#include "search7_shared.h"  // search7_blocks(), channels_7_sse, depth_7
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator>  // std::next

constexpr std::size_t CHANNELS = channels_7_sse;
constexpr std::size_t CDEPTH = depth_7;

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

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

}  // anonymous namespace

auto search7(BYTE * * const q_start,
	BYTE const gap_open_penalty,
	BYTE const gap_extend_penalty,
	BYTE const * const score_matrix,
	BYTE * const dprofile,
	BYTE * const hearray,
	db_thread_s & dbt,
	long const sequences,
	long const * const seqnos,
	long * const scores,
	long const qlen) -> void
{
  auto const build_profile = [](BYTE * const profile,
                                BYTE const * const matrix,
                                BYTE const * const dseq) -> void
  {
    dprofile_fill7(profile, matrix, dseq);
  };
  search7_blocks(build_profile, q_start, gap_open_penalty, gap_extend_penalty,
                 score_matrix, dprofile, hearray, dbt, sequences, seqnos,
                 scores, qlen);
}
