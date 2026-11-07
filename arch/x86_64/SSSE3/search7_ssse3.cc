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

// The 7-bit kernel with SSSE3: the score profile built by table
// lookups (pshufb) in the transposed score matrix, compiled with
// -mssse3 and selected at run time (arch/x86_64/search_dispatch.cc).
// The search loop is shared with the SSE2 version (search7_shared.h).
// [moved from search7.cc, formerly compiled with -DSWIPE_SSSE3]

#include "swipe.h"
#include "search_kernels.h"  // the prototype of the kernel
#include "arch/x86_64/intrinsics_to_functions.h"  // v128_t, v_load, v_store, v_shuffle_8, ...
#include "search7_shared.h"  // search7_blocks(), depth_7, byte_0x80
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator>  // std::next

constexpr std::size_t CDEPTH = depth_7;

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

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

}  // anonymous namespace

auto search7_ssse3(BYTE * * const q_start,
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
    dprofile_shuffle7(profile, matrix, dseq);
  };
  search7_blocks(build_profile, q_start, gap_open_penalty, gap_extend_penalty,
                 score_matrix, dprofile, hearray, dbt, sequences, seqnos,
                 scores, qlen);
}
