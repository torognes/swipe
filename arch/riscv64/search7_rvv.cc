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

// The 7-bit kernel with RVV: the score profile built by table lookups
// (vrgather) in the transposed score matrix, as search7_ssse3() does
// with pshufb (arch/x86_64/SSSE3/search7_ssse3.cc), and the cells of
// arch/riscv64/align_rvv.cc. The search loop is shared with the other
// 128-bit kernels (search7_shared.h).

#include "swipe.h"
#include "search_kernels.h"  // the prototype of the kernel

#ifdef __riscv_vector

#include "search7_shared.h"  // search7_blocks(), channels_7_sse, depth_7
#include <riscv_vector.h>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator>  // std::next

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

constexpr std::size_t lanes = channels_7_sse;
constexpr std::size_t row_bytes = score_matrix_width;  // the table of a lookup
constexpr std::size_t profile_row_bytes = depth_7 * channels_7_sse;
constexpr auto row_stride = static_cast<std::ptrdiff_t>(row_bytes);
constexpr auto profile_row_stride = static_cast<std::ptrdiff_t>(profile_row_bytes);

// The profile vector of query symbol 'row' and database position j is
// the row of the transposed score matrix (the scores of that query
// symbol against the 32 database symbols), indexed by the database
// symbols of position j. LMUL = 2: a table of 32 bytes needs two
// registers when VLEN = 128 (the minimum of the V extension); vl is 32
// for the table, 16 (the channels) for the lookups.
inline auto dprofile_gather7(BYTE * const dprofile,
                             BYTE const * const score_matrix,
                             BYTE const * const dseq) -> void
{
  static_assert(depth_7 == 4, "four database positions per block");
  auto const position_0 = __riscv_vle8_v_u8m2(dseq, lanes);
  auto const position_1 = __riscv_vle8_v_u8m2(std::next(dseq, lanes), lanes);
  auto const position_2 = __riscv_vle8_v_u8m2(std::next(dseq, 2 * lanes), lanes);
  auto const position_3 = __riscv_vle8_v_u8m2(std::next(dseq, 3 * lanes), lanes);

  for (std::ptrdiff_t row = 0; row < static_cast<std::ptrdiff_t>(score_matrix_width); ++row)
  {
    auto const table = __riscv_vle8_v_u8m2(std::next(score_matrix, row * row_stride), row_bytes);
    auto * const line = std::next(dprofile, row * profile_row_stride);
    __riscv_vse8_v_u8m2(line, __riscv_vrgather_vv_u8m2(table, position_0, lanes), lanes);
    __riscv_vse8_v_u8m2(std::next(line, lanes), __riscv_vrgather_vv_u8m2(table, position_1, lanes), lanes);
    __riscv_vse8_v_u8m2(std::next(line, 2 * lanes), __riscv_vrgather_vv_u8m2(table, position_2, lanes), lanes);
    __riscv_vse8_v_u8m2(std::next(line, 3 * lanes), __riscv_vrgather_vv_u8m2(table, position_3, lanes), lanes);
  }
}

}  // anonymous namespace

auto search7_rvv(BYTE * * const q_start,
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
    dprofile_gather7(profile, matrix, dseq);
  };
  search7_blocks(build_profile, q_start, gap_open_penalty, gap_extend_penalty,
                 score_matrix, dprofile, hearray, dbt, sequences, seqnos,
                 scores, qlen);
}

#endif  // __riscv_vector
