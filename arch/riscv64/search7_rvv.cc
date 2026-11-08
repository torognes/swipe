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
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator>  // std::next

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

constexpr std::size_t row_bytes = score_matrix_width;  // the table of a lookup
// the profile vectors of a query symbol (its row of the profile), and
// the database symbols of a block: depth_7 positions of 16 channels
constexpr std::size_t block_bytes = depth_7 * channels_7_sse;
constexpr auto row_stride = static_cast<std::ptrdiff_t>(row_bytes);
constexpr auto block_stride = static_cast<std::ptrdiff_t>(block_bytes);
// register group sizes (VLEN is at least 128 bits: LMUL = 4 holds 64
// bytes)
constexpr std::size_t lmul_2 = 2;
constexpr std::size_t lmul_4 = 4;

// The RVV operations of the profile builder, at one register group
// size (LMUL = 1, 2 or 4)
struct Lmul_1 {
  static auto load(BYTE const * const source, std::size_t const length) -> vuint8m1_t {
    return __riscv_vle8_v_u8m1(source, length);
  }
  static auto lookup(vuint8m1_t const table, vuint8m1_t const indices) -> vuint8m1_t {
    return __riscv_vrgather_vv_u8m1(table, indices, block_bytes);
  }
  static auto store(BYTE * const target, vuint8m1_t const vector) -> void {
    __riscv_vse8_v_u8m1(target, vector, block_bytes);
  }
};

struct Lmul_2 {
  static auto load(BYTE const * const source, std::size_t const length) -> vuint8m2_t {
    return __riscv_vle8_v_u8m2(source, length);
  }
  static auto lookup(vuint8m2_t const table, vuint8m2_t const indices) -> vuint8m2_t {
    return __riscv_vrgather_vv_u8m2(table, indices, block_bytes);
  }
  static auto store(BYTE * const target, vuint8m2_t const vector) -> void {
    __riscv_vse8_v_u8m2(target, vector, block_bytes);
  }
};

struct Lmul_4 {
  static auto load(BYTE const * const source, std::size_t const length) -> vuint8m4_t {
    return __riscv_vle8_v_u8m4(source, length);
  }
  static auto lookup(vuint8m4_t const table, vuint8m4_t const indices) -> vuint8m4_t {
    return __riscv_vrgather_vv_u8m4(table, indices, block_bytes);
  }
  static auto store(BYTE * const target, vuint8m4_t const vector) -> void {
    __riscv_vse8_v_u8m4(target, vector, block_bytes);
  }
};

// The profile vectors of query symbol 'row' (one per database position
// of the block, 64 bytes in all) are the row of the transposed score
// matrix (the scores of that query symbol against the 32 database
// symbols), indexed by the 64 database symbols of the block (dseq, in
// the same order): one lookup per row. vl is 32 for the table, 64 for
// the lookup.
template <typename Lmul>
inline auto dprofile_gather7(BYTE * const dprofile,
                             BYTE const * const score_matrix,
                             BYTE const * const dseq) -> void
{
  auto const symbols = Lmul::load(dseq, block_bytes);
  for (std::ptrdiff_t row = 0; row < static_cast<std::ptrdiff_t>(score_matrix_width); ++row)
  {
    auto const table = Lmul::load(std::next(score_matrix, row * row_stride), row_bytes);
    Lmul::store(std::next(dprofile, row * block_stride), Lmul::lookup(table, symbols));
  }
}

// The smallest register group that holds the 64 bytes of a block: a
// lookup costs as much as its register group, not as its vl (BPI-F3,
// VLEN = 256: LMUL = 2 is 3.2 times faster than LMUL = 4, for the same
// result)
inline auto dprofile_build7(BYTE * const dprofile,
                            BYTE const * const score_matrix,
                            BYTE const * const dseq) -> void
{
  static_assert(block_bytes == 64, "four positions of 16 channels");
  auto const register_bytes = __riscv_vlenb();  // VLEN / 8, at least 16
  if (register_bytes >= block_bytes)
  {
    dprofile_gather7<Lmul_1>(dprofile, score_matrix, dseq);
  }
  else if (lmul_2 * register_bytes >= block_bytes)
  {
    dprofile_gather7<Lmul_2>(dprofile, score_matrix, dseq);
  }
  else
  {
    assert(lmul_4 * register_bytes >= block_bytes);
    dprofile_gather7<Lmul_4>(dprofile, score_matrix, dseq);
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
    dprofile_build7(profile, matrix, dseq);
  };
  search7_blocks(build_profile, q_start, gap_open_penalty, gap_extend_penalty,
                 score_matrix, dprofile, hearray, dbt, sequences, seqnos,
                 scores, qlen);
}

#endif  // __riscv_vector
