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

#include "arch/search_dispatch.h"
#include "search_kernels.h"  // search7, search7_rvv, search7_rvv256, search16, ...
#include "swipe.h"  // BYTE, WORD, ScoreMatrices, vector_bytes, vector_bytes_avx2
#include <cstddef>  // std::size_t
#ifdef __riscv_vector
#include <riscv_vector.h>
#endif


// riscv64: with the V extension, the RVV kernels; with 32-byte vectors
// when VLEN >= 256 (one register holds them: the search7_rvv256(),
// search16_rvv256() and search16s_rvv256() kernels), with 16-byte
// vectors otherwise (search7_rvv(), and search16() and search16s()
// with the RVV cells). Without V (make SWIPE_RISCV_BASELINE=1): the
// kernels of the other machines (arch/generic/search_dispatch.cc).
// Build-time test switch: -DSWIPE_NO_RVV256 (the 16-byte kernels
// whatever VLEN, as -DSWIPE_NO_AVX2 on x86-64).

namespace {

// the 16-bit score matrix, as the 16-bit kernels read it
auto matrix_16(ScoreMatrices const & matrices) noexcept -> WORD const *
{
  return reinterpret_cast<WORD const *>(matrices.score_16.data());
}

// the 32-byte kernels: VLEN >= 256 (vlenb: VLEN / 8, in bytes)
auto wide_vectors() noexcept -> bool
{
#if defined(__riscv_vector) && not defined(SWIPE_NO_RVV256)
  return __riscv_vlenb() >= vector_bytes_avx2;
#else
  return false;
#endif
}

}  // anonymous namespace


auto kernel_vector_bytes() noexcept -> std::size_t
{
  return wide_vectors() ? vector_bytes_avx2 : vector_bytes;
}


auto search7_dispatch(BYTE * * const q_start,
                      BYTE const gap_open_penalty,
                      BYTE const gap_extend_penalty,
                      ScoreMatrices const & matrices,
                      BYTE * const dprofile,
                      BYTE * const hearray,
                      db_thread_s & dbt,
                      long const sequences,
                      long const * const seqnos,
                      long * const scores,
                      long const qlen) -> void
{
#ifdef __riscv_vector
  auto const * const transposed = reinterpret_cast<BYTE const *>(matrices.score_7t.data());
  if (wide_vectors())
  {
    search7_rvv256(q_start, gap_open_penalty, gap_extend_penalty, transposed,
                   dprofile, hearray, dbt, sequences, seqnos, scores, qlen);
    return;
  }
  search7_rvv(q_start, gap_open_penalty, gap_extend_penalty, transposed,
              dprofile, hearray, dbt, sequences, seqnos, scores, qlen);
#else
  search7(q_start, gap_open_penalty, gap_extend_penalty,
          reinterpret_cast<BYTE const *>(matrices.score_7.data()),
          dprofile, hearray, dbt, sequences, seqnos, scores, qlen);
#endif
}


auto search16_dispatch(WORD * * const q_start,
                       WORD const gap_open_penalty,
                       WORD const gap_extend_penalty,
                       ScoreMatrices const & matrices,
                       WORD * const dprofile,
                       WORD * const hearray,
                       db_thread_s & dbt,
                       long const sequences,
                       long const * const seqnos,
                       long * const scores,
                       long * const bestpos,
                       int const qlen) -> void
{
#ifdef __riscv_vector
  if (wide_vectors())
  {
    search16_rvv256(q_start, gap_open_penalty, gap_extend_penalty, matrix_16(matrices),
                    dprofile, hearray, dbt, sequences, seqnos, scores, bestpos, qlen);
    return;
  }
#endif
  search16(q_start, gap_open_penalty, gap_extend_penalty, matrix_16(matrices),
           dprofile, hearray, dbt, sequences, seqnos, scores, bestpos, qlen);
}


auto search16s_dispatch(WORD * * const q_start,
                        WORD const gap_open_penalty,
                        WORD const gap_extend_penalty,
                        ScoreMatrices const & matrices,
                        WORD * const dprofile,
                        WORD * const hearray,
                        DbThread const * const dbta,
                        long const sequences,
                        long const * const seqnos,
                        long * const scores,
                        long * const bestpos,
                        long * const bestq,
                        int const qlen) -> void
{
#ifdef __riscv_vector
  if (wide_vectors())
  {
    search16s_rvv256(q_start, gap_open_penalty, gap_extend_penalty, matrix_16(matrices),
                     dprofile, hearray, dbta, sequences, seqnos, scores, bestpos, bestq, qlen);
    return;
  }
#endif
  search16s(q_start, gap_open_penalty, gap_extend_penalty, matrix_16(matrices),
            dprofile, hearray, dbta, sequences, seqnos, scores, bestpos, bestq, qlen);
}
