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
#include "arch/x86_64/cpu_features.h"  // cpu_features
#include "search_kernels.h"  // search7, search7_ssse3, search7_avx2, ...
#include "swipe.h"  // BYTE, WORD, ScoreMatrices, vector_bytes, vector_bytes_avx2
#include <cassert>
#include <cstddef>  // std::size_t


// The x86-64 kernels: AVX2 when the processor has it (32-byte
// vectors), else SSSE3 for the 7-bit profile (pshufb, which reads the
// transposed score matrix), else SSE2 (16-byte vectors)

namespace {

// the 16-bit score matrix, as the 16-bit kernels read it
auto matrix_16(ScoreMatrices const & matrices) noexcept -> WORD const *
{
  return reinterpret_cast<WORD const *>(matrices.score_16.data());
}

}  // anonymous namespace


// whatever the query and the threads
auto kernel_vector_bytes(long /* longest_query */, long /* threads */) noexcept -> std::size_t
{
  return cpu_features.avx2 ? vector_bytes_avx2 : vector_bytes;
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
                      long const qlen,
                      std::size_t const kernel_bytes) -> void
{
  if (kernel_bytes == vector_bytes_avx2)
  {
    assert(cpu_features.avx2);
    search7_avx2(q_start, gap_open_penalty, gap_extend_penalty,
                 reinterpret_cast<BYTE const *>(matrices.score_7t.data()),
                 dprofile, hearray, dbt, sequences, seqnos, scores, qlen);
    return;
  }
  if (cpu_features.ssse3)
  {
    search7_ssse3(q_start, gap_open_penalty, gap_extend_penalty,
                  reinterpret_cast<BYTE const *>(matrices.score_7t.data()),
                  dprofile, hearray, dbt, sequences, seqnos, scores, qlen);
    return;
  }
  search7(q_start, gap_open_penalty, gap_extend_penalty,
          reinterpret_cast<BYTE const *>(matrices.score_7.data()),
          dprofile, hearray, dbt, sequences, seqnos, scores, qlen);
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
                       int const qlen,
                       std::size_t const kernel_bytes) -> void
{
  if (kernel_bytes == vector_bytes_avx2)
  {
    assert(cpu_features.avx2);
    search16_avx2(q_start, gap_open_penalty, gap_extend_penalty, matrix_16(matrices),
                  dprofile, hearray, dbt, sequences, seqnos, scores, bestpos, qlen);
    return;
  }
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
                        int const qlen,
                        std::size_t const kernel_bytes) -> void
{
  if (kernel_bytes == vector_bytes_avx2)
  {
    assert(cpu_features.avx2);
    search16s_avx2(q_start, gap_open_penalty, gap_extend_penalty, matrix_16(matrices),
                   dprofile, hearray, dbta, sequences, seqnos, scores, bestpos, bestq, qlen);
    return;
  }
  search16s(q_start, gap_open_penalty, gap_extend_penalty, matrix_16(matrices),
            dprofile, hearray, dbta, sequences, seqnos, scores, bestpos, bestq, qlen);
}
