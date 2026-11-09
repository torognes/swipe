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
#include "search_kernels.h"  // search7, search16, search16s
#include "swipe.h"  // BYTE, WORD, ScoreMatrices, vector_bytes
#include <cassert>
#include <cstddef>  // std::size_t


// Any machine other than x86-64: one set of kernels, with 16-byte
// vectors (search7() builds its profile from the score matrix as is,
// not transposed)

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
  return vector_bytes;
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
  assert(kernel_bytes == vector_bytes);
  static_cast<void>(kernel_bytes);  // without assert() (NDEBUG)
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
  assert(kernel_bytes == vector_bytes);
  static_cast<void>(kernel_bytes);  // without assert() (NDEBUG)
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
  assert(kernel_bytes == vector_bytes);
  static_cast<void>(kernel_bytes);  // without assert() (NDEBUG)
  search16s(q_start, gap_open_penalty, gap_extend_penalty, matrix_16(matrices),
            dprofile, hearray, dbta, sequences, seqnos, scores, bestpos, bestq, qlen);
}
