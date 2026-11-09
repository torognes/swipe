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

// The choice of the search kernels, common interface (as swarm's
// src/arch/*/search_dispatch.hpp): one implementation per
// architecture, chosen by the Makefile (ARCH). arch/x86_64/ picks the
// AVX2, SSSE3 or SSE2 kernels from the CPU features detected at run
// time; arch/generic/ has one set of kernels. The callers
// (search_threads.cc, align_threads.cc) name no instruction set.
// The parameters are those of the kernels (search_kernels.h), except
// the score matrix: the dispatch picks the layout its kernel reads.

#ifndef SWIPE_ARCH_SEARCH_DISPATCH_H
#define SWIPE_ARCH_SEARCH_DISPATCH_H

#include "swipe.h"  // BYTE, WORD, ScoreMatrices, db_thread_s, DbThread
#include <cstddef>  // std::size_t


// the bytes of a vector of the kernels chosen for this processor (16
// or 32): the channels of the 7-bit kernel (a byte each) and of the
// 16-bit kernels (a word each), and the strides of the score profile
// rows and of the H/E array that follow from them. Chosen once per
// query, from the longest strand or frame searched and the number of
// threads (riscv64: see arch/riscv64/search_dispatch.cc), and passed
// to the kernel dispatch below (kernel_bytes)
auto kernel_vector_bytes(long longest_query, long threads) noexcept -> std::size_t;

// 7-bit search, 4 database residues per block
auto search7_dispatch(BYTE * * q_start,
                      BYTE gap_open_penalty,
                      BYTE gap_extend_penalty,
                      ScoreMatrices const & matrices,
                      BYTE * dprofile,
                      BYTE * hearray,
                      db_thread_s & dbt,
                      long sequences,
                      long const * seqnos,
                      long * scores,
                      long qlen,
                      std::size_t kernel_bytes) -> void;

// 16-bit search, 4 database residues per block, with the alignment end
auto search16_dispatch(WORD * * q_start,
                       WORD gap_open_penalty,
                       WORD gap_extend_penalty,
                       ScoreMatrices const & matrices,
                       WORD * dprofile,
                       WORD * hearray,
                       db_thread_s & dbt,
                       long sequences,
                       long const * seqnos,
                       long * scores,
                       long * bestpos,
                       int qlen,
                       std::size_t kernel_bytes) -> void;

// 16-bit search, 1 database residue per block, with the alignment end
// on both sequences (the alignment threads)
auto search16s_dispatch(WORD * * q_start,
                        WORD gap_open_penalty,
                        WORD gap_extend_penalty,
                        ScoreMatrices const & matrices,
                        WORD * dprofile,
                        WORD * hearray,
                        DbThread const * dbta,
                        long sequences,
                        long const * seqnos,
                        long * scores,
                        long * bestpos,
                        long * bestq,
                        int qlen,
                        std::size_t kernel_bytes) -> void;

#endif  // SWIPE_ARCH_SEARCH_DISPATCH_H
