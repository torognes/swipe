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

// The SIMD search kernels (search7.cc, search16.cc, search16s.cc, and
// their x86-64 variants: search7_ssse3, compiled from search7.cc with
// -mssse3, and the AVX2 kernels of search7_avx2.cc and
// search16_avx2.cc). Only the kernel dispatch (arch/*/search_dispatch.cc)
// calls them: the other files call the dispatch functions of
// arch/search_dispatch.h.
// [moved from swipe.h]

#ifndef SWIPE_SEARCH_KERNELS_H
#define SWIPE_SEARCH_KERNELS_H

#include "swipe.h"  // BYTE, WORD, db_thread_s, DbThread


auto search7(BYTE * * q_start,
	     BYTE gap_open_penalty,
	     BYTE gap_extend_penalty,
	     BYTE const * score_matrix,
	     BYTE * dprofile,
	     BYTE * hearray,
	     db_thread_s & dbt,
	     long sequences,
	     long const * seqnos,
	     long * scores,
	     long qlen) -> void;

auto search7_ssse3(BYTE * * q_start,
		   BYTE gap_open_penalty,
		   BYTE gap_extend_penalty,
		   BYTE const * score_matrix,
		   BYTE * dprofile,
		   BYTE * hearray,
		   db_thread_s & dbt,
		   long sequences,
		   long const * seqnos,
		   long * scores,
		   long qlen) -> void;

auto search7_avx2(BYTE * * q_start,
		  BYTE gap_open_penalty,
		  BYTE gap_extend_penalty,
		  BYTE const * score_matrix,
		  BYTE * dprofile,
		  BYTE * hearray,
		  db_thread_s & dbt,
		  long sequences,
		  long const * seqnos,
		  long * scores,
		  long qlen) -> void;

auto search16(WORD * * q_start,
	      WORD gap_open_penalty,
	      WORD gap_extend_penalty,
	      WORD const * score_matrix,
	      WORD * dprofile,
	      WORD * hearray,
	      db_thread_s & dbt,
	      long sequences,
	      long const * seqnos,
	      long * scores,
	      long * bestpos,
	      int qlen) -> void;

auto search16s(WORD * * q_start,
	       WORD gap_open_penalty,
	       WORD gap_extend_penalty,
	       WORD const * score_matrix,
	       WORD * dprofile,
	       WORD * hearray,
	       DbThread const * dbta,
	       long sequences,
	       long const * seqnos,
	       long * scores,
	       long * bestpos,
	       long * bestq,
	       int qlen) -> void;

auto search16_avx2(WORD * * q_start,
		   WORD gap_open_penalty,
		   WORD gap_extend_penalty,
		   WORD const * score_matrix,
		   WORD * dprofile,
		   WORD * hearray,
		   db_thread_s & dbt,
		   long sequences,
		   long const * seqnos,
		   long * scores,
		   long * bestpos,
		   int qlen) -> void;

auto search16s_avx2(WORD * * q_start,
		    WORD gap_open_penalty,
		    WORD gap_extend_penalty,
		    WORD const * score_matrix,
		    WORD * dprofile,
		    WORD * hearray,
		    DbThread const * dbta,
		    long sequences,
		    long const * seqnos,
		    long * scores,
		    long * bestpos,
		    long * bestq,
		    int qlen) -> void;

#endif  // SWIPE_SEARCH_KERNELS_H
