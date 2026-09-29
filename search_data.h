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

// the data of a search or alignment thread, and the functions shared
// by the search threads (search_threads.cc) and the alignment threads
// (align_threads.cc)

#ifndef SWIPE_SEARCH_DATA_H
#define SWIPE_SEARCH_DATA_H

#include "swipe.h"
#include "fatal_allocator.h"  // Buffer
#include <cstddef>  // std::size_t

// the score profile of the SIMD kernels: 32 symbols x 64 bytes (4
// database residues x 16 bytes of lanes)
constexpr std::size_t profile_bytes = std::size_t{32} * 64;

struct search_data
{
  struct db_thread_s * dbt;
  struct db_thread_s * dbta[8];

  Buffer<BYTE> dprofile;  // profile_bytes
  Buffer<BYTE> hearray;
  Buffer<BYTE *> qtable[6];  // empty: tables not allocated

  Buffer<long> scores;
  Buffer<long> bestpos;
  Buffer<long> bestq;
  Buffer<long> start_list;
  Buffer<long> in_list;
  Buffer<long> out_list;

  long qlen[6];

  std::size_t start_count;
  std::size_t in_count;
  std::size_t out_count;

  Buffer<long> start_hits;

  long seqfirst, seqlast;

  long qstrand1, qstrand2, qframe1, qframe2;
  long dstrand1, dstrand2, dframe1, dframe2;
};

auto calc_chunks(long volcount,
		 long par,
		 long channels,
		 long const * volume_sequences,
		 long * volume_chunks,
		 long * totalchunks,
		 long * biggestchunk) -> void;

// search_threads.cc: the search of a query by parameters.threads threads
auto prepare_search(long par) -> void;
auto run_threads(Parameters const & parameters) -> void;

// align_threads.cc: the alignment of the hits of a query
auto align_threads(Parameters const & parameters) -> void;

#endif  // SWIPE_SEARCH_DATA_H
