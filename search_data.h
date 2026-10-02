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
#include <algorithm>  // std::find_if
#include <array>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator>  // std::distance, std::next

// the score profile of the SIMD kernels: 32 symbols x 64 bytes (4
// database residues x 16 bytes of lanes)
constexpr std::size_t profile_bytes = std::size_t{32} * 64;
// the profile of the AVX2 kernel: 32 symbols x 128 bytes (4 database
// residues x 32 bytes of lanes)
constexpr std::size_t profile_bytes_avx2 = std::size_t{32} * 128;

struct search_data
{
  DbThread dbt;
  std::array<DbThread, 8> dbta;

  Buffer<BYTE> dprofile;  // profile_bytes
  Buffer<BYTE> hearray;
  std::array<Buffer<BYTE *>, frame_count> qtable;  // empty: tables not allocated
  std::array<Buffer<BYTE *>, frame_count> qtable_avx2;  // rows of profile_bytes_avx2

  Buffer<long> scores;
  Buffer<long> bestpos;
  Buffer<long> bestq;
  Buffer<long> start_list;
  Buffer<long> in_list;
  Buffer<long> out_list;

  std::array<long, frame_count> qlen;

  std::size_t start_count;
  std::size_t in_count;
  std::size_t out_count;

  Buffer<long> start_hits;

  long seqfirst, seqlast;

  long qstrand1, qstrand2, qframe1, qframe2;
  long dstrand1, dstrand2, dframe1, dframe2;
};

// the first bin (volume, or bin of hits) from 'first' on with chunks
// left, or chunks.size() when none is left; first is at most
// chunks.size()
inline auto next_bin_with_chunks(View<long> const chunks, std::size_t const first) -> std::size_t
{
  assert(first <= chunks.size());
  auto const found = std::find_if(std::next(chunks.begin(), static_cast<std::ptrdiff_t>(first)),
                                  chunks.end(),
                                  [](long const count) -> bool { return count != 0; });
  return static_cast<std::size_t>(std::distance(chunks.begin(), found));
}

// the threads that share the work, and the channels of their kernel
struct Chunking
{
  long threads;
  long channels;
};

// the chunks of each volume (or bin) of volume_sequences, written to
// volume_chunks (as many entries); returns the size of the largest one
auto calc_chunks(View<long> volume_sequences,
		 long * volume_chunks,
		 Chunking chunking) -> long;

// the query tables (qtables: data.qtable or data.qtable_avx2) and
// lengths (data.qlen) of the strands or frames searched, for a profile
// with rows of row_bytes; returns the longest query length
auto query_tables_init(Parameters const & parameters,
		       search_data & data,
		       std::ptrdiff_t row_bytes,
		       std::array<Buffer<BYTE *>, frame_count> & qtables) -> long;

// search_threads.cc: the search of a query by parameters.threads threads
auto prepare_search(long par) -> void;
auto run_threads(Parameters const & parameters) -> void;

// align_threads.cc: the alignment of the hits of a query
auto align_threads(Parameters const & parameters) -> void;

#endif  // SWIPE_SEARCH_DATA_H
