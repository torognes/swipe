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

#include "search_data.h"
#include <algorithm>  // std::generate
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <functional>  // std::cref
#include <iterator>  // std::begin, std::end, std::next
#include <mutex>  // std::mutex, std::lock_guard
#include <thread>
#include <vector>

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

// the alignment work is distributed in 7 bins: one per query strand
// and frame (3 x 2), and one for the hits that are not aligned
constexpr std::size_t unaligned_bin = frame_count;
constexpr std::size_t align_bins = frame_count + 1;

// the distribution of the hits among the alignment threads: chunks of
// each bin, handed out by align_getwork() under the mutex
struct AlignWork
{
  std::mutex mutex;
  long maxchunksize = 0;  // the largest chunk: the size of the lists
  long alignedhits = 0;  // the first hit of the next chunk
  Buffer<long> hits_sorted;  // the hits, in alignment order
  std::size_t volnext = 0;  // the next bin with chunks left
  std::array<long, align_bins> volseqs {{}};  // the hits left in each bin
  std::array<long, align_bins> volchunks {{}};  // the chunks left in each bin
};

AlignWork align_work;

// the bytes of a symbol's row in the score profile of search16s(): 1
// database residue x 16 bytes of lanes
constexpr std::ptrdiff_t profile_row_bytes = 16;

auto align_init(Parameters const & parameters, search_data & data) -> void
{
  data.dbt = db_thread_create();

  std::generate(std::begin(data.dbta), std::end(data.dbta), db_thread_create);

  data.dprofile.resize(profile_bytes);
  long const hearraylen = query_tables_init(parameters, data, profile_row_bytes);

  //  fprintf(out, "hearray length = %ld\n", hearraylen);

  data.hearray.resize(static_cast<std::size_t>(hearraylen) * hearray_row_bytes);

  auto const listsize = static_cast<std::size_t>(align_work.maxchunksize);
  //  if ((symtype == 3) || (symtype == 4))
  //    listsize *= 6;

  data.start_list.resize(listsize);
  data.start_hits.resize(listsize);
  data.in_list.resize(listsize);
  data.out_list.resize(listsize);
  data.scores.resize(listsize);
  data.bestpos.resize(listsize);
  data.bestq.resize(listsize);

  if (parameters.symtype == SymbolType::blastn)
  {
    data.qstrand1 = parameters.querystrands == QueryStrands::minus ? 1 : 0;
    data.qframe1 = 0;
    data.qstrand2 = parameters.querystrands == QueryStrands::plus ? 0 : 1;
    data.qframe2 = 0;

    data.dstrand1 = 0;
    data.dframe1 = 0;
    data.dstrand2 = 0;
    data.dframe2 = 0;
  }
  else if (parameters.symtype == SymbolType::blastx)
  {
    data.qstrand1 = parameters.querystrands == QueryStrands::minus ? 1 : 0;
    data.qframe1 = 0;
    data.qstrand2 = parameters.querystrands == QueryStrands::plus ? 0 : 1;
    data.qframe2 = 2;

    data.dstrand1 = 0;
    data.dframe1 = 0;
    data.dstrand2 = 0;
    data.dframe2 = 0;
  }
  else if (parameters.symtype == SymbolType::tblastn)
  {
    data.qstrand1 = 0;
    data.qframe1 = 0;
    data.qstrand2 = 0;
    data.qframe2 = 0;

    data.dstrand1 = 0;
    data.dframe1 = 0;
    data.dstrand2 = 1;
    data.dframe2 = 2;
  }
  else if (parameters.symtype == SymbolType::tblastx)
  {
    data.qstrand1 = parameters.querystrands == QueryStrands::minus ? 1 : 0;
    data.qframe1 = 0;
    data.qstrand2 = parameters.querystrands == QueryStrands::plus ? 0 : 1;
    data.qframe2 = 2;

    data.dstrand1 = 0;
    data.dframe1 = 0;
    data.dstrand2 = 1;
    data.dframe2 = 2;
  }
  else
  {
    data.qstrand1 = 0;
    data.qframe1 = 0;
    data.qstrand2 = 0;
    data.qframe2 = 0;

    data.dstrand1 = 0;
    data.dframe1 = 0;
    data.dstrand2 = 0;
    data.dframe2 = 0;
  }
}

// the hits hitfirst to hitlast of the list, for align_chunk()
struct HitChunk
{
  long first;
  long last;
};

auto align_chunk(Parameters const & parameters, search_data & data, HitChunk const chunk) -> void
{
  long const hitfirst = chunk.first;
  long const hitlast = chunk.last;
  auto const chunk_hits = make_view(align_work.hits_sorted)
    .subspan(static_cast<std::size_t>(hitfirst), static_cast<std::size_t>(hitlast - hitfirst + 1));
  if (hitlast < parameters.alignments)
  {

    for (long qstrand = data.qstrand1; qstrand <= data.qstrand2; qstrand++)
    {
      for(long qframe = data.qframe1; qframe <= data.qframe2; qframe++)
      {
	data.start_count = 0;

	for (auto const hs : chunk_hits)
	{
	  auto const hit = hits_gethit(hs);

	  if ((qstrand == hit.strands.qstrand) && (qframe == hit.strands.qframe))
	  {
	    data.start_hits[data.start_count] = hs;
	    data.start_list[data.start_count] = 
	      (hit.seqno << 3) | (hit.strands.dstrand << 2) | hit.strands.dframe;
	    data.start_count++;
	  }
	}

	if (data.start_count != 0)
	{
	  //	  printf("Aligning %ld sequences.\n", data.start_count);


	  BYTE ** qtable = data.qtable[frame_index(qstrand, qframe)].data();
	  long const qlen = data.qlen[frame_index(qstrand, qframe)];
      
	  /* 16-bit search, 8x1 db symbols, with alignment end */
	  
	  // the 16-bit penalties are only used when they fit (KI-13:
	  // otherwise no 16-bit result is accepted)
	  search16s(reinterpret_cast<WORD**>(qtable),
		    static_cast<WORD>(parameters.gapopenextend),
		    static_cast<WORD>(parameters.gapextend),
		    reinterpret_cast<WORD*>(score_matrices.score_16.data()),
		    reinterpret_cast<WORD*>(data.dprofile.data()),
		    reinterpret_cast<WORD*>(data.hearray.data()),
		    data.dbta.data(),
		    static_cast<long>(data.start_count),
		    data.start_list.data(),
		    data.scores.data(),
		    data.bestpos.data(),
		    data.bestq.data(),
		    static_cast<int>(qlen));
	
	  for (std::size_t i = 0; i < data.start_count; i++)
	  {
	    long const pos = data.bestpos[i];
	    long const bestq = data.bestq[i];
	  
	    //	  fprintf(out, "seqno=%ld score=%ld bestpos=%ld\n", seqno, score, pos);
	  
	    long const hitno = data.start_hits[i];

	    if (data.scores[i] < score_matrices.limit_16)
	    {
	      hits_enter_align_hint(hitno, bestq, pos);
	    }
	  }
	}
      }
    }
  }

  for (auto const hitno : chunk_hits)
  {
    hits_align(parameters, *data.dbt, hitno);
  }
}

auto align_done(search_data & data) -> void
{

  for (auto & db_thread : data.dbta)
  {
    db_thread.reset();
  }

  data.dbt.reset();
}



auto align_threads_init(Parameters const & parameters) -> void
{
  auto const hits = hits_getcount();

  align_work.hits_sorted = hits_sort();

  align_work.volseqs.fill(0);

  for(long i = 0; i<hits; i++)
  {
    if (i >= parameters.alignments)
    {
      align_work.volseqs[unaligned_bin]++;
    }
    else
    {
      auto const strands = hits_gethit(i).strands;
      align_work.volseqs[frame_index(strands.qstrand, strands.qframe)]++;
    }
  }

  align_work.maxchunksize = calc_chunks(make_view(align_work.volseqs),
                                        align_work.volchunks.data(),
                                        {parameters.threads, static_cast<long>(channels_16)});

  align_work.alignedhits = 0;
  align_work.volnext = 0;

  align_work.volnext = next_bin_with_chunks(make_view(align_work.volchunks), align_work.volnext);
}

auto align_threads_done() -> void
{
  align_work.hits_sorted = Buffer<long>();
}

auto align_getwork(long * first, long * last) -> int
{
  int status = 0;

  std::lock_guard<std::mutex> const lock(align_work.mutex);
  if (align_work.volnext < align_bins)
  {
    long const seqcount = align_work.volseqs[align_work.volnext];
    long const chunks = align_work.volchunks[align_work.volnext];
    long const chunksize = ((seqcount+chunks-1) / chunks);

    * first = align_work.alignedhits;
    * last = align_work.alignedhits + chunksize - 1;

    align_work.alignedhits += chunksize;
    status = 1;

    align_work.volseqs[align_work.volnext] -= chunksize;
    align_work.volchunks[align_work.volnext]--;

    align_work.volnext = next_bin_with_chunks(make_view(align_work.volchunks), align_work.volnext);
  }
  return status;
}

auto align_worker(Parameters const & parameters) -> void
{
  search_data sd;
  align_init(parameters, sd);

  long i = 0;
  long j = 0;
  while (align_getwork(&i, &j) != 0)
  {
    align_chunk(parameters, sd, {i, j});
  }

  align_done(sd);
}

}  // anonymous namespace

auto align_threads(Parameters const & parameters) -> void
{
  align_threads_init(parameters);

  std::vector<std::thread> workers;
  workers.reserve(static_cast<std::size_t>(parameters.threads));
  for (long t = 0; t < parameters.threads; t++)
  {
    workers.emplace_back(align_worker, std::cref(parameters));
  }
  for (auto & worker_thread : workers)
  {
    worker_thread.join();
  }

  align_threads_done();
}
