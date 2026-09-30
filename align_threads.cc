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

auto align_init(Parameters const & parameters, struct search_data * sdp) -> void
{
  sdp->dbt = db_thread_create();

  std::generate(std::begin(sdp->dbta), std::end(sdp->dbta), db_thread_create);

  sdp->dprofile.resize(profile_bytes);
  long hearraylen = 0;

  if (parameters.symtype == SymbolType::blastn)
  {
    for (long s = 0; s < 2; s++)
    {
      if (searches_strand(parameters.querystrands, s))
      {
	long const qlen = query.nt[strand_index(s)].len;
	sdp->qlen[frame_index(s, 0)] = qlen;
	sdp->qtable[frame_index(s, 0)].resize(static_cast<std::size_t>(qlen));
	for (std::size_t i = 0; i < sdp->qtable[frame_index(s, 0)].size(); i++)
	{
	  sdp->qtable[frame_index(s, 0)][i] = std::next(sdp->dprofile.data(), profile_row_bytes * query.nt[strand_index(s)].seq[i]);
	}
	hearraylen = qlen > hearraylen ? qlen : hearraylen;
      }
    }
  }
  else if ((parameters.symtype == SymbolType::blastp) || (parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::sound))
  {
    long const qlen = query.aa[0].len;
    sdp->qlen[0] = qlen;
    sdp->qtable[0].resize(static_cast<std::size_t>(qlen));
    for (std::size_t i = 0; i < sdp->qtable[0].size(); i++)
    {
      sdp->qtable[0][i] = std::next(sdp->dprofile.data(), profile_row_bytes * query.aa[0].seq[i]);
    }
    hearraylen = qlen > hearraylen ? qlen : hearraylen;
  }
  else if ((parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastx))
  {
    for (long s = 0; s < 2; s++)
    {
      if (searches_strand(parameters.querystrands, s))
      {
	for(long f=0; f<3; f++)
	{
	  long const qlen = query.aa[frame_index(s, f)].len;
	  sdp->qlen[frame_index(s, f)] = qlen;
	  sdp->qtable[frame_index(s, f)].resize(static_cast<std::size_t>(qlen));
	  for (std::size_t i = 0; i < sdp->qtable[frame_index(s, f)].size(); i++)
	  {
	    sdp->qtable[frame_index(s, f)][i] = std::next(sdp->dprofile.data(), profile_row_bytes * query.aa[frame_index(s, f)].seq[i]);
	  }
	  hearraylen = qlen > hearraylen ? qlen : hearraylen;
	}
      }
    }
  }
  
  //  fprintf(out, "hearray length = %ld\n", hearraylen);

  sdp->hearray.resize(static_cast<std::size_t>(hearraylen) * hearray_row_bytes);

  auto const listsize = static_cast<std::size_t>(align_work.maxchunksize);
  //  if ((symtype == 3) || (symtype == 4))
  //    listsize *= 6;

  sdp->start_list.resize(listsize);
  sdp->start_hits.resize(listsize);
  sdp->in_list.resize(listsize);
  sdp->out_list.resize(listsize);
  sdp->scores.resize(listsize);
  sdp->bestpos.resize(listsize);
  sdp->bestq.resize(listsize);

  if (parameters.symtype == SymbolType::blastn)
  {
    sdp->qstrand1 = parameters.querystrands == QueryStrands::minus ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = parameters.querystrands == QueryStrands::plus ? 0 : 1;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (parameters.symtype == SymbolType::blastx)
  {
    sdp->qstrand1 = parameters.querystrands == QueryStrands::minus ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = parameters.querystrands == QueryStrands::plus ? 0 : 1;
    sdp->qframe2 = 2;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (parameters.symtype == SymbolType::tblastn)
  {
    sdp->qstrand1 = 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = 0;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 1;
    sdp->dframe2 = 2;
  }
  else if (parameters.symtype == SymbolType::tblastx)
  {
    sdp->qstrand1 = parameters.querystrands == QueryStrands::minus ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = parameters.querystrands == QueryStrands::plus ? 0 : 1;
    sdp->qframe2 = 2;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 1;
    sdp->dframe2 = 2;
  }
  else
  {
    sdp->qstrand1 = 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = 0;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
}

auto align_chunk(Parameters const & parameters, struct search_data * sdp, long hitfirst, long hitlast) -> void
{
  if (hitlast < parameters.alignments)
  {

    for (long qstrand = sdp->qstrand1; qstrand <= sdp->qstrand2; qstrand++)
    {
      for(long qframe = sdp->qframe1; qframe <= sdp->qframe2; qframe++)
      {
	sdp->start_count = 0;

	for(long hitno = hitfirst; hitno <= hitlast; hitno++)
	{
	  long const hs = align_work.hits_sorted[static_cast<std::size_t>(hitno)];
	  auto const hit = hits_gethit(hs);

	  if ((qstrand == hit.strands.qstrand) && (qframe == hit.strands.qframe))
	  {
	    sdp->start_hits[sdp->start_count] = hs;
	    sdp->start_list[sdp->start_count] = 
	      (hit.seqno << 3) | (hit.strands.dstrand << 2) | hit.strands.dframe;
	    sdp->start_count++;
	  }
	}

	if (sdp->start_count != 0)
	{
	  //	  printf("Aligning %ld sequences.\n", sdp->start_count);


	  BYTE ** qtable = sdp->qtable[frame_index(qstrand, qframe)].data();
	  long const qlen = sdp->qlen[frame_index(qstrand, qframe)];
      
	  /* 16-bit search, 8x1 db symbols, with alignment end */
	  
	  // the 16-bit penalties are only used when they fit (KI-13:
	  // otherwise no 16-bit result is accepted)
	  search16s(reinterpret_cast<WORD**>(qtable),
		    static_cast<WORD>(parameters.gapopenextend),
		    static_cast<WORD>(parameters.gapextend),
		    reinterpret_cast<WORD*>(score_matrices.score_16.data()),
		    reinterpret_cast<WORD*>(sdp->dprofile.data()),
		    reinterpret_cast<WORD*>(sdp->hearray.data()),
		    sdp->dbta.data(),
		    static_cast<long>(sdp->start_count),
		    sdp->start_list.data(),
		    sdp->scores.data(),
		    sdp->bestpos.data(),
		    sdp->bestq.data(),
		    static_cast<int>(qlen));
	
	  for (std::size_t i = 0; i < sdp->start_count; i++)
	  {
	    long const pos = sdp->bestpos[i];
	    long const bestq = sdp->bestq[i];
	  
	    //	  fprintf(out, "seqno=%ld score=%ld bestpos=%ld\n", seqno, score, pos);
	  
	    long const hitno = sdp->start_hits[i];

	    if (sdp->scores[i] < score_matrices.limit_16)
	    {
	      hits_enter_align_hint(hitno, bestq, pos);
	    }
	  }
	}
      }
    }
  }

  for (long hitno = hitfirst; hitno <= hitlast; hitno++)
  {
    hits_align(parameters, sdp->dbt, align_work.hits_sorted[static_cast<std::size_t>(hitno)]);
  }
}

auto align_done(struct search_data * sdp) -> void
{

  for (auto * db_thread : sdp->dbta)
  {
    db_thread_destruct(db_thread);
  }

  db_thread_destruct(sdp->dbt);
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

  long totalchunks = 0;

  calc_chunks(static_cast<long>(align_bins),
	      parameters.threads,
	      static_cast<long>(channels_16),
	      align_work.volseqs.data(),
	      align_work.volchunks.data(),
	      & totalchunks,
	      & align_work.maxchunksize);

  align_work.alignedhits = 0;
  align_work.volnext = 0;

  while ((align_work.volnext < align_bins) && (align_work.volchunks[align_work.volnext] == 0))
  {
    align_work.volnext++;
  }
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

    while ((align_work.volnext < align_bins) && (align_work.volchunks[align_work.volnext] == 0))
    {
      align_work.volnext++;
    }
  }
  return status;
}

auto align_worker(Parameters const & parameters) -> void
{
  search_data sd;
  align_init(parameters, &sd);

  long i = 0;
  long j = 0;
  while (align_getwork(&i, &j) != 0)
  {
    align_chunk(parameters, &sd, i, j);
  }

  align_done(&sd);
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
