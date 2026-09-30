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
#include <algorithm>  // std::copy_n, std::max, std::min
#include <cmath>  // std::floor, std::sqrt
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <functional>  // std::cref
#include <iterator>  // std::next
#include <limits>
#include <mutex>  // std::mutex, std::lock_guard
#include <thread>
#include <utility>  // std::swap
#include <vector>

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

// the distribution of the database sequences among the search threads:
// chunks of each volume, handed out by search_getwork() under the mutex
struct SearchWork
{
  std::mutex mutex;
  std::mutex count_mutex;  // for run.compute7 (swipe.cc)
  long maxchunksize = 0;  // the largest chunk: the size of the lists
  std::size_t volnext = 0;  // the next volume with chunks left
  long seqnext = 0;  // the first sequence of the next chunk
  Buffer<long> volchunks;  // the chunks left in each volume
  Buffer<long> volseqs;  // the sequences left in each volume
};

SearchWork search_work;

// the bytes of a symbol's row in the score profile of search7() and
// search16(): 4 database residues x 16 bytes of lanes
constexpr std::ptrdiff_t profile_row_bytes = 64;

auto search_init(Parameters const & parameters, struct search_data * sdp) -> void
{
  sdp->dbt = db_thread_create();
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

  // at least one row: the kernels memset() the array, and an empty
  // Buffer has no storage (a null data(), for an empty query)
  sdp->hearray.resize(static_cast<std::size_t>(std::max(hearraylen, 1L)) * hearray_row_bytes);

  auto listsize = static_cast<std::size_t>(search_work.maxchunksize);
  if ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx))
  {
    listsize *= frame_count;  // the database frames
  }

  sdp->start_list.resize(listsize);
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

auto search_done(struct search_data * sdp) -> void
{
  db_thread_destruct(sdp->dbt);
}

auto search_getwork(long * first, long * last) -> int
{
  int status = 0;
  auto const volcount = static_cast<std::size_t>(db_getvolumecount());
  
  std::lock_guard<std::mutex> const lock(search_work.mutex);
  if (search_work.volnext < volcount)
  {
    long const seqcount = search_work.volseqs[search_work.volnext];
    long const chunks = search_work.volchunks[search_work.volnext];
    long const chunksize = ((seqcount+chunks-1) / chunks);

    * first = search_work.seqnext;
    * last = search_work.seqnext + chunksize - 1;
    search_work.seqnext += chunksize;
    status = 1;

    //    fprintf(out, "Processing sequences %d to %d (%d sequences) in volume %ld.\n", *first, *last, *last - * first + 1, volnext);

    search_work.volseqs[search_work.volnext] -= chunksize;
    search_work.volchunks[search_work.volnext]--;

    while ((search_work.volnext < volcount) && (search_work.volchunks[search_work.volnext] == 0))
    {
      search_work.volnext++;
    }
  }
  return status;
}


// blastn: a hit of the reverse complement of the query is entered as
// a hit of the query on the reverse strand of the database sequence
auto reported_strands(SymbolType const symbol_type, HitStrands const & strands) -> HitStrands
{
  if ((symbol_type == SymbolType::blastn) && (strands.qstrand != 0))
  {
    return {0, 0, 1, 0};
  }
  return strands;
}

auto search_chunk(Parameters const & parameters, struct search_data * sdp) -> void
{
  // the 7-bit engine uses signed bytes: gap penalties are clamped to
  // 127 (KI-11). This is exact: 7-bit scores are in [0, 127], so a
  // penalty of 127 already takes any score down to zero. The 16-bit
  // and 63-bit engines, and the alignments, use the real penalties
  long const max_7 = std::numeric_limits<signed char>::max();
  auto const gapopenextend_7 = static_cast<BYTE>(std::min(parameters.gapopenextend, max_7));
  auto const gapextend_7 = static_cast<BYTE>(std::min(parameters.gapextend, max_7));

  //  fprintf(out, "Searching seqnos %ld to %ld\n", sdp->seqfirst, sdp->seqlast);

  if (parameters.taxidfilename != nullptr)
  {
    db_mapheaders(sdp->dbt, sdp->seqfirst, sdp->seqlast);
  }

  sdp->start_count = 0;
  for(long seqno = sdp->seqfirst; seqno <= sdp->seqlast; seqno++)
  {
    if (db_check_inclusion(sdp->dbt, seqno) != 0)
    {
      if ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx))
      {
	for (long dstrand = sdp->dstrand1; dstrand <= sdp->dstrand2; dstrand++)
	{
	  for(long dframe = sdp->dframe1; dframe <= sdp->dframe2; dframe++)
	  {
	    sdp->start_list[sdp->start_count++] =
	      (seqno << 3) | (dstrand << 2) | dframe;
	  }
	}
      }
      else
      {
	sdp->start_list[sdp->start_count++] = seqno << 3;
      }
    }
  }

  if (sdp->start_count == 0)
  {
    return;
  }

  long const s1 = sdp->start_list[0] >> 3;
  long const s2 = sdp->start_list[sdp->start_count-1] >> 3;
  
  // fprintf(out, "Mapping seqnos %ld to %ld\n", s1, s2);

  db_mapsequences(sdp->dbt, s1, s2);

  for (long qstrand = sdp->qstrand1; qstrand <= sdp->qstrand2; qstrand++)
  {
    for(long qframe = sdp->qframe1; qframe <= sdp->qframe2; qframe++)
    {
      sdp->out_count = sdp->start_count;
      std::copy_n(sdp->start_list.begin(), sdp->start_count, sdp->out_list.begin());
      
      BYTE ** qtable = sdp->qtable[frame_index(qstrand, qframe)].data();
      long const qlen = sdp->qlen[frame_index(qstrand, qframe)];
      
      /* 7-bit search */
	  
      std::swap(sdp->in_list, sdp->out_list);
      sdp->in_count = sdp->out_count;
	  
      if (sdp->in_count > 0)
      {
	{
	  std::lock_guard<std::mutex> const lock(search_work.count_mutex);
	  run.compute7 += static_cast<long>(sdp->in_count);
	}
	    
	// fprintf(out, "Searching seqnos %ld to %ld\n", sdp->in_list[0], sdp->in_list[sdp->in_count-1]);

	if (cpu_features.ssse3)
	{
	  search7_ssse3(qtable,
			gapopenextend_7,
			gapextend_7,
			reinterpret_cast<BYTE const *>(score_matrices.score_7t.data()),
			sdp->dprofile.data(),
			sdp->hearray.data(),
			sdp->dbt,
			static_cast<long>(sdp->in_count),
			sdp->in_list.data(),
			sdp->scores.data(),
			qlen);
	}
	else
	{
	  search7(qtable,
		  gapopenextend_7,
		  gapextend_7,
		  reinterpret_cast<BYTE const *>(score_matrices.score_7.data()),
		  sdp->dprofile.data(),
		  sdp->hearray.data(),
		  sdp->dbt,
		  static_cast<long>(sdp->in_count),
		  sdp->in_list.data(),
		  sdp->scores.data(),
		  qlen);
	}

	sdp->out_count = 0;
    
	for (std::size_t i = 0; i < sdp->in_count; i++)
	{
	  long const seqnosf = sdp->in_list[i];
	  long const score = sdp->scores[i];
      
	  if (score < score_matrices.limit_7)
	  {
	    long const seqno = seqnosf >> 3;
	    long const dstrand = (seqnosf >> 2) & 1;
	    long const dframe = seqnosf & 3;

	    hits_enter(seqno, score,
		       reported_strands(parameters.symtype, {qstrand, qframe, dstrand, dframe}));
	  }
	  else
	  {
	    sdp->out_list[sdp->out_count++] = seqnosf;
	  }
	}
      }

      /* 16-bit search */
	  
      std::swap(sdp->in_list, sdp->out_list);
      sdp->in_count = sdp->out_count;
  
      if (sdp->in_count > 0)
      {
	  
	// the 16-bit penalties are only used when they fit (KI-13:
	// otherwise no 16-bit result is accepted)
	search16(reinterpret_cast<WORD**>(qtable),
		 static_cast<WORD>(parameters.gapopenextend),
		 static_cast<WORD>(parameters.gapextend),
		 reinterpret_cast<WORD*>(score_matrices.score_16.data()),
		 reinterpret_cast<WORD*>(sdp->dprofile.data()),
		 reinterpret_cast<WORD*>(sdp->hearray.data()),
		 sdp->dbt,
		 static_cast<long>(sdp->in_count),
		 sdp->in_list.data(),
		 sdp->scores.data(),
		 sdp->bestpos.data(),
		 static_cast<int>(qlen));
    
	sdp->out_count = 0;
    
	for (std::size_t i = 0; i < sdp->in_count; i++)
	{
	  long const seqnosf = sdp->in_list[i];
	  long const score = sdp->scores[i];
	  if (score < score_matrices.limit_16)
	  {
	    long const seqno = seqnosf >> 3;
	    long const dstrand = (seqnosf >> 2) & 1;
	    long const dframe = seqnosf & 3;

	    hits_enter(seqno, score,
		       reported_strands(parameters.symtype, {qstrand, qframe, dstrand, dframe}));
	  }
	  else
	  {
	    sdp->out_list[sdp->out_count++] = seqnosf;
	  }
	}
      }
      
      /* 63-bit search */

      std::swap(sdp->in_list, sdp->out_list);
      sdp->in_count = sdp->out_count;
  
      if (sdp->in_count > 0)
      {
    
	for (std::size_t i = 0; i < sdp->in_count; i++)
	{
	  long const seqnosf = sdp->in_list[i];
	  long const seqno = seqnosf >> 3;
	  long const dstrand = (seqnosf >> 2) & 1;
	  long const dframe = seqnosf & 3;
      
	  long ntlen = 0;
	  View<char> const sequence = db_getsequence(sdp->dbt, seqno, {dstrand, dframe},
						     & ntlen, 0);
	  auto const * dbegin = sequence.begin();
	  auto const * dend = sequence.end();
      
	  char const * q = nullptr;
	  if (parameters.symtype == SymbolType::blastn)
	  {
	    q = query.nt[strand_index(qstrand)].seq;
	  }
	  else
	  {
	    q = query.aa[frame_index(qstrand, qframe)].seq;
	  }

	  long const score = fullsw(dbegin,
			      dend,
			      q, 
			      std::next(q, qlen),
			      reinterpret_cast<long*>(sdp->hearray.data()),
			      score_matrices.score_63.data(),
			      parameters.gapopenextend,
			      parameters.gapextend);

	  hits_enter(seqno, score,
		     reported_strands(parameters.symtype, {qstrand, qframe, dstrand, dframe}));
	}
      }
  
    }
  }
}


auto worker(Parameters const & parameters) -> void
{
  struct search_data sd;
  search_init(parameters, &sd);

  while (search_getwork(&sd.seqfirst, &sd.seqlast) != 0)
  {
    search_chunk(parameters, &sd);
  }

  search_done(&sd);
}

}  // anonymous namespace

auto calc_chunks(View<long> const volume_sequences,
		 long * volume_chunks,
		 Chunking const chunking) -> long
{
  long const par = chunking.threads;
  long const channels = chunking.channels;

  long volsused = 0;
  auto const volumes = volume_sequences.size();
  std::vector<long> chunksizes(volumes);
  long totalseqs = 0;
  long biggest_chunk_size = 0;
  std::size_t vv = 0;
  for(std::size_t v = 0; v < volumes; v++)
  {
    if (volume_sequences[v] != 0)
    {
      totalseqs += volume_sequences[v];
      volsused++;
      volume_chunks[v] = 1;
      chunksizes[v] = volume_sequences[v];
      if (chunksizes[v] > biggest_chunk_size)
      {
	biggest_chunk_size = chunksizes[v];
	vv = v;
      }
    }
    else
    {
      chunksizes[v] = 0;
      volume_chunks[v] = 0;
    }
  }

  long upper = channels;
  if (totalseqs >= 4 * channels * par)
  {
    upper *= static_cast<long>(floor(sqrt((1.0 * static_cast<double>(totalseqs)) / static_cast<double>(channels * par))));
  }

  long chunks = volsused;
  long const minchunks = totalseqs < par ? totalseqs : par;

  while((biggest_chunk_size > upper) || (chunks < minchunks))
  {
    volume_chunks[vv]++;
    chunks++;
    chunksizes[vv] = (volume_sequences[vv] + volume_chunks[vv] - 1) / volume_chunks[vv];

    biggest_chunk_size = 0;
    vv = 0;
    for(std::size_t v = 0; v < volumes; v++)
    {
      if (chunksizes[v] > biggest_chunk_size)
      {
	vv = v;
	biggest_chunk_size = chunksizes[v];
      }
    }
  }
  
  return biggest_chunk_size;
}

auto prepare_search(long par) -> void
{
  search_work.volnext = 0;
  search_work.seqnext = 0;

  auto const volcount = static_cast<std::size_t>(db_getvolumecount());
  search_work.volseqs.resize(volcount);
  search_work.volchunks.resize(volcount);
  for (std::size_t v = 0; v < volcount; v++)
  {
    search_work.volseqs[v] = db_getseqcount_volume(static_cast<long>(v));
  }

  search_work.maxchunksize = calc_chunks(make_view(search_work.volseqs),
                                         search_work.volchunks.data(),
                                         {par, static_cast<long>(channels_7)});

  while ((search_work.volnext < volcount) && (search_work.volchunks[search_work.volnext] == 0))
  {
    search_work.volnext++;
  }
}

auto run_threads(Parameters const & parameters) -> void
{
  std::vector<std::thread> workers;
  workers.reserve(static_cast<std::size_t>(parameters.threads));
  for (long t = 0; t < parameters.threads; t++)
  {
    workers.emplace_back(worker, std::cref(parameters));
  }
  for (auto & worker_thread : workers)
  {
    worker_thread.join();
  }
}
