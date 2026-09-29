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

#include "swipe.h"
#include "fatal_allocator.h"  // Buffer
#include <algorithm>  // std::copy_n, std::generate, std::min
#include <array>
#include <cassert>
#include <cmath>  // std::floor, std::sqrt
#include <cstddef>  // std::size_t
#include <cstdlib>  // exit, posix_memalign
#include <functional>  // std::cref
#include <iterator>  // std::begin, std::end
#include <limits>
#include <mutex>  // std::mutex, std::lock_guard
#include <string>  // std::string (fatal)
#include <thread>
#include <utility>  // std::swap
#include <vector>


// the version number is read from the file VERSION by the Makefile
#ifndef SWIPE_VERSION
#ifdef __CPPCHECK__
// static analysis with cppcheck, run without the Makefile's flags
#define SWIPE_VERSION "0.0.0"
#else
#error "SWIPE_VERSION is not defined: build swipe with make"
#endif
#endif

extern char const swipe_name_and_version[] = "SWIPE " SWIPE_VERSION;

/* Other variables */

long queryno;

long cpu_feature_ssse3;
long cpu_feature_sse41;

long compute7;

long totalhits;

FILE * out = stdout;  // default output: stdout (--out FILE)

struct time_info ti;

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

long cpu_feature_sse2;
std::mutex countmutex;
std::mutex workmutex;
long maxchunksize;
std::size_t volnext;
long seqnext;
Buffer<long> volchunks;
Buffer<long> volseqs;

struct search_data
{
  struct db_thread_s * dbt;
  struct db_thread_s * dbta[8];

  Buffer<BYTE> dprofile;
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

}  // anonymous namespace

[[noreturn]] auto fatal(char const * message) noexcept -> void
{
  if (message != nullptr)
  {
    fprintf(stderr, "%s\n", message);
  }
  exit(1);
}

[[noreturn]] auto fatal(std::string const & message) noexcept -> void
{
  fatal(message.c_str());
}

auto xmalloc(size_t size) -> void *
{
  size_t const alignment = 16;
  void * t = nullptr;
  if (posix_memalign(&t, alignment, size) != 0)
  {
    t = nullptr;
  }

  if (t == nullptr)
  {
    fatal("Unable to allocate enough memory.");
  }

  return t;
}

namespace {

long alignedhits;
Buffer<long> hits_sorted;

std::size_t align_volnext;

// the alignment work is distributed in 7 bins: one per query strand
// and frame (3 x 2), and one for the hits that are not aligned
constexpr std::size_t align_bins = 7;
std::array<long, align_bins> align_volseqs {};
std::array<long, align_bins> align_volchunks {};

auto align_init(Parameters const & parameters, struct search_data * sdp) -> void
{
  sdp->dbt = db_thread_create();

  std::generate(std::begin(sdp->dbta), std::end(sdp->dbta), db_thread_create);

  sdp->dprofile.resize(4 * 16 * 32);
  long qlen = 0;
  long hearraylen = 0;

  if (parameters.symtype == SymbolType::blastn)
  {
    for (int s = 0; s < 2; s++)
    {
      if (searches_strand(parameters.querystrands, s))
      {
	qlen = query.nt[s].len;
	sdp->qlen[3*s] = qlen;
	sdp->qtable[3*s].resize(static_cast<std::size_t>(qlen));
	for (std::size_t i = 0; i < sdp->qtable[3*s].size(); i++)
	{
	  sdp->qtable[3*s][i] = std::next(sdp->dprofile.data(), 16 * query.nt[s].seq[i]);
	}
	hearraylen = qlen > hearraylen ? qlen : hearraylen;
      }
    }
  }
  else if ((parameters.symtype == SymbolType::blastp) || (parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::sound))
  {
    qlen = query.aa[0].len;
    sdp->qlen[0] = qlen;
    sdp->qtable[0].resize(static_cast<std::size_t>(qlen));
    for (std::size_t i = 0; i < sdp->qtable[0].size(); i++)
    {
      sdp->qtable[0][i] = std::next(sdp->dprofile.data(), 16 * query.aa[0].seq[i]);
    }
    hearraylen = qlen > hearraylen ? qlen : hearraylen;
  }
  else if ((parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastx))
  {
    for (int s = 0; s < 2; s++)
    {
      if (searches_strand(parameters.querystrands, s))
      {
	for(int f=0; f<3; f++)
	{
	  qlen = query.aa[(3*s)+f].len;
	  sdp->qlen[(3*s)+f] = qlen;
	  sdp->qtable[(3*s)+f].resize(static_cast<std::size_t>(qlen));
	  for (std::size_t i = 0; i < sdp->qtable[(3*s)+f].size(); i++)
	  {
	    sdp->qtable[(3*s)+f][i] = std::next(sdp->dprofile.data(), 16 * query.aa[(3*s)+f].seq[i]);
	  }
	  hearraylen = qlen > hearraylen ? qlen : hearraylen;
	}
      }
    }
  }
  
  //  fprintf(out, "hearray length = %ld\n", hearraylen);

  sdp->hearray.resize(static_cast<std::size_t>(hearraylen) * 32);

  auto const listsize = static_cast<std::size_t>(maxchunksize);
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
	  long const hs = hits_sorted[static_cast<std::size_t>(hitno)];
	  long seqno = 0;
	  long score = 0;
	  long hqstrand = 0;
	  long hqframe = 0;
	  long hdstrand = 0;
	  long hdframe = 0;
	
	  hits_gethit(hs, & seqno, & score, & hqstrand, & hqframe, 
		      & hdstrand, & hdframe);

	  if ((qstrand == hqstrand) && (qframe == hqframe))
	  {
	    sdp->start_hits[sdp->start_count] = hs;
	    sdp->start_list[sdp->start_count] = 
	      (seqno << 3) | (hdstrand << 2) | hdframe;
	    sdp->start_count++;
	  }
	}

	if (sdp->start_count != 0)
	{
	  //	  printf("Aligning %ld sequences.\n", sdp->start_count);


	  BYTE ** qtable = sdp->qtable[(3*qstrand)+qframe].data();
	  long const qlen = sdp->qlen[(3*qstrand)+qframe];
      
	  /* 16-bit search, 8x1 db symbols, with alignment end */
	  
	  // the 16-bit penalties are only used when they fit (KI-13:
	  // otherwise no 16-bit result is accepted)
	  search16s(reinterpret_cast<WORD**>(qtable),
		    static_cast<WORD>(parameters.gapopenextend),
		    static_cast<WORD>(parameters.gapextend),
		    reinterpret_cast<WORD*>(score_matrix_16),
		    reinterpret_cast<WORD*>(sdp->dprofile.data()),
		    reinterpret_cast<WORD*>(sdp->hearray.data()),
		    sdp->dbta,
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

	    if (sdp->scores[i] < SCORELIMIT_16)
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
    hits_align(parameters, sdp->dbt, hits_sorted[static_cast<std::size_t>(hitno)]);
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


auto calc_chunks(long volcount, 
		 long par,
		 long channels,
		 long const * volume_sequences,
		 long * volume_chunks,
		 long * totalchunks,
		 long * biggestchunk) -> void
{

  long volsused = 0;
  auto const volumes = static_cast<std::size_t>(volcount);
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
  
  *biggestchunk = biggest_chunk_size;
  *totalchunks = chunks;
}

auto align_threads_init(Parameters const & parameters) -> void
{
  long const hits = hits_getcount();

  hits_sorted = hits_sort();

  align_volseqs.fill(0);

  for(long i = 0; i<hits; i++)
  {
    long seqno = 0;
    long score = 0;
    long qstrand = 0;
    long qframe = 0;
    long dstrand = 0;
    long dframe = 0;

    if (i >= parameters.alignments)
    {
      align_volseqs[6]++;
    }
    else
    {
      hits_gethit(i, & seqno, & score,
		  & qstrand, & qframe,
		  & dstrand, & dframe);
      
      align_volseqs[static_cast<std::size_t>((3*qstrand)+qframe)]++;
    }
  }

  long totalchunks = 0;

  calc_chunks(static_cast<long>(align_bins),
	      parameters.threads,
	      8,
	      align_volseqs.data(),
	      align_volchunks.data(),
	      & totalchunks,
	      & maxchunksize);

  alignedhits = 0;
  align_volnext = 0;

  while ((align_volnext < align_bins) && (align_volchunks[align_volnext] == 0))
  {
    align_volnext++;
  }
}

auto align_threads_done() -> void
{
  hits_sorted = Buffer<long>();
}

auto align_getwork(long * first, long * last) -> int
{
  int status = 0;

  std::lock_guard<std::mutex> const lock(workmutex);
  if (align_volnext < align_bins)
  {
    long const seqcount = align_volseqs[align_volnext];
    long const chunks = align_volchunks[align_volnext];
    long const chunksize = ((seqcount+chunks-1) / chunks);

    * first = alignedhits;
    * last = alignedhits + chunksize - 1;

    alignedhits += chunksize;
    status = 1;

    align_volseqs[align_volnext] -= chunksize;
    align_volchunks[align_volnext]--;

    while ((align_volnext < align_bins) && (align_volchunks[align_volnext] == 0))
    {
      align_volnext++;
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


auto search_init(Parameters const & parameters, struct search_data * sdp) -> void
{
  sdp->dbt = db_thread_create();
  sdp->dprofile.resize(4 * 16 * 32);
  long qlen = 0;
  long hearraylen = 0;

  if (parameters.symtype == SymbolType::blastn)
  {
    for (int s = 0; s < 2; s++)
    {
      if (searches_strand(parameters.querystrands, s))
      {
	qlen = query.nt[s].len;
	sdp->qlen[3*s] = qlen;
	sdp->qtable[3*s].resize(static_cast<std::size_t>(qlen));
	for (std::size_t i = 0; i < sdp->qtable[3*s].size(); i++)
	{
	  sdp->qtable[3*s][i] = std::next(sdp->dprofile.data(), 64 * query.nt[s].seq[i]);
	}
	hearraylen = qlen > hearraylen ? qlen : hearraylen;
      }
    }
  }
  else if ((parameters.symtype == SymbolType::blastp) || (parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::sound))
  {
    qlen = query.aa[0].len;
    sdp->qlen[0] = qlen;
    sdp->qtable[0].resize(static_cast<std::size_t>(qlen));
    for (std::size_t i = 0; i < sdp->qtable[0].size(); i++)
    {
      sdp->qtable[0][i] = std::next(sdp->dprofile.data(), 64 * query.aa[0].seq[i]);
    }
    hearraylen = qlen > hearraylen ? qlen : hearraylen;
  }
  else if ((parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastx))
  {
    for (int s = 0; s < 2; s++)
    {
      if (searches_strand(parameters.querystrands, s))
      {
	for(int f=0; f<3; f++)
	{
	  qlen = query.aa[(3*s)+f].len;
	  sdp->qlen[(3*s)+f] = qlen;
	  sdp->qtable[(3*s)+f].resize(static_cast<std::size_t>(qlen));
	  for (std::size_t i = 0; i < sdp->qtable[(3*s)+f].size(); i++)
	  {
	    sdp->qtable[(3*s)+f][i] = std::next(sdp->dprofile.data(), 64 * query.aa[(3*s)+f].seq[i]);
	  }
	  hearraylen = qlen > hearraylen ? qlen : hearraylen;
	}
      }
    }
  }
  
  //  fprintf(out, "hearray length = %ld\n", hearraylen);

  sdp->hearray.resize(static_cast<std::size_t>(hearraylen) * 32);

  auto listsize = static_cast<std::size_t>(maxchunksize);
  if ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx))
  {
    listsize *= 6;
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
  
  std::lock_guard<std::mutex> const lock(workmutex);
  if (volnext < volcount)
  {
    long const seqcount = volseqs[volnext];
    long const chunks = volchunks[volnext];
    long const chunksize = ((seqcount+chunks-1) / chunks);

    * first = seqnext;
    * last = seqnext + chunksize - 1;
    seqnext += chunksize;
    status = 1;

    //    fprintf(out, "Processing sequences %d to %d (%d sequences) in volume %ld.\n", *first, *last, *last - * first + 1, volnext);

    volseqs[volnext] -= chunksize;
    volchunks[volnext]--;

    while ((volnext < volcount) && (volchunks[volnext] == 0))
    {
      volnext++;
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
  BYTE const gapopenextend_7 = static_cast<BYTE>(std::min(parameters.gapopenextend, max_7));
  BYTE const gapextend_7 = static_cast<BYTE>(std::min(parameters.gapextend, max_7));

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
      
      BYTE ** qtable = sdp->qtable[(3*qstrand)+qframe].data();
      long const qlen = sdp->qlen[(3*qstrand)+qframe];
      
      /* 7-bit search */
	  
      std::swap(sdp->in_list, sdp->out_list);
      sdp->in_count = sdp->out_count;
	  
      if (sdp->in_count > 0)
      {
	{
	  std::lock_guard<std::mutex> const lock(countmutex);
	  compute7 += static_cast<long>(sdp->in_count);
	}
	    
	// fprintf(out, "Searching seqnos %ld to %ld\n", sdp->in_list[0], sdp->in_list[sdp->in_count-1]);

	if (cpu_feature_ssse3 != 0)
	{
	  search7_ssse3(qtable,
			gapopenextend_7,
			gapextend_7,
			reinterpret_cast<BYTE*>(score_matrix_7t),
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
		  reinterpret_cast<BYTE *>(score_matrix_7),
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
      
	  if (score < SCORELIMIT_7)
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
		 reinterpret_cast<WORD*>(score_matrix_16),
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
	  if (score < SCORELIMIT_16)
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
      
	  char * address = nullptr;
	  long length = 0;
	  long ntlen = 0;
	  db_getsequence(sdp->dbt, seqno, dstrand, dframe, 
			 & address, & length, & ntlen, 0);
	  char * dbegin = address;
	  char const * dend = address + length - 1;
      
	  char * q = nullptr;
	  if (parameters.symtype == SymbolType::blastn)
	  {
	    q = query.nt[qstrand].seq;
	  }
	  else
	  {
	    q = query.aa[(3 * qstrand) + qframe].seq;
	  }

	  long const score = fullsw(dbegin,
			      dend,
			      q, 
			      q + qlen,
			      reinterpret_cast<long*>(sdp->hearray.data()),
			      score_matrix_63,
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


auto prepare_search(long par) -> void
{
  volnext = 0;
  seqnext = 0;

  auto const volcount = static_cast<std::size_t>(db_getvolumecount());
  for (std::size_t v = 0; v < volcount; v++)
  {
    volseqs[v] = db_getseqcount_volume(static_cast<long>(v));
  }

  long totalchunks = 0;

  calc_chunks(static_cast<long>(volcount),
	      par,
	      16,
	      volseqs.data(),
	      volchunks.data(),
	      & totalchunks,
	      & maxchunksize);

  while ((volnext < volcount) && (volchunks[volnext] == 0))
  {
    volnext++;
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

}  // anonymous namespace

#define cpuid(l1,l2,a,b,c,d)						\
  __asm__ __volatile__							\
    ("cpuid": "=a" (a), "=b" (b), "=c" (c), "=d" (d) : "a" (l1), "c" (l2));

namespace {

auto cpu_features() -> void
{
  unsigned int a __attribute__ ((unused)) = 0;
  unsigned int b __attribute__ ((unused)) = 0;
  unsigned int c = 0;
  unsigned int d = 0;
  cpuid(1,0,a,b,c,d);
  cpu_feature_sse2  = (d >> 26) & 1;
  cpu_feature_ssse3 = (c >>  9) & 1;
  cpu_feature_sse41 = (c >> 19) & 1;
}

auto clock_start(struct time_info * tip) -> void
{
  time(& tip->t1);                 /* time(2)   */
  tip->clock1 = std::chrono::steady_clock::now();
}

auto clock_stop(Parameters const & parameters, struct time_info * tip) -> void
{
  struct tm tms;
  char const timeformat[] = "%a, %e %b %Y %T UTC";

  tip->clock2 = std::chrono::steady_clock::now();
  time(& tip->t2);

  gmtime_r(&tip->t1, & tms);
  strftime(tip->starttime.data(), tip->starttime.size(), timeformat, & tms);
  
  gmtime_r(&tip->t2, & tms);
  strftime(tip->endtime.data(), tip->endtime.size(), timeformat, & tms);

  tip->elapsed = std::chrono::duration<double>(tip->clock2 - tip->clock1).count();
  
  double speed = (static_cast<double>(db_getsymcount_masked()));

  if (parameters.symtype == SymbolType::blastn)
  {
    speed *= static_cast<double>(query.nt[0].len);
    if (parameters.querystrands == QueryStrands::both)
    {
      speed *= 2;
    }
  }
  else if ((parameters.symtype == SymbolType::blastp) || (parameters.symtype == SymbolType::sound))
  {
    /* sound queries are stored as amino acid queries (KI-33) */
    speed *= static_cast<double>(query.aa[0].len);
  }
  else if (parameters.symtype == SymbolType::blastx)
  {
    speed *= static_cast<double>(query.nt[0].len);
    if (parameters.querystrands == QueryStrands::both)
    {
      speed *= 2;
    }
  }
  else if (parameters.symtype == SymbolType::tblastn)
  {
    speed *= 2;
    speed *= static_cast<double>(query.aa[0].len);
  }
  else if (parameters.symtype == SymbolType::tblastx)
  {
    speed *= 2;
    speed *= static_cast<double>(query.nt[0].len);
    if (parameters.querystrands == QueryStrands::both)
    {
      speed *= 2;
    }
  }
  /* the speed is unknown when no time elapsed (KI-33) */
  tip->speed = (tip->elapsed > 0.0) ? speed / tip->elapsed : 0.0;
  
  if (parameters.view == OutputFormat::plain)
  {
    fprintf(out, "Search started:    %s\n", tip->starttime.data());
    fprintf(out, "Search completed:  %s\n", tip->endtime.data());
    fprintf(out, "Elapsed:           %.2fs\n", tip->elapsed);
    if (tip->elapsed > 0.0)
    {
      fprintf(out, "Speed:             %.3f GCUPS\n", tip->speed / 1e9);
    }
    else
    {
      fprintf(out, "Speed:             n/a\n");
    }
    fprintf(out, "\n");
  }
}



auto work(Parameters const & parameters) -> void
{
  args_show(parameters);
  hits_init(parameters);

  compute7 = 0;

  //  totalhits = 0;

  prepare_search(parameters.threads);

  if (parameters.view==OutputFormat::plain)
  {
    fprintf(out, "Searching...");
    fflush(out);
  }

  clock_start(&ti);
  
  run_threads(parameters);
 
  if (parameters.view == OutputFormat::plain)
  {
    fprintf(out, "...............................................done\n\n");
  }
 
  clock_stop(parameters, &ti);

  //  if (view == 0)
  //    clock_start(&ti);

  align_threads(parameters);
  
  //  if (view == 0)
  //    clock_stop(&ti);

  hits_show(parameters);
  hits_exit();
}

}  // anonymous namespace

auto main(int argc, char**argv) -> int
{

  cpu_features();

  if (cpu_feature_sse2 == 0)
  {
    fatal("This program requires a processor with SSE2.");
  }

  auto const parameters = args_init(argc, argv);

  db_open(parameters);
  
  volchunks.resize(static_cast<std::size_t>(db_getvolumecount()));
  volseqs.resize(static_cast<std::size_t>(db_getvolumecount()));

  if(parameters.dump != 0)
  {
    struct db_thread_s * t = db_thread_create();
    long const seqcount = db_getseqcount();
    for (long i = 0; i < seqcount; i++)
    {
      db_show_fasta(t, i, 0, 0, parameters.dump - 1);
    }
    db_thread_destruct(t);
  }
  else
  {
    score_matrix_init(parameters);

    queryno = 0;
    
    query_init(parameters.queryname, parameters.symtype, parameters.querystrands);
    
    {
      hits_show_begin(parameters.view);
    }
    
    while (query_read() != 0)
    {
      
      work(parameters);
      
      queryno++;
    }
    
    {
      hits_show_end(parameters.view);
    }
    
    query_exit();
  }
  
  db_close();

  if (parameters.outfile != nullptr)
  {
    fclose(out);
  }
}
