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
#include "search_data.h"
#include <algorithm>  // std::generate
#include <array>
#include <cstddef>  // std::size_t
#include <cstdlib>  // exit, posix_memalign
#include <functional>  // std::cref
#include <iterator>  // std::begin, std::end, std::next
#include <mutex>  // std::mutex, std::lock_guard
#include <string>  // std::string (fatal)
#include <thread>
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
std::mutex workmutex;
long maxchunksize;

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
