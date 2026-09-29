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
#include "search_data.h"  // prepare_search, run_threads, align_threads
#include "print_view.h"  // as_c_string, fprint
#include <cassert>
#include <cstddef>  // size_t
#include <cstdio>  // std::fclose, std::ferror, std::fflush
#include <cstdlib>  // exit, posix_memalign
#include <string>  // std::string (fatal)


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

// the four registers set by the cpuid instruction
struct CpuidRegisters
{
  unsigned int eax = 0;
  unsigned int ebx = 0;
  unsigned int ecx = 0;
  unsigned int edx = 0;
};

auto cpuid(unsigned int const leaf, unsigned int const subleaf) -> CpuidRegisters
{
  CpuidRegisters registers;
  __asm__ __volatile__
    ("cpuid" : "=a" (registers.eax), "=b" (registers.ebx),
     "=c" (registers.ecx), "=d" (registers.edx) : "a" (leaf), "c" (subleaf));
  return registers;
}

auto cpu_features() -> void
{
  CpuidRegisters const registers = cpuid(1, 0);
  cpu_feature_sse2  = (registers.edx >> 26) & 1;
  cpu_feature_ssse3 = (registers.ecx >>  9) & 1;
  cpu_feature_sse41 = (registers.ecx >> 19) & 1;
}

auto clock_start(struct time_info * tip) -> void
{
  static_cast<void>(time(& tip->t1));  /* time(2)   */
  tip->clock1 = std::chrono::steady_clock::now();
}

auto clock_stop(Parameters const & parameters, struct time_info * tip) -> void
{
  struct tm tms;
  char const timeformat[] = "%a, %e %b %Y %T UTC";

  tip->clock2 = std::chrono::steady_clock::now();
  static_cast<void>(time(& tip->t2));

  // strftime() returns 0 when the buffer is too small (its contents
  // are then undefined): the buffers hold the longest date
  gmtime_r(&tip->t1, & tms);
  auto const start_length = strftime(tip->starttime.data(), tip->starttime.size(), timeformat, & tms);
  assert(start_length != 0);
  static_cast<void>(start_length);
  
  gmtime_r(&tip->t2, & tms);
  auto const end_length = strftime(tip->endtime.data(), tip->endtime.size(), timeformat, & tms);
  assert(end_length != 0);
  static_cast<void>(end_length);

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
    fprint(out, "Search started:    ");
    fprint(out, as_c_string(tip->starttime.data()));
    fprint(out, '\n');
    fprint(out, "Search completed:  ");
    fprint(out, as_c_string(tip->endtime.data()));
    fprint(out, '\n');
    fprintf(out, "Elapsed:           %.2fs\n", tip->elapsed);
    if (tip->elapsed > 0.0)
    {
      fprintf(out, "Speed:             %.3f GCUPS\n", tip->speed / 1e9);
    }
    else
    {
      fprint(out, "Speed:             n/a\n");
    }
    fprint(out, '\n');
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
    fprint(out, "Searching...");
    static_cast<void>(fflush(out));  // a write error is reported at the end (main())
  }

  clock_start(&ti);
  
  run_threads(parameters);
 
  if (parameters.view == OutputFormat::plain)
  {
    fprint(out, "...............................................done\n\n");
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

  /* A write error (full disk, quota, broken pipe) is often deferred by
     stdio until the buffer is flushed, so check fflush and the error
     flag before closing; fclose also flushes and can report the same
     error (as vsearch's CheckedCloseOutputHandle) */
  if ((std::fflush(out) != 0) or (std::ferror(out) != 0))
  {
    fatal("Unable to write to output file (disk full, quota exceeded, or broken pipe?)");
  }
  if ((parameters.outfile != nullptr) and (std::fclose(out) != 0))
  {
    fatal("Unable to close output file (disk full or quota exceeded?)");
  }
}
