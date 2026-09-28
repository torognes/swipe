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
#include <algorithm>  // std::min
#include <cassert>
#include <cerrno>  // errno, ERANGE
#include <cmath>  // std::floor, std::isfinite
#include <cstdlib>  // std::strtol, std::strtod
#include <limits>
#include <string>  // std::string (fatal)
#include <vector>

/* ARGUMENTS AND THEIR DEFAULTS */

#define DEFAULT_MAXMATCHES 250
#define DEFAULT_ALIGNMENTS 100
#define DEFAULT_MINSCORE 1
#define DEFAULT_MAXSCORE (LONG_MAX)
#define DEFAULT_QUERYNAME "-"
#define DEFAULT_DATABASENAME ""
#define DEFAULT_GAPOPEN 0
#define DEFAULT_GAPEXTEND 0
#define DEFAULT_MATRIXNAME "BLOSUM62"
#define DEFAULT_MATCHSCORE 1
#define DEFAULT_MISMATCHSCORE (-3)
#define DEFAULT_THREADS 1
#define DEFAULT_VIEW 0
#define DEFAULT_SYMTYPE 1
#define DEFAULT_SHOW_GIS 0
#define DEFAULT_SHOW_TAXID 0
#define DEFAULT_EXPECT 10.0
#define DEFAULT_MINEXPECT 0.0
#define DEFAULT_QUERYSTRANDS 3
#define DEFAULT_QUERY_GENCODE 1
#define DEFAULT_DB_GENCODE 1
#define DEFAULT_SUBALIGNMENTS 1
#define DEFAULT_DUMP 0
#define DEFAULT_OUT stdout
#define DEFAULT_EFFDBSIZE 0

char * progname;
const char * matrixname;
const char * databasename;
const char * queryname;
char * taxidfilename;
char * outfile = nullptr;

double expect;
double minexpect;
long minscore;
long maxscore;
long alignments;
long maxmatches;
long gapopen;
long gapextend;
long threads;
long view;
long symtype;
long show_gis;
long show_taxid;
long matchscore;
long mismatchscore;
long gapopenextend;
long querystrands;
long query_gencode;
long db_gencode;
long subalignments;
long dump;
long effdbsize;

/* Other variables */

long queryno;

long cpu_feature_sse2;
long cpu_feature_ssse3;
long cpu_feature_sse41;

#define MAX_THREADS 256

pthread_mutex_t countmutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t workmutex = PTHREAD_MUTEX_INITIALIZER;

pthread_t pthread_id[MAX_THREADS];


long maxchunksize;
long volnext;
long seqnext;
long * volchunks;
long * volseqs;

long compute7;
long compute16;
long compute32;
long compute63;

long rounds7;
long rounds16;
long rounds32;
long rounds63;

long totalhits;

FILE * out = DEFAULT_OUT;

struct time_info ti;

struct search_data
{
  struct db_thread_s * dbt;
  struct db_thread_s * dbta[8];

  BYTE * dprofile;
  BYTE * hearray;
  BYTE ** qtable[6];

  long * scores;
  long * bestpos;
  long * bestq;
  long * start_list;
  long * in_list;
  long * out_list;
  long * tmp_list;

  long qlen[6];

  long start_count;
  long in_count;
  long out_count;

  long * start_hits;

  long seqfirst, seqlast;

  long qstrand1, qstrand2, qframe1, qframe2;
  long dstrand1, dstrand2, dframe1, dframe2;
};

void fatal(const char * message)
{
  if (message)
    fprintf(stderr, "%s\n", message);
  exit(1);
}

void fatal(std::string const & message)
{
  fatal(message.c_str());
}

void * xmalloc(size_t size)
{
  const size_t alignment = 16;
  void * t = nullptr;
  if (posix_memalign(& t, alignment, size) != 0)
    t = nullptr;
  
  if (t==nullptr)
    fatal("Unable to allocate enough memory.");

  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  void * t = realloc(ptr, size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

long alignedhits;
long alignedall;
long * hits_sorted;

long align_volnext;
long align_seqnext;

long * align_volseqs;
long * align_volchunks;

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

void align_init(struct search_data * sdp)
{
  sdp->dbt = db_thread_create();

  for(int i=0; i<8; i++)
    sdp->dbta[i] = db_thread_create();

  sdp->dprofile = (BYTE*) xmalloc(4*16*32);
  long qlen = 0;
  long hearraylen = 0;
  
  for(int i = 0; i < 6; i++)
  {
    sdp->qtable[i] = nullptr;
  }

  if (symtype == 0)
  {
    for(int s=0; s<2; s++)
      if ((s+1) & querystrands)
      {
	qlen = query.nt[s].len;
	sdp->qlen[3*s] = qlen;
	sdp->qtable[3*s] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
	for(int i=0; i<qlen; i++)
	{
	  sdp->qtable[3*s][i] = sdp->dprofile + 16*query.nt[s].seq[i];
	}
	hearraylen = qlen > hearraylen ? qlen : hearraylen;
      }
  }
  else if ((symtype == 1) || (symtype == 3) || (symtype == 5))
  {
    qlen = query.aa[0].len;
    sdp->qlen[0] = qlen;
    sdp->qtable[0] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
    for(int i=0; i<qlen; i++)
    {
      sdp->qtable[0][i] = sdp->dprofile + 16*query.aa[0].seq[i];
    }
    hearraylen = qlen > hearraylen ? qlen : hearraylen;
  }
  else if ((symtype == 2) || (symtype == 4))
  {
    for(int s=0; s<2; s++)
      if ((s+1) & querystrands)
	for(int f=0; f<3; f++)
	{
	  qlen = query.aa[3*s+f].len;
	  sdp->qlen[3*s+f] = qlen;
	  sdp->qtable[3*s+f] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
	  for(int i=0; i<qlen; i++)
	  {
	    sdp->qtable[3*s+f][i] = sdp->dprofile + 16*query.aa[3*s+f].seq[i];
	  }
	  hearraylen = qlen > hearraylen ? qlen : hearraylen;
	}
  }
  
  //  fprintf(out, "hearray length = %ld\n", hearraylen);

  sdp->hearray = (BYTE*) xmalloc(hearraylen*32);

  long listsize = maxchunksize * sizeof(long);
  //  if ((symtype == 3) || (symtype == 4))
  //    listsize *= 6;

  sdp->start_list = (long*) xmalloc(listsize);
  sdp->start_hits = (long*) xmalloc(listsize);
  sdp->in_list = (long*) xmalloc(listsize);
  sdp->out_list = (long*) xmalloc(listsize);
  sdp->scores = (long*) xmalloc(listsize);
  sdp->bestpos = (long*) xmalloc(listsize);
  sdp->bestq = (long*) xmalloc(listsize);

  if (symtype == 0)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (symtype == 2)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 2;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (symtype == 3)
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
  else if (symtype == 4)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
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

void align_chunk(struct search_data * sdp, long hitfirst, long hitlast)
{
  if (hitlast < alignments)
  {

    for (long qstrand = sdp->qstrand1; qstrand <= sdp->qstrand2; qstrand++)
      for(long qframe = sdp->qframe1; qframe <= sdp->qframe2; qframe++)
      {
	sdp->start_count = 0;

	for(long hitno = hitfirst; hitno <= hitlast; hitno++)
	{
	  long hs = hits_sorted[hitno];
	  long seqno, score, hqstrand, hqframe, hdstrand, hdframe;
	
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

	if (sdp->start_count)
	{
	  //	  printf("Aligning %ld sequences.\n", sdp->start_count);


	  BYTE ** qtable = sdp->qtable[3*qstrand+qframe];
	  long qlen = sdp->qlen[3*qstrand+qframe];
      
	  /* 16-bit search, 8x1 db symbols, with alignment end */
	  
	  pthread_mutex_lock(&countmutex);
	  compute32 += sdp->in_count;
	  rounds32++;
	  pthread_mutex_unlock(&countmutex);
	
	  search16s((WORD**)qtable,
		    gapopenextend,
		    gapextend,
		    (WORD*)(score_matrix_16),
		    (WORD*)(sdp->dprofile),
		    (WORD*)(sdp->hearray),
		    sdp->dbta,
		    sdp->start_count,
		    sdp->start_list,
		    sdp->scores,
		    sdp->bestpos,
		    sdp->bestq,
		    qlen);
	
	  for (int i=0; i<sdp->start_count; i++)
	  {
	    long pos = sdp->bestpos[i];
	    long bestq = sdp->bestq[i];
	  
	    //	  fprintf(out, "seqno=%ld score=%ld bestpos=%ld\n", seqno, score, pos);
	  
	    long hitno = sdp->start_hits[i];
	  
	    if (sdp->scores[i] < SCORELIMIT_16)
	      hits_enter_align_hint(hitno, bestq, pos);
	  }
	}
      }

  }

  for(long hitno = hitfirst; hitno <= hitlast; hitno++)
    hits_align(sdp->dbt, hits_sorted[hitno]);
}

void align_done(struct search_data * sdp)
{
  for(int i = 0; i < 6; i++)
  {
    if (sdp->qtable[i])
      free(sdp->qtable[i]);
  }

  free(sdp->dprofile);
  free(sdp->hearray);
  free(sdp->scores);
  free(sdp->bestpos);
  free(sdp->bestq);
  free(sdp->start_list);
  free(sdp->start_hits);
  free(sdp->in_list);
  free(sdp->out_list);
  
  for(int i=0; i<8; i++)
    db_thread_destruct(sdp->dbta[i]);

  db_thread_destruct(sdp->dbt);
}


void calc_chunks(long volcount, 
		 long par,
		 long channels,
		 long * volume_sequences,
		 long * volume_chunks,
		 long * totalchunks,
		 long * biggestchunk)
{
#ifdef DEBUG
  fprintf(out, "Calculating chunk distribution among volumes\n");
#endif

  long volsused = 0;
  std::vector<long> chunksizes(volcount);
  long totalseqs = 0;
  long biggest_chunk_size = 0;
  long vv = 0;
  for(long v = 0; v < volcount; v++)
  {
    if (volume_sequences[v])
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
    upper *= (long) (floor(sqrt((1.0 * totalseqs) / (channels * par))));

  long chunks = volsused;
  long minchunks = totalseqs < par ? totalseqs : par;

  while((biggest_chunk_size > upper) || (chunks < minchunks))
  {
    volume_chunks[vv]++;
    chunks++;
    chunksizes[vv] = (volume_sequences[vv] + volume_chunks[vv] - 1) / volume_chunks[vv];

    biggest_chunk_size = 0;
    vv = 0;
    for(long v=0; v < volcount; v++)
    {
      if (chunksizes[v] > biggest_chunk_size)
      {
	vv = v;
	biggest_chunk_size = chunksizes[v];
      }
    }
  }
  
#ifdef DEBUG
  {
    fprintf(out, "\n");
    fprintf(out, "Chunk distribution:\n");
    fprintf(out, "Total seqs: %ld\n", totalseqs);
    fprintf(out, "Total chunks: %ld  Upper size limit: %ld\n", chunks, upper);
    fprintf(out, "Vol Sequences Chunks Maxchunksize\n");
    for(long v = 0; v < volcount; v++)
    {
      fprintf(out, "%3ld %9ld %6ld %12ld\n",
	      v, volume_sequences[v], volume_chunks[v], chunksizes[v]);
    }
    fprintf(out, "\n");
  }
#endif
  
  *biggestchunk = biggest_chunk_size;
  *totalchunks = chunks;
}

void align_threads_init()
{
  long hits = hits_getcount();

  hits_sorted = hits_sort();

  long bins = 7;

  align_volseqs = (long*) xmalloc(bins*sizeof(long));
  align_volchunks = (long*) xmalloc(bins*sizeof(long));

  for(long i = 0; i<bins; i++)
    align_volseqs[i] = 0;

  for(long i = 0; i<hits; i++)
  {
    long seqno;
    long score;
    long qstrand;
    long qframe;
    long dstrand;
    long dframe;
    
    if (i>=alignments)
      align_volseqs[6]++;
    else
    {
      hits_gethit(i, & seqno, & score,
		  & qstrand, & qframe,
		  & dstrand, & dframe);
      
      align_volseqs[3*qstrand+qframe]++;
    }
  }

  long totalchunks;

  calc_chunks(bins,
	      threads,
	      8,
	      align_volseqs,
	      align_volchunks,
	      & totalchunks,
	      & maxchunksize);

  alignedhits = 0;
  align_volnext = 0;

  while ((align_volnext < bins) && (align_volchunks[align_volnext] == 0))
    align_volnext++;
}

void align_threads_done()
{
  free(hits_sorted);
  free(align_volchunks);
  free(align_volseqs);
}

int align_getwork(long * first, long * last)
{
  int status = 0;
  long bins = 7;
  long volcount = bins;

  pthread_mutex_lock(&workmutex);
  if (align_volnext < volcount)
  {
    long seqcount = align_volseqs[align_volnext];
    long chunks = align_volchunks[align_volnext];
    long chunksize = ((seqcount+chunks-1) / chunks);

    * first = alignedhits;
    * last = alignedhits + chunksize - 1;

    alignedhits += chunksize;
    status = 1;

    align_volseqs[align_volnext] -= chunksize;
    align_volchunks[align_volnext]--;
    
    while ((align_volnext < bins) && (align_volchunks[align_volnext] == 0))
      align_volnext++;
  }
  pthread_mutex_unlock(&workmutex);
  return status;
}

void * align_worker(void *)
{
  search_data sd;
  align_init(&sd);

  long i, j;
  while(align_getwork(&i, &j))
    align_chunk(&sd, i, j);
  
  align_done(&sd);
  return nullptr;
}

void align_threads()
{
  long t;
  void * status;

  align_threads_init();
  
  for(t=0; t<threads; t++)
    {
      if (pthread_create(pthread_id + t, nullptr, align_worker, (void *)t))
	fatal("Cannot create thread.");
    }
  
  for(t=0; t<threads; t++) {
    if (pthread_join(pthread_id[t], &status))
      fatal("Cannot join thread.");
  }

  align_threads_done();
}

void args_show()
{
  if (view == 0)
  {
    
    if (! cpu_feature_ssse3)
    {
      fprintf(out, "The performance is reduced because this CPU lacks SSSE3.\n\n");
    }
    
    const char * symtypestring[] = { "Nucleotide", "Amino acid", "Translated query", "Translated database", "Both translated", "Sound" };
    
    //      char * viewtypestring[] = { "plain", 0, 0, 0, 0, 0, 0, "xml",
    //			  "tab-separated", "tab-separated with comments" };
    
    fprintf(out, "Database file:     %s\n", databasename);
    fprintf(out, "Database title:    %s\n", db_gettitle());
    fprintf(out, "Database time:     %s\n", db_gettime());
    
    if (db_ismasked())
      {
	fprintf(out, "Database size:     %ld residues", db_getsymcount_masked());
	fprintf(out, " in %ld sequences\n", db_getseqcount_masked());
      }
      else
      {
	fprintf(out, "Database size:     %ld residues", db_getsymcount());
	fprintf(out, " in %ld sequences\n", db_getseqcount());
      }

      fprintf(out, "Longest db seq:    %ld residues\n", db_getlongest());

      if (effdbsize > 0)
	fprintf(out, "Effective db size: %ld\n", effdbsize);

      fprintf(out, "Query file name:   %s\n", queryname);

      long qlen;
      if ((symtype == 0) || (symtype == 2) || (symtype == 4))
	qlen = query.nt[0].len;
      else
	qlen = query.aa[0].len;
      
      fprintf(out, "Query length:      %ld residues\n", qlen);

      query_show();

      if (symtype == 0)
      {
	fprintf(out, "Query strands:     ");
	switch (querystrands)
	{
	case 1:
	  fprintf(out, "Plus");
	  break;
	case 2:
	  fprintf(out, "Minus");
	  break;
	case 3:
	  fprintf(out, "Plus and minus");
	  break;
	default:
	  break;
	}
	fprintf(out, "\n");
	fprintf(out, "Score matrix:      %ld/%ld\n", matchscore, mismatchscore);
      }
      else
      {
	fprintf(out, "Score matrix:      %s\n", matrixname);
      }

      fprintf(out, "Gap penalty:       %ld+%ldk\n", gapopen, gapextend);
      fprintf(out, "Max expect shown:  %-g\n", expect);
      fprintf(out, "Min score shown:   %ld\n", minscore);
      fprintf(out, "Max matches shown: %ld\n", maxmatches);
      fprintf(out, "Alignments shown:  %ld\n", alignments);
      fprintf(out, "Show gi's:         %ld\n", show_gis);
      fprintf(out, "Show taxid's:      %ld\n", show_taxid);
      fprintf(out, "Threads:           %ld\n", threads);
      fprintf(out, "Symbol type:       %s\n", symtypestring[symtype]);
      if ((symtype == 2) || (symtype == 4))
	fprintf(out, "Query genetic code:%s (%ld)\n", gencode_names[query_gencode-1], query_gencode);
      if ((symtype == 3) || (symtype == 4))
	fprintf(out, "DB genetic code:   %s (%ld)\n", gencode_names[db_gencode-1], db_gencode);
      
#if 0
      if ((symtype == 2) || (symtype == 4))
      {
	for(long s=0; s<2; s++)
	  for(long f=0; f<3; f++)
	  {
	    fprintf(out, "Translation of query, frame %c%ld:\n", s?'-':'+', f+1);
	    //	    translate(query_, qlen, s, f, 0, & prot, & plen);
	    
	    long plen = query.aa[3*s+f].len;
	    for (int j=0; j<plen; j+=60)
	    {
	      for(int k=0; k<60; k++)
	      {
		if (j+k < plen)
		  putc(sym_ncbi_aa[query.aa[3*s+f].seq[j+k]], out);
		else
		  break;
	      }
	      fprintf(out, "\n");
	    }
	    fprintf(out, "\n");
	  }
      }
#endif

      // fprintf(out, "View:              %s\n", viewtypestring[view]);
      if (taxidfilename)
	fprintf(out, "Taxid filename:    %s\n", taxidfilename);
      fprintf(out, "\n");
    }
}
  
void args_usage()
{
  /* options unused by BLAST: chkuxHN */
  /* options used by SWIPE:   chkuxHN  */

  fprintf(out, "Usage: %s [OPTIONS]\n", progname);
  fprintf(out, "  -h, --help                 show help\n");
  fprintf(out, "      --version              show version\n");
  fprintf(out, "  -d, --db=FILE              sequence database base name (required)\n");
  fprintf(out, "  -i, --query=FILE           query sequence filename (stdin)\n");
  fprintf(out, "  -M, --matrix=NAME/FILE     score matrix name or filename (BLOSUM62)\n");
  fprintf(out, "  -q, --penalty=NUM          penalty for nucleotide mismatch (-3)\n");
  fprintf(out, "  -r, --reward=NUM           reward for nucleotide match (1)\n");
  fprintf(out, "  -G, --gapopen=NUM          gap open penalty (11)\n");
  fprintf(out, "  -E, --gapextend=NUM        gap extension penalty (1)\n");
  fprintf(out, "  -v, --num_descriptions=NUM sequence descriptions to show (250)\n");
  fprintf(out, "  -b, --num_alignments=NUM   sequence alignments to show (100)\n");
  fprintf(out, "  -e, --evalue=REAL          maximum expect value of sequences to show (10.0)\n");
  fprintf(out, "  -k, --minevalue=REAL       minimum expect value of sequences to show (0.0)\n");
  fprintf(out, "  -c, --min_score=NUM        minimum score of sequences to show (1)\n");
  fprintf(out, "  -u, --max_score=NUM        maximum score of sequences to show (inf.)\n");
  fprintf(out, "  -a, --num_threads=NUM      number of threads to use [1-%d] (1)\n", MAX_THREADS);
  fprintf(out, "  -m, --outfmt=NUM           output format [0,7-9=plain,xml,tsv,tsv+] (0)\n");
  fprintf(out, "  -I, --show_gis             show gi numbers in results (no)\n");
  fprintf(out, "  -p, --symtype=NAME/NUM     symbol type/translation [0-4] (1)\n");
  fprintf(out, "  -S, --strand=NAME/NUM      query strands to search [1-3] (3)\n");
  fprintf(out, "  -Q, --query_gencode=NUM    query genetic code [1-23] (1)\n");
  fprintf(out, "  -D, --db_gencode=NUM       database genetic code [1-23] (1)\n");
  fprintf(out, "  -x, --taxidlist=FILE       taxid list filename (none)\n");
  fprintf(out, "  -N, --dump=NUM             dump database [0-2=no,yes,split headers] (0)\n");
  fprintf(out, "  -H, --show_taxid           show taxid etc in results (no)\n");
  fprintf(out, "  -o, --out=FILE             output file (stdout)\n");
  fprintf(out, "  -z, --dbsize=NUM           set effective database size (0)\n");
}

void args_version()
{
  char const title[] = "SWIPE " SWIPE_VERSION;
  char const ref[] = "Reference: T. Rognes (2011) Faster Smith-Waterman database searches\nwith inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.";
  fprintf(out, "%s\n\n%s\n", title, ref);
}

void args_help()
{
  args_version();
  fprintf(out, "\n");
  
  args_usage();
}

// strict conversions of option values (KI-9): the whole value must be
// a number, without trailing characters, and within the range of the
// type; otherwise, swipe stops with the error message of the option
auto parse_long(char const * const text, char const * const message) -> long
{
  assert(text != nullptr);
  char * end = nullptr;
  errno = 0;
  auto const value = std::strtol(text, &end, 10);
  if ((end == text) or (*end != '\0') or (errno == ERANGE))
    fatal(message);
  return value;
}

auto parse_double(char const * const text, char const * const message) -> double
{
  assert(text != nullptr);
  char * end = nullptr;
  errno = 0;
  auto const value = std::strtod(text, &end);
  if ((end == text) or (*end != '\0') or (errno == ERANGE) or
      (not std::isfinite(value)))
    fatal(message);
  return value;
}

// the effective database size accepts the real notation of blastall's
// -z (e.g. 7.06e+06, GitHub #9), but must be a non-negative integer
auto parse_dbsize(char const * const text) -> long
{
  static char const message[] = "Illegal effective db size specified";
  constexpr auto upper_limit = static_cast<double>(std::numeric_limits<long>::max());
  auto const value = parse_double(text, message);
  if ((value < 0.0) or (std::floor(value) < value) or (value >= upper_limit))
    fatal(message);
  return static_cast<long>(value);
}

void args_init(int argc, char **argv)
{
  /* Set defaults */
  gapopen = DEFAULT_GAPOPEN;
  gapextend = DEFAULT_GAPEXTEND;
  matrixname = "";
  queryname = DEFAULT_QUERYNAME;
  databasename = DEFAULT_DATABASENAME;
  minscore = DEFAULT_MINSCORE;
  maxscore = DEFAULT_MAXSCORE;
  maxmatches = DEFAULT_MAXMATCHES;
  alignments = DEFAULT_ALIGNMENTS;
  threads = DEFAULT_THREADS;
  view = DEFAULT_VIEW;
  symtype = DEFAULT_SYMTYPE;
  show_gis = DEFAULT_SHOW_GIS;
  show_taxid = DEFAULT_SHOW_TAXID;
  expect = DEFAULT_EXPECT;
  minexpect = DEFAULT_MINEXPECT;
  taxidfilename = nullptr;
  matchscore = DEFAULT_MATCHSCORE;
  mismatchscore = DEFAULT_MISMATCHSCORE;
  querystrands = DEFAULT_QUERYSTRANDS;
  query_gencode = DEFAULT_QUERY_GENCODE;
  db_gencode = DEFAULT_DB_GENCODE;
  subalignments = DEFAULT_SUBALIGNMENTS;
  dump = DEFAULT_DUMP;
  effdbsize = DEFAULT_EFFDBSIZE;

  progname = argv[0];

  opterr = 1;
  char short_options[] = "d:i:M:q:r:G:E:S:v:b:c:u:e:k:a:m:p:x:C:Q:D:F:K:N:o:z:IHh";

  static struct option long_options[] =
  {
    {"db",               required_argument, nullptr, 'd' },
    {"query",            required_argument, nullptr, 'i' },
    {"matrix",           required_argument, nullptr, 'M' },
    {"penalty",          required_argument, nullptr, 'q' },
    {"reward",           required_argument, nullptr, 'r' },
    {"gapopen",          required_argument, nullptr, 'G' },
    {"gapextend",        required_argument, nullptr, 'E' },
    {"strand",           required_argument, nullptr, 'S' },
    {"num_descriptions", required_argument, nullptr, 'v' },
    {"num_alignments",   required_argument, nullptr, 'b' },
    {"min_score",        required_argument, nullptr, 'c' },
    {"max_score",        required_argument, nullptr, 'u' },
    {"evalue",           required_argument, nullptr, 'e' },
    {"minevalue",        required_argument, nullptr, 'k' },
    {"num_threads",      required_argument, nullptr, 'a' },
    {"outfmt",           required_argument, nullptr, 'm' },
    {"symtype",          required_argument, nullptr, 'p' },
    {"taxidlist",        required_argument, nullptr, 'x' },
    {"taxid",            required_argument, nullptr, 'x' },  /* alias (2.1.1 and older) */
    {"comp_based_stats", required_argument, nullptr, 'C' },
    {"query_gencode",    required_argument, nullptr, 'Q' },
    {"db_gencode",       required_argument, nullptr, 'D' },
    {"filter",           required_argument, nullptr, 'F' },
    {"subalignments",    required_argument, nullptr, 'K' },
    {"dump",             required_argument, nullptr, 'N' },
    {"out",              required_argument, nullptr, 'o' },
    {"dbsize",           required_argument, nullptr, 'z' },
    {"show_gis",         no_argument,       nullptr, 'I' },
    {"show_taxid",       no_argument,       nullptr, 'H' },
    {"help",             no_argument,       nullptr, 'h' },
    {"version",          no_argument,       nullptr, 'V' },
    { nullptr, 0, nullptr, 0 }
  };
  
  int option_index = 0;
  int c;

  // gap penalties not given on the command line take the default
  // values of the score matrix or of the symbol type; a penalty of
  // zero is a valid value (KI-6)
  auto gapopen_given = false;
  auto gapextend_given = false;
  
  while (true)
    {
      c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1)
	break;

      switch(c)
	{
	case 'a':
	  /* threads */
	  threads = parse_long(optarg, "Illegal number of threads specified");
	  break;
	  
	case 'b':
	  /* alignments */
	  alignments = parse_long(optarg, "Illegal number of alignments specified.");
	  break;
	  
	case 'c':
	  /* min score threshold */
	  minscore = parse_long(optarg, "Illegal minimum score specified.");
	  break;
	  
	case 'C':
	  /* composition-based adjustments */
	  if ( (strcasecmp(optarg, "F") != 0) && (strcmp(optarg, "0") != 0) )
	    fatal("Composition-based score adjustments not supported.");
	  break;

	case 'd':
	  /* database */
	  databasename = optarg;
	  break;
	  
	case 'D':
	  /* database genetic code */
	  db_gencode = parse_long(optarg, "Illegal database genetic code specified.");
	  break;
	  
	case 'e':
	  /* evalue */
	  expect = parse_double(optarg, "Illegal expect value specified.");
	  break;
	  
	case 'E':
	  /* gap extend */
	  gapextend = parse_long(optarg, "Illegal gap penalties.");
	  gapextend_given = true;
	  break;
	  
	case 'F':
	  /* filter */
	  if ( (strlen(optarg) != 0) && (strcasecmp(optarg, "F") != 0) )
	    fatal("Query sequence filtering not supported.");
	  break;
	  
	case 'G':
	  /* gap open */
	  gapopen = parse_long(optarg, "Illegal gap penalties.");
	  gapopen_given = true;
	  break;
	  
	case 'h':
	  args_help();
	  exit(0);
	  break;

	case 'V':
	  /* long option only: -v is --num_descriptions */
	  args_version();
	  exit(0);
	  break;
	  
	case 'H':
	  /* show_taxid */
	  show_taxid = 1;
	  break;
	  
	case 'i':
	  /* query */
	  queryname = optarg;
	  break;
	  
	case 'I':
	  /* show_gis */
	  show_gis = 1;
	  break;
	  
	case 'k':
	  /* min evalue threshold */
	  minexpect = parse_double(optarg, "Illegal minimum expect value specified.");
	  break;
	  
	case 'K':
	  /* subalignments */
	  subalignments = parse_long(optarg, "Illegal number of subalignments specified.");
	  break;
	  
	case 'm':
	  /* view */
	  view = parse_long(optarg, "Illegal view type.");
	  break;
	  
	case 'M':
	  /* matrix */
	  matrixname = optarg;
	  break;
	  
	case 'N':
	  /* dump */
	  dump = parse_long(optarg, "Illegal dump mode.");
	  break;
	  
	case 'o':
	  /* output file */
	  outfile = optarg;
	  break;
	  
	case 'p':
	  /* symtype */
	  if (strcmp(optarg, "blastn") == 0)
	    symtype = 0;
	  else if (strcmp(optarg, "blastp") == 0)
	    symtype = 1;
	  else if (strcmp(optarg, "blastx") == 0)
	    symtype = 2;
	  else if (strcmp(optarg, "tblastn") == 0)
	    symtype = 3;
	  else if (strcmp(optarg, "tblastx") == 0)
	    symtype = 4;
	  else if (strcmp(optarg, "sound") == 0)
	    symtype = 5;
	  else
	    symtype = parse_long(optarg, "Illegal symbol type.");
	  break;
	  
	case 'q':
	  /* penalty */
	  mismatchscore = parse_long(optarg, "Illegal mismatch penalty specified.");
	  break;
	  
	case 'Q':
	  /* query genetic code */
	  query_gencode = parse_long(optarg, "Illegal query genetic code specified.");
	  break;
	  
	case 'r':
	  /* reward */
	  matchscore = parse_long(optarg, "Illegal match reward specified.");
	  break;
	  
	case 'S':
	  if (strcmp(optarg, "plus") == 0)
	    querystrands = 1;
	  else if (strcmp(optarg, "minus") == 0)
	    querystrands = 2;
	  else if (strcmp(optarg, "both") == 0)
	    querystrands = 3;
	  else
	    querystrands = parse_long(optarg, "Illegal query strands specified.");
	  break;

	case 'u':
	  /* maxscore */
	  maxscore = parse_long(optarg, "Illegal maximum score specified.");
	  break;
	  
	case 'v':
	  /* max matches shown */
	  maxmatches = parse_long(optarg, "Illegal number of descriptions specified.");
	  break;
	  
	case 'x':
	  /* taxid filename */
	  taxidfilename = optarg;
	  break;
	  
	case 'z':
	  /* effective db size */
	  effdbsize = parse_dbsize(optarg);
	  break;
	  
	case '?':
	default:
	  args_usage();
	  exit(1);
	  break;
	}
    }
  
  long gopen_default;
  long gextend_default;

  if (symtype == 0)
  {
    if (not gapopen_given)
      gapopen = 5;
    if (not gapextend_given)
      gapextend = 2;
  }
  else if (symtype < 5)
  {
    if (strlen(matrixname) == 0)
      matrixname = DEFAULT_MATRIXNAME;

    if (stats_getprefs(matrixname, & gopen_default, & gextend_default))
    {
      if (not gapopen_given)
	gapopen = gopen_default;
      if (not gapextend_given)
	gapextend = gextend_default;
    }
    else
    {
      // no default for this matrix: a penalty not given is zero
      if ((not gapopen_given) && (not gapextend_given))
	fatal("Unknown score matrix. Gap penalties must be specified (-G and -E).");
    }
  }
  else if (symtype == 5)
  {
    if (strlen(matrixname) == 0)
      matrixname = "IDENTITY_5_1";
    if (not gapopen_given)
      gapopen = 15;
    if (not gapextend_given)
      gapextend = 5;
  }

  gapopenextend = gapopen + gapextend;

  if (effdbsize < 0)
    fatal("Illegal effective db size specified");

  if ((threads < 1) || (threads > MAX_THREADS))
    fatal("Illegal number of threads specified");

  if (strlen(databasename) == 0)
    fatal("No database specified.");
  
  if (!((view==0)||(view==7)||(view==8)||(view==9)||(view==99)))
    fatal("Illegal view type.");
  
  if ((symtype < 0) || (symtype > 5))
    fatal("Illegal symbol type.");

  if ((gapopen < 0) || (gapextend < 0) || ((gapopen + gapextend) < 1))
    fatal("Illegal gap penalties.");

  if ((querystrands < 1) || (querystrands > 3))
    fatal("Illegal query strands specified.");

  if ((querystrands == 2) && ((symtype == 1) || (symtype == 3)))
    fatal("Illegal strand specified for protein query.");

  if ((query_gencode < 1)  || (query_gencode > 23) || (! gencode_names[query_gencode-1]))
    fatal("Illegal query genetic code specified.");

  if ((db_gencode < 1) || (db_gencode > 23) || (! gencode_names[db_gencode-1]))
    fatal("Illegal database genetic code specified.");

  if ((dump<0) || (dump>2))
    fatal("Illegal dump mode.");

  /* ranges of the result limits (KI-7, KI-9) */
  if (maxmatches < 0)
    fatal("Illegal number of descriptions specified.");

  if (alignments < 0)
    fatal("Illegal number of alignments specified.");

  /* scores below 1 are not alignments ("Internal error in align
     function.") */
  if (minscore < 1)
    fatal("Illegal minimum score specified.");

  if (maxscore < 0)
    fatal("Illegal maximum score specified.");

  if (expect <= 0.0)
    fatal("Illegal expect value specified.");

  if (minexpect < 0.0)
    fatal("Illegal minimum expect value specified.");

  /* the output file is opened (and truncated) only once all the
     options are checked (KI-8) */
  if (outfile)
  {
    FILE * f = fopen(outfile, "w");
    if (! f)
      fatal("Unable to open output file for writing.");
    out = f;
  }
  
  translate_init(query_gencode, db_gencode);
}

}  // anonymous namespace


void vector_print(BYTE * vector)
{
  for(int i=0; i<16; i++)
    fprintf(out, " %02x", vector[i]);
}

void vector_print_word(WORD * vector)
{
  for(int i=0; i<8; i++)
    fprintf(out, " %04x", vector[i]);
}


namespace {

void search_init(struct search_data * sdp)
{
  sdp->dbt = db_thread_create();
  sdp->dprofile = (BYTE*) xmalloc(4*16*32);
  long qlen = 0;
  long hearraylen = 0;
  
  for(int i = 0; i < 6; i++)
  {
    sdp->qtable[i] = nullptr;
  }

  if (symtype == 0)
  {
    for(int s=0; s<2; s++)
      if ((s+1) & querystrands)
      {
	qlen = query.nt[s].len;
	sdp->qlen[3*s] = qlen;
	sdp->qtable[3*s] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
	for(int i=0; i<qlen; i++)
	{
	  sdp->qtable[3*s][i] = sdp->dprofile + 64*query.nt[s].seq[i];
	}
	hearraylen = qlen > hearraylen ? qlen : hearraylen;
      }
  }
  else if ((symtype == 1) || (symtype == 3) || (symtype == 5))
  {
    qlen = query.aa[0].len;
    sdp->qlen[0] = qlen;
    sdp->qtable[0] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
    for(int i=0; i<qlen; i++)
    {
      sdp->qtable[0][i] = sdp->dprofile + 64*query.aa[0].seq[i];
    }
    hearraylen = qlen > hearraylen ? qlen : hearraylen;
  }
  else if ((symtype == 2) || (symtype == 4))
  {
    for(int s=0; s<2; s++)
      if ((s+1) & querystrands)
	for(int f=0; f<3; f++)
	{
	  qlen = query.aa[3*s+f].len;
	  sdp->qlen[3*s+f] = qlen;
	  sdp->qtable[3*s+f] = (BYTE**) xmalloc(qlen*sizeof(BYTE*));
	  for(int i=0; i<qlen; i++)
	  {
	    sdp->qtable[3*s+f][i] = sdp->dprofile + 64*query.aa[3*s+f].seq[i];
	  }
	  hearraylen = qlen > hearraylen ? qlen : hearraylen;
	}
  }
  
  //  fprintf(out, "hearray length = %ld\n", hearraylen);

  sdp->hearray = (BYTE*) xmalloc(hearraylen*32);

  long listsize = maxchunksize * sizeof(long);
  if ((symtype == 3) || (symtype == 4))
    listsize *= 6;

  sdp->start_list = (long*) xmalloc(listsize);
  sdp->in_list = (long*) xmalloc(listsize);
  sdp->out_list = (long*) xmalloc(listsize);
  sdp->scores = (long*) xmalloc(listsize);
  sdp->bestpos = (long*) xmalloc(listsize);
  sdp->bestq = (long*) xmalloc(listsize);

  if (symtype == 0)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 0;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (symtype == 2)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
    sdp->qframe2 = 2;

    sdp->dstrand1 = 0;
    sdp->dframe1 = 0;
    sdp->dstrand2 = 0;
    sdp->dframe2 = 0;
  }
  else if (symtype == 3)
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
  else if (symtype == 4)
  {
    sdp->qstrand1 = querystrands == 2 ? 1 : 0;
    sdp->qframe1 = 0;
    sdp->qstrand2 = querystrands == 1 ? 0 : 1;
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

void search_done(struct search_data * sdp)
{
  for(int i = 0; i < 6; i++)
  {
    if (sdp->qtable[i])
      free(sdp->qtable[i]);
  }

  free(sdp->dprofile);
  free(sdp->hearray);
  free(sdp->scores);
  free(sdp->bestpos);
  free(sdp->bestq);
  free(sdp->start_list);
  free(sdp->in_list);
  free(sdp->out_list);
  db_thread_destruct(sdp->dbt);
}

int search_getwork(long * first, long * last)
{
  int status = 0;
  long volcount = db_getvolumecount();
  
  pthread_mutex_lock(&workmutex);
  if (volnext < volcount)
  {
    long seqcount = volseqs[volnext];
    long chunks = volchunks[volnext];
    long chunksize = ((seqcount+chunks-1) / chunks);

    * first = seqnext;
    * last = seqnext + chunksize - 1;
    seqnext += chunksize;
    status = 1;

    //    fprintf(out, "Processing sequences %d to %d (%d sequences) in volume %ld.\n", *first, *last, *last - * first + 1, volnext);

    volseqs[volnext] -= chunksize;
    volchunks[volnext]--;
    
    while ((volnext < volcount) && (volchunks[volnext] == 0))
      volnext++;
  }
  pthread_mutex_unlock(&workmutex);
  return status;
}


void search_chunk(struct search_data * sdp)
{
  // the 7-bit engine uses signed bytes: gap penalties are clamped to
  // 127 (KI-11). This is exact: 7-bit scores are in [0, 127], so a
  // penalty of 127 already takes any score down to zero. The 16-bit
  // and 63-bit engines, and the alignments, use the real penalties
  long const max_7 = std::numeric_limits<signed char>::max();
  BYTE const gapopenextend_7 = static_cast<BYTE>(std::min(gapopenextend, max_7));
  BYTE const gapextend_7 = static_cast<BYTE>(std::min(gapextend, max_7));

  //  fprintf(out, "Searching seqnos %ld to %ld\n", sdp->seqfirst, sdp->seqlast);
  
  if(taxidfilename)
    db_mapheaders(sdp->dbt, sdp->seqfirst, sdp->seqlast);
  
  sdp->start_count = 0;
  for(long seqno = sdp->seqfirst; seqno <= sdp->seqlast; seqno++)
  {
    if (db_check_inclusion(sdp->dbt, seqno))
    {
      if ((symtype == 3) || (symtype == 4))
      {
	for(long dstrand = sdp->dstrand1; dstrand <= sdp->dstrand2; dstrand++)
	  for(long dframe = sdp->dframe1; dframe <= sdp->dframe2; dframe++)
	  {
	    sdp->start_list[sdp->start_count++] =
	      (seqno << 3) | (dstrand << 2) | dframe;
	  }
      }
      else
      {
	sdp->start_list[sdp->start_count++] = seqno << 3;
      }
    }
  }
  
  if (sdp->start_count == 0)
    return;
  
  long s1 = sdp->start_list[0] >> 3;
  long s2 = sdp->start_list[sdp->start_count-1] >> 3;
  
  // fprintf(out, "Mapping seqnos %ld to %ld\n", s1, s2);

  db_mapsequences(sdp->dbt, s1, s2);

  for (long qstrand = sdp->qstrand1; qstrand <= sdp->qstrand2; qstrand++)
    for(long qframe = sdp->qframe1; qframe <= sdp->qframe2; qframe++)
    {
      long dstrand, dframe;
      
      sdp->out_count = sdp->start_count;
      memcpy(sdp->out_list, sdp->start_list, sdp->start_count * sizeof(long));
      
      BYTE ** qtable = sdp->qtable[3*qstrand+qframe];
      long qlen = sdp->qlen[3*qstrand+qframe];
      
#if 1
      
      /* 7-bit search */
	  
      sdp->tmp_list = sdp->in_list;
      sdp->in_list = sdp->out_list;
      sdp->out_list = sdp->tmp_list;
      sdp->in_count = sdp->out_count;
	  
      if (sdp->in_count > 0)
      {
	pthread_mutex_lock(&countmutex);
	compute7 += sdp->in_count;
	rounds7++;
	pthread_mutex_unlock(&countmutex);
	    
	// fprintf(out, "Searching seqnos %ld to %ld\n", sdp->in_list[0], sdp->in_list[sdp->in_count-1]);
	    
	if (cpu_feature_ssse3)
	  search7_ssse3(qtable,
			gapopenextend_7,
			gapextend_7,
			(BYTE*) score_matrix_7t,
			sdp->dprofile,
			sdp->hearray,
			sdp->dbt,
			sdp->in_count,
			sdp->in_list,
			sdp->scores,
			qlen);
	else
	  search7(qtable,
		  gapopenextend_7,
		  gapextend_7,
		  (BYTE*) score_matrix_7,
		  sdp->dprofile,
		  sdp->hearray,
		  sdp->dbt,
		  sdp->in_count,
		  sdp->in_list,
		  sdp->scores,
		  qlen);
	
	sdp->out_count = 0;
    
	for (int i=0; i<sdp->in_count; i++)
	{
	  long seqnosf = sdp->in_list[i];
	  long score = sdp->scores[i];
      
	  if (score < SCORELIMIT_7)
	  {
	    long seqno = seqnosf >> 3;
	    dstrand = (seqnosf >> 2) & 1;
	    dframe = seqnosf & 3;

	    if ((symtype == 0) && qstrand)
	      hits_enter(seqno, score, 0, 0, 1, 0, -1, -1);
	    else
	      hits_enter(seqno, score, qstrand, qframe, dstrand, dframe, -1, -1);
	  }
	  else
	  {
	    sdp->out_list[sdp->out_count++] = seqnosf;
	  }
	}
      }
#endif

#if 1

      /* 16-bit search */
	  
      sdp->tmp_list = sdp->in_list;
      sdp->in_list = sdp->out_list;
      sdp->out_list = sdp->tmp_list;
      sdp->in_count = sdp->out_count;
  
      if (sdp->in_count > 0)
      {
	pthread_mutex_lock(&countmutex);
	compute16 += sdp->in_count;
	rounds16++;
	pthread_mutex_unlock(&countmutex);
	  
	search16((WORD**)qtable,
		 gapopenextend,
		 gapextend,
		 (WORD*)(score_matrix_16),
		 (WORD*)(sdp->dprofile),
		 (WORD*)(sdp->hearray),
		 sdp->dbt,
		 sdp->in_count,
		 sdp->in_list,
		 sdp->scores,
		 sdp->bestpos,
		 qlen);
    
	sdp->out_count = 0;
    
	for (int i=0; i<sdp->in_count; i++)
	{
	  long seqnosf = sdp->in_list[i];
	  long score = sdp->scores[i];
	  if (score < SCORELIMIT_16)
	  {
	    long seqno = seqnosf >> 3;
	    dstrand = (seqnosf >> 2) & 1;
	    dframe = seqnosf & 3;
		
	    long pos = sdp->bestpos[i];
	    
	    //	    fprintf(out, "seqno=%ld score=%ld bestpos=%ld\n", seqno, score, pos);

	    if ((symtype == 0) && qstrand)
	      hits_enter(seqno, score, 0, 0, 1, 0, pos, -1);
	    else
	      hits_enter(seqno, score, qstrand, qframe, dstrand, dframe, pos, -1);
	  }
	  else
	  {
	    sdp->out_list[sdp->out_count++] = seqnosf;
	  }
	}
      }
#endif
      
#if 1

      /* 63-bit search */

      sdp->tmp_list = sdp->in_list;
      sdp->in_list = sdp->out_list;
      sdp->out_list = sdp->tmp_list;
      sdp->in_count = sdp->out_count;
  
      if (sdp->in_count > 0)
      {
	pthread_mutex_lock(&countmutex);
	compute63 += sdp->in_count;
	rounds63++;
	pthread_mutex_unlock(&countmutex);
    
	for (int i=0; i<sdp->in_count; i++)
	{
	  long seqnosf = sdp->in_list[i];
	  long seqno = seqnosf >> 3;
	  dstrand = (seqnosf >> 2) & 1;
	  dframe = seqnosf & 3;
      
	  char * address;
	  long length;
	  long ntlen;
	  db_getsequence(sdp->dbt, seqno, dstrand, dframe, 
			 & address, & length, & ntlen, 0);
	  char * dbegin = address;
	  char * dend = address + length - 1;
      
	  char * q;
	  if (symtype == 0)
	    q = query.nt[qstrand].seq;
	  else
	    q = query.aa[3*qstrand+qframe].seq;

	  long score = fullsw(dbegin,
			      dend,
			      (char*) q, 
			      (char*) q + qlen,
			      (long*) sdp->hearray,
			      score_matrix_63,
			      gapopenextend,
			      gapextend);

	  if ((symtype == 0) && qstrand)
	    hits_enter(seqno, score, 0, 0, 1, 0, -1, -1);
	  else
	    hits_enter(seqno, score, qstrand, qframe, dstrand, dframe, -1, -1);
	}
      }
  
#endif
    }
}


void * worker(void *)
{
  struct search_data sd;
  search_init(&sd);

  while(search_getwork(&sd.seqfirst, &sd.seqlast))
    search_chunk(&sd);

  search_done(&sd);
  return nullptr;
}


void prepare_search(long par)
{
  volnext = 0;
  seqnext = 0;

#if 1

  long volcount = db_getvolumecount();
  for(long v = 0; v < volcount; v++)
    volseqs[v] = db_getseqcount_volume(v);

  long totalchunks;

  calc_chunks(volcount,
	      par,
	      16,
	      volseqs,
	      volchunks,
	      & totalchunks,
	      & maxchunksize);

  while ((volnext < volcount) && (volchunks[volnext] == 0))
    volnext++;
  
#else
  long seqcount = db_getseqcount();

  long chunkcount;
  if (seqcount >= 16*par)
  {
    chunkcount = par * (long) (floor(sqrt((1.0 * seqcount) / (16 * par))));
  }
  else if (seqcount >= par)
  {
    chunkcount = par;
  }
  else
  {
    chunkcount = seqcount;
  }

//  fprintf(out, "Chunkcount: %ld\n", chunkcount);

  long volcount = db_getvolumecount();
  long rest_chunks = chunkcount;
  long rest_seqcount = seqcount;
  
  //  fprintf(out, "Seqcount, chunkcount, ratio (total): %ld, %ld, %f\n", seqcount, chunkcount, 1.0 * seqcount / chunkcount);

  maxchunksize = 0;

  for(long i = 0; i < volcount; i++)
  {
    long volsize = db_getseqcount_volume(i);
    long volchunk = ((rest_chunks * volsize) + rest_seqcount - 1) / rest_seqcount;
    
    volseqs[i] = volsize;
    volchunks[i] = volchunk;

    long chunksize = (volsize + volchunk -1) / volchunk;

    if (chunksize > maxchunksize)
      maxchunksize = chunksize;

    // fprintf(out, "Seqcount, chunkcount, ratio (%ld): %ld, %ld, %f\n", i, volsize, volchunk, 1.0 * volsize / volchunk);

    rest_chunks -= volchunk;
    rest_seqcount -= volsize;
  }
#endif
}

void run_threads()
{
  long t;
  void * status;

  for(t=0; t<threads; t++)
    {
      if (pthread_create(pthread_id + t, nullptr, worker, (void *)t))
	fatal("Cannot create thread.");
    }
  
  for(t=0; t<threads; t++) {
    if (pthread_join(pthread_id[t], &status))
      fatal("Cannot join thread.");
  }
}

}  // anonymous namespace

#define cpuid(l1,l2,a,b,c,d)						\
  __asm__ __volatile__							\
    ("cpuid": "=a" (a), "=b" (b), "=c" (c), "=d" (d) : "a" (l1), "c" (l2));

namespace {

void cpu_features()
{
  unsigned int a __attribute__ ((unused));
  unsigned int b __attribute__ ((unused));
  unsigned int c,d;
  cpuid(1,0,a,b,c,d);
  cpu_feature_sse2  = (d >> 26) & 1;
  cpu_feature_ssse3 = (c >>  9) & 1;
  cpu_feature_sse41 = (c >> 19) & 1;
}

void clock_start(struct time_info * tip)
{
  time(& tip->t1);                 /* time(2)   */
  tip->clock1 = std::chrono::steady_clock::now();
}

void clock_stop(struct time_info * tip)
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
  
  double speed = ((double)db_getsymcount_masked());

  if (symtype == 0)
  {
    speed *= query.nt[0].len;
    if (querystrands == 3)
      speed *= 2;
  }
  else if ((symtype == 1) || (symtype == 5))
  {
    /* sound queries are stored as amino acid queries (KI-33) */
    speed *= query.aa[0].len;
  }
  else if (symtype == 2)
  {
    speed *= query.nt[0].len;
    if (querystrands == 3)
      speed *= 2;
  }
  else if (symtype == 3)
  {
    speed *= 2;
    speed *= query.aa[0].len;
  }
  else if (symtype == 4)
  {
    speed *= 2;
    speed *= query.nt[0].len;
    if (querystrands == 3)
      speed *= 2;
  }
  /* the speed is unknown when no time elapsed (KI-33) */
  tip->speed = (tip->elapsed > 0.0) ? speed / tip->elapsed : 0.0;
  
  if (view == 0)
  {
    fprintf(out, "Search started:    %s\n", tip->starttime.data());
    fprintf(out, "Search completed:  %s\n", tip->endtime.data());
    fprintf(out, "Elapsed:           %.2fs\n", tip->elapsed);
    if (tip->elapsed > 0.0)
      fprintf(out, "Speed:             %.3f GCUPS\n", tip->speed / 1e9);
    else
      fprintf(out, "Speed:             n/a\n");
    fprintf(out, "\n");
  }
}



void work()
{
  args_show();
  hits_init(maxmatches, alignments, minscore, maxscore, minexpect, expect, view==0);

  compute7 = 0;
  compute16 = 0;
  compute32 = 0;
  compute63 = 0;
  rounds7 = 0;
  rounds16 = 0;
  rounds32 = 0;
  rounds63 = 0;

  //  totalhits = 0;

  prepare_search(threads);

  if (view==0)
  {
    fprintf(out, "Searching...");
    fflush(out);
  }

  clock_start(&ti);
  
  run_threads();
 
#if 1
  if (view == 0)
    fprintf(out, "...............................................done\n\n");
#endif
 
  clock_stop(&ti);

#if 0
  if (view == 0)
  {
    fprintf(out, "Computed (7bit):   %ld sequences in %ld rounds\n", compute7, rounds7);
    fprintf(out, "Computed (16bit):  %ld sequences in %ld rounds\n", compute16, rounds16);
    fprintf(out, "Computed (63bit):  %ld sequences in %ld rounds\n", compute63, rounds63);
    //    fprintf(out, "Total hits:        %ld\n", totalhits);
    fprintf(out, "\n");
  }
#endif

#if 0
  if (view==0)
  {
    fprintf(out, "Aligning...");
    fflush(out);
  }
#endif

  //  if (view == 0)
  //    clock_start(&ti);

  align_threads();
  
#if 0
  if (view == 0)
    fprintf(out, "...............................................done\n\n");
#endif
 
  //  if (view == 0)
  //    clock_stop(&ti);

  hits_show(view, show_gis);
  hits_exit();
}

}  // anonymous namespace

int main(int argc, char**argv)
{

  cpu_features();

  if (! cpu_feature_sse2)
    fatal("This program requires a processor with SSE2.");

  args_init(argc,argv);

  db_open(symtype, databasename, taxidfilename);
  
  volchunks = (long*) xmalloc(db_getvolumecount() * sizeof(long));
  volseqs   = (long*) xmalloc(db_getvolumecount() * sizeof(long));

  if(dump)
  {
    struct db_thread_s * t = db_thread_create();
    long seqcount = db_getseqcount();
    for (long i=0; i < seqcount; i++)
      db_show_fasta(t, i, 0, 0, dump-1);
    db_thread_destruct(t);
  }
  else
  {
    score_matrix_init();

    queryno = 0;
    
    query_init(queryname, symtype, querystrands);
    
    {
      hits_show_begin(view);
    }
    
    while (query_read())
    {
      
      work();
      
      queryno++;
    }
    
    {
      hits_show_end(view);
    }
    
    query_exit();

#ifdef DEBUG
    fprintf(out, "Freeing score matrix\n");
#endif
    score_matrix_free();
  }
  
  free(volchunks);
  free(volseqs);

#ifdef DEBUG
  fprintf(out, "Closing db\n");
#endif
  db_close();

  
  if (outfile)
    fclose(out);
}
