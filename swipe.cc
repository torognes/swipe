/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2012 Torbjorn Rognes, University of Oslo, 
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

char * progname;
char * matrixname;
char * databasename;
char * queryname;
char * taxidfilename;
char * outfile = NULL;

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

/* Other variables */

long queryno;

long sse2_present;
long ssse3_present;
long sse41_present;

#define MAX_THREADS 256

pthread_mutex_t countmutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t workmutex = PTHREAD_MUTEX_INITIALIZER;

pthread_t pthread_id[MAX_THREADS];

#ifdef MPISWIPE
int mpirank, mpisize;
#endif

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

FILE * out = stdout;

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

void fatal(char * message)
{
  if (message)
    fprintf(stderr, "%s\n", message);
  exit(1);
}

void fatal(char * format, char * message)
{
  fprintf(stderr, format, message);
  fprintf(stderr, "\n");
  exit(1);
}

void * xmalloc(size_t size)
{
  const size_t alignment = 16;
  void * t;
  posix_memalign(& t, alignment, size);
  
  if (t==NULL)
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
    sdp->qtable[i] = 0;
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
  free(sdp->in_list);
  free(sdp->out_list);
  
  for(int i=0; i<8; i++)
    db_thread_destruct(sdp->dbta[i]);

  db_thread_destruct(sdp->dbt);
}


void calc_chunks(long volcount, 
		 long par,
		 long channels,
		 long * volseqs,
		 long * * volchunksp,
		 long * totalchunks,
		 long * biggestchunk)
{
#ifdef DEBUG
  fprintf(out, "Calculating chunk distribution among volumes\n");
#endif

  long volsused = 0;
  long chunksizes[volcount];
  long * volchunks = (long*) xmalloc(volcount * sizeof(long));
  long totalseqs = 0;
  long maxchunksize = 0;
  long vv = 0;
  for(long v = 0; v < volcount; v++)
  {
    if (volseqs[v])
    {
      totalseqs += volseqs[v];
      volsused++;
      volchunks[v] = 1;
      chunksizes[v] = volseqs[v];
      if (chunksizes[v] > maxchunksize)
      {
	maxchunksize = chunksizes[v];
	vv = v;
      }
    }
    else
    {
      chunksizes[v] = 0;
      volchunks[v] = 0;
    }
  }

  long upper = channels;
  if (totalseqs >= 4 * channels * par)
    upper *= (long) (floor(sqrt((1.0 * totalseqs) / (channels * par))));

  long chunks = volsused;
  long minchunks = totalseqs < par ? totalseqs : par;

  while((maxchunksize > upper) || (chunks < minchunks))
  {
    volchunks[vv]++;
    chunks++;
    chunksizes[vv] = (volseqs[vv] + volchunks[vv] - 1) / volchunks[vv];

    maxchunksize = 0;
    vv = 0;
    for(long v=0; v < volcount; v++)
    {
      if (chunksizes[v] > maxchunksize)
      {
	vv = v;
	maxchunksize = chunksizes[v];
      }
    }
  }
  
#ifdef DEBUG
#ifdef MPISWIPE
  if (!mpirank)
#endif
  {
    fprintf(out, "\n");
    fprintf(out, "Chunk distribution:\n");
    fprintf(out, "Total seqs: %ld\n", totalseqs);
    fprintf(out, "Total chunks: %ld  Upper size limit: %ld\n", chunks, upper);
    fprintf(out, "Vol Sequences Chunks Maxchunksize\n");
    for(long v = 0; v < volcount; v++)
    {
      fprintf(out, "%3ld %9ld %6ld %12ld\n",
	      v, volseqs[v], volchunks[v], chunksizes[v]);
    }
    fprintf(out, "\n");
  }
#endif
  
  *volchunksp = volchunks;
  *biggestchunk = maxchunksize;
  *totalchunks = chunks;
}

void align_threads_init()
{
  long hits = hits_getcount();

  hits_sorted = hits_sort();

  long bins = 7;

  align_volseqs = (long*) xmalloc(bins*sizeof(long));

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
  long maxchunksize;

  calc_chunks(bins,
	      threads,
	      8,
	      align_volseqs,
	      & align_volchunks,
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
  return 0;
}

void align_threads()
{
  long t;
  void * status;

  align_threads_init();
  
  for(t=0; t<threads; t++)
    {
      if (pthread_create(pthread_id + t, 0, align_worker, (void *)t))
	fatal("Cannot create thread.");
    }
  
  for(t=0; t<threads; t++) {
    if (pthread_join(pthread_id[t], &status))
      fatal("Cannot join thread.");
  }

  align_threads_done();
}

void args_getstring(int i, int argc, char **argv, char ** result, char * error)
{
  if (i+1 < argc)
    *result = argv[i+1];
  else
    fatal(error);
}

void args_getnum(int i, int argc, char **argv, long * result, char * error)
{
  if (i+1 < argc)
    *result = atol(argv[i+1]);
  else
    fatal(error);
}

void args_show()
{
  if (view == 0)
  {
    
#ifndef MPISWIPE
    if (! ssse3_present)
    {
      fprintf(out, "The performance is reduced because this CPU lacks SSSE3.\n\n");
    }
#endif
    
    char * symtypestring[] = { "Nucleotide", "Amino acid", "Translated query", "Translated database", "Both translated", "Sound" };
    
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
  fprintf(out, "  -k  --minevalue=REAL       minimum expect value of sequences to show (0.0)\n");
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
}

void args_help()
{
  char title[] = "SWIPE " SWIPE_VERSION;
  char ref[] = "Reference: T. Rognes (2011) Faster Smith-Waterman database searches\nwith inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.";
  fprintf(out, "%s [%s %s]\n\n%s\n\n", title, __DATE__, __TIME__, ref);
  
  args_usage();
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
  taxidfilename = NULL;
  matchscore = DEFAULT_MATCHSCORE;
  mismatchscore = DEFAULT_MISMATCHSCORE;
  querystrands = DEFAULT_QUERYSTRANDS;
  query_gencode = DEFAULT_QUERY_GENCODE;
  db_gencode = DEFAULT_DB_GENCODE;
  subalignments = DEFAULT_SUBALIGNMENTS;
  dump = DEFAULT_DUMP;

  progname = argv[0];

  opterr = 1;
  char short_options[] = "d:i:M:q:r:G:E:S:v:b:c:u:e:k:a:m:p:x:C:Q:D:F:K:N:o:IHh";

  static struct option long_options[] =
  {
    {"db",               required_argument, NULL, 'd' },
    {"query",            required_argument, NULL, 'i' },
    {"matrix",           required_argument, NULL, 'M' },
    {"penalty",          required_argument, NULL, 'q' },
    {"reward",           required_argument, NULL, 'r' },
    {"gapopen",          required_argument, NULL, 'G' },
    {"gapextend",        required_argument, NULL, 'E' },
    {"strand",           required_argument, NULL, 'S' },
    {"num_descriptions", required_argument, NULL, 'v' },
    {"num_alignments",   required_argument, NULL, 'b' },
    {"min_score",        required_argument, NULL, 'c' },
    {"max_score",        required_argument, NULL, 'u' },
    {"evalue",           required_argument, NULL, 'e' },
    {"minevalue",        required_argument, NULL, 'k' },
    {"num_threads",      required_argument, NULL, 'a' },
    {"outfmt",           required_argument, NULL, 'm' },
    {"symtype",          required_argument, NULL, 'p' },
    {"taxid",            required_argument, NULL, 'x' },
    {"comp_based_stats", required_argument, NULL, 'C' },
    {"query_gencode",    required_argument, NULL, 'Q' },
    {"db_gencode",       required_argument, NULL, 'D' },
    {"filter",           required_argument, NULL, 'F' },
    {"subalignments",    required_argument, NULL, 'K' },
    {"dump",             required_argument, NULL, 'N' },
    {"out",              required_argument, NULL, 'o' },
    {"show_gis",         no_argument,       NULL, 'I' },
    {"show_taxid",       no_argument,       NULL, 'H' },
    {"help",             no_argument,       NULL, 'h' },
    { 0, 0, 0, 0 }
  };
  
  int option_index = 0;
  int c;
  
  while (1)
    {
      c = getopt_long(argc, argv, short_options, long_options, &option_index);
      if (c == -1)
	break;

      switch(c)
	{
	case 'a':
	  /* threads */
	  threads = atol(optarg);
	  break;
	  
	case 'b':
	  /* alignments */
	  alignments = atol(optarg);
	  break;
	  
	case 'c':
	  /* min score threshold */
	  minscore = atol(optarg);
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
	  db_gencode = atol(optarg);
	  break;
	  
	case 'e':
	  /* evalue */
	  expect = atof(optarg);
	  break;
	  
	case 'E':
	  /* gap extend */
	  gapextend = atol(optarg);
	  break;
	  
	case 'F':
	  /* filter */
	  if ( (strlen(optarg) != 0) && (strcasecmp(optarg, "F") != 0) )
	    fatal("Query sequence filtering not supported.");
	  break;
	  
	case 'G':
	  /* gap open */
	  gapopen = atol(optarg);
	  break;
	  
	case 'h':
	  args_help();
	  exit(1);
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
	  minexpect = atof(optarg);
	  break;
	  
	case 'K':
	  /* subalignments */
	  subalignments = atol(optarg);
	  break;
	  
	case 'm':
	  /* view */
	  view = atol(optarg);
	  break;
	  
	case 'M':
	  /* matrix */
	  matrixname = optarg;
	  break;
	  
	case 'N':
	  /* dump */
	  dump = atol(optarg);
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
	    symtype = atol(optarg);
	  break;
	  
	case 'q':
	  /* penalty */
	  mismatchscore = atol(optarg);
	  break;
	  
	case 'Q':
	  /* query genetic code */
	  query_gencode = atol(optarg);
	  break;
	  
	case 'r':
	  /* reward */
	  matchscore = atol(optarg);
	  break;
	  
	case 'S':
	  if (strcmp(optarg, "plus") == 0)
	    querystrands = 1;
	  else if (strcmp(optarg, "minus") == 0)
	    querystrands = 2;
	  else if (strcmp(optarg, "both") == 0)
	    querystrands = 3;
	  else
	    querystrands = atol(optarg);
	  break;

	case 'u':
	  /* maxscore */
	  maxscore = atol(optarg);
	  break;
	  
	case 'v':
	  /* max matches shown */
	  maxmatches = atol(optarg);
	  break;
	  
	case 'x':
	  /* taxid filename */
	  taxidfilename = optarg;
	  break;
	  
	case '?':
	default:
	  args_usage();
	  exit(1);
	  break;
	}
    }
  
  if (outfile)
  {
    FILE * f = fopen(outfile, "w");
    if (! f)
      fatal("Unable to open output file for writing.");
    out = f;
  }

  long gopen_default;
  long gextend_default;

  if (symtype == 0)
  {
    if (gapopen == 0)
      gapopen = 5;
    if (gapextend == 0)
      gapextend = 2;
  }
  else if (symtype < 5)
  {
    if (strlen(matrixname) == 0)
      matrixname = DEFAULT_MATRIXNAME;

    if (stats_getprefs(matrixname, & gopen_default, & gextend_default))
    {
      if (gapopen == 0)
	gapopen = gopen_default;
      if (gapextend == 0)
	gapextend = gextend_default;
    }
    else
    {
      if ((gapopen == 0) && (gapextend == 0))
	fatal("Unknown score matrix. Gap penalties must be specified (-G and -E).");
    }
  }
  else if (symtype == 5)
  {
    if (strlen(matrixname) == 0)
      matrixname = "IDENTITY_5_1";
    if (gapopen == 0)
      gapopen = 15;
    if (gapextend == 0)
      gapextend = 5;
  }

  gapopenextend = gapopen + gapextend;

  if ((threads < 1) || (threads > MAX_THREADS))
    fatal("Illegal number of threads specified");

  if (strlen(databasename) == 0)
    fatal("No database specified.");
  
  if (!((view==0)||(view==7)||(view==8)||(view==9)||(view==99)))
    fatal("Illegal view type.");
  
  if ((gapopen < 0) || (gapextend < 0) || ((gapopen + gapextend) < 1))
    fatal("Illegal gap penalties.");
  
  if ((symtype < 0) || (symtype > 5))
    fatal("Illegal symbol type.");

  if ((querystrands < 1) || (querystrands > 3))
    fatal("Illegal query strands specified.");

  if ((querystrands == 2) && ((symtype == 1) || (symtype == 3) || (symtype == 4)))
    fatal("Illegal strand specified for protein query.");

  if ((query_gencode < 1)  || (query_gencode > 23) || (! gencode_names[query_gencode-1]))
    fatal("Illegal query genetic code specified.");

  if ((db_gencode < 1) || (db_gencode > 23) || (! gencode_names[db_gencode-1]))
    fatal("Illegal database genetic code specified.");

  if ((dump<0) || (dump>2))
    fatal("Illegal dump mode.");
  
  translate_init(query_gencode, db_gencode);
}


void vector_fill(BYTE * vector, BYTE value)
{
  memset(vector, value, 16);
}

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


void search_init(struct search_data * sdp)
{
  sdp->dbt = db_thread_create();
  sdp->dprofile = (BYTE*) xmalloc(4*16*32);
  long qlen = 0;
  long hearraylen = 0;
  
  for(int i = 0; i < 6; i++)
  {
    sdp->qtable[i] = 0;
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
	    
	search7(qtable,
		gapopenextend,
		gapextend,
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
  return 0;
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
	      & volchunks,
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
      if (pthread_create(pthread_id + t, 0, worker, (void *)t))
	fatal("Cannot create thread.");
    }
  
  for(t=0; t<threads; t++) {
    if (pthread_join(pthread_id[t], &status))
      fatal("Cannot join thread.");
  }
}

#define cpuid(f,a,b,c,d) asm("cpuid": "=a" (a), "=b" (b), "=c" (c), "=d" (d) : "a" (f));


void cpu_features()
{
  unsigned int a __attribute__ ((unused));
  unsigned int b __attribute__ ((unused));
  unsigned int c,d;
  cpuid(1,a,b,c,d);
  //  fprintf(out, "cpuid: %08x %08x %08x %08x\n", a, b, c, d);

  sse2_present  = (d >> 26) & 1;
  ssse3_present = (c >>  9) & 1;
  sse41_present = (c >> 19) & 1;

#if 0
  int has_mmx   = (d >> 23) & 1;
  int has_sse   = (d >> 25) & 1;
  int has_sse3  = (c      ) & 1;
  int has_sse42 = (c >> 20) & 1;
  
  fprintf(out, "MMX:    %d\n", has_mmx);
  fprintf(out, "SSE:    %d\n", has_sse);
  fprintf(out, "SSE2:   %d\n", has_sse2);
  fprintf(out, "SSE3:   %d\n", has_sse3);
  fprintf(out, "SSSE3:  %d\n", has_ssse3);
  fprintf(out, "SSE4.1: %d\n", has_sse41);
  fprintf(out, "SSE4.2: %d\n", has_sse42);
#endif
}

struct time_info ti;

void clock_start(struct time_info * tip)
{
  tip->clk_tck = sysconf(_SC_CLK_TCK);
  time(& tip->t1);                 /* time(2)   */
  tip->wc1 = times(& tip->times1); /* times (2) */
}

void clock_stop(struct time_info * tip)
{
  struct tm tms;
  char buf[30];
  char timeformat[] = "%a, %e %b %Y %T UTC";

  tip->wc2 = times(& tip->times2);
  time(& tip->t2);

  gmtime_r(&tip->t1, & tms);
  strftime(buf, 30, timeformat, & tms);
  tip->starttime = (char*) xmalloc(30);
  strcpy(tip->starttime, buf);
  
  gmtime_r(&tip->t2, & tms);
  strftime(buf, 30, timeformat, & tms);
  tip->endtime = (char*) xmalloc(30);
  strcpy(tip->endtime, buf);

  tip->elapsed = ((double)(tip->wc2 - tip->wc1)) / tip->clk_tck;
  
  double speed = ((double)db_getsymcount_masked());

  if (symtype == 0)
  {
    speed *= query.nt[0].len;
    if (querystrands == 3)
      speed *= 2;
  }
  else if (symtype == 1)
  {
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
  speed /= tip->elapsed;
  tip->speed = speed;
  
  if (view == 0)
  {
    fprintf(out, "Search started:    %s\n", tip->starttime);
    fprintf(out, "Search completed:  %s\n", tip->endtime);
    fprintf(out, "Elapsed:           %.2fs\n", tip->elapsed);
    fprintf(out, "Speed:             %.3f GCUPS\n", tip->speed / 1e9);
    fprintf(out, "\n");
  }
}


#ifdef MPISWIPE

#define tag_die 0
#define tag_search 1
#define tag_search_done 2
#define tag_search_get 3
#define tag_search_report 4
#define tag_align 5
#define tag_align_init 11
#define tag_align_done 6
#define tag_align_get 7
#define tag_align_report 8
#define tag_features 9
#define tag_stats 10
#define tag_align_header 11
#define tag_align_seq 12
#define tag_align_coord 13
#define tag_align_string 14

void master(int size)
{
  //  fprintf(out, "sizeof(long) = %lu\n", sizeof(long));

  long nodes_with_sse2 = sse2_present;
  long nodes_with_ssse3 = ssse3_present;

  for(int i=1; i < size; i++)
  {
    long features;
    MPI_Status s;
    MPI_Recv(&features, 1, MPI_LONG, MPI_ANY_SOURCE, tag_features, 
	     MPI_COMM_WORLD, &s);
    nodes_with_sse2 += features >> 1;
    nodes_with_ssse3 += features & 1;
  }

  if (nodes_with_sse2 < size)
    fatal("Sorry, some nodes lack SSE2.");
  
  if (nodes_with_ssse3 < size)
  {
    fprintf(out, "Performance reduced because %ld nodes lack SSSE3.\n",
	   size - nodes_with_ssse3);
  }
  
  if (size < 2)
    fatal("Cannot run on a single node.");

  args_show();


  /* master node, job distributor */

  //fprintf(out, "This is the master.\n");

  /* init */
  
  hits_init(maxmatches, alignments, minscore, maxscore, minexpect, expect, view==0);

  prepare_search(size-1);

  if (view==0)
    {
      fprintf(out, "Searching...");
      fflush(out);
    }

  long slaves_active = 0;
  long first, last;
  
  long buffersize = 1024;

  //  fprintf(out, "master buffersize: %ld.\n", buffersize);
  
  char * buffer = (char*) xmalloc(buffersize);
  long * longbuffer = (long*) buffer;

  compute7 = 0;
  compute16 = 0;
  compute32 = 0;
  compute63 = 0;
  rounds7 = 0;
  rounds16 = 0;
  rounds32 = 0;
  rounds63 = 0;
  totalhits = 0;

  clock_start(&ti);
  
  /* phase 1 - searching */

  for (int i = 1; i < size; i++)
  {
    if (search_getwork(&first, &last))
    {
      longbuffer[0] = first;
      longbuffer[1] = last;
      /* To be made non-blocking with mpi_isend.
	 Need individual buffers. 
	 Need request list. */
      MPI_Send(longbuffer, 2, MPI_LONG, i, tag_search, MPI_COMM_WORLD);
      //      fprintf(out, "Asking slave %d to search from %ld to %ld.\n", i, first, last);
      slaves_active++;
    }	
    else
      break;
  }
  
  while(slaves_active)
  {
    long morework = search_getwork(&first, &last);

    MPI_Status status;
    int source;
    int tag;
    int len;
    
    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, & status);

    source = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    MPI_Get_count(& status, MPI_CHAR, & len);
    
    if (buffersize < len)
    {
      buffersize = len;
#ifdef DEBUG
      fprintf(out, "resizing buffer to length %ld.\n", buffersize);
#endif
      buffer = (char *) xrealloc(buffer, buffersize);
      longbuffer = (long*) buffer;
    }
    
    MPI_Recv(buffer, buffersize, MPI_CHAR, source, tag, MPI_COMM_WORLD, NULL);
    
#ifdef DEBUG
    fprintf(out, "Received message with tag %d and length %d from node %d.\n",
	   tag, len, source);
#endif

    switch(tag)
    {

    case tag_search_done:

      if (morework)
      {
	longbuffer[0] = first;
	longbuffer[1] = last;
	MPI_Send(longbuffer, 2, MPI_LONG, source, tag_search, MPI_COMM_WORLD);
	//	fprintf(out, "Asking slave %d to search from %ld to %ld.\n", source, first, last);
      }	
      else
      {
	MPI_Send(0, 0, MPI_CHAR, source, tag_search_get, MPI_COMM_WORLD);
	//	fprintf(out, "Asking slave %d to send results\n", source);
      }
      break;

    case tag_search_report:
      
#ifdef DEBUG
      fprintf(out, "Master receiving %d hits from slave %d.\n", len/56, source);
#endif

      for (int i = 0; i < (long)(len/sizeof(long)); i+=8)
      {

	long seqno = longbuffer[i+0];
	long score = longbuffer[i+1];
	long qstrand = longbuffer[i+2];
	long qframe = longbuffer[i+3];
	long dstrand = longbuffer[i+4];
	long dframe = longbuffer[i+5];
	long pos = longbuffer[i+6];
	long bestq = longbuffer[i+7];

	//	fprintf(out, "seqno=%ld\n", seqno);
	if ((symtype == 0) && qstrand)
	  hits_enter(seqno, score, 0, 0, 1, 0, pos, bestq);
	else
	  hits_enter(seqno, score, qstrand, qframe, dstrand, dframe, pos, bestq);
      }
      
      break;

    case tag_stats:
      
      rounds7 += longbuffer[0];
      rounds16 += longbuffer[1];
      rounds32 += longbuffer[2];
      rounds63 += longbuffer[3];
      compute7 += longbuffer[4];
      compute16 += longbuffer[5];
      compute32 += longbuffer[6];
      compute63 += longbuffer[7];
      totalhits += longbuffer[8];
      
      slaves_active--;

      break;
    }
  }

  if (view == 0)
    fprintf(out, "...............................................done\n\n");

  clock_stop(&ti);

#if 0
  if (view == 0)
  {
    fprintf(out, "Computed (7bit):   %ld sequences in %ld rounds\n", compute7, rounds7);
    fprintf(out, "Computed (16bit):  %ld sequences in %ld rounds\n", compute16, rounds16);
    //    fprintf(out, "Computed (32bit):  %ld sequences in %ld rounds\n", compute32, rounds32);
    fprintf(out, "Computed (63bit):  %ld sequences in %ld rounds\n", compute63, rounds63);
    fprintf(out, "\n");
  }
#endif
  
  /* phase 2 - aligning */

  //OK master: align_threads_init();
  // slaves: align_init(&sd);
  // master: align_getwork -> slaves (i,j)
  // slaves: align_chunk(&sd, i, j)
  // slaves: align_done(&sd);
  //OK master: align_threads_done();
  

  long hitno = 0;
  long n = hits_getcount();
  long showalignments = n < alignments ? n : alignments;
  long * node_job = (long *) xmalloc(size * sizeof(long));
  long seqno, score, qstrand, qframe, dstrand, dframe;
  
  align_threads_init();

  //  if (view == 0)
  //    clock_start(&ti);

  while ((hitno < n) && (hitno < size - 1))
  {
    hits_gethit(hitno, & seqno, & score, & qstrand, & qframe, & dstrand, & dframe);
    longbuffer[0] = seqno;
    longbuffer[1] = hitno < showalignments;
    longbuffer[2] = qstrand;
    longbuffer[3] = qframe;
    longbuffer[4] = dstrand;
    longbuffer[5] = dframe;
    node_job[hitno+1] = hitno;
#ifdef DEBUG
    fprintf(out, "Asking node %ld to align hitno %ld seqno %ld.\n", hitno+1, hitno, seqno);
#endif
    MPI_Send(longbuffer, 6, MPI_LONG, hitno+1, tag_align, MPI_COMM_WORLD);
    slaves_active++;
    hitno++;
  }

  while(slaves_active)
  {
    MPI_Status status;
    int source;
    int tag;
    int len;

    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, & status);

    source = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    MPI_Get_count(& status, MPI_CHAR, & len);
    
    if (buffersize < len)
    {
      buffersize = len;
#ifdef DEBUG
      fprintf(out, "resizing buffer to length %ld.\n", buffersize);
#endif
      buffer = (char *) xrealloc(buffer, buffersize);
      longbuffer = (long*) buffer;
    }
    
    MPI_Recv(buffer, buffersize, MPI_CHAR, source, tag, MPI_COMM_WORLD, NULL);
    
#ifdef DEBUG
    fprintf(out, "Received message with tag %d and length %d from node %d.\n",
	   tag, len, source);
#endif

    switch(tag)
    {

    case tag_align_seq:

#ifdef DEBUG
      fprintf(out, "Receiving sequence from node %d\n", source);
#endif
      hits_enter_seq(node_job[source], buffer, len);
      break;

    case tag_align_coord:
      
#ifdef DEBUG
      fprintf(out, "Receiving coord from node %d\n", source);
#endif
      hits_enter_align_coord(node_job[source],
			     longbuffer[0],
			     longbuffer[1],
			     longbuffer[2],
			     longbuffer[3],
			     longbuffer[4]);
      break;

    case tag_align_string:

#ifdef DEBUG
      fprintf(out, "Receiving alignment from node %d\n", source);
#endif
      hits_enter_align_string(node_job[source], buffer, len);
      break;

    case tag_align_header:

#ifdef DEBUG
      fprintf(out, "Receiving header from node %d (hitno %ld)\n", source, node_job[source]);
#endif
      hits_enter_header(node_job[source], buffer, len);
      
      if (hitno < n)
      {
	hits_gethit(hitno, & seqno, & score, & qstrand, & qframe, & dstrand, & dframe);
	longbuffer[0] = seqno;
	longbuffer[1] = hitno < showalignments;
	longbuffer[2] = qstrand;
	longbuffer[3] = qframe;
	longbuffer[4] = dstrand;
	longbuffer[5] = dframe;
	node_job[source] = hitno;
#ifdef DEBUG
	fprintf(out, "Asking node %d to align hitno %ld seqno %ld.\n", source, hitno, seqno);
#endif
	MPI_Send(longbuffer, 6, MPI_LONG, source, tag_align, MPI_COMM_WORLD);
	hitno++;
      }
      else
	slaves_active--;
      break;

    }
  }

  //  if (view == 0)
  //    clock_stop(&ti);

#ifdef DEBUG
  fprintf(out, "Showing hits.\n");
#endif

  hits_show(view, show_gis);
  
#ifdef DEBUG
  fprintf(out, "Shown hits.\n");
#endif

  /* ask all slaves to terminate. */
  for (int i=1; i < size; i++)
    MPI_Send(0, 0, MPI_CHAR, i, tag_die, MPI_COMM_WORLD);
    
  hits_exit();

  align_threads_done();

  free(node_job);
  free(buffer);

#ifdef DEBUG
  fprintf(out, "Master completed.\n");
#endif
}

void slave(int rank, int size)
{
  (void) rank; /* avoid warning */
  
  long features = sse2_present * 2 + ssse3_present;
  MPI_Send(&features, 1, MPI_LONG, 0, tag_features, MPI_COMM_WORLD);

  //  fprintf(out, "Node %d: SSE2: %ld SSSE3: %ld\n", rank, sse2_present, ssse3_present);
  //  fprintf(out, "This is slave no %d.\n", rank);

  hits_init(maxmatches, alignments, minscore, maxscore, minexpect, expect, 0);

  prepare_search(size-1);

  long morework = 1;

  long buffersize = 1024;
  char * buffer = (char*) xmalloc(buffersize);
  long * longbuffer = (long*) buffer;
  
  //  fprintf(out, "buffersize: %ld.\n", buffersize);

  //OK master: align_threads_init();
  // slaves: align_init(&sd);
  // master: align_getwork -> slaves (i,j)
  // slaves: align_chunk(&sd, i, j)
  // slaves: align_done(&sd);
  //OK master: align_threads_done();

  struct search_data sd;

  search_init(&sd);

  totalhits = 0;
  compute7 = 0;
  compute16 = 0;
  compute32 = 0;
  compute63 = 0;
  rounds7 = 0;
  rounds16 = 0;
  rounds32 = 0;
  rounds63 = 0;

  struct db_thread_s * t = db_thread_create();

  while (morework)
  {
    MPI_Status status;
    int source;
    int tag;
    int len;

    MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, & status);

    source = status.MPI_SOURCE;
    tag = status.MPI_TAG;
    MPI_Get_count(& status, MPI_CHAR, & len);
    
    if (buffersize < len)
    {
      buffersize = len;
#ifdef DEBUG
      fprintf(out, "resizing buffer to length %ld.\n", buffersize);
#endif
      buffer = (char *) xrealloc(buffer, buffersize);
      longbuffer = (long*) buffer;
    }
    
    MPI_Recv(buffer, buffersize, MPI_CHAR, source, tag, MPI_COMM_WORLD, NULL);
    
#ifdef DEBUG
    fprintf(out, "Slave %d receiving message with tag %d from source %d of length %d.\n",
    	   rank, tag, source, len);
#endif

    long n;

    long seqno;
    long qframe;
    long qstrand;
    long dframe;
    long dstrand;
    long score;
    long do_align;

    long align_q_start;
    long align_q_end;
    long align_d_start;
    long align_d_end;

    long header_len;
    long seq_len;
    long ntlen;
    long alignment_len;
    char * header;
    char * seq;
    char * alignment;
    
    switch(tag)
    {

    case tag_search:
      
      sd.seqfirst = longbuffer[0];
      sd.seqlast = longbuffer[1];

      // fprintf(out, "Slave %d doing search from %d to %d.\n",
      //	     rank, sd.seqfirst, sd.seqlast);

      search_chunk(&sd);
      
      MPI_Send(0, 0, MPI_CHAR, 0, tag_search_done, MPI_COMM_WORLD);
      
      //      fprintf(out, "Slave %d reporting finished search.\n", rank);
      break;

    case tag_search_get:
      
      //      fprintf(out, "Slave %d sending %d results...\n", rank, hits_count);

      n = hits_getcount();
      
      if (buffersize < (long)(8*n*sizeof(long)))
      {
	buffersize = 8*n*sizeof(long);
#ifdef DEBUG
	fprintf(out, "resizing buffer to length %ld.\n", buffersize);
#endif
	buffer = (char *) xrealloc(buffer, buffersize);
	longbuffer = (long*) buffer;
      }
      
      for(long i=0; i < n; i++)
      {
	hits_gethit(i, & seqno, & score, & qstrand, & qframe, & dstrand, & dframe);

	//	fprintf(out, "slave sending seqno=%ld\n", seqno);

	longbuffer[8*i+0] = seqno;
	longbuffer[8*i+1] = score;
	longbuffer[8*i+2] = qstrand;
	longbuffer[8*i+3] = qframe;
	longbuffer[8*i+4] = dstrand;
	longbuffer[8*i+5] = dframe;
	longbuffer[8*i+6] = -1;
	longbuffer[8*i+7] = -1;
      }

      MPI_Send(longbuffer, 8*n, MPI_LONG, 0, tag_search_report, MPI_COMM_WORLD);

      longbuffer[0] = rounds7;
      longbuffer[1] = rounds16;
      longbuffer[2] = rounds32;
      longbuffer[3] = rounds63;
      longbuffer[4] = compute7;
      longbuffer[5] = compute16;
      longbuffer[6] = compute32;
      longbuffer[7] = compute63;
      longbuffer[8] = totalhits;

      MPI_Send(longbuffer, 9, MPI_LONG, 0, tag_stats, MPI_COMM_WORLD);

      break;

    case tag_align:

      seqno = longbuffer[0];
      do_align = longbuffer[1];
      qstrand = longbuffer[2];
      qframe = longbuffer[3];
      dstrand = longbuffer[4];
      dframe = longbuffer[5];

      if (do_align)
      {
	db_mapsequences(t, seqno, seqno);
	db_getsequence(t, seqno, dstrand, dframe, 
		       & seq, & seq_len, & ntlen, 0);

	char * dseq = seq;
	long dlen = seq_len - 1;

	char * qseq;
	long qlen;

	if (symtype == 0)
	{
	  qseq = query.nt[0].seq;
	  qlen = query.nt[0].len;
	}
	else
	{
	  qseq = query.aa[3*qstrand + qframe].seq;
	  qlen = query.aa[3*qstrand + qframe].len;
	}

	long align_score;

	// give no hint of alignment end
	
	align_q_end = 0;
	align_d_end = 0;
	align_score = 0;

	align(qseq,
	      dseq,
	      qlen,
	      dlen,
	      score_matrix_63,
	      gapopen,
	      gapextend,
	      & align_q_start,
	      & align_d_start,
	      & align_q_end,
	      & align_d_end,
	      & alignment,
	      & align_score);
	
	MPI_Send(dseq, dlen, MPI_CHAR, 0, tag_align_seq, MPI_COMM_WORLD);

	longbuffer[0] = align_q_start;
	longbuffer[1] = align_q_end;
	longbuffer[2] = align_d_start;
	longbuffer[3] = align_d_end;
	longbuffer[4] = ntlen;
	MPI_Send(longbuffer, 5, MPI_LONG, 0, tag_align_coord, MPI_COMM_WORLD);
	
 	alignment_len = strlen(alignment);
	MPI_Send(alignment, alignment_len+1, MPI_CHAR, 0, tag_align_string, MPI_COMM_WORLD);
	free(alignment);
      }

      db_mapheaders(t, seqno, seqno);
      db_getheader(t, seqno, & header, & header_len);
#ifdef DEBUG
      fprintf(out, "Node %d sending header of length %ld.\n", rank, header_len);
#endif
      MPI_Send(header, header_len, MPI_CHAR, 0, tag_align_header, MPI_COMM_WORLD);
      
      break;

    case tag_die:
      //      fprintf(out, "Slave %d asked to die.\n", rank);
      morework = 0;
      break;

    }
  }


  db_thread_destruct(t);
  
  //  fprintf(out, "Slave %d dies.\n", rank);

  search_done(&sd);

  free(buffer);

  hits_exit();

}

#endif

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

int main(int argc, char**argv)
{
#ifdef MPISWIPE
  int rc = MPI_Init(&argc, &argv);
  if (rc != MPI_SUCCESS)
    fatal("Unable to initialize MPI.");
  
  int rank;
  int size;
  MPI_Comm_rank(MPI_COMM_WORLD, & rank);
  MPI_Comm_size(MPI_COMM_WORLD, & size);

  mpirank = rank;
  mpisize = size;

#endif

  cpu_features();

#ifndef MPISWIPE
  if (! sse2_present)
    fatal("This program requires a processor with SSE2.\n");
#endif

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
    
#ifdef MPISWIPE
    if (!mpirank)
#endif
    {
      hits_show_begin(view);
    }
    
    while (query_read())
    {
      
#ifdef MPISWIPE
      if (rank == 0)
	master(size);
      else
	slave(rank,size);
#else
      work();
#endif
      
      queryno++;
    }
    
#ifdef MPISWIPE
    if (!mpirank)
#endif
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

#ifdef MPISWIPE
#ifdef DEBUG
  fprintf(out, "Finalizing (rank %d).\n", rank);
#endif
  MPI_Finalize();
#ifdef DEBUG
  fprintf(out, "Finalized (rank %d).\n", rank);
#endif
#endif
  
  if (outfile)
    fclose(out);
}
