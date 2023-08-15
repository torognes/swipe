/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2021 Torbjorn Rognes, University of Oslo,
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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/times.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <getopt.h>
#include <math.h>
#include <x86intrin.h>

#ifdef MPISWIPE
#include <mpi.h>
#endif

#ifdef __APPLE__
#include <libkern/OSByteOrder.h>
#define bswap_32 OSSwapInt32
#define bswap_64 OSSwapInt64
#else
#include <byteswap.h>
#endif

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#define SWIPE_VERSION "2.1.1"

// Should be 32bits integer
typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

#define WIDTH 32
#define WIDTH_SHIFT 5
#define BLOCKWIDTH 32

#define ext1 ".ssq"
#define ext2 ".ssi"
#define ext3 ".shd"
#define ext4 ".shi"

extern char BIAS;


//#define BIASED

#ifdef BIASED
#define ZERO 0x00
#else
#define ZERO 0x80
#endif

void vector_print(BYTE * vector);
void vector_print_word(WORD * vector);

void * xmalloc(size_t size);
void * xrealloc(void *ptr, size_t size);


extern long cpu_feature_ssse3;
extern long cpu_feature_sse41;

extern const char * queryname;
extern const char * matrixname;
extern long gapopen;
extern long gapextend;
extern long gapopenextend;
extern long * score_matrix_63;
extern long symtype;
extern long matchscore;
extern long mismatchscore;
extern long totalhits;
extern const char * gencode_names[];
extern long querystrands;
extern double minexpect;
extern double expect;
extern long maxmatches;
extern long threads;
extern const char * databasename;
extern long alignments;
extern long queryno;
extern long compute7;
extern long show_taxid;
extern long effdbsize;

extern char map_ncbi_nt4[];
extern char map_ncbi_nt16[];
extern char map_ncbi_aa[];
extern char map_sound[];

extern const char * sym_ncbi_nt4;
extern const char * sym_ncbi_nt16;
extern const char * sym_ncbi_nt16u;
extern const char * sym_ncbi_aa;
extern const char * sym_sound;

extern char ntcompl[];
extern char d_translate[];

extern FILE * out;

extern const char mat_blosum45[];
extern const char mat_blosum50[];
extern const char mat_blosum62[];
extern const char mat_blosum80[];
extern const char mat_blosum90[];
extern const char mat_pam30[];
extern const char mat_pam70[];
extern const char mat_pam250[];

extern long SCORELIMIT_7;
extern long SCORELIMIT_8;
extern long SCORELIMIT_16;
extern long SCORELIMIT_32;
extern long SCORELIMIT_63;
extern char BIAS;

extern char * score_matrix_7;
extern char * score_matrix_7t;
extern unsigned char * score_matrix_8;
extern short * score_matrix_16;
extern unsigned int * score_matrix_32;
extern long * score_matrix_63;

struct sequence
{
  char * seq;
  long len;
};

struct query_s
{
  struct sequence nt[2]; /* 2 strands */
  struct sequence aa[6]; /* 6 frames */
  char * description;
  long dlen;
  long symtype;
  long strands;
  char * map;
  const char * sym;
};

extern struct query_s query;
//extern char * qseq;
//extern long qlen;

struct db_thread_s;

struct time_info
{
  time_t t1, t2;
  struct tms times1, times2;
  clock_t wc1, wc2;
  long clk_tck;

  char * starttime;
  char * endtime;
  double elapsed;
  double speed;
};

extern struct time_info ti;

void fatal(const char * message);
void fatal(const char * format, const char * message);

void search7(BYTE * * q_start,
	     BYTE gap_open_penalty,
	     BYTE gap_extend_penalty,
	     BYTE * score_matrix,
	     BYTE * dprofile,
	     BYTE * hearray,
	     struct db_thread_s * dbt,
	     long sequences,
	     long * seqnos,
	     long * scores,
	     long qlen);

void search7_ssse3(BYTE * * q_start,
		   BYTE gap_open_penalty,
		   BYTE gap_extend_penalty,
		   BYTE * score_matrix,
		   BYTE * dprofile,
		   BYTE * hearray,
		   struct db_thread_s * dbt,
		   long sequences,
		   long * seqnos,
		   long * scores,
		   long qlen);

void search16(WORD * * q_start,
	      WORD gap_open_penalty,
	      WORD gap_extend_penalty,
	      WORD * score_matrix,
	      WORD * dprofile,
	      WORD * hearray,
	      struct db_thread_s * dbt,
	      long sequences,
	      long * seqnos,
	      long * scores,
	      long * bestpos,
	      int qlen);

void search16s(WORD * * q_start,
	       WORD gap_open_penalty,
	       WORD gap_extend_penalty,
	       WORD * score_matrix,
	       WORD * dprofile,
	       WORD * hearray,
	       struct db_thread_s * * dbta,
	       long sequences,
	       long * seqnos,
	       long * scores,
	       long * bestpos,
	       long * bestq,
	       int qlen);

long fullsw(char * dseq,
	    char * dend,
	    char * qseq,
	    char * qend,
	    long * hearray, 
	    long * score_matrix,
	    BYTE gap_open_penalty,
	    BYTE gap_extend_penalty);

void align(char * a_seq,
	   char * b_seq,
	   long M,
	   long N,
	   long * scorematrix,
	   long q,
	   long r,
	   long * a_begin,
	   long * b_begin,
	   long * a_end,
	   long * b_end,
	   char ** alignment,
	   long * s);

void query_init(const char * queryname, long symtype, long strands);
void query_exit();
int query_read();
void query_show();

void score_matrix_init();
void score_matrix_free();

void translate_init(long qtableno, long dtableno);
char * revcompl(char * seq, long len);
void translate(char * dna, long dlen,
               long strand, long frame, long table,
               char ** protp, long * plenp);

struct asnparse_info;
typedef struct asnparse_info * apt;

apt parser_create();
void parser_destruct(apt p);

long parse_header(apt p, unsigned char * buf, long len, long memb, long (*f)(long),
		  long show_gis, long indent, long maxlen, 
		  long linelen, long maxdeflines, long show_descr);

void parse_getdeflines(apt p, unsigned char* buf, long len, long memb, long (*f_checktaxid)(long), long show_gis, long * deflines, char *** deflinetable);

long parse_getdeflinecount(apt p, unsigned char * buf, long len,
                           long memb, long(*f_checktaxid)(long));

void db_open(long symtype, const char * basename, char * taxidfilename);
void db_close();
long db_getseqcount();
long db_getseqcount_masked();
long db_getsymcount();
long db_getsymcount_masked();
long db_getlongest();
char* db_gettitle();
char* db_gettime();
long db_getvolumecount();
long db_getseqcount_volume(long v);
long db_getseqcount_volume_masked(long v);
long db_ismasked();
long db_getversion();

long db_getvolume(long seqno);

struct db_thread_s * db_thread_create();
void db_thread_destruct(struct db_thread_s * t);

long db_check_taxid(long taxid);

void db_parse_header(struct db_thread_s * t, char * address, long length,
		     long show_gis,
		     long * deflines, char *** deflinetable);

void db_showheader(struct db_thread_s * t, char * address, long length, 
		   long show_gis, long indent,
		   long maxlen, long linelen, long maxdeflines, long show_descr);
void db_getshowheader(struct db_thread_s * t, long seqno,
		      long show_gis, long indent,
		      long maxlen, long linelen, long maxdeflines);

void db_show_fasta(struct db_thread_s * t, long seqno,
		   long strand, long frame, long split);

long db_check_inclusion(struct db_thread_s * t, long seqno);

void db_mapsequences(struct db_thread_s * t, long firstseqno, long lastseqno);
void db_mapheaders(struct db_thread_s * t, long firstseqno, long lastseqno);

void db_getsequence(struct db_thread_s * t, long seqno, long strand, long frame, 
		    char ** addressp, long * lengthp, long * ntlenp, int c);
void db_getheader(struct db_thread_s * t, long seqno, char ** address, 
		  long * length);

void hits_init(long descriptions, long alignments, long minscore, 
	       long maxscore, double minexpect, double expect, int show_nostats);
void hits_enter(long seqno, long score, long qstrand, long qframe,
		long dstrand, long dframe, long align_hint, long bestq);
long * hits_sort();
long hits_getcount();
void hits_align(struct db_thread_s * t, long i);
void hits_show_begin(long view);
void hits_show_end(long view);
void hits_show(long view, long show_gis, long show_best);
void hits_empty();
void hits_exit();
void hits_gethit(long i, long * seqno, long * score, 
		 long * qstrand, long * qframe,
		 long * dstrand, long * dframe);
void hits_getfull(long i, 
		  long * seqno, 
		  long * score,
		  long * align_q_start,
		  long * align_q_end,
		  long * align_d_start,
		  long * align_d_end,
		  char ** header, long * header_len,
		  char ** seq, long * seq_len,
		  char ** align, long * align_len);
void hits_enter_align_hint(long i, long q_end, long d_end);
void hits_enter_header(long i, char * header, long header_len);
void hits_enter_seq(long hitno, char* buffer, long len);
void hits_enter_align_coord(long i,
			    long align_q_start,
			    long align_q_end,
			    long align_d_start,
			    long align_d_end,
			    long dlennt);
void hits_enter_align_string(long hitno, char * align, long align_len);


long stats_getparams_nt(long matchscore,
			long mismatchscore, 
			long gopen,
			long gextend,
			double * lambda,
			double * K,
			double * H,
			double * alpha,
			double * beta);

long stats_getparams(const char * matrix,
		     long gopen,
		     long gextend,
		     double * lambda,
		     double * K,
		     double * H,
		     double * alpha,
		     double * beta);

long stats_getprefs(const char * matrix,
		    long * gopen,
		    long * gextend);


typedef int Int4;
typedef long Int8;
typedef double Nlm_FloatHi;

#include "blastkar_partial.h"
