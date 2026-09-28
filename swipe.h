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

#ifndef SWIPE_H
#define SWIPE_H

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <cctype>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <getopt.h>
#include <cmath>
#include <x86intrin.h>
#include <array>
#include <chrono>
#include <ctime>
#include <string>


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

// the version number is read from the file VERSION by the Makefile
#ifndef SWIPE_VERSION
#ifdef __CPPCHECK__
// static analysis with cppcheck, run without the Makefile's flags
#define SWIPE_VERSION "0.0.0"
#else
#error "SWIPE_VERSION is not defined: build swipe with make"
#endif
#endif

// Should be 32bits integer
using UINT32 = unsigned int;
using WORD = unsigned short;
using BYTE = unsigned char;

extern char BIAS;

auto xmalloc(size_t size) -> void *;
auto xrealloc(void *ptr, size_t size) -> void *;


extern long cpu_feature_ssse3;
extern long cpu_feature_sse41;

extern char const * queryname;
extern char const * matrixname;
extern long gapopen;
extern long gapextend;
extern long gapopenextend;
extern long * score_matrix_63;
extern long symtype;
extern long matchscore;
extern long mismatchscore;
extern long totalhits;
extern char const * gencode_names[];
extern long querystrands;
extern double minexpect;
extern double expect;
extern long maxmatches;
extern long threads;
extern char const * databasename;
extern long alignments;
extern long queryno;
extern long compute7;
extern long show_taxid;
extern long effdbsize;

extern char map_ncbi_nt4[];
extern char map_ncbi_nt16[];
extern char map_ncbi_aa[];
extern char map_sound[];

extern char const * sym_ncbi_nt4;
extern char const * sym_ncbi_nt16;
extern char const * sym_ncbi_nt16u;
extern char const * sym_ncbi_aa;
extern char const * sym_sound;

extern char ntcompl[];
extern char d_translate[];

extern FILE * out;

extern char const mat_blosum45[];
extern char const mat_blosum50[];
extern char const mat_blosum62[];
extern char const mat_blosum80[];
extern char const mat_blosum90[];
extern char const mat_pam30[];
extern char const mat_pam70[];
extern char const mat_pam250[];

extern long SCORELIMIT_7;
extern long SCORELIMIT_8;
extern long SCORELIMIT_16;
extern long SCORELIMIT_32;
extern long SCORELIMIT_63;

extern char * score_matrix_7;
extern char * score_matrix_7t;
extern unsigned char * score_matrix_8;
extern short * score_matrix_16;
extern unsigned int * score_matrix_32;

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
  char const * sym;
};

extern struct query_s query;
//extern char * qseq;
//extern long qlen;

struct db_thread_s;

struct time_info
{
  time_t t1, t2;
  // monotonic clock for the elapsed time (KI-33)
  std::chrono::steady_clock::time_point clock1;
  std::chrono::steady_clock::time_point clock2;

  // kept until the results are shown (-m 99, KI-28)
  std::array<char, 30> starttime;
  std::array<char, 30> endtime;
  double elapsed;
  double speed;
};

extern struct time_info ti;

auto fatal(char const * message) -> void;
auto fatal(std::string const & message) -> void;

auto search7(BYTE * * q_start,
	     BYTE gap_open_penalty,
	     BYTE gap_extend_penalty,
	     BYTE * score_matrix,
	     BYTE * dprofile,
	     BYTE * hearray,
	     struct db_thread_s * dbt,
	     long sequences,
	     long const * seqnos,
	     long * scores,
	     long qlen) -> void;

auto search7_ssse3(BYTE * * q_start,
		   BYTE gap_open_penalty,
		   BYTE gap_extend_penalty,
		   BYTE * score_matrix,
		   BYTE * dprofile,
		   BYTE * hearray,
		   struct db_thread_s * dbt,
		   long sequences,
		   long const * seqnos,
		   long * scores,
		   long qlen) -> void;

auto search16(WORD * * q_start,
	      WORD gap_open_penalty,
	      WORD gap_extend_penalty,
	      WORD * score_matrix,
	      WORD * dprofile,
	      WORD * hearray,
	      struct db_thread_s * dbt,
	      long sequences,
	      long const * seqnos,
	      long * scores,
	      long * bestpos,
	      int qlen) -> void;

auto search16s(WORD * * q_start,
	       WORD gap_open_penalty,
	       WORD gap_extend_penalty,
	       WORD * score_matrix,
	       WORD * dprofile,
	       WORD * hearray,
	       struct db_thread_s * const * dbta,
	       long sequences,
	       long const * seqnos,
	       long * scores,
	       long * bestpos,
	       long * bestq,
	       int qlen) -> void;

auto fullsw(char * dseq,
	    char const * dend,
	    char * qseq,
	    char const * qend,
	    long * hearray, 
	    long * score_matrix,
	    long gap_open_extend,
	    long gap_extend_penalty) -> long;

auto align(char * a_seq,
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
	   long * s) -> void;

auto query_init(char const * query_filename, long symbol_type, long strands) -> void;
auto query_exit() -> void;
auto query_read() -> int;
auto query_show() -> void;

auto score_matrix_init() -> void;
auto score_matrix_free() -> void;

auto translate_init(long qtableno, long dtableno) -> void;
auto revcompl(char const * seq, long len) -> char *;
auto translate(char const * dna, long dlen,
               long strand, long frame, long table,
               char ** protp, long * plenp) -> void;

struct asnparse_info;
using apt = asnparse_info *;

auto parser_create() -> apt;
auto parser_destruct(apt p) -> void;

// XML outputs: the five special characters are escaped (KI-27)
enum struct Escaping : int { none, xml };

// print a character to out, escaped as XML (&amp; &lt; &gt; &quot; &apos;)
auto xml_putc(char symbol) noexcept -> void;

auto parse_header(apt p, unsigned char * buf, long len, long memb, long (*f)(long),
		  long show_gis, long indent, long maxlen, 
		  long linelen, long maxdeflines, long show_descr,
		  Escaping escaping = Escaping::none) -> long;

auto parse_getdeflines(apt p, unsigned char* buf, long len, long memb, long (*f_checktaxid)(long), long show_gis, long * deflines, char *** deflinetable) -> void;

auto parse_getdeflinecount(apt p, unsigned char * buf, long len,
                           long memb, long(*f_checktaxid)(long)) -> long;

auto db_open(long symbol_type, char const * basename, char * taxidfilename) -> void;
auto db_close() -> void;
auto db_getseqcount() -> long;
auto db_getseqcount_masked() -> long;
auto db_getsymcount() -> long;
auto db_getsymcount_masked() -> long;
auto db_getlongest() -> long;
auto db_gettitle() -> char*;
auto db_gettime() -> char*;
auto db_getvolumecount() -> long;
auto db_getseqcount_volume(long v) -> long;
auto db_getseqcount_volume_masked(long v) -> long;
auto db_ismasked() -> long;
auto db_getversion() -> long;

auto db_getvolume(long seqno) -> long;

auto db_thread_create() -> struct db_thread_s *;
auto db_thread_destruct(struct db_thread_s * t) -> void;

auto db_check_taxid(long taxid) -> long;

auto db_parse_header(struct db_thread_s const * t, char * address, long length,
		     long show_gis,
		     long * deflines, char *** deflinetable) -> void;

auto db_showheader(struct db_thread_s const * t, char * address, long length, 
		   long show_gis, long indent,
		   long maxlen, long linelen, long maxdeflines, long show_descr,
		   Escaping escaping = Escaping::none) -> void;
auto db_getshowheader(struct db_thread_s * t, long seqno,
		      long show_gis, long indent,
		      long maxlen, long linelen, long maxdeflines) -> void;

auto db_show_fasta(struct db_thread_s * t, long seqno,
		   long strand, long frame, long split) -> void;

auto db_check_inclusion(struct db_thread_s * t, long seqno) -> long;

auto db_mapsequences(struct db_thread_s const * t, long firstseqno, long lastseqno) -> void;
auto db_mapheaders(struct db_thread_s const * t, long firstseqno, long lastseqno) -> void;

// frame value asking db_getsequence() for the nucleotide sequence of
// a translated database (symtypes 3 and 4), without translation
constexpr long untranslated_frame = -1;
auto db_getsequence(struct db_thread_s * t, long seqno, long strand, long frame, 
		    char ** addressp, long * lengthp, long * ntlenp, int c) -> void;
auto db_getheader(struct db_thread_s const * t, long seqno, char ** address, 
		  long * length) -> void;

auto hits_init(long descriptions, long max_alignments, long minscore, 
	       long maxscore, double min_expect, double max_expect, int show_nostats) -> void;
auto hits_enter(long seqno, long score, long qstrand, long qframe,
		long dstrand, long dframe, long align_hint, long bestq) -> void;
auto hits_sort() -> long *;
auto hits_getcount() -> long;
auto hits_align(struct db_thread_s * t, long i) -> void;
auto hits_show_begin(long view) -> void;
auto hits_show_end(long view) -> void;
auto hits_show(long view, long show_gis) -> void;
auto hits_empty() -> void;
auto hits_exit() -> void;
auto hits_gethit(long i, long * seqno, long * score, 
		 long * qstrand, long * qframe,
		 long * dstrand, long * dframe) -> void;
auto hits_getfull(long i, 
		  long * seqno, 
		  long * score,
		  long * align_q_start,
		  long * align_q_end,
		  long * align_d_start,
		  long * align_d_end,
		  char ** header, long * header_len,
		  char ** seq, long * seq_len,
		  char ** align, long * align_len) -> void;
auto hits_enter_align_hint(long i, long q_end, long d_end) -> void;
auto hits_enter_header(long i, char const * header, long header_len) -> void;
auto hits_enter_seq(long hitno, char const * seq, long seq_len) -> void;
auto hits_enter_align_coord(long i,
			    long align_q_start,
			    long align_q_end,
			    long align_d_start,
			    long align_d_end,
			    long dlennt) -> void;
auto hits_enter_align_string(long hitno, char const * align, long align_len) -> void;


auto stats_getparams_nt(long match_score,
			long mismatch_score, 
			long gopen,
			long gextend,
			double * lambda,
			double * K,
			double * H,
			double * alpha,
			double * beta) -> long;

auto stats_getparams(char const * matrix,
		     long gopen,
		     long gextend,
		     double * lambda,
		     double * K,
		     double * H,
		     double * alpha,
		     double * beta) -> long;

auto stats_getprefs(char const * matrix,
		    long * gopen,
		    long * gextend) -> long;


using Int4 = int;
using Int8 = long;
using Nlm_FloatHi = double;

#include "blastkar_partial.h"

#endif  // SWIPE_H
