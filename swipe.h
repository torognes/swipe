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
#include <cstdint>  // std::int32_t, std::int64_t
#include <cstdlib>
#include <climits>
#include <cctype>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include <arpa/inet.h>
#include <getopt.h>
#include <cmath>
#include <x86intrin.h>
#include <array>
#include <chrono>
#include <ctime>
#include <string>
#include "fatal_allocator.h"  // Buffer, xmalloc


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

// "SWIPE X.Y.Z": the program name and its version (defined in swipe.cc)
extern char const swipe_name_and_version[];

// Should be 32bits integer
using UINT32 = unsigned int;
using WORD = unsigned short;
using BYTE = unsigned char;

// symbol type of the search (option -p, --symtype): the query and
// database sequence types, and their translations
enum struct SymbolType : long
{
  blastn = 0,   // nucleotide query, nucleotide database
  blastp = 1,   // amino acid query, amino acid database
  blastx = 2,   // translated nucleotide query, amino acid database
  tblastn = 3,  // amino acid query, translated nucleotide database
  tblastx = 4,  // translated query, translated database
  sound = 5     // sound codes
};

// output format of the results (option -m, --outfmt)
enum struct OutputFormat : long
{
  plain = 0,                  // BLAST-like plain text
  xml = 7,                    // simple XML
  tabular = 8,                // tabular (BLAST -m 8)
  tabular_with_comments = 9,  // tabular with comment lines (BLAST -m 9)
  paralign_xml = 99           // ParAlign XML
};

// query strands to search (option -S, --strand): a bit mask of the
// plus (1) and minus (2) strands
enum struct QueryStrands : long
{
  plus = 1,
  minus = 2,
  both = 3
};

// true when the query strand of index strand (0: plus, 1: minus) is
// searched
inline auto searches_strand(QueryStrands const strands, long const strand) -> bool
{
  return ((strand + 1) & static_cast<long>(strands)) != 0;
}

/* ARGUMENTS AND THEIR DEFAULTS */

constexpr long default_maxmatches = 250;
constexpr long default_alignments = 100;
constexpr long default_minscore = 1;
constexpr long default_maxscore = LONG_MAX;
constexpr char const * default_queryname = "-";
constexpr char const * default_databasename = "";
constexpr long default_gapopen = 0;
constexpr long default_gapextend = 0;
constexpr char const * default_matrixname = "BLOSUM62";
constexpr long default_matchscore = 1;
constexpr long default_mismatchscore = -3;
constexpr long default_threads = 1;
constexpr OutputFormat default_view = OutputFormat::plain;
constexpr SymbolType default_symtype = SymbolType::blastp;
constexpr long default_show_gis = 0;
constexpr long default_show_taxid = 0;
constexpr double default_expect = 10.0;
constexpr double default_minexpect = 0.0;
constexpr QueryStrands default_querystrands = QueryStrands::both;
constexpr long default_query_gencode = 1;
constexpr long default_db_gencode = 1;
constexpr long default_subalignments = 1;
constexpr long default_dump = 0;
constexpr std::int64_t default_effdbsize = 0;

// the command-line options, final once args_init() has parsed and
// checked them (including the gap penalties that default to those of
// the score matrix or of the symbol type); gapopenextend is derived
struct Parameters
{
  char * progname = nullptr;
  char const * matrixname = "";
  char const * databasename = default_databasename;
  char const * queryname = default_queryname;
  char * taxidfilename = nullptr;
  char * outfile = nullptr;
  double expect = default_expect;
  double minexpect = default_minexpect;
  long alignments = default_alignments;
  long maxmatches = default_maxmatches;
  long minscore = default_minscore;
  long maxscore = default_maxscore;
  long gapopen = default_gapopen;
  long gapextend = default_gapextend;
  long gapopenextend = 0;
  long matchscore = default_matchscore;
  long mismatchscore = default_mismatchscore;
  long threads = default_threads;
  SymbolType symtype = default_symtype;
  QueryStrands querystrands = default_querystrands;
  OutputFormat view = default_view;
  long show_gis = default_show_gis;
  long show_taxid = default_show_taxid;
  long query_gencode = default_query_gencode;
  long db_gencode = default_db_gencode;
  long subalignments = default_subalignments;
  long dump = default_dump;
  std::int64_t effdbsize = default_effdbsize;
};

auto xrealloc(void *ptr, size_t size) -> void *;


extern long cpu_feature_ssse3;
extern long cpu_feature_sse41;

extern long * const score_matrix_63;
extern long totalhits;
extern char const * gencode_names[];
extern long queryno;
extern long compute7;

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
extern long SCORELIMIT_16;

extern char * const score_matrix_7;
extern char * const score_matrix_7t;
extern short * const score_matrix_16;

struct sequence
{
  char * seq;  // storage.data(), or nullptr
  long len;
  Buffer<char> storage;  // owns seq
};

struct query_s
{
  struct sequence nt[2]; /* 2 strands */
  struct sequence aa[6]; /* 6 frames */
  std::string description;
  long dlen;
  SymbolType symtype;
  QueryStrands strands;
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

// print the message to stderr and exit with status 1; [[noreturn]]
// belongs on the declarations: callers know that fatal() never returns
[[noreturn]] auto fatal(char const * message) noexcept -> void;
[[noreturn]] auto fatal(std::string const & message) noexcept -> void;

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

auto query_init(char const * query_filename, SymbolType symbol_type, QueryStrands strands) -> void;
auto query_exit() -> void;
auto query_read() -> int;
auto query_show() -> void;

auto score_matrix_init(Parameters const & parameters) -> void;

auto translate_init(long qtableno, long dtableno) -> void;
auto revcompl(char const * seq, long len) -> Buffer<char>;
auto translate(char const * dna, long dlen,
               long strand, long frame, long table,
               Buffer<char> & protein, long * plenp) -> void;

struct asnparse_info;
using apt = asnparse_info *;

auto parser_create(long show_taxid) -> apt;
auto parser_destruct(apt p) -> void;

// XML outputs: the five special characters are escaped (KI-27)
enum struct Escaping : int { none, xml };

// print a character to out, escaped as XML (&amp; &lt; &gt; &quot; &apos;)
auto xml_putc(char symbol) noexcept -> void;

// deflines: whole, or only their identifier (up to the first space)
enum struct DeflineText : int { identifier, full };

// how parse_header() and db_showheader() print the deflines of a
// database sequence (the defaults: the first defline, whole, on one
// line, neither truncated nor padded)
struct HeaderLayout
{
  long show_gis = 0;  // non-zero: show the gi numbers (-I)
  long indent = 0;  // continuation lines, when maxdeflines > 1
  long maxlen = 0;  // truncated after maxlen characters (0: never)
  long linelen = LONG_MAX;  // wrapped and padded to linelen (LONG_MAX: never)
  long maxdeflines = 1;  // more than one: one defline per line
  DeflineText text = DeflineText::full;
  Escaping escaping = Escaping::none;
};

auto parse_header(apt p, unsigned char * buf, long len, long memb, long (*f)(long),
		  HeaderLayout const & layout) -> long;

auto parse_getdeflines(apt p, unsigned char* buf, long len, long memb, long (*f_checktaxid)(long), long show_gis, long * deflines, char *** deflinetable) -> void;

auto parse_getdeflinecount(apt p, unsigned char * buf, long len,
                           long memb, long(*f_checktaxid)(long)) -> long;

auto db_open(Parameters const & parameters) -> void;
auto db_close() -> void;
auto db_getseqcount() -> std::int64_t;
auto db_getseqcount_masked() -> std::int64_t;
auto db_getsymcount() -> std::int64_t;
auto db_getsymcount_masked() -> std::int64_t;
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
		   HeaderLayout const & layout) -> void;

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

auto hits_init(Parameters const & parameters) -> void;
// strands and frames of a hit: query and database sequence
struct HitStrands
{
  long qstrand;
  long qframe;
  long dstrand;
  long dframe;
};

auto hits_enter(long seqno, long score, HitStrands const & strands) -> void;
auto hits_sort() -> Buffer<long>;
auto hits_getcount() -> long;
auto hits_align(Parameters const & parameters, struct db_thread_s * t, long i) -> void;
auto hits_show_begin(OutputFormat view) -> void;
auto hits_show_end(OutputFormat view) -> void;
auto hits_show(Parameters const & parameters) -> void;
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


// the NCBI integer types of blastkar_partial.c: 4 and 8 bytes
using Int4 = std::int32_t;
using Int8 = std::int64_t;
using Nlm_FloatHi = double;

#include "blastkar_partial.h"

#endif  // SWIPE_H
