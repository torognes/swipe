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
#include <array>
#include <cassert>
#include <chrono>
#include <ctime>
#include <memory>  // std::unique_ptr
#include <string>
#include <vector>
#include "fatal_allocator.h"  // Buffer, xmalloc
#include "view.h"  // View


#include "os_byteswap.h"  // bswap_32, bswap_64

// the size of the line buffers: a line of a score matrix, and the
// first allocation (then the growth step) of a query sequence
constexpr std::size_t line_buffer_size = 2048;

// "SWIPE X.Y.Z": the program name and its version (defined in swipe.cc)
extern char const * const swipe_name_and_version;

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
  sound = 5,     // sound codes
};

// output format of the results (option -m, --outfmt)
enum struct OutputFormat : long
{
  plain = 0,                  // BLAST-like plain text
  xml = 7,                    // simple XML
  tabular = 8,                // tabular (BLAST -m 8)
  tabular_with_comments = 9,  // tabular with comment lines (BLAST -m 9)
  paralign_xml = 99,           // ParAlign XML
};

// query strands to search (option -S, --strand): a bit mask of the
// plus (1) and minus (2) strands
enum struct QueryStrands : long
{
  plus = 1,
  minus = 2,
  both = 3,
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
// the gap penalties of blastn and of sound searches when not given
constexpr long default_blastn_gapopen = 5;
constexpr long default_blastn_gapextend = 2;
constexpr long default_sound_gapopen = 15;
constexpr long default_sound_gapextend = 5;
constexpr char const * default_sound_matrixname = "IDENTITY_5_1";
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

// options.cc
auto args_init(int argc, char * const * argv) -> Parameters;
auto args_show(Parameters const & parameters) -> void;



// the SIMD instruction sets of the processor, detected once (cpuid)
struct CpuFeatures
{
  bool sse2;
  bool ssse3;
};
extern CpuFeatures const cpu_features;

// the score matrices are 32 x 32, row-major (symbol codes 0 to 31)
constexpr std::size_t score_matrix_width = 32;

inline auto score_matrix_cell(std::size_t const row, std::size_t const column) -> std::size_t
{
  assert((row < score_matrix_width) and (column < score_matrix_width));
  return (row * score_matrix_width) + column;
}

// the genetic codes 1 to 23 (nullptr: no code of that number)
constexpr std::size_t gencode_count = 23;
extern std::array<char const *, gencode_count> const gencode_names;

// the tables indexed by a byte (an unsigned char)
constexpr std::size_t byte_values = 256;

// the base of the numbers read by strtol() and its siblings
constexpr int decimal_base = 10;

extern std::array<char, byte_values> const map_ncbi_nt16;
extern std::array<char, byte_values> const map_ncbi_aa;
extern std::array<char, byte_values> const map_sound;

extern char const * const sym_ncbi_nt16;
extern char const * const sym_ncbi_nt16u;
extern char const * const sym_ncbi_aa;
extern char const * const sym_sound;

// the 4-bit nucleotide codes: a bit per base (A, C, G, T), 16 values
// (0: none, 15: any base)
constexpr std::size_t nucleotide_codes = 16;

extern std::array<char, nucleotide_codes> const ntcompl;
constexpr std::size_t translation_table_size = nucleotide_codes * nucleotide_codes * nucleotide_codes;

// the codon translation tables of the genetic codes of the query (-Q)
// and of the database (-D), in query.cc
struct TranslationTables
{
  std::array<char, translation_table_size> query {{}};
  // the codon translation table of the database (16 x 16 x 16 codes of
  // nucleotides), filled by translate_init()
  std::array<char, translation_table_size> database {{}};
};

extern TranslationTables translation_tables;

extern FILE * out;


struct sequence
{
  char * seq;  // storage.data(), or nullptr
  long len;
  Buffer<char> storage;  // owns seq

  // the residues, as a read-only view
  auto view() const -> View<char>
  {
    assert(len >= 0);
    return View<char>(seq, static_cast<std::size_t>(len));
  }
};

// a nucleotide sequence has two strands, each translated in three
// reading frames
constexpr std::size_t strand_count = 2;
constexpr std::size_t frames_per_strand = 3;
constexpr std::size_t frame_count = strand_count * frames_per_strand;

// the index of a query strand (0: plus, 1: minus), and of a frame of a
// translated query or of its search tables: (3 x strand) + frame
inline auto strand_index(long const strand) -> std::size_t
{
  assert((strand >= 0) and (static_cast<std::size_t>(strand) < strand_count));
  return static_cast<std::size_t>(strand);
}

inline auto frame_index(long const strand, long const frame) -> std::size_t
{
  assert((frame >= 0) and (static_cast<std::size_t>(frame) < frames_per_strand));
  return (frames_per_strand * strand_index(strand)) + static_cast<std::size_t>(frame);
}

struct query_s
{
  std::array<struct sequence, strand_count> nt; /* 2 strands */
  std::array<struct sequence, frame_count> aa; /* 6 frames */
  std::string description;
  long dlen;
  SymbolType symtype;
  QueryStrands strands;
  char const * map;
  char const * sym;

  FILE * input;  // the query file (stdin with "-")
  // next line of the query file, with its end-of-line character (an
  // empty string means the end of the file)
  std::string line;
};

extern struct query_s query;
//extern char * qseq;
//extern long qlen;

struct db_thread_s;

// a date of the outputs ("%a, %e %b %Y %T UTC"): at most 29 characters,
// "Wed, 30 Sep 2026 07:10:01 UTC", and the terminating NUL
constexpr std::size_t date_string_size = 30;

// a speed in GCUPS: billions of cell updates per second
inline auto gcups(double const cell_updates_per_second) -> double
{
  constexpr double billion = 1e9;
  return cell_updates_per_second / billion;
}

struct time_info
{
  time_t t1, t2;
  // monotonic clock for the elapsed time (KI-33)
  std::chrono::steady_clock::time_point clock1;
  std::chrono::steady_clock::time_point clock2;

  // kept until the results are shown (-m 99, KI-28)
  std::array<char, date_string_size> starttime;
  std::array<char, date_string_size> endtime;
  double elapsed;
  double speed;
};

// the state of the run (swipe.cc): the query being searched, the
// counts of the -m 99 output, and the timing of the search
struct SearchRun
{
  long queryno = 0;  // the number of the query, from 0
  long compute7 = 0;  // sequences searched by the 7-bit stage
  long totalhits = 0;  // hits at or above the initial score threshold
  struct time_info ti;
};

extern SearchRun run;

// print the message to stderr and exit with status 1; [[noreturn]]
// belongs on the declarations: callers know that fatal() never returns
[[noreturn]] auto fatal(char const * message) noexcept -> void;
[[noreturn]] auto fatal(std::string const & message) noexcept -> void;

// the state of a thread reading the database (database.cc: its maps of
// the files, header parser and sequence buffers), owned by a DbThread
struct db_thread_s;

struct DbThreadDelete
{
  auto operator()(db_thread_s * thread) const noexcept -> void;
};

using DbThread = std::unique_ptr<db_thread_s, DbThreadDelete>;

auto search7(BYTE * * q_start,
	     BYTE gap_open_penalty,
	     BYTE gap_extend_penalty,
	     BYTE const * score_matrix,
	     BYTE * dprofile,
	     BYTE * hearray,
	     db_thread_s & dbt,
	     long sequences,
	     long const * seqnos,
	     long * scores,
	     long qlen) -> void;

auto search7_ssse3(BYTE * * q_start,
		   BYTE gap_open_penalty,
		   BYTE gap_extend_penalty,
		   BYTE const * score_matrix,
		   BYTE * dprofile,
		   BYTE * hearray,
		   db_thread_s & dbt,
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
	      db_thread_s & dbt,
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
	       DbThread const * dbta,
	       long sequences,
	       long const * seqnos,
	       long * scores,
	       long * bestpos,
	       long * bestq,
	       int qlen) -> void;

auto fullsw(char const * dseq,
	    char const * dend,
	    char const * qseq,
	    char const * qend,
	    long * hearray, 
	    long const * score_matrix,
	    long gap_open_extend,
	    long gap_extend_penalty) -> long;

auto query_init(char const * query_filename, SymbolType symbol_type, QueryStrands strands) -> void;
auto query_exit() -> void;
auto query_read() -> int;
auto query_show() -> void;

auto score_matrix_init(Parameters const & parameters) -> void;

auto translate_init(long qtableno, long dtableno) -> void;
// the reverse complement of a nucleotide sequence, NUL-terminated, into
// complement (at least sequence.size() + 1 bytes)
auto reverse_complement(View<char> sequence, char * complement) -> void;
auto revcompl(View<char> sequence) -> Buffer<char>;
// a strand (0: plus, 1: minus) and a reading frame (0 to 2, or
// untranslated_frame) of a nucleotide sequence
struct StrandFrame
{
  long strand;
  long frame;
};

// the translation of one strand and frame of a nucleotide sequence
// with a genetic code table (translation_tables.query or .database),
// NUL-terminated, into prot (at least length / 3 + 1 bytes); returns
// the protein length
auto translate_codons(View<char> sequence, StrandFrame where,
                      std::array<char, translation_table_size> const & table,
                      char * prot) -> long;

// the same for the query (-Q), into a buffer resized to fit
auto translate(View<char> sequence, StrandFrame where,
               Buffer<char> & protein) -> long;

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

auto parse_header(apt p, View<char> header, long memb, long (*f)(long),
		  HeaderLayout const & layout) -> long;

// the deflines of a header that pass the membership and taxid filters
auto parse_getdeflines(apt p, View<char> header, long memb, long (*f_checktaxid)(long), long show_gis) -> std::vector<std::string>;

auto parse_getdeflinecount(apt p, View<char> header,
                           long memb, long(*f_checktaxid)(long)) -> long;

auto db_open(Parameters const & parameters) -> void;
auto db_close() -> void;
auto db_getseqcount() -> std::int64_t;
auto db_getseqcount_masked() -> std::int64_t;
auto db_getsymcount() -> std::int64_t;
auto db_getsymcount_masked() -> std::int64_t;
auto db_getlongest() -> long;
auto db_gettitle() -> char const *;
auto db_gettime() -> char const *;
auto db_getvolumecount() -> long;
auto db_getseqcount_volume(long v) -> long;
auto db_getseqcount_volume_masked(long v) -> long;
auto db_ismasked() -> long;
auto db_getversion() -> long;

auto db_getvolume(long seqno) -> long;

auto db_thread_create() -> DbThread;

auto db_check_taxid(long taxid) -> long;

auto db_parse_header(db_thread_s const & t, View<char> header,
                     long show_gis) -> std::vector<std::string>;

auto db_showheader(db_thread_s const & t, View<char> header,
		   HeaderLayout const & layout) -> void;

auto db_show_fasta(db_thread_s & t, long seqno,
		   StrandFrame where, long split) -> void;

auto db_check_inclusion(db_thread_s & t, long seqno) -> long;

auto db_mapsequences(db_thread_s const & t, long firstseqno, long lastseqno) -> void;
auto db_mapheaders(db_thread_s const & t, long firstseqno, long lastseqno) -> void;

// frame value asking db_getsequence() for the nucleotide sequence of
// a translated database (symtypes 3 and 4), without translation
constexpr long untranslated_frame = -1;
// the channels (database sequences searched at once) of the kernels:
// 16 bytes in the 7-bit kernel, 8 words in the 16-bit kernels; c, the
// channel of db_getsequence(), is below max_channels
constexpr std::size_t channels_7 = 16;
constexpr std::size_t channels_16 = 8;
constexpr std::size_t max_channels = channels_7;
static_assert(channels_16 <= max_channels, "a buffer per channel");

// the SIMD vectors of the kernels (SSE, __m128i) are 16 bytes, and
// the buffers they load from and store to are aligned on them
constexpr std::size_t vector_bytes = 16;

// the H/E array of the kernels: per query position, H and E, one
// vector each
constexpr std::size_t hearray_row_bytes = 2 * vector_bytes;

// the score matrices of the search (matrices.cc): 32 x 32, aligned for
// the SIMD kernels, in the four score widths of the search stages, and
// the limits below which a 7-bit or 16-bit score is accepted
constexpr std::size_t score_matrix_size = score_matrix_width * score_matrix_width;

struct ScoreMatrices
{
  alignas(vector_bytes) std::array<char, score_matrix_size> score_7 {{}};
  alignas(vector_bytes) std::array<char, score_matrix_size> score_7t {{}};  // transposed
  alignas(vector_bytes) std::array<short, score_matrix_size> score_16 {{}};
  alignas(vector_bytes) std::array<long, score_matrix_size> score_63 {{}};
  long limit_7 = 0;  // SCORELIMIT_7
  long limit_16 = 0;  // SCORELIMIT_16
};

extern ScoreMatrices score_matrices;

// the residues of a sequence (GitHub #27: without the separator that
// follows it); ntlenp receives its length in nucleotides
auto db_getsequence(db_thread_s & t, long seqno, StrandFrame where,
		    long * ntlenp, std::size_t c) -> View<char>;
// the header of a sequence, as stored: binary ASN.1 (a Blast-def-line-set)
auto db_getheader(db_thread_s const & t, long seqno) -> View<char>;

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
auto hits_align(Parameters const & parameters, db_thread_s & t, long i) -> void;
auto hits_show_begin(OutputFormat view) -> void;
auto hits_show_end(OutputFormat view) -> void;
auto hits_show(Parameters const & parameters) -> void;
auto hits_empty() -> void;
auto hits_exit() -> void;
// a hit of the list: its database sequence, score, strands and frames
struct Hit
{
  long seqno;
  long score;
  HitStrands strands;
};

auto hits_gethit(long i) -> Hit;

auto hits_enter_align_hint(long i, long q_end, long d_end) -> void;


// the gap penalties of a scoring system: opening and extension
struct GapPenalties
{
  long open;
  long extend;
};

// the cells where a local alignment begins and ends (a: in the query,
// b: in the database sequence) and its score; as a hint to align(), a
// non-zero score with the end cell (the beginning is then ignored)
struct AlignmentRegion
{
  long a_begin;
  long b_begin;
  long a_end;
  long b_end;
  long score;
};

// the optimal local alignment of the two sequences (align.cc), as a
// string of operations (e.g. M12D1M5) in alignment, and its region
auto align(View<char> query_sequence,
           View<char> database_sequence,
           long const * scorematrix,
           GapPenalties gaps,
           AlignmentRegion const & hint,
           std::string & alignment) -> AlignmentRegion;

// the scores of a nucleotide match and mismatch (blastn)
struct BlastnScores
{
  long match;
  long mismatch;
};

// the Karlin-Altschul parameters of a scoring system (NCBI tables)
struct KarlinAltschul
{
  double lambda;
  double K;
  double H;
  double alpha;
  double beta;
};

// the parameters of a scoring system, if the NCBI tables have them
struct StatisticsLookup
{
  bool found;
  KarlinAltschul values;
};

// the default gap penalties of a score matrix, if it has some
struct DefaultGaps
{
  bool found;
  GapPenalties penalties;
};

auto stats_getparams_nt(BlastnScores scores, GapPenalties gaps) -> StatisticsLookup;
auto stats_getparams(char const * matrix, GapPenalties gaps) -> StatisticsLookup;
auto stats_getprefs(char const * matrix) -> DefaultGaps;


// the NCBI integer types of blastkar_partial.ccc: 4 and 8 bytes
using Int4 = std::int32_t;
using Int8 = std::int64_t;
using Nlm_FloatHi = double;

#include "blastkar_partial.h"

#endif  // SWIPE_H
