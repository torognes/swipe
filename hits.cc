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
#include "decimal_digits.h"  // decimal::Buffer, decimal::to_decimal
#include "print_view.h"  // as_c_string, fprint, fprint_integer, fprint_spaces
#include <algorithm>  // std::find_if, std::min, std::move_backward, std::sort
#include <array>
#include <cassert>
#include <cctype>  // std::isspace
#include <cmath>  // std::isnan
#include <cstddef>  // std::size_t
#include <cstdint>  // std::int64_t, INT64_C
#include <cstdlib>  // std::strtol
#include <cstring>  // std::strncmp
#include <initializer_list>
#include <iterator>  // std::distance, std::next
#include <limits>
#include <mutex>  // std::mutex, std::lock_guard
#include <numeric>  // std::iota
#include <string>
#include <tuple>  // std::make_tuple
#include <utility>  // std::move
#include <vector>

/* parameters for bit scores and expect values */

namespace {

// the statistical parameters of the search (Karlin-Altschul), set by
// hits_init()
struct Statistics
{
  long available = 0;

  double alpha = 0;
  double beta = 0;
  double lambda = 0;
  double K = 0;
  double H = 0;
  double Kmn = 0;

  /* ungapped statistical parameters, only shown with -m 99 (KI-31) */
  double ungapped_lambda = 0;
  double ungapped_K = 0;
  double ungapped_H = 0;

  double logK = 0;
  double lambda_d_log2 = 0;
  double logK_d_log2 = 0;
};

Statistics statistics;

// the parameters of a lookup, if found: 1 (statistics available), or 0
auto take_statistics(StatisticsLookup const & lookup) -> long
{
  if (not lookup.found)
  {
    return 0;
  }
  statistics.lambda = lookup.values.lambda;
  statistics.K = lookup.values.K;
  statistics.H = lookup.values.H;
  statistics.alpha = lookup.values.alpha;
  statistics.beta = lookup.values.beta;
  return 1;
}

}  // anonymous namespace

/* gap penalties of the ungapped rows of the NCBI score matrix tables
   (INT2_MAX, see blastkar_partial.cc) */
constexpr long ungapped_penalty = 32767;

// the natural logarithm of 2: bit scores are in base 2
constexpr double ln_2 = 0.693147180559945309417;

namespace {

// E-value and bit score of a raw score, and a percentage (the same
// expressions as at their former call sites: identical results)
auto expect_value_of(long const score) -> double
{
  return statistics.Kmn * exp(- statistics.lambda * static_cast<double>(score));
}

auto bit_score_of(long const score) -> double
{
  return (statistics.lambda_d_log2 * static_cast<double>(score)) - statistics.logK_d_log2;
}

// the most hits of a database sequence: one per searched pair of a
// query strand or frame and a database frame
auto hits_per_sequence(Parameters const & parameters) -> std::int64_t
{
  auto const strands = static_cast<std::int64_t>(strand_count);
  auto const frames = static_cast<std::int64_t>(frames_per_strand);
  auto const all_frames = static_cast<std::int64_t>(frame_count);
  auto const query_strands =
    (parameters.querystrands == QueryStrands::both) ? strands : 1;

  if (parameters.symtype == SymbolType::blastn)
  {
    return query_strands;
  }
  if (parameters.symtype == SymbolType::blastx)
  {
    return query_strands * frames;
  }
  if (parameters.symtype == SymbolType::tblastn)
  {
    return all_frames;
  }
  if (parameters.symtype == SymbolType::tblastx)
  {
    return query_strands * frames * all_frames;
  }
  return 1;  // blastp, sound
}

// the score of an aligned pair of residues (symbol codes)
auto pair_score(char const query_symbol, char const db_symbol) -> long
{
  return score_matrices.score_63[score_matrix_cell(static_cast<unsigned char>(query_symbol),
                                           static_cast<unsigned char>(db_symbol))];
}

// the plain output: the list of hits has a description column, then
// the strand or frames of the hit, then the score column(s); the
// header of an alignment is indented and wrapped
constexpr long description_width = 67;
constexpr std::size_t score_width = 5;
constexpr long alignment_header_indent = 10;
constexpr long alignment_header_width = 79;

// the width of the strand or frames of a hit, after its description:
// " +" (blastn), " +1" (blastx, tblastn), " +1/-2" (tblastx)
auto frame_mark_width(SymbolType const symbol_type) -> long
{
  constexpr long strand_mark = 2;
  constexpr long frame_mark = 3;
  constexpr long frame_pair_mark = 6;
  if (symbol_type == SymbolType::blastn)
  {
    return strand_mark;
  }
  if ((symbol_type == SymbolType::blastx) or (symbol_type == SymbolType::tblastn))
  {
    return frame_mark;
  }
  if (symbol_type == SymbolType::tblastx)
  {
    return frame_pair_mark;
  }
  return 0;
}

// a whole percentage, rounded down, as in "Identities = 9/12 (75%)"
auto whole_percentage(long const part, long const whole) -> long
{
  constexpr long hundred = 100;
  return part * hundred / whole;
}

auto percentage(long const part, long const whole) -> double
{
  return 100.0 * static_cast<double>(part) / static_cast<double>(whole);
}

// effective length of the database in the search space: its length
// minus the length adjustment of each of its sequences, computed in
// 64 bits (KI-40: an int product overflowed for large databases)
constexpr auto effective_db_length(std::int64_t const db_length,
                                   std::int64_t const sequence_count,
                                   std::int64_t const length_adjustment) -> std::int64_t
{
  return db_length - (sequence_count * length_adjustment);
}

static_assert(effective_db_length(1000, 10, 3) == 970, "search space");
// 50 million sequences, adjustment of 100: the product (5e9) does not
// fit in an int (KI-40)
static_assert(effective_db_length(INT64_C(20000000000), 50000000, 100) == INT64_C(15000000000),
              "search space of a large database (KI-40)");

struct hits_entry
{
  std::string alignment;
  Buffer<char> dseq;
  Buffer<char> header_address;
  long seqno;
  long qstrand;
  long qframe;
  long dstrand;
  long dframe;
  long dlen;
  long dlennt;
  long score;
  long score_align;
  long align_hint;
  long bestq;
  long align_q_start;
  long align_q_end;
  long align_d_start;
  long align_d_end;
};

// the best hits of the query, sorted by decreasing score, entered by
// the search threads (hits_enter(), under the mutex)
struct HitList
{
  Buffer<hits_entry> entries;
  int count = 0;
  long keep = 0;  // the size of the list: the most descriptions or alignments
  long score_threshold = 0;  // a hit below it is not kept
  long upper_score_threshold = 0;
  long init_threshold = 0;
  long obvious = 0;  // the hits above upper_score_threshold
  long descriptions = 0;  // -v
  long alignments = 0;  // -b
  std::mutex mutex;
};

HitList hit_list;

// the hit of rank i in the list (a long, as the hit counts)
auto hit_entry(long const i) -> struct hits_entry &
{
  assert((i >= 0) and (static_cast<std::size_t>(i) < hit_list.entries.size()));
  return hit_list.entries[static_cast<std::size_t>(i)];
}


// the order of the hits for the alignment step: the hits to align
// first, then by query strand and frame, sequence, database strand and
// frame
auto hits_less(long const lhs, long const rhs) -> bool
{
  auto const & a = hit_entry(lhs);
  auto const & b = hit_entry(rhs);
  bool const a_unaligned = lhs >= hit_list.alignments;
  bool const b_unaligned = rhs >= hit_list.alignments;
  return std::make_tuple(a_unaligned, a.qstrand, a.qframe, a.seqno, a.dstrand, a.dframe) <
    std::make_tuple(b_unaligned, b.qstrand, b.qframe, b.seqno, b.dstrand, b.dframe);
}

}  // anonymous namespace

auto hits_sort() -> Buffer<long>
{
  Buffer<long> hits_sorted(static_cast<std::size_t>(hit_list.count));
  std::iota(hits_sorted.begin(), hits_sorted.end(), 0L);
  std::sort(hits_sorted.begin(), hits_sorted.end(), hits_less);
  return hits_sorted;
}


auto hits_enter(long seqno, long score, HitStrands const & strands) -> void
{
  // show_progress();

  //  fprintf(out, "Entering score %u from sequence no %u.\n", score, seqno);
  
  // find correct place

  std::lock_guard<std::mutex> const lock(hit_list.mutex);

  if (score > hit_list.upper_score_threshold)
  {
    hit_list.obvious++;
  }

  if (score >= hit_list.init_threshold)
  {
    run.totalhits++;
  }

  if ((score < hit_list.score_threshold) || (score > hit_list.upper_score_threshold))
  {
    return;
  }

  long place = hit_list.count;

  while ((place > 0) && ((score > hit_entry(place-1).score) ||
			 ((score == hit_entry(place-1).score) &&
			  (seqno > hit_entry(place-1).seqno))))
  {
    place--;
  }

  // move entries down
  
  long const move = (hit_list.count < hit_list.keep ? hit_list.count : hit_list.keep - 1) - place;

  //  fprintf(out, "Inserting at place %d, moving %d.\n", place, move);

  if (move > 0)
  {
    // entries place to place + move - 1 shift to place + 1 to place + move
    assert(static_cast<std::size_t>(place + move) < hit_list.entries.size());
    auto const first = std::next(hit_list.entries.begin(), place);
    auto const last = std::next(first, move);
    std::move_backward(first, last, std::next(last));
  }

  // fill new entry

  if (place < hit_list.keep)
  {
    hit_entry(place).seqno = seqno;
    hit_entry(place).qstrand = strands.qstrand;
    hit_entry(place).qframe = strands.qframe;
    hit_entry(place).dstrand = strands.dstrand;
    hit_entry(place).dframe = strands.dframe;
    hit_entry(place).score = score;
    // set by hits_enter_align_hint(), for the hits to align
    hit_entry(place).align_hint = -1;
    hit_entry(place).bestq = -1;
    if (hit_list.count < hit_list.keep)
    {
      hit_list.count++;
    }
  }
  
  // no hit is kept with -v 0 -b 0: the list is empty (KI-10)
  if ((hit_list.keep > 0) and (hit_list.count == hit_list.keep))
  {
    hit_list.score_threshold = hit_entry(hit_list.keep - 1).score;
  }

}

auto hits_getcount() -> long
{
  return hit_list.count;
}

auto hits_gethit(long i) -> Hit
{
  auto const & h = hit_entry(i);
  return {h.seqno, h.score, {h.qstrand, h.qframe, h.dstrand, h.dframe}};
}

auto hits_enter_align_hint(long i, long q_end, long d_end) -> void
{
  hit_entry(i).bestq = q_end;
  hit_entry(i).align_hint = d_end;
}

namespace {

// score thresholds computed from E-values can be infinite (e.g. an
// empty database, Kmn = 0) or beyond the range of long: converting
// them with a cast is undefined behaviour (KI-7)
auto threshold_to_long(double const value) -> long
{
  assert(not std::isnan(value));
  constexpr auto upper_limit = static_cast<double>(std::numeric_limits<long>::max());
  constexpr auto lower_limit = static_cast<double>(std::numeric_limits<long>::min());
  if (value >= upper_limit)
  {
    return std::numeric_limits<long>::max();
  }
  if (value <= lower_limit)
  {
    return std::numeric_limits<long>::min();
  }
  return static_cast<long>(value);
}

}  // anonymous namespace

auto hits_init(Parameters const & parameters) -> void
{
  long const descriptions = parameters.maxmatches;
  long const max_alignments = parameters.alignments;
  long const minscore = parameters.minscore;
  long const maxscore = parameters.maxscore;
  double const min_expect = parameters.minexpect;
  double const max_expect = parameters.expect;
  auto const show_nostats = static_cast<int>(parameters.view == OutputFormat::plain);

  hit_list.descriptions = descriptions;
  hit_list.alignments = max_alignments;
  hit_list.keep = descriptions > max_alignments ? descriptions : max_alignments;
  
  auto const maxhits = db_getseqcount_masked() * hits_per_sequence(parameters);

  hit_list.keep = static_cast<long>(std::min<std::int64_t>(hit_list.keep, maxhits));

  hit_list.obvious = 0;
  hit_list.count = 0;
  hit_list.entries.clear();
  hit_list.entries.resize(static_cast<std::size_t>(hit_list.keep));

  std::int64_t seqcount = 0;
  std::int64_t symcount = 0;

  if (db_ismasked() != 0)
  {
    seqcount = db_getseqcount_masked();
    symcount = db_getsymcount_masked();
  }
  else
  {
    seqcount = db_getseqcount();
    symcount = db_getsymcount();
  }

  //fprintf(out, "matrix=%s, go=%ld, ge=%ld\n", matrixname, gapopen, gapextend);

  statistics.available = 0;

  if (parameters.symtype == SymbolType::blastn)
  {
    statistics.available = take_statistics(stats_getparams_nt({parameters.matchscore, parameters.mismatchscore},
                                                              {parameters.gapopen, parameters.gapextend}));
    if (statistics.available != 0)
    {

      /*
      fprintf(out, "Params: lambda=%6.3g K=%6.3g H=%6.3g alpha=%6.3g beta=%6.3g\n",
	      lambda, K, H, alpha, beta);
      */

      statistics.logK = log(statistics.K);
      statistics.lambda_d_log2 = statistics.lambda / ln_2;
      statistics.logK_d_log2 = statistics.logK / ln_2;
      
      long const qlen = query.nt[0].len;

      std::int64_t dlen = 0;
      if (parameters.effdbsize > 0)
      {
	dlen = parameters.effdbsize;
      }
      else
      {
	dlen = symcount;
      }

      int const lenadj = length_adjustment(statistics.K, statistics.logK, statistics.alpha / statistics.lambda, statistics.beta,
					   qlen, dlen, seqcount);
    
      //      fprintf(out, "lenadj: %d\n", lenadj);

      std::int64_t const m = qlen - lenadj;

      std::int64_t const n = (parameters.effdbsize > 0)
	? parameters.effdbsize
	: effective_db_length(dlen, seqcount, lenadj);

      statistics.Kmn = statistics.K * static_cast<double>(m) * static_cast<double>(n);
    }
  }
  else if (parameters.symtype < SymbolType::sound)
  {
    if (parameters.symtype == SymbolType::tblastx)
    {
      statistics.available = take_statistics(stats_getparams(parameters.matrixname,
                                                             {ungapped_penalty, ungapped_penalty}));
    }
    else
    {
      statistics.available = take_statistics(stats_getparams(parameters.matrixname,
                                                             {parameters.gapopen, parameters.gapextend}));
    }


    if (statistics.available != 0)
    {
      
      statistics.logK = log(statistics.K);
      statistics.lambda_d_log2 = statistics.lambda / ln_2;
      statistics.logK_d_log2 = statistics.logK / ln_2;
      
      long qlen = query.aa[0].len;
      if ((parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastx))
      {
	qlen = query.nt[0].len / 3;
      }

      std::int64_t dlen = 0;
      if (parameters.effdbsize > 0)
      {
	dlen = parameters.effdbsize;
      }
      else
      {
	if ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx))
	{
	  dlen = symcount / 3;
	}
	else
	{
	  dlen = symcount;
	}
      }

      int const lenadj = length_adjustment(statistics.K, statistics.logK, statistics.alpha / statistics.lambda, statistics.beta,
					   qlen, dlen, seqcount);

      std::int64_t const m = qlen - lenadj;

      std::int64_t const n = (parameters.effdbsize > 0)
	? parameters.effdbsize
	: effective_db_length(dlen, seqcount, lenadj);

      statistics.Kmn = statistics.K * static_cast<double>(m) * static_cast<double>(n);
    }
  }

  /* ungapped statistical parameters (-m 99): the (0, 0) rows of the
     nucleotide tables, the ungapped rows of the matrix tables; the
     gapped values when there are none (KI-31) */
  statistics.ungapped_lambda = statistics.lambda;
  statistics.ungapped_K = statistics.K;
  statistics.ungapped_H = statistics.H;
  StatisticsLookup ungapped {false, {0, 0, 0, 0, 0}};
  if (parameters.symtype == SymbolType::blastn)
  {
    ungapped = stats_getparams_nt({parameters.matchscore, parameters.mismatchscore}, {0, 0});
  }
  else if (parameters.symtype < SymbolType::sound)
  {
    ungapped = stats_getparams(parameters.matrixname, {ungapped_penalty, ungapped_penalty});
  }
  if (ungapped.found)
  {
    statistics.ungapped_lambda = ungapped.values.lambda;
    statistics.ungapped_K = ungapped.values.K;
    statistics.ungapped_H = ungapped.values.H;
  }

  hit_list.score_threshold = minscore;
  hit_list.upper_score_threshold = maxscore;
  
  if (statistics.available != 0)
  {
    auto const minscore_expect = threshold_to_long(ceil(- log(max_expect / statistics.Kmn) / statistics.lambda));
    if (minscore_expect > minscore)
    {
      hit_list.score_threshold = minscore_expect;
    }

    if (min_expect > 0.0)
    {
      auto const maxscore_expect = threshold_to_long(floor(- log(min_expect / statistics.Kmn) / statistics.lambda));
      if (maxscore_expect < maxscore)
      {
	hit_list.upper_score_threshold = maxscore_expect;
      }
    }
  }
  else
  {
    if (show_nostats != 0)
    {
      fprint(out, "Statistical parameters are not available for the scoring system specified.\nBit scores and E-values will not be computed.\n\n");
    }
  }

  hit_list.init_threshold = hit_list.score_threshold;

  //  fprintf(out, "scorethreshold: %ld\n", scorethreshold);
}

auto hits_empty() -> void
{
  for (long i=0; i<hit_list.count; i++)
  {
    struct hits_entry * h = &hit_entry(i);

    h->header_address = Buffer<char>();
    h->dseq = Buffer<char>();
    h->alignment = std::string();
  }
}

auto hits_exit() -> void
{
  hits_empty();
  hit_list.entries = Buffer<hits_entry>();
}

auto hits_align(Parameters const & parameters, db_thread_s & t, long i) -> void
{
  long ntlen = 0;

  struct hits_entry * h = &hit_entry(i);

  db_mapheaders(t, h->seqno, h->seqno);

  auto const header = db_getheader(t, h->seqno);
  h->header_address.assign(header.begin(), header.end());

  // the sequence length is needed for every hit shown (-m 7 <len>,
  // KI-37), the sequence itself only for hits with an alignment
  db_mapsequences(t, h->seqno, h->seqno);

  View<char> const sequence = db_getsequence(t, h->seqno, {h->dstrand, h->dframe},
					     & ntlen, 0);
  h->dlen = static_cast<long>(sequence.size());
  h->dlennt = ntlen;

  if (i < hit_list.alignments)
  {
    h->dseq.assign(sequence.begin(), sequence.end());
    
    auto const & query_sequence = (parameters.symtype == SymbolType::blastn) ?
      query.nt[0] : query.aa[frame_index(h->qstrand, h->qframe)];

    // give hint of alignment end

    if ((h->bestq > 0) && (h->align_hint != 0))
    {
      h->score_align = h->score;
      h->align_q_end = h->bestq;
      h->align_d_end = h->align_hint;
    }
    else
    {
      h->score_align = 0;
      h->align_q_end = 0;
      h->align_d_end = 0;
    }

    auto const region = align(query_sequence.view(),
			      make_view(h->dseq),
			      score_matrices.score_63.data(),
			      {parameters.gapopen, parameters.gapextend},
			      {0, 0, h->align_q_end, h->align_d_end, h->score_align},
			      h->alignment);
    h->align_q_start = region.a_begin;
    h->align_d_start = region.b_begin;
    h->align_q_end = region.a_end;
    h->align_d_end = region.b_end;
    h->score_align = region.score;
  }
}


constexpr std::size_t ALIGNLEN = 60;

namespace {

// a hit, as seen by the alignment printers: its sequences, strands
// and frames, and the first and last positions of its alignment, as
// displayed (1-based, in nucleotides for translated sequences)
struct AlignedHit
{
  SymbolType symtype = default_symtype;
  char const * alignment = nullptr;
  char const * sym = nullptr;
  char const * q_seq = nullptr;
  char const * d_seq = nullptr;
  long q_align_start = 0;
  long d_align_start = 0;
  long q_len = 0;
  long q_len_nt = 0;
  long d_len = 0;
  long d_len_nt = 0;
  long q_strand = 0;
  long q_frame = 0;
  long d_strand = 0;
  long d_frame = 0;
  long q_first = 0;
  long q_last = 0;
  long d_first = 0;
  long d_last = 0;
  int poswidth = 1;
};

auto aligned_hit(Parameters const & parameters, long const i) -> AlignedHit
{
  auto const & entry = hit_entry(i);
  AlignedHit hit;
  hit.symtype = parameters.symtype;
  hit.alignment = entry.alignment.c_str();
  hit.q_align_start = entry.align_q_start;
  hit.d_align_start = entry.align_d_start;
  hit.q_strand = entry.qstrand;
  hit.q_frame = entry.qframe;
  hit.d_strand = entry.dstrand;
  hit.d_frame = entry.dframe;
  
  if (parameters.symtype == SymbolType::blastn)
  {
    // hits of the reverse complement of the query are entered on the
    // reverse strand of the database sequence (reported_strands())
    assert(hit.q_strand == 0);
    hit.sym = sym_ncbi_nt16;
    hit.q_seq = query.nt[strand_index(hit.q_strand)].seq;
    hit.q_len = query.nt[strand_index(hit.q_strand)].len;
  }
  else if (parameters.symtype == SymbolType::sound)
  {
    hit.sym = sym_sound;
    hit.q_seq = query.aa[0].seq;
    hit.q_len = query.aa[0].len;
  }
  else
  {
    hit.sym = sym_ncbi_aa;
    hit.q_seq = query.aa[frame_index(hit.q_strand, hit.q_frame)].seq;
    hit.q_len = query.aa[frame_index(hit.q_strand, hit.q_frame)].len;
    hit.q_len_nt = query.nt[0].len;
    hit.d_len_nt = entry.dlennt;
  }

  hit.d_seq = entry.dseq.data();
  hit.d_len = entry.dlen;

  /* calculate first and last alignment positions for display */

  hit.q_first = entry.align_q_start;
  hit.q_last = entry.align_q_end;
  hit.d_first = entry.align_d_start;
  hit.d_last = entry.align_d_end;
  
  if (parameters.symtype == SymbolType::blastn)
  {
    if (hit.q_strand != 0)
    {
      hit.q_first = hit.q_len - 1 - hit.q_first;
      hit.q_last = hit.q_len - 1 - hit.q_last;
    }

    if (hit.d_strand != 0)
    {
      hit.d_first = hit.d_len - 1 - hit.d_first;
      hit.d_last = hit.d_len - 1 - hit.d_last;
    }
  }
  
  if ((parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastx))
  {
    if (hit.q_strand != 0)
    {
      hit.q_first = query.nt[0].len - 1 - (3 * hit.q_first) - hit.q_frame;
      hit.q_last = query.nt[0].len - 1 - (3 * hit.q_last) - hit.q_frame - 2;
    }
    else
    {
      hit.q_first = (3 * hit.q_first) + hit.q_frame;
      hit.q_last = (3 * hit.q_last) + hit.q_frame + 2;
    }
  }
  
  if ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx))
  {
    if (hit.d_strand != 0)
    {
      hit.d_first = hit.d_len_nt - 1 - (3 * hit.d_first) - hit.d_frame;
      hit.d_last = hit.d_len_nt - 1 - (3 * hit.d_last) - hit.d_frame - 2;
    }
    else
    {
      hit.d_first = (3 * hit.d_first) + hit.d_frame;
      hit.d_last = (3 * hit.d_last) + hit.d_frame + 2;
    }
  }

  hit.q_first++;
  hit.q_last++;
  hit.d_first++;
  hit.d_last++;

  long const maxqpos = hit.q_first > hit.q_last ? hit.q_first : hit.q_last; 
  long const maxdpos = hit.d_first > hit.d_last ? hit.d_first : hit.d_last; 
  long const maxpos = maxqpos > maxdpos ? maxqpos : maxdpos;
  decimal::Buffer buffer {{}};
  hit.poswidth = static_cast<int>(decimal::to_decimal(buffer, maxpos).size());

  return hit;
}

// an operation of an alignment string: a letter (M, D, I; 0 ends the
// alignment) and its count
struct AlignmentOperation
{
  char op;
  long len;
};

// show_align(): the alignment of a hit, in lines of ALIGNLEN columns
class AlignmentLines
{
public:
  explicit AlignmentLines(AlignedHit const & aligned) :
    hit(aligned), q_pos(aligned.q_align_start), d_pos(aligned.d_align_start) {}

  auto putalignop(AlignmentOperation const & operation) -> void;

private:
  AlignedHit hit;
  std::size_t line_pos = 0;
  long q_start = 0;
  long d_start = 0;
  long q_pos = 0;
  long d_pos = 0;
  std::array<char, ALIGNLEN + 1> q_line {{}};
  std::array<char, ALIGNLEN + 1> a_line {{}};
  std::array<char, ALIGNLEN + 1> d_line {{}};
};

auto AlignmentLines::putalignop(AlignmentOperation const & operation) -> void
{
  char const c = operation.op;
  long const len = operation.len;

  long count = len;
  while(count != 0)
  {
    if (line_pos == 0)
    {
      q_start = q_pos;
      d_start = d_pos;
    }

    char qs = 0;
    char ds = 0;

    switch(c)
    {
    case 'M':
      qs = hit.q_seq[q_pos++];
      ds = hit.d_seq[d_pos++];
      q_line[line_pos] = hit.sym[static_cast<int>(qs)];
      if (hit.symtype == SymbolType::blastn)
      {
	a_line[line_pos] = (qs == ds) ? '|' : ' ';
      }
      else
      {
	a_line[line_pos] = (qs == ds) ? hit.sym[static_cast<int>(qs)] : 
	  (pair_score(qs, ds) > 0 ? '+' : ' ');
      }
      d_line[line_pos] = hit.sym[static_cast<int>(ds)];
      line_pos++;
      break;

    case 'D':
      qs = hit.q_seq[q_pos++];
      q_line[line_pos] = hit.sym[static_cast<int>(qs)];
      a_line[line_pos] = ' ';
      d_line[line_pos] = '-';
      line_pos++;
      break;

    case 'I':
      ds = hit.d_seq[d_pos++];
      q_line[line_pos] = '-';
      a_line[line_pos] = ' ';
      d_line[line_pos] = hit.sym[static_cast<int>(ds)];
      line_pos++;
      break;
    default:
      break;
    }

    if ((line_pos == ALIGNLEN) || ((c == 0) && (line_pos > 0)))
    {
      // print alignment lines

      q_line[line_pos] = 0;
      a_line[line_pos] = 0;
      d_line[line_pos] = 0;

      long q1 = q_start + 1;
      long q2 = q_pos;

      long d1 = d_start + 1;
      long d2 = d_pos;

      if ((hit.symtype == SymbolType::blastn) && (hit.d_strand != 0))
      {
	d1 = hit.d_len - d1 + 1;
	d2 = hit.d_len - d2 + 1;
      }

      if ((hit.symtype == SymbolType::blastx) || (hit.symtype == SymbolType::tblastx))
      {
	if (hit.q_strand != 0)
	{
	  q1 = hit.q_len_nt - (3*q_start) - hit.q_frame;
	  q2 = hit.q_len_nt - (3*q_pos) - hit.q_frame + 1;
	}
	else
	{
	  q1 = (3*q_start) + hit.q_frame + 1;
	  q2 = (3*q_pos) + hit.q_frame;
	}
      }
      
      if ((hit.symtype == SymbolType::tblastn) || (hit.symtype == SymbolType::tblastx))
      {
	if (hit.d_strand != 0)
	{
	  d1 = hit.d_len_nt - (3*d_start) - hit.d_frame;
	  d2 = hit.d_len_nt - (3*d_pos) - hit.d_frame + 1;
	}
	else
	{
	  d1 = (3*d_start) + hit.d_frame + 1;
	  d2 = (3*d_pos) + hit.d_frame;
	}
      }


      fprint(out, "\n");
      // positions right-aligned on poswidth columns (was "%*ld")
      assert(hit.poswidth > 0);
      auto const width = static_cast<std::size_t>(hit.poswidth);
      fprint(out, "Query: ");
      fprint_integer(out, q1, width);
      fprint(out, ' ');
      fprint(out, as_c_string(q_line.data()));
      fprint(out, ' ');
      fprint_integer(out, q2);
      fprint(out, '\n');
      fprint(out, "       ");
      fprint_spaces(out, width);
      fprint(out, ' ');
      fprint(out, as_c_string(a_line.data()));
      fprint(out, '\n');
      fprint(out, "Sbjct: ");
      fprint_integer(out, d1, width);
      fprint(out, ' ');
      fprint(out, as_c_string(d_line.data()));
      fprint(out, ' ');
      fprint_integer(out, d2);
      fprint(out, '\n');

      line_pos = 0;
    }

    count--;
  }
}

// one operation of an alignment string ("M12D3...", align.cc): its
// letter and its count; the cursor moves past the count's digits
auto next_operation(char const * & cursor) -> AlignmentOperation
{
  AlignmentOperation operation {0, 0};
  operation.op = *cursor;
  cursor = std::next(cursor);
  char * end = nullptr;
  operation.len = std::strtol(cursor, & end, decimal_base);
  cursor = end;
  return operation;
}

auto show_align(AlignedHit const & hit) -> void
{
  AlignmentLines lines(hit);
  
  char const * p = hit.alignment;
  auto const * e = std::next(hit.alignment, static_cast<std::ptrdiff_t>(strlen(hit.alignment)));
  
  while(p < e)
  {
    auto const operation = next_operation(p);
    lines.putalignop(operation);
  }
  
  lines.putalignop({0, 1});
}

// the counts of an alignment, for its summaries
struct AlignmentCounts
{
  long identities = 0;
  long positives = 0;
  long indels = 0;  // gap positions
  long aligned = 0;  // alignment columns
  long gaps = 0;  // gap openings
};

// an alignment: its counts, and its query, middle and database lines
struct WholeAlignment
{
  AlignmentCounts counts;
  std::string qline;
  std::string aline;
  std::string dline;
};

auto whole_align(AlignedHit const & hit) -> WholeAlignment
{
  WholeAlignment alignment;
  auto & counts = alignment.counts;
  auto & qline = alignment.qline;
  auto & aline = alignment.aline;
  auto & dline = alignment.dline;

  long al = 0;
  char const * p = hit.alignment;
  while((*p) != 0)
  {
    al += next_operation(p).len;
  }

  for (auto * line : {&qline, &aline, &dline})
  {
    line->reserve(static_cast<std::size_t>(al));
  }
  
  long q_pos = hit.q_align_start;
  long d_pos = hit.d_align_start;
  
  p = hit.alignment;

  while((*p) != 0)
  {
    auto const operation = next_operation(p);
    char const op = operation.op;
    long const len = operation.len;
    
    counts.aligned += len;
    if (op == 'D')
    {
      for(long j=0; j<len; j++)
      {
	char const qs = hit.q_seq[q_pos++];
	qline += hit.sym[static_cast<int>(qs)];
	aline += ' ';
	dline += '-';
      }
      counts.gaps += 1;
      counts.indels += len;
    }
    else if (op == 'I')
    {
      for(long j=0; j<len; j++)
      {
	char const ds = hit.d_seq[d_pos++];
	qline += '-';
	aline += ' ';
	dline += hit.sym[static_cast<int>(ds)];
      }
      counts.gaps += 1;
      counts.indels += len;
    }
    else if (op == 'M')
    {
      for(long j=0; j<len; j++)
      {
	char const qs = hit.q_seq[q_pos++];
	char const ds = hit.d_seq[d_pos++];
	qline += hit.sym[static_cast<int>(qs)];
	if (qs == ds)
	{
	  aline += '|';
	  counts.identities++;
	  counts.positives++;
	}
	else if (pair_score(qs, ds) > 0)
	{
	  aline += '+';
	  counts.positives++;
	}
	else
	{
	  aline += ' ';
	}
	dline += hit.sym[static_cast<int>(ds)];
      }
    }
    else
    {
      fatal("Illegal alignment string.");
    }
  }

  return alignment;
}

auto count_align(AlignedHit const & hit) -> AlignmentCounts
{
  AlignmentCounts counts;
  
  long q_pos = hit.q_align_start;
  long d_pos = hit.d_align_start;
  
  char const * p = hit.alignment;
  auto const * e = std::next(hit.alignment, static_cast<std::ptrdiff_t>(strlen(hit.alignment)));

  while(p < e)
  {
    auto const operation = next_operation(p);
    char const op = operation.op;
    long const len = operation.len;
    
    counts.aligned += len;
    if (op == 'D')
    {
      counts.gaps += 1;
      counts.indels += len;
      q_pos += len;
    }
    else if (op == 'I')
    {
      counts.gaps += 1;
      counts.indels += len;
      d_pos += len;
    }
    else
    {
      for(long j=0; j<len; j++)
      {
	char const qs = hit.q_seq[q_pos++];
	char const ds = hit.d_seq[d_pos++];
	if (qs == ds)
	{
	  counts.identities++;
	  counts.positives++;
	}
	else if (pair_score(qs, ds) > 0)
	{
	  counts.positives++;
	}
      }
    }
  }
  return counts;
}

auto hits_show_expect(double expect_value) -> void
{
  // the format of an expect value depends on its range: each bound is
  // where the rounded value would need one more character
  constexpr double zero_below = 1e-180;
  constexpr double three_digit_exponent_below = 9.5e-100;
  constexpr double exponent_below = 0.00095;
  constexpr double three_decimals_below = 0.0995;
  constexpr double two_decimals_below = 0.95;
  constexpr double one_decimal_below = 9.5;
  // "%-6.0e" of a value from 1e-180 on: at most 6 characters and the NUL
  constexpr std::size_t exponent_text_size = 10;

  std::array<char, exponent_text_size> temp {{}};
  if (expect_value < zero_below)
  {
    fprint(out, "0.0  ");
  }
  else if (expect_value < three_digit_exponent_below)
  {
    // C++17 refactoring: replace the printf() formats of this function with std::to_chars
    snprintf(temp.data(), temp.size(), "%-6.0e", expect_value);
    fprint(out, as_c_string(std::next(temp.data())));  // without the first character
  }
  else if (expect_value < exponent_below)
  {
    fprintf(out, "%-5.0e", expect_value);
  }
  else if (expect_value < three_decimals_below)
  {
    fprintf(out, "%-5.3f", expect_value);
  }
  else if (expect_value < two_decimals_below)
  {
    fprintf(out, "%-5.2f", expect_value);
  }
  else if (expect_value < one_decimal_below)
  {
    fprintf(out, "%-5.1f", expect_value);
  }
  else
  {
    fprintf(out, "%5.0f", expect_value);
  }
}

}  // anonymous namespace

auto xml_putc(char const symbol) noexcept -> void
{
  switch (symbol)
    {
    case '&':
      fprint(out, "&amp;");
      break;
    case '<':
      fprint(out, "&lt;");
      break;
    case '>':
      fprint(out, "&gt;");
      break;
    case '"':
      fprint(out, "&quot;");
      break;
    case '\'':
      fprint(out, "&apos;");
      break;
    default:
      fprint(out, symbol);
      break;
    }
}

namespace {

// print at most max_length characters of text, escaped as XML (KI-27);
// the text is truncated before it is escaped
auto xml_print(View<char> const text,
               std::size_t const max_length = std::numeric_limits<std::size_t>::max()) noexcept -> void
{
  for (auto const symbol : text.first(std::min(max_length, text.size())))
  {
    xml_putc(symbol);
  }
}

// ParAlign XML (-m 99): the length of the short name of a hit (the
// start of its title)
constexpr std::size_t short_name_length = 35;

// the anchor of a hit: "query_hit_frame_strand..." (numbers and marks)
auto make_anchor(SymbolType symbol_type, long query_index, long i) -> std::string
{
  auto const & hit = hit_entry(i);
  auto const sign = [](long const strand) -> char
  {
    return (strand != 0) ? '-' : '+';
  };
  auto const prefix = std::to_string(query_index) + "_" + std::to_string(hit.seqno);

  switch(symbol_type)
  {
  case SymbolType::blastn:
    // blastn: the strand of a hit is stored as its database strand
    // (KI-29)
    return prefix + "__" + sign(hit.dstrand) + "__+";
  case SymbolType::blastx:
    return prefix + "_" + std::to_string(hit.qframe + 1) + "_" + sign(hit.qstrand) + "__";
  case SymbolType::tblastn:
    return prefix + "___" + std::to_string(hit.dframe + 1) + "_" + sign(hit.dstrand);
  case SymbolType::tblastx:
    return prefix + "_" + std::to_string(hit.qframe + 1) + "_" + sign(hit.qstrand)
      + "_" + std::to_string(hit.dframe + 1) + "_" + sign(hit.dstrand);
  default:
    return prefix + "____";
  }
}

// the parts of a defline: its gi (0: none), the link (the identifier
// before the first space, empty when there is no space), and the rest
// (the title)
struct DeflineParts
{
  long gi;
  View<char> link;
  View<char> rest;
};

auto hits_defline_split(char const * defline) -> DeflineParts
{
  char const * p = defline;

  // no gi (KI-42: it kept the gi of the previous defline)
  DeflineParts parts {0, View<char>{}, View<char>{}};
  
  // "gi|" and a number, as the header parser writes them (set_id(),
  // asnparse.cc)
  constexpr std::size_t gi_prefix_length = 3;
  if (std::strncmp(p, "gi|", gi_prefix_length) == 0)
  {
    auto * const number = std::next(p, gi_prefix_length);
    char * end = nullptr;
    auto const value = std::strtol(number, & end, decimal_base);
    if (end != number)
    {
      parts.gi = value;
      p = end;
    }
  }

  if (*p == '|')
  {
    p = std::next(p);
  }

  auto const * const r = strchr(p, ' ');
  if (r != nullptr)
  {
    parts.link = View<char>{p, static_cast<std::size_t>(std::distance(p, r))};
    parts.rest = as_c_string(std::next(r));
  }
  else
  {
    parts.rest = as_c_string(p);
  }

  return parts;
}

// the numbers of hits shown: descriptions (-v) and alignments (-b),
// at most the hits found
struct ShownHits
{
  long descriptions;
  long alignments;
};

// the tabular output (-m 8, -m 9): with or without comment lines
enum struct TabularComments : bool { without, with };

auto hits_show_xml_paralign(Parameters const & parameters,
			    ShownHits const & shown,
			    db_thread_s const & t) -> void
{
  /* ParAlign XML */
  
  fprint(out, "\t<paralignOutput>\n");
  
  char const * qseqtypedescr = nullptr;
  bool const protein_query = (query.symtype == SymbolType::blastp) || (query.symtype == SymbolType::tblastn) || (query.symtype == SymbolType::sound);
  auto const & q = protein_query ? query.aa[0] : query.nt[0];
  if ((query.symtype == SymbolType::blastp) || (query.symtype == SymbolType::tblastn))
  {
    qseqtypedescr = "Amino Acid";
  }
  else if (query.symtype == SymbolType::sound)
  {
    /* sound queries are stored as amino acid queries (KI-30) */
    qseqtypedescr = "Sound";
  }
  else
  {
    qseqtypedescr = "Nucleotide";
  }
  
  fprint(out, "\t\t<queryInformation>\n");
  fprint(out, "\t\t\t<queryFilename>");
  xml_print(as_c_string(parameters.queryname));
  fprint(out, "</queryFilename>\n");
  fprint(out, "\t\t\t<querySequencetype>");
  fprint(out, as_c_string(qseqtypedescr));
  fprint(out, "</querySequencetype>\n");
  fprint(out, "\t\t\t<queryDescription>");
  xml_print(as_c_string(query.description));
  fprint(out, "</queryDescription>\n");
  fprint(out, "\t\t\t<queryLength>");
  fprint_integer(out, q.len);
  fprint(out, "</queryLength>\n");
  fprint(out, "\t\t\t<querySequence>");
  for (auto const residue : q.view())
  {
    fprint(out, query.sym[static_cast<int>(residue)]);
  }
  fprint(out, "</querySequence>\n");
  fprint(out, "\t\t</queryInformation>\n");
  
  char const * dbseqtypedescr = nullptr;
  char const * ncbidb = nullptr;
  char const * ncbiopt = nullptr;
  if ((query.symtype == SymbolType::blastn) || (query.symtype == SymbolType::tblastn) || (query.symtype == SymbolType::tblastx))
  {
    dbseqtypedescr = "Nucleotide";
    ncbidb = "Nucleotide";
    ncbiopt = "GenBank";
  }
  else
  {
    dbseqtypedescr = "Amino Acid";
    ncbidb = "Protein";
    ncbiopt = "GenPept";
  }

  /* sound databases are stored as amino acid databases (KI-30) */
  if (query.symtype == SymbolType::sound)
  {
    dbseqtypedescr = "Sound";
  }
  fprint(out, "\t\t<databaseInformation>\n");
  fprint(out, "\t\t\t<databaseFilename>");
  xml_print(as_c_string(parameters.databasename));
  fprint(out, "</databaseFilename>\n");
  fprint(out, "\t\t\t<databaseSequencetype>");
  fprint(out, as_c_string(dbseqtypedescr));
  fprint(out, "</databaseSequencetype>\n");
  fprint(out, "\t\t\t<databaseDescription>");
  xml_print(as_c_string(db_gettitle()));
  fprint(out, "</databaseDescription>\n");
  fprint(out, "\t\t\t<databaseVersion>");
  fprint_integer(out, db_getversion());
  fprint(out, "</databaseVersion>\n");
  fprint(out, "\t\t\t<databaseDate>");
  xml_print(as_c_string(db_gettime()));
  fprint(out, "</databaseDate>\n");
  fprint(out, "\t\t\t<residueCount>");
  fprint_integer(out, db_getsymcount_masked());
  fprint(out, "</residueCount>\n");
  fprint(out, "\t\t\t<sequenceCount>");
  fprint_integer(out, db_getseqcount_masked());
  fprint(out, "</sequenceCount>\n");
  fprint(out, "\t\t\t<longestSequenceLength>");
  fprint_integer(out, db_getlongest());
  fprint(out, "</longestSequenceLength>\n");
  fprint(out, "\t\t</databaseInformation>\n");
  
  char const * strands = "";
  switch(parameters.querystrands)
  {
  case QueryStrands::plus:
    strands = "Plus";
    break;
  case QueryStrands::minus:
    strands = "Minus";
    break;
  case QueryStrands::both:
    strands = "Both";
    break;
  default:
    break;
  }

  fprint(out, "\t\t<options>\n");
  fprint(out, "\t\t\t<algorithm>Smith-Waterman</algorithm>\n");

  if ((parameters.symtype == SymbolType::blastn) || (parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastx))
  {
    fprint(out, "\t\t\t<queryStrands>");
    fprint(out, as_c_string(strands));
    fprint(out, "</queryStrands>\n");
  }

  if (parameters.symtype == SymbolType::blastn)
  {
    fprint(out, "\t\t\t<scoreMatrix>NT</scoreMatrix>\n");
  }
  else
  {
    fprint(out, "\t\t\t<scoreMatrix>");
    xml_print(as_c_string(parameters.matrixname));
    fprint(out, "</scoreMatrix>\n");
  }

  fprint(out, "\t\t\t<gapPenalties>\n");
  fprint(out, "\t\t\t\t<gapPenaltyOpen>");
  fprint_integer(out, parameters.gapopen);
  fprint(out, "</gapPenaltyOpen>\n");
  fprint(out, "\t\t\t\t<gapPenaltyExtension>");
  fprint_integer(out, parameters.gapextend);
  fprint(out, "</gapPenaltyExtension>\n");
  fprint(out, "\t\t\t\t<ungapped>\n");
  fprintf(out, "\t\t\t\t\t<ungappedLambda>%.4g</ungappedLambda>\n", statistics.ungapped_lambda);
  fprintf(out, "\t\t\t\t\t<ungappedKappa>%.4g</ungappedKappa>\n", statistics.ungapped_K);
  fprintf(out, "\t\t\t\t\t<ungappedEta>%.4g</ungappedEta>\n", statistics.ungapped_H);
  fprint(out, "\t\t\t\t</ungapped>\n");
  fprint(out, "\t\t\t\t<gapped>\n");
  fprintf(out, "\t\t\t\t\t<gappedLambda>%.4g</gappedLambda>\n", statistics.lambda);
  fprintf(out, "\t\t\t\t\t<gappedKappa>%.4g</gappedKappa>\n", statistics.K);
  fprintf(out, "\t\t\t\t\t<gappedEta>%.4g</gappedEta>\n", statistics.H);
  fprint(out, "\t\t\t\t</gapped>\n");

  fprint(out, "\t\t\t</gapPenalties>\n");
  fprint(out, "\t\t\t<expectRange>\n");
  fprintf(out, "\t\t\t\t<expectRangeFrom>%.2g</expectRangeFrom>\n", parameters.minexpect);
  fprintf(out, "\t\t\t\t<expectRangeTo>%.2g</expectRangeTo>\n", parameters.expect);
  fprint(out, "\t\t\t</expectRange>\n");
  fprint(out, "\t\t\t<displayLimits>\n");
  fprint(out, "\t\t\t\t<hitLimit>");
  fprint_integer(out, parameters.maxmatches);
  fprint(out, "</hitLimit>\n");
  fprint(out, "\t\t\t\t<alignmentLimit>");
  fprint_integer(out, parameters.alignments);
  fprint(out, "</alignmentLimit>\n");
  fprint(out, "\t\t\t\t<subalignmentLimit>");
  fprint_integer(out, static_cast<long>(1));
  fprint(out, "</subalignmentLimit>\n");
  fprint(out, "\t\t\t</displayLimits>\n");
  fprint(out, "\t\t\t<threads>");
  fprint_integer(out, parameters.threads);
  fprint(out, "</threads>\n");
  fprint(out, "\t\t</options>\n");

  fprint(out, "\t\t\t<searchInformation>\n");
  fprint(out, "\t\t\t\t<searchStarted>");
  fprint(out, as_c_string(run.ti.starttime.data()));
  fprint(out, "</searchStarted>\n");
  fprint(out, "\t\t\t\t<searchCompleted>");
  fprint(out, as_c_string(run.ti.endtime.data()));
  fprint(out, "</searchCompleted>\n");
  fprintf(out, "\t\t\t\t<searchElapsedTime>%.2fs</searchElapsedTime>\n", run.ti.elapsed);
  if (run.ti.elapsed > 0.0)
  {
    fprintf(out, "\t\t\t\t<searchSpeed>%.3f GCUPS</searchSpeed>\n", gcups(run.ti.speed));
  }
  else
  {
    fprint(out, "\t\t\t\t<searchSpeed>n/a</searchSpeed>\n");
  }
  fprint(out, "\t\t\t\t<searchSWAlignments>\n");
  fprint(out, "\t\t\t\t\t<SWAbsolute>");
  fprint_integer(out, run.compute7);
  fprint(out, "</SWAbsolute>\n");
  fprint(out, "\t\t\t\t\t<SWPercent>100</SWPercent>\n");
  fprint(out, "\t\t\t\t</searchSWAlignments>\n");
  fprint(out, "\t\t\t</searchInformation>\n");

  fprint(out, "\t\t<resultInformation>\n");
  fprint(out, "\t\t\t<resultHits>\n");
  fprint(out, "\t\t\t\t<totalCount>");
  fprint_integer(out, run.totalhits);
  fprint(out, "</totalCount>\n");
  fprint(out, "\t\t\t\t<obviousCount>");
  fprint_integer(out, hit_list.obvious);
  fprint(out, "</obviousCount>\n");
  fprint(out, "\t\t\t\t<shownCount>");
  fprint_integer(out, shown.descriptions);
  fprint(out, "</shownCount>\n");
  fprint(out, "\t\t\t</resultHits>\n");
  fprint(out, "\t\t\t<alignmentCount>");
  fprint_integer(out, shown.alignments);
  fprint(out, "</alignmentCount>\n");
  fprint(out, "\t\t</resultInformation>\n");
  
  fprint(out, "\t\t<shortVersionHits>\n");
  
  for(long i=0; i<shown.descriptions; i++)
  {
    auto const score = hit_entry(i).score;
    auto const e = expect_value_of(score);

    auto const anchor = make_anchor(query.symtype, run.queryno, i);

    auto const deflinetable = db_parse_header(t, make_view(hit_entry(i).header_address), 1);
    auto const parts = hits_defline_split(deflinetable[0].c_str());

    fprint(out, "\t\t\t<shortVersionHit>\n");
    fprint(out, "\t\t\t\t<shortVersionAnchor>");
    fprint(out, make_view(anchor));
    fprint(out, "</shortVersionAnchor>\n");
    if (parts.gi != 0)
      {
    fprint(out, "\t\t\t\t<shortVersionLink>\n");
    fprint(out, "\t\t\t\t\t<shortVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=");
    fprint(out, as_c_string(ncbidb));
    fprint(out, "&amp;list_uids=");
    fprint_integer(out, parts.gi);
    fprint(out, "&amp;dopt=");
    fprint(out, as_c_string(ncbiopt));
    fprint(out, "</shortVersionLinkDestination>\n");
    fprint(out, "\t\t\t\t\t<shortVersionLinkText>gi|");
    fprint_integer(out, parts.gi);
    fprint(out, "</shortVersionLinkText>\n");
    fprint(out, "\t\t\t\t</shortVersionLink>\n");
      }
    fprint(out, "\t\t\t\t<shortVersionLink>\n");
    fprint(out, "\t\t\t\t\t<shortVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&amp;db=");
    fprint(out, as_c_string(ncbidb));
    fprint(out, "&amp;term=");
    xml_print(parts.link);
    fprint(out, "&amp;doptcmdl=");
    fprint(out, as_c_string(ncbiopt));
    fprint(out, "</shortVersionLinkDestination>\n");
    fprint(out, "\t\t\t\t\t<shortVersionLinkText>");
    xml_print(parts.link);
    fprint(out, "</shortVersionLinkText>\n");
    fprint(out, "\t\t\t\t</shortVersionLink>\n");
    fprint(out, "\t\t\t\t<shortVersionName>");
    xml_print(parts.rest, short_name_length);
    fprint(out, "</shortVersionName>\n");
    if (parameters.symtype == SymbolType::blastn)
    {
      fprint(out, "\t\t\t\t<shortVersionStrand>");
      fprint(out, (hit_entry(i).dstrand != 0) ? '-' : '+');
      fprint(out, "</shortVersionStrand>\n");
    }
    else if (parameters.symtype == SymbolType::blastx)
    {
      fprint(out, "\t\t\t\t<shortVersionFrame>");
      fprint(out, (hit_entry(i).qstrand != 0) ? '-' : '+');
      fprint_integer(out, hit_entry(i).qframe+1);
      fprint(out, "</shortVersionFrame>\n");
    }
    else if (parameters.symtype == SymbolType::tblastn)
    {
      fprint(out, "\t\t\t\t<shortVersionFrame>");
      fprint(out, (hit_entry(i).dstrand != 0) ? '-' : '+');
      fprint_integer(out, hit_entry(i).dframe+1);
      fprint(out, "</shortVersionFrame>\n");
    }
    else if (parameters.symtype == SymbolType::tblastx)
    {
      fprint(out, "\t\t\t\t<shortVersionFrame>");
      fprint(out, (hit_entry(i).qstrand != 0) ? '-' : '+');
      fprint_integer(out, hit_entry(i).qframe+1);
      fprint(out, '/');
      fprint(out, (hit_entry(i).dstrand != 0) ? '-' : '+');
      fprint_integer(out, hit_entry(i).dframe+1);
      fprint(out, "</shortVersionFrame>\n");
    }
    fprint(out, "\t\t\t\t<shortVersionScore>");
    fprint_integer(out, score);
    fprint(out, "</shortVersionScore>\n");
    fprintf(out, "\t\t\t\t<shortVersionEValue>%.2g</shortVersionEValue>\n", e);
    fprint(out, "\t\t\t</shortVersionHit>\n");

  }

  fprint(out, "\t\t</shortVersionHits>\n");

  if (shown.alignments != 0)
  {
    fprint(out, "\t\t<longVersionHits>\n");
    
    for(long i=0; i<shown.alignments; i++)
    {
      
      auto const anchor = make_anchor(query.symtype, run.queryno, i);
      
      fprint(out, "\t\t\t<longVersionHit>\n");
      fprint(out, "\t\t\t\t<longVersionAnchor>");
      fprint(out, make_view(anchor));
      fprint(out, "</longVersionAnchor>\n");
      
      auto const deflinetable = db_parse_header(t, make_view(hit_entry(i).header_address), 1);
      fprint(out, "\t\t\t\t<linkContainer>\n");
      
      for (auto const & defline : deflinetable)
      {
	auto const parts = hits_defline_split(defline.c_str());
  
        if (parts.gi != 0)
	{
          fprint(out, "\t\t\t\t\t<longVersionLink>\n");
	  fprint(out, "\t\t\t\t\t\t<longVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=");
	  fprint(out, as_c_string(ncbidb));
	  fprint(out, "&amp;list_uids=");
	  fprint_integer(out, parts.gi);
	  fprint(out, "&amp;dopt=");
	  fprint(out, as_c_string(ncbiopt));
	  fprint(out, "</longVersionLinkDestination>\n");
	  fprint(out, "\t\t\t\t\t\t<longVersionLinkText>gi|");
	  fprint_integer(out, parts.gi);
	  fprint(out, "</longVersionLinkText>\n");
	  fprint(out, "\t\t\t\t\t</longVersionLink>\n");
	}
      
	fprint(out, "\t\t\t\t\t<longVersionLink>\n");
	fprint(out, "\t\t\t\t\t\t<longVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&amp;db=");
	fprint(out, as_c_string(ncbidb));
	fprint(out, "&amp;term=");
	xml_print(parts.link);
	fprint(out, "&amp;doptcmdl=");
	fprint(out, as_c_string(ncbiopt));
	fprint(out, "</longVersionLinkDestination>\n");
	fprint(out, "\t\t\t\t\t\t<longVersionLinkText>");
	xml_print(parts.link);
	fprint(out, "</longVersionLinkText>\n");
	fprint(out, "\t\t\t\t\t</longVersionLink>\n");
      
	fprint(out, "\t\t\t\t\t<longVersionName>");
	xml_print(parts.rest);
	fprint(out, "</longVersionName>\n");
      }
        
      fprint(out, "\t\t\t\t</linkContainer>\n");
    
      auto const dlen = hit_entry(i).dlen;
      auto const dlennt = hit_entry(i).dlennt;

      if (parameters.symtype == SymbolType::blastn)
      {
	fprint(out, "\t\t\t\t<databaseSequenceLength>");
	fprint_integer(out, dlen);
	fprint(out, " nt</databaseSequenceLength>\n");
      }
      else if ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx))
      {
	fprint(out, "\t\t\t\t<databaseSequenceLength>");
	fprint_integer(out, dlennt);
	fprint(out, " nt</databaseSequenceLength>\n");
      }
      else
      {
	fprint(out, "\t\t\t\t<databaseSequenceLength>");
	fprint_integer(out, dlen);
	fprint(out, " aa</databaseSequenceLength>\n");
      }

      if (parameters.symtype == SymbolType::blastn)
      {
	fprint(out, "\t\t\t\t<alignmentMatchLocation>");
	fprint(out, as_c_string((hit_entry(i).dstrand != 0) ? "Matches on complementary strands." : "Matches on same strands."));
	fprint(out, "</alignmentMatchLocation>\n");
      }
      else if ((parameters.symtype>=SymbolType::blastx) && (parameters.symtype<=SymbolType::tblastx))
      {
	fprint(out, "\t\t\t\t<longVersionFrames>\n");

	if ((parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastx))
	{
	  fprint(out, "\t\t\t\t\t<longVersionQueryFrame>\n");
	  fprint(out, "\t\t\t\t\t\t<queryStrand>");
	  fprint(out, (hit_entry(i).qstrand != 0) ? '-' : '+');
	  fprint(out, "</queryStrand>\n");
	  fprint(out, "\t\t\t\t\t\t<queryFrame>");
	  fprint_integer(out, hit_entry(i).qframe+1);
	  fprint(out, "</queryFrame>\n");
	  fprint(out, "\t\t\t\t\t</longVersionQueryFrame>\n");
	}
	
	if ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx))
	{
	  fprint(out, "\t\t\t\t\t<longVersionDatabaseFrame>\n");
	  fprint(out, "\t\t\t\t\t\t<databaseStrand>");
	  fprint(out, (hit_entry(i).dstrand != 0) ? '-' : '+');
	  fprint(out, "</databaseStrand>\n");
	  fprint(out, "\t\t\t\t\t\t<databaseFrame>");
	  fprint_integer(out, hit_entry(i).dframe+1);
	  fprint(out, "</databaseFrame>\n");
	  fprint(out, "\t\t\t\t\t</longVersionDatabaseFrame>\n");
	}

	fprint(out, "\t\t\t\t</longVersionFrames>\n");
      }

      auto const score = hit_entry(i).score;
      auto const e = expect_value_of(score);

    
        
      auto const hit = aligned_hit(parameters, i);
      auto const alignment = whole_align(hit);
      auto const & counts = alignment.counts;

      fprint(out, "\t\t\t\t<alignment>\n");
      fprint(out, "\t\t\t\t\t<subalignment>\n");
      fprint(out, "\t\t\t\t\t\t<longVersionScore>");
      fprint_integer(out, score);
      fprint(out, "</longVersionScore>\n");
      fprintf(out, "\t\t\t\t\t\t<longVersionEValue>%.2g</longVersionEValue>\n", e);
      fprint(out, "\t\t\t\t\t\t<identical>\n");
      fprint(out, "\t\t\t\t\t\t\t<identicalNominator>");
      fprint_integer(out, counts.identities);
      fprint(out, "</identicalNominator>\n");
      fprint(out, "\t\t\t\t\t\t\t<identicalDenominator>");
      fprint_integer(out, counts.aligned);
      fprint(out, "</identicalDenominator>\n");
      fprintf(out, "\t\t\t\t\t\t\t<identicalPercentage>%.1f</identicalPercentage>\n", percentage(counts.identities, counts.aligned));
      fprint(out, "\t\t\t\t\t\t</identical>\n");

      if (parameters.symtype != SymbolType::blastn)
      {
	fprint(out, "\t\t\t\t\t\t<positive>\n");
	fprint(out, "\t\t\t\t\t\t\t<positiveNominator>");
	fprint_integer(out, counts.positives);
	fprint(out, "</positiveNominator>\n");
	fprint(out, "\t\t\t\t\t\t\t<positiveDenominator>");
	fprint_integer(out, counts.aligned);
	fprint(out, "</positiveDenominator>\n");
	fprintf(out, "\t\t\t\t\t\t\t<positivePercentage>%.1f</positivePercentage>\n", percentage(counts.positives, counts.aligned));
	fprint(out, "\t\t\t\t\t\t</positive>\n");
      }

      fprint(out, "\t\t\t\t\t\t<indels>\n");
      fprint(out, "\t\t\t\t\t\t\t<indelsNominator>");
      fprint_integer(out, counts.indels);
      fprint(out, "</indelsNominator>\n");
      fprint(out, "\t\t\t\t\t\t\t<indelsDenominator>");
      fprint_integer(out, counts.aligned);
      fprint(out, "</indelsDenominator>\n");
      fprintf(out, "\t\t\t\t\t\t\t<indelsPercentage>%.1f</indelsPercentage>\n", percentage(counts.indels, counts.aligned));
      fprint(out, "\t\t\t\t\t\t</indels>\n");
      fprint(out, "\t\t\t\t\t\t<gaps>");
      fprint_integer(out, counts.gaps);
      fprint(out, "</gaps>\n");
      fprint(out, "\t\t\t\t\t\t<alignmentQuery>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentQueryStart>");
      fprint_integer(out, hit.q_first);
      fprint(out, "</alignmentQueryStart>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentQueryLine>");
      fprint(out, as_c_string(alignment.qline.c_str()));
      fprint(out, "</alignmentQueryLine>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentQueryEnd>");
      fprint_integer(out, hit.q_last);
      fprint(out, "</alignmentQueryEnd>\n");
      fprint(out, "\t\t\t\t\t\t</alignmentQuery>\n");
      fprint(out, "\t\t\t\t\t\t<alignmentLine>");
      fprint(out, as_c_string(alignment.aline.c_str()));
      fprint(out, "</alignmentLine>\n");
      fprint(out, "\t\t\t\t\t\t<alignmentDatabase>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentDatabaseStart>");
      fprint_integer(out, hit.d_first);
      fprint(out, "</alignmentDatabaseStart>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentDatabaseLine>");
      fprint(out, as_c_string(alignment.dline.c_str()));
      fprint(out, "</alignmentDatabaseLine>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentDatabaseEnd>");
      fprint_integer(out, hit.d_last);
      fprint(out, "</alignmentDatabaseEnd>\n");
      fprint(out, "\t\t\t\t\t\t</alignmentDatabase>\n");
      fprint(out, "\t\t\t\t\t</subalignment>\n");
      fprint(out, "\t\t\t\t</alignment>\n");
      fprint(out, "\t\t\t</longVersionHit>\n");


    }

    fprint(out, "\t\t</longVersionHits>\n");
  }

  fprint(out, "\t</paralignOutput>\n");
}

// the query id ends at the first whitespace character (space, tab,
// ...), as in BLAST (KI-20)
auto ends_query_id(char const symbol) -> bool
{
  return (symbol == '\0') or
    (std::isspace(static_cast<unsigned char>(symbol)) != 0);
}

// the query id: the description up to its first whitespace character
auto query_id(std::string const & description) -> View<char>
{
  auto const end = std::find_if(description.begin(), description.end(), ends_query_id);
  return make_view(description).first(static_cast<std::size_t>(std::distance(description.begin(), end)));
}

auto show_description(std::string const & description) -> void
{
  fprint(out, query_id(description));
}

// query id (the description up to its first whitespace character),
// escaped as XML
// (KI-27)
auto show_description_xml(std::string const & description) -> void
{
  xml_print(query_id(description));
}

auto hits_show_xml(Parameters const & parameters,
		   ShownHits const & shown,
		   db_thread_s const & t) -> void
{
  /* Simple XML */
  
  fprint(out, "<result>\n");
  fprint(out, "  <general>\n");
  fprint(out, "    <hitcount>");
  fprint_integer(out, hit_list.count);
  fprint(out, "</hitcount>\n");
  fprint(out, "  </general>\n");
  fprint(out, "  <hits>\n");
  
  for(long i=0; i<shown.descriptions; i++)
  {
    auto const seqno = hit_entry(i).seqno;
    auto const score = hit_entry(i).score;
    // the database sequence length in nucleotides for tblastn and
    // tblastx, as in the other outputs (KI-36)
    long const dlen = ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx)) ?
      hit_entry(i).dlennt : hit_entry(i).dlen;
    
    fprint(out, "    <hit>\n");
    fprint(out, "      <hitno>");
    fprint_integer(out, i+1);
    fprint(out, "</hitno>\n");
    fprint(out, "      <track>");
    fprint_integer(out, seqno);
    fprint(out, "</track>\n");
    fprint(out, "      <query>");
    show_description_xml(query.description);
    fprint(out, "</query>\n");
    fprint(out, "      <name>");
    HeaderLayout layout;
    layout.show_gis = parameters.show_gis;
    layout.escaping = Escaping::xml;
    db_showheader(t, make_view(hit_entry(i).header_address), layout);
    fprint(out, "</name>\n");
    fprint(out, "      <len>");
    fprint_integer(out, dlen);
    fprint(out, "</len>\n");
    fprint(out, "      <score>");
    fprint_integer(out, score);
    fprint(out, "</score>\n");
    
    if (i < shown.alignments)
    {

        
      auto const hit = aligned_hit(parameters, i);
      auto const alignment = whole_align(hit);

      fprint(out, "      <alignment>");
      fprint(out, as_c_string(hit_entry(i).alignment.c_str()));
      fprint(out, "</alignment>\n");

      fprint(out, "      <qpos>");
      fprint_integer(out, hit.q_first);
      fprint(out, ',');
      fprint_integer(out, hit.q_last);
      fprint(out, "</qpos>\n");
      fprint(out, "      <dpos>");
      fprint_integer(out, hit.d_first);
      fprint(out, ',');
      fprint_integer(out, hit.d_last);
      fprint(out, "</dpos>\n");
      
      fprint(out, "      <qseq>");
      fprint(out, as_c_string(alignment.qline.c_str()));
      fprint(out, "</qseq>\n");
      fprint(out, "      <aseq>");
      fprint(out, as_c_string(alignment.aline.c_str()));
      fprint(out, "</aseq>\n");
      fprint(out, "      <dseq>");
      fprint(out, as_c_string(alignment.dline.c_str()));
      fprint(out, "</dseq>\n");

    }
    fprint(out, "    </hit>\n");
  }
  fprint(out, "  </hits>\n");
  fprint(out, "</result>\n");
}

auto hits_show_tsv(Parameters const & parameters,
		   ShownHits const & shown,
		   TabularComments const comments,
		   db_thread_s const & t) -> void
{
  if (comments == TabularComments::with)
    {
      constexpr char const * ref = "Reference: T. Rognes (2011) Faster Smith-Waterman database searches with inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.";
      fprint(out, "# ");
      fprint(out, as_c_string(swipe_name_and_version));
      fprint(out, " - ");
      fprint(out, as_c_string(ref));
      fprint(out, '\n');
      fprint(out, "# Query: ");
      fprint(out, as_c_string(query.description.c_str()));
      fprint(out, '\n');
      fprint(out, "# Database: ");
      fprint(out, as_c_string(parameters.databasename));
      fprint(out, '\n');
      if (statistics.available != 0)
      {
	fprint(out, "# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n");
      }
      else
      {
	fprint(out, "# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, score\n");
      }
    }

  for(long i=0; i<shown.alignments; i++)
  {
    show_description(query.description);
    fprint(out, '\t');
    HeaderLayout layout;
    layout.show_gis = 1;
    layout.text = DeflineText::identifier;
    db_showheader(t, make_view(hit_entry(i).header_address), layout);
    
    
    auto const hit = aligned_hit(parameters, i);
    auto const counts = count_align(hit);
    
    auto const score = hit_entry(i).score;
    
    fprintf(out, "\t%.2f\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", 
	    percentage(counts.identities, counts.aligned),
	    counts.aligned,
	    counts.aligned - counts.identities - counts.indels,
	    counts.gaps,
	    hit.q_first,
	    hit.q_last,
	    hit.d_first,
	    hit.d_last);
    
    if (statistics.available != 0)
    {
      auto const expect_value = expect_value_of(score);
      fprintf(out, "\t%.2g", expect_value);
      auto const bits = bit_score_of(score);
      fprintf(out, "\t%.1f", bits);
    }
    else
    {
      fprint(out, "\t");
      fprint_integer(out, score);
    }

    fprint(out, "\n");
  }
}

auto hits_show_plain(Parameters const & parameters,
		     ShownHits const & shown,
		     db_thread_s const & t) -> void
{
    if (hit_list.count == 0)
    {
      fprint(out, "\nNo hits.\n");
    }
    else
    {
      if (statistics.available != 0)
      {
	fprint(out, "                                                                 Score    E\n");
	fprint(out, "Sequences producing significant alignments:                      (bits) Value\n\n");
      }
      else
      {
	fprint(out, "Sequences producing significant alignments:                         Score\n\n");
      }
	  
      assert(shown.descriptions <= hit_list.count);
      for (auto const & hit : make_view(hit_list.entries).first(static_cast<std::size_t>(shown.descriptions)))
      {
	long const headerlen = description_width - frame_mark_width(parameters.symtype);

	HeaderLayout layout;
	layout.show_gis = parameters.show_gis;
	layout.maxlen = headerlen;
	layout.linelen = headerlen;
	db_showheader(t, 
		      make_view(hit.header_address), layout);

	auto const score = hit.score;

	if (parameters.symtype == SymbolType::blastn)
	{
	  fprint(out, ' ');
	  fprint(out, (hit.dstrand != 0) ? '-' : '+');
	}
	else if (parameters.symtype == SymbolType::blastx)
	{
	  fprint(out, ' ');
	  fprint(out, (hit.qstrand != 0) ? '-' : '+');
	  fprint_integer(out, hit.qframe+1);
	}
	else if (parameters.symtype == SymbolType::tblastn)
	{
	  fprint(out, ' ');
	  fprint(out, (hit.dstrand != 0) ? '-' : '+');
	  fprint_integer(out, hit.dframe+1);
	}
	else if (parameters.symtype == SymbolType::tblastx)
	{
	  fprint(out, ' ');
	  fprint(out, (hit.qstrand != 0) ? '-' : '+');
	  fprint_integer(out, hit.qframe + 1);
	  fprint(out, '/');
	  fprint(out, (hit.dstrand != 0) ? '-' : '+');
	  fprint_integer(out, hit.dframe + 1);
	}

	if (statistics.available != 0)
	{
	  auto const bits = static_cast<long>(floor(bit_score_of(score) + 0.5));
	  auto const expect_value = expect_value_of(score);
		
	  fprint(out, ' ');
	  fprint_integer(out, bits, score_width);
		
	  fprint(out, "   ");
		
	  hits_show_expect(expect_value);
	}
	else
	{
	  fprint(out, ' ');
	  fprint_integer(out, score, score_width);
	}

	fprint(out, '\n');
      }

      for(long i=0; i<shown.alignments; i++)
      {
	fprint(out, "\n");
	HeaderLayout layout;
	layout.show_gis = parameters.show_gis;
	layout.indent = alignment_header_indent;
	layout.linelen = alignment_header_width;
	layout.maxdeflines = LONG_MAX;
	db_showheader(t, make_view(hit_entry(i).header_address), layout);
	if ((parameters.symtype == SymbolType::tblastn) || (parameters.symtype == SymbolType::tblastx))
	{
	  fprint(out, "          Length = ");
	  fprint_integer(out, hit_entry(i).dlennt);
	  fprint(out, '\n');
	}
	else
	{
	  fprint(out, "          Length = ");
	  fprint_integer(out, hit_entry(i).dlen);
	  fprint(out, '\n');
	}
	fprint(out, "\n");
	      
	auto const score = hit_entry(i).score;

	if (statistics.available != 0)
	{
	  auto const bits = bit_score_of(score);
	  auto const expect_value = expect_value_of(score);
		
	  fprintf(out, " Score = %.1lf bits (%ld), Expect = ", bits, score);
	  hits_show_expect(expect_value);
	}
	else
	{
	  fprint(out, " Score = ");
	  fprint_integer(out, score);
	}

	fprint(out, '\n');


	auto const hit = aligned_hit(parameters, i);
	auto const counts = count_align(hit);
	      
	fprint(out, " Identities = ");
	fprint_integer(out, counts.identities);
	fprint(out, '/');
	fprint_integer(out, counts.aligned);
	fprint(out, " (");
	fprint_integer(out, whole_percentage(counts.identities, counts.aligned));
	fprint(out, "%)");
	if (parameters.symtype > SymbolType::blastn)
	{
	  fprint(out, ", Positives = ");
	  fprint_integer(out, counts.positives);
	  fprint(out, '/');
	  fprint_integer(out, counts.aligned);
	  fprint(out, " (");
	  fprint_integer(out, whole_percentage(counts.positives, counts.aligned));
	  fprint(out, "%)");
	}
	if (counts.indels != 0)
	{
	  fprint(out, ", Gaps = ");
	  fprint_integer(out, counts.indels);
	  fprint(out, '/');
	  fprint_integer(out, counts.aligned);
	  fprint(out, " (");
	  fprint_integer(out, whole_percentage(counts.indels, counts.aligned));
	  fprint(out, "%)");
	}
	fprint(out, "\n");

	if (parameters.symtype == SymbolType::blastn)
	{
	  fprint(out, " Strand = ");
	  fprint(out, as_c_string((hit_entry(i).dstrand != 0) ? "Plus / Minus" : "Plus / Plus"));
	  fprint(out, '\n');
	}
	else if (parameters.symtype == SymbolType::blastx)
	{
	  fprint(out, " Frame = ");
	  fprint(out, (hit_entry(i).qstrand != 0) ? '-':'+');
	  fprint_integer(out, hit_entry(i).qframe+1);
	  fprint(out, '\n');
	}
	else if (parameters.symtype == SymbolType::tblastn)
	{
	  fprint(out, " Frame = ");
	  fprint(out, (hit_entry(i).dstrand != 0) ? '-':'+');
	  fprint_integer(out, hit_entry(i).dframe+1);
	  fprint(out, '\n');
	}
	else if (parameters.symtype == SymbolType::tblastx)
	{
	  fprint(out, " Frame = ");
	  fprint(out, (hit_entry(i).qstrand != 0) ? '-' : '+');
	  fprint_integer(out, hit_entry(i).qframe + 1);
	  fprint(out, " / ");
	  fprint(out, (hit_entry(i).dstrand != 0) ? '-' : '+');
	  fprint_integer(out, hit_entry(i).dframe + 1);
	  fprint(out, '\n');
	}

	show_align(hit);
	fprint(out, "\n");
      }
	  
    }
    //      fprintf(out, "\n");
}

}  // anonymous namespace

auto hits_show_begin(OutputFormat view) -> void
{
  if (view==OutputFormat::plain)
    {
      fprint(out, as_c_string(swipe_name_and_version));
      fprint(out, "\n\n");
      fprint(out, as_c_string("Reference: T. Rognes (2011) Faster Smith-Waterman database searches\nwith inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221."));
      fprint(out, "\n\n");
    }
  else if (view==OutputFormat::xml)
    {
      // one root element around the results of all queries (KI-26)
      fprint(out, "<?xml version=\"1.0\"?>\n");
      fprint(out, "<results>\n");
    }
  else if (view==OutputFormat::paralign_xml)
    {
      constexpr char const * url1 = "http://www.w3.org/2001/XMLSchema-instance";
      constexpr char const * url2 = "http://www.paralign.org/ParalignXML.xsd";

      fprint(out, "<?xml version=\"1.0\"?>\n");
      fprint(out, "<ParalignXML xmlns:xsi=\"");
      fprint(out, as_c_string(url1));
      fprint(out, "\" xsi:noNamespaceSchemaLocation=\"");
      fprint(out, as_c_string(url2));
      fprint(out, "\">\n");
      fprint(out, "\t<programInformation>\n");
      fprint(out, "\t\t<programName>swipe</programName>\n");
      fprint(out, "\t\t<programVersion>");
      fprint(out, as_c_string(swipe_name_and_version));
      fprint(out, "</programVersion>\n");
      fprint(out, "\t\t<programDescription>Smith-Waterman database searches with inter-sequence SIMD parallelisation</programDescription>\n");
      fprint(out, "\t\t<articleReferences>\n");
      fprint(out, "\t\t\t<reference>T. Rognes (2011) Faster Smith-Waterman database searches with inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.</reference>\n");
      fprint(out, "\t\t</articleReferences>\n");
      fprint(out, "\t\t<license>SWIPE is available under the GNU Affero General Public License, version 3</license>\n");
      fprint(out, "\t</programInformation>\n");
    }
}

auto hits_show_end(OutputFormat view) -> void
{
  if (view==OutputFormat::xml)
  {
    fprint(out, "</results>\n");
  }
  else if (view==OutputFormat::paralign_xml)
  {
    fprint(out, "</ParalignXML>\n");
  }
}

auto hits_show(Parameters const & parameters) -> void
{
  OutputFormat const view = parameters.view;

  // compute number of hits and alignments to actually show

  long const count = hit_list.count;
  ShownHits const shown {std::min(count, hit_list.descriptions),
                         std::min(count, hit_list.alignments)};

  auto const t = db_thread_create();

  if(view == OutputFormat::plain)
  {
    hits_show_plain(parameters, shown, *t);
  }
  else if (view==OutputFormat::xml)
  {
    hits_show_xml(parameters, shown, *t);
  }
  else if ((view==OutputFormat::tabular)||(view==OutputFormat::tabular_with_comments))
  {
    auto const comments = (view == OutputFormat::tabular_with_comments) ?
      TabularComments::with : TabularComments::without;
    hits_show_tsv(parameters, shown, comments, *t);
  }
  else if (view==OutputFormat::paralign_xml)
  {
    hits_show_xml_paralign(parameters, shown, *t);
  }
}

