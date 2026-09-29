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
#include "print_view.h"  // as_c_string, fprint, fprint_integer, fprint_spaces
#include <algorithm>  // std::min, std::sort
#include <array>
#include <cassert>
#include <cctype>  // std::isspace
#include <cmath>  // std::isnan
#include <cinttypes>  // PRId64
#include <cstddef>  // std::size_t
#include <cstdint>  // std::int64_t, INT64_C
#include <cstdlib>  // std::strtol
#include <cstring>  // std::strncmp
#include <initializer_list>
#include <iterator>  // std::next
#include <limits>
#include <mutex>  // std::mutex, std::lock_guard
#include <numeric>  // std::iota
#include <string>
#include <utility>  // std::move
#include <vector>

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

long keephits;
long scorethreshold;
long upperscorethreshold;
int hits_count;
long init_threshold;
long obvious;

long opt_descriptions;
long opt_alignments;

}  // anonymous namespace

/* parameters for bit scores and expect values */

namespace {

long stats_available = 0;

double alpha;
double beta;
double lambda;
double K;
double H;
double Kmn = 0;

/* ungapped statistical parameters, only shown with -m 99 (KI-31) */
double ungapped_lambda = 0;
double ungapped_K = 0;
double ungapped_H = 0;

}  // anonymous namespace

/* gap penalties of the ungapped rows of the NCBI score matrix tables
   (INT2_MAX, see blastkar_partial.c) */
constexpr long ungapped_penalty = 32767;

namespace {

double logK;
double lambda_d_log2;
double logK_d_log2;

// E-value and bit score of a raw score, and a percentage (the same
// expressions as at their former call sites: identical results)
auto expect_value_of(long const score) -> double
{
  return Kmn * exp(- lambda * static_cast<double>(score));
}

auto bit_score_of(long const score) -> double
{
  return (lambda_d_log2 * static_cast<double>(score)) - logK_d_log2;
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

Buffer<hits_entry> hits_list;

// the hit of rank i in the list (a long, as the hit counts)
auto hit_entry(long const i) -> struct hits_entry &
{
  assert((i >= 0) and (static_cast<std::size_t>(i) < hits_list.size()));
  return hits_list[static_cast<std::size_t>(i)];
}


std::mutex hitsmutex;

auto hits_compare(void const * a, void const * b) -> int
{
  auto const index_a = *static_cast<long const *>(a);
  auto const index_b = *static_cast<long const *>(b);
  struct hits_entry const * ap = &hit_entry(index_a);
  struct hits_entry const * bp = &hit_entry(index_b);
  
  if ( static_cast<int>(index_a >= opt_alignments) < static_cast<int>(index_b >= opt_alignments) )
  {
    return -1;
  }
  if ( static_cast<int>(index_a >= opt_alignments) > static_cast<int>(index_b >= opt_alignments) )
  {
    return +1;
  }
  if (ap->qstrand < bp->qstrand)
  {
    return -1;
  }
  if (ap->qstrand > bp->qstrand)
  {
    return +1;
  }
  if (ap->qframe < bp->qframe)
  {
    return -1;
  }
  if (ap->qframe > bp->qframe)
  {
    return +1;
  }
  if (ap->seqno < bp->seqno)
  {
    return -1;
  }
  if (ap->seqno > bp->seqno)
  {
    return +1;
  }
  if (ap->dstrand < bp->dstrand)
  {
    return -1;
  }
  if (ap->dstrand > bp->dstrand)
  {
    return +1;
  }
  if (ap->dframe < bp->dframe)
  {
    return -1;
  }
  if (ap->dframe > bp->dframe)
  {
    return +1;
  }

  return 0;
}

}  // anonymous namespace

auto hits_sort() -> Buffer<long>
{
  Buffer<long> hits_sorted(static_cast<std::size_t>(hits_count));
  std::iota(hits_sorted.begin(), hits_sorted.end(), 0L);
  std::sort(hits_sorted.begin(), hits_sorted.end(),
            [](long const lhs, long const rhs) -> bool {
              return hits_compare(&lhs, &rhs) < 0;
            });
  return hits_sorted;
}


auto hits_enter(long seqno, long score, HitStrands const & strands) -> void
{
  // show_progress();

  //  fprintf(out, "Entering score %u from sequence no %u.\n", score, seqno);
  
  // find correct place

  std::lock_guard<std::mutex> const lock(hitsmutex);

  if (score > upperscorethreshold)
  {
    obvious++;
  }

  if (score >= init_threshold)
  {
    totalhits++;
  }

  if ((score < scorethreshold) || (score > upperscorethreshold))
  {
    return;
  }

  long place = hits_count;

  while ((place > 0) && ((score > hit_entry(place-1).score) ||
			 ((score == hit_entry(place-1).score) &&
			  (seqno > hit_entry(place-1).seqno))))
  {
    place--;
  }

  // move entries down
  
  long const move = (hits_count < keephits ? hits_count : keephits - 1) - place;

  //  fprintf(out, "Inserting at place %d, moving %d.\n", place, move);

  for (long j = move; j > 0; j--)
  {
    hit_entry(place + j) = std::move(hit_entry(place + j - 1));
  }

  // fill new entry

  if (place < keephits)
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
    if (hits_count < keephits)
    {
      hits_count++;
    }
  }
  
  // no hit is kept with -v 0 -b 0: the list is empty (KI-10)
  if ((keephits > 0) and (hits_count == keephits))
  {
    scorethreshold = hit_entry(keephits - 1).score;
  }

}

auto hits_getcount() -> long
{
  return hits_count;
}

auto hits_gethit(long i, long * seqno, long * score, 
		 long * qstrand, long * qframe,
		 long * dstrand, long * dframe) -> void
{
  struct hits_entry const * h = &hit_entry(i);
  *seqno = h->seqno;
  *score = h->score;
  *qstrand = h->qstrand;
  *qframe = h->qframe;
  *dstrand = h->dstrand;
  *dframe = h->dframe;
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
  int const show_nostats = static_cast<int>(parameters.view == OutputFormat::plain);

  opt_descriptions = descriptions;
  opt_alignments = max_alignments;
  keephits = descriptions > max_alignments ? descriptions : max_alignments;
  
  std::int64_t maxhits = db_getseqcount_masked();
  if (parameters.symtype == SymbolType::blastn)
    {
      if (parameters.querystrands == QueryStrands::both)
      {
	maxhits *= 2;
      }
    }
  else if (parameters.symtype == SymbolType::blastx)
    {
      if (parameters.querystrands == QueryStrands::both)
      {
	maxhits *= 6;
      }
      else
      {
	maxhits *= 3;
      }
    }
  else if (parameters.symtype == SymbolType::tblastn)
    {
      maxhits *= 6;
    }
  else if (parameters.symtype == SymbolType::tblastx)
    {
      if (parameters.querystrands == QueryStrands::both)
      {
	maxhits *= 36;
      }
      else
      {
	maxhits *= 18;
      }
    }

  keephits = static_cast<long>(std::min<std::int64_t>(keephits, maxhits));

  obvious = 0;
  hits_count = 0;
  hits_list.clear();
  hits_list.resize(static_cast<std::size_t>(keephits));

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

  stats_available = 0;

  std::int64_t m = 0;
  std::int64_t n = 0;
  int lenadj = 0;
  if (parameters.symtype == SymbolType::blastn)
  {
    if (stats_getparams_nt(parameters.matchscore,
			   parameters.mismatchscore,
			   parameters.gapopen,
			   parameters.gapextend,
			   & lambda,
			   & K,
			   & H,
			   & alpha,
			   & beta) != 0)
    {
      stats_available = 1;

      /*
      fprintf(out, "Params: lambda=%6.3g K=%6.3g H=%6.3g alpha=%6.3g beta=%6.3g\n",
	      lambda, K, H, alpha, beta);
      */

      logK = log(K);
      lambda_d_log2 = lambda / log(2.0);
      logK_d_log2 = logK / log(2.0);
      
      lenadj = 0;
      
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

      BlastComputeLengthAdjustment(K,
				   logK,
				   alpha / lambda,
				   beta,
				   static_cast<Int4>(qlen),
				   dlen,
				   static_cast<Int4>(seqcount),
				   & lenadj);
    
      //      fprintf(out, "lenadj: %d\n", lenadj);

      m = qlen - lenadj;

      if (parameters.effdbsize > 0)
      {
	n = parameters.effdbsize;
      }
      else
      {
	n = effective_db_length(dlen, seqcount, lenadj);
      }

      Kmn = K * static_cast<double>(m) * static_cast<double>(n);
    }
  }
  else if (parameters.symtype < SymbolType::sound)
  {
    if (parameters.symtype == SymbolType::tblastx)
    {
      stats_available = stats_getparams(parameters.matrixname,
					32767,
					32767,
					& lambda,
					& K,
					& H,
					& alpha,
					& beta);
    }
    else
    {
      stats_available = stats_getparams(parameters.matrixname,
					parameters.gapopen,
					parameters.gapextend,
					& lambda,
					& K,
					& H,
					& alpha,
					& beta);
    }


    if (stats_available != 0)
    {
      
      logK = log(K);
      lambda_d_log2 = lambda / log(2.0);
      logK_d_log2 = logK / log(2.0);
      
      lenadj = 0;
      
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

      BlastComputeLengthAdjustment(K,
				   logK,
				   alpha / lambda,
				   beta,
				   static_cast<Int4>(qlen),
				   dlen,
				   static_cast<Int4>(seqcount),
				   & lenadj);

      m = qlen - lenadj;

      if (parameters.effdbsize > 0)
      {
	n = parameters.effdbsize;
      }
      else
      {
	n = effective_db_length(dlen, seqcount, lenadj);
      }

      Kmn = K * static_cast<double>(m) * static_cast<double>(n);
    }
  }

  /* ungapped statistical parameters (-m 99): the (0, 0) rows of the
     nucleotide tables, the ungapped rows of the matrix tables; the
     gapped values when there are none (KI-31) */
  ungapped_lambda = lambda;
  ungapped_K = K;
  ungapped_H = H;
  double ungapped_alpha = 0;
  double ungapped_beta = 0;
  if (parameters.symtype == SymbolType::blastn)
  {
    stats_getparams_nt(parameters.matchscore, parameters.mismatchscore, 0, 0,
                       & ungapped_lambda, & ungapped_K, & ungapped_H,
                       & ungapped_alpha, & ungapped_beta);
  }
  else if (parameters.symtype < SymbolType::sound)
  {
    stats_getparams(parameters.matrixname, ungapped_penalty, ungapped_penalty,
		    &ungapped_lambda, &ungapped_K, &ungapped_H,
		    &ungapped_alpha, &ungapped_beta);
  }

  scorethreshold = minscore;
  upperscorethreshold = maxscore;
  
  if (stats_available != 0)
  {
    long const minscore_expect = threshold_to_long(ceil(- log(max_expect / Kmn) / lambda));
    if (minscore_expect > minscore)
    {
      scorethreshold = minscore_expect;
    }

    if (min_expect > 0.0)
    {
      long const maxscore_expect = threshold_to_long(floor(- log(min_expect / Kmn) / lambda));
      if (maxscore_expect < maxscore)
      {
	upperscorethreshold = maxscore_expect;
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

  init_threshold = scorethreshold;

  //  fprintf(out, "scorethreshold: %ld\n", scorethreshold);
}

auto hits_empty() -> void
{
  for (long i=0; i<hits_count; i++)
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
  hits_list = Buffer<hits_entry>();
}

auto hits_align(Parameters const & parameters, struct db_thread_s * t, long i) -> void
{
  long ntlen = 0;

  struct hits_entry * h = &hit_entry(i);

  db_mapheaders(t, h->seqno, h->seqno);

  View<char> const header = db_getheader(t, h->seqno);
  h->header_address.assign(header.begin(), header.end());

  // the sequence length is needed for every hit shown (-m 7 <len>,
  // KI-37), the sequence itself only for hits with an alignment
  db_mapsequences(t, h->seqno, h->seqno);

  View<char> const sequence = db_getsequence(t, h->seqno, h->dstrand,
					     h->dframe, & ntlen, 0);
  h->dlen = static_cast<long>(sequence.size());
  h->dlennt = ntlen;

  if (i < opt_alignments)
  {
    h->dseq.assign(sequence.begin(), sequence.end());
    
    char * qseq = nullptr;
    long qlen = 0;
    if (parameters.symtype == SymbolType::blastn)
    {
      qseq = query.nt[0].seq;
      qlen = query.nt[0].len;
    }
    else
    {
      qseq = query.aa[(3*h->qstrand) + h->qframe].seq;
      qlen = query.aa[(3*h->qstrand) + h->qframe].len;
    }

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

    align(qseq,
	  h->dseq.data(),
	  qlen,
	  h->dlen,
	  score_matrix_63,
	  parameters.gapopen,
	  parameters.gapextend,
	  & h->align_q_start,
	  & h->align_d_start,
	  & h->align_q_end,
	  & h->align_d_end,
	  h->alignment,
	  & h->score_align);
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
  struct hits_entry const & entry = hit_entry(i);
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
    hit.q_seq = query.nt[hit.q_strand].seq;
    hit.q_len = query.nt[hit.q_strand].len;
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
    hit.q_seq = query.aa[(3*hit.q_strand)+hit.q_frame].seq;
    hit.q_len = query.aa[(3*hit.q_strand)+hit.q_frame].len;
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
  long maxpos = maxqpos > maxdpos ? maxqpos : maxdpos;
  hit.poswidth = 1;
  while (maxpos > 9)
  {
    maxpos /= 10;
    hit.poswidth++;
  }

  return hit;
}

// show_align(): the alignment of a hit, in lines of ALIGNLEN columns
class AlignmentLines
{
public:
  explicit AlignmentLines(AlignedHit const & aligned) :
    hit(aligned), q_pos(aligned.q_align_start), d_pos(aligned.d_align_start) {}

  auto putalignop(char c, long len) -> void;

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

auto AlignmentLines::putalignop(char c, long len) -> void
{

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
	  (score_matrix_63[(32*qs)+ds] > 0 ? '+' : ' ');
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
struct AlignmentOperation
{
  char op = 0;
  long len = 0;
};

auto next_operation(char const * & cursor) -> AlignmentOperation
{
  AlignmentOperation operation;
  operation.op = *cursor;
  cursor = std::next(cursor);
  char * end = nullptr;
  operation.len = std::strtol(cursor, & end, 10);
  cursor = end;
  return operation;
}

auto show_align(AlignedHit const & hit) -> void
{
  AlignmentLines lines(hit);
  
  char const * p = hit.alignment;
  char const * e = hit.alignment + strlen(hit.alignment);
  
  while(p < e)
  {
    AlignmentOperation const operation = next_operation(p);
    lines.putalignop(operation.op, operation.len);
  }
  
  lines.putalignop(0, 1);
}

auto whole_align(AlignedHit const & hit,
		 long * identities,
		 long * positives,
		 long * indels,
		 long * aligned,
		 long * gaps,
		 std::string & qline,
		 std::string & aline,
		 std::string & dline) -> void
{

  long al = 0;
  char const * p = hit.alignment;
  while((*p) != 0)
  {
    al += next_operation(p).len;
  }

  for (auto * line : {&qline, &aline, &dline})
  {
    line->clear();
    line->reserve(static_cast<std::size_t>(al));
  }

  *identities = 0;
  *positives = 0;
  *indels = 0;
  *gaps = 0;
  *aligned = 0;
  
  long q_pos = hit.q_align_start;
  long d_pos = hit.d_align_start;
  
  p = hit.alignment;

  while((*p) != 0)
  {
    AlignmentOperation const operation = next_operation(p);
    char const op = operation.op;
    long const len = operation.len;
    
    *aligned += len;
    if (op == 'D')
    {
      for(long j=0; j<len; j++)
      {
	char const qs = hit.q_seq[q_pos++];
	qline += hit.sym[static_cast<int>(qs)];
	aline += ' ';
	dline += '-';
      }
      *gaps += 1;
      *indels += len;
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
      *gaps += 1;
      *indels += len;
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
	  (*identities)++;
	  (*positives)++;
	}
	else if (score_matrix_63[(32*qs)+ds] > 0)
	{
	  aline += '+';
	  (*positives)++;
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

}

auto count_align(AlignedHit const & hit,
		 long * identities,
		 long * positives,
		 long * indels,
		 long * aligned,
		 long * gaps) -> void
{
  *identities = 0;
  *positives = 0;
  *indels = 0;
  *gaps = 0;
  *aligned = 0;
  
  long q_pos = hit.q_align_start;
  long d_pos = hit.d_align_start;
  
  char const * p = hit.alignment;
  char const * e = hit.alignment + strlen(hit.alignment);

  while(p < e)
  {
    AlignmentOperation const operation = next_operation(p);
    char const op = operation.op;
    long const len = operation.len;
    
    *aligned += len;
    if (op == 'D')
    {
      *gaps += 1;
      *indels += len;
      q_pos += len;
    }
    else if (op == 'I')
    {
      *gaps += 1;
      *indels += len;
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
	  (*identities)++;
	  (*positives)++;
	}
	else if (score_matrix_63[(32 * qs) + ds] > 0)
	{
	  (*positives)++;
	}
      }
    }
  }
}

auto hits_show_expect(double expect_value) -> void
{
  char temp[10];
  if (expect_value < 1e-180)
  {
    fprint(out, "0.0  ");
  }
  else if (expect_value < 9.5e-100)
  {
    snprintf(temp, sizeof(temp), "%-6.0e", expect_value);
    fputs(temp+1, out);
  }
  else if (expect_value < 0.00095)
  {
    fprintf(out, "%-5.0e", expect_value);
  }
  else if (expect_value < 0.0995)
  {
    fprintf(out, "%-5.3f", expect_value);
  }
  else if (expect_value < 0.95)
  {
    fprintf(out, "%-5.2f", expect_value);
  }
  else if (expect_value < 9.5)
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
auto xml_print(char const * const text,
               std::size_t const max_length = std::numeric_limits<std::size_t>::max()) noexcept -> void
{
  for (std::size_t i = 0; (i < max_length) and (text[i] != '\0'); ++i)
  {
    xml_putc(text[i]);
  }
}

auto make_anchor(char * anchor, std::size_t const size, SymbolType symbol_type, long query_index, long i) -> void
{
  switch(symbol_type)
  {
  case SymbolType::blastn:
    // blastn: the strand of a hit is stored as its database strand
    // (KI-29)
    snprintf(anchor, size, "%ld_%ld__%c__+",
	     query_index,
	     hit_entry(i).seqno,
	     (hit_entry(i).dstrand != 0) ? '-' : '+');
    break;
  case SymbolType::blastx:
    snprintf(anchor, size, "%ld_%ld_%ld_%c__",
	     query_index,
	     hit_entry(i).seqno,
	     hit_entry(i).qframe+1,
	     (hit_entry(i).qstrand != 0) ? '-' : '+');
    break;
  case SymbolType::tblastn:
    snprintf(anchor, size, "%ld_%ld___%ld_%c",
	     query_index,
	     hit_entry(i).seqno,
	     hit_entry(i).dframe+1,
	     (hit_entry(i).dstrand != 0) ? '-' : '+');
    break;
  case SymbolType::tblastx:
    snprintf(anchor, size, "%ld_%ld_%ld_%c_%ld_%c",
	     query_index,
	     hit_entry(i).seqno,
	     hit_entry(i).qframe+1,
	     (hit_entry(i).qstrand != 0) ? '-' : '+',
	     hit_entry(i).dframe+1,
	     (hit_entry(i).dstrand != 0) ? '-' : '+');
    break;
  default:
    snprintf(anchor, size, "%ld_%ld____",
	     query_index,
	     hit_entry(i).seqno);
    break;
  }
}

auto hits_defline_split(char * defline, 
			long * gi,
			char ** link, std::size_t * linklen, 
			char ** rest) -> void
{
  char * p = defline;

  *gi = 0;  // no gi (KI-42: it kept the gi of the previous defline)
  *link = nullptr;
  *linklen = 0;
  *rest = nullptr;
  
  // "gi|" and a number, as the header parser writes them (set_id(),
  // asnparse.cc)
  constexpr std::size_t gi_prefix_length = 3;
  if (std::strncmp(p, "gi|", gi_prefix_length) == 0)
  {
    char * const number = std::next(p, gi_prefix_length);
    char * end = nullptr;
    long const value = std::strtol(number, & end, 10);
    if (end != number)
    {
      *gi = value;
      p = end;
    }
  }

  if (*p == '|')
  {
    p++;
  }

  char * r = strchr(p, ' ');
  if (r != nullptr)
  {
    *linklen = static_cast<std::size_t>(r - p);
    *link = p;
    *rest = r+1;
  }
  else
  {
    * rest = p;
  }
}

auto hits_show_xml_paralign(Parameters const & parameters,
			    long showalignments,
			    long showhits,
			    struct db_thread_s const * t) -> void
{
  /* ParAlign XML */
  
  fprint(out, "\t<paralignOutput>\n");
  
  char const * qseqtypedescr = nullptr;
  struct sequence q;
  if ((query.symtype == SymbolType::blastp) || (query.symtype == SymbolType::tblastn))
  {
    qseqtypedescr = "Amino Acid";
    q = query.aa[0];
  }
  else if (query.symtype == SymbolType::sound)
  {
    /* sound queries are stored as amino acid queries (KI-30) */
    qseqtypedescr = "Sound";
    q = query.aa[0];
  }
  else
  {
    qseqtypedescr = "Nucleotide";
    q = query.nt[0];
  }
  
  fprint(out, "\t\t<queryInformation>\n");
  fprint(out, "\t\t\t<queryFilename>");
  xml_print(parameters.queryname);
  fprint(out, "</queryFilename>\n");
  fprint(out, "\t\t\t<querySequencetype>");
  fprint(out, as_c_string(qseqtypedescr));
  fprint(out, "</querySequencetype>\n");
  fprint(out, "\t\t\t<queryDescription>");
  xml_print(query.description.c_str());
  fprint(out, "</queryDescription>\n");
  fprint(out, "\t\t\t<queryLength>");
  fprint_integer(out, q.len);
  fprint(out, "</queryLength>\n");
  fprint(out, "\t\t\t<querySequence>");
  for (int i = 0; i < q.len; i++)
  {
    fprint(out, query.sym[static_cast<int>(q.seq[i])]);
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
  xml_print(parameters.databasename);
  fprint(out, "</databaseFilename>\n");
  fprint(out, "\t\t\t<databaseSequencetype>");
  fprint(out, as_c_string(dbseqtypedescr));
  fprint(out, "</databaseSequencetype>\n");
  fprint(out, "\t\t\t<databaseDescription>");
  xml_print(db_gettitle());
  fprint(out, "</databaseDescription>\n");
  fprint(out, "\t\t\t<databaseVersion>");
  fprint_integer(out, db_getversion());
  fprint(out, "</databaseVersion>\n");
  fprint(out, "\t\t\t<databaseDate>");
  xml_print(db_gettime());
  fprint(out, "</databaseDate>\n");
  fprintf(out, "\t\t\t<residueCount>%" PRId64 "</residueCount>\n", db_getsymcount_masked());
  fprintf(out, "\t\t\t<sequenceCount>%" PRId64 "</sequenceCount>\n", db_getseqcount_masked());
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
    xml_print(parameters.matrixname);
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
  fprintf(out, "\t\t\t\t\t<ungappedLambda>%.4g</ungappedLambda>\n", ungapped_lambda);
  fprintf(out, "\t\t\t\t\t<ungappedKappa>%.4g</ungappedKappa>\n", ungapped_K);
  fprintf(out, "\t\t\t\t\t<ungappedEta>%.4g</ungappedEta>\n", ungapped_H);
  fprint(out, "\t\t\t\t</ungapped>\n");
  fprint(out, "\t\t\t\t<gapped>\n");
  fprintf(out, "\t\t\t\t\t<gappedLambda>%.4g</gappedLambda>\n", lambda);
  fprintf(out, "\t\t\t\t\t<gappedKappa>%.4g</gappedKappa>\n", K);
  fprintf(out, "\t\t\t\t\t<gappedEta>%.4g</gappedEta>\n", H);
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
  fprint(out, as_c_string(ti.starttime.data()));
  fprint(out, "</searchStarted>\n");
  fprint(out, "\t\t\t\t<searchCompleted>");
  fprint(out, as_c_string(ti.endtime.data()));
  fprint(out, "</searchCompleted>\n");
  fprintf(out, "\t\t\t\t<searchElapsedTime>%.2fs</searchElapsedTime>\n", ti.elapsed);
  if (ti.elapsed > 0.0)
  {
    fprintf(out, "\t\t\t\t<searchSpeed>%.3f GCUPS</searchSpeed>\n", ti.speed / 1e9);
  }
  else
  {
    fprint(out, "\t\t\t\t<searchSpeed>n/a</searchSpeed>\n");
  }
  fprint(out, "\t\t\t\t<searchSWAlignments>\n");
  fprint(out, "\t\t\t\t\t<SWAbsolute>");
  fprint_integer(out, compute7);
  fprint(out, "</SWAbsolute>\n");
  fprint(out, "\t\t\t\t\t<SWPercent>100</SWPercent>\n");
  fprint(out, "\t\t\t\t</searchSWAlignments>\n");
  fprint(out, "\t\t\t</searchInformation>\n");

  fprint(out, "\t\t<resultInformation>\n");
  fprint(out, "\t\t\t<resultHits>\n");
  fprint(out, "\t\t\t\t<totalCount>");
  fprint_integer(out, totalhits);
  fprint(out, "</totalCount>\n");
  fprint(out, "\t\t\t\t<obviousCount>");
  fprint_integer(out, obvious);
  fprint(out, "</obviousCount>\n");
  fprint(out, "\t\t\t\t<shownCount>");
  fprint_integer(out, showhits);
  fprint(out, "</shownCount>\n");
  fprint(out, "\t\t\t</resultHits>\n");
  fprint(out, "\t\t\t<alignmentCount>");
  fprint_integer(out, showalignments);
  fprint(out, "</alignmentCount>\n");
  fprint(out, "\t\t</resultInformation>\n");
  
  fprint(out, "\t\t<shortVersionHits>\n");
  
  for(long i=0; i<showhits; i++)
  {
    long const score = hit_entry(i).score;
    double const e = expect_value_of(score);

    char anchor[200];
    make_anchor(anchor, 200, query.symtype, queryno, i);

    long deflines = 0;
    std::vector<std::string> deflinetable;
    long gi = 0;
    char * link = nullptr;
    char * title = nullptr;
    std::size_t linklen = 0;
    db_parse_header(t, make_view(hit_entry(i).header_address),
		    1, & deflines, & deflinetable);
    hits_defline_split(&deflinetable[0][0], 
		       & gi,
		       & link, & linklen,
		       & title);

    fprint(out, "\t\t\t<shortVersionHit>\n");
    fprint(out, "\t\t\t\t<shortVersionAnchor>");
    fprint(out, as_c_string(anchor));
    fprint(out, "</shortVersionAnchor>\n");
    if (gi != 0)
      {
    fprint(out, "\t\t\t\t<shortVersionLink>\n");
    fprint(out, "\t\t\t\t\t<shortVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=");
    fprint(out, as_c_string(ncbidb));
    fprint(out, "&amp;list_uids=");
    fprint_integer(out, gi);
    fprint(out, "&amp;dopt=");
    fprint(out, as_c_string(ncbiopt));
    fprint(out, "</shortVersionLinkDestination>\n");
    fprint(out, "\t\t\t\t\t<shortVersionLinkText>gi|");
    fprint_integer(out, gi);
    fprint(out, "</shortVersionLinkText>\n");
    fprint(out, "\t\t\t\t</shortVersionLink>\n");
      }
    fprint(out, "\t\t\t\t<shortVersionLink>\n");
    fprint(out, "\t\t\t\t\t<shortVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&amp;db=");
    fprint(out, as_c_string(ncbidb));
    fprint(out, "&amp;term=");
    xml_print(link, linklen);
    fprint(out, "&amp;doptcmdl=");
    fprint(out, as_c_string(ncbiopt));
    fprint(out, "</shortVersionLinkDestination>\n");
    fprint(out, "\t\t\t\t\t<shortVersionLinkText>");
    xml_print(link, linklen);
    fprint(out, "</shortVersionLinkText>\n");
    fprint(out, "\t\t\t\t</shortVersionLink>\n");
    fprint(out, "\t\t\t\t<shortVersionName>");
    xml_print(title, 35);
    fprint(out, "</shortVersionName>\n");
    if (parameters.symtype == SymbolType::blastn)
    {
      fprint(out, "\t\t\t\t<shortVersionStrand>");
      fprint(out, (hit_entry(i).dstrand != 0) ? '-' : '+');
      fprint(out, "</shortVersionStrand>\n");
    }
    else if (parameters.symtype == SymbolType::blastx)
    {
      fprintf(out, "\t\t\t\t<shortVersionFrame>%c%ld</shortVersionFrame>\n", 
	     (hit_entry(i).qstrand != 0) ? '-' : '+', 
	     hit_entry(i).qframe+1);
    }
    else if (parameters.symtype == SymbolType::tblastn)
    {
      fprintf(out, "\t\t\t\t<shortVersionFrame>%c%ld</shortVersionFrame>\n", 
	     (hit_entry(i).dstrand != 0) ? '-' : '+', 
	     hit_entry(i).dframe+1);
    }
    else if (parameters.symtype == SymbolType::tblastx)
    {
      fprintf(out, "\t\t\t\t<shortVersionFrame>%c%ld/%c%ld</shortVersionFrame>\n", 
	     (hit_entry(i).qstrand != 0) ? '-' : '+', 
	     hit_entry(i).qframe+1,
	     (hit_entry(i).dstrand != 0) ? '-' : '+', 
	     hit_entry(i).dframe+1);
    }
    fprint(out, "\t\t\t\t<shortVersionScore>");
    fprint_integer(out, score);
    fprint(out, "</shortVersionScore>\n");
    fprintf(out, "\t\t\t\t<shortVersionEValue>%.2g</shortVersionEValue>\n", e);
    fprint(out, "\t\t\t</shortVersionHit>\n");

  }

  fprint(out, "\t\t</shortVersionHits>\n");

  if (showalignments != 0)
  {
    fprint(out, "\t\t<longVersionHits>\n");
    
    for(long i=0; i<showalignments; i++)
    {
      
      char anchor[200];
      make_anchor(anchor, 200, query.symtype, queryno, i);
      
      fprint(out, "\t\t\t<longVersionHit>\n");
      fprint(out, "\t\t\t\t<longVersionAnchor>");
      fprint(out, as_c_string(anchor));
      fprint(out, "</longVersionAnchor>\n");
      
      long deflines = 0;
      std::vector<std::string> deflinetable;
      long gi = 0;
      char * link = nullptr;
      char * title = nullptr;
      std::size_t linklen = 0;
      db_parse_header(t, make_view(hit_entry(i).header_address),
		      1, & deflines, & deflinetable);
      fprint(out, "\t\t\t\t<linkContainer>\n");
      
      for (int d=0; d < deflines; d++)
      {
	hits_defline_split(&deflinetable[static_cast<std::size_t>(d)][0], 
			   & gi,
			   & link, & linklen,
			   & title);
  
        if (gi != 0)
	{
          fprint(out, "\t\t\t\t\t<longVersionLink>\n");
	  fprint(out, "\t\t\t\t\t\t<longVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=");
	  fprint(out, as_c_string(ncbidb));
	  fprint(out, "&amp;list_uids=");
	  fprint_integer(out, gi);
	  fprint(out, "&amp;dopt=");
	  fprint(out, as_c_string(ncbiopt));
	  fprint(out, "</longVersionLinkDestination>\n");
	  fprint(out, "\t\t\t\t\t\t<longVersionLinkText>gi|");
	  fprint_integer(out, gi);
	  fprint(out, "</longVersionLinkText>\n");
	  fprint(out, "\t\t\t\t\t</longVersionLink>\n");
	}
      
	fprint(out, "\t\t\t\t\t<longVersionLink>\n");
	fprint(out, "\t\t\t\t\t\t<longVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&amp;db=");
	fprint(out, as_c_string(ncbidb));
	fprint(out, "&amp;term=");
	xml_print(link, linklen);
	fprint(out, "&amp;doptcmdl=");
	fprint(out, as_c_string(ncbiopt));
	fprint(out, "</longVersionLinkDestination>\n");
	fprint(out, "\t\t\t\t\t\t<longVersionLinkText>");
	xml_print(link, linklen);
	fprint(out, "</longVersionLinkText>\n");
	fprint(out, "\t\t\t\t\t</longVersionLink>\n");
      
	fprint(out, "\t\t\t\t\t<longVersionName>");
	xml_print(title);
	fprint(out, "</longVersionName>\n");
      }
        
      fprint(out, "\t\t\t\t</linkContainer>\n");
    
      long const dlen = hit_entry(i).dlen;
      long const dlennt = hit_entry(i).dlennt;

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

      long const score = hit_entry(i).score;
      double const e = expect_value_of(score);

      long identities = 0;
      long positives = 0;
      long indels = 0;
      long gaps = 0;
      long aligned = 0;
    
      std::string qline;
      std::string aline;
      std::string dline;
        
      AlignedHit const hit = aligned_hit(parameters, i);
      whole_align(hit, & identities, & positives, & indels, & aligned, & gaps,
		  qline, aline, dline);

      fprint(out, "\t\t\t\t<alignment>\n");
      fprint(out, "\t\t\t\t\t<subalignment>\n");
      fprint(out, "\t\t\t\t\t\t<longVersionScore>");
      fprint_integer(out, score);
      fprint(out, "</longVersionScore>\n");
      fprintf(out, "\t\t\t\t\t\t<longVersionEValue>%.2g</longVersionEValue>\n", e);
      fprint(out, "\t\t\t\t\t\t<identical>\n");
      fprint(out, "\t\t\t\t\t\t\t<identicalNominator>");
      fprint_integer(out, identities);
      fprint(out, "</identicalNominator>\n");
      fprint(out, "\t\t\t\t\t\t\t<identicalDenominator>");
      fprint_integer(out, aligned);
      fprint(out, "</identicalDenominator>\n");
      fprintf(out, "\t\t\t\t\t\t\t<identicalPercentage>%.1f</identicalPercentage>\n", percentage(identities, aligned));
      fprint(out, "\t\t\t\t\t\t</identical>\n");

      if (parameters.symtype != SymbolType::blastn)
      {
	fprint(out, "\t\t\t\t\t\t<positive>\n");
	fprint(out, "\t\t\t\t\t\t\t<positiveNominator>");
	fprint_integer(out, positives);
	fprint(out, "</positiveNominator>\n");
	fprint(out, "\t\t\t\t\t\t\t<positiveDenominator>");
	fprint_integer(out, aligned);
	fprint(out, "</positiveDenominator>\n");
	fprintf(out, "\t\t\t\t\t\t\t<positivePercentage>%.1f</positivePercentage>\n", percentage(positives, aligned));
	fprint(out, "\t\t\t\t\t\t</positive>\n");
      }

      fprint(out, "\t\t\t\t\t\t<indels>\n");
      fprint(out, "\t\t\t\t\t\t\t<indelsNominator>");
      fprint_integer(out, indels);
      fprint(out, "</indelsNominator>\n");
      fprint(out, "\t\t\t\t\t\t\t<indelsDenominator>");
      fprint_integer(out, aligned);
      fprint(out, "</indelsDenominator>\n");
      fprintf(out, "\t\t\t\t\t\t\t<indelsPercentage>%.1f</indelsPercentage>\n", percentage(indels, aligned));
      fprint(out, "\t\t\t\t\t\t</indels>\n");
      fprint(out, "\t\t\t\t\t\t<gaps>");
      fprint_integer(out, gaps);
      fprint(out, "</gaps>\n");
      fprint(out, "\t\t\t\t\t\t<alignmentQuery>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentQueryStart>");
      fprint_integer(out, hit.q_first);
      fprint(out, "</alignmentQueryStart>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentQueryLine>");
      fprint(out, as_c_string(qline.c_str()));
      fprint(out, "</alignmentQueryLine>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentQueryEnd>");
      fprint_integer(out, hit.q_last);
      fprint(out, "</alignmentQueryEnd>\n");
      fprint(out, "\t\t\t\t\t\t</alignmentQuery>\n");
      fprint(out, "\t\t\t\t\t\t<alignmentLine>");
      fprint(out, as_c_string(aline.c_str()));
      fprint(out, "</alignmentLine>\n");
      fprint(out, "\t\t\t\t\t\t<alignmentDatabase>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentDatabaseStart>");
      fprint_integer(out, hit.d_first);
      fprint(out, "</alignmentDatabaseStart>\n");
      fprint(out, "\t\t\t\t\t\t\t<alignmentDatabaseLine>");
      fprint(out, as_c_string(dline.c_str()));
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

auto show_description(char const *desc) -> void
{
  char const *dptr = nullptr;

  for (dptr = desc; not ends_query_id(*dptr); dptr++)
  {
    fprint(out, *dptr);
  }
}

// query id (the description up to its first whitespace character),
// escaped as XML
// (KI-27)
auto show_description_xml(char const * const desc) -> void
{
  for (auto const * dptr = desc; not ends_query_id(*dptr); ++dptr)
  {
    xml_putc(*dptr);
  }
}

auto hits_show_xml(Parameters const & parameters,
		   long show_gis,
		   long showalignments,
		   long showhits,
		   struct db_thread_s const * t) -> void
{
  /* Simple XML */
  
  fprint(out, "<result>\n");
  fprint(out, "  <general>\n");
  fprint(out, "    <hitcount>");
  fprint_integer(out, hits_count);
  fprint(out, "</hitcount>\n");
  fprint(out, "  </general>\n");
  fprint(out, "  <hits>\n");
  
  for(long i=0; i<showhits; i++)
  {
    long const seqno = hit_entry(i).seqno;
    long const score = hit_entry(i).score;
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
    show_description_xml(query.description.c_str());
    fprintf(out,"</query>\n");
    fprint(out, "      <name>");
    HeaderLayout layout;
    layout.show_gis = show_gis;
    layout.escaping = Escaping::xml;
    db_showheader(t, make_view(hit_entry(i).header_address), layout);
    fprint(out, "</name>\n");
    fprint(out, "      <len>");
    fprint_integer(out, dlen);
    fprint(out, "</len>\n");
    fprint(out, "      <score>");
    fprint_integer(out, score);
    fprint(out, "</score>\n");
    
    if (i < showalignments)
    {
      long identities = 0;
      long positives = 0;
      long gaps = 0;
      long aligned = 0;
      long indels = 0;

      std::string qline;
      std::string aline;
      std::string dline;
        
      AlignedHit const hit = aligned_hit(parameters, i);
      whole_align(hit, & identities, & positives, & indels, & aligned, & gaps,
		  qline, aline, dline);

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
      fprint(out, as_c_string(qline.c_str()));
      fprint(out, "</qseq>\n");
      fprint(out, "      <aseq>");
      fprint(out, as_c_string(aline.c_str()));
      fprint(out, "</aseq>\n");
      fprint(out, "      <dseq>");
      fprint(out, as_c_string(dline.c_str()));
      fprint(out, "</dseq>\n");

    }
    fprint(out, "    </hit>\n");
  }
  fprint(out, "  </hits>\n");
  fprint(out, "</result>\n");
}

auto hits_show_tsv(Parameters const & parameters,
		   long showalignments,
		   long showcomments,
		   struct db_thread_s const * t) -> void
{
  char ref[] = "Reference: T. Rognes (2011) Faster Smith-Waterman database searches with inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.";
  
  if (showcomments != 0)
    {
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
      if (stats_available != 0)
      {
	fprint(out, "# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n");
      }
      else
      {
	fprint(out, "# Fields: Query id, Subject id, % identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, score\n");
      }
    }

  for(long i=0; i<showalignments; i++)
  {
    show_description(query.description.c_str());
    fprint(out, '\t');
    HeaderLayout layout;
    layout.show_gis = 1;
    layout.text = DeflineText::identifier;
    db_showheader(t, make_view(hit_entry(i).header_address), layout);
    
    long identities = 0;
    long positives = 0;
    long gaps = 0;
    long aligned = 0;
    long indels = 0;
    
    AlignedHit const hit = aligned_hit(parameters, i);
    count_align(hit, & identities, & positives, & indels, & aligned, & gaps);
    
    long const score = hit_entry(i).score;
    
    fprintf(out, "\t%.2f\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", 
	    percentage(identities, aligned),
	    aligned,
	    aligned - identities - indels,
	    gaps,
	    hit.q_first,
	    hit.q_last,
	    hit.d_first,
	    hit.d_last);
    
    if (stats_available != 0)
    {
      double const expect_value = expect_value_of(score);
      fprintf(out, "\t%.2g", expect_value);
      double const bits = bit_score_of(score);
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
		     long show_gis,
		     long showalignments,
		     long showhits,
		     struct db_thread_s const * t) -> void
{
    if (hits_count == 0)
    {
      fprint(out, "\nNo hits.\n");
    }
    else
    {
      if (stats_available != 0)
      {
	fprint(out, "                                                                 Score    E\n");
	fprint(out, "Sequences producing significant alignments:                      (bits) Value\n\n");
      }
      else
      {
	fprint(out, "Sequences producing significant alignments:                         Score\n\n");
      }
	  
      for(long i=0; i<showhits; i++)
      {
	long headerlen = 67;
	if (parameters.symtype == SymbolType::blastn)
	{
	  headerlen = 65;
	}
	else if ((parameters.symtype == SymbolType::blastx) || (parameters.symtype == SymbolType::tblastn))
	{
	  headerlen = 64;
	}
	else if (parameters.symtype == SymbolType::tblastx)
	{
	  headerlen = 61;
	}

	HeaderLayout layout;
	layout.show_gis = show_gis;
	layout.maxlen = headerlen;
	layout.linelen = headerlen;
	db_showheader(t, 
		      make_view(hit_entry(i).header_address), layout);

	long const score = hit_entry(i).score;

	if (parameters.symtype == SymbolType::blastn)
	{
	  fprint(out, ' ');
	  fprint(out, (hit_entry(i).dstrand != 0) ? '-' : '+');
	}
	else if (parameters.symtype == SymbolType::blastx)
	{
	  fprintf(out, " %c%ld", (hit_entry(i).qstrand != 0) ? '-' : '+',
		 hit_entry(i).qframe+1);
	}
	else if (parameters.symtype == SymbolType::tblastn)
	{
	  fprintf(out, " %c%ld", (hit_entry(i).dstrand != 0) ? '-' : '+',
		 hit_entry(i).dframe+1);
	}
	else if (parameters.symtype == SymbolType::tblastx)
	{
	  fprintf(out, " %c%ld/%c%ld",
		  (hit_entry(i).qstrand != 0) ? '-' : '+',
		  hit_entry(i).qframe + 1,
		  (hit_entry(i).dstrand != 0) ? '-' : '+',
		  hit_entry(i).dframe + 1);
	}

	if (stats_available != 0)
	{
	  long const bits = static_cast<long>(floor(bit_score_of(score) + 0.5));
	  double const expect_value = expect_value_of(score);
		
	  fprint(out, ' ');
	  fprint_integer(out, bits, 5);
		
	  fprint(out, "   ");
		
	  hits_show_expect(expect_value);
	}
	else
	{
	  fprint(out, ' ');
	  fprint_integer(out, score, 5);
	}

	fprint(out, '\n');
      }

      for(long i=0; i<showalignments; i++)
      {
	fprint(out, "\n");
	HeaderLayout layout;
	layout.show_gis = show_gis;
	layout.indent = 10;
	layout.linelen = 79;
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
	      
	long const score = hit_entry(i).score;

	if (stats_available != 0)
	{
	  double const bits = bit_score_of(score);
	  double const expect_value = expect_value_of(score);
		
	  fprintf(out, " Score = %.1lf bits (%ld), Expect = ", bits, score);
	  hits_show_expect(expect_value);
	}
	else
	{
	  fprint(out, " Score = ");
	  fprint_integer(out, score);
	}

	fprint(out, '\n');

	long identities = 0;
	long positives = 0;
	long gaps = 0;
	long aligned = 0;
	long indels = 0;

	AlignedHit const hit = aligned_hit(parameters, i);
	count_align(hit, & identities, & positives, & indels, & aligned, & gaps);
	      
	fprintf(out, " Identities = %ld/%ld (%ld%%)",
	       identities, aligned, identities * 100 / aligned);
	if (parameters.symtype > SymbolType::blastn)
	{
	  fprintf(out, ", Positives = %ld/%ld (%ld%%)",
		  positives, aligned, positives * 100 / aligned);
	}
	if (indels != 0)
	{
	  fprint(out, ", Gaps = ");
	  fprint_integer(out, indels);
	  fprint(out, '/');
	  fprint_integer(out, aligned);
	  fprint(out, " (");
	  fprint_integer(out, indels * 100 / aligned);
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
	  fprintf(out, " Frame = %c%ld / %c%ld\n",
		  (hit_entry(i).qstrand != 0) ? '-' : '+',
		  hit_entry(i).qframe + 1,
		  (hit_entry(i).dstrand != 0) ? '-' : '+',
		  hit_entry(i).dframe + 1);
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
      fprintf(out, "%s\n\n%s\n\n", 
	      swipe_name_and_version, 
	      "Reference: T. Rognes (2011) Faster Smith-Waterman database searches\nwith inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.");
    }
  else if (view==OutputFormat::xml)
    {
      // one root element around the results of all queries (KI-26)
      fprint(out, "<?xml version=\"1.0\"?>\n");
      fprint(out, "<results>\n");
    }
  else if (view==OutputFormat::paralign_xml)
    {
      char url1[] = "http://www.w3.org/2001/XMLSchema-instance";
      char url2[] = "http://www.paralign.org/ParalignXML.xsd";

      fprint(out, "<?xml version=\"1.0\"?>\n");
      fprintf(out, "<ParalignXML xmlns:xsi=\"%s\" xsi:noNamespaceSchemaLocation=\"%s\">\n",
	      url1, url2);
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
  long const show_gis = parameters.show_gis;

  // compute number of hits and alignments to actually show

  long showalignments = 0;
  long showhits = 0;

  if (hits_count < opt_descriptions)
  {
    showhits = hits_count;
  }
  else
  {
    showhits = opt_descriptions;
  }

  if (hits_count < opt_alignments)
  {
    showalignments = hits_count;
  }
  else
  {
    showalignments = opt_alignments;
  }

  struct db_thread_s * t = db_thread_create();

  if(view == OutputFormat::plain)
  {
    hits_show_plain(parameters, show_gis, showalignments, showhits, t);
  }
  else if (view==OutputFormat::xml)
  {
    hits_show_xml(parameters, show_gis, showalignments, showhits, t);
  }
  else if ((view==OutputFormat::tabular)||(view==OutputFormat::tabular_with_comments))
  {
    hits_show_tsv(parameters, showalignments, static_cast<long>(view == OutputFormat::tabular_with_comments), t);
  }
  else if (view==OutputFormat::paralign_xml)
  {
    hits_show_xml_paralign(parameters, showalignments, showhits, t);
  }
  db_thread_destruct(t);
}

