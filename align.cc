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
#include <algorithm>  // std::fill_n, std::max
#include <cassert>
#include <initializer_list>  // std::max({...})
#include <cstddef>  // std::size_t
#include <string>  // std::string, std::to_string
#include <utility>  // std::move

// These functions are based on the following articles:
// - Huang, Hardison & Miller (1990) CABIOS 6:373-381
// - Myers & Miller (1988) CABIOS 4:11-17

// For consistency with non-symmetric score matrices:
// a (of length M) is the query sequence
// b (of length N) is the database sequence

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

// a cell of the alignment matrix: position a in the query sequence,
// position b in the database sequence
struct Cell
{
  long a;
  long b;
};

// Reverse pass of region(): from the end cell of the best local
// alignment, whose score is known, find the cell where it begins. HH
// and EE are work arrays of at least end.b + 1 elements.
auto region_begin(char const * a_seq,
		  char const * b_seq,
		  long const * scorematrix,
		  long const q,
		  long const r,
		  Cell const end,
		  long const score,
		  long * HH,
		  long * EE) -> Cell
{
  assert(end.b >= 0);
  std::fill_n(HH, end.b + 1, -1L);
  std::fill_n(EE, end.b + 1, -1L);

  long Cost = 0;

  for (long i = end.a; i >= 0; i--)
    {
      long h = -1;
      long f = -1;
      long p = 0;
      if (i == end.a)
      {
	p = 0;
      }
      else
      {
	p = -1;
      }
      for (long j = end.b; j >= 0; j--)
	{
	  f = std::max(f, h - q) - r;
	  EE[j] = std::max(EE[j], HH[j] - q) - r;

	  h = p + (scorematrix + (b_seq[j]<<5))[static_cast<int>(a_seq[i])];

	  h = std::max({h, f, EE[j]});


	  p = HH[j];

	  HH[j] = h;

	  if (h > Cost)
	    {
	      Cost = h;
	      if (Cost >= score)
	      {
		return {i, j};
	      }
	    }
	}
    }

  fatal("Internal error in align function.");
}

// the region of the best local alignment of a_sequence and
// b_sequence; with a non-zero hint.score, its end cell is that of the
// hint, and only the reverse pass is run
auto region(View<char> const a_sequence,
	    View<char> const b_sequence,
	    long const * scorematrix,
	    GapPenalties const gaps,
	    AlignmentRegion const & hint) -> AlignmentRegion
{
  auto const * const a_seq = a_sequence.data();
  auto const * const b_seq = b_sequence.data();
  auto const M = static_cast<long>(a_sequence.size());
  auto const N = static_cast<long>(b_sequence.size());
  long const q = gaps.open;
  long const r = gaps.extend;
  long a_end = hint.a_end;
  long b_end = hint.b_end;

  Buffer<long> hh_buffer(static_cast<std::size_t>(N));
  Buffer<long> ee_buffer(static_cast<std::size_t>(N));
  long * HH = hh_buffer.data();
  long * EE = ee_buffer.data();

  long score = 0;

  // Forward pass

  if (hint.score != 0)
  {
    score = hint.score;
  }
  else
  {

    std::fill_n(HH, N, 0L);
    std::fill_n(EE, N, - q);
    
    for (long i = 0; i < M; i++)
    {
      long h = 0;
      long p = 0;
      long f = - q;
      for (long j = 0; j < N; j++)
      {
	f = std::max(f, h - q) - r;
	EE[j] = std::max(EE[j], HH[j] - q) - r;
	
	h = p + (scorematrix + (b_seq[j]<<5))[static_cast<int>(a_seq[i])];
	
	h = std::max({h, 0L, f, EE[j]});
	
	p = HH[j];
	
	HH[j] = h;
	
	if (h > score)
	{
	  score = h;
	  a_end = i;
	  b_end = j;
	}
      }
    }
  }

  // Reverse pass

  auto const begin = region_begin(a_seq, b_seq, scorematrix, q, r,
                                  {a_end, b_end}, score, HH, EE);
  return {begin.a, begin.b, a_end, b_end, score};
}

struct aligner_info
{
  char op;
  long count;
  std::string alignment;
};

auto push(aligner_info & info) -> void
{
  if (info.count > 0)
  {
    // the operation and its length, as "%c%ld" (e.g. M12)
    info.alignment += info.op;
    info.alignment += std::to_string(info.count);
  }
}

auto newop(aligner_info & info, char op, long len) -> void
{
  if (info.op == op)
  {
    info.count += len;
  }
  else
  {
    push(info);
    info.op = op;
    info.count = len;
  }
}

auto delete_a(aligner_info & info, long len) -> void
{
  newop(info, 'D', len);
}

auto insert_b(aligner_info & info, long len) -> void
{
  newop(info, 'I', len);
}

auto match(aligner_info & info) -> void
{
  newop(info, 'M', 1);
}

// what the recursion of diff() does not change: the two sequences,
// the score matrix and the gap penalties
struct DiffInput
{
  char const * a_seq;
  char const * b_seq;
  long const * scorematrix;
  GapPenalties gaps;
};

// a block of the alignment matrix: M residues of a from a_pos, N
// residues of b from b_pos
struct DiffBlock
{
  long a_pos;
  long b_pos;
  long M;
  long N;
};

// the gap open penalties at the left (tb) and right (te) ends of a
// block: 0 when a gap is already open there, q otherwise
struct EndGaps
{
  long tb;
  long te;
};

auto diff(aligner_info & info,
	  DiffInput const & input,
	  DiffBlock const & block,
	  EndGaps const ends) -> void
{
  auto const * const a_seq = input.a_seq;
  auto const * const b_seq = input.b_seq;
  auto const * const scorematrix = input.scorematrix;
  long const q = input.gaps.open;
  long const r = input.gaps.extend;
  long const a_pos = block.a_pos;
  long const b_pos = block.b_pos;
  long const M = block.M;
  long const N = block.N;
  long const tb = ends.tb;
  long const te = ends.te;

  if (N == 0)
    {
      if (M > 0)
      {
	delete_a(info, M);
      }
    }
  else if (M == 0)
    {
      insert_b(info, N);
    }
  else if (M == 1)
    {
      // Conversion (1 char from A, N chars from B)

      // tb = gap open penalty on extreme left
      // te = gap open penalty on extreme right
      // tb = 0 or q depending on whether a gap is already open on left of B
      // te = 0 or q depending on whether a gap is already open on right of B

      long MaxScore = 0;
      long J = 0;

      if (tb <= te)
	{
	  // Delete 1 from A, Insert N from B
	  // A----
	  // -BBBB

	  MaxScore = - tb - ((1 + N) * r) - q;
	  J = -1;
	}
      else
	{
	  // Insert N from B, Delete 1 from A
	  // ----A
	  // BBBB-

	  MaxScore = - q - ((1 + N) * r) - te;
	  J = N;
	}

      for (long j = 0; j < N; j++)
	{
	  // Insert J from B, replace 1, insert rest of B
	  // -A--
	  // BBBB

	  long Score = (scorematrix + (b_seq[b_pos+j]<<5))[static_cast<int>(a_seq[a_pos])] - (r * (N-1));

	  if (j > 0)
	  {
	    Score -= q;
	  }
	  if (j < N - 1)
	  {
	    Score -= q;
	  }

	  if (Score > MaxScore)
	    {
	      MaxScore = Score;
	      J = j;
	    }
	}

      if (J == -1)
	{
	  delete_a(info, 1);
	  insert_b(info, N);
	}
      else if (J == N)
	{
	  insert_b(info, N);
	  delete_a(info, 1);
	}
      else
	{
	  if (J > 0)
	  {
	    insert_b(info, J);
	  }
	  match(info);
	  if (J < N - 1)
	  {
	    insert_b(info, N - 1 - J);
	  }
	}
    }
  else
    {

      long const I = M/2;
      long i = 0;
      long j = 0;
      long t = 0;

      // Compute HH & EE in forward phase with tb

      Buffer<long> hh_buffer(static_cast<std::size_t>(N) + 1);
      Buffer<long> ee_buffer(static_cast<std::size_t>(N) + 1);
      long * HH = hh_buffer.data();
      long * EE = ee_buffer.data();

      HH[0] = 0;
      t = -q;
      for (j = 1; j <= N; j++)
	{
	  t -= r;
	  HH[j] = t;
	  EE[j] = t - q;
	}
      t = -tb;
      for (i = 1; i <= I; i++)
	{
	  long p = HH[0];
	  t -= r;
	  long h = t;
	  HH[0] = t;
	  long f = t - q;

	  for (j = 1; j <= N; j++)
	    {
	      f = std::max(f, h - q) - r;
	      EE[j] = std::max(EE[j], HH[j] - q) - r;

	      h = p + (scorematrix + (b_seq[b_pos+j-1]<<5))[static_cast<int>(a_seq[a_pos+i-1])];

	      h = std::max({h, f, EE[j]});
	      p = HH[j];
	      HH[j] = h;
	    }
	}
      EE[0] = HH[0];


      // Compute XX & YY in reverse phase with te

      Buffer<long> xx_buffer(static_cast<std::size_t>(N) + 1);
      Buffer<long> yy_buffer(static_cast<std::size_t>(N) + 1);
      long * XX = xx_buffer.data();
      long * YY = yy_buffer.data();

      XX[0] = 0;
      t = -q;
      for (j = 1; j <= N; j++)
	{
	  t -= r;
	  XX[j] = t;
	  YY[j] = t - q;
	}

      t = -te;
      for (i = 1; i <= M-I; i++)
	{
	  long p = XX[0];
	  t -= r;
	  long h = t;
	  XX[0] = t;
	  long f = t - q;

	  for (j = 1; j <= N; j++)
	    {
	      f = std::max(f, h - q) - r;
	      YY[j] = std::max(YY[j], XX[j] - q) - r;

	      h = p + (scorematrix + (b_seq[b_pos+N-j]<<5))[static_cast<int>(a_seq[a_pos+M-i])];

	      h = std::max({h, f, YY[j]});
	      p = XX[j];
	      XX[j] = h;
	    }
	}
      YY[0] = XX[0];




      long MaxScore = LONG_MIN;
      long P = -1;
      long J = -1;

      for (j=0; j <= N; j++)
	{
	  long const Score = HH[j] + XX[N-j];
	  if (Score > MaxScore)
	    {
	      MaxScore = Score;
	      P = 0;
	      J = j;
	    }
	}

      // released before the recursive calls (peak memory: one level)
      Buffer<long>().swap(hh_buffer);
      Buffer<long>().swap(xx_buffer);

      for (j=0; j <= N; j++)
	{
	  long const Score = EE[j] + YY[N-j] + q;
	  if (Score >= MaxScore)
	    {
	      MaxScore = Score;
	      P = 1;
	      J = j;
	    }
	}

      Buffer<long>().swap(ee_buffer);
      Buffer<long>().swap(yy_buffer);

      if (P == 0)
	{
	  diff(info, input, {a_pos, b_pos, I, J}, {tb, q});
	  diff(info, input, {a_pos+I, b_pos+J, M-I, N-J}, {q, te});
	}
      else if (P == 1)
	{
	  diff(info, input, {a_pos, b_pos, I-1, J}, {tb, 0});
	  delete_a(info, 2);
	  diff(info, input, {a_pos+I+1, b_pos+J, M-I-1, N-J}, {0, te});
	}
    }
}

}  // anonymous namespace

auto align(View<char> const query_sequence,
	   View<char> const database_sequence,
	   long const * scorematrix,
	   GapPenalties const gaps,
	   AlignmentRegion const & hint,
	   std::string & alignment) -> AlignmentRegion
{
  aligner_info ai {0, 0, std::string()};

  auto const result = region(query_sequence, database_sequence, scorematrix, gaps, hint);

  diff(ai,
       {query_sequence.data(), database_sequence.data(), scorematrix, gaps},
       {result.a_begin, result.b_begin,
        result.a_end - result.a_begin + 1, result.b_end - result.b_begin + 1},
       {gaps.open, gaps.open});

  push(ai);

  alignment = std::move(ai.alignment);
  return result;
}
