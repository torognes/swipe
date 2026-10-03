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
#include "fatal_allocator.h"  // Buffer
#include "view.h"  // View
#include <algorithm>  // std::fill_n, std::max
#include <array>
#include <cassert>
#include <initializer_list>  // std::max({...})
#include <cstddef>  // std::size_t
#include <limits>  // std::numeric_limits
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

// the two sequences, the score matrix and the gap penalties of an
// alignment: what region_begin() and the recursion of diff() do not
// change
struct AlignmentInput
{
  View<char> a_seq;
  View<char> b_seq;
  ScoreMatrix63 const & scorematrix;
  GapPenalties gaps;

  // the score of query position a against database position b
  auto score(long const a, long const b) const -> long
  {
    assert((a >= 0) and (b >= 0));
    auto const query_symbol = static_cast<unsigned char>(a_seq[static_cast<std::size_t>(a)]);
    auto const database_symbol = static_cast<unsigned char>(b_seq[static_cast<std::size_t>(b)]);
    return scorematrix[score_matrix_cell(database_symbol, query_symbol)];
  }
};

// a work array of the alignment (HH, EE, XX, YY): one score per
// position of the database sequence, indexed as the positions are
// (long)
class Row
{
public:
  explicit Row(long const length) : cells_(static_cast<std::size_t>(length)) {}

  auto operator[](long const j) -> long &
  {
    assert((j >= 0) and (static_cast<std::size_t>(j) < cells_.size()));
    return cells_[static_cast<std::size_t>(j)];
  }

  // the first 'length' cells set to 'value'
  auto fill(long const length, long const value) -> void
  {
    assert((length >= 0) and (static_cast<std::size_t>(length) <= cells_.size()));
    std::fill_n(cells_.begin(), length, value);
  }

  // the memory given back
  auto release() -> void
  {
    Buffer<long>().swap(cells_);
  }

private:
  Buffer<long> cells_;
};

// Reverse pass of region(): from the end cell of the best local
// alignment, whose score is known, find the cell where it begins. HH
// and EE are work arrays of at least end.b + 1 elements.
auto region_begin(AlignmentInput const & input,
		  Cell const end,
		  long const score,
		  Row & HH,
		  Row & EE) -> Cell
{
  long const q = input.gaps.open;
  long const r = input.gaps.extend;
  assert(end.b >= 0);
  HH.fill(end.b + 1, -1L);
  EE.fill(end.b + 1, -1L);

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

	  h = p + input.score(i, j);

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
auto region(AlignmentInput const & input,
	    AlignmentRegion const & hint) -> AlignmentRegion
{
  auto const M = static_cast<long>(input.a_seq.size());
  auto const N = static_cast<long>(input.b_seq.size());
  long const q = input.gaps.open;
  long const r = input.gaps.extend;
  long a_end = hint.a_end;
  long b_end = hint.b_end;

  Row HH(N);
  Row EE(N);

  long score = 0;

  // Forward pass

  if (hint.score != 0)
  {
    score = hint.score;
  }
  else
  {

    HH.fill(N, 0L);
    EE.fill(N, - q);
    
    for (long i = 0; i < M; i++)
    {
      long h = 0;
      long p = 0;
      long f = - q;
      for (long j = 0; j < N; j++)
      {
	f = std::max(f, h - q) - r;
	EE[j] = std::max(EE[j], HH[j] - q) - r;
	
	h = p + input.score(i, j);
	
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

  auto const begin = region_begin(input, {a_end, b_end}, score, HH, EE);
  return {begin.a, begin.b, a_end, b_end, score};
}

// the operations of an alignment, as they are written in it
enum struct Operation : char
{
  none = 0,
  deletion = 'D',
  insertion = 'I',
  match = 'M'
};

struct aligner_info
{
  Operation op;
  long count;
  std::string alignment;
};

auto push(aligner_info & info) -> void
{
  if (info.count > 0)
  {
    // the operation and its length, as "%c%ld" (e.g. M12)
    info.alignment += static_cast<char>(info.op);
    info.alignment += std::to_string(info.count);
  }
}

auto newop(aligner_info & info, Operation const op, long const len) -> void
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
  newop(info, Operation::deletion, len);
}

auto insert_b(aligner_info & info, long len) -> void
{
  newop(info, Operation::insertion, len);
}

auto match(aligner_info & info) -> void
{
  newop(info, Operation::match, 1);
}

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
	  AlignmentInput const & input,
	  DiffBlock const & block,
	  EndGaps const ends) -> void
{
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

	  long Score = input.score(a_pos, b_pos + j) - (r * (N-1));

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

      Row HH(N + 1);
      Row EE(N + 1);

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

	      h = p + input.score(a_pos + i - 1, b_pos + j - 1);

	      h = std::max({h, f, EE[j]});
	      p = HH[j];
	      HH[j] = h;
	    }
	}
      EE[0] = HH[0];


      // Compute XX & YY in reverse phase with te

      Row XX(N + 1);
      Row YY(N + 1);

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

	      h = p + input.score(a_pos + M - i, b_pos + N - j);

	      h = std::max({h, f, YY[j]});
	      p = XX[j];
	      XX[j] = h;
	    }
	}
      YY[0] = XX[0];




      long MaxScore = std::numeric_limits<long>::min();
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
      HH.release();
      XX.release();

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

      EE.release();
      YY.release();

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
	   ScoreMatrix63 const & scorematrix,
	   GapPenalties const gaps,
	   AlignmentRegion const & hint,
	   std::string & alignment) -> AlignmentRegion
{
  aligner_info ai {Operation::none, 0, std::string()};
  AlignmentInput const input {query_sequence, database_sequence, scorematrix, gaps};

  auto const result = region(input, hint);

  diff(ai,
       input,
       {result.a_begin, result.b_begin,
        result.a_end - result.a_begin + 1, result.b_end - result.b_begin + 1},
       {gaps.open, gaps.open});

  push(ai);

  alignment = std::move(ai.alignment);
  return result;
}
