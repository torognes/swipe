/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2013 Torbjorn Rognes, University of Oslo, 
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

// These functions are based on the following articles:
// - Huang, Hardison & Miller (1990) CABIOS 6:373-381
// - Myers & Miller (1988) CABIOS 4:11-17

#define MAX(a,b) (a > b ? a : b)

void region(char * a_seq,
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
	    long * s)
{
  
  long * HH = (long *) xmalloc(N * sizeof(long));
  long * EE = (long *) xmalloc(N * sizeof(long));

  long i, j;

  long score = 0;

  // Forward pass

#if 1
  if (*s)
  {
    score = *s;
  }
  else
#endif
  {

    for (j = 0; j < N; j++)
    {
      HH[j] = 0;
      EE[j] = - q;
    }
    
    for (i = 0; i < M; i++)
    {
      long h = 0;
      long p = 0;
      long f = - q;
      for (j = 0; j < N; j++)
      {
	f = MAX(f, h - q) - r;
	EE[j] = MAX(EE[j], HH[j] - q) - r;
	
	h = p + (scorematrix + (a_seq[i]<<5))[(int)(b_seq[j])];
	
	if (h < 0)
	  h = 0;
	if (f > h)
	  h = f;
	if (EE[j] > h)
	  h = EE[j];
	
	p = HH[j];
	
	HH[j] = h;
	
	if (h > score)
	{
	  score = h;
	  *a_end = i;
	  *b_end = j;
	}
      }
    }
  }

  // Reverse pass

  for (j = *b_end; j >= 0; j--)
    {
      HH[j] = -1;
      EE[j] = -1;
    }

  long Cost = 0;

  for (i = *a_end; i >= 0; i--)
    {
      long h = -1;
      long f = -1;
      long p;
      if (i == *a_end)
	p = 0;
      else
	p = -1;
      for (j = *b_end; j >= 0; j--)
	{
	  f = MAX(f, h - q) - r;
	  EE[j] = MAX(EE[j], HH[j] - q) - r;

	  h = p + (scorematrix + (a_seq[i]<<5))[(int)(b_seq[j])];

	  if (f > h)
	    h = f;
	  if (EE[j] > h)
	    h = EE[j];


	  p = HH[j];

	  HH[j] = h;

	  if (h > Cost)
	    {
	      Cost = h;
	      *a_begin = i;
	      *b_begin = j;
	      if (Cost >= score)
		goto Found;
	    }
	}
    }

  fatal("Internal error in align function.");

 Found:

  free(EE);
  free(HH);

  *s = score;
}

struct aligner_info
{
  char op;
  long count;
  char * alignment;
  long length;
  long size;
};

void init(struct aligner_info * aip)
{
  aip->op = 0;
  aip->count = 0;
  aip->size = 64;
  aip->alignment = (char*) xmalloc(aip->size);
  aip->length = 0;
}

void push(struct aligner_info * aip)
{
  if (aip->count > 0)
  {
    while (1)
    {
      long rest = aip->size - aip->length;
      int n = snprintf(aip->alignment + aip->length,
		       rest,
		       "%c%ld", aip->op, aip->count);
      if ((n < 0) || (n >= rest))
      {
	aip->size += 64;
	aip->alignment = (char*) xrealloc(aip->alignment, aip->size);
	//	fprintf(stderr, "Reallocating memory for alignment: %ld\n", aip->size);
      }
      else
      {
	aip->length += n;
	break;
      }
    }
  }
}

void newop(struct aligner_info * aip, char op, long len)
{
  if (aip->op == op)
    aip->count += len;
  else
    {
      push(aip);
      aip->op = op;
      aip->count = len;
    }
}

void delete_a(struct aligner_info * aip, long len)
{
  newop(aip, 'D', len);
}

void insert_b(struct aligner_info * aip, long len)
{
  newop(aip, 'I', len);
}

void match(struct aligner_info * aip)
{
  newop(aip, 'M', 1);
}

void diff(struct aligner_info * aip,
	  char * a_seq,
	  char * b_seq,
	  long M,
	  long N,
	  long a_pos,
	  long b_pos,
	  long * scorematrix,
	  long q,
	  long r,
	  long tb,
	  long te)
{
  long MaxScore = 0;

  if (N == 0)
    {
      if (M > 0)
	delete_a(aip, M);
    }
  else if (M == 0)
    {
      insert_b(aip, N);
    }
  else if (M == 1)
    {
      // Conversion (1 char from A, N chars from B)

      // tb = gap open penalty on extreme left
      // te = gap open penalty on extreme right
      // tb = 0 or q depending on whether a gap is already open on left of B
      // te = 0 or q depending on whether a gap is already open on right of B

      long J;

      if (tb <= te)
	{
	  // Delete 1 from A, Insert N from B
	  // A----
	  // -BBBB

	  MaxScore = - tb - (1 + N) * r - q;
	  J = -1;
	}
      else
	{
	  // Insert N from B, Delete 1 from A
	  // ----A
	  // BBBB-

	  MaxScore = - q - (1 + N) * r - te;
	  J = N;
	}

      for (long j = 0; j < N; j++)
	{
	  // Insert J from B, replace 1, insert rest of B
	  // -A--
	  // BBBB

	  long Score = (scorematrix + (a_seq[a_pos]<<5))[(int)(b_seq[b_pos+j])] - r * (N-1);

	  if (j > 0)
	    Score -= q;
	  if (j < N-1)
	    Score -= q;

	  if (Score > MaxScore)
	    {
	      MaxScore = Score;
	      J = j;
	    }
	}

      if (J == -1)
	{
	  delete_a(aip, 1);
	  insert_b(aip, N);
	}
      else if (J == N)
	{
	  insert_b(aip, N);
	  delete_a(aip, 1);
	}
      else
	{
	  if (J > 0)
	    insert_b(aip, J);
	  match(aip);
	  if (J < N-1)
	    insert_b(aip, N-1-J);
	}
    }
  else
    {

      long I = M/2;
      long i, j;
      long t;

      // Compute HH & EE in forward phase with tb

      long * HH = (long *) xmalloc((N+1) * sizeof(long));
      long * EE = (long *) xmalloc((N+1) * sizeof(long));

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
	      f = MAX(f, h - q) - r;
	      EE[j] = MAX(EE[j], HH[j] - q) - r;

	      h = p + (scorematrix + (a_seq[a_pos+i-1]<<5))[(int)(b_seq[b_pos+j-1])];

	      if (f > h)
		h = f;
	      if (EE[j] > h)
		h = EE[j];
	      p = HH[j];
	      HH[j] = h;
	    }
	}
      EE[0] = HH[0];


      // Compute XX & YY in reverse phase with te

      long * XX = (long *) xmalloc((N+1) * sizeof(long));
      long * YY = (long *) xmalloc((N+1) * sizeof(long));

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
	      f = MAX(f, h - q) - r;
	      YY[j] = MAX(YY[j], XX[j] - q) - r;

	      h = p + (scorematrix + (a_seq[a_pos+M-i]<<5))[(int)(b_seq[b_pos+N-j])];

	      if (f > h)
		h = f;
	      if (YY[j] > h)
		h = YY[j];
	      p = XX[j];
	      XX[j] = h;
	    }
	}
      YY[0] = XX[0];




      MaxScore = LONG_MIN;
      long P = -1;
      long J = -1;

      for (j=0; j <= N; j++)
	{
	  long Score = HH[j] + XX[N-j];
	  if (Score > MaxScore)
	    {
	      MaxScore = Score;
	      P = 0;
	      J = j;
	    }
	}

      free(HH);
      free(XX);

      for (j=0; j <= N; j++)
	{
	  long Score = EE[j] + YY[N-j] + q;
	  if (Score >= MaxScore)
	    {
	      MaxScore = Score;
	      P = 1;
	      J = j;
	    }
	}

      free(EE);
      free(YY);

      if (P == 0)
	{
	  diff(aip, a_seq, b_seq, I, J, a_pos, b_pos,
	       scorematrix, q, r, tb, q);
	  diff(aip, a_seq, b_seq, M-I, N-J, a_pos+I, b_pos+J, 
	       scorematrix, q, r, q, te);
	}
      else if (P == 1)
	{
	  diff(aip, a_seq, b_seq, I-1, J, a_pos, b_pos,
	       scorematrix, q, r, tb, 0);
	  delete_a(aip, 2);
	  diff(aip, a_seq, b_seq, M-I-1, N-J, a_pos+I+1, b_pos+J,
	       scorematrix, q, r, 0, te);
	}
    }
}

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
	   long * s)
{
  struct aligner_info ai;

  long score = *s;

  init(&ai);

  region(a_seq,
	 b_seq,
	 M,
	 N,
	 scorematrix,
	 q,
	 r,
	 a_begin,
	 b_begin,
	 a_end,
	 b_end,
	 & score);
  
  diff(& ai,
       a_seq,
       b_seq,
       *a_end - *a_begin + 1,
       *b_end - *b_begin + 1,
       *a_begin, 
       *b_begin, 
       scorematrix,
       q,
       r,
       q,
       q);

  push(& ai);

  *alignment = ai.alignment;
  *s = score;
}
