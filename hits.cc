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

long keephits;
long scorethreshold;
long upperscorethreshold;
int hits_count;
long init_threshold;
long obvious;

long opt_descriptions;
long opt_alignments;
long opt_show_best;

/* parameters for bit scores and expect values */

long stats_available = 0;

double alpha;
double beta;
double lambda;
double K;
double H;
double Kmn = 0;
double logK;
double lambda_d_log2;
double logK_d_log2;

struct hits_entry
{
  char * alignment;
  char * dseq;
  char * header_address;
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
  long header_length;
} * hits_list;


pthread_mutex_t hitsmutex = PTHREAD_MUTEX_INITIALIZER;

int hits_compare(const void * a, const void * b)
{
  struct hits_entry * ap = hits_list + *((long *) a);
  struct hits_entry * bp = hits_list + *((long *) b);
  
  if ( (*((long*)a) >= opt_alignments) < (*((long*)b) >= opt_alignments) )
  {
    return -1;
  }
  else if ( (*((long*)a) >= opt_alignments) > (*((long*)b) >= opt_alignments) )
  {
    return +1;
  }
  else if (ap->qstrand < bp->qstrand)
  {
    return -1;
  }
  else if (ap->qstrand > bp->qstrand)
  {
    return +1;
  }
  else if (ap->qframe < bp->qframe)
  {
    return -1;
  }
  else if (ap->qframe > bp->qframe)
  {
    return +1;
  }
  else if (ap->seqno < bp->seqno)
  {
    return -1;
  }
  else if (ap->seqno > bp->seqno)
  {
    return +1;
  }
  else if (ap->dstrand < bp->dstrand)
  {
    return -1;
  }
  else if (ap->dstrand > bp->dstrand)
  {
    return +1;
  }
  else if (ap->dframe < bp->dframe)
  {
    return -1;
  }
  else if (ap->dframe > bp->dframe)
  {
    return +1;
  }

  return 0;
}

long * hits_sort()
{
  long * hits_sorted = (long *) xmalloc(hits_count * sizeof(long));
  
  for(long i=0; i<hits_count; i++)
    hits_sorted[i] = i;
  
  qsort(hits_sorted, hits_count, sizeof(long), hits_compare);

#if 0
  fprintf(out, "Hits sorted:\n");
  fprintf(out, "Hitno Order Vol    Seqno S F S F Score\n");
  for(long i=0; i<hits_count; i++)
  {
    struct hits_entry * h = hits_list + hits_sorted[i];
    fprintf(out, "%5ld %5ld %3ld %8ld %ld %ld %ld %ld %5ld\n", 
	    h - hits_list + 1,
	    i+1,
	    db_getvolume(h->seqno),
	    h->seqno, h->qstrand, h->qframe, h->dstrand, h->dframe, h->score);
  }
  fprintf(out, "\n");
#endif

  return hits_sorted;
}


void hits_enter(long seqno, long score, long qstrand, long qframe,
		long dstrand, long dframe, long align_hint, long bestq)
{
  // show_progress();

  //  fprintf(out, "Entering score %u from sequence no %u.\n", score, seqno);
  
  // find correct place

  pthread_mutex_lock(&hitsmutex);

  if (score > upperscorethreshold)
    obvious++;
  
  if (score >= init_threshold)
    totalhits++;

  if ((score < scorethreshold) || (score > upperscorethreshold))
  {
    pthread_mutex_unlock(&hitsmutex);
    return;
  }

  long place = hits_count;

  while ((place > 0) && ((score > hits_list[place-1].score) || 
			 ((score == hits_list[place-1].score) && 
			  (seqno > hits_list[place-1].seqno))))
    place--;

  long Nbest = 1;
  if (place > 0)
  {
    if (score < hits_list[place-1].score)
      Nbest++;
    
    for (long j=(place-1); j>0; j--)
    {
      if (hits_list[j].score < hits_list[j-1].score)
        Nbest++;
    }
  }

  if (Nbest > opt_show_best)
  {
    pthread_mutex_unlock(&hitsmutex);
    return;
  }

  // move entries down
  
  long move = (hits_count < keephits ? hits_count : keephits - 1) - place;

  //  fprintf(out, "Inserting at place %d, moving %d.\n", place, move);

  for (long j=move; j>0; j--)
    hits_list[place+j] = hits_list[place+j-1];
  
  // fill new entry

  if (place < keephits)
  {
    hits_list[place].seqno = seqno;
    hits_list[place].qstrand = qstrand;
    hits_list[place].qframe = qframe;
    hits_list[place].dstrand = dstrand;
    hits_list[place].dframe = dframe;
    hits_list[place].score = score;
    hits_list[place].align_hint = align_hint;
    hits_list[place].bestq = bestq;
    if (hits_count < keephits)
      hits_count++;
  }
  
  if (hits_count == keephits)
    scorethreshold = hits_list[keephits-1].score;
  
  pthread_mutex_unlock(&hitsmutex);
}

long hits_getcount()
{
  return hits_count;
}

void hits_gethit(long i, long * seqno, long * score, 
		 long * qstrand, long * qframe,
		 long * dstrand, long * dframe)
{
  struct hits_entry * h = hits_list + i;
  *seqno = h->seqno;
  *score = h->score;
  *qstrand = h->qstrand;
  *qframe = h->qframe;
  *dstrand = h->dstrand;
  *dframe = h->dframe;
}

void hits_enter_seq(long i, char* seq, long seq_len)
{
  hits_list[i].dseq = (char*) xmalloc(seq_len);
  memcpy(hits_list[i].dseq, seq, seq_len);
  hits_list[i].dlen = seq_len;
}

void hits_enter_align_hint(long i, long q_end, long d_end)
{
  hits_list[i].bestq = q_end;
  hits_list[i].align_hint = d_end;
}

void hits_enter_align_coord(long i,
			    long align_q_start,
			    long align_q_end,
			    long align_d_start,
			    long align_d_end,
			    long dlennt)
{
  hits_list[i].align_q_start = align_q_start;
  hits_list[i].align_q_end = align_q_end;
  hits_list[i].align_d_start = align_d_start;
  hits_list[i].align_d_end = align_d_end;
  hits_list[i].dlennt = dlennt;
}

void hits_enter_header(long i, char * header, long header_len)
{
  hits_list[i].header_address = (char*) xmalloc(header_len);
  memcpy(hits_list[i].header_address, header, header_len);
  hits_list[i].header_length = header_len;
}

void hits_enter_align_string(long i, char * align, long align_len)
{
  hits_list[i].alignment = (char*) xmalloc(align_len);
  memcpy(hits_list[i].alignment, align, align_len);
  //  hits_list[i].alignment[align_len] = 0;
}

void hits_init(long descriptions, long alignments, long show_best, long minscore, long maxscore, double minexpect, double expect, int show_nostats)
{
  opt_descriptions = descriptions;
  opt_alignments = alignments;
  opt_show_best = show_best;
  keephits = descriptions > alignments ? descriptions : alignments;
  
  long maxhits = db_getseqcount_masked();
  if (symtype == 0)
    {
      if (querystrands == 3)
	maxhits *= 2;
    }
  else if (symtype == 2)
    {
      if (querystrands == 3)
	maxhits *= 6;
      else
	maxhits *= 3;
    }
  else if (symtype == 3)
    {
      maxhits *= 6;
    }
  else if (symtype == 4)
    {
      if (querystrands == 3)
	maxhits *= 36;
      else
	maxhits *= 18;
    }

  if (keephits > maxhits)
    keephits = maxhits;

  obvious = 0;
  hits_count = 0;
  hits_list = (struct hits_entry *) xmalloc(keephits * sizeof(struct hits_entry));

  for(int i=0; i<keephits; i++)
  {
    struct hits_entry * h = hits_list + i;
    h->header_address = 0;
    h->header_length = 0;
    h->dseq = 0;
    h->alignment = 0;
  }

  int seqcount;
  long symcount;

  if (db_ismasked())
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

  long m = 0;
  long n = 0;
  int lenadj = 0;
  if (symtype == 0)
  {
    if (stats_getparams_nt(matchscore,
			   mismatchscore,
			   gapopen,
			   gapextend,
			   & lambda,
			   & K,
			   & H,
			   & alpha,
			   & beta))
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
      
      long qlen = query.nt[0].len;

      long dlen;
      if (effdbsize > 0)
	dlen = effdbsize;
      else
	dlen = symcount;

      BlastComputeLengthAdjustment(K,
				   logK,
				   alpha / lambda,
				   beta,
				   qlen,
				   dlen,
				   seqcount,
				   & lenadj);
    
      //      fprintf(out, "lenadj: %d\n", lenadj);

      m = qlen - lenadj;

      if (effdbsize > 0)
	n = effdbsize;
      else
	n = dlen - seqcount * lenadj;

      Kmn = K * (double)m * (double)n;
    }
  }
  else if (symtype < 5)
  {
    if (symtype == 4)
    {
      stats_available = stats_getparams(matrixname,
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
      stats_available = stats_getparams(matrixname,
					gapopen,
					gapextend,
					& lambda,
					& K,
					& H,
					& alpha,
					& beta);
    }


    if (stats_available)
    {
      
#ifdef DEBUG
      fprintf(out, "Params: lambda=%6.3g K=%6.3g H=%6.3g alpha=%6.3g beta=%6.3g\n",
	      lambda, K, H, alpha, beta);
#endif

      logK = log(K);
      lambda_d_log2 = lambda / log(2.0);
      logK_d_log2 = logK / log(2.0);
      
      lenadj = 0;
      
      long qlen = query.aa[0].len;
      if ((symtype == 2) || (symtype == 4))
	qlen = query.nt[0].len / 3;

      long dlen;
      if (effdbsize > 0)
      {
	dlen = effdbsize;
      }
      else
      {
	if ((symtype == 3) || (symtype == 4))
	  dlen = symcount / 3;
	else
	  dlen = symcount;
      }

      BlastComputeLengthAdjustment(K,
				   logK,
				   alpha / lambda,
				   beta,
				   qlen,
				   dlen,
				   seqcount,
				   & lenadj);

      m = qlen - lenadj;

      if (effdbsize > 0)
	n = effdbsize;
      else
	n = dlen - seqcount * lenadj;

      Kmn = K * (double)m * (double)n;
    }
  }

#ifdef DEBUG
  fprintf(out, "lenadj=%d m=%ld n=%ld mn=%.1f\n", lenadj, m, n, (double)m * (double)n);
#endif

  scorethreshold = minscore;
  upperscorethreshold = maxscore;
  
  if (stats_available)
  {
    long minscore_expect = (long)(ceil(- log(expect / Kmn) / lambda));
    if (minscore_expect > minscore)
      scorethreshold = minscore_expect;
    
    if (minexpect > 0.0)
    {
      long maxscore_expect = (long)(floor(- log(minexpect / Kmn) / lambda));
      if (maxscore_expect < maxscore)
	upperscorethreshold = maxscore_expect;
    }
  }
  else
  {
    if (show_nostats)
      fprintf(out, "Statistical parameters are not available for the scoring system specified.\nBit scores and E-values will not be computed.\n\n");
  }

  init_threshold = scorethreshold;

  //  fprintf(out, "scorethreshold: %ld\n", scorethreshold);
}

void hits_empty()
{
  for (long i=0; i<hits_count; i++)
  {
    struct hits_entry * h = hits_list + i;

    if (h->header_address)
    {
      free(h->header_address);
      h->header_address = 0;
    }
    
    if (h->dseq)
    {
      free(h->dseq);
      h->dseq = 0;
    }

    if (h->alignment)
    {
      free(h->alignment);
      h->alignment = 0;
    }
  }
}

void hits_exit()
{
  hits_empty();
  free(hits_list);
  hits_list = 0;
}

void hits_align(struct db_thread_s * t, long i)
{
  char * address;
  long length, ntlen;

  struct hits_entry * h = hits_list + i;

  db_mapheaders(t, h->seqno, h->seqno);

  db_getheader(t, h->seqno, & address, & length);
  h->header_length = length;
  h->header_address = (char*) xmalloc(length);
  memcpy(h->header_address, address, length);

  if (i < opt_alignments)
  {
    db_mapsequences(t, h->seqno, h->seqno);
    
    db_getsequence(t, h->seqno, h->dstrand, h->dframe,
		   & address, & length, & ntlen, 0);
    h->dlen = length - 1;
    h->dlennt = ntlen;

    h->dseq = (char*)xmalloc(h->dlen);
    memcpy(h->dseq, address, h->dlen);
    
    char * qseq;
    long qlen;
    if (symtype == 0)
    {
      qseq = query.nt[0].seq;
      qlen = query.nt[0].len;
    }
    else
    {
      qseq = query.aa[3*h->qstrand + h->qframe].seq;
      qlen = query.aa[3*h->qstrand + h->qframe].len;
    }

    // give hint of alignment end

    if ((h->bestq > 0) && (h->align_hint))
    {
      h->score_align = h->score;
      h->align_q_end = h->bestq;
      h->align_d_end = h->align_hint;
#if 0
      printf("Align with hints: score=%ld, q_end=%ld, d_end=%ld\n", 
	     h->score_align, h->align_q_end, h->align_d_end);
#endif
    }
    else
    {
      h->score_align = 0;
      h->align_q_end = 0;
      h->align_d_end = 0;
    }

    align(qseq,
	  h->dseq,
	  qlen,
	  h->dlen,
	  score_matrix_63,
	  gapopen,
	  gapextend,
	  & h->align_q_start,
	  & h->align_d_start,
	  & h->align_q_end,
	  & h->align_d_end,
	  & h->alignment,
	  & h->score_align);
  }
}


#define ALIGNLEN 60
long line_pos;
long q_start;
long d_start;
long q_pos;
long d_pos;
long q_len;
long q_len_nt;
long d_len;
long d_len_nt;
long d_strand;
long d_frame;
long q_strand;
long q_frame;
long q_first;
long q_last;
long d_first;
long d_last;
char * q_seq;
char * d_seq;
char q_line[ALIGNLEN+1];
char a_line[ALIGNLEN+1];
char d_line[ALIGNLEN+1];
const char * sym;
int poswidth;

void putalignop(char c, long len)
{

  long count = len;
  while(count)
  {
    if (line_pos == 0)
    {
      q_start = q_pos;
      d_start = d_pos;
    }

    char qs;
    char ds;

    switch(c)
    {
    case 'M':
      qs = q_seq[q_pos++];
      ds = d_seq[d_pos++];
      q_line[line_pos] = sym[(int)(qs)];
      if (symtype == 0)
      {
	a_line[line_pos] = (qs == ds) ? '|' : ' ';
      }
      else
      {
	a_line[line_pos] = (qs == ds) ? sym[(int)(qs)] : 
	  (score_matrix_63[32*qs+ds] > 0 ? '+' : ' ');
      }
      d_line[line_pos] = sym[(int)(ds)];
      line_pos++;
      break;

    case 'D':
      qs = q_seq[q_pos++];
      q_line[line_pos] = sym[(int)(qs)];
      a_line[line_pos] = ' ';
      d_line[line_pos] = '-';
      line_pos++;
      break;

    case 'I':
      ds = d_seq[d_pos++];
      q_line[line_pos] = '-';
      a_line[line_pos] = ' ';
      d_line[line_pos] = sym[(int)(ds)];
      line_pos++;
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

      if (symtype == 0)
      {
	if (d_strand)
	{
	  d1 = d_len - d1 + 1;
	  d2 = d_len - d2 + 1;
	}
      }

      if ((symtype == 2) || (symtype == 4))
      {
	if (q_strand)
	{
	  q1 = q_len_nt - 3*q_start - q_frame;
	  q2 = q_len_nt - 3*q_pos - q_frame + 1;
	}
	else
	{
	  q1 = 3*q_start + q_frame + 1;
	  q2 = 3*q_pos + q_frame;
	}
      }
      
      if ((symtype == 3) || (symtype == 4))
      {
	if (d_strand)
	{
	  d1 = d_len_nt - 3*d_start - d_frame;
	  d2 = d_len_nt - 3*d_pos - d_frame + 1;
	}
	else
	{
	  d1 = 3*d_start + d_frame + 1;
	  d2 = 3*d_pos + d_frame;
	}
      }


      fprintf(out, "\n");
      fprintf(out, "Query: %*ld %s %ld\n", poswidth, q1, q_line, q2);
      fprintf(out, "       %*s %s\n", poswidth, "", a_line);
      fprintf(out, "Sbjct: %*ld %s %ld\n", poswidth, d1, d_line, d2);

      line_pos = 0;
    }

    count--;
  }
}

void show_align(long i)
{
  long q_align_start = hits_list[i].align_q_start;
  long d_align_start = hits_list[i].align_d_start;
  q_strand = hits_list[i].qstrand;
  q_frame = hits_list[i].qframe;
  d_strand = hits_list[i].dstrand;
  d_frame = hits_list[i].dframe;
  char * alignment = hits_list[i].alignment;
  
  if (symtype == 0)
  {
    sym = sym_ncbi_nt16;
    q_seq = query.nt[0].seq;
    q_len = query.nt[0].len;
  }
  else if (symtype == 5)
  {
    sym = sym_sound;
    q_seq = query.aa[0].seq;
    q_len = query.aa[0].len;
  }
  else
  {
    sym = sym_ncbi_aa;
    q_seq = query.aa[3*q_strand+q_frame].seq;
    q_len = query.aa[3*q_strand+q_frame].len;
    q_len_nt = query.nt[0].len;
    d_len_nt = hits_list[i].dlennt;
  }

  d_seq = hits_list[i].dseq;
  d_len = hits_list[i].dlen;
  
  q_pos = q_align_start;
  d_pos = d_align_start;
  
  char * p = alignment;
  char * e = alignment + strlen(alignment);
  
  while(p < e)
  {
    char op = *p++;
    long len;
    int n;
    sscanf(p, "%ld%n", & len, & n);
    p += n;
    putalignop(op, len);
  }
  
  putalignop(0, 1);
}

void whole_align(long i,
		 long * identities,
		 long * positives,
		 long * indels,
		 long * aligned,
		 long * gaps,
		 char ** qline,
		 char ** aline,
		 char ** dline)
{

  long al = 0;
  char * alignment = hits_list[i].alignment;
  char * p = alignment;
  while(*p)
  {
    long len;
    int n;
    sscanf(p, "%*c%ld%n", & len, & n);
    p += n;
    al += len;
  }

  char * qlinep = (char*) xmalloc(al+1);
  char * alinep = (char*) xmalloc(al+1);
  char * dlinep = (char*) xmalloc(al+1);

  *qline = qlinep;
  *aline = alinep;
  *dline = dlinep;

  long q_align_start = hits_list[i].align_q_start;
  long d_align_start = hits_list[i].align_d_start;

  d_strand = hits_list[i].dstrand;
  d_frame = hits_list[i].dframe;
  q_strand = hits_list[i].qstrand;
  q_frame = hits_list[i].qframe;
  
  if (symtype == 0)
  {
    sym = sym_ncbi_nt16;
    q_seq = query.nt[q_strand].seq;
    q_len = query.nt[q_strand].len;
  }
  else if (symtype == 5)
  {
    sym = sym_sound;
    q_seq = query.aa[0].seq;
    q_len = query.aa[0].len;
  }
  else
  {
    sym = sym_ncbi_aa;
    q_seq = query.aa[3*q_strand+q_frame].seq;
    q_len = query.aa[3*q_strand+q_frame].len;
    q_len_nt = query.nt[0].len;
    d_len_nt = hits_list[i].dlennt;
  }

  *identities = 0;
  *positives = 0;
  *indels = 0;
  *gaps = 0;
  *aligned = 0;
  
  d_seq = hits_list[i].dseq;
  d_len = hits_list[i].dlen;
  
  q_pos = q_align_start;
  d_pos = d_align_start;
  
  p = alignment;

  while(*p)
  {
    char op = *p++;
    long len;
    int n;
    sscanf(p, "%ld%n", & len, & n);
    p += n;
    
    *aligned += len;
    if (op == 'D')
    {
      for(long j=0; j<len; j++)
      {
	char qs = q_seq[q_pos++];
	*qlinep++ = sym[(int)(qs)];
	*alinep++ = ' ';
	*dlinep++ = '-';
      }
      *gaps += 1;
      *indels += len;
    }
    else if (op == 'I')
    {
      for(long j=0; j<len; j++)
      {
	char ds = d_seq[d_pos++];
	*qlinep++ = '-';
	*alinep++ = ' ';
	*dlinep++ = sym[(int)(ds)];
      }
      *gaps += 1;
      *indels += len;
    }
    else if (op == 'M')
    {
      for(long j=0; j<len; j++)
      {
	char qs = q_seq[q_pos++];
	char ds = d_seq[d_pos++];
	*qlinep++ = sym[(int)(qs)];
	if (qs == ds)
	{
	  *alinep++ = '|';
	  (*identities)++;
	  (*positives)++;
	}
	else if (score_matrix_63[32*qs+ds] > 0)
	{
	  *alinep++ = '+';
	  (*positives)++;
	}
	else
	{
	  *alinep++ = ' ';
	}
	*dlinep++ = sym[(int)(ds)];
      }
    }
    else
      fatal("Illegal alignment string.");
  }

  *qlinep = 0;
  *alinep = 0;
  *dlinep = 0;

  /* calculate first and last alignment positions for display */

  q_first = hits_list[i].align_q_start;
  q_last = hits_list[i].align_q_end;
  d_first = hits_list[i].align_d_start;
  d_last = hits_list[i].align_d_end;
  
  if (symtype == 0)
  {
    if (q_strand)
    {
      q_first = q_len - 1 - q_first;
      q_last = q_len - 1 - q_last;
    }

    if (d_strand)
    {
      d_first = d_len - 1 - d_first;
      d_last = d_len - 1 - d_last;
    }
  }
  
  if ((symtype == 2) || (symtype == 4))
  {
    if (q_strand)
    {
      q_first = query.nt[0].len - 1 - 3 * q_first - q_frame;
      q_last = query.nt[0].len - 1 - 3 * q_last - q_frame - 2;
    }
    else
    {
      q_first = 3 * q_first + q_frame;
      q_last = 3 * q_last + q_frame + 2;
    }
  }
  
  if ((symtype == 3) || (symtype == 4))
  {
    if (d_strand)
    {
      d_first = d_len_nt - 1 - 3 * d_first - d_frame;
      d_last = d_len_nt - 1 - 3 * d_last - d_frame - 2;
    }
    else
    {
      d_first = 3 * d_first + d_frame;
      d_last = 3 * d_last + d_frame + 2;
    }
  }

  q_first++;
  q_last++;
  d_first++;
  d_last++;

  long maxqpos = q_first > q_last ? q_first : q_last; 
  long maxdpos = d_first > d_last ? d_first : d_last; 
  long maxpos = maxqpos > maxdpos ? maxqpos : maxdpos;
  poswidth = 1;
  while (maxpos > 9)
  {
    maxpos /= 10;
    poswidth++;
  }
}

void count_align(long i,
		 long * identities,
		 long * positives,
		 long * indels,
		 long * aligned,
		 long * gaps)
{
  long q_align_start = hits_list[i].align_q_start;
  long d_align_start = hits_list[i].align_d_start;
  char * alignment = hits_list[i].alignment;

  d_strand = hits_list[i].dstrand;
  d_frame = hits_list[i].dframe;
  q_strand = hits_list[i].qstrand;
  q_frame = hits_list[i].qframe;
  
  if (symtype == 0)
  {
    sym = sym_ncbi_nt16;
    q_seq = query.nt[q_strand].seq;
    q_len = query.nt[q_strand].len;
  }
  else if (symtype == 5)
  {
    sym = sym_sound;
    q_seq = query.aa[0].seq;
    q_len = query.aa[0].len;
  }
  else
  {
    sym = sym_ncbi_aa;
    q_seq = query.aa[3*q_strand+q_frame].seq;
    q_len = query.aa[3*q_strand+q_frame].len;
    q_len_nt = query.nt[0].len;
    d_len_nt = hits_list[i].dlennt;
  }

  *identities = 0;
  *positives = 0;
  *indels = 0;
  *gaps = 0;
  *aligned = 0;
  
  d_seq = hits_list[i].dseq;
  d_len = hits_list[i].dlen;
  
  q_pos = q_align_start;
  d_pos = d_align_start;
  
  char * p = alignment;
  char * e = alignment + strlen(alignment);

  while(p < e)
  {
    char op = *p++;
    long len;
    int n;
    sscanf(p, "%ld%n", & len, & n);
    p += n;
    
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
	char qs = q_seq[q_pos++];
	char ds = d_seq[d_pos++];
	if (qs == ds)
	{
	  (*identities)++;
	  (*positives)++;
	}
	else if (score_matrix_63[32*qs+ds] > 0)
	  (*positives)++;
      }
    }
  }

  /* calculate first and last alignment positions for display */

  q_first = hits_list[i].align_q_start;
  q_last = hits_list[i].align_q_end;
  d_first = hits_list[i].align_d_start;
  d_last = hits_list[i].align_d_end;
  
  if (symtype == 0)
  {
    if (q_strand)
    {
      q_first = q_len - 1 - q_first;
      q_last = q_len - 1 - q_last;
    }

    if (d_strand)
    {
      d_first = d_len - 1 - d_first;
      d_last = d_len - 1 - d_last;
    }
  }
  
  if ((symtype == 2) || (symtype == 4))
  {
    if (q_strand)
    {
      q_first = query.nt[0].len - 1 - 3 * q_first - q_frame;
      q_last = query.nt[0].len - 1 - 3 * q_last - q_frame - 2;
    }
    else
    {
      q_first = 3 * q_first + q_frame;
      q_last = 3 * q_last + q_frame + 2;
    }
  }
  
  if ((symtype == 3) || (symtype == 4))
  {
    if (d_strand)
    {
      d_first = d_len_nt - 1 - 3 * d_first - d_frame;
      d_last = d_len_nt - 1 - 3 * d_last - d_frame - 2;
    }
    else
    {
      d_first = 3 * d_first + d_frame;
      d_last = 3 * d_last + d_frame + 2;
    }
  }

  q_first++;
  q_last++;
  d_first++;
  d_last++;

  long maxqpos = q_first > q_last ? q_first : q_last; 
  long maxdpos = d_first > d_last ? d_first : d_last; 
  long maxpos = maxqpos > maxdpos ? maxqpos : maxdpos;
  poswidth = 1;
  while (maxpos > 9)
  {
    maxpos /= 10;
    poswidth++;
  }
}

void hits_show_expect(double expect)
{
  char temp[10];
  if (expect < 1e-180)
    fprintf(out, "0.0  ");
  else if (expect < 9.5e-100)
  {
    sprintf(temp, "%-6.0e", expect);
    fputs(temp+1, out);
  }
  else if (expect < 0.00095)
    fprintf(out, "%-5.0e", expect);
  else if (expect < 0.0995)
    fprintf(out, "%-5.3f", expect);
  else if (expect < 0.95)
    fprintf(out, "%-5.2f", expect);
  else if (expect < 9.5)
    fprintf(out, "%-5.1f", expect);
  else
    fprintf(out, "%5.0f", expect);
}

void hits_show_expect_nospace(double expect)
{
  if (expect < 1e-180)
    fprintf(out, "0.0");
  else if (expect < 9.5e-100)
    fprintf(out, "%.0e", expect);
  else if (expect < 0.0995)
    fprintf(out, "%.3f", expect);
  else if (expect < 0.95)
    fprintf(out, "%.2f", expect);
  else if (expect < 9.5)
    fprintf(out, "%.1f", expect);
  else
    fprintf(out, "%.0f", expect);
}

void make_anchor(char * anchor, long size, long symtype, long queryno, long i)
{
  switch(symtype)
  {
  case 0:
    snprintf(anchor, size, "%ld_%ld__%c__+",
	     queryno,
	     hits_list[i].seqno,
	     hits_list[i].qstrand ? '-' : '+');
    break;
  case 2:
    snprintf(anchor, size, "%ld_%ld_%ld_%c__",
	     queryno,
	     hits_list[i].seqno,
	     hits_list[i].qframe+1,
	     hits_list[i].qstrand ? '-' : '+');
    break;
  case 3:
    snprintf(anchor, size, "%ld_%ld___%ld_%c",
	     queryno,
	     hits_list[i].seqno,
	     hits_list[i].dframe+1,
	     hits_list[i].dstrand ? '-' : '+');
    break;
  case 4:
    snprintf(anchor, size, "%ld_%ld_%ld_%c_%ld_%c",
	     queryno,
	     hits_list[i].seqno,
	     hits_list[i].qframe+1,
	     hits_list[i].qstrand ? '-' : '+',
	     hits_list[i].dframe+1,
	     hits_list[i].dstrand ? '-' : '+');
    break;
  default:
    snprintf(anchor, size, "%ld_%ld____",
	     queryno,
	     hits_list[i].seqno);
    break;
  }
}

void hits_defline_split(char * defline, 
			long * gi,
			char ** link, int * linklen, 
			char ** rest)
{
  char * p = defline;
  int len;

  *link = 0;
  *linklen = 0;
  *rest = 0;
  
  int m = sscanf(p, "gi|%ld%n", gi, & len);
  if (m > 0)
    //  if (len > 0)
    p += len;
  
  if (*p == '|')
    p++;
  
  char * r = strchr(p, ' ');
  if (r)
  {
    *linklen = r - p;
    *link = p;
    *rest = r+1;
  }
  else
  {
    * rest = p;
  }
}

void hits_show_xml_paralign(long showalignments,
			    long showhits,
			    struct db_thread_s * t)
{
  /* ParAlign XML */
  
  fprintf(out, "\t<paralignOutput>\n");
  
  const char * qseqtypedescr;
  struct sequence q;
  if ((query.symtype == 1) || (query.symtype == 3))
  {
    qseqtypedescr = "Amino Acid";
    q = query.aa[0];
  }
  else
  {
    qseqtypedescr = "Nucleotide";
    q = query.nt[0];
  }
  
  fprintf(out, "\t\t<queryInformation>\n");
  fprintf(out, "\t\t\t<queryFilename>./%s</queryFilename>\n", queryname);
  fprintf(out, "\t\t\t<querySequencetype>%s</querySequencetype>\n", qseqtypedescr);
  fprintf(out, "\t\t\t<queryDescription>%s</queryDescription>\n", query.description);
  fprintf(out, "\t\t\t<queryLength>%ld</queryLength>\n", q.len);
  fprintf(out, "\t\t\t<querySequence>");
  for(int i=0; i<q.len; i++)
    putc(query.sym[(int)(q.seq[i])], out);
  fprintf(out, "</querySequence>\n");
  fprintf(out, "\t\t</queryInformation>\n");
  
  const char * dbseqtypedescr;
  const char * ncbidb;
  const char * ncbiopt;
  if ((query.symtype == 0) || (query.symtype == 3) || (query.symtype == 4))
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
  fprintf(out, "\t\t<databaseInformation>\n");
  fprintf(out, "\t\t\t<databaseFilename>%s</databaseFilename>\n", databasename);
  fprintf(out, "\t\t\t<databaseSequencetype>%s</databaseSequencetype>\n", dbseqtypedescr);
  fprintf(out, "\t\t\t<databaseDescription>%s</databaseDescription>\n", db_gettitle());
  fprintf(out, "\t\t\t<databaseVersion>%ld</databaseVersion>\n", db_getversion());
  fprintf(out, "\t\t\t<databaseDate>%s</databaseDate>\n", db_gettime());
  fprintf(out, "\t\t\t<residueCount>%ld</residueCount>\n", db_getsymcount_masked());
  fprintf(out, "\t\t\t<sequenceCount>%ld</sequenceCount>\n", db_getseqcount_masked());
  fprintf(out, "\t\t\t<longestSequenceLength>%ld</longestSequenceLength>\n", db_getlongest());
  fprintf(out, "\t\t</databaseInformation>\n");
  
  const char * strands = "";
  switch(querystrands)
  {
  case 1:
    strands = "Plus";
    break;
  case 2:
    strands = "Minus";
    break;
  case 3:
    strands = "Both";
    break;
  }

  fprintf(out, "\t\t<options>\n");
  fprintf(out, "\t\t\t<algorithm>Smith-Waterman</algorithm>\n");

  if ((symtype == 0) || (symtype == 2) || (symtype == 4))
  {
    fprintf(out, "\t\t\t<queryStrands>%s</queryStrands>\n", strands);
  }

  if (symtype == 0)
    fprintf(out, "\t\t\t<scoreMatrix>NT</scoreMatrix>\n");
  else
    fprintf(out, "\t\t\t<scoreMatrix>%s</scoreMatrix>\n", matrixname);

  fprintf(out, "\t\t\t<gapPenalties>\n");
  fprintf(out, "\t\t\t\t<gapPenaltyOpen>%ld</gapPenaltyOpen>\n", gapopen);
  fprintf(out, "\t\t\t\t<gapPenaltyExtension>%ld</gapPenaltyExtension>\n", gapextend);
  fprintf(out, "\t\t\t\t<ungapped>\n");
  fprintf(out, "\t\t\t\t\t<ungappedLambda>%.4g</ungappedLambda>\n", lambda);
  fprintf(out, "\t\t\t\t\t<ungappedKappa>%.4g</ungappedKappa>\n", K);
  fprintf(out, "\t\t\t\t\t<ungappedEta>%.4g</ungappedEta>\n", H);
  fprintf(out, "\t\t\t\t</ungapped>\n");
  fprintf(out, "\t\t\t\t<gapped>\n");
  fprintf(out, "\t\t\t\t\t<gappedLambda>%.4g</gappedLambda>\n", lambda);
  fprintf(out, "\t\t\t\t\t<gappedKappa>%.4g</gappedKappa>\n", K);
  fprintf(out, "\t\t\t\t\t<gappedEta>%.4g</gappedEta>\n", H);
  fprintf(out, "\t\t\t\t</gapped>\n");

  fprintf(out, "\t\t\t</gapPenalties>\n");
  fprintf(out, "\t\t\t<expectRange>\n");
  fprintf(out, "\t\t\t\t<expectRangeFrom>%.2g</expectRangeFrom>\n", minexpect);
  fprintf(out, "\t\t\t\t<expectRangeTo>%.2g</expectRangeTo>\n", expect);
  fprintf(out, "\t\t\t</expectRange>\n");
  fprintf(out, "\t\t\t<displayLimits>\n");
  fprintf(out, "\t\t\t\t<hitLimit>%ld</hitLimit>\n", maxmatches);
  fprintf(out, "\t\t\t\t<alignmentLimit>%ld</alignmentLimit>\n", alignments);
  fprintf(out, "\t\t\t\t<subalignmentLimit>%ld</subalignmentLimit>\n", (long)1);
  fprintf(out, "\t\t\t</displayLimits>\n");
  fprintf(out, "\t\t\t<threads>%ld</threads>\n", threads);
  fprintf(out, "\t\t</options>\n");

  fprintf(out, "\t\t\t<searchInformation>\n");
  fprintf(out, "\t\t\t\t<searchStarted>%s</searchStarted>\n", ti.starttime);
  fprintf(out, "\t\t\t\t<searchCompleted>%s</searchCompleted>\n", ti.endtime);
  fprintf(out, "\t\t\t\t<searchElapsedTime>%.2fs</searchElapsedTime>\n", ti.elapsed);
  fprintf(out, "\t\t\t\t<searchSpeed>%.3f GCUPS</searchSpeed>\n", ti.speed / 1e9);
  fprintf(out, "\t\t\t\t<searchSWAlignments>\n");
  fprintf(out, "\t\t\t\t\t<SWAbsolute>%ld</SWAbsolute>\n", compute7);
  fprintf(out, "\t\t\t\t\t<SWPercent>100</SWPercent>\n");
  fprintf(out, "\t\t\t\t</searchSWAlignments>\n");
  fprintf(out, "\t\t\t</searchInformation>\n");

  fprintf(out, "\t\t<resultInformation>\n");
  fprintf(out, "\t\t\t<resultHits>\n");
  fprintf(out, "\t\t\t\t<totalCount>%ld</totalCount>\n", totalhits);
  fprintf(out, "\t\t\t\t<obviousCount>%ld</obviousCount>\n", obvious);
  fprintf(out, "\t\t\t\t<shownCount>%ld</shownCount>\n", showhits);
  fprintf(out, "\t\t\t</resultHits>\n");
  fprintf(out, "\t\t\t<alignmentCount>%ld</alignmentCount>\n", showalignments);
  fprintf(out, "\t\t</resultInformation>\n");
  
  fprintf(out, "\t\t<shortVersionHits>\n");
  
  for(long i=0; i<showhits; i++)
  {
    long score = hits_list[i].score;
    double e = Kmn * exp(- lambda * score);

    char anchor[200];
    make_anchor(anchor, 200, query.symtype, queryno, i);

    long deflines;
    char ** deflinetable;
    long gi = 0;
    char * link;
    char * title;
    int linklen;
    db_parse_header(t, hits_list[i].header_address, hits_list[i].header_length,
		    1, & deflines, & deflinetable);
    hits_defline_split(deflinetable[0], 
		       & gi,
		       & link, & linklen,
		       & title);

    fprintf(out, "\t\t\t<shortVersionHit>\n");
    fprintf(out, "\t\t\t\t<shortVersionAnchor>%s</shortVersionAnchor>\n", anchor);
    if (gi)
      {
    fprintf(out, "\t\t\t\t<shortVersionLink>\n");
    fprintf(out, "\t\t\t\t\t<shortVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=%s&amp;list_uids=%ld&amp;dopt=%s</shortVersionLinkDestination>\n", ncbidb, gi, ncbiopt);
    fprintf(out, "\t\t\t\t\t<shortVersionLinkText>gi|%ld</shortVersionLinkText>\n", gi);
    fprintf(out, "\t\t\t\t</shortVersionLink>\n");
      }
    fprintf(out, "\t\t\t\t<shortVersionLink>\n");
    fprintf(out, "\t\t\t\t\t<shortVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&amp;db=%s&amp;term=%.*s&amp;doptcmdl=%s</shortVersionLinkDestination>\n", ncbidb, linklen, link, ncbiopt);
    fprintf(out, "\t\t\t\t\t<shortVersionLinkText>%.*s</shortVersionLinkText>\n", linklen, link);
    fprintf(out, "\t\t\t\t</shortVersionLink>\n");
    fprintf(out, "\t\t\t\t<shortVersionName>%.35s</shortVersionName>\n", title);
    if (symtype == 0)
    {
      fprintf(out, "\t\t\t\t<shortVersionStrand>%c</shortVersionStrand>\n", hits_list[i].qstrand ? '-' : '+');
    }
    else if (symtype == 2)
    {
      fprintf(out, "\t\t\t\t<shortVersionFrame>%c%ld</shortVersionFrame>\n", 
	     hits_list[i].qstrand ? '-' : '+', 
	     hits_list[i].qframe+1);
    }
    else if (symtype == 3)
    {
      fprintf(out, "\t\t\t\t<shortVersionFrame>%c%ld</shortVersionFrame>\n", 
	     hits_list[i].dstrand ? '-' : '+', 
	     hits_list[i].dframe+1);
    }
    else if (symtype == 4)
    {
      fprintf(out, "\t\t\t\t<shortVersionFrame>%c%ld/%c%ld</shortVersionFrame>\n", 
	     hits_list[i].qstrand ? '-' : '+', 
	     hits_list[i].qframe+1,
	     hits_list[i].dstrand ? '-' : '+', 
	     hits_list[i].dframe+1);
    }
    fprintf(out, "\t\t\t\t<shortVersionScore>%ld</shortVersionScore>\n", score);
    fprintf(out, "\t\t\t\t<shortVersionEValue>%.2g</shortVersionEValue>\n", e);
    fprintf(out, "\t\t\t</shortVersionHit>\n");

    for (int d=0; d<deflines; d++)
      free(deflinetable[d]);
    free(deflinetable);
  }

  fprintf(out, "\t\t</shortVersionHits>\n");

  if (showalignments)
  {
    fprintf(out, "\t\t<longVersionHits>\n");
    
    for(long i=0; i<showalignments; i++)
    {
      
      char anchor[200];
      make_anchor(anchor, 200, query.symtype, queryno, i);
      
      fprintf(out, "\t\t\t<longVersionHit>\n");
      fprintf(out, "\t\t\t\t<longVersionAnchor>%s</longVersionAnchor>\n", anchor);
      
      long deflines;
      char ** deflinetable;
      long gi = 0;
      char * link;
      char * title;
      int linklen;
      db_parse_header(t, hits_list[i].header_address, hits_list[i].header_length,
		      1, & deflines, & deflinetable);
      fprintf(out, "\t\t\t\t<linkContainer>\n");
      
      for (int d=0; d < deflines; d++)
      {
	hits_defline_split(deflinetable[d], 
			   & gi,
			   & link, & linklen,
			   & title);
  
        if (gi)
	{
          fprintf(out, "\t\t\t\t\t<longVersionLink>\n");
	  fprintf(out, "\t\t\t\t\t\t<longVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Retrieve&amp;db=%s&amp;list_uids=%ld&amp;dopt=%s</longVersionLinkDestination>\n", ncbidb, gi, ncbiopt);
	  fprintf(out, "\t\t\t\t\t\t<longVersionLinkText>gi|%ld</longVersionLinkText>\n", gi);
	  fprintf(out, "\t\t\t\t\t</longVersionLink>\n");
	}
      
	fprintf(out, "\t\t\t\t\t<longVersionLink>\n");
	fprintf(out, "\t\t\t\t\t\t<longVersionLinkDestination>http://www.ncbi.nlm.nih.gov/entrez/query.fcgi?cmd=Search&amp;db=%s&amp;term=%.*s&amp;doptcmdl=%s</longVersionLinkDestination>\n", ncbidb, linklen, link, ncbiopt);
	fprintf(out, "\t\t\t\t\t\t<longVersionLinkText>%.*s</longVersionLinkText>\n", linklen, link);
	fprintf(out, "\t\t\t\t\t</longVersionLink>\n");
      
	fprintf(out, "\t\t\t\t\t<longVersionName>%s</longVersionName>\n", title);
      
	free(deflinetable[d]);
      }
    
      free(deflinetable);
        
      fprintf(out, "\t\t\t\t</linkContainer>\n");
    
      long dlen = hits_list[i].dlen;
      long dlennt = hits_list[i].dlennt;

      if (symtype == 0)
	fprintf(out, "\t\t\t\t<databaseSequenceLength>%ld nt</databaseSequenceLength>\n", dlen);
      else if ((symtype == 3) || (symtype == 4))
	fprintf(out, "\t\t\t\t<databaseSequenceLength>%ld nt</databaseSequenceLength>\n", dlennt);
      else
	fprintf(out, "\t\t\t\t<databaseSequenceLength>%ld aa</databaseSequenceLength>\n", dlen);
      
      if (symtype == 0)
      {
	fprintf(out, "\t\t\t\t<alignmentMatchLocation>%s</alignmentMatchLocation>\n", hits_list[i].qstrand ? "Matches on complementary strands." : "Matches on same strands.");
      }
      else if ((symtype>=2) && (symtype<=4))
      {
	fprintf(out, "\t\t\t\t<longVersionFrames>\n");

	if ((symtype == 2) || (symtype == 4))
	{
	  fprintf(out, "\t\t\t\t\t<longVersionQueryFrame>\n");
	  fprintf(out, "\t\t\t\t\t\t<queryStrand>%c</queryStrand>\n", hits_list[i].qstrand ? '-' : '+');
	  fprintf(out, "\t\t\t\t\t\t<queryFrame>%ld</queryFrame>\n", hits_list[i].qframe+1);
	  fprintf(out, "\t\t\t\t\t</longVersionQueryFrame>\n");
	}
	
	if ((symtype == 3) || (symtype == 4))
	{
	  fprintf(out, "\t\t\t\t\t<longVersionDatabaseFrame>\n");
	  fprintf(out, "\t\t\t\t\t\t<databaseStrand>%c</databaseStrand>\n", hits_list[i].dstrand ? '-' : '+');
	  fprintf(out, "\t\t\t\t\t\t<databaseFrame>%ld</databaseFrame>\n", hits_list[i].dframe+1);
	  fprintf(out, "\t\t\t\t\t</longVersionDatabaseFrame>\n");
	}

	fprintf(out, "\t\t\t\t</longVersionFrames>\n");
      }

      long score = hits_list[i].score;
      double e = Kmn * exp(- lambda * score);

      long identities;
      long positives;
      long indels;
      long gaps;
      long aligned;
    
      char *qline;
      char *aline;
      char *dline;
        
      whole_align(i, & identities, & positives, & indels, & aligned, & gaps,
		  & qline, & aline, & dline);

      fprintf(out, "\t\t\t\t<alignment>\n");
      fprintf(out, "\t\t\t\t\t<subalignment>\n");
      fprintf(out, "\t\t\t\t\t\t<longVersionScore>%ld</longVersionScore>\n", score);
      fprintf(out, "\t\t\t\t\t\t<longVersionEValue>%.2g</longVersionEValue>\n", e);
      fprintf(out, "\t\t\t\t\t\t<identical>\n");
      fprintf(out, "\t\t\t\t\t\t\t<identicalNominator>%ld</identicalNominator>\n", identities);
      fprintf(out, "\t\t\t\t\t\t\t<identicalDenominator>%ld</identicalDenominator>\n", aligned);
      fprintf(out, "\t\t\t\t\t\t\t<identicalPercentage>%.1f</identicalPercentage>\n", 100.0*identities/aligned);
      fprintf(out, "\t\t\t\t\t\t</identical>\n");

      if (symtype != 0)
      {
	fprintf(out, "\t\t\t\t\t\t<positive>\n");
	fprintf(out, "\t\t\t\t\t\t\t<positiveNominator>%ld</positiveNominator>\n", positives);
	fprintf(out, "\t\t\t\t\t\t\t<positiveDenominator>%ld</positiveDenominator>\n", aligned);
	fprintf(out, "\t\t\t\t\t\t\t<positivePercentage>%.1f</positivePercentage>\n", 100.0*positives/aligned);
	fprintf(out, "\t\t\t\t\t\t</positive>\n");
      }

      fprintf(out, "\t\t\t\t\t\t<indels>\n");
      fprintf(out, "\t\t\t\t\t\t\t<indelsNominator>%ld</indelsNominator>\n", indels);
      fprintf(out, "\t\t\t\t\t\t\t<indelsDenominator>%ld</indelsDenominator>\n", aligned);
      fprintf(out, "\t\t\t\t\t\t\t<indelsPercentage>%.1f</indelsPercentage>\n", 100.0*indels/aligned);
      fprintf(out, "\t\t\t\t\t\t</indels>\n");
      fprintf(out, "\t\t\t\t\t\t<gaps>%ld</gaps>\n", gaps);
      fprintf(out, "\t\t\t\t\t\t<alignmentQuery>\n");
      fprintf(out, "\t\t\t\t\t\t\t<alignmentQueryStart>%ld</alignmentQueryStart>\n", q_first);
      fprintf(out, "\t\t\t\t\t\t\t<alignmentQueryLine>%s</alignmentQueryLine>\n", qline);
      fprintf(out, "\t\t\t\t\t\t\t<alignmentQueryEnd>%ld</alignmentQueryEnd>\n", q_last);
      fprintf(out, "\t\t\t\t\t\t</alignmentQuery>\n");
      fprintf(out, "\t\t\t\t\t\t<alignmentLine>%s</alignmentLine>\n", aline);
      fprintf(out, "\t\t\t\t\t\t<alignmentDatabase>\n");
      fprintf(out, "\t\t\t\t\t\t\t<alignmentDatabaseStart>%ld</alignmentDatabaseStart>\n", d_first);
      fprintf(out, "\t\t\t\t\t\t\t<alignmentDatabaseLine>%s</alignmentDatabaseLine>\n", dline);
      fprintf(out, "\t\t\t\t\t\t\t<alignmentDatabaseEnd>%ld</alignmentDatabaseEnd>\n", d_last);
      fprintf(out, "\t\t\t\t\t\t</alignmentDatabase>\n");
      fprintf(out, "\t\t\t\t\t</subalignment>\n");
      fprintf(out, "\t\t\t\t</alignment>\n");
      fprintf(out, "\t\t\t</longVersionHit>\n");

      free(qline);
      free(aline);
      free(dline);

    }

    fprintf(out, "\t\t</longVersionHits>\n");
  }

  fprintf(out, "\t</paralignOutput>\n");
}

static void show_description(const char *desc)
{
  const char *dptr;

  for (dptr = desc; *dptr != '\0' && *dptr != ' '; dptr++)
  {
    putc(*dptr, out);
  }
}

void hits_show_xml(long show_gis,
		   long showalignments,
		   long showhits,
		   struct db_thread_s * t)
{
  /* Simple XML */
  
  fprintf(out, "<result>\n");
  fprintf(out, "  <general>\n");
  fprintf(out, "    <hitcount>%d</hitcount>\n", hits_count);
  fprintf(out, "  </general>\n");
  fprintf(out, "  <hits>\n");
  
  for(long i=0; i<showhits; i++)
  {
    long seqno = hits_list[i].seqno;
    long score = hits_list[i].score;
    long dlen = hits_list[i].dlen;
    
    fprintf(out, "    <hit>\n");
    fprintf(out, "      <hitno>%ld</hitno>\n", i+1);
    fprintf(out, "      <track>%ld</track>\n", seqno);
    fprintf(out, "      <query>");
    show_description(query.description);
    fprintf(out,"</query>\n");
    fprintf(out, "      <name>");
    db_showheader(t, hits_list[i].header_address,
		  hits_list[i].header_length,
		  show_gis, 0, 0, LONG_MAX, 1, 1);
    fprintf(out, "</name>\n");
    fprintf(out, "      <len>%ld</len>\n", dlen);
    fprintf(out, "      <score>%ld</score>\n", score);
    
    if (i < showalignments)
    {
      long identities;
      long positives;
      long gaps;
      long aligned;
      long indels;

      char *qline;
      char *aline;
      char *dline;
        
      whole_align(i, & identities, & positives, & indels, & aligned, & gaps,
		  & qline, & aline, & dline);

      fprintf(out, "      <alignment>");
      fprintf(out, "%s", hits_list[i].alignment);
      fprintf(out, "</alignment>\n");

      fprintf(out, "      <qpos>%ld,%ld</qpos>\n", q_first, q_last);
      fprintf(out, "      <dpos>%ld,%ld</dpos>\n", d_first, d_last);
      
      fprintf(out, "      <qseq>%s</qseq>\n", qline);
      fprintf(out, "      <aseq>%s</aseq>\n", aline);
      fprintf(out, "      <dseq>%s</dseq>\n", dline);

      free(qline);
      free(aline);
      free(dline);
    }
    fprintf(out, "    </hit>\n");
  }
  fprintf(out, "  </hits>\n");
  fprintf(out, "</result>\n");
}

void hits_show_tsv(long showalignments,
		   long showcomments,
		   struct db_thread_s * t)
{
  char title[] = "SWIPE " SWIPE_VERSION;
  char ref[] = "Reference: T. Rognes (2011) Faster Smith-Waterman database searches with inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.";
  
  if (showcomments)
    {
      fprintf(out, "# %s - Compiled %s %s - %s\n", title, __DATE__, __TIME__, ref);
      fprintf(out, "# Query: %s\n", query.description);
      fprintf(out, "# Database: %s\n", databasename);
      if (stats_available)
	fprintf(out, "# Fields: Query id, Subject id, %% identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, e-value, bit score\n");
      else
	fprintf(out, "# Fields: Query id, Subject id, %% identity, alignment length, mismatches, gap openings, q. start, q. end, s. start, s. end, score\n");
    }

  for(long i=0; i<showalignments; i++)
  {
    show_description(query.description);
    putc('\t', out);
    db_showheader(t, hits_list[i].header_address,
		  hits_list[i].header_length,
		  1, 0, 0, LONG_MAX, 1, 0);
    
    long identities;
    long positives;
    long gaps;
    long aligned;
    long indels;
    
    count_align(i, & identities, & positives, & indels, & aligned, & gaps);
    
    long score = hits_list[i].score;
    
    fprintf(out, "\t%.2f\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld\t%ld", 
	    100.0 * identities / aligned,
	    aligned,
	    aligned - identities - indels,
	    gaps,
	    q_first,
	    q_last,
	    d_first,
	    d_last);
    
    if (stats_available)
    {
      double expect = Kmn * exp(- lambda * score);
      fprintf(out, "\t%.2g", expect);
      double bits = lambda_d_log2 * score - logK_d_log2;
      fprintf(out, "\t%.1f", bits);
    }
    else
    {
      fprintf(out, "\t%ld", score);
    }

    fprintf(out, "\n");
  }
}

void hits_show_plain(long show_gis,
		     long showalignments,
		     long showhits,
		     struct db_thread_s * t)
{
    if (hits_count == 0)
    {
      fprintf(out, "\nNo hits.\n");
    }
    else
    {
      if (stats_available)
      {
	fprintf(out, "                                                                 Score    E\n");
	fprintf(out, "Sequences producing significant alignments:                      (bits) Value\n\n");
      }
      else
      {
	fprintf(out, "Sequences producing significant alignments:                         Score\n\n");
      }
	  
      for(long i=0; i<showhits; i++)
      {
	long headerlen = 67;
	if (symtype == 0)
	  headerlen = 65;
	else if ((symtype == 2) || (symtype == 3))
	  headerlen = 64;
	else if (symtype == 4)
	  headerlen = 61;
	      
	db_showheader(t, 
		      hits_list[i].header_address,
		      hits_list[i].header_length, 
		      show_gis, 0, headerlen, headerlen, 1, 1);

	long score = hits_list[i].score;
	      
	if (symtype == 0)
	  fprintf(out, " %c", hits_list[i].dstrand ? '-' : '+');
	else if (symtype == 2)
	  fprintf(out, " %c%ld", hits_list[i].qstrand ? '-' : '+',
		 hits_list[i].qframe+1);
	else if (symtype == 3)
	  fprintf(out, " %c%ld", hits_list[i].dstrand ? '-' : '+',
		 hits_list[i].dframe+1);
	else if (symtype == 4)
	  fprintf(out, " %c%ld/%c%ld", 
		 hits_list[i].qstrand ? '-' : '+',
		 hits_list[i].qframe+1,
		 hits_list[i].dstrand ? '-' : '+',
		 hits_list[i].dframe+1);
	      
	if (stats_available)
	{
	  long bits = (long) floor(lambda_d_log2 * score - logK_d_log2 + 0.5);
	  double expect = Kmn * exp(- lambda * score);
		
	  fprintf(out, " %5ld", bits);
		
	  fprintf(out, "   ");
		
	  hits_show_expect(expect);
	}
	else
	  fprintf(out, " %5ld", score);

	putc('\n', out);
      }

#if 0
      fprintf(out, "\n");
      if (showalignments)
	fprintf(out, "Alignments (%ld)\n", showalignments);
#endif

      for(long i=0; i<showalignments; i++)
      {
	fprintf(out, "\n");
	db_showheader(t, hits_list[i].header_address,
		      hits_list[i].header_length,
		      show_gis, 10, 0, 79, LONG_MAX, 1);
	if ((symtype == 3) || (symtype == 4))
	  fprintf(out, "          Length = %ld\n", hits_list[i].dlennt);
	else
	  fprintf(out, "          Length = %ld\n", hits_list[i].dlen);
	fprintf(out, "\n");
	      
	long score = hits_list[i].score;

	if (stats_available)
	{
	  double bits = lambda_d_log2 * score - logK_d_log2;
	  double expect = Kmn * exp(- lambda * score);
		
	  fprintf(out, " Score = %.1lf bits (%ld), Expect = ", bits, score);
	  hits_show_expect(expect);
	}
	else
	{
	  fprintf(out, " Score = %ld", score);
	}

	putc('\n', out);

	long identities;
	long positives;
	long gaps;
	long aligned;
	long indels;

	count_align(i, & identities, & positives, & indels, & aligned, & gaps);
	      
	fprintf(out, " Identities = %ld/%ld (%ld%%)",
	       identities, aligned, identities * 100 / aligned);
	if (symtype > 0)
	  fprintf(out, ", Positives = %ld/%ld (%ld%%)",
		 positives, aligned, positives * 100 / aligned);
	if (indels)
	  fprintf(out, ", Gaps = %ld/%ld (%ld%%)", indels, aligned, indels * 100 / aligned);
	fprintf(out, "\n");

	if (symtype == 0)
	  fprintf(out, " Strand = %s\n", hits_list[i].dstrand ? "Plus / Minus" : "Plus / Plus");
	else if (symtype == 2)
	  fprintf(out, " Frame = %c%ld\n", hits_list[i].qstrand ? '-':'+', hits_list[i].qframe+1);
	else if (symtype == 3)
	  fprintf(out, " Frame = %c%ld\n", hits_list[i].dstrand ? '-':'+', hits_list[i].dframe+1);
	else if (symtype == 4)
	  fprintf(out, " Frame = %c%ld / %c%ld\n", 
		 hits_list[i].qstrand ? '-' : '+',
		 hits_list[i].qframe+1,
		 hits_list[i].dstrand ? '-' : '+',
		 hits_list[i].dframe+1);

#if 0
	// fprintf(out, "String: %s\n", hits_list[i].alignment);

	fprintf(out, "\nAlignment end: %ld, %ld\n", 
	       hits_list[i].bestq+1,
	       hits_list[i].align_hint+1);
	fprintf(out, "Hint offset: %ld, %ld\n", 
	       hits_list[i].align_q_end - hits_list[i].bestq,
	       hits_list[i].align_d_end - hits_list[i].align_hint);

	fprintf(out, "Alignment score: %ld\n",hits_list[i].score_align); 
#endif

	show_align(i);
	fprintf(out, "\n");
      }
	  
    }
    //      fprintf(out, "\n");
}

void hits_show_begin(long view)
{
  if (view==0)
    {
      fprintf(out, "%s [%s %s]\n\n%s\n\n", 
	      "SWIPE " SWIPE_VERSION, 
	      __DATE__, 
	      __TIME__, 
	      "Reference: T. Rognes (2011) Faster Smith-Waterman database searches\nwith inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.");
    }
  else if (view==7)
    {
      fprintf(out, "<?xml version=\"1.0\"?>\n");
    }
  else if (view==99)
    {
      char url1[] = "http://www.w3.org/2001/XMLSchema-instance";
      char url2[] = "http://www.paralign.org/ParalignXML.xsd";

      fprintf(out, "<?xml version=\"1.0\"?>\n");
      fprintf(out, "<ParalignXML xmlns:xsi=\"%s\" xsi:noNamespaceSchemaLocation=\"%s\">\n",
	      url1, url2);
      fprintf(out, "\t<programInformation>\n");
      fprintf(out, "\t\t<programName>swipe</programName>\n");
      fprintf(out, "\t\t<programVersion>SWIPE " SWIPE_VERSION "</programVersion>\n");
      fprintf(out, "\t\t<programDescription>Smith-Waterman database searches with inter-sequence SIMD parallelisation</programDescription>\n");
      fprintf(out, "\t\t<articleReferences>\n");
      fprintf(out, "\t\t\t<reference>T. Rognes (2011) Faster Smith-Waterman database searches with inter-sequence SIMD parallelisation, BMC Bioinformatics, 12:221.</reference>\n");
      fprintf(out, "\t\t</articleReferences>\n");
      fprintf(out, "\t\t<license>SWIPE is available under the GNU Affero General Public License, version 3</license>\n");
      fprintf(out, "\t</programInformation>\n");
    }
}

void hits_show_end(long view)
{
  if (view==99)
  {
    fprintf(out, "</ParalignXML>\n");
  }
}

void hits_show(long view, long show_gis)
{
  // compute number of hits and alignments to actually show

  long showalignments;
  long showhits;

  if (hits_count < opt_descriptions)
    showhits = hits_count;
  else
    showhits = opt_descriptions;

  if (hits_count < opt_alignments)
    showalignments = hits_count;
  else
    showalignments = opt_alignments;

  long Nbest = 1;
  long place = 1;
  while (place < hits_count)
  {
    if (hits_list[place-1].score > hits_list[place].score)
      Nbest++;

    if (Nbest <= opt_show_best)
      place++;
    else
      break;
  }

  if (showhits > place)
    showhits = place;
  
  if (showalignments > place)
    showalignments = place;

  struct db_thread_s * t = db_thread_create();

  if(view == 0)
  {
    hits_show_plain(show_gis, showalignments, showhits, t);
  }
  else if (view==7)
  {
    hits_show_xml(show_gis, showalignments, showhits, t);
  }
  else if ((view==8)||(view==9))
  {
    hits_show_tsv(showalignments, view == 9, t);
  }
  else if (view==99)
  {
    hits_show_xml_paralign(showalignments, showhits, t);
  }
  db_thread_destruct(t);
}

