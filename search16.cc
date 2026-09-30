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

constexpr std::size_t CHANNELS = channels_16;
constexpr std::size_t CDEPTH = 4;

// the word 0x8000 (the lanes of _mm_set_epi16() are short: 0x8000
// does not fit in a signed short, -32768 has the same bits)
constexpr short word_0x8000 = static_cast<short>(-32768);

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

// C++26 refactoring: std::simd, with std::add_sat and std::sub_sat

// One cell of a block (the ONESTEP macro of the former inline
// assembly): H is the score of the diagonal cell, N receives the score
// of this cell (the diagonal of the next column), F and E are the
// vertical and horizontal gap scores, S is the running maximum
inline auto onestep16(__m128i const H,
                      __m128i & N,
                      __m128i & F,
                      __m128i const V,
                      __m128i & E,
                      __m128i & S,
                      __m128i const Q,
                      __m128i const R) -> void
{
  auto cell = _mm_adds_epi16(H, V);
  cell = _mm_max_epi16(cell, F);
  cell = _mm_max_epi16(cell, E);
  S = _mm_max_epi16(cell, S);
  F = _mm_subs_epi16(F, R);
  E = _mm_subs_epi16(E, R);
  N = cell;
  cell = _mm_subs_epi16(cell, Q);
  E = _mm_max_epi16(cell, E);
  F = _mm_max_epi16(cell, F);
}

inline auto donormal16(volatile __m128i * Sm,
                       __m128i * hep,
                       __m128i * const * qp,
                       __m128i const * Qm,
                       __m128i const * Rm,
                       long ql,
                       __m128i const * Zm) -> void
{
  auto S = *Sm;
  auto const Q = *Qm;
  auto const R = *Rm;
  auto H0 = *Zm;
  auto H1 = H0;
  auto H2 = H0;
  auto H3 = H0;
  auto F0 = H0;
  auto F1 = H0;
  auto F2 = H0;
  auto F3 = H0;
  __m128i N1;
  __m128i N2;
  __m128i N3;

  for (long qi = 0; qi < ql; ++qi)
  {
    __m128i const * const x = qp[qi];  // load x from qp[qi]
    auto const N0 = hep[2 * qi];  // load N0
    auto E = hep[(2 * qi) + 1];  // load E

    onestep16(H0, N1, F0, x[0], E, S, Q, R);
    onestep16(H1, N2, F1, x[1], E, S, Q, R);
    onestep16(H2, N3, F2, x[2], E, S, Q, R);
    onestep16(H3, hep[2 * qi], F3, x[3], E, S, Q, R);

    hep[(2 * qi) + 1] = E;  // save E
    H0 = N0;
    H1 = N1;
    H2 = N2;
    H3 = N3;
  }

  *Sm = S;  // save S
}

inline auto domasked16(volatile __m128i * Sm,
                       __m128i * hep,
                       __m128i * const * qp,
                       __m128i const * Qm,
                       __m128i const * Rm,
                       long ql,
                       __m128i const * Zm,
                       __m128i const * Mm) -> void
{
  auto const M = *Mm;
  auto S = _mm_adds_epi16(_mm_adds_epi16(*Sm, M), M);  // add M
  auto const Q = *Qm;
  auto const R = *Rm;
  auto H0 = *Zm;
  auto H1 = H0;
  auto H2 = H0;
  auto H3 = H0;
  auto F0 = H0;
  auto F1 = H0;
  auto F2 = H0;
  auto F3 = H0;
  __m128i N1;
  __m128i N2;
  __m128i N3;

  for (long qi = 0; qi < ql; ++qi)
  {
    __m128i const * const x = qp[qi];  // load x from qp[qi]
    auto const N0 = _mm_adds_epi16(_mm_adds_epi16(hep[2 * qi], M), M);  // load N0, add M
    auto E = _mm_adds_epi16(_mm_adds_epi16(hep[(2 * qi) + 1], M), M);  // load E, add M

    onestep16(H0, N1, F0, x[0], E, S, Q, R);
    onestep16(H1, N2, F1, x[1], E, S, Q, R);
    onestep16(H2, N3, F2, x[2], E, S, Q, R);
    onestep16(H3, hep[2 * qi], F3, x[3], E, S, Q, R);

    hep[(2 * qi) + 1] = E;  // save E
    H0 = N0;
    H1 = N1;
    H2 = N2;
    H3 = N3;
  }

  *Sm = S;  // save S
}

inline auto dprofile_fill16(WORD * dprofile_word,
			    WORD * score_matrix_word,
			    BYTE const * dseq) -> void
{
  __m128i xmm0;
  __m128i xmm1;
  __m128i xmm2;
  __m128i xmm3;
  __m128i xmm4;
  __m128i xmm5;
  __m128i xmm6;
  __m128i xmm7;
  __m128i xmm8;
  __m128i xmm9;
  __m128i xmm10;
  __m128i xmm11;
  __m128i xmm12;
  __m128i xmm13;
  __m128i xmm14;
  __m128i xmm15;
  __m128i xmm16;
  __m128i xmm17;
  __m128i xmm18;
  __m128i xmm19;
  __m128i xmm20;
  __m128i xmm21;
  __m128i xmm22;
  __m128i xmm23;
  __m128i xmm24;
  __m128i xmm25;
  __m128i xmm26;
  __m128i xmm27;
  __m128i xmm28;
  __m128i xmm29;
  __m128i xmm30;
  __m128i xmm31;
  
  for (std::size_t j = 0; j < CDEPTH; j++)
  {
    std::array<int, CHANNELS> d;
    for (std::size_t z = 0; z < CHANNELS; z++)
    {
      d[z] = dseq[(j * CHANNELS) + z] << 5;
    }

    //      for(int i=0; i<24; i += 8)
    for(std::size_t i=0; i<32; i += 8)
    {
      xmm0  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix_word + d[0] + i));
      xmm1  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix_word + d[1] + i));
      xmm2  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix_word + d[2] + i));
      xmm3  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix_word + d[3] + i));
      xmm4  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix_word + d[4] + i));
      xmm5  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix_word + d[5] + i));
      xmm6  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix_word + d[6] + i));
      xmm7  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix_word + d[7] + i));
      
      xmm8  = _mm_unpacklo_epi16(xmm0,  xmm1);
      xmm9  = _mm_unpackhi_epi16(xmm0,  xmm1);
      xmm10 = _mm_unpacklo_epi16(xmm2,  xmm3);
      xmm11 = _mm_unpackhi_epi16(xmm2,  xmm3);
      xmm12 = _mm_unpacklo_epi16(xmm4,  xmm5);
      xmm13 = _mm_unpackhi_epi16(xmm4,  xmm5);
      xmm14 = _mm_unpacklo_epi16(xmm6,  xmm7);
      xmm15 = _mm_unpackhi_epi16(xmm6,  xmm7);
      
      xmm16 = _mm_unpacklo_epi32(xmm8,  xmm10);
      xmm17 = _mm_unpackhi_epi32(xmm8,  xmm10);
      xmm18 = _mm_unpacklo_epi32(xmm12, xmm14);
      xmm19 = _mm_unpackhi_epi32(xmm12, xmm14);
      xmm20 = _mm_unpacklo_epi32(xmm9,  xmm11);
      xmm21 = _mm_unpackhi_epi32(xmm9,  xmm11);
      xmm22 = _mm_unpacklo_epi32(xmm13, xmm15);
      xmm23 = _mm_unpackhi_epi32(xmm13, xmm15);
      
      xmm24 = _mm_unpacklo_epi64(xmm16, xmm18);
      xmm25 = _mm_unpackhi_epi64(xmm16, xmm18);
      xmm26 = _mm_unpacklo_epi64(xmm17, xmm19);
      xmm27 = _mm_unpackhi_epi64(xmm17, xmm19);
      xmm28 = _mm_unpacklo_epi64(xmm20, xmm22);
      xmm29 = _mm_unpackhi_epi64(xmm20, xmm22);
      xmm30 = _mm_unpacklo_epi64(xmm21, xmm23);
      xmm31 = _mm_unpackhi_epi64(xmm21, xmm23);
      
      _mm_store_si128(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+0)) + (CHANNELS*j)), xmm24);
      _mm_store_si128(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+1)) + (CHANNELS*j)), xmm25);
      _mm_store_si128(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+2)) + (CHANNELS*j)), xmm26);
      _mm_store_si128(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+3)) + (CHANNELS*j)), xmm27);
      _mm_store_si128(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+4)) + (CHANNELS*j)), xmm28);
      _mm_store_si128(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+5)) + (CHANNELS*j)), xmm29);
      _mm_store_si128(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+6)) + (CHANNELS*j)), xmm30);
      _mm_store_si128(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+7)) + (CHANNELS*j)), xmm31);
    }
  }
}

}  // anonymous namespace

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
	      int qlen) -> void
{
  
  volatile __m128i S;
  __m128i SL;
  __m128i Q;
  __m128i R;
  __m128i T;
  __m128i M;
  __m128i Z;
  __m128i T0;
  auto * const hep = reinterpret_cast<__m128i*>(hearray);
  __m128i ** const qp = reinterpret_cast<__m128i**>(q_start);
  std::array<BYTE const *, CHANNELS> d_begin;
  std::array<BYTE const *, CHANNELS> d_pos;
  std::array<BYTE const *, CHANNELS> d_best;
  std::array<BYTE const *, CHANNELS> d_end;

  // the database residues of the channels, 16-byte aligned for the loads
  alignas(__m128i) std::array<BYTE, CDEPTH * sizeof(__m128i)> dseqalloc;

  auto * dseq = dseqalloc.data();
  BYTE const zero = 0;

  std::array<long, CHANNELS> seq_id;
  long next_id = 0;
  unsigned done = 0;
  
  Z = _mm_set1_epi16(word_0x8000);
  T0 = _mm_set_epi16(0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, 0x0000, word_0x8000);
  Q  = _mm_set1_epi16(static_cast<short>(gap_open_penalty));
  R  = _mm_set1_epi16(static_cast<short>(gap_extend_penalty));
  
  S = Z;
  SL = Z;
      
  for(long a=0; a < qlen; a++)
  {
    hep[2*a] = Z;
    hep[(2*a)+1] = Z;
  }

  for (std::size_t c = 0; c < CHANNELS; c++)
  {
    d_begin[c] = &zero;
    d_pos[c] = d_begin[c];
    d_best[c] = d_begin[c];
    d_end[c] = d_begin[c];
    seq_id[c] = -1;
  }

  int easy = 0;

  while(true)
  {
    if (easy != 0)
    {
      for (std::size_t c = 0; c < CHANNELS; c++)
      {
	for (std::size_t j = 0; j < CDEPTH; j++)
	{
	  if (d_pos[c] < d_end[c])
	  {
	    dseq[(CHANNELS*j)+c] = *(d_pos[c]++);
	  }
	  else
	  {
	    dseq[(CHANNELS * j) + c] = 0;
	  }
	}
	if ((d_pos[c] == d_end[c]) && (seq_id[c] > -1))
	{
	  easy = 0;
	}
      }
	
      dprofile_fill16(dprofile, score_matrix, dseq);
      	  
      donormal16(&S, hep, qp, &Q, &R, qlen, &Z);

      /* save column address if new highscore */
      
      auto const mask = _mm_movemask_epi8(_mm_cmpgt_epi16(S, SL));
      for (std::size_t c = 0; c < CHANNELS; c++)
      {
	if ((mask & (3 << 2 * c)) != 0)
	{
	  d_best[c] = d_pos[c];
	}
      }

      SL = S;
    }	  
    else
    {

      easy = 1;
 
      M = _mm_setzero_si128();
      T = T0;

      for (std::size_t c = 0; c < CHANNELS; c++)
      {
	if (d_pos[c] < d_end[c])
	{
	  for (std::size_t j = 0; j < CDEPTH; j++)
	  {
	    if (d_pos[c] < d_end[c])
	    {
	      dseq[(CHANNELS*j)+c] = *(d_pos[c]++);
	    }
	    else
	    {
	      dseq[(CHANNELS * j) + c] = 0;
	    }
	  }

	  if (d_pos[c] == d_end[c])
	  {
	    easy = 0;
	  }
	}
	else
	{
	  M = _mm_xor_si128(M, T);
		  
	  long const cand_id = seq_id[c];
		  
	  if (cand_id >= 0)
	  {
	    long const score = reinterpret_cast<WORD *>(const_cast<__m128i *>(&S))[c] ^ 0x8000;
	    scores[cand_id] = score;
	    bestpos[cand_id] = d_best[c] - d_begin[c];
	    done++;
	  }
		  
	  if (next_id < sequences)
	  {
	    seq_id[c] = next_id;
	    long const seqnosf = seqnos[next_id];
	    long ntlen = 0;

	    long const strand = (seqnosf >> 2) & 1;
	    long const frame = seqnosf & 3;
	    long const seqno = seqnosf >> 3;

	    View<char> const sequence =
	      db_getsequence(dbt, seqno, {strand, frame}, & ntlen, c);
		      
	    d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
	    d_pos[c] = d_begin[c];
	    d_best[c] = d_begin[c];
	    d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
	    next_id++;
		      
	    for (std::size_t j = 0; j < CDEPTH; j++)
	    {
	      if (d_pos[c] < d_end[c])
	      {
		dseq[(CHANNELS*j)+c] = *(d_pos[c]++);
	      }
	      else
	      {
		dseq[(CHANNELS * j) + c] = 0;
	      }
	    }
	    if (d_pos[c] == d_end[c])
	    {
	      easy = 0;
	    }
	  }
	  else
	  {
	    seq_id[c] = -1;
	    d_begin[c] = &zero;
	    d_pos[c] = d_begin[c];
	    d_best[c] = d_begin[c];
	    d_end[c] = d_begin[c];
	    for (std::size_t j = 0; j < CDEPTH; j++)
	    {
	      dseq[(CHANNELS * j) + c] = 0;
	    }
	  }
	}
	T = _mm_slli_si128(T, 2);
      }

      if (done == sequences)
      {
	break;
      }

      dprofile_fill16(dprofile, score_matrix, dseq);
      	  
      domasked16(&S, hep, qp, &Q, &R, qlen, &Z, &M);

      /* save column address if new highscore */
      
      SL = _mm_adds_epi16(SL, M);
      SL = _mm_adds_epi16(SL, M);
      auto const mask = _mm_movemask_epi8(_mm_cmpgt_epi16(S, SL));
      for (std::size_t c = 0; c < CHANNELS; c++)
      {
	if ((mask & (3 << 2 * c)) != 0)
	{
	  d_best[c] = d_pos[c];
	}
      }

      SL = S;
    }
  }
}
