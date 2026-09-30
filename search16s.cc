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
#include "intrinsics_to_functions.h"  // v_load, v_store, v_merge_*, v_dup_*, ...
#include "align_cells.h"  // Ops_16, onestep(), No_mask, Mask

constexpr std::size_t CHANNELS = channels_16;
constexpr std::size_t CDEPTH = 1;

// the word 0x8000 (the lanes of _mm_set_epi16() are short: 0x8000
// does not fit in a signed short, -32768 has the same bits)
constexpr short word_0x8000 = static_cast<short>(-32768);

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

// One pass over the query for a block of database residues (the
// former donormal and domasked kernels, selected by Masking)
template <typename Masking>
inline auto align_cells16s(__m128i & S,
                           __m128i * hep,
                           __m128i * const * qp,
                           __m128i const Q,
                           __m128i const R,
                           long ql,
                           __m128i const Z,
                           Masking const & masking) -> void
{
  auto score = apply_mask<Ops_16>(S, masking);  // mask
  auto H0 = Z;
  auto F0 = H0;

  for (long qi = 0; qi < ql; ++qi)
  {
    __m128i const * const x = qp[qi];  // load x from qp[qi]
    auto const N0 = apply_mask<Ops_16>(hep[2 * qi], masking);  // load N0, mask
    auto E = apply_mask<Ops_16>(hep[(2 * qi) + 1], masking);  // load E, mask

    onestep<Ops_16>(H0, hep[2 * qi], F0, x[0], E, score, Q, R);

    hep[(2 * qi) + 1] = E;  // save E
    H0 = N0;
  }

  S = score;  // save S
}

inline auto dprofile_fill16s(WORD * dprofile_word,
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
      xmm0  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[0] + i));
      xmm1  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[1] + i));
      xmm2  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[2] + i));
      xmm3  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[3] + i));
      xmm4  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[4] + i));
      xmm5  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[5] + i));
      xmm6  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[6] + i));
      xmm7  = v_load(reinterpret_cast<__m128i*>(score_matrix_word + d[7] + i));
      
      xmm8  = v_merge_lo_16(xmm0,  xmm1);
      xmm9  = v_merge_hi_16(xmm0,  xmm1);
      xmm10 = v_merge_lo_16(xmm2,  xmm3);
      xmm11 = v_merge_hi_16(xmm2,  xmm3);
      xmm12 = v_merge_lo_16(xmm4,  xmm5);
      xmm13 = v_merge_hi_16(xmm4,  xmm5);
      xmm14 = v_merge_lo_16(xmm6,  xmm7);
      xmm15 = v_merge_hi_16(xmm6,  xmm7);
      
      xmm16 = v_merge_lo_32(xmm8,  xmm10);
      xmm17 = v_merge_hi_32(xmm8,  xmm10);
      xmm18 = v_merge_lo_32(xmm12, xmm14);
      xmm19 = v_merge_hi_32(xmm12, xmm14);
      xmm20 = v_merge_lo_32(xmm9,  xmm11);
      xmm21 = v_merge_hi_32(xmm9,  xmm11);
      xmm22 = v_merge_lo_32(xmm13, xmm15);
      xmm23 = v_merge_hi_32(xmm13, xmm15);
      
      xmm24 = v_merge_lo_64(xmm16, xmm18);
      xmm25 = v_merge_hi_64(xmm16, xmm18);
      xmm26 = v_merge_lo_64(xmm17, xmm19);
      xmm27 = v_merge_hi_64(xmm17, xmm19);
      xmm28 = v_merge_lo_64(xmm20, xmm22);
      xmm29 = v_merge_hi_64(xmm20, xmm22);
      xmm30 = v_merge_lo_64(xmm21, xmm23);
      xmm31 = v_merge_hi_64(xmm21, xmm23);
      
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+0)) + (CHANNELS*j)), xmm24);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+1)) + (CHANNELS*j)), xmm25);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+2)) + (CHANNELS*j)), xmm26);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+3)) + (CHANNELS*j)), xmm27);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+4)) + (CHANNELS*j)), xmm28);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+5)) + (CHANNELS*j)), xmm29);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+6)) + (CHANNELS*j)), xmm30);
      v_store(reinterpret_cast<__m128i*>(dprofile_word + (CDEPTH*CHANNELS*(i+7)) + (CHANNELS*j)), xmm31);
    }
  }
}

}  // anonymous namespace

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
	       int qlen) -> void
{
  __m128i S;
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
  std::array<BYTE const *, CHANNELS> d_end;
  std::array<BYTE const *, CHANNELS> d_best;
  std::array<long, CHANNELS> q_best;

  // the database residues of the channels, 16-byte aligned for the loads
  alignas(__m128i) std::array<BYTE, CDEPTH * sizeof(__m128i)> dseqalloc;

  auto * dseq = dseqalloc.data();
  BYTE const zero = 0;

  std::array<long, CHANNELS> seq_id;
  long next_id = 0;
  unsigned done = 0;
  
  Z = v_dup_i16(word_0x8000);
  T0 = v_first_lane_i16(word_0x8000);
  Q  = v_dup_i16(static_cast<short>(gap_open_penalty));
  R  = v_dup_i16(static_cast<short>(gap_extend_penalty));
  
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
    d_end[c] = d_begin[c];
    d_best[c] = d_begin[c];
    q_best[c] = -1;
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
	
      dprofile_fill16s(dprofile, score_matrix, dseq);
      	  
      align_cells16s(S, hep, qp, Q, R, qlen, Z, No_mask{});
      
      /* save column address if new highscore */

      auto const mask = v_mask_gt_i16(S, SL);
      if (mask != 0)
      {
	for (std::size_t c = 0; c < CHANNELS; c++)
	{
	  if ((mask & (3 << 2 * c)) != 0)
	  {
	    d_best[c] = d_pos[c] - 1;
	  }
	}

	for(long i = qlen-1; i >= 0; i--)
	{
	  int const m2 = mask & v_mask_eq_i16(hep[2*i], S);
	  if (m2 != 0)
	  {
	    for (std::size_t c = 0; c < CHANNELS; c++)
	    {
	      if ((m2 & (3 << 2 * c)) != 0)
	      {
		q_best[c] = i;
	      }
	    }
	  }
	}
      }

      SL = S;
    }	  
    else
    {

      easy = 1;
 
      M = v_zero();
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
	  M = v_xor(M, T);
		  
	  long const cand_id = seq_id[c];
		  
	  if (cand_id >= 0)
	  {
	    long const score = reinterpret_cast<WORD *>(&S)[c] ^ 0x8000;
	    scores[cand_id] = score;
	    bestpos[cand_id] = d_best[c] - d_begin[c];	    
	    bestq[cand_id] = q_best[c];
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

	    db_mapsequences(dbta[static_cast<std::size_t>(c)], seqno, seqno);

	    View<char> const sequence =
	      db_getsequence(dbta[static_cast<std::size_t>(c)], seqno, {strand, frame}, & ntlen, c);
		      
	    d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
	    d_pos[c] = d_begin[c];
	    d_best[c] = d_begin[c];
	    d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
	    q_best[c] = -1;
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
	    d_pos[c] = &zero;
	    d_end[c] = d_pos[c];
	    for (std::size_t j = 0; j < CDEPTH; j++)
	    {
	      dseq[(CHANNELS * j) + c] = 0;
	    }
	  }
	}
	T = v_shift_bytes_left<2>(T);
      }

      if (done == sequences)
      {
	break;
      }

      dprofile_fill16s(dprofile, score_matrix, dseq);
      	  
      align_cells16s(S, hep, qp, Q, R, qlen, Z, Mask{M});

      /* save column address if new highscore */

      SL = v_adds_i16(SL, M);
      SL = v_adds_i16(SL, M);
      auto const mask = v_mask_gt_i16(S, SL);
      if (mask != 0)
      {
	for (std::size_t c = 0; c < CHANNELS; c++)
	{
	  if ((mask & (3 << 2 * c)) != 0)
	  {
	    d_best[c] = d_pos[c] - 1;
	  }
	}

	for(long i = qlen-1; i >= 0; i--)
	{
	  int const m2 = mask & v_mask_eq_i16(hep[2*i], S);
	  if (m2 != 0)
	  {
	    for (std::size_t c = 0; c < CHANNELS; c++)
	    {
	      if ((m2 & (3 << 2 * c)) != 0)
	      {
		q_best[c] = i;
	      }
	    }
	  }
	}
      }

      SL = S;
    }
  }
}
