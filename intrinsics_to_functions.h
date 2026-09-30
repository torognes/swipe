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

// The SSE2 operations of the alignment kernels and of the score
// profile builders, as named functions. Same file name as in swarm (src/arch/x86_64/), but the
// functions are defined inline here (swarm defines them out of line
// and relies on -flto), and their names carry the signedness and the
// lane width, as swipe mixes signed saturated arithmetic and unsigned
// maxima on the same bytes: another architecture (NEON, AltiVec) can
// provide the same functions (phase 10)

#ifndef SWIPE_INTRINSICS_TO_FUNCTIONS_H
#define SWIPE_INTRINSICS_TO_FUNCTIONS_H

#include <emmintrin.h>  // SSE2
#ifdef __SSSE3__
#include <tmmintrin.h>  // SSSE3
#endif


// 16 lanes of 8 bits

// signed saturated addition (paddsb)
inline auto v_adds_i8(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_adds_epi8(lhs, rhs);
}

// signed saturated subtraction (psubsb)
inline auto v_subs_i8(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_subs_epi8(lhs, rhs);
}

// unsigned maximum (pmaxub)
inline auto v_max_u8(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_max_epu8(lhs, rhs);
}


// 8 lanes of 16 bits

// signed saturated addition (paddsw)
inline auto v_adds_i16(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_adds_epi16(lhs, rhs);
}

// signed saturated subtraction (psubsw)
inline auto v_subs_i16(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_subs_epi16(lhs, rhs);
}

// signed maximum (pmaxsw)
inline auto v_max_i16(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_max_epi16(lhs, rhs);
}


// whole vectors

// aligned load of 16 bytes (movdqa)
inline auto v_load(__m128i const * const ptr) -> __m128i
{
  return _mm_load_si128(ptr);
}

// aligned store of 16 bytes (movdqa)
inline auto v_store(__m128i * const ptr, __m128i const vector) -> void
{
  _mm_store_si128(ptr, vector);
}

// load of 8 bytes into the low half, the high half set to zero (movq)
inline auto v_load_64(__m128i const * const ptr) -> __m128i
{
  return _mm_loadl_epi64(ptr);
}


// interleaving of the low (merge_lo) or high (merge_hi) halves of two
// vectors, by lanes of 8, 16, 32 or 64 bits (punpckl*, punpckh*): the
// transpositions of the score profile builders

inline auto v_merge_lo_8(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_unpacklo_epi8(lhs, rhs);
}

inline auto v_merge_lo_16(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_unpacklo_epi16(lhs, rhs);
}

inline auto v_merge_hi_16(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_unpackhi_epi16(lhs, rhs);
}

inline auto v_merge_lo_32(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_unpacklo_epi32(lhs, rhs);
}

inline auto v_merge_hi_32(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_unpackhi_epi32(lhs, rhs);
}

inline auto v_merge_lo_64(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_unpacklo_epi64(lhs, rhs);
}

inline auto v_merge_hi_64(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_unpackhi_epi64(lhs, rhs);
}


// bitwise operations (whole vectors)

inline auto v_and(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_and_si128(lhs, rhs);
}

inline auto v_or(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_or_si128(lhs, rhs);
}

inline auto v_xor(__m128i const lhs, __m128i const rhs) -> __m128i
{
  return _mm_xor_si128(lhs, rhs);
}


// all 16 lanes of 8 bits set to value
inline auto v_dup_i8(char const value) -> __m128i
{
  return _mm_set1_epi8(value);
}

// each lane of 16 bits shifted left by 'count' bits (psllw); the count
// is a template parameter, as the instruction takes an immediate
template <int count>
inline auto v_shift_left_i16(__m128i const vector) -> __m128i
{
  return _mm_slli_epi16(vector, count);
}

#ifdef __SSSE3__
// table lookup (pshufb): lane i of the result is the lane
// (indices[i] & 0x0f) of table, or zero when bit 7 of indices[i] is
// set. Porting note: NEON's vqtbl1q_u8 returns zero for every index
// from 16, and dprofile_shuffle7() uses indices 16 to 31 (bit 7 clear)
// to read the low 4 bits
inline auto v_shuffle_8(__m128i const table, __m128i const indices) -> __m128i
{
  return _mm_shuffle_epi8(table, indices);
}
#endif

#endif  // SWIPE_INTRINSICS_TO_FUNCTIONS_H
