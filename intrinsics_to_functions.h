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

// The SSE2 lane operations of the alignment kernels, as named
// functions. Same file name as in swarm (src/arch/x86_64/), but the
// functions are defined inline here (swarm defines them out of line
// and relies on -flto), and their names carry the signedness and the
// lane width, as swipe mixes signed saturated arithmetic and unsigned
// maxima on the same bytes: another architecture (NEON, AltiVec) can
// provide the same functions (phase 10)

#ifndef SWIPE_INTRINSICS_TO_FUNCTIONS_H
#define SWIPE_INTRINSICS_TO_FUNCTIONS_H

#include <emmintrin.h>  // SSE2


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

#endif  // SWIPE_INTRINSICS_TO_FUNCTIONS_H
