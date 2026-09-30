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

// The cell of the alignment kernels, shared by search7.cc (7-bit
// scores), search16.cc and search16s.cc (16-bit scores). As in swarm
// (src/utils/align_cells.hpp), the score width arrives as a policy
// type, Ops, with three static members: add and sub (saturated) and
// max. A class rather than function pointers, so that the calls
// disappear at instantiation.

#ifndef SWIPE_ALIGN_CELLS_H
#define SWIPE_ALIGN_CELLS_H

#include "intrinsics_to_functions.h"
#include <emmintrin.h>  // __m128i


// 16 lanes of 7-bit scores: signed saturated arithmetic, unsigned maxima
struct Ops_7 {
  static auto add(__m128i const lhs, __m128i const rhs) -> __m128i { return v_adds_i8(lhs, rhs); }
  static auto sub(__m128i const lhs, __m128i const rhs) -> __m128i { return v_subs_i8(lhs, rhs); }
  static auto max(__m128i const lhs, __m128i const rhs) -> __m128i { return v_max_u8(lhs, rhs); }
};

// 8 lanes of 16-bit scores: signed saturated arithmetic, signed maxima
struct Ops_16 {
  static auto add(__m128i const lhs, __m128i const rhs) -> __m128i { return v_adds_i16(lhs, rhs); }
  static auto sub(__m128i const lhs, __m128i const rhs) -> __m128i { return v_subs_i16(lhs, rhs); }
  static auto max(__m128i const lhs, __m128i const rhs) -> __m128i { return v_max_i16(lhs, rhs); }
};


// C++26 refactoring: std::simd, with std::add_sat and std::sub_sat

// One cell of a block (the ONESTEP macro of the former inline
// assembly): H is the score of the diagonal cell, N receives the score
// of this cell (the diagonal of the next column), F and E are the
// vertical and horizontal gap scores, S is the running maximum
template <typename Ops>
inline auto onestep(__m128i const H,
                    __m128i & N,
                    __m128i & F,
                    __m128i const V,
                    __m128i & E,
                    __m128i & S,
                    __m128i const Q,
                    __m128i const R) -> void
{
  auto cell = Ops::add(H, V);
  cell = Ops::max(cell, F);
  cell = Ops::max(cell, E);
  S = Ops::max(cell, S);
  F = Ops::sub(F, R);
  E = Ops::sub(E, R);
  N = cell;
  cell = Ops::sub(cell, Q);
  E = Ops::max(cell, E);
  F = Ops::max(cell, F);
}

#endif  // SWIPE_ALIGN_CELLS_H
