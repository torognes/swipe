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
  using Vector = __m128i;
  // the zero score is 0x80, and scores are at most 0xff: as signed
  // bytes, at most -1, so that one saturated addition of a -128 lane
  // resets the score to zero
  static auto mask(__m128i const lhs, __m128i const rhs) -> __m128i { return v_adds_i8(lhs, rhs); }
  static auto add(__m128i const lhs, __m128i const rhs) -> __m128i { return v_adds_i8(lhs, rhs); }
  static auto sub(__m128i const lhs, __m128i const rhs) -> __m128i { return v_subs_i8(lhs, rhs); }
  static auto max(__m128i const lhs, __m128i const rhs) -> __m128i { return v_max_u8(lhs, rhs); }
};

// 8 lanes of 16-bit scores: signed saturated arithmetic, signed maxima
struct Ops_16 {
  using Vector = __m128i;
  // the zero score is 0x8000 (-32768), and scores span the whole signed
  // range: two saturated additions of a -32768 lane reset any score to
  // zero
  static auto mask(__m128i const lhs, __m128i const rhs) -> __m128i { return v_adds_i16(v_adds_i16(lhs, rhs), rhs); }
  static auto add(__m128i const lhs, __m128i const rhs) -> __m128i { return v_adds_i16(lhs, rhs); }
  static auto sub(__m128i const lhs, __m128i const rhs) -> __m128i { return v_subs_i16(lhs, rhs); }
  static auto max(__m128i const lhs, __m128i const rhs) -> __m128i { return v_max_i16(lhs, rhs); }
};


#ifdef __AVX2__
// 32 lanes of 7-bit scores (AVX2): the operations of Ops_7
struct Ops_7_avx2 {
  using Vector = __m256i;
  static auto mask(__m256i const lhs, __m256i const rhs) -> __m256i { return v256_adds_i8(lhs, rhs); }
  static auto add(__m256i const lhs, __m256i const rhs) -> __m256i { return v256_adds_i8(lhs, rhs); }
  static auto sub(__m256i const lhs, __m256i const rhs) -> __m256i { return v256_subs_i8(lhs, rhs); }
  static auto max(__m256i const lhs, __m256i const rhs) -> __m256i { return v256_max_u8(lhs, rhs); }
};

// 16 lanes of 16-bit scores (AVX2): the operations of Ops_16
struct Ops_16_avx2 {
  using Vector = __m256i;
  static auto mask(__m256i const lhs, __m256i const rhs) -> __m256i { return v256_adds_i16(v256_adds_i16(lhs, rhs), rhs); }
  static auto add(__m256i const lhs, __m256i const rhs) -> __m256i { return v256_adds_i16(lhs, rhs); }
  static auto sub(__m256i const lhs, __m256i const rhs) -> __m256i { return v256_subs_i16(lhs, rhs); }
  static auto max(__m256i const lhs, __m256i const rhs) -> __m256i { return v256_max_i16(lhs, rhs); }
};
#endif


// The masking of a kernel pass, selected by the type of its last
// argument (as swarm's src/utils/mask_vectors.hpp): No_mask when every
// channel continues its database sequence, Mask when some channels
// start a new one (their lanes of 'lanes' are set, the others are
// zero), so that their scores are reset first. The No_mask overload
// returns its argument: the regular kernel contains no masking code.
struct No_mask {};

// (the vector type comes from Ops: a vector type such as __m128i
// cannot be a template argument, its attributes would be dropped)
template <typename Ops>
struct Mask {
  typename Ops::Vector lanes;
};

template <typename Ops>
inline auto apply_mask(typename Ops::Vector const vector, No_mask const & /*masking*/) -> typename Ops::Vector
{
  return vector;
}

template <typename Ops>
inline auto apply_mask(typename Ops::Vector const vector, Mask<Ops> const & masking) -> typename Ops::Vector
{
  return Ops::mask(vector, masking.lanes);
}


// C++26 refactoring: std::simd, with std::add_sat and std::sub_sat

// One cell of a block (the ONESTEP macro of the former inline
// assembly): H is the score of the diagonal cell, N receives the score
// of this cell (the diagonal of the next column), F and E are the
// vertical and horizontal gap scores, S is the running maximum
template <typename Ops>
inline auto onestep(typename Ops::Vector const H,
                    typename Ops::Vector & N,
                    typename Ops::Vector & F,
                    typename Ops::Vector const V,
                    typename Ops::Vector & E,
                    typename Ops::Vector & S,
                    typename Ops::Vector const Q,
                    typename Ops::Vector const R) -> void
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


// One pass over the query for a block of four database residues per
// channel (the former donormal and domasked kernels of search7.cc and
// search16.cc, selected by Masking)
template <typename Ops, typename Masking>
inline auto align_cells(typename Ops::Vector & S,
                        typename Ops::Vector * hep,
                        typename Ops::Vector * const * qp,
                        typename Ops::Vector const Q,
                        typename Ops::Vector const R,
                        long ql,
                        typename Ops::Vector const Z,
                        Masking const & masking) -> void
{
  using Vector = typename Ops::Vector;
  auto score = apply_mask<Ops>(S, masking);  // mask
  auto H0 = Z;
  auto H1 = H0;
  auto H2 = H0;
  auto H3 = H0;
  auto F0 = H0;
  auto F1 = H0;
  auto F2 = H0;
  auto F3 = H0;
  Vector N1;
  Vector N2;
  Vector N3;

  for (long qi = 0; qi < ql; ++qi)
  {
    Vector const * const x = qp[qi];  // load x from qp[qi]
    auto const N0 = apply_mask<Ops>(hep[2 * qi], masking);  // load N0, mask
    auto E = apply_mask<Ops>(hep[(2 * qi) + 1], masking);  // load E, mask

    onestep<Ops>(H0, N1, F0, x[0], E, score, Q, R);
    onestep<Ops>(H1, N2, F1, x[1], E, score, Q, R);
    onestep<Ops>(H2, N3, F2, x[2], E, score, Q, R);
    onestep<Ops>(H3, hep[2 * qi], F3, x[3], E, score, Q, R);

    hep[(2 * qi) + 1] = E;  // save E
    H0 = N0;
    H1 = N1;
    H2 = N2;
    H3 = N3;
  }

  S = score;  // save S
}


// One pass over the query for a block of database residues (the
// former donormal and domasked kernels, selected by Masking)
// [moved from search16s.cc (align_cells16s), generic over Ops: one
// database residue per channel, for search16s() and its AVX2 version]
template <typename Ops, typename Masking>
inline auto align_cells_single(typename Ops::Vector & S,
                               typename Ops::Vector * hep,
                               typename Ops::Vector * const * qp,
                               typename Ops::Vector const Q,
                               typename Ops::Vector const R,
                               long ql,
                               typename Ops::Vector const Z,
                               Masking const & masking) -> void
{
  using Vector = typename Ops::Vector;
  auto score = apply_mask<Ops>(S, masking);  // mask
  auto H0 = Z;
  auto F0 = H0;

  for (long qi = 0; qi < ql; ++qi)
  {
    Vector const * const x = qp[qi];  // load x from qp[qi]
    auto const N0 = apply_mask<Ops>(hep[2 * qi], masking);  // load N0, mask
    auto E = apply_mask<Ops>(hep[(2 * qi) + 1], masking);  // load E, mask

    onestep<Ops>(H0, hep[2 * qi], F0, x[0], E, score, Q, R);

    hep[(2 * qi) + 1] = E;  // save E
    H0 = N0;
  }

  S = score;  // save S
}

#endif  // SWIPE_ALIGN_CELLS_H
