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

// The cells of align_cells() and align_cells_single() (align_cells.h)
// with RVV intrinsics: one instruction for each saturated addition,
// subtraction and maximum, where the generic vectors need five to
// seven (arch/generic/intrinsics_to_functions.h). As in swarm
// (src/arch/riscv64/align_rvv.cpp), the RVV types are sizeless: the
// vectors live in registers here, and the buffers (H and E scores,
// score profile) keep their v128_t layout.
//
// vl is fixed at 128 bits, the width the buffers are laid out for (16
// lanes of 8 bits, or 8 lanes of 16 bits): the V extension guarantees
// VLEN >= 128.

#include "arch/riscv64/align_cells_rvv.h"

#ifdef __riscv_vector

#include "arch/intrinsics_to_functions.h"  // v128_t, v128_bytes
#include <riscv_vector.h>
#include <cstddef>  // std::size_t
#include <cstdint>  // std::int8_t, std::int16_t
#include <iterator>  // std::next


namespace {

// 16 lanes of 7-bit scores: signed saturated arithmetic, unsigned
// maxima (Ops_7)
struct Rvv_7 {
  using Vector = vint8m1_t;
  static constexpr std::size_t lanes = v128_bytes / sizeof(std::int8_t);
  static auto load(v128_t const * const source) -> Vector {
    return __riscv_vle8_v_i8m1(reinterpret_cast<std::int8_t const *>(source), lanes);
  }
  static auto store(v128_t * const target, Vector const vector) -> void {
    __riscv_vse8_v_i8m1(reinterpret_cast<std::int8_t *>(target), vector, lanes);
  }
  static auto add(Vector const lhs, Vector const rhs) -> Vector {
    return __riscv_vsadd_vv_i8m1(lhs, rhs, lanes);
  }
  static auto sub(Vector const lhs, Vector const rhs) -> Vector {
    return __riscv_vssub_vv_i8m1(lhs, rhs, lanes);
  }
  static auto max(Vector const lhs, Vector const rhs) -> Vector {
    return __riscv_vreinterpret_v_u8m1_i8m1(
        __riscv_vmaxu_vv_u8m1(__riscv_vreinterpret_v_i8m1_u8m1(lhs),
                              __riscv_vreinterpret_v_i8m1_u8m1(rhs), lanes));
  }
  // as Ops_7::mask()
  static auto mask(Vector const lhs, Vector const rhs) -> Vector {
    return add(lhs, rhs);
  }
};

// 8 lanes of 16-bit scores: signed saturated arithmetic, signed maxima
// (Ops_16)
struct Rvv_16 {
  using Vector = vint16m1_t;
  static constexpr std::size_t lanes = v128_bytes / sizeof(std::int16_t);
  static auto load(v128_t const * const source) -> Vector {
    return __riscv_vle16_v_i16m1(reinterpret_cast<std::int16_t const *>(source), lanes);
  }
  static auto store(v128_t * const target, Vector const vector) -> void {
    __riscv_vse16_v_i16m1(reinterpret_cast<std::int16_t *>(target), vector, lanes);
  }
  static auto add(Vector const lhs, Vector const rhs) -> Vector {
    return __riscv_vsadd_vv_i16m1(lhs, rhs, lanes);
  }
  static auto sub(Vector const lhs, Vector const rhs) -> Vector {
    return __riscv_vssub_vv_i16m1(lhs, rhs, lanes);
  }
  static auto max(Vector const lhs, Vector const rhs) -> Vector {
    return __riscv_vmax_vv_i16m1(lhs, rhs, lanes);
  }
  // as Ops_16::mask()
  static auto mask(Vector const lhs, Vector const rhs) -> Vector {
    return add(add(lhs, rhs), rhs);
  }
};


// The gap penalties, in registers
template <typename Rvv>
struct Gaps {
  typename Rvv::Vector open;
  typename Rvv::Vector extend;
};


// One cell, as onestep() of align_cells.h
template <typename Rvv>
inline auto onestep(typename Rvv::Vector & S,
                    typename Rvv::Vector const H,
                    typename Rvv::Vector & N,
                    typename Rvv::Vector const V,
                    typename Rvv::Vector & F,
                    typename Rvv::Vector const gap_open,
                    typename Rvv::Vector const gap_extend,
                    typename Rvv::Vector & E) -> void
{
  auto cell = Rvv::add(H, V);
  cell = Rvv::max(cell, F);
  cell = Rvv::max(cell, E);
  S = Rvv::max(cell, S);
  F = Rvv::sub(F, gap_extend);
  E = Rvv::sub(E, gap_extend);
  N = cell;
  cell = Rvv::sub(cell, gap_open);
  E = Rvv::max(cell, E);
  F = Rvv::max(cell, F);
}


// align_cells(): a block of four database residues per channel; the
// masking of the lanes of 'mask' when 'masked' (a compile-time
// constant: the regular kernel contains no masking code)
template <typename Rvv, bool masked>
inline auto align_four(v128_t & S_vector,
                       v128_t * const hep,
                       v128_t const * const * const qp,
                       v128_t const & gap_open_vector,
                       v128_t const & gap_extend_vector,
                       long const ql,
                       v128_t const & Z_vector,
                       v128_t const * const mask) -> void
{
  auto const gap_open = Rvv::load(&gap_open_vector);
  auto const gap_extend = Rvv::load(&gap_extend_vector);
  auto const Z = Rvv::load(&Z_vector);
  auto const M = masked ? Rvv::load(mask) : Z;  // unused when not masked
  auto score = Rvv::load(&S_vector);
  if (masked)
  {
    score = Rvv::mask(score, M);
  }
  auto H0 = Z;
  auto H1 = Z;
  auto H2 = Z;
  auto H3 = Z;
  auto F0 = Z;
  auto F1 = Z;
  auto F2 = Z;
  auto F3 = Z;
  auto N1 = Z;
  auto N2 = Z;
  auto N3 = Z;
  auto N4 = Z;

  for (long qi = 0; qi < ql; ++qi)
  {
    v128_t const * const x = *std::next(qp, qi);
    auto * const h_cell = std::next(hep, 2 * qi);
    auto * const e_cell = std::next(h_cell);
    auto N0 = Rvv::load(h_cell);
    auto E = Rvv::load(e_cell);
    if (masked)
    {
      N0 = Rvv::mask(N0, M);
      E = Rvv::mask(E, M);
    }

    onestep<Rvv>(score, H0, N1, Rvv::load(x), F0, gap_open, gap_extend, E);
    onestep<Rvv>(score, H1, N2, Rvv::load(std::next(x, 1)), F1, gap_open, gap_extend, E);
    onestep<Rvv>(score, H2, N3, Rvv::load(std::next(x, 2)), F2, gap_open, gap_extend, E);
    onestep<Rvv>(score, H3, N4, Rvv::load(std::next(x, 3)), F3, gap_open, gap_extend, E);

    Rvv::store(h_cell, N4);
    Rvv::store(e_cell, E);
    H0 = N0;
    H1 = N1;
    H2 = N2;
    H3 = N3;
  }

  Rvv::store(&S_vector, score);
}


// align_cells_single(): one database residue per channel
template <typename Rvv, bool masked>
inline auto align_single(v128_t & S_vector,
                         v128_t * const hep,
                         v128_t const * const * const qp,
                         v128_t const & gap_open_vector,
                         v128_t const & gap_extend_vector,
                         long const ql,
                         v128_t const & Z_vector,
                         v128_t const * const mask) -> void
{
  auto const gap_open = Rvv::load(&gap_open_vector);
  auto const gap_extend = Rvv::load(&gap_extend_vector);
  auto const Z = Rvv::load(&Z_vector);
  auto const M = masked ? Rvv::load(mask) : Z;  // unused when not masked
  auto score = Rvv::load(&S_vector);
  if (masked)
  {
    score = Rvv::mask(score, M);
  }
  auto H0 = Z;
  auto F0 = Z;
  auto N1 = Z;

  for (long qi = 0; qi < ql; ++qi)
  {
    v128_t const * const x = *std::next(qp, qi);
    auto * const h_cell = std::next(hep, 2 * qi);
    auto * const e_cell = std::next(h_cell);
    auto N0 = Rvv::load(h_cell);
    auto E = Rvv::load(e_cell);
    if (masked)
    {
      N0 = Rvv::mask(N0, M);
      E = Rvv::mask(E, M);
    }

    onestep<Rvv>(score, H0, N1, Rvv::load(x), F0, gap_open, gap_extend, E);

    Rvv::store(h_cell, N1);
    Rvv::store(e_cell, E);
    H0 = N0;
  }

  Rvv::store(&S_vector, score);
}

}  // anonymous namespace


auto align_rvv_7(v128_t & S, v128_t * const hep, v128_t const * const * const qp,
                 v128_t const & gap_open, v128_t const & gap_extend,
                 long const ql, v128_t const & Z, v128_t const * const mask) -> void
{
  if (mask == nullptr)
  {
    align_four<Rvv_7, false>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
    return;
  }
  align_four<Rvv_7, true>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
}


auto align_rvv_16(v128_t & S, v128_t * const hep, v128_t const * const * const qp,
                  v128_t const & gap_open, v128_t const & gap_extend,
                  long const ql, v128_t const & Z, v128_t const * const mask) -> void
{
  if (mask == nullptr)
  {
    align_four<Rvv_16, false>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
    return;
  }
  align_four<Rvv_16, true>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
}


auto align_rvv_16_single(v128_t & S, v128_t * const hep, v128_t const * const * const qp,
                         v128_t const & gap_open, v128_t const & gap_extend,
                         long const ql, v128_t const & Z, v128_t const * const mask) -> void
{
  if (mask == nullptr)
  {
    align_single<Rvv_16, false>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
    return;
  }
  align_single<Rvv_16, true>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
}

#endif  // __riscv_vector
