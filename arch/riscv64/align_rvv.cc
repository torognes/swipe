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
// vl is fixed at the width the buffers are laid out for: 128 bits (16
// lanes of 8 bits, or 8 lanes of 16 bits), which the V extension
// guarantees (VLEN >= 128), or 256 bits (the align_rvv256_* functions,
// for the 32-byte kernels, only called when VLEN >= 256: one register
// at LMUL = 1 in both cases).

#include "arch/riscv64/align_cells_rvv.h"

#ifdef __riscv_vector

#include "arch/intrinsics_to_functions.h"  // v128_t, v128_bytes
#include <riscv_vector.h>
#include <cstddef>  // std::size_t
#include <array>
#include <cassert>
#include <cstdint>  // std::int8_t, std::int16_t, std::uint8_t, std::uintptr_t
#include <type_traits>  // std::is_same
#include <iterator>  // std::next


namespace {

// The vectors of a kernel are stored in buffers of 16 bytes (v128_t),
// or 32 bytes (v256_bytes_t, when VLEN >= 256), seen here as byte
// arrays (Storage): a vector type such as v128_t cannot be a template
// argument, its attributes would be dropped

// 7-bit scores, one byte lane per channel: signed saturated
// arithmetic, unsigned maxima (Ops_7)
template <std::size_t bytes>
struct Rvv_7 {
  using Vector = vint8m1_t;
  using Storage = std::array<std::uint8_t, bytes>;
  static constexpr std::size_t lanes = bytes / sizeof(std::int8_t);
  static auto load(Storage const * const source) -> Vector {
    return __riscv_vle8_v_i8m1(reinterpret_cast<std::int8_t const *>(source), lanes);
  }
  static auto store(Storage * const target, Vector const vector) -> void {
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

// the words of a buffer (aligned on a word: a cast through void *, as
// a direct cast from a byte array warns, -Wcast-align)
template <typename Storage>
inline auto as_words(Storage * const buffer) -> std::int16_t *
{
  assert(reinterpret_cast<std::uintptr_t>(buffer) % alignof(std::int16_t) == 0);
  return static_cast<std::int16_t *>(static_cast<void *>(buffer));
}

template <typename Storage>
inline auto as_words(Storage const * const buffer) -> std::int16_t const *
{
  assert(reinterpret_cast<std::uintptr_t>(buffer) % alignof(std::int16_t) == 0);
  return static_cast<std::int16_t const *>(static_cast<void const *>(buffer));
}

// 16-bit scores, one word lane per channel: signed saturated
// arithmetic, signed maxima (Ops_16)
template <std::size_t bytes>
struct Rvv_16 {
  using Vector = vint16m1_t;
  using Storage = std::array<std::uint8_t, bytes>;
  static constexpr std::size_t lanes = bytes / sizeof(std::int16_t);
  static auto load(Storage const * const source) -> Vector {
    return __riscv_vle16_v_i16m1(as_words(source), lanes);
  }
  static auto store(Storage * const target, Vector const vector) -> void {
    __riscv_vse16_v_i16m1(as_words(target), vector, lanes);
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

using Rvv_7_128 = Rvv_7<v128_bytes>;
using Rvv_16_128 = Rvv_16<v128_bytes>;
using Rvv_7_256 = Rvv_7<v256_bytes>;
using Rvv_16_256 = Rvv_16<v256_bytes>;
static_assert(std::is_same<Rvv_7_256::Storage, v256_bytes_t>::value, "32-byte buffers");

// a 16-byte buffer, as the 128-bit templates see it
using Bytes_128 = Rvv_7_128::Storage;
static_assert(sizeof(Bytes_128) == sizeof(v128_t), "16-byte buffers");

inline auto bytes_128(v128_t & vector) -> Bytes_128 &
{
  return *reinterpret_cast<Bytes_128 *>(&vector);
}

inline auto bytes_128(v128_t const & vector) -> Bytes_128 const &
{
  return *reinterpret_cast<Bytes_128 const *>(&vector);
}

inline auto bytes_128(v128_t * const vectors) -> Bytes_128 *
{
  return reinterpret_cast<Bytes_128 *>(vectors);
}

inline auto bytes_128(v128_t const * const vectors) -> Bytes_128 const *
{
  return reinterpret_cast<Bytes_128 const *>(vectors);
}

inline auto bytes_128(v128_t const * const * const vectors) -> Bytes_128 const * const *
{
  return reinterpret_cast<Bytes_128 const * const *>(vectors);
}


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
inline auto align_four(typename Rvv::Storage & S_vector,
                       typename Rvv::Storage * const hep,
                       typename Rvv::Storage const * const * const qp,
                       typename Rvv::Storage const & gap_open_vector,
                       typename Rvv::Storage const & gap_extend_vector,
                       long const ql,
                       typename Rvv::Storage const & Z_vector,
                       typename Rvv::Storage const * const mask) -> void
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
    auto const * const x = *std::next(qp, qi);
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
inline auto align_single(typename Rvv::Storage & S_vector,
                         typename Rvv::Storage * const hep,
                         typename Rvv::Storage const * const * const qp,
                         typename Rvv::Storage const & gap_open_vector,
                         typename Rvv::Storage const & gap_extend_vector,
                         long const ql,
                         typename Rvv::Storage const & Z_vector,
                         typename Rvv::Storage const * const mask) -> void
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
    auto const * const x = *std::next(qp, qi);
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
    align_four<Rvv_7_128, false>(bytes_128(S), bytes_128(hep), bytes_128(qp), bytes_128(gap_open),
        bytes_128(gap_extend), ql, bytes_128(Z), bytes_128(mask));
    return;
  }
  align_four<Rvv_7_128, true>(bytes_128(S), bytes_128(hep), bytes_128(qp), bytes_128(gap_open),
        bytes_128(gap_extend), ql, bytes_128(Z), bytes_128(mask));
}


auto align_rvv_16(v128_t & S, v128_t * const hep, v128_t const * const * const qp,
                  v128_t const & gap_open, v128_t const & gap_extend,
                  long const ql, v128_t const & Z, v128_t const * const mask) -> void
{
  if (mask == nullptr)
  {
    align_four<Rvv_16_128, false>(bytes_128(S), bytes_128(hep), bytes_128(qp), bytes_128(gap_open),
        bytes_128(gap_extend), ql, bytes_128(Z), bytes_128(mask));
    return;
  }
  align_four<Rvv_16_128, true>(bytes_128(S), bytes_128(hep), bytes_128(qp), bytes_128(gap_open),
        bytes_128(gap_extend), ql, bytes_128(Z), bytes_128(mask));
}


auto align_rvv_16_single(v128_t & S, v128_t * const hep, v128_t const * const * const qp,
                         v128_t const & gap_open, v128_t const & gap_extend,
                         long const ql, v128_t const & Z, v128_t const * const mask) -> void
{
  if (mask == nullptr)
  {
    align_single<Rvv_16_128, false>(bytes_128(S), bytes_128(hep), bytes_128(qp), bytes_128(gap_open),
        bytes_128(gap_extend), ql, bytes_128(Z), bytes_128(mask));
    return;
  }
  align_single<Rvv_16_128, true>(bytes_128(S), bytes_128(hep), bytes_128(qp), bytes_128(gap_open),
        bytes_128(gap_extend), ql, bytes_128(Z), bytes_128(mask));
}



auto align_rvv256_7(v256_bytes_t & S, v256_bytes_t * const hep,
                    v256_bytes_t const * const * const qp,
                    v256_bytes_t const & gap_open, v256_bytes_t const & gap_extend,
                    long const ql, v256_bytes_t const & Z,
                    v256_bytes_t const * const mask) -> void
{
  if (mask == nullptr)
  {
    align_four<Rvv_7_256, false>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
    return;
  }
  align_four<Rvv_7_256, true>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
}


auto align_rvv256_16(v256_bytes_t & S, v256_bytes_t * const hep,
                     v256_bytes_t const * const * const qp,
                     v256_bytes_t const & gap_open, v256_bytes_t const & gap_extend,
                     long const ql, v256_bytes_t const & Z,
                     v256_bytes_t const * const mask) -> void
{
  if (mask == nullptr)
  {
    align_four<Rvv_16_256, false>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
    return;
  }
  align_four<Rvv_16_256, true>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
}


auto align_rvv256_16_single(v256_bytes_t & S, v256_bytes_t * const hep,
                            v256_bytes_t const * const * const qp,
                            v256_bytes_t const & gap_open, v256_bytes_t const & gap_extend,
                            long const ql, v256_bytes_t const & Z,
                            v256_bytes_t const * const mask) -> void
{
  if (mask == nullptr)
  {
    align_single<Rvv_16_256, false>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
    return;
  }
  align_single<Rvv_16_256, true>(S, hep, qp, gap_open, gap_extend, ql, Z, mask);
}

#endif  // __riscv_vector
