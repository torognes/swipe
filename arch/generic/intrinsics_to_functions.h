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

// The vector operations of the 128-bit kernels, without intrinsics:
// the same names and the same results as
// arch/x86_64/intrinsics_to_functions.h (the SSE2 operations), on the
// generic vectors of GCC and clang (vector_size(16)), that each target
// compiles to its own SIMD instructions (SSE2, NEON, AltiVec, RVV), or
// to scalar code without SIMD (GitHub #14). For any machine without a
// backend of its own (riscv64 first), and for x86-64 with "make
// ARCH=generic", so that the whole test-suite checks these functions
// on x86-64 too.
//
// The saturated additions and subtractions are written without
// widening (wrapping sum, then the lanes that overflowed replaced by
// the bound of their sign), and the maxima with comparison masks: a
// first version, with one loop over the lanes per operation, was 48
// times slower than the SSE2 kernels on x86-64 (GCC 15 unrolls these
// loops in the kernel before it can vectorize them).
//
// Lanes are numbered in memory order, as on x86: lane i of 8 bits is
// byte i, lane i of 16 bits is bytes 2i and 2i + 1, in host byte order
// (as the WORD arrays the kernels store to and read from). The
// interleaves (v_merge_*) move groups of bytes, whatever the byte
// order of the host.
//
// Not C++11: the vector_size attribute is a GCC and clang extension
// (GCC 4.8 has it, with the operators used here; this file is only
// compiled for ARCH=generic, with a recent compiler).

#ifndef SWIPE_ARCH_GENERIC_INTRINSICS_TO_FUNCTIONS_H
#define SWIPE_ARCH_GENERIC_INTRINSICS_TO_FUNCTIONS_H

#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // std::int8_t, std::int16_t, std::uint8_t, std::uint16_t
#include <cstring>  // std::memcpy


constexpr std::size_t v128_bytes = 16;

// 16 bytes: 16 lanes of 8 bits, or 8 lanes of 16 bits. The kernels
// cast their byte buffers to v128_t pointers and dereference them, as
// they do with __m128i: may_alias, as __m128i has, tells the compiler
// that a v128_t can be read and written over other types
using v128_t [[gnu::vector_size(v128_bytes), gnu::may_alias]] = std::uint8_t;

static_assert(sizeof(v128_t) == v128_bytes, "16 bytes");


namespace generic_lanes {

// the same 16 bytes, seen as signed bytes, or as words
using i8x16 [[gnu::vector_size(v128_bytes), gnu::may_alias]] = std::int8_t;
using i16x8 [[gnu::vector_size(v128_bytes), gnu::may_alias]] = std::int16_t;
using u16x8 [[gnu::vector_size(v128_bytes), gnu::may_alias]] = std::uint16_t;

constexpr std::int8_t i8_max = INT8_MAX;
constexpr std::int16_t i16_max = INT16_MAX;
constexpr int i8_sign_shift = 7;  // all bits of a lane set to its sign bit
constexpr int i16_sign_shift = 15;

// the lanes of sum (a wrapping sum or difference of lhs and rhs) that
// overflowed, replaced by the bound of the sign of lhs: lanes of type
// Vector, sign_shift = bits per lane - 1. The sum overflowed when its
// sign differs from the sign of lhs and from the sign of rhs (sum), or
// of -rhs (difference): 'other' is rhs for a sum, ~rhs for a
// difference
template <typename Vector, typename Lane>
inline auto saturated(Vector const lhs, Vector const other, Vector const sum,
                      Lane const lane_max, int const sign_shift) -> Vector
{
  Vector const overflow = ((lhs ^ sum) & (other ^ sum)) >> sign_shift;
  Vector const bound = (lhs >> sign_shift) ^ lane_max;
  return (sum & ~overflow) | (bound & overflow);
}

// the larger of each pair of lanes (a comparison is a mask: all bits
// of a lane set where it is true)
template <typename Vector>
inline auto maximum(Vector const lhs, Vector const rhs) -> Vector
{
  auto const larger = reinterpret_cast<Vector>(lhs > rhs);
  return (lhs & larger) | (rhs & ~larger);
}

// the interleave of groups of 'group' bytes of the low (half = 0) or
// high (half = 1) halves of two vectors: lhs group, rhs group, lhs
// group, ... (punpckl* and punpckh*)
template <std::size_t group, std::size_t half>
inline auto merge(v128_t const lhs, v128_t const rhs) -> v128_t
{
  static_assert(v128_bytes % (2 * group) == 0, "groups of a power of two bytes");
  constexpr std::size_t groups = v128_bytes / 2 / group;  // per half
  constexpr std::size_t first = half * (v128_bytes / 2);
  std::array<std::uint8_t, v128_bytes> left {{}};
  std::array<std::uint8_t, v128_bytes> right {{}};
  std::array<std::uint8_t, v128_bytes> result {{}};
  std::memcpy(left.data(), &lhs, v128_bytes);
  std::memcpy(right.data(), &rhs, v128_bytes);
  for (std::size_t index = 0; index < groups; ++index)
  {
    auto const source = first + (index * group);
    auto const target = 2 * index * group;
    std::memcpy(&result[target], &left[source], group);
    std::memcpy(&result[target + group], &right[source], group);
  }
  v128_t vector;
  std::memcpy(&vector, result.data(), v128_bytes);
  return vector;
}

// the lanewise comparison of 16-bit signed lanes (a mask), packed as
// pmovmskb does: two bits per lane, bits 2i and 2i + 1 for lane i
inline auto byte_mask(i16x8 const comparison) -> unsigned int
{
  constexpr unsigned int lane_bits = 3;  // the two bits of a 16-bit lane
  constexpr std::size_t lanes = v128_bytes / sizeof(std::int16_t);
  unsigned int mask = 0;
  for (std::size_t lane = 0; lane < lanes; ++lane)
  {
    if (comparison[lane] != 0)
    {
      mask |= lane_bits << (2 * lane);
    }
  }
  return mask;
}

}  // namespace generic_lanes


// 16 lanes of 8 bits

// signed saturated addition (paddsb)
inline auto v_adds_i8(v128_t const lhs, v128_t const rhs) -> v128_t
{
  using generic_lanes::i8x16;
  auto const sum = reinterpret_cast<i8x16>(lhs + rhs);  // wraps
  return reinterpret_cast<v128_t>(
      generic_lanes::saturated(reinterpret_cast<i8x16>(lhs), reinterpret_cast<i8x16>(rhs),
                               sum, generic_lanes::i8_max, generic_lanes::i8_sign_shift));
}

// signed saturated subtraction (psubsb)
inline auto v_subs_i8(v128_t const lhs, v128_t const rhs) -> v128_t
{
  using generic_lanes::i8x16;
  auto const difference = reinterpret_cast<i8x16>(lhs - rhs);  // wraps
  return reinterpret_cast<v128_t>(
      generic_lanes::saturated(reinterpret_cast<i8x16>(lhs), reinterpret_cast<i8x16>(~rhs),
                               difference, generic_lanes::i8_max, generic_lanes::i8_sign_shift));
}

// unsigned maximum (pmaxub)
inline auto v_max_u8(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return generic_lanes::maximum(lhs, rhs);
}


// 8 lanes of 16 bits

// signed saturated addition (paddsw)
inline auto v_adds_i16(v128_t const lhs, v128_t const rhs) -> v128_t
{
  using generic_lanes::i16x8;
  using generic_lanes::u16x8;
  auto const sum = reinterpret_cast<i16x8>(reinterpret_cast<u16x8>(lhs) + reinterpret_cast<u16x8>(rhs));
  return reinterpret_cast<v128_t>(
      generic_lanes::saturated(reinterpret_cast<i16x8>(lhs), reinterpret_cast<i16x8>(rhs),
                               sum, generic_lanes::i16_max, generic_lanes::i16_sign_shift));
}

// signed saturated subtraction (psubsw)
inline auto v_subs_i16(v128_t const lhs, v128_t const rhs) -> v128_t
{
  using generic_lanes::i16x8;
  using generic_lanes::u16x8;
  auto const difference = reinterpret_cast<i16x8>(reinterpret_cast<u16x8>(lhs) - reinterpret_cast<u16x8>(rhs));
  return reinterpret_cast<v128_t>(
      generic_lanes::saturated(reinterpret_cast<i16x8>(lhs), reinterpret_cast<i16x8>(~rhs),
                               difference, generic_lanes::i16_max, generic_lanes::i16_sign_shift));
}

// signed maximum (pmaxsw)
inline auto v_max_i16(v128_t const lhs, v128_t const rhs) -> v128_t
{
  using generic_lanes::i16x8;
  return reinterpret_cast<v128_t>(
      generic_lanes::maximum(reinterpret_cast<i16x8>(lhs), reinterpret_cast<i16x8>(rhs)));
}


// whole vectors

// aligned load of 16 bytes
inline auto v_load(v128_t const * const ptr) -> v128_t
{
  return *ptr;
}

// aligned store of 16 bytes
inline auto v_store(v128_t * const ptr, v128_t const vector) -> void
{
  *ptr = vector;
}

// load of 8 bytes into the low half, the high half set to zero (movq):
// only 8 bytes are read, ptr needs no alignment
inline auto v_load_64(void const * const ptr) -> v128_t
{
  v128_t vector {};
  std::memcpy(&vector, ptr, v128_bytes / 2);
  return vector;
}


// interleaving of the low (merge_lo) or high (merge_hi) halves of two
// vectors, by lanes of 8, 16, 32 or 64 bits: the transpositions of the
// score profile builders

inline auto v_merge_lo_8(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return generic_lanes::merge<1, 0>(lhs, rhs);
}

inline auto v_merge_lo_16(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return generic_lanes::merge<2, 0>(lhs, rhs);
}

inline auto v_merge_hi_16(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return generic_lanes::merge<2, 1>(lhs, rhs);
}

inline auto v_merge_lo_32(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return generic_lanes::merge<4, 0>(lhs, rhs);
}

inline auto v_merge_hi_32(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return generic_lanes::merge<4, 1>(lhs, rhs);
}

inline auto v_merge_lo_64(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return generic_lanes::merge<8, 0>(lhs, rhs);
}

inline auto v_merge_hi_64(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return generic_lanes::merge<8, 1>(lhs, rhs);
}


// bitwise operations (whole vectors)

inline auto v_and(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return lhs & rhs;
}

inline auto v_or(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return lhs | rhs;
}

inline auto v_xor(v128_t const lhs, v128_t const rhs) -> v128_t
{
  return lhs ^ rhs;
}


// all 16 lanes of 8 bits set to value
inline auto v_dup_i8(char const value) -> v128_t
{
  return v128_t {} + static_cast<std::uint8_t>(value);
}

// each lane of 16 bits shifted left by 'count' bits (psllw)
template <int count>
inline auto v_shift_left_i16(v128_t const vector) -> v128_t
{
  static_assert(count >= 0 and count < 16, "a shift within the lane");
  return reinterpret_cast<v128_t>(reinterpret_cast<generic_lanes::u16x8>(vector) << count);
}


// all 8 lanes of 16 bits set to value
inline auto v_dup_i16(short const value) -> v128_t
{
  return reinterpret_cast<v128_t>(generic_lanes::u16x8 {} + static_cast<std::uint16_t>(value));
}

// lanewise comparisons of 16-bit signed lanes, packed into a byte mask
// (as pcmpgtw or pcmpeqw, then pmovmskb): two bits per lane, bits 2i
// and 2i + 1 for lane i
inline auto v_mask_gt_i16(v128_t const lhs, v128_t const rhs) -> unsigned int
{
  using generic_lanes::i16x8;
  return generic_lanes::byte_mask(reinterpret_cast<i16x8>(lhs) > reinterpret_cast<i16x8>(rhs));
}

inline auto v_mask_eq_i16(v128_t const lhs, v128_t const rhs) -> unsigned int
{
  using generic_lanes::i16x8;
  return generic_lanes::byte_mask(reinterpret_cast<i16x8>(lhs) == reinterpret_cast<i16x8>(rhs));
}

#endif  // SWIPE_ARCH_GENERIC_INTRINSICS_TO_FUNCTIONS_H
