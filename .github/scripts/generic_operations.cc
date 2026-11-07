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

// Prints the results of every vector operation of the 128-bit kernels
// (arch/intrinsics_to_functions.h) on many inputs: every pair of
// 8-bit lanes, edge and pseudo-random 16-bit lanes. Built twice, once
// with the x86-64 intrinsics and once with the portable operations
// (-DSWIPE_ARCH_GENERIC), the two outputs must be identical: see
// generic_operations.sh.

#include "arch/intrinsics_to_functions.h"
#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // std::uint8_t, std::uint16_t, std::uint32_t
#include <cstdio>  // std::printf
#include <cstring>  // std::memcpy


namespace {

using Bytes = std::array<std::uint8_t, sizeof(v128_t)>;

auto to_vector(Bytes const & bytes) -> v128_t
{
  v128_t vector;
  std::memcpy(&vector, bytes.data(), sizeof(v128_t));
  return vector;
}

auto print(char const * name, v128_t const vector) -> void
{
  Bytes bytes {{}};
  std::memcpy(bytes.data(), &vector, sizeof(v128_t));
  std::printf("%s", name);
  for (auto const byte : bytes)
  {
    std::printf(" %02x", static_cast<unsigned int>(byte));
  }
  std::printf("\n");
}

// a xorshift generator: the same numbers everywhere
auto next_random(std::uint32_t & state) -> std::uint32_t
{
  state ^= state << 13U;
  state ^= state >> 17U;
  state ^= state << 5U;
  return state;
}

auto check_pair(v128_t const lhs, v128_t const rhs) -> void
{
  print("adds_i8", v_adds_i8(lhs, rhs));
  print("subs_i8", v_subs_i8(lhs, rhs));
  print("max_u8", v_max_u8(lhs, rhs));
  print("adds_i16", v_adds_i16(lhs, rhs));
  print("subs_i16", v_subs_i16(lhs, rhs));
  print("max_i16", v_max_i16(lhs, rhs));
  print("merge_lo_8", v_merge_lo_8(lhs, rhs));
  print("merge_lo_16", v_merge_lo_16(lhs, rhs));
  print("merge_hi_16", v_merge_hi_16(lhs, rhs));
  print("merge_lo_32", v_merge_lo_32(lhs, rhs));
  print("merge_hi_32", v_merge_hi_32(lhs, rhs));
  print("merge_lo_64", v_merge_lo_64(lhs, rhs));
  print("merge_hi_64", v_merge_hi_64(lhs, rhs));
  print("and", v_and(lhs, rhs));
  print("or", v_or(lhs, rhs));
  print("xor", v_xor(lhs, rhs));
  print("shift_left_i16<3>", v_shift_left_i16<3>(lhs));
  std::printf("mask_gt_i16 %04x\n", v_mask_gt_i16(lhs, rhs));
  std::printf("mask_eq_i16 %04x\n", v_mask_eq_i16(lhs, rhs));
  std::printf("mask_eq_i16 (self) %04x\n", v_mask_eq_i16(lhs, lhs));
  // (a byte array: a vector type cannot be a template argument, its
  // attributes would be dropped)
  alignas(v128_t) Bytes memory {{}};
  auto * const vector = reinterpret_cast<v128_t *>(memory.data());
  v_store(vector, lhs);
  print("load", v_load(vector));
  print("load_64", v_load_64(vector));
  v_store(vector, rhs);
  print("store", to_vector(memory));
}

}  // anonymous namespace


auto main() -> int
{
  constexpr std::size_t byte_values = 256;
  constexpr std::size_t lanes = sizeof(v128_t);

  // every pair of 8-bit values: 16 values of lhs against one value of
  // rhs per pair of vectors
  for (std::size_t first = 0; first < byte_values; first += lanes)
  {
    Bytes left {{}};
    for (std::size_t lane = 0; lane < lanes; ++lane)
    {
      left[lane] = static_cast<std::uint8_t>(first + lane);
    }
    for (std::size_t value = 0; value < byte_values; ++value)
    {
      Bytes right {{}};
      right.fill(static_cast<std::uint8_t>(value));
      check_pair(to_vector(left), to_vector(right));
      check_pair(to_vector(right), to_vector(left));
    }
  }

  // 16-bit edge values, then pseudo-random vectors
  constexpr std::array<std::uint16_t, 8> edges {{0x0000, 0x0001, 0x7ffe, 0x7fff,
                                                 0x8000, 0x8001, 0xfffe, 0xffff}};
  for (auto const edge : edges)
  {
    std::array<std::uint16_t, lanes / 2> words {{}};
    std::array<std::uint16_t, lanes / 2> others {{}};
    words.fill(edge);
    others = edges;
    v128_t lhs;
    v128_t rhs;
    std::memcpy(&lhs, words.data(), sizeof(v128_t));
    std::memcpy(&rhs, others.data(), sizeof(v128_t));
    check_pair(lhs, rhs);
    check_pair(rhs, lhs);
  }
  std::uint32_t state = 2463534242U;
  constexpr int random_pairs = 100000;
  for (int pair = 0; pair < random_pairs; ++pair)
  {
    Bytes left {{}};
    Bytes right {{}};
    for (std::size_t lane = 0; lane < lanes; ++lane)
    {
      left[lane] = static_cast<std::uint8_t>(next_random(state));
      right[lane] = static_cast<std::uint8_t>(next_random(state));
    }
    check_pair(to_vector(left), to_vector(right));
  }

  // duplicates
  for (std::size_t value = 0; value < byte_values; ++value)
  {
    print("dup_i8", v_dup_i8(static_cast<char>(value)));
    print("dup_i16", v_dup_i16(static_cast<short>(value * 257)));
  }
  return 0;
}
