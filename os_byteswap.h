/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2021 Torbjorn Rognes, University of Oslo,
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

// imported from vsearch (src/utils/os_byteswap.hpp and os_byteswap.cpp,
// commit 234522e1); swipe change: the definitions are inline, in this
// header (out of line, the calls of db_getsequence() cost +0.17%
// instructions: callgrind, blastp vs pdbaa)

#ifndef SWIPE_OS_BYTESWAP_H
#define SWIPE_OS_BYTESWAP_H

// Operating System specific commands to swap bytes
// C++23 refactoring: replace with std::byteswap()


#if defined(_MSC_VER) || defined(_WIN32)

#include <cstdint>  // uint16_t, uint32_t, uint64_t
#include <cstdlib>  // _byteswap_ushort, _byteswap_ulong, _byteswap_uint64

inline auto bswap_16(uint16_t const bsx) noexcept -> uint16_t {
  return _byteswap_ushort(bsx);
}

inline auto bswap_32(uint32_t const bsx) noexcept -> uint32_t {
  return _byteswap_ulong(bsx);
}

inline auto bswap_64(uint64_t const bsx) noexcept -> uint64_t {
  return _byteswap_uint64(bsx);
}


#elif defined(__APPLE__)

// Mac OS X / Darwin features
#include <cstdint>  // uint16_t, uint32_t, uint64_t
#include <libkern/OSByteOrder.h>

inline auto bswap_16(uint16_t const bsx) noexcept -> uint16_t {
  return OSSwapInt16(bsx);
}

inline auto bswap_32(uint32_t const bsx) noexcept -> uint32_t {
  return OSSwapInt32(bsx);
}

inline auto bswap_64(uint64_t const bsx) noexcept -> uint64_t {
  return OSSwapInt64(bsx);
}


#elif defined(__FreeBSD__) || defined(__NetBSD__)

// FreeBSD and NetBSD spell the three functions the same way and differ only
// in which header supplies them.
#include <cstdint>  // uint16_t, uint32_t, uint64_t
#if defined(__FreeBSD__)
#include <sys/endian.h>  // bswap16, bswap32, bswap64
#else
#include <sys/types.h>
#include <machine/bswap.h>  // bswap16, bswap32, bswap64
#endif

inline auto bswap_16(uint16_t const bsx) noexcept -> uint16_t {
  return bswap16(bsx);
}

inline auto bswap_32(uint32_t const bsx) noexcept -> uint32_t {
  return bswap32(bsx);
}

inline auto bswap_64(uint64_t const bsx) noexcept -> uint64_t {
  return bswap64(bsx);
}


#else

// Linux and other operating systems. Use the compiler byteswap builtins
// (GCC 4.8+ / Clang) when available — identical to what glibc's <byteswap.h>
// macros expand to — and fall back to a portable shift-based implementation
// for any other compiler, so the build no longer depends on <byteswap.h>
// being present.
#include <cstdint>  // uint16_t, uint32_t, uint64_t

#if ! (defined(__GNUC__) || defined(__clang__))
// the fallback moves whole byte lanes, so every shift distance below is a
// multiple of this
static constexpr auto bits_per_byte = 8U;
#endif

inline auto bswap_16(uint16_t const bsx) noexcept -> uint16_t {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap16(bsx);
#else
  return static_cast<uint16_t>((bsx >> bits_per_byte) | (bsx << bits_per_byte));
#endif
}

inline auto bswap_32(uint32_t const bsx) noexcept -> uint32_t {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap32(bsx);
#else
  return ((bsx & UINT32_C(0x000000FF)) << (3U * bits_per_byte)) |
         ((bsx & UINT32_C(0x0000FF00)) << (1U * bits_per_byte)) |
         ((bsx & UINT32_C(0x00FF0000)) >> (1U * bits_per_byte)) |
         ((bsx & UINT32_C(0xFF000000)) >> (3U * bits_per_byte));
#endif
}

inline auto bswap_64(uint64_t const bsx) noexcept -> uint64_t {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap64(bsx);
#else
  return ((bsx & UINT64_C(0x00000000000000FF)) << (7U * bits_per_byte)) |
         ((bsx & UINT64_C(0x000000000000FF00)) << (5U * bits_per_byte)) |
         ((bsx & UINT64_C(0x0000000000FF0000)) << (3U * bits_per_byte)) |
         ((bsx & UINT64_C(0x00000000FF000000)) << (1U * bits_per_byte)) |
         ((bsx & UINT64_C(0x000000FF00000000)) >> (1U * bits_per_byte)) |
         ((bsx & UINT64_C(0x0000FF0000000000)) >> (3U * bits_per_byte)) |
         ((bsx & UINT64_C(0x00FF000000000000)) >> (5U * bits_per_byte)) |
         ((bsx & UINT64_C(0xFF00000000000000)) >> (7U * bits_per_byte));
#endif
}

#endif

#endif  // SWIPE_OS_BYTESWAP_H
