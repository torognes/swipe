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

// The vector operations of the 128-bit kernels (v128_t, v_*): the x86
// intrinsics on x86-64, the portable lane loops of arch/generic/ on any
// other machine, or on x86-64 when the Makefile builds ARCH=generic
// (it then defines SWIPE_ARCH_GENERIC). The AVX2 kernels include
// arch/x86_64/intrinsics_to_functions.h directly.

#ifndef SWIPE_ARCH_INTRINSICS_TO_FUNCTIONS_H
#define SWIPE_ARCH_INTRINSICS_TO_FUNCTIONS_H

#if defined(__x86_64__) && not defined(SWIPE_ARCH_GENERIC)
#include "arch/x86_64/intrinsics_to_functions.h"
#else
#include "arch/generic/intrinsics_to_functions.h"
#endif

#include <cassert>
#include <cstdint>  // std::uintptr_t


// The vectors at 'bytes', a buffer aligned on a vector (alignas(v128_t),
// or allocated by xmalloc()): a cast through void *, as a direct cast
// from a byte or word pointer warns on the targets that require
// aligned accesses (-Wcast-align: mips64, and riscv64 by default)
inline auto as_v128(void * const bytes) -> v128_t *
{
  assert(reinterpret_cast<std::uintptr_t>(bytes) % alignof(v128_t) == 0);
  return static_cast<v128_t *>(bytes);
}

inline auto as_v128(void const * const bytes) -> v128_t const *
{
  assert(reinterpret_cast<std::uintptr_t>(bytes) % alignof(v128_t) == 0);
  return static_cast<v128_t const *>(bytes);
}

#endif  // SWIPE_ARCH_INTRINSICS_TO_FUNCTIONS_H
