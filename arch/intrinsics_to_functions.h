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

#endif  // SWIPE_ARCH_INTRINSICS_TO_FUNCTIONS_H
