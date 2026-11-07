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

// The SIMD instruction sets of an x86-64 processor, detected once
// (cpuid): read by arch/x86_64/cpu_features.cc and by the kernel
// dispatch, never outside arch/x86_64/

#ifndef SWIPE_ARCH_X86_64_CPU_FEATURES_H
#define SWIPE_ARCH_X86_64_CPU_FEATURES_H


// the SIMD instruction sets of the processor, detected once (cpuid)
struct CpuFeatures
{
  bool sse2;
  bool ssse3;
  bool avx2;
};
extern CpuFeatures const cpu_features;

#endif  // SWIPE_ARCH_X86_64_CPU_FEATURES_H
