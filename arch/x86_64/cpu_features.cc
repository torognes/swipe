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

#include "arch/cpu_features.h"
#include "arch/x86_64/cpu_features.h"  // CpuFeatures, cpu_features
#include "print_view.h"  // fprint
#include "swipe.h"  // fatal
#include <cpuid.h>  // __get_cpuid, bit_SSE2, bit_SSSE3
#include <cstdio>  // std::FILE


namespace {

auto detect_cpu_features() noexcept -> CpuFeatures
{
  // the feature bits of cpuid leaf 1 (bit_SSE2 in edx, bit_SSSE3 in
  // ecx); __get_cpuid() returns 0, and no feature is reported, when the
  // processor does not support that leaf
  unsigned int eax = 0;
  unsigned int ebx = 0;
  unsigned int ecx = 0;
  unsigned int edx = 0;
  if (__get_cpuid(1, &eax, &ebx, &ecx, &edx) == 0)
  {
    return {false, false, false};
  }
  bool const sse2 = (edx & bit_SSE2) != 0;
#ifdef SWIPE_NO_SSSE3
  bool const ssse3 = false;
#else
  bool const ssse3 = (ecx & bit_SSSE3) != 0;
#endif
  // AVX2: the processor (cpuid leaf 7) and the operating system (the
  // 256-bit registers saved on context switches, XCR0): GCC and clang
  // check both (sketch: check that GCC 4.8's libgcc does the XCR0 part)
  __builtin_cpu_init();
#ifdef SWIPE_NO_AVX2
  bool const avx2 = false;
#else
  bool const avx2 = __builtin_cpu_supports("avx2") != 0;
#endif
  return {sse2, ssse3, avx2};
}

}  // anonymous namespace


// detected before main() (dynamic initialization: cpuid has no dependency)
extern CpuFeatures const cpu_features = detect_cpu_features();


auto cpu_features_check() noexcept -> void
{
  if (not cpu_features.sse2)
  {
    fatal("This program requires a processor with SSE2.");
  }
}


auto cpu_features_show(std::FILE * const output_handle) -> void
{
  if (not cpu_features.ssse3)
  {
    fprint(output_handle, "The performance is reduced because this CPU lacks SSSE3.\n\n");
  }
}
