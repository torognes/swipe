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

// The processor checks of swipe, common interface (as vsearch's
// src/arch/cpu_features.hpp): one implementation per architecture,
// arch/x86_64/cpu_features.cc (cpuid) and arch/generic/cpu_features.cc
// (any other machine: the kernels need no run-time detection), chosen
// by the Makefile (ARCH)

#ifndef SWIPE_ARCH_CPU_FEATURES_H
#define SWIPE_ARCH_CPU_FEATURES_H

#include <cstdio>  // std::FILE


// fatal when the processor cannot run swipe (x86-64: no SSE2)
auto cpu_features_check() noexcept -> void;

// the notice of the plain output about the processor, if any (x86-64:
// the performance is reduced without SSSE3)
auto cpu_features_show(std::FILE * output_handle) -> void;

#endif  // SWIPE_ARCH_CPU_FEATURES_H
