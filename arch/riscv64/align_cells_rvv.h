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

// The native RVV versions of align_cells() and align_cells_single()
// (align_cells.h), for riscv64 with the V extension: the explicit
// specializations below replace the generic-vector kernels with the
// functions of arch/riscv64/align_rvv.cc. Same arguments, same
// results. Included at the end of align_cells.h.

#ifndef SWIPE_ARCH_RISCV64_ALIGN_CELLS_RVV_H
#define SWIPE_ARCH_RISCV64_ALIGN_CELLS_RVV_H

#ifdef __riscv_vector

#include "align_cells.h"  // align_cells(), align_cells_single(), Ops_7, Ops_16, ...
#include "arch/intrinsics_to_functions.h"  // v128_t
#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // std::uint8_t


// The masking of a pass: no_mask (nullptr) for the regular kernels
// [the vectors are passed by reference: a generic vector argument
// would be passed through memory anyway]
auto align_rvv_7(v128_t & S, v128_t * hep, v128_t const * const * qp,
                 v128_t const & gap_open, v128_t const & gap_extend,
                 long ql, v128_t const & Z, v128_t const * mask) -> void;

auto align_rvv_16(v128_t & S, v128_t * hep, v128_t const * const * qp,
                  v128_t const & gap_open, v128_t const & gap_extend,
                  long ql, v128_t const & Z, v128_t const * mask) -> void;

auto align_rvv_16_single(v128_t & S, v128_t * hep, v128_t const * const * qp,
                         v128_t const & gap_open, v128_t const & gap_extend,
                         long ql, v128_t const & Z, v128_t const * mask) -> void;


// The cells of the 32-byte kernels (arch/riscv64/search7_rvv.cc,
// arch/riscv64/search16_rvv256.cc), when VLEN >= 256: 32 lanes of 8
// bits, or 16 lanes of 16 bits, stored in 32-byte buffers
constexpr std::size_t v256_bytes = 32;
using v256_bytes_t = std::array<std::uint8_t, v256_bytes>;

auto align_rvv256_7(v256_bytes_t & S, v256_bytes_t * hep, v256_bytes_t const * const * qp,
                    v256_bytes_t const & gap_open, v256_bytes_t const & gap_extend,
                    long ql, v256_bytes_t const & Z, v256_bytes_t const * mask) -> void;

auto align_rvv256_16(v256_bytes_t & S, v256_bytes_t * hep, v256_bytes_t const * const * qp,
                     v256_bytes_t const & gap_open, v256_bytes_t const & gap_extend,
                     long ql, v256_bytes_t const & Z, v256_bytes_t const * mask) -> void;

auto align_rvv256_16_single(v256_bytes_t & S, v256_bytes_t * hep, v256_bytes_t const * const * qp,
                            v256_bytes_t const & gap_open, v256_bytes_t const & gap_extend,
                            long ql, v256_bytes_t const & Z, v256_bytes_t const * mask) -> void;


template <>
inline auto align_cells<Ops_7, No_mask>(v128_t & S,
                                        v128_t * const hep,
                                        v128_t const * const * const qp,
                                        GapVectors<Ops_7> const gaps,
                                        long const ql,
                                        v128_t const Z,
                                        No_mask const & /*masking*/) -> void
{
  align_rvv_7(S, hep, qp, gaps.open, gaps.extend, ql, Z, nullptr);
}

template <>
inline auto align_cells<Ops_7, Mask<Ops_7>>(v128_t & S,
                                            v128_t * const hep,
                                            v128_t const * const * const qp,
                                            GapVectors<Ops_7> const gaps,
                                            long const ql,
                                            v128_t const Z,
                                            Mask<Ops_7> const & masking) -> void
{
  align_rvv_7(S, hep, qp, gaps.open, gaps.extend, ql, Z, &masking.lanes);
}

template <>
inline auto align_cells<Ops_16, No_mask>(v128_t & S,
                                         v128_t * const hep,
                                         v128_t const * const * const qp,
                                         GapVectors<Ops_16> const gaps,
                                         long const ql,
                                         v128_t const Z,
                                         No_mask const & /*masking*/) -> void
{
  align_rvv_16(S, hep, qp, gaps.open, gaps.extend, ql, Z, nullptr);
}

template <>
inline auto align_cells<Ops_16, Mask<Ops_16>>(v128_t & S,
                                              v128_t * const hep,
                                              v128_t const * const * const qp,
                                              GapVectors<Ops_16> const gaps,
                                              long const ql,
                                              v128_t const Z,
                                              Mask<Ops_16> const & masking) -> void
{
  align_rvv_16(S, hep, qp, gaps.open, gaps.extend, ql, Z, &masking.lanes);
}

template <>
inline auto align_cells_single<Ops_16, No_mask>(v128_t & S,
                                                v128_t * const hep,
                                                v128_t const * const * const qp,
                                                GapVectors<Ops_16> const gaps,
                                                long const ql,
                                                v128_t const Z,
                                                No_mask const & /*masking*/) -> void
{
  align_rvv_16_single(S, hep, qp, gaps.open, gaps.extend, ql, Z, nullptr);
}

template <>
inline auto align_cells_single<Ops_16, Mask<Ops_16>>(v128_t & S,
                                                     v128_t * const hep,
                                                     v128_t const * const * const qp,
                                                     GapVectors<Ops_16> const gaps,
                                                     long const ql,
                                                     v128_t const Z,
                                                     Mask<Ops_16> const & masking) -> void
{
  align_rvv_16_single(S, hep, qp, gaps.open, gaps.extend, ql, Z, &masking.lanes);
}

#endif  // __riscv_vector

#endif  // SWIPE_ARCH_RISCV64_ALIGN_CELLS_RVV_H
