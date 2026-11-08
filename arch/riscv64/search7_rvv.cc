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

// The 7-bit kernel with RVV: the score profile built by table lookups
// (vrgather) in the transposed score matrix, as search7_ssse3() does
// with pshufb (arch/x86_64/SSSE3/search7_ssse3.cc), and the cells of
// arch/riscv64/align_rvv.cc. The search loop is shared with the other
// 128-bit kernels (search7_shared.h).

#include "swipe.h"
#include "search_kernels.h"  // the prototype of the kernel

#ifdef __riscv_vector

#include "search7_shared.h"  // search7_blocks(), channels_7_sse, depth_7, zero_score_7
#include "arch/riscv64/align_cells_rvv.h"  // align_rvv256_7(), v256_bytes_t
#include "view.h"  // View
#include <riscv_vector.h>
#include <array>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstring>  // std::memset
#include <iterator>  // std::next

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

constexpr std::size_t row_bytes = score_matrix_width;  // the table of a lookup
// the profile vectors of a query symbol (its row of the profile), and
// the database symbols of a block: depth_7 positions of 16 channels
// (32 channels in the 32-byte kernel: two such blocks)
constexpr std::size_t block_bytes = depth_7 * channels_7_sse;
constexpr auto row_stride = static_cast<std::ptrdiff_t>(row_bytes);
constexpr auto block_stride = static_cast<std::ptrdiff_t>(block_bytes);
// register group sizes (VLEN is at least 128 bits: LMUL = 4 holds 64
// bytes)
constexpr std::size_t lmul_2 = 2;
constexpr std::size_t lmul_4 = 4;

// The RVV operations of the profile builder, at one register group
// size (LMUL = 1, 2 or 4)
struct Lmul_1 {
  static auto load(BYTE const * const source, std::size_t const length) -> vuint8m1_t {
    return __riscv_vle8_v_u8m1(source, length);
  }
  static auto lookup(vuint8m1_t const table, vuint8m1_t const indices) -> vuint8m1_t {
    return __riscv_vrgather_vv_u8m1(table, indices, block_bytes);
  }
  static auto store(BYTE * const target, vuint8m1_t const vector) -> void {
    __riscv_vse8_v_u8m1(target, vector, block_bytes);
  }
};

struct Lmul_2 {
  static auto load(BYTE const * const source, std::size_t const length) -> vuint8m2_t {
    return __riscv_vle8_v_u8m2(source, length);
  }
  static auto lookup(vuint8m2_t const table, vuint8m2_t const indices) -> vuint8m2_t {
    return __riscv_vrgather_vv_u8m2(table, indices, block_bytes);
  }
  static auto store(BYTE * const target, vuint8m2_t const vector) -> void {
    __riscv_vse8_v_u8m2(target, vector, block_bytes);
  }
};

struct Lmul_4 {
  static auto load(BYTE const * const source, std::size_t const length) -> vuint8m4_t {
    return __riscv_vle8_v_u8m4(source, length);
  }
  static auto lookup(vuint8m4_t const table, vuint8m4_t const indices) -> vuint8m4_t {
    return __riscv_vrgather_vv_u8m4(table, indices, block_bytes);
  }
  static auto store(BYTE * const target, vuint8m4_t const vector) -> void {
    __riscv_vse8_v_u8m4(target, vector, block_bytes);
  }
};

// The profile vectors of query symbol 'row' (one per database position
// of the block, 64 bytes in all) are the row of the transposed score
// matrix (the scores of that query symbol against the 32 database
// symbols), indexed by the 64 database symbols of the block (dseq, in
// the same order): one lookup per row. vl is 32 for the table, 64 for
// the lookup.
//
// With 32 channels (blocks = 2), a profile row is 128 bytes, looked up
// in two halves of 64 bytes (the database symbols of positions 0 and 1,
// then of positions 2 and 3): a 128-byte lookup would need a register
// group twice as large, which costs as much as two lookups.
template <typename Lmul, std::ptrdiff_t blocks>
inline auto dprofile_gather7(BYTE * const dprofile,
                             BYTE const * const score_matrix,
                             BYTE const * const dseq) -> void
{
  static_assert(blocks == 1 or blocks == 2, "16 or 32 channels");
  constexpr std::ptrdiff_t profile_row_stride = blocks * block_stride;
  auto const symbols_0 = Lmul::load(dseq, block_bytes);
  auto const symbols_1 = Lmul::load(std::next(dseq, (blocks - 1) * block_stride), block_bytes);
  for (std::ptrdiff_t row = 0; row < static_cast<std::ptrdiff_t>(score_matrix_width); ++row)
  {
    auto const table = Lmul::load(std::next(score_matrix, row * row_stride), row_bytes);
    auto * const line = std::next(dprofile, row * profile_row_stride);
    Lmul::store(line, Lmul::lookup(table, symbols_0));
    if (blocks == 2)
    {
      Lmul::store(std::next(line, block_stride), Lmul::lookup(table, symbols_1));
    }
  }
}

// The smallest register group that holds the 64 bytes of a block: a
// lookup costs as much as its register group, not as its vl (BPI-F3,
// VLEN = 256: LMUL = 2 is 3.2 times faster than LMUL = 4, for the same
// result)
template <std::ptrdiff_t blocks>
inline auto dprofile_build7(BYTE * const dprofile,
                            BYTE const * const score_matrix,
                            BYTE const * const dseq) -> void
{
  static_assert(block_bytes == 64, "four positions of 16 channels");
  auto const register_bytes = __riscv_vlenb();  // VLEN / 8, at least 16
  if (register_bytes >= block_bytes)
  {
    dprofile_gather7<Lmul_1, blocks>(dprofile, score_matrix, dseq);
  }
  else if (lmul_2 * register_bytes >= block_bytes)
  {
    dprofile_gather7<Lmul_2, blocks>(dprofile, score_matrix, dseq);
  }
  else
  {
    assert(lmul_4 * register_bytes >= block_bytes);
    dprofile_gather7<Lmul_4, blocks>(dprofile, score_matrix, dseq);
  }
}

}  // anonymous namespace

auto search7_rvv(BYTE * * const q_start,
                 BYTE const gap_open_penalty,
                 BYTE const gap_extend_penalty,
                 BYTE const * const score_matrix,
                 BYTE * const dprofile,
                 BYTE * const hearray,
                 db_thread_s & dbt,
                 long const sequences,
                 long const * const seqnos,
                 long * const scores,
                 long const qlen) -> void
{
  auto const build_profile = [](BYTE * const profile,
                                BYTE const * const matrix,
                                BYTE const * const dseq) -> void
  {
    dprofile_build7<1>(profile, matrix, dseq);
  };
  search7_blocks(build_profile, q_start, gap_open_penalty, gap_extend_penalty,
                 score_matrix, dprofile, hearray, dbt, sequences, seqnos,
                 scores, qlen);
}



// The 32-channel version, when VLEN >= 256 (arch/riscv64/search_dispatch.cc):
// the search loop of search7_avx2() (arch/x86_64/AVX2/search7_avx2.cc),
// the vectors stored in 32-byte buffers between the calls of the cells
auto search7_rvv256(BYTE * * const q_start,
                    BYTE const gap_open_penalty,
                    BYTE const gap_extend_penalty,
                    BYTE const * const score_matrix,
                    BYTE * const dprofile,
                    BYTE * const hearray,
                    db_thread_s & dbt,
                    long const sequences,
                    long const * const seqnos,
                    long * const scores,
                    long const qlen) -> void
{
  constexpr std::size_t channels = channels_7_avx2;
  static_assert(channels == v256_bytes, "one byte lane per channel");
  static_assert(channels == 2 * channels_7_sse, "two 16-channel blocks");
  assert(__riscv_vlenb() >= v256_bytes);

  auto * const hep = reinterpret_cast<v256_bytes_t *>(hearray);
  auto const * const * const qp = reinterpret_cast<v256_bytes_t const * const *>(q_start);
  std::array<BYTE const *, channels> d_begin {{}};
  std::array<BYTE const *, channels> d_end {{}};

  // the database residues of the channels
  std::array<BYTE, depth_7 * channels> dseq {{}};
  // the lanes of the channels that start a new sequence (0x80)
  v256_bytes_t restart {{}};
  BYTE const zero = 0;

  std::array<long, channels> seq_id {{}};
  long next_id = 0;
  long done = 0;

  std::memset(hearray, zero_score_7, static_cast<std::size_t>(qlen) * hearray_row_bytes_avx2);

  v256_bytes_t Z {{}};
  Z.fill(zero_score_7);
  v256_bytes_t gap_open {{}};
  gap_open.fill(gap_open_penalty);
  v256_bytes_t gap_extend {{}};
  gap_extend.fill(gap_extend_penalty);

  // the scores of the channels, read lane by lane
  auto S = Z;

  d_begin.fill(&zero);
  d_end.fill(&zero);
  seq_id.fill(-1);

  // the next depth_7 residues of channel c (zeros after its end);
  // false when the channel reached the end of its sequence
  auto const fill_channel = [&](std::size_t const c) -> bool
  {
    for (std::size_t j = 0; j < depth_7; j++)
    {
      if (d_begin[c] < d_end[c])
      {
        dseq[(channels * j) + c] = *d_begin[c];
        d_begin[c] = std::next(d_begin[c]);
      }
      else
      {
        dseq[(channels * j) + c] = 0;
      }
    }
    return d_begin[c] != d_end[c];
  };

  bool easy = false;

  while (true)
  {
    if (easy)
    {
      // fill all channels
      for (std::size_t c = 0; c < channels; c++)
      {
        if (not fill_channel(c))
        {
          easy = false;
        }
      }

      dprofile_build7<2>(dprofile, score_matrix, dseq.data());
      align_rvv256_7(S, hep, qp, gap_open, gap_extend, qlen, Z, nullptr);
      continue;
    }

    // One or more sequences ended in the previous block: switch over
    // to new sequences
    easy = true;
    restart.fill(0);

    for (std::size_t c = 0; c < channels; c++)
    {
      if (d_begin[c] < d_end[c])
      {
        // this channel has more sequence
        if (not fill_channel(c))
        {
          easy = false;
        }
        continue;
      }

      // sequence in channel c ended: change of sequence
      restart[c] = zero_score_7;

      long const cand_id = seq_id[c];
      if (cand_id >= 0)
      {
        // save score
        scores[cand_id] = S[c] - zero_score_7;
        done++;
      }

      if (next_id < sequences)
      {
        // get next sequence
        seq_id[c] = next_id;
        long const seqnosf = seqnos[next_id];
        View<char> const sequence =
          db_getsequence(dbt, entry_seqno(seqnosf), entry_where(seqnosf), c);
        d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
        d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
        next_id++;
        if (not fill_channel(c))
        {
          easy = false;
        }
      }
      else
      {
        // no more sequences, empty channel
        seq_id[c] = -1;
        d_begin[c] = &zero;
        d_end[c] = d_begin[c];
        static_cast<void>(fill_channel(c));
      }
    }

    if (done == sequences)
    {
      break;
    }

    dprofile_build7<2>(dprofile, score_matrix, dseq.data());
    align_rvv256_7(S, hep, qp, gap_open, gap_extend, qlen, Z, &restart);
  }
}

#endif  // __riscv_vector
