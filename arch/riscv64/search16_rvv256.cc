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

// The 16-bit kernels with 32-byte vectors (16 database sequences at
// once), for riscv64 when VLEN >= 256 (arch/riscv64/search_dispatch.cc):
// search16_rvv256() and search16s_rvv256(), the search loops of
// search16_avx2() and search16s_avx2() (arch/x86_64/AVX2/search16_avx2.cc),
// with the cells of arch/riscv64/align_rvv.cc. The vectors are stored
// in 32-byte buffers between the calls of the cells; the comparison
// masks have one bit per lane (RVV mask registers stored with vsm). The
// score profile is built by table lookups (vrgather) in a transposed
// copy of the score matrix, as in arch/riscv64/search7_rvv.cc.

#include "swipe.h"
#include "search_kernels.h"  // the prototypes of the kernels

#ifdef __riscv_vector

#include "arch/riscv64/align_cells_rvv.h"  // align_rvv256_16(), align_rvv256_16_single(), v256_bytes_t
#include "view.h"  // View
#include <riscv_vector.h>
#include <array>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // std::int16_t, std::uint8_t, std::uintptr_t
#include <cstring>  // std::memcpy
#include <iterator>  // std::next, std::prev

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

constexpr std::size_t channels = channels_16_avx2;
static_assert(channels * sizeof(WORD) == v256_bytes, "one word lane per channel");

// the zero score of a 16-bit lane (the scores are stored with a bias
// of 0x8000), which is also the lane value that resets a score to zero
// (Ops_16::mask(), align_cells.h)
constexpr WORD zero_score = 0x8000;

// the lanes of a 32-byte vector, as words; aligned on 32 bytes, as the
// word loads of the cells need (at least 2 bytes)
struct alignas(v256_bytes) Lanes16
{
  v256_bytes_t bytes {{}};

  auto lane(std::size_t const c) const -> WORD
  {
    WORD word = 0;
    std::memcpy(&word, &bytes[c * sizeof(WORD)], sizeof(WORD));
    return word;
  }

  auto set_lane(std::size_t const c, WORD const word) -> void
  {
    std::memcpy(&bytes[c * sizeof(WORD)], &word, sizeof(WORD));
  }

  auto fill(WORD const word) -> void
  {
    for (std::size_t c = 0; c < channels; c++)
    {
      set_lane(c, word);
    }
  }
};

// the words of a 32-byte buffer (a cast through void *: -Wcast-align)
inline auto as_words(v256_bytes_t & vector) -> std::int16_t *
{
  assert(reinterpret_cast<std::uintptr_t>(vector.data()) % alignof(std::int16_t) == 0);
  return static_cast<std::int16_t *>(static_cast<void *>(vector.data()));
}

inline auto as_words(v256_bytes_t const & vector) -> std::int16_t const *
{
  assert(reinterpret_cast<std::uintptr_t>(vector.data()) % alignof(std::int16_t) == 0);
  return static_cast<std::int16_t const *>(static_cast<void const *>(vector.data()));
}

inline auto load16(v256_bytes_t const & vector) -> vint16m1_t
{
  return __riscv_vle16_v_i16m1(as_words(vector), channels);
}

// the lanes set in an RVV mask, one bit per lane (lane c: bit c)
inline auto mask_bits(vbool16_t const mask) -> unsigned int
{
  std::array<std::uint8_t, channels / 8> bits {{}};
  __riscv_vsm_v_b16(bits.data(), mask, channels);
  return bits[0] | (static_cast<unsigned int>(bits[1]) << 8U);
}

// the lanes where lhs > rhs (signed words)
inline auto mask_gt(v256_bytes_t const & lhs, v256_bytes_t const & rhs) -> unsigned int
{
  return mask_bits(__riscv_vmsgt_vv_i16m1_b16(load16(lhs), load16(rhs), channels));
}

// the lanes where lhs == rhs
inline auto mask_eq(v256_bytes_t const & lhs, v256_bytes_t const & rhs) -> unsigned int
{
  return mask_bits(__riscv_vmseq_vv_i16m1_b16(load16(lhs), load16(rhs), channels));
}

inline auto lane_is_set(unsigned int const mask, std::size_t const c) -> bool
{
  return ((mask >> c) & 1U) != 0;
}

// lhs + rhs, signed saturated (v256_adds_i16)
inline auto adds16(Lanes16 const & lhs, Lanes16 const & rhs) -> Lanes16
{
  Lanes16 sum;
  __riscv_vse16_v_i16m1(as_words(sum.bytes),
                        __riscv_vsadd_vv_i16m1(load16(lhs.bytes), load16(rhs.bytes), channels),
                        channels);
  return sum;
}


// The score profile, by table lookups: the transposed score matrix
// (row q: the scores of query symbol q against the 32 database
// symbols), built once per call of a kernel
using Transposed16 = std::array<WORD, score_matrix_size>;

inline auto transpose16(WORD const * const score_matrix, Transposed16 & transposed) -> void
{
  for (std::size_t row = 0; row < score_matrix_width; row++)
  {
    for (std::size_t column = 0; column < score_matrix_width; column++)
    {
      transposed[score_matrix_cell(column, row)] = score_matrix[score_matrix_cell(row, column)];
    }
  }
}

// the lookups of 32 words (a table of 64 bytes), at LMUL = 1 (VLEN >=
// 512) or 2 (VLEN = 256)
constexpr std::size_t table_words = score_matrix_width;

struct Lmul16_1 {
  using Indices = vuint16m1_t;
  static auto indices(BYTE const * const symbols, std::size_t const length) -> Indices {
    return __riscv_vzext_vf2_u16m1(__riscv_vle8_v_u8mf2(symbols, length), length);
  }
  static auto table(WORD const * const row) -> vuint16m1_t {
    return __riscv_vle16_v_u16m1(row, table_words);
  }
  static auto lookup(WORD * const target, vuint16m1_t const table, Indices const indices,
                     std::size_t const length) -> void {
    __riscv_vse16_v_u16m1(target, __riscv_vrgather_vv_u16m1(table, indices, length), length);
  }
};

struct Lmul16_2 {
  using Indices = vuint16m2_t;
  static auto indices(BYTE const * const symbols, std::size_t const length) -> Indices {
    return __riscv_vzext_vf2_u16m2(__riscv_vle8_v_u8m1(symbols, length), length);
  }
  static auto table(WORD const * const row) -> vuint16m2_t {
    return __riscv_vle16_v_u16m2(row, table_words);
  }
  static auto lookup(WORD * const target, vuint16m2_t const table, Indices const indices,
                     std::size_t const length) -> void {
    __riscv_vse16_v_u16m2(target, __riscv_vrgather_vv_u16m2(table, indices, length), length);
  }
};

// The profile vectors of query symbol 'row' (depth x 16 words: one
// vector per database position) are the row of the transposed matrix,
// indexed by the database symbols of the block (dseq, in the same
// order), in lookups of at most 32 words
template <typename Lmul, std::size_t depth>
inline auto dprofile_gather16(WORD * const dprofile,
                              Transposed16 const & transposed,
                              BYTE const * const dseq) -> void
{
  static_assert(depth == 1 or depth == 4, "the depths of search16s() and search16()");
  constexpr std::size_t row_words = depth * channels;
  constexpr std::size_t lookup_words = depth == 1 ? channels : table_words;
  constexpr auto lookup_stride = static_cast<std::ptrdiff_t>(lookup_words);
  static_assert(row_words % lookup_words == 0, "whole lookups");
  auto const indices_0 = Lmul::indices(dseq, lookup_words);
  auto const indices_1 = Lmul::indices(std::next(dseq, (depth == 1) ? 0 : lookup_stride), lookup_words);
  for (std::size_t row = 0; row < score_matrix_width; row++)
  {
    auto const table = Lmul::table(&transposed[row * table_words]);
    auto * const line = std::next(dprofile, static_cast<std::ptrdiff_t>(row * row_words));
    Lmul::lookup(line, table, indices_0, lookup_words);
    if (depth == 4)
    {
      Lmul::lookup(std::next(line, lookup_stride), table, indices_1, lookup_words);
    }
  }
}

template <std::size_t depth>
inline auto dprofile_build16(WORD * const dprofile,
                             Transposed16 const & transposed,
                             BYTE const * const dseq) -> void
{
  constexpr std::size_t table_bytes = table_words * sizeof(WORD);
  if (__riscv_vlenb() >= table_bytes)
  {
    dprofile_gather16<Lmul16_1, depth>(dprofile, transposed, dseq);
    return;
  }
  dprofile_gather16<Lmul16_2, depth>(dprofile, transposed, dseq);
}

}  // anonymous namespace


auto search16_rvv256(WORD * * const q_start,
                     WORD const gap_open_penalty,
                     WORD const gap_extend_penalty,
                     WORD const * const score_matrix,
                     WORD * const dprofile,
                     WORD * const hearray,
                     db_thread_s & dbt,
                     long const sequences,
                     long const * const seqnos,
                     long * const scores,
                     long * const bestpos,
                     int const qlen) -> void
{
  constexpr std::size_t CDEPTH = 4;
  assert(__riscv_vlenb() >= v256_bytes);
  auto * const hep = reinterpret_cast<v256_bytes_t *>(hearray);
  auto const * const * const qp = reinterpret_cast<v256_bytes_t const * const *>(q_start);
  std::array<BYTE const *, channels> d_begin {{}};
  std::array<BYTE const *, channels> d_pos {{}};
  std::array<BYTE const *, channels> d_best {{}};
  std::array<BYTE const *, channels> d_end {{}};

  // the database residues of the channels
  std::array<BYTE, CDEPTH * channels> dseq {{}};
  BYTE const zero = 0;

  std::array<long, channels> seq_id {{}};
  long next_id = 0;
  long done = 0;

  Transposed16 transposed {{}};
  transpose16(score_matrix, transposed);

  Lanes16 Z;
  Z.fill(zero_score);
  Lanes16 Q;
  Q.fill(gap_open_penalty);
  Lanes16 R;
  R.fill(gap_extend_penalty);

  auto S = Z;
  auto SL = Z;

  // the H and E scores of each query residue
  for (long i = 0; i < 2L * qlen; i++)
  {
    *std::next(hep, i) = Z.bytes;
  }

  d_begin.fill(&zero);
  d_pos.fill(&zero);
  d_best.fill(&zero);
  d_end.fill(&zero);
  seq_id.fill(-1);

  // the next CDEPTH residues of channel c (zeros after its end)
  auto const fill_channel = [&](std::size_t const c) -> void
  {
    for (std::size_t j = 0; j < CDEPTH; j++)
    {
      if (d_pos[c] < d_end[c])
      {
        dseq[(channels * j) + c] = *d_pos[c];
        d_pos[c] = std::next(d_pos[c]);
      }
      else
      {
        dseq[(channels * j) + c] = 0;
      }
    }
  };

  // save column address if new highscore
  auto const save_best_positions = [&](unsigned int const mask) -> void
  {
    for (std::size_t c = 0; c < channels; c++)
    {
      if (lane_is_set(mask, c))
      {
        d_best[c] = d_pos[c];
      }
    }
  };

  bool easy = false;

  while (true)
  {
    if (easy)
    {
      for (std::size_t c = 0; c < channels; c++)
      {
        fill_channel(c);
        if ((d_pos[c] == d_end[c]) && (seq_id[c] > -1))
        {
          easy = false;
        }
      }

      dprofile_build16<CDEPTH>(dprofile, transposed, dseq.data());
      align_rvv256_16(S.bytes, hep, qp, Q.bytes, R.bytes, qlen, Z.bytes, nullptr);
      save_best_positions(mask_gt(S.bytes, SL.bytes));
      SL = S;
      continue;
    }

    easy = true;
    Lanes16 restart;
    auto const lane_scores = S;

    for (std::size_t c = 0; c < channels; c++)
    {
      if (d_pos[c] < d_end[c])
      {
        fill_channel(c);
        if (d_pos[c] == d_end[c])
        {
          easy = false;
        }
        continue;
      }

      restart.set_lane(c, zero_score);

      long const cand_id = seq_id[c];
      if (cand_id >= 0)
      {
        long const score = lane_scores.lane(c) ^ zero_score;
        scores[cand_id] = score;
        bestpos[cand_id] = d_best[c] - d_begin[c];
        done++;
      }

      if (next_id < sequences)
      {
        seq_id[c] = next_id;
        long const seqnosf = seqnos[next_id];
        View<char> const sequence =
          db_getsequence(dbt, entry_seqno(seqnosf), entry_where(seqnosf), c);
        d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
        d_pos[c] = d_begin[c];
        d_best[c] = d_begin[c];
        d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
        next_id++;
        fill_channel(c);
        if (d_pos[c] == d_end[c])
        {
          easy = false;
        }
      }
      else
      {
        seq_id[c] = -1;
        d_begin[c] = &zero;
        d_pos[c] = d_begin[c];
        d_best[c] = d_begin[c];
        d_end[c] = d_begin[c];
        fill_channel(c);
      }
    }

    if (done == sequences)
    {
      break;
    }

    dprofile_build16<CDEPTH>(dprofile, transposed, dseq.data());
    align_rvv256_16(S.bytes, hep, qp, Q.bytes, R.bytes, qlen, Z.bytes, &restart.bytes);

    SL = adds16(SL, restart);
    SL = adds16(SL, restart);
    save_best_positions(mask_gt(S.bytes, SL.bytes));
    SL = S;
  }
}


auto search16s_rvv256(WORD * * const q_start,
                      WORD const gap_open_penalty,
                      WORD const gap_extend_penalty,
                      WORD const * const score_matrix,
                      WORD * const dprofile,
                      WORD * const hearray,
                      DbThread const * const dbta,
                      long const sequences,
                      long const * const seqnos,
                      long * const scores,
                      long * const bestpos,
                      long * const bestq,
                      int const qlen) -> void
{
  constexpr std::size_t CDEPTH = 1;
  assert(__riscv_vlenb() >= v256_bytes);
  auto * const hep = reinterpret_cast<v256_bytes_t *>(hearray);
  auto const * const * const qp = reinterpret_cast<v256_bytes_t const * const *>(q_start);
  std::array<BYTE const *, channels> d_begin {{}};
  std::array<BYTE const *, channels> d_pos {{}};
  std::array<BYTE const *, channels> d_end {{}};
  std::array<BYTE const *, channels> d_best {{}};
  std::array<long, channels> q_best {{}};

  // the database residues of the channels
  std::array<BYTE, CDEPTH * channels> dseq {{}};
  BYTE const zero = 0;

  std::array<long, channels> seq_id {{}};
  long next_id = 0;
  long done = 0;

  Transposed16 transposed {{}};
  transpose16(score_matrix, transposed);

  Lanes16 Z;
  Z.fill(zero_score);
  Lanes16 Q;
  Q.fill(gap_open_penalty);
  Lanes16 R;
  R.fill(gap_extend_penalty);

  auto S = Z;
  auto SL = Z;

  // the H and E scores of each query residue
  for (long i = 0; i < 2L * qlen; i++)
  {
    *std::next(hep, i) = Z.bytes;
  }

  d_begin.fill(&zero);
  d_pos.fill(&zero);
  d_best.fill(&zero);
  d_end.fill(&zero);
  q_best.fill(-1);
  seq_id.fill(-1);

  // the next residue of channel c (zero after its end)
  auto const fill_channel = [&](std::size_t const c) -> void
  {
    for (std::size_t j = 0; j < CDEPTH; j++)
    {
      if (d_pos[c] < d_end[c])
      {
        dseq[(channels * j) + c] = *d_pos[c];
        d_pos[c] = std::next(d_pos[c]);
      }
      else
      {
        dseq[(channels * j) + c] = 0;
      }
    }
  };

  // save column address if new highscore, and the query position of
  // that score (the first one, scanning backwards)
  auto const save_best_positions = [&](unsigned int const mask) -> void
  {
    if (mask == 0)
    {
      return;
    }
    for (std::size_t c = 0; c < channels; c++)
    {
      if (lane_is_set(mask, c))
      {
        d_best[c] = std::prev(d_pos[c]);
      }
    }
    for (long i = qlen - 1; i >= 0; i--)
    {
      auto const m2 = mask & mask_eq(*std::next(hep, 2 * i), S.bytes);
      if (m2 == 0)
      {
        continue;
      }
      for (std::size_t c = 0; c < channels; c++)
      {
        if (lane_is_set(m2, c))
        {
          q_best[c] = i;
        }
      }
    }
  };

  bool easy = false;

  while (true)
  {
    if (easy)
    {
      for (std::size_t c = 0; c < channels; c++)
      {
        fill_channel(c);
        if ((d_pos[c] == d_end[c]) && (seq_id[c] > -1))
        {
          easy = false;
        }
      }

      dprofile_build16<CDEPTH>(dprofile, transposed, dseq.data());
      align_rvv256_16_single(S.bytes, hep, qp, Q.bytes, R.bytes, qlen, Z.bytes, nullptr);
      save_best_positions(mask_gt(S.bytes, SL.bytes));
      SL = S;
      continue;
    }

    easy = true;
    Lanes16 restart;
    auto const lane_scores = S;

    for (std::size_t c = 0; c < channels; c++)
    {
      if (d_pos[c] < d_end[c])
      {
        fill_channel(c);
        if (d_pos[c] == d_end[c])
        {
          easy = false;
        }
        continue;
      }

      restart.set_lane(c, zero_score);

      long const cand_id = seq_id[c];
      if (cand_id >= 0)
      {
        long const score = lane_scores.lane(c) ^ zero_score;
        scores[cand_id] = score;
        bestpos[cand_id] = d_best[c] - d_begin[c];
        bestq[cand_id] = q_best[c];
        done++;
      }

      if (next_id < sequences)
      {
        seq_id[c] = next_id;
        long const seqnosf = seqnos[next_id];
        long const seqno = entry_seqno(seqnosf);
        db_mapsequences(*dbta[c], seqno, seqno);
        View<char> const sequence =
          db_getsequence(*dbta[c], seqno, entry_where(seqnosf), c);
        d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
        d_pos[c] = d_begin[c];
        d_best[c] = d_begin[c];
        d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
        q_best[c] = -1;
        next_id++;
        fill_channel(c);
        if (d_pos[c] == d_end[c])
        {
          easy = false;
        }
      }
      else
      {
        seq_id[c] = -1;
        d_pos[c] = &zero;
        d_end[c] = d_pos[c];
        fill_channel(c);
      }
    }

    if (done == sequences)
    {
      break;
    }

    dprofile_build16<CDEPTH>(dprofile, transposed, dseq.data());
    align_rvv256_16_single(S.bytes, hep, qp, Q.bytes, R.bytes, qlen, Z.bytes, &restart.bytes);

    SL = adds16(SL, restart);
    SL = adds16(SL, restart);
    save_best_positions(mask_gt(S.bytes, SL.bytes));
    SL = S;
  }
}

#endif  // __riscv_vector
