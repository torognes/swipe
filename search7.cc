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

#include "swipe.h"
#include "align_cells.h"  // Ops_7, onestep(), No_mask, Mask
#include <array>
#include <cstddef>  // std::ptrdiff_t, std::size_t

constexpr std::size_t CHANNELS = channels_7;
static_assert(sizeof(__m128i) == vector_bytes, "an SSE vector");
constexpr std::size_t CDEPTH = 4;

// the byte 0x80 (the lanes of _mm_set_epi8() are char: 0x80 does not
// fit in a signed char, -128 has the same bits)
constexpr char byte_0x80 = static_cast<char>(-128);

#ifdef SWIPE_SSSE3

// profline(j) strides, in 16-byte vectors: a row of the 32 x 32 score
// matrix is two vectors, a row of the profile one vector per CDEPTH
constexpr std::ptrdiff_t matrix_row_vectors = 2;
constexpr auto profile_row_vectors = static_cast<std::ptrdiff_t>(CDEPTH);

inline auto dprofile_shuffle7(BYTE * dprofile,
			      BYTE * score_matrix,
			      BYTE * dseq_byte) -> void
{
  __m128i a;
  __m128i b;
  __m128i c;
  __m128i d;
  __m128i x;
  __m128i y;
  __m128i m0;
  __m128i m1;
  __m128i m2;
  __m128i m3;
  __m128i m4;
  __m128i m5;
  __m128i m6;
  __m128i m7;
  __m128i t0;
  __m128i t1;
  __m128i t2;
  __m128i t3;
  __m128i t4;
  __m128i t5;
  __m128i t6;
  __m128i t7;
  __m128i t8;
  __m128i t9;
  __m128i t10;
  __m128i t11;
  __m128i t12;
  __m128i t13;
  __m128i u0, u1, u2, u3, u4, u5,         u8, u9, u10, u11, u12, u13;

  auto * dseq = reinterpret_cast<__m128i*>(dseq_byte);
  
  // 16 x 4 = 64 db symbols
  // ca 458 instructions

  // make masks

  /* Note: pshufb only on modern Intel cpus (SSSE3), not AMD */
  /* SSSE3: Supplemental SSE3 */

  x = _mm_set_epi8(0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10,
                   0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10);

  y = _mm_set1_epi8(byte_0x80);

  a  = _mm_load_si128(dseq);
  t0 = _mm_and_si128(a, x);
  t1 = _mm_slli_epi16(t0, 3);
  t2 = _mm_xor_si128(t1, y);
  m0 = _mm_or_si128(a, t1);
  m1 = _mm_or_si128(a, t2);

  b  = _mm_load_si128(dseq+1);
  t3 = _mm_and_si128(b, x);
  t4 = _mm_slli_epi16(t3, 3);
  t5 = _mm_xor_si128(t4, y);
  m2 = _mm_or_si128(b, t4);
  m3 = _mm_or_si128(b, t5);

  c  = _mm_load_si128(dseq+2);
  u0 = _mm_and_si128(c, x);
  u1 = _mm_slli_epi16(u0, 3);
  u2 = _mm_xor_si128(u1, y);
  m4 = _mm_or_si128(c, u1);
  m5 = _mm_or_si128(c, u2);

  d  = _mm_load_si128(dseq+3);
  u3 = _mm_and_si128(d, x);
  u4 = _mm_slli_epi16(u3, 3);
  u5 = _mm_xor_si128(u4, y);
  m6 = _mm_or_si128(d, u4);
  m7 = _mm_or_si128(d, u5);

#define profline(j)					\
  t6  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix)+(matrix_row_vectors*(j)));   \
  t7  = _mm_load_si128(reinterpret_cast<__m128i*>(score_matrix)+(matrix_row_vectors*(j))+1); \
  t8  = _mm_shuffle_epi8(t6, m0);			\
  t9  = _mm_shuffle_epi8(t7, m1);			\
  t10 = _mm_shuffle_epi8(t6, m2);			\
  t11 = _mm_shuffle_epi8(t7, m3);			\
  u8  = _mm_shuffle_epi8(t6, m4);			\
  u9  = _mm_shuffle_epi8(t7, m5);			\
  u10 = _mm_shuffle_epi8(t6, m6);			\
  u11 = _mm_shuffle_epi8(t7, m7);			\
  t12 = _mm_or_si128(t8,  t9);				\
  t13 = _mm_or_si128(t10, t11);				\
  u12 = _mm_or_si128(u8,  u9);				\
  u13 = _mm_or_si128(u10, u11);				\
  _mm_store_si128(reinterpret_cast<__m128i*>(dprofile)+(profile_row_vectors*(j)),   t12);	\
  _mm_store_si128(reinterpret_cast<__m128i*>(dprofile)+(profile_row_vectors*(j))+1, t13);	\
  _mm_store_si128(reinterpret_cast<__m128i*>(dprofile)+(profile_row_vectors*(j))+2, u12);	\
  _mm_store_si128(reinterpret_cast<__m128i*>(dprofile)+(profile_row_vectors*(j))+3, u13)

  profline(0);
  profline(1);
  profline(2);
  profline(3);
  profline(4);
  profline(5);
  profline(6);
  profline(7);
  profline(8);
  profline(9);
  profline(10);
  profline(11);
  profline(12);
  profline(13);
  profline(14);
  profline(15);
  profline(16);
  profline(17);
  profline(18);
  profline(19);
  profline(20);
  profline(21);
  profline(22);
  profline(23);
  profline(24);
  profline(25);
  profline(26);
  profline(27);
  profline(28);
  profline(29);
  profline(30);
  profline(31);
}

#else

inline auto dprofile_fill7(BYTE * dprofile,
			   BYTE * score_matrix,
			   BYTE const * dseq) -> void
{
  __m128i xmm0;
  __m128i xmm1;
  __m128i xmm2;
  __m128i xmm3;
  __m128i xmm4;
  __m128i xmm5;
  __m128i xmm6;
  __m128i xmm7;
  __m128i xmm8;
  __m128i xmm9;
  __m128i xmm10;
  __m128i xmm11;
  __m128i xmm12;
  __m128i xmm13;
  __m128i xmm14;
  __m128i xmm15;
  
  // 4 x 16 db symbols
  // ca (60x2+68x2)x4 = 976 instructions

  for (std::size_t j = 0; j < CDEPTH; j++)
  {
    std::array<unsigned, CHANNELS> d;
    for (std::size_t i = 0; i < CHANNELS; i++)
    {
      d[i] = static_cast<unsigned>(dseq[(j * CHANNELS) + i]) << 5;
    }

    xmm0  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + d[0] ));
    xmm2  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + d[2] ));
    xmm4  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + d[4] ));
    xmm6  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + d[6] ));
    xmm8  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + d[8] ));
    xmm10 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + d[10]));
    xmm12 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + d[12]));
    xmm14 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + d[14]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  *reinterpret_cast<__m128i*>(score_matrix + d[1] ));
    xmm2  = _mm_unpacklo_epi8(xmm2,  *reinterpret_cast<__m128i*>(score_matrix + d[3] ));
    xmm4  = _mm_unpacklo_epi8(xmm4,  *reinterpret_cast<__m128i*>(score_matrix + d[5] ));
    xmm6  = _mm_unpacklo_epi8(xmm6,  *reinterpret_cast<__m128i*>(score_matrix + d[7] ));
    xmm8  = _mm_unpacklo_epi8(xmm8,  *reinterpret_cast<__m128i*>(score_matrix + d[9] ));
    xmm10 = _mm_unpacklo_epi8(xmm10, *reinterpret_cast<__m128i*>(score_matrix + d[11]));
    xmm12 = _mm_unpacklo_epi8(xmm12, *reinterpret_cast<__m128i*>(score_matrix + d[13]));
    xmm14 = _mm_unpacklo_epi8(xmm14, *reinterpret_cast<__m128i*>(score_matrix + d[15]));
      
    xmm1 = xmm0;
    xmm0 = _mm_unpacklo_epi16(xmm0, xmm2);
    xmm1 = _mm_unpackhi_epi16(xmm1, xmm2);
    xmm5 = xmm4;
    xmm4 = _mm_unpacklo_epi16(xmm4, xmm6);
    xmm5 = _mm_unpackhi_epi16(xmm5, xmm6);
    xmm9 = xmm8;
    xmm8 = _mm_unpacklo_epi16(xmm8, xmm10);
    xmm9 = _mm_unpackhi_epi16(xmm9, xmm10);
    xmm13 = xmm12;
    xmm12 = _mm_unpacklo_epi16(xmm12, xmm14);
    xmm13 = _mm_unpackhi_epi16(xmm13, xmm14);

    xmm2  = xmm0;
    xmm0  = _mm_unpacklo_epi32(xmm0, xmm4);
    xmm2  = _mm_unpackhi_epi32(xmm2, xmm4);
    xmm6  = xmm1;
    xmm1  = _mm_unpacklo_epi32(xmm1, xmm5);
    xmm6  = _mm_unpackhi_epi32(xmm6, xmm5);
    xmm10 = xmm8;
    xmm8  = _mm_unpacklo_epi32(xmm8, xmm12);
    xmm10 = _mm_unpackhi_epi32(xmm10, xmm12);
    xmm14 = xmm9;
    xmm9  = _mm_unpacklo_epi32(xmm9, xmm13);
    xmm14 = _mm_unpackhi_epi32(xmm14, xmm13);
      
    xmm3  = xmm0;
    xmm0  = _mm_unpacklo_epi64(xmm0, xmm8);
    xmm3  = _mm_unpackhi_epi64(xmm3, xmm8);
    xmm7  = xmm2;
    xmm2  = _mm_unpacklo_epi64(xmm2, xmm10);
    xmm7  = _mm_unpackhi_epi64(xmm7, xmm10);
    xmm11 = xmm1;
    xmm1  = _mm_unpacklo_epi64(xmm1, xmm9);
    xmm11 = _mm_unpackhi_epi64(xmm11, xmm9);
    xmm15 = xmm6;
    xmm6  = _mm_unpacklo_epi64(xmm6, xmm14);
    xmm15 = _mm_unpackhi_epi64(xmm15, xmm14);

    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+  0), xmm0);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+ 64), xmm3);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+128), xmm2);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+192), xmm7);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+256), xmm1);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+320), xmm11);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+384), xmm6);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+448), xmm15);


    // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

    xmm0  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[0 ]));
    xmm1  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[1 ]));
    xmm2  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[2 ]));
    xmm3  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[3 ]));
    xmm4  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[4 ]));
    xmm5  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[5 ]));
    xmm6  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[6 ]));
    xmm7  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[7 ]));
    xmm8  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[8 ]));
    xmm9  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[9 ]));
    xmm10 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[10]));
    xmm11 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[11]));
    xmm12 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[12]));
    xmm13 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[13]));
    xmm14 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[14]));
    xmm15 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 8 + d[15]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  xmm1);
    xmm2  = _mm_unpacklo_epi8(xmm2,  xmm3);
    xmm4  = _mm_unpacklo_epi8(xmm4,  xmm5);
    xmm6  = _mm_unpacklo_epi8(xmm6,  xmm7);
    xmm8  = _mm_unpacklo_epi8(xmm8,  xmm9);
    xmm10 = _mm_unpacklo_epi8(xmm10, xmm11);
    xmm12 = _mm_unpacklo_epi8(xmm12, xmm13);
    xmm14 = _mm_unpacklo_epi8(xmm14, xmm15);
      
    xmm1 = xmm0;
    xmm0 = _mm_unpacklo_epi16(xmm0, xmm2);
    xmm1 = _mm_unpackhi_epi16(xmm1, xmm2);
    xmm5 = xmm4;
    xmm4 = _mm_unpacklo_epi16(xmm4, xmm6);
    xmm5 = _mm_unpackhi_epi16(xmm5, xmm6);
    xmm9 = xmm8;
    xmm8 = _mm_unpacklo_epi16(xmm8, xmm10);
    xmm9 = _mm_unpackhi_epi16(xmm9, xmm10);
    xmm13 = xmm12;
    xmm12 = _mm_unpacklo_epi16(xmm12, xmm14);
    xmm13 = _mm_unpackhi_epi16(xmm13, xmm14);

    xmm2  = xmm0;
    xmm0  = _mm_unpacklo_epi32(xmm0, xmm4);
    xmm2  = _mm_unpackhi_epi32(xmm2, xmm4);
    xmm6  = xmm1;
    xmm1  = _mm_unpacklo_epi32(xmm1, xmm5);
    xmm6  = _mm_unpackhi_epi32(xmm6, xmm5);
    xmm10 = xmm8;
    xmm8  = _mm_unpacklo_epi32(xmm8, xmm12);
    xmm10 = _mm_unpackhi_epi32(xmm10, xmm12);
    xmm14 = xmm9;
    xmm9  = _mm_unpacklo_epi32(xmm9, xmm13);
    xmm14 = _mm_unpackhi_epi32(xmm14, xmm13);
      
    xmm3  = xmm0;
    xmm0  = _mm_unpacklo_epi64(xmm0, xmm8);
    xmm3  = _mm_unpackhi_epi64(xmm3, xmm8);
    xmm7  = xmm2;
    xmm2  = _mm_unpacklo_epi64(xmm2, xmm10);
    xmm7  = _mm_unpackhi_epi64(xmm7, xmm10);
    xmm11 = xmm1;
    xmm1  = _mm_unpacklo_epi64(xmm1, xmm9);
    xmm11 = _mm_unpackhi_epi64(xmm11, xmm9);
    xmm15 = xmm6;
    xmm6  = _mm_unpacklo_epi64(xmm6, xmm14);
    xmm15 = _mm_unpackhi_epi64(xmm15, xmm14);

    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+512+  0), xmm0);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+512+ 64), xmm3);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+512+128), xmm2);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+512+192), xmm7);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+512+256), xmm1);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+512+320), xmm11);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+512+384), xmm6);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+512+448), xmm15);


    xmm0  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 16 + d[0 ]));
    xmm2  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 16 + d[2 ]));
    xmm4  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 16 + d[4 ]));
    xmm6  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 16 + d[6 ]));
    xmm8  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 16 + d[8 ]));
    xmm10 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 16 + d[10]));
    xmm12 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 16 + d[12]));
    xmm14 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 16 + d[14]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  *reinterpret_cast<__m128i*>(score_matrix + 16 + d[1 ]));
    xmm2  = _mm_unpacklo_epi8(xmm2,  *reinterpret_cast<__m128i*>(score_matrix + 16 + d[3 ]));
    xmm4  = _mm_unpacklo_epi8(xmm4,  *reinterpret_cast<__m128i*>(score_matrix + 16 + d[5 ]));
    xmm6  = _mm_unpacklo_epi8(xmm6,  *reinterpret_cast<__m128i*>(score_matrix + 16 + d[7 ]));
    xmm8  = _mm_unpacklo_epi8(xmm8,  *reinterpret_cast<__m128i*>(score_matrix + 16 + d[9 ]));
    xmm10 = _mm_unpacklo_epi8(xmm10, *reinterpret_cast<__m128i*>(score_matrix + 16 + d[11 ]));
    xmm12 = _mm_unpacklo_epi8(xmm12, *reinterpret_cast<__m128i*>(score_matrix + 16 + d[13 ]));
    xmm14 = _mm_unpacklo_epi8(xmm14, *reinterpret_cast<__m128i*>(score_matrix + 16 + d[15 ]));
      
    xmm1 = xmm0;
    xmm0 = _mm_unpacklo_epi16(xmm0, xmm2);
    xmm1 = _mm_unpackhi_epi16(xmm1, xmm2);
    xmm5 = xmm4;
    xmm4 = _mm_unpacklo_epi16(xmm4, xmm6);
    xmm5 = _mm_unpackhi_epi16(xmm5, xmm6);
    xmm9 = xmm8;
    xmm8 = _mm_unpacklo_epi16(xmm8, xmm10);
    xmm9 = _mm_unpackhi_epi16(xmm9, xmm10);
    xmm13 = xmm12;
    xmm12 = _mm_unpacklo_epi16(xmm12, xmm14);
    xmm13 = _mm_unpackhi_epi16(xmm13, xmm14);

    xmm2  = xmm0;
    xmm0  = _mm_unpacklo_epi32(xmm0, xmm4);
    xmm2  = _mm_unpackhi_epi32(xmm2, xmm4);
    xmm6  = xmm1;
    xmm1  = _mm_unpacklo_epi32(xmm1, xmm5);
    xmm6  = _mm_unpackhi_epi32(xmm6, xmm5);
    xmm10 = xmm8;
    xmm8  = _mm_unpacklo_epi32(xmm8, xmm12);
    xmm10 = _mm_unpackhi_epi32(xmm10, xmm12);
    xmm14 = xmm9;
    xmm9  = _mm_unpacklo_epi32(xmm9, xmm13);
    xmm14 = _mm_unpackhi_epi32(xmm14, xmm13);
      
    xmm3  = xmm0;
    xmm0  = _mm_unpacklo_epi64(xmm0, xmm8);
    xmm3  = _mm_unpackhi_epi64(xmm3, xmm8);
    xmm7  = xmm2;
    xmm2  = _mm_unpacklo_epi64(xmm2, xmm10);
    xmm7  = _mm_unpackhi_epi64(xmm7, xmm10);
    xmm11 = xmm1;
    xmm1  = _mm_unpacklo_epi64(xmm1, xmm9);
    xmm11 = _mm_unpackhi_epi64(xmm11, xmm9);
    xmm15 = xmm6;
    xmm6  = _mm_unpacklo_epi64(xmm6, xmm14);
    xmm15 = _mm_unpackhi_epi64(xmm15, xmm14);

    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1024+  0), xmm0);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1024+ 64), xmm3);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1024+128), xmm2);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1024+192), xmm7);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1024+256), xmm1);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1024+320), xmm11);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1024+384), xmm6);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1024+448), xmm15);


    // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

    xmm0  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[0 ]));
    xmm1  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[1 ]));
    xmm2  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[2 ]));
    xmm3  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[3 ]));
    xmm4  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[4 ]));
    xmm5  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[5 ]));
    xmm6  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[6 ]));
    xmm7  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[7 ]));
    xmm8  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[8 ]));
    xmm9  = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[9 ]));
    xmm10 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[10]));
    xmm11 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[11]));
    xmm12 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[12]));
    xmm13 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[13]));
    xmm14 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[14]));
    xmm15 = _mm_loadl_epi64(reinterpret_cast<__m128i*>(score_matrix + 24 + d[15]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  xmm1);
    xmm2  = _mm_unpacklo_epi8(xmm2,  xmm3);
    xmm4  = _mm_unpacklo_epi8(xmm4,  xmm5);
    xmm6  = _mm_unpacklo_epi8(xmm6,  xmm7);
    xmm8  = _mm_unpacklo_epi8(xmm8,  xmm9);
    xmm10 = _mm_unpacklo_epi8(xmm10, xmm11);
    xmm12 = _mm_unpacklo_epi8(xmm12, xmm13);
    xmm14 = _mm_unpacklo_epi8(xmm14, xmm15);
      
    xmm1 = xmm0;
    xmm0 = _mm_unpacklo_epi16(xmm0, xmm2);
    xmm1 = _mm_unpackhi_epi16(xmm1, xmm2);
    xmm5 = xmm4;
    xmm4 = _mm_unpacklo_epi16(xmm4, xmm6);
    xmm5 = _mm_unpackhi_epi16(xmm5, xmm6);
    xmm9 = xmm8;
    xmm8 = _mm_unpacklo_epi16(xmm8, xmm10);
    xmm9 = _mm_unpackhi_epi16(xmm9, xmm10);
    xmm13 = xmm12;
    xmm12 = _mm_unpacklo_epi16(xmm12, xmm14);
    xmm13 = _mm_unpackhi_epi16(xmm13, xmm14);

    xmm2  = xmm0;
    xmm0  = _mm_unpacklo_epi32(xmm0, xmm4);
    xmm2  = _mm_unpackhi_epi32(xmm2, xmm4);
    xmm6  = xmm1;
    xmm1  = _mm_unpacklo_epi32(xmm1, xmm5);
    xmm6  = _mm_unpackhi_epi32(xmm6, xmm5);
    xmm10 = xmm8;
    xmm8  = _mm_unpacklo_epi32(xmm8, xmm12);
    xmm10 = _mm_unpackhi_epi32(xmm10, xmm12);
    xmm14 = xmm9;
    xmm9  = _mm_unpacklo_epi32(xmm9, xmm13);
    xmm14 = _mm_unpackhi_epi32(xmm14, xmm13);
      
    xmm3  = xmm0;
    xmm0  = _mm_unpacklo_epi64(xmm0, xmm8);
    xmm3  = _mm_unpackhi_epi64(xmm3, xmm8);
    xmm7  = xmm2;
    xmm2  = _mm_unpacklo_epi64(xmm2, xmm10);
    xmm7  = _mm_unpackhi_epi64(xmm7, xmm10);
    xmm11 = xmm1;
    xmm1  = _mm_unpacklo_epi64(xmm1, xmm9);
    xmm11 = _mm_unpackhi_epi64(xmm11, xmm9);
    xmm15 = xmm6;
    xmm6  = _mm_unpacklo_epi64(xmm6, xmm14);
    xmm15 = _mm_unpackhi_epi64(xmm15, xmm14);

    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1536+  0), xmm0);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1536+ 64), xmm3);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1536+128), xmm2);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1536+192), xmm7);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1536+256), xmm1);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1536+320), xmm11);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1536+384), xmm6);
    _mm_store_si128(reinterpret_cast<__m128i*>(dprofile+(16*j)+1536+448), xmm15);
  }

}

#endif

// One pass over the query for a block of database residues (the
// former donormal and domasked kernels, selected by Masking)
template <typename Masking>
inline auto align_cells7(__m128i & S,
                         __m128i * hep,
                         __m128i * const * qp,
                         __m128i const Q,
                         __m128i const R,
                         long ql,
                         __m128i const Z,
                         Masking const & masking) -> void
{
  auto score = apply_mask<Ops_7>(S, masking);  // mask
  auto H0 = Z;
  auto H1 = H0;
  auto H2 = H0;
  auto H3 = H0;
  auto F0 = H0;
  auto F1 = H0;
  auto F2 = H0;
  auto F3 = H0;
  __m128i N1;
  __m128i N2;
  __m128i N3;

  for (long qi = 0; qi < ql; ++qi)
  {
    __m128i const * const x = qp[qi];  // load x from qp[qi]
    auto const N0 = apply_mask<Ops_7>(hep[2 * qi], masking);  // load N0, mask
    auto E = apply_mask<Ops_7>(hep[(2 * qi) + 1], masking);  // load E, mask

    onestep<Ops_7>(H0, N1, F0, x[0], E, score, Q, R);
    onestep<Ops_7>(H1, N2, F1, x[1], E, score, Q, R);
    onestep<Ops_7>(H2, N3, F2, x[2], E, score, Q, R);
    onestep<Ops_7>(H3, hep[2 * qi], F3, x[3], E, score, Q, R);

    hep[(2 * qi) + 1] = E;  // save E
    H0 = N0;
    H1 = N1;
    H2 = N2;
    H3 = N3;
  }

  S = score;  // save S
}

void
#ifdef SWIPE_SSSE3
search7_ssse3
#else
search7
#endif
       (BYTE * * q_start,
	BYTE gap_open_penalty,
	BYTE gap_extend_penalty,
	BYTE * score_matrix,
	BYTE * dprofile,
	BYTE * hearray,
	struct db_thread_s * dbt,
	long sequences,
	long const * seqnos,
	long * scores,
	long qlen)
{
  __m128i S;
  __m128i Q;
  __m128i R;
  __m128i T;
  __m128i M;
  __m128i Z;
  __m128i T0;
  auto * const hep = reinterpret_cast<__m128i*>(hearray);
  __m128i ** const qp = reinterpret_cast<__m128i**>(q_start);
  std::array<BYTE const *, CHANNELS> d_begin;
  std::array<BYTE const *, CHANNELS> d_end;
  
  // the database residues of the channels, 16-byte aligned for the loads
  alignas(__m128i) std::array<BYTE, CDEPTH * sizeof(__m128i)> dseqalloc;
  
  auto * dseq = dseqalloc.data();
  BYTE const zero = 0;

  std::array<long, CHANNELS> seq_id;
  long next_id = 0;
  unsigned done = 0;
  
  memset(hearray, 0x80, static_cast<std::size_t>(qlen) * hearray_row_bytes);

  Z  = _mm_set1_epi8(byte_0x80);
  T0 = _mm_set_epi8(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
		    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, byte_0x80);
  Q  = _mm_set1_epi8(static_cast<char>(gap_open_penalty));
  R  = _mm_set1_epi8(static_cast<char>(gap_extend_penalty));

  S = Z;

  for (std::size_t c = 0; c < CHANNELS; c++)
  {
    d_begin[c] = &zero;
    d_end[c] = d_begin[c];
    seq_id[c] = -1;
  }

  int easy = 0;

  while(true)
  {
    if (easy != 0)
    {
      // fill all channels

      for (std::size_t c = 0; c < CHANNELS; c++)
      {
	for (std::size_t j = 0; j < CDEPTH; j++)
	{
	  if (d_begin[c] < d_end[c])
	  {
	    dseq[(CHANNELS*j)+c] = *(d_begin[c]++);
	  }
	  else
	  {
	    dseq[(CHANNELS * j) + c] = 0;
	  }
	}
	if (d_begin[c] == d_end[c])
	{
	  easy = 0;
	}
      }

#ifdef SWIPE_SSSE3
      dprofile_shuffle7(dprofile, score_matrix, dseq);
#else
      dprofile_fill7(dprofile, score_matrix, dseq);
#endif

      align_cells7(S, hep, qp, Q, R, qlen, Z, No_mask{});
    }
    else
    {
      // One or more sequences ended in the previous block 
      // We have to switch over to a new sequence

      easy = 1;

      M = _mm_setzero_si128();
      T = T0;
      for (std::size_t c = 0; c < CHANNELS; c++)
      {
	if (d_begin[c] < d_end[c])
	{
	  // this channel has more sequence

	  for (std::size_t j = 0; j < CDEPTH; j++)
	  {
	    if (d_begin[c] < d_end[c])
	    {
	      dseq[(CHANNELS*j)+c] = *(d_begin[c]++);
	    }
	    else
	    {
	      dseq[(CHANNELS * j) + c] = 0;
	    }
	  }
	  if (d_begin[c] == d_end[c])
	  {
	    easy = 0;
	  }
	}
	else
	{
	  // sequence in channel c ended
	  // change of sequence

	  M = _mm_xor_si128(M, T);

	  long const cand_id = seq_id[c];
		  
	  if (cand_id >= 0)
	  {
	    // save score
	    long const score = (reinterpret_cast<BYTE*>(&S))[c] - 0x80;
	    scores[cand_id] = score;
	    done++;
	  }

	  if (next_id < sequences)
	  {
	    // get next sequence
	    seq_id[c] = next_id;
	    long const seqnosf = seqnos[next_id];

	    long ntlen = 0;
	    long const strand = (seqnosf >> 2) & 1;
	    long const frame = seqnosf & 3;
	    long const seqno = seqnosf >> 3;

	    View<char> const sequence =
	      db_getsequence(dbt, seqno, {strand, frame}, &ntlen, c);
		      
	    // printf("Seqno: %ld Address: %p\n", seqno, address);
	    d_begin[c] = reinterpret_cast<BYTE const *>(sequence.begin());
	    d_end[c] = reinterpret_cast<BYTE const *>(sequence.end());
	    next_id++;
		      
	    // fill channel
	    for (std::size_t j = 0; j < CDEPTH; j++)
	    {
	      if (d_begin[c] < d_end[c])
	      {
		dseq[(CHANNELS*j)+c] = *(d_begin[c]++);
	      }
	      else
	      {
		dseq[(CHANNELS * j) + c] = 0;
	      }
	    }
	    if (d_begin[c] == d_end[c])
	    {
	      easy = 0;
	    }
	  }
	  else
	  {
	    // no more sequences, empty channel
	    seq_id[c] = -1;
	    d_begin[c] = &zero;
	    d_end[c] = d_begin[c];
	    for (std::size_t j = 0; j < CDEPTH; j++)
	    {
	      dseq[(CHANNELS * j) + c] = 0;
	    }
	  }


	}

	T = _mm_slli_si128(T, 1);
      }

      if (done == sequences)
      {
	break;
      }

#ifdef SWIPE_SSSE3
      dprofile_shuffle7(dprofile, score_matrix, dseq);
#else
      dprofile_fill7(dprofile, score_matrix, dseq);
#endif
	  
      align_cells7(S, hep, qp, Q, R, qlen, Z, Mask{M});
    }
  }
}
