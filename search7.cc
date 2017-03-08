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

// #define DEBUG

#define CHANNELS 16
#define CDEPTH 4
#define MATRIXWIDTH 32

#ifdef SWIPE_SSSE3

inline void dprofile_shuffle7(BYTE * dprofile,
			      BYTE * score_matrix,
			      BYTE * dseq_byte)
{
#if MATRIXWIDTH > 16
  __m128i a, b, c, d, x, y, m0, m1, m2, m3, m4, m5, m6, m7;
  __m128i t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13;
  __m128i u0, u1, u2, u3, u4, u5,         u8, u9, u10, u11, u12, u13;
#else
  __m128i m0, m1, m2, m3, t0, t1, t2, t3, t4;
#endif

  __m128i * dseq = (__m128i*) dseq_byte;
  
  // 16 x 4 = 64 db symbols
  // ca 458 instructions

  // make masks

  /* Note: pshufb only on modern Intel cpus (SSSE3), not AMD */
  /* SSSE3: Supplemental SSE3 */

#if MATRIXWIDTH > 16
  x = _mm_set_epi8(0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10,
                   0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10, 0x10);

  y = _mm_set_epi8(0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80,
                   0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80);

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
  t6  = _mm_load_si128((__m128i*)(score_matrix)+2*j);   \
  t7  = _mm_load_si128((__m128i*)(score_matrix)+2*j+1); \
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
  _mm_store_si128((__m128i*)(dprofile)+4*j,   t12);	\
  _mm_store_si128((__m128i*)(dprofile)+4*j+1, t13);	\
  _mm_store_si128((__m128i*)(dprofile)+4*j+2, u12);	\
  _mm_store_si128((__m128i*)(dprofile)+4*j+3, u13)

#else

  m0 = _mm_load_si128(dseq);
  m1 = _mm_load_si128(dseq+1);
  m2 = _mm_load_si128(dseq+2);
  m3 = _mm_load_si128(dseq+3);

#define profline(j)					\
  t0 = _mm_load_si128((__m128i*)(score_matrix)+2*j);	\
  t1 = _mm_shuffle_epi8(t0, m0);			\
  t2 = _mm_shuffle_epi8(t0, m1);			\
  t3 = _mm_shuffle_epi8(t0, m2);			\
  t4 = _mm_shuffle_epi8(t0, m3);			\
  _mm_store_si128((__m128i*)(dprofile)+4*j+0, t1);	\
  _mm_store_si128((__m128i*)(dprofile)+4*j+1, t2);	\
  _mm_store_si128((__m128i*)(dprofile)+4*j+2, t3);	\
  _mm_store_si128((__m128i*)(dprofile)+4*j+3, t4)

#endif

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

#if MATRIXWIDTH > 16
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

#if MATRIXWIDTH > 28
  profline(28);
  profline(29);
  profline(30);
  profline(31);
#endif

#endif

  //  dprofile_dump7(dprofile);
}

#else

inline void dprofile_fill7(BYTE * dprofile,
			   BYTE * score_matrix,
			   BYTE * dseq)
{
  __m128i xmm0,  xmm1, xmm2,  xmm3,  xmm4,  xmm5,  xmm6,  xmm7;
  __m128i xmm8,  xmm9, xmm10, xmm11, xmm12, xmm13, xmm14, xmm15;
  
  // 4 x 16 db symbols
  // ca (60x2+68x2)x4 = 976 instructions

  for(int j=0; j<CDEPTH; j++)
  {
    unsigned d[CHANNELS];
    for(int i=0; i<CHANNELS; i++)
      d[i] = dseq[j*CHANNELS+i] << 5;
      
    xmm0  = _mm_loadl_epi64((__m128i*)(score_matrix + d[0] ));
    xmm2  = _mm_loadl_epi64((__m128i*)(score_matrix + d[2] ));
    xmm4  = _mm_loadl_epi64((__m128i*)(score_matrix + d[4] ));
    xmm6  = _mm_loadl_epi64((__m128i*)(score_matrix + d[6] ));
    xmm8  = _mm_loadl_epi64((__m128i*)(score_matrix + d[8] ));
    xmm10 = _mm_loadl_epi64((__m128i*)(score_matrix + d[10]));
    xmm12 = _mm_loadl_epi64((__m128i*)(score_matrix + d[12]));
    xmm14 = _mm_loadl_epi64((__m128i*)(score_matrix + d[14]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  *(__m128i*)(score_matrix + d[1] ));
    xmm2  = _mm_unpacklo_epi8(xmm2,  *(__m128i*)(score_matrix + d[3] ));
    xmm4  = _mm_unpacklo_epi8(xmm4,  *(__m128i*)(score_matrix + d[5] ));
    xmm6  = _mm_unpacklo_epi8(xmm6,  *(__m128i*)(score_matrix + d[7] ));
    xmm8  = _mm_unpacklo_epi8(xmm8,  *(__m128i*)(score_matrix + d[9] ));
    xmm10 = _mm_unpacklo_epi8(xmm10, *(__m128i*)(score_matrix + d[11]));
    xmm12 = _mm_unpacklo_epi8(xmm12, *(__m128i*)(score_matrix + d[13]));
    xmm14 = _mm_unpacklo_epi8(xmm14, *(__m128i*)(score_matrix + d[15]));
      
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

    _mm_store_si128((__m128i*)(dprofile+16*j+  0), xmm0);
    _mm_store_si128((__m128i*)(dprofile+16*j+ 64), xmm3);
    _mm_store_si128((__m128i*)(dprofile+16*j+128), xmm2);
    _mm_store_si128((__m128i*)(dprofile+16*j+192), xmm7);
    _mm_store_si128((__m128i*)(dprofile+16*j+256), xmm1);
    _mm_store_si128((__m128i*)(dprofile+16*j+320), xmm11);
    _mm_store_si128((__m128i*)(dprofile+16*j+384), xmm6);
    _mm_store_si128((__m128i*)(dprofile+16*j+448), xmm15);


    // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

    xmm0  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[0 ]));
    xmm1  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[1 ]));
    xmm2  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[2 ]));
    xmm3  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[3 ]));
    xmm4  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[4 ]));
    xmm5  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[5 ]));
    xmm6  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[6 ]));
    xmm7  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[7 ]));
    xmm8  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[8 ]));
    xmm9  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[9 ]));
    xmm10 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[10]));
    xmm11 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[11]));
    xmm12 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[12]));
    xmm13 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[13]));
    xmm14 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[14]));
    xmm15 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[15]));

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

    _mm_store_si128((__m128i*)(dprofile+16*j+512+  0), xmm0);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+ 64), xmm3);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+128), xmm2);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+192), xmm7);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+256), xmm1);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+320), xmm11);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+384), xmm6);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+448), xmm15);


    xmm0  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[0 ]));
    xmm2  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[2 ]));
    xmm4  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[4 ]));
    xmm6  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[6 ]));
    xmm8  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[8 ]));
    xmm10 = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[10]));
    xmm12 = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[12]));
    xmm14 = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[14]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  *(__m128i*)(score_matrix + 16 + d[1 ]));
    xmm2  = _mm_unpacklo_epi8(xmm2,  *(__m128i*)(score_matrix + 16 + d[3 ]));
    xmm4  = _mm_unpacklo_epi8(xmm4,  *(__m128i*)(score_matrix + 16 + d[5 ]));
    xmm6  = _mm_unpacklo_epi8(xmm6,  *(__m128i*)(score_matrix + 16 + d[7 ]));
    xmm8  = _mm_unpacklo_epi8(xmm8,  *(__m128i*)(score_matrix + 16 + d[9 ]));
    xmm10 = _mm_unpacklo_epi8(xmm10, *(__m128i*)(score_matrix + 16 + d[11 ]));
    xmm12 = _mm_unpacklo_epi8(xmm12, *(__m128i*)(score_matrix + 16 + d[13 ]));
    xmm14 = _mm_unpacklo_epi8(xmm14, *(__m128i*)(score_matrix + 16 + d[15 ]));
      
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

    _mm_store_si128((__m128i*)(dprofile+16*j+1024+  0), xmm0);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+ 64), xmm3);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+128), xmm2);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+192), xmm7);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+256), xmm1);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+320), xmm11);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+384), xmm6);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+448), xmm15);


    // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

    xmm0  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[0 ]));
    xmm1  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[1 ]));
    xmm2  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[2 ]));
    xmm3  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[3 ]));
    xmm4  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[4 ]));
    xmm5  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[5 ]));
    xmm6  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[6 ]));
    xmm7  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[7 ]));
    xmm8  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[8 ]));
    xmm9  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[9 ]));
    xmm10 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[10]));
    xmm11 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[11]));
    xmm12 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[12]));
    xmm13 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[13]));
    xmm14 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[14]));
    xmm15 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[15]));

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

    _mm_store_si128((__m128i*)(dprofile+16*j+1536+  0), xmm0);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+ 64), xmm3);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+128), xmm2);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+192), xmm7);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+256), xmm1);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+320), xmm11);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+384), xmm6);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+448), xmm15);
  }

  //  dprofile_dump7(dprofile);
}

void dprofile_dump7(BYTE * dprofile)
{
  const char * ss = sym_ncbi_aa;
  //  char * ss = sym_sound;

  printf("\ndprofile:\n");
  for(int k=0; k<4; k++)
  {
    printf("k=%d 0 1 2 3 4 5 6 7 8 9 a b c d e f\n", k);
    for(int i=0; i<32; i++)
    {
      printf("%c: ",ss[i]);
      for(int j=0; j<16; j++)
	printf("%2d", (char) dprofile[i*64+16*k+j]);
      printf("\n");
    }
  }
  printf("\n");
  exit(1);
}

int dumpcounter = 0;
char lines[4*16*1000];

void dseq_dump7(BYTE * dseq)
{
  const char * s = sym_ncbi_aa;

  if (dumpcounter < 21)
  {
    for(int i=0; i<CHANNELS; i++)
    {
      for(int j=0; j<CDEPTH; j++)
      {
	lines[4000*i+4*dumpcounter+j] = s[dseq[j*CHANNELS+i]];
      }
    }
    dumpcounter++;
  }
  else
  {
    for(int i=0; i<16; i++)
    {
      printf("%.1000s\n", lines+4000*i);
    }
    exit(1);
  }
}

#endif

// Register usage
// rdi:   hep
// rsi:   qp
// rdx:   Qm
// rcx:   Rm
// r8:    ql
// r9:    Sm/Mm

// rax:   x, temp
// r10:   ql2
// r11:   qi
// xmm0:  H0
// xmm1:  H1
// xmm2:  H2
// xmm3:  H3
// xmm4:  F0
// xmm5:  F1
// xmm6:  F2
// xmm7:  F3
// xmm8:  N0
// xmm9:  N1
// xmm10: N2
// xmm11: N3
// xmm12: E
// xmm13: S
// xmm14: Q 
// xmm15: R


#define INITIALIZE					\
  "        movq    %0, %%rax        \n"			\
  "        movdqa  (%%rax), %%xmm13 \n"			\
  "        movdqa  (%3), %%xmm14    \n"			\
  "        movdqa  (%4), %%xmm15    \n"			\
  "        movq    %6, %%rax        \n"			\
  "        movdqa  (%%rax), %%xmm0  \n"			\
  "        movdqa  %%xmm0, %%xmm1   \n"			\
  "        movdqa  %%xmm0, %%xmm2   \n"			\
  "        movdqa  %%xmm0, %%xmm3   \n"			\
  "        movdqa  %%xmm0, %%xmm4   \n"			\
  "        movdqa  %%xmm0, %%xmm5   \n"			\
  "        movdqa  %%xmm0, %%xmm6   \n"			\
  "        movdqa  %%xmm0, %%xmm7   \n"			\
  "        movq    %5, %%r12        \n"			\
  "        shlq    $3, %%r12        \n"			\
  "        movq    %%r12, %%r10     \n"			\
  "        andq    $-16, %%r10      \n"			\
  "        xorq    %%r11, %%r11     \n"

#define ONESTEP(H, N, F, V)                             \
  "        paddsb  " V "(%%rax), " H "\n"               \
  "        pmaxub  " F ", " H "       \n"               \
  "        pmaxub  %%xmm12, " H "     \n"               \
  "        pmaxub  " H ", %%xmm13     \n"               \
  "        psubsb  %%xmm15, " F "     \n"               \
  "        psubsb  %%xmm15, %%xmm12   \n"               \
  "        movdqa  " H ", " N "       \n"               \
  "        psubsb  %%xmm14, " H "     \n"               \
  "        pmaxub  " H ", %%xmm12     \n"               \
  "        pmaxub  " H ", " F "       \n"

inline void donormal7(__m128i * Sm,
		      __m128i * hep,
		      __m128i ** qp,
		      __m128i * Qm,
		      __m128i * Rm,
		      long ql,
		      __m128i * Zm)
{
#ifdef DEBUG
  printf("donormal\n");
  printf("Sm=%p\n", Sm);
  printf("hep=%p\n", hep);
  printf("qp=%p\n", qp);
  printf("Qm=%p\n", Qm);
  printf("Rm=%p\n", Rm);
  printf("qlen=%ld\n", ql);
  printf("Zm=%p\n", Zm);
#endif
  
  __asm__
    __volatile__
    (
     "## donormal7                             \n"
     INITIALIZE
     "        jmp     2f                      \n"
     
     "1:      movq    0(%2,%%r11,1), %%rax    \n" // load x from qp[qi]
     "        movdqa  0(%1,%%r11,4), %%xmm8   \n" // load N0
     "        movdqa  16(%1,%%r11,4), %%xmm12 \n" // load E
     
     ONESTEP("%%xmm0",  "%%xmm9",          "%%xmm4", "0" )
     ONESTEP("%%xmm1",  "%%xmm10",         "%%xmm5", "16")
     ONESTEP("%%xmm2",  "%%xmm11",         "%%xmm6", "32")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",   "%%xmm7", "48")
     
     "        movdqa  %%xmm12, 16(%1,%%r11,4) \n" // save E
     "        movq    8(%2,%%r11,1), %%rax    \n" // load x from qp[qi+1]
     "        movdqa  32(%1,%%r11,4), %%xmm0  \n" // load H0
     "        movdqa  48(%1,%%r11,4), %%xmm12 \n" // load E
     
     ONESTEP("%%xmm8",  "%%xmm1",           "%%xmm4", "0" )
     ONESTEP("%%xmm9",  "%%xmm2",           "%%xmm5", "16")
     ONESTEP("%%xmm10", "%%xmm3",           "%%xmm6", "32")
     ONESTEP("%%xmm11", "32(%1,%%r11,4)",   "%%xmm7", "48")
     
     "        movdqa  %%xmm12, 48(%1,%%r11,4) \n" // save E
     "        addq    $16, %%r11              \n" // qi++
     "2:      cmpq    %%r11, %%r10            \n" // qi = ql4 ?
     "        jne     1b                      \n" // loop
     
     "        cmpq    %%r11, %%r12            \n" 
     "        je      3f                      \n"
     "        movq    0(%2,%%r11,1), %%rax    \n" // load x from qp[qi]
     "        movdqa  16(%1,%%r11,4), %%xmm12 \n" // load E
     
     ONESTEP("%%xmm0",  "%%xmm9",          "%%xmm4", "0" )
     ONESTEP("%%xmm1",  "%%xmm10",         "%%xmm5", "16")
     ONESTEP("%%xmm2",  "%%xmm11",         "%%xmm6", "32")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",   "%%xmm7", "48")
     
     "        movdqa  %%xmm12, 16(%1,%%r11,4) \n" // save E
     "3:      movq    %0, %%rax               \n" // save S
     "        movdqa  %%xmm13, (%%rax)          "
     :
     : "m"(Sm), "r"(hep),"r"(qp), "r"(Qm), "r"(Rm), "r"(ql), "m"(Zm)
     : "xmm0",  "xmm1",  "xmm2",  "xmm3",
       "xmm4",  "xmm5",  "xmm6",  "xmm7",
       "xmm8",  "xmm9",  "xmm10", "xmm11", 
       "xmm12", "xmm13", "xmm14", "xmm15",
       "rax",   "r10",   "r11",   "r12",
       "cc"
     );
}

inline void domasked7(__m128i * Sm,
		      __m128i * hep,
		      __m128i ** qp,
		      __m128i * Qm, 
		      __m128i * Rm, 
		      long ql,      
		      __m128i * Zm,
		      __m128i * Mm)
{
  
#ifdef DEBUG
  printf("domasked\n");
  printf("Sm=%p\n", Sm);
  printf("hep=%p\n", hep);
  printf("qp=%p\n", qp);
  printf("Qm=%p\n", Qm);
  printf("Rm=%p\n", Rm);
  printf("qlen=%ld\n", ql);
  printf("Zm=%p\n", Zm);
  printf("Mm=%p\n", Mm);
#endif
  
  __asm__
    __volatile__
    (
     "## domasked7                             \n"
     INITIALIZE
     "        paddsb  (%7), %%xmm13            \n" // mask
     "        jmp     2f                       \n"
     
     "1:      movq    0(%2,%%r11,1), %%rax     \n" // load x from qp[qi]
     "        movdqa  0(%1,%%r11,4), %%xmm8    \n" // load N0
     "        paddsb  (%7), %%xmm8             \n" // mask
     "        movdqa  16(%1,%%r11,4), %%xmm12  \n" // load E
     "        paddsb  (%7), %%xmm12            \n" // mask
     
     ONESTEP("%%xmm0",  "%%xmm9",          "%%xmm4", "0" )
     ONESTEP("%%xmm1",  "%%xmm10",         "%%xmm5", "16")
     ONESTEP("%%xmm2",  "%%xmm11",         "%%xmm6", "32")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",   "%%xmm7", "48")
     
     "        movdqa  %%xmm12, 16(%1,%%r11,4)  \n" // save E
     "        movq    8(%2,%%r11,1), %%rax     \n" // load x from qp[qi+1]
     "        movdqa  32(%1,%%r11,4), %%xmm0   \n" // load H0
     "        paddsb  (%7), %%xmm0             \n" // mask
     "        movdqa  48(%1,%%r11,4), %%xmm12  \n" // load E
     "        paddsb  (%7), %%xmm12            \n" // mask
     
     ONESTEP("%%xmm8",  "%%xmm1",           "%%xmm4", "0" )
     ONESTEP("%%xmm9",  "%%xmm2",           "%%xmm5", "16")
     ONESTEP("%%xmm10", "%%xmm3",           "%%xmm6", "32")
     ONESTEP("%%xmm11", "32(%1,%%r11,4)",   "%%xmm7", "48")
     
     "        movdqa  %%xmm12, 48(%1,%%r11,4)  \n" // save E
     "        addq    $16, %%r11               \n" // qi++
     "2:      cmpq    %%r11, %%r10             \n" // qi = ql4 ?
     "        jne     1b                       \n" // loop
     
     "        cmpq    %%r11, %%r12             \n" 
     "        je      3f                       \n"
     "        movq    0(%2,%%r11,1), %%rax     \n" // load x from qp[qi]
     "        movdqa  16(%1,%%r11,4), %%xmm12  \n" // load E
     "        paddsb  (%7), %%xmm12            \n" // mask
     
     ONESTEP("%%xmm0",  "%%xmm9",          "%%xmm4", "0" )
     ONESTEP("%%xmm1",  "%%xmm10",         "%%xmm5", "16")
     ONESTEP("%%xmm2",  "%%xmm11",         "%%xmm6", "32")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",   "%%xmm7", "48")
     
     "        movdqa  %%xmm12, 16(%1,%%r11,4)  \n" // save E
     "3:      movq    %0, %%rax                \n" // save S
     "        movdqa  %%xmm13, (%%rax)           "
     : 
     : "m"(Sm), "r"(hep),"r"(qp), "r"(Qm), "r"(Rm), "r"(ql), "m"(Zm),
       "r"(Mm)
     : "xmm0",  "xmm1",  "xmm2",  "xmm3",
       "xmm4",  "xmm5",  "xmm6",  "xmm7",
       "xmm8",  "xmm9",  "xmm10", "xmm11", 
       "xmm12", "xmm13", "xmm14", "xmm15",
       "rax",   "r10",   "r11",   "r12",
       "cc"
     );
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
	long * seqnos,
	long * scores,
	long qlen)
{
  __m128i S, Q, R, T, M, Z, T0;
  __m128i *hep, **qp;
  BYTE * d_begin[CHANNELS];
  BYTE * d_end[CHANNELS];
  
  __m128i dseqalloc[CDEPTH];
  
  BYTE * dseq = (BYTE*) & dseqalloc;
  BYTE zero;

  long seq_id[CHANNELS];
  long next_id = 0;
  unsigned done;
  
  memset(hearray, 0x80, qlen * 32);

  Z  = _mm_set_epi8(0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 
		    0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80, 0x80);
  T0 = _mm_set_epi8(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
		    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x80);
  Q  = _mm_set_epi8(gap_open_penalty, gap_open_penalty,
		    gap_open_penalty, gap_open_penalty,
		    gap_open_penalty, gap_open_penalty,
		    gap_open_penalty, gap_open_penalty,
		    gap_open_penalty, gap_open_penalty,
		    gap_open_penalty, gap_open_penalty,
		    gap_open_penalty, gap_open_penalty,
		    gap_open_penalty, gap_open_penalty);
  R  = _mm_set_epi8(gap_extend_penalty, gap_extend_penalty,
		    gap_extend_penalty, gap_extend_penalty,
		    gap_extend_penalty, gap_extend_penalty,
		    gap_extend_penalty, gap_extend_penalty,
		    gap_extend_penalty, gap_extend_penalty,
		    gap_extend_penalty, gap_extend_penalty,
		    gap_extend_penalty, gap_extend_penalty,
		    gap_extend_penalty, gap_extend_penalty);
  zero = 0;
  done = 0;

  S = Z;

  hep = (__m128i*) hearray;
  qp = (__m128i**) q_start;

#ifdef DEBUG
  //  printf("Searching %ld sequences...\n", sequences);
#endif

  for (int c=0; c<CHANNELS; c++)
  {
    d_begin[c] = &zero;
    d_end[c] = d_begin[c];
    seq_id[c] = -1;
  }

  int easy = 0;

  while(1)
  {
    if (easy)
    {
      // fill all channels

      for(int c=0; c<CHANNELS; c++)
      {
	for(int j=0; j<CDEPTH; j++)
	{
	  if (d_begin[c] < d_end[c])
	    dseq[CHANNELS*j+c] = *(d_begin[c]++);
	  else
	    dseq[CHANNELS*j+c] = 0;
	}
	if (d_begin[c] == d_end[c])
	  easy = 0;
      }

#ifdef SWIPE_SSSE3
      dprofile_shuffle7(dprofile, score_matrix, dseq);
#else
      dprofile_fill7(dprofile, score_matrix, dseq);
#endif

      donormal7(&S, hep, qp, &Q, &R, qlen, &Z);
    }
    else
    {
      // One or more sequences ended in the previous block 
      // We have to switch over to a new sequence

      easy = 1;

      M = _mm_setzero_si128();
      T = T0;
      for (int c=0; c<CHANNELS; c++)
      {
	if (d_begin[c] < d_end[c])
	{
	  // this channel has more sequence

	  for(int j=0; j<CDEPTH; j++)
	  {
	    if (d_begin[c] < d_end[c])
	      dseq[CHANNELS*j+c] = *(d_begin[c]++);
	    else
	      dseq[CHANNELS*j+c] = 0;
	  }
	  if (d_begin[c] == d_end[c])
	    easy = 0;
	}
	else
	{
	  // sequence in channel c ended
	  // change of sequence

	  M = _mm_xor_si128(M, T);

	  long cand_id = seq_id[c];
		  
	  if (cand_id >= 0)
	  {
	    // save score
	    long score = ((BYTE*)&S)[c] - 0x80;
	    scores[cand_id] = score;
	    done++;
	  }

	  if (next_id < sequences)
	  {
	    // get next sequence
	    seq_id[c] = next_id;
	    long seqnosf = seqnos[next_id];

	    char* address;
	    long length, ntlen;
	    long strand = (seqnosf >> 2) & 1;
	    long frame = seqnosf & 3;
	    long seqno = seqnosf >> 3;

	    db_getsequence(dbt, seqno, strand, frame, 
			   & address, & length, &ntlen, c);
		      
	    // printf("Seqno: %ld Address: %p\n", seqno, address);
	    d_begin[c] = (unsigned char*) address;
	    d_end[c] = (unsigned char*) address + length - 1;
	    next_id++;
		      
	    // fill channel
	    for(int j=0; j<CDEPTH; j++)
	    {
	      if (d_begin[c] < d_end[c])
		dseq[CHANNELS*j+c] = *(d_begin[c]++);
	      else
		dseq[CHANNELS*j+c] = 0;
	    }
	    if (d_begin[c] == d_end[c])
	      easy = 0;
	  }
	  else
	  {
	    // no more sequences, empty channel
	    seq_id[c] = -1;
	    d_begin[c] = &zero;
	    d_end[c] = d_begin[c];
	    for (int j=0; j<CDEPTH; j++)
	      dseq[CHANNELS*j+c] = 0;
	  }


	}

	T = _mm_slli_si128(T, 1);
      }

      if (done == sequences)
	break;
	  
#ifdef SWIPE_SSSE3
      dprofile_shuffle7(dprofile, score_matrix, dseq);
#else
      dprofile_fill7(dprofile, score_matrix, dseq);
#endif
	  
      domasked7(&S, hep, qp, &Q, &R, qlen, &Z, &M);
    }
  }
}
