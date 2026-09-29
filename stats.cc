/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2013 Torbjorn Rognes, University of Oslo,
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
#include <algorithm>  // std::max
#include <array>
#include <cmath>
#include <cstddef>  // std::size_t
#include <cstring>
#include <cstdio>

// the NCBI names used by blastkar_partial.cc
using array_of_8 = std::array<double, 8>;  // a row of statistical parameters
constexpr Int4 BLAST_MATRIX_NOMINAL = 0;
constexpr Int4 BLAST_MATRIX_BEST = 1;
constexpr Int4 INT2_MAX = 32767;

#include "blastkar_partial.cc"

auto stats_getparams_nt(long match_score,
			long mismatch_score, 
			long gopen,
			long gextend,
			double * lambda,
			double * K,
			double * H,
			double * alpha,
			double * beta) -> long
{
  View<array_of_8> bv;
  long gomax = 0;
  long gemax = 0;

  if      ((match_score == 1) && (mismatch_score == -5))
  {
    bv = make_view(blastn_values_1_5);
    gomax = 3;
    gemax = 3;
  }
  else if ((match_score == 1) && (mismatch_score == -4))
  {
    bv = make_view(blastn_values_1_4);
    gomax = 2;
    gemax = 2;
  }
  else if ((match_score == 2) && (mismatch_score == -7))
  {
    bv = make_view(blastn_values_2_7);
    gomax = 4;
    gemax = 4;
  }
  else if ((match_score == 1) && (mismatch_score == -3))
  {
    bv = make_view(blastn_values_1_3);
    gomax = 2;
    gemax = 2;
  }
  else if ((match_score == 2) && (mismatch_score == -5))
  {
    bv = make_view(blastn_values_2_5);
    gomax = 4;
    gemax = 4;
  }
  else if ((match_score == 1) && (mismatch_score == -2))
  {
    bv = make_view(blastn_values_1_2);
    gomax = 2;
    gemax = 2;
  }
  else if ((match_score == 2) && (mismatch_score == -3))
  {
    bv = make_view(blastn_values_2_3);
    gomax = 6;
    gemax = 4;
  }
  else if ((match_score == 3) && (mismatch_score == -4))
  {
    bv = make_view(blastn_values_3_4);
    gomax = 6;
    gemax = 3;
  }
  else if ((match_score == 4) && (mismatch_score == -5))
  {
    bv = make_view(blastn_values_4_5);
    gomax = 4;
    gemax = 2;
  }
  else if ((match_score == 1) && (mismatch_score == -1))
  {
    bv = make_view(blastn_values_1_1);
    gomax = 5;
    gemax = 5;
  }
  else if ((match_score == 3) && (mismatch_score == -2))
  {
    bv = make_view(blastn_values_3_2);
    gomax = 12;
    gemax = 8;
  }
  else if ((match_score == 5) && (mismatch_score == -4))
  {
    bv = make_view(blastn_values_5_4);
    gomax = 25;
    gemax = 10;
  }
  else
  {
    return 0;
  }

  if ((gopen >= gomax) && (gextend >= gemax))
  {
    gopen = 0;
    gextend = 0;
  }

  for(std::size_t i = 0; i < bv.size(); i++)
  {
    if ( (fabs(bv[i][0] - (static_cast<double>(gopen))) < 0.1) &&
	 (fabs(bv[i][1] - (static_cast<double>(gextend))) < 0.1) )
    {
      * lambda = bv[i][2];
      * K = bv[i][3];
      * H = bv[i][4];
      * alpha = bv[i][5];
      * beta = bv[i][6];
      return 1;
    }
  }

  return 0;
}

auto stats_getparams(char const * matrix,
		     long gopen,
		     long gextend,
		     double * lambda,
		     double * K,
		     double * H,
		     double * alpha,
		     double * beta) -> long
{
  View<array_of_8> mat;

  if (strcasecmp(matrix, "BLOSUM45") == 0)
  {
    mat = make_view(blosum45_values);
  }
  else if (strcasecmp(matrix, "BLOSUM50") == 0)
  {
    mat = make_view(blosum50_values);
  }
  else if (strcasecmp(matrix, "BLOSUM62") == 0)
  {
    mat = make_view(blosum62_values);
  }
  else if (strcasecmp(matrix, "BLOSUM80") == 0)
  {
    mat = make_view(blosum80_values);
  }
  else if (strcasecmp(matrix, "BLOSUM90") == 0)
  {
    mat = make_view(blosum90_values);
  }
  else if (strcasecmp(matrix, "PAM30") == 0)
  {
    mat = make_view(pam30_values);
  }
  else if (strcasecmp(matrix, "PAM70") == 0)
  {
    mat = make_view(pam70_values);
  }
  else if (strcasecmp(matrix, "PAM250") == 0)
  {
    mat = make_view(pam250_values);
  }
  else
  {
    return 0;
  }

  for (std::size_t i = 0; i < mat.size(); i++)
  {
    if ( (fabs(mat[i][0] - (static_cast<double>(gopen))) < 0.1) &&
	 (fabs(mat[i][1] - (static_cast<double>(gextend))) < 0.1) )
    {
      * lambda = mat[i][3];
      * K = mat[i][4];
      * H = mat[i][5];
      * alpha = mat[i][6];
      * beta = mat[i][7];

      //      printf("m=%s go=%ld ge=%ld: Chose index %ld: %-g %-g\n", matrix, gopen, gextend, i, mat[i][0], mat[i][1]);
            
      return 1;
    }
  }

  return 0;
}

auto stats_getprefs(char const * matrix,
		    long * gopen,
		    long * gextend) -> long
{
  View<array_of_8> mat;
  View<Int4> prefs;

  if (strcasecmp(matrix, "BLOSUM45") == 0)
  {
    mat = make_view(blosum45_values);
    prefs = make_view(blosum45_prefs);
  }
  else if (strcasecmp(matrix, "BLOSUM50") == 0)
  {
    mat = make_view(blosum50_values);
    prefs = make_view(blosum50_prefs);
  }
  else if (strcasecmp(matrix, "BLOSUM62") == 0)
  {
    mat = make_view(blosum62_values);
    prefs = make_view(blosum62_prefs);
  }
  else if (strcasecmp(matrix, "BLOSUM80") == 0)
  {
    mat = make_view(blosum80_values);
    prefs = make_view(blosum80_prefs);
  }
  else if (strcasecmp(matrix, "BLOSUM90") == 0)
  {
    mat = make_view(blosum90_values);
    prefs = make_view(blosum90_prefs);
  }
  else if (strcasecmp(matrix, "PAM30") == 0)
  {
    mat = make_view(pam30_values);
    prefs = make_view(pam30_prefs);
  }
  else if (strcasecmp(matrix, "PAM70") == 0)
  {
    mat = make_view(pam70_values);
    prefs = make_view(pam70_prefs);
  }
  else if (strcasecmp(matrix, "PAM250") == 0)
  {
    mat = make_view(pam250_values);
    prefs = make_view(pam250_prefs);
  }
  else
  {
    return 0;
  }

  for (std::size_t i = 0; i < mat.size(); i++)
  {
    if (prefs[i] != 0)
    {
      * gopen = static_cast<long>(mat[i][0]);
      * gextend = static_cast<long>(mat[i][1]);
      return 1;
    }
  }

  return 0;
}
