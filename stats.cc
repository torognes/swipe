/*
    SWIPE
    Smith-Waterman database searches with Inter-sequence Parallel Execution

    Copyright (C) 2008-2012 Torbjorn Rognes, University of Oslo,
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

typedef int Int4;
typedef long Int8;
typedef double Nlm_FloatHi;
typedef int Boolean;
typedef double array_of_8[8];
#define FALSE 0
#define TRUE 1
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define BLAST_MATRIX_NOMINAL 0
#define BLAST_MATRIX_BEST 1
#define INT2_MAX 32767

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "blastkar_partial.c"

long stats_getparams_nt(long matchscore,
			long mismatchscore, 
			long gopen,
			long gextend,
			double * lambda,
			double * K,
			double * H,
			double * alpha,
			double * beta)
{
  const array_of_8 * bv;
  long bm;
  long gomax;
  long gemax;

  if      ((matchscore == 1) && (mismatchscore == -5))
  {
    bv = blastn_values_1_5;
    bm = sizeof(blastn_values_1_5) / sizeof(array_of_8);
    gomax = 3;
    gemax = 3;
  }
  else if ((matchscore == 1) && (mismatchscore == -4))
  {
    bv = blastn_values_1_4;
    bm = sizeof(blastn_values_1_4) / sizeof(array_of_8);
    gomax = 2;
    gemax = 2;
  }
  else if ((matchscore == 2) && (mismatchscore == -7))
  {
    bv = blastn_values_2_7;
    bm = sizeof(blastn_values_2_7) / sizeof(array_of_8);
    gomax = 4;
    gemax = 4;
  }
  else if ((matchscore == 1) && (mismatchscore == -3))
  {
    bv = blastn_values_1_3;
    bm = sizeof(blastn_values_1_3) / sizeof(array_of_8);
    gomax = 2;
    gemax = 2;
  }
  else if ((matchscore == 2) && (mismatchscore == -5))
  {
    bv = blastn_values_2_5;
    bm = sizeof(blastn_values_2_5) / sizeof(array_of_8);
    gomax = 4;
    gemax = 4;
  }
  else if ((matchscore == 1) && (mismatchscore == -2))
  {
    bv = blastn_values_1_2;
    bm = sizeof(blastn_values_1_2) / sizeof(array_of_8);
    gomax = 2;
    gemax = 2;
  }
  else if ((matchscore == 2) && (mismatchscore == -3))
  {
    bv = blastn_values_2_3;
    bm = sizeof(blastn_values_2_3) / sizeof(array_of_8);
    gomax = 6;
    gemax = 4;
  }
  else if ((matchscore == 3) && (mismatchscore == -4))
  {
    bv = blastn_values_3_4;
    bm = sizeof(blastn_values_3_4) / sizeof(array_of_8);
    gomax = 6;
    gemax = 3;
  }
  else if ((matchscore == 4) && (mismatchscore == -5))
  {
    bv = blastn_values_4_5;
    bm = sizeof(blastn_values_4_5) / sizeof(array_of_8);
    gomax = 4;
    gemax = 2;
  }
  else if ((matchscore == 1) && (mismatchscore == -1))
  {
    bv = blastn_values_1_1;
    bm = sizeof(blastn_values_1_1) / sizeof(array_of_8);
    gomax = 5;
    gemax = 5;
  }
  else if ((matchscore == 3) && (mismatchscore == -2))
  {
    bv = blastn_values_3_2;
    bm = sizeof(blastn_values_3_2) / sizeof(array_of_8);
    gomax = 12;
    gemax = 8;
  }
  else if ((matchscore == 5) && (mismatchscore == -4))
  {
    bv = blastn_values_5_4;
    bm = sizeof(blastn_values_5_4) / sizeof(array_of_8);
    gomax = 25;
    gemax = 10;
  }
  else
    return 0;

  if ((gopen >= gomax) && (gextend >= gemax))
  {
    gopen = 0;
    gextend = 0;
  }

  for(long i = 0; i < bm; i++)
  {
    if ( (fabs(bv[i][0] - ((double) gopen)) < 0.1) &&
	 (fabs(bv[i][1] - ((double) gextend)) < 0.1) )
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

long stats_getparams(char * matrix,
		     long gopen,
		     long gextend,
		     double * lambda,
		     double * K,
		     double * H,
		     double * alpha,
		     double * beta)
{
  double (*mat)[8]; 
  long val;

  if (strcasecmp(matrix, "BLOSUM45") == 0)
  {
    mat = blosum45_values;
    val = BLOSUM45_VALUES_MAX;
  }
  else if (strcasecmp(matrix, "BLOSUM50") == 0)
  {
    mat = blosum50_values;
    val = BLOSUM50_VALUES_MAX;
  }
  else if (strcasecmp(matrix, "BLOSUM62") == 0)
  {
    mat = blosum62_values;
    val = BLOSUM62_VALUES_MAX;
  }
  else if (strcasecmp(matrix, "BLOSUM62_20") == 0)
  {
    mat = blosum62_20_values;
    val = BLOSUM62_20_VALUES_MAX;
  }
  else if (strcasecmp(matrix, "BLOSUM80") == 0)
  {
    mat = blosum80_values;
    val = BLOSUM80_VALUES_MAX;
  }
  else if (strcasecmp(matrix, "BLOSUM90") == 0)
  {
    mat = blosum90_values;
    val = BLOSUM90_VALUES_MAX;
  }
  else if (strcasecmp(matrix, "PAM30") == 0)
  {
    mat = pam30_values;
    val = PAM30_VALUES_MAX;
  }
  else if (strcasecmp(matrix, "PAM70") == 0)
  {
    mat = pam70_values;
    val = PAM70_VALUES_MAX;
  }
  else if (strcasecmp(matrix, "PAM250") == 0)
  {
    mat = pam250_values;
    val = PAM250_VALUES_MAX;
  }
  else
    return 0;

  for (long i=0; i<val; i++)
  {
    if ( (fabs(mat[i][0] - ((double) gopen)) < 0.1) &&
	 (fabs(mat[i][1] - ((double) gextend)) < 0.1) )
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

long stats_getprefs(char * matrix,
		    long * gopen,
		    long * gextend)
{
  double (*mat)[8]; 
  long val;
  Int4 *prefs;

  if (strcasecmp(matrix, "BLOSUM45") == 0)
  {
    mat = blosum45_values;
    val = BLOSUM45_VALUES_MAX;
    prefs = blosum45_prefs;
  }
  else if (strcasecmp(matrix, "BLOSUM50") == 0)
  {
    mat = blosum50_values;
    val = BLOSUM50_VALUES_MAX;
    prefs = blosum50_prefs;
  }
  else if (strcasecmp(matrix, "BLOSUM62") == 0)
  {
    mat = blosum62_values;
    val = BLOSUM62_VALUES_MAX;
    prefs = blosum62_prefs;
  }
  else if (strcasecmp(matrix, "BLOSUM62_20") == 0)
  {
    mat = blosum62_20_values;
    val = BLOSUM62_20_VALUES_MAX;
    prefs = blosum62_20_prefs;
  }
  else if (strcasecmp(matrix, "BLOSUM80") == 0)
  {
    mat = blosum80_values;
    val = BLOSUM80_VALUES_MAX;
    prefs = blosum80_prefs;
  }
  else if (strcasecmp(matrix, "BLOSUM90") == 0)
  {
    mat = blosum90_values;
    val = BLOSUM90_VALUES_MAX;
    prefs = blosum90_prefs;
  }
  else if (strcasecmp(matrix, "PAM30") == 0)
  {
    mat = pam30_values;
    val = PAM30_VALUES_MAX;
    prefs = pam30_prefs;
  }
  else if (strcasecmp(matrix, "PAM70") == 0)
  {
    mat = pam70_values;
    val = PAM70_VALUES_MAX;
    prefs = pam70_prefs;
  }
  else if (strcasecmp(matrix, "PAM250") == 0)
  {
    mat = pam250_values;
    val = PAM250_VALUES_MAX;
    prefs = pam250_prefs;
  }
  else
    return 0;

  for (long i=0; i<val; i++)
  {
    if (prefs[i])
    {
      * gopen = (long) mat[i][0];
      * gextend = (long) mat[i][1];
      return 1;
    }
  }

  return 0;
}
