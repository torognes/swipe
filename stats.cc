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
#include <cmath>
#include <cstddef>  // std::size_t

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
  auto const tables = blastn_tables(match_score, mismatch_score);
  auto const bv = tables.values;
  if (bv.empty())
  {
    return 0;
  }

  if ((gopen >= tables.gap_open_max) && (gextend >= tables.gap_extend_max))
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
  auto const mat = blast_matrix_values(matrix);
  if (mat.empty())
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
  auto const mat = blast_matrix_values(matrix);
  auto const prefs = blast_matrix_prefs(matrix);
  if (mat.empty())
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
