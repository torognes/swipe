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

// anonymous namespace: limit visibility and usage to this translation unit
namespace {

// the row of a table for these gap penalties (the tables store the
// integer penalties as doubles); the first two columns are the same in
// both kinds of tables
auto has_penalties(array_of_8 const & row, long const gopen, long const gextend) -> bool
{
  constexpr double tolerance = 0.1;
  return (fabs(value_of(row, MatrixColumn::gap_open) - static_cast<double>(gopen)) < tolerance) and
    (fabs(value_of(row, MatrixColumn::gap_extend) - static_cast<double>(gextend)) < tolerance);
}

}  // anonymous namespace

auto stats_getparams_nt(BlastnScores const scores, GapPenalties const gaps) -> StatisticsLookup
{
  StatisticsLookup const not_found {false, {0, 0, 0, 0, 0}};
  long gopen = gaps.open;
  long gextend = gaps.extend;
  auto const tables = blastn_tables(scores.match, scores.mismatch);
  auto const bv = tables.values;
  if (bv.empty())
  {
    return not_found;
  }

  if ((gopen >= tables.gap_open_max) && (gextend >= tables.gap_extend_max))
  {
    gopen = 0;
    gextend = 0;
  }

  for(std::size_t i = 0; i < bv.size(); i++)
  {
    if (has_penalties(bv[i], gopen, gextend))
    {
      return {true, {value_of(bv[i], BlastnColumn::lambda),
                     value_of(bv[i], BlastnColumn::K),
                     value_of(bv[i], BlastnColumn::H),
                     value_of(bv[i], BlastnColumn::alpha),
                     value_of(bv[i], BlastnColumn::beta)}};
    }
  }

  return not_found;
}

auto stats_getparams(char const * matrix, GapPenalties const gaps) -> StatisticsLookup
{
  StatisticsLookup const not_found {false, {0, 0, 0, 0, 0}};
  long const gopen = gaps.open;
  long const gextend = gaps.extend;
  auto const mat = blast_matrix_values(matrix);
  if (mat.empty())
  {
    return not_found;
  }

  for (std::size_t i = 0; i < mat.size(); i++)
  {
    if (has_penalties(mat[i], gopen, gextend))
    {
      //      printf("m=%s go=%ld ge=%ld: Chose index %ld: %-g %-g\n", matrix, gopen, gextend, i, mat[i][0], mat[i][1]);
            
      return {true, {value_of(mat[i], MatrixColumn::lambda),
                     value_of(mat[i], MatrixColumn::K),
                     value_of(mat[i], MatrixColumn::H),
                     value_of(mat[i], MatrixColumn::alpha),
                     value_of(mat[i], MatrixColumn::beta)}};
    }
  }

  return not_found;
}

auto stats_getprefs(char const * matrix) -> DefaultGaps
{
  DefaultGaps const not_found {false, {0, 0}};
  auto const mat = blast_matrix_values(matrix);
  auto const prefs = blast_matrix_prefs(matrix);
  if (mat.empty())
  {
    return not_found;
  }

  for (std::size_t i = 0; i < mat.size(); i++)
  {
    if (prefs[i] != 0)
    {
      return {true, {static_cast<long>(value_of(mat[i], MatrixColumn::gap_open)),
                     static_cast<long>(value_of(mat[i], MatrixColumn::gap_extend))}};
    }
  }

  return not_found;
}
