/* ===========================================================================
 *
 *                            PUBLIC DOMAIN NOTICE
 *               National Center for Biotechnology Information
 *
 *  This software/database is a "United States Government Work" under the
 *  terms of the United States Copyright Act.  It was written as part of
 *  the author's official duties as a United States Government employee and
 *  thus cannot be copyrighted.  This software/database is freely available
 *  to the public for use. The National Library of Medicine and the U.S.
 *  Government have not placed any restriction on its use or reproduction.
 *
 *  Although all reasonable efforts have been taken to ensure the accuracy
 *  and reliability of the software and data, the NLM and the U.S.
 *  Government do not and cannot warrant the performance or results that
 *  may be obtained by using this software or data. The NLM and the U.S.
 *  Government disclaim all warranties, express or implied, including
 *  warranties of performance, merchantability or fitness for any particular
 *  purpose.
 *
 *  Please cite the author in any work or product based on this material.
 *
 * ===========================================================================*/

#ifndef SWIPE_BLASTKAR_PARTIAL_H
#define SWIPE_BLASTKAR_PARTIAL_H

// a row of statistical parameters (NCBI's tables: gap open, gap
// extension, then the Karlin-Altschul parameters)
using array_of_8 = std::array<double, 8>;

// the tables of a score matrix and of a blastn score pair (swipe
// additions, in blastkar_partial.cc): empty views when unknown
struct BlastnTables
{
  View<array_of_8> values;
  long gap_open_max;  // from these gap costs on, the ungapped row
  long gap_extend_max;
};

auto blast_matrix_values(char const * matrix) -> View<array_of_8>;
auto blast_matrix_prefs(char const * matrix) -> View<Int4>;
auto blastn_tables(long match_score, long mismatch_score) -> BlastnTables;

// NCBI's length adjustment (the NCBI types are std::int32_t, std::int64_t
// and double, see swipe.h); documented in blastkar_partial.cc
auto BlastComputeLengthAdjustment(double K,
                                  double logK,
                                  double alpha_d_lambda,
                                  double beta,
                                  std::int32_t query_length,
                                  std::int64_t db_length,
                                  std::int32_t db_num_seqs,
                                  std::int32_t * length_adjustment) -> std::int32_t;

// swipe addition: the length adjustment itself (the NCBI function's
// status, 1 when the iteration did not converge, is not used)
auto length_adjustment(double K,
                       double logK,
                       double alpha_d_lambda,
                       double beta,
                       long query_length,
                       std::int64_t db_length,
                       std::int64_t db_sequences) -> std::int32_t;

#endif  // SWIPE_BLASTKAR_PARTIAL_H
