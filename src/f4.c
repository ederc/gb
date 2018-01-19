/* gb: Gröbner Basis
 * Copyright (C) 2018 Christian Eder <ederc@mathematik.uni-kl.de>
 * This file is part of gb.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file f4.c
 * \brief Implementation of the F4 Algorithm
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

/* we get from julia the generators as three arrays:
 * 1.  an array of the lengths of each generator
 * 2.  an array of all coefficients of all generators in the order:
 *     first all coefficients of generator 1, then all of generator 2, ...
 * 3.  an array of all exponents of all generators in the order:
 *     first all exponents of generator 1, then all of generator 2, ...
 *     MORE TODO */
int32_t *f4_julia(
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const int32_t field_char,
    const int32_t nr_vars,
    const int32_t nr_gens,
    const int32_t ht_size,
    const int32_t nr_threads,
    const int32_t la_option
    )
{
  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  int32_t i, round;
  val_t **mat;
  len_t *hcm; /* hash-column-map */

  /* checks and set all meta data. if a nonzero value is returned then
   * some of the input data is corrupted. */
  if (check_and_set_meta_data(lens, cfs, exps, field_char, nr_vars,
      nr_gens, ht_size, nr_threads, la_option)) {
    return NULL;
  }

  /* initialize stuff */
  initialize_statistics();
  GB_DEBUG(GBDBG, "-------------------------------------------------\n");
  GB_DEBUG(GBDBG, "#variables             %15d\n", nvars);
  GB_DEBUG(GBDBG, "#equations             %15d\n", nr_gens);
  GB_DEBUG(GBDBG, "field characteristic   %15d\n", fc);
  GB_DEBUG(GBDBG, "linear algebra option  %15d\n", laopt);
  GB_DEBUG(GBDBG, "intial hash table size %15d (2^%d)\n",
      (int32_t)pow(2,htes), htes);
  GB_DEBUG(GBDBG, "#threads               %15d\n", nthrds);
  GB_DEBUG(GBDBG, "-------------------------------------------------\n");

  initialize_basis(nr_gens);
  initialize_pairset();
  initialize_global_hash_table();
  initialize_local_hash_table();

  mat = import_julia_data(lens, cfs, exps, nr_gens);

  /* for faster divisibility checks, needs to be done after we have
   * read some input data for applying heuristics */
  calculate_divmask();

  /* sort initial elements, smallest lead term first */
  qsort(mat, (unsigned long)nrows, sizeof(val_t *),
      matrix_row_initial_input_cmp);
  /* normalize input generators */
  for (i = 0; i < nrows; ++i) {
    normalize_matrix_row(mat[i]);
  }

  /* move input generators to basis and generate first spairs */
  update_basis(mat);

  free(mat);
  mat = NULL;

  /* let's start the f4 rounds,  we are done when no more spairs
   * are left in the pairset */
  for (round = 1; pload > 0; ++round) {
    GB_DEBUG(GBDBG, "%3d", round);

    /* preprocess data for next reduction round */
    mat = select_spairs();
    mat = symbolic_preprocessing(mat);
    /* exponent hashes mapped to column indices for linear algebra */
    hcm = convert_hashes_to_columns(mat);
    /* sort matrix rows by decreasing pivots */
    mat = sort_matrix_rows(mat);

    /* here starts the linear algebra part depending on
     * the chosen options */
    /* switch (laopt) {
     *   case 1:
     *     mat = sparse_linear_algebra(mat);
     *     break;
     *   case 42:
     *     mat = probabilistic_sparse_linear_algebra(mat);
     *     break;
     *   default:
     *     mat = sparse_linear_algebra(mat);
     * } */
    mat = linear_algebra(mat);

    mat = convert_columns_to_hashes(mat, hcm);

    free(hcm);
    hcm = NULL;

    update_basis(mat);

    free(mat);
    mat = NULL;
    GB_DEBUG(GBDBG, "\n");
  }


  int32_t *basis  = export_julia_data();

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  GB_DEBUG(GBDBG, "-------------------------------------------------\n");
  GB_DEBUG(GBDBG, "overall                %15.3f sec\n", rt1-rt0);
  GB_DEBUG(GBDBG, "symbol                 %15.3f sec\n", symbol_rtime);
  GB_DEBUG(GBDBG, "update                 %15.3f sec\n", update_rtime);
  GB_DEBUG(GBDBG, "convert                %15.3f sec\n", convert_rtime);
  GB_DEBUG(GBDBG, "la                     %15.3f sec\n", la_rtime);
  GB_DEBUG(GBDBG, "-------------------------------------------------\n");
  GB_DEBUG(GBDBG, "size of basis          %15d\n", basis[1]);
  GB_DEBUG(GBDBG, "#pairs reduced         %15ld\n", num_pairsred);
  GB_DEBUG(GBDBG, "#GM criterion          %15ld\n", num_gb_crit);
  GB_DEBUG(GBDBG, "#redundant             %15ld\n", num_redundant);
  GB_DEBUG(GBDBG, "#rows reduced          %15ld\n", num_rowsred);
  GB_DEBUG(GBDBG, "#zero reductions       %15ld\n", num_zerored);
  GB_DEBUG(GBDBG, "global hash table size %15d (2^%d)\n",
      msize, (int32_t)(log(msize)/log(2)));
  GB_DEBUG(GBDBG, "local hash table size  %15d (2^%d)\n",
      mlsize, (int32_t)(log(mlsize)/log(2)));
  if (laopt == 42) {
    GB_DEBUG(GBDBG, "#probable blocks fail  %15ld\n", num_prob_fail);
  }
  GB_DEBUG(GBDBG, "-------------------------------------------------\n");

  /* free and clean up */
  free_local_hash_table();
  free_global_hash_table();
  free_pairset();
  /* note that all rows kept from mat during the overall computation are
   * basis elements and thus we do not need to free the rows itself, but
   * just the matrix structure */
  free(mat);
  free_basis();

  return basis;
}
