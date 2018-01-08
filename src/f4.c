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
 * 2.  an array of all exponents of all generators in the order:
 *     first all exponents of generator 1, then all of generator 2, ... */
int32_t *f4_julia(
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const int32_t field_char,
    const int32_t nr_vars,
    const int32_t nr_gens,
    const int32_t ht_size
    )
{
  int32_t i, hts_safe, round;
  val_t **mat;

  if (nr_gens == 0
    || nr_vars == 0
    || field_char == 0
    || lens == NULL
    || cfs == NULL
    || exps == NULL) {
    return NULL;
  }

  hts_safe  = ht_size;
  if (hts_safe == 0) {
    hts_safe  = 12;
  }
  
  /* initialize stuff */
  initialize_statistics();

  initialize_basis(nr_gens);
  initialize_pairset();
  initialize_global_hash_table(nr_vars, hts_safe, field_char);
  initialize_local_hash_table(hts_safe);

  mat = import_julia_data(lens, cfs, exps, nr_gens);

  /* for faster divisibility checks, needs to be done after we have
   * read some input data for applying heuristics */
  calculate_divmask();

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
    DEBUG(MATDBG, "rd %3d", round);

    /* preprocess data for next reduction round */
    mat = select_spairs();
    mat = symbolic_preprocessing(mat);
    /* exponent hashes mapped to column indices for linear algebra */
    mat = convert_hashes_to_columns(&mat);
    /* sort matrix rows by decreasing pivots */
    mat = sort_matrix_rows(mat);

    /* here starts the linear algebra part depending on
     * the chosen options */
    switch (laopt) {
      case 1:
        mat = sparse_linear_algebra(mat);
        break;
      default:
        mat = sparse_linear_algebra(mat);
    }

    convert_columns_to_hashes(mat);

  }

  /* free and clean up */
  free_local_hash_table();
  free_global_hash_table();
  free_pairset();
  /* note that all rows kept from mat during the overall computation are
   * basis elements and thus we do not need to free the rows itself, but
   * just the matrix structure */
  free(mat);
  free_basis();

  return mat[0];
}
