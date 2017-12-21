/* gb: Gröbner Basis
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
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
    int32_t *lens,
    int32_t *cfs,
    int32_t *exps,
    int32_t nr_vars,
    int32_t nr_gens,
    int32_t ht_size
    )
{
  int i, j;
  val_t **mat;
  
  /* initialize stuff */
  initialize_basis(nr_gens);
  initialize_global_hash_table(nr_vars, ht_size);
  initialize_local_hash_table(ht_size);

  mat = import_julia_data(lens, cfs, exps, nr_gens);

  for (i = 0; i < nr_gens; ++i) {
    for (j = 1; j < mat[i][0]; ++j) {
      printf("%d ", mat[i][j]);
    }
    printf("\n");
  }

  /* free and clean up */
  free_local_hash_table();
  free_global_hash_table();
  free_basis();

  return lens;
}
