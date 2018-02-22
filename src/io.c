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
 * \file io.c
 * \brief Input and output handling
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

static inline void set_function_pointers(
    void
    )
{
  /* todo: this needs to be generalized for different monomial orders */
  matrix_row_initial_input_cmp  =
    matrix_row_initial_input_cmp_drl;
  monomial_cmp  = monomial_cmp_drl;
  hcm_cmp       = hcm_cmp_pivots_drl;
  select_spairs = select_spairs_by_minimal_degree;

  switch (laopt) {
    case 1:
      linear_algebra  = sparse_linear_algebra;
      break;
    case 42:
      linear_algebra  = probabilistic_sparse_linear_algebra;
      break;
    default:
      linear_algebra  = sparse_linear_algebra;
  }
  
  if (fc < pow(2, 16)) {
    reduce_dense_row_by_known_pivots = reduce_dense_row_by_known_pivots_16_bit;
  } else {
    /* TODO: 32 bit implementation */
    reduce_dense_row_by_known_pivots = reduce_dense_row_by_known_pivots_32_bit;
  }
}

static inline int32_t check_and_set_meta_data(
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const int32_t field_char,
    const int32_t nr_vars,
    const int32_t nr_gens,
    const int32_t ht_size,
    const int32_t nr_threads,
    const int32_t max_nr_pairs,
    const int32_t la_option
    )
{
  if (nr_gens <= 0
    || nr_vars <= 0
    || field_char <= 0
    || lens == NULL
    || cfs == NULL
    || exps == NULL) {
    return 1;
  }

  nvars = nr_vars;
  /* TODO: add prime check */
  fc    = field_char;
  /* set hash table size */
  htes  = ht_size;
  if (htes <= 0) {
    htes  = 12;
  }

  /* set number of threads */
  if (nr_threads <= 0) {
    nthrds  = 1;
  } else {
    nthrds  = nr_threads;
  }

  if (max_nr_pairs <= 0) {
    mnsel = 2147483647; /* 2^31-1 */
  } else {
    mnsel = max_nr_pairs;
  }

  /* set linear algebra option */
  if (la_option <= 0) {
    laopt = 1;
  } else {
    laopt = la_option;
  }

  set_function_pointers();

  return 0;
}

/* note that depending on the input data we set the corresponding
 * function pointers for monomial resp. spair comparisons, taking
 * spairs by a given minimal property for symbolic preprocessing, etc. */
static val_t **import_julia_data(
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const int32_t nr_gens
    )
{
  int32_t i, j;
  int32_t off = 0; /* offset in arrays */
  val_t **mat = malloc((unsigned long)nr_gens * sizeof(val_t *));
  
  for (i = 0; i < nr_gens; ++i) {
    /* each matrix row has the following structure:
     * [length | offset | eh1 | cf1 | eh2 | cf2 | .. | ehl | cfl]
     * where piv? is a label for being a known or an unknown pivot */
    mat[i]    = malloc((2*(unsigned long)lens[i]+2) * sizeof(val_t));
    mat[i][0] = 2*lens[i]+2;
    mat[i][1] = 2*(lens[i] % UNROLL)+2; /* offset for starting loop unrolling */
    for (j = off; j < off+lens[i]; ++j) {
      mat[i][2*(j+1-off)]   = insert_in_global_hash_table(exps+(nvars*j));
      mat[i][2*(j+1-off)+1] = cfs[j];
    }
    /* mark initial generators, they have to be added to the basis first */
    off +=  lens[i];
  }
  npivs = nrows = nrall = nr_gens;

  return mat;
}

static int64_t export_julia_data(
    int32_t **bp
    )
{
  int32_t i, j, k, ctr_lengths, ctr_elements;
  int32_t *basis  = *bp;

  int64_t len = 0; /* complete length of exported array */
  int64_t nb  = 0; /* # elemnts in basis */

  const int32_t lterm = 1 + nvars; /* length of a term */

  /* compute number of terms */
  for (i = 0; i < bload; ++i) {
    if ((long)bs[i] & bred) {
      continue;
    } else {
      len +=  (int64_t)((bs[i][0]-2)/2);
      nb++;
    }
  }

  /* compute the length considering the number of variables per exponent */
  len = len * (int64_t)lterm;
  /* add storage for length of each element */
  len = len + nb;
  /* add storage for length of complete array */
  len++;
  /* add storage for number of generators in basis */
  len++;

  basis  = (int32_t *)malloc((unsigned long)len * sizeof(int32_t));

  if (nb > (int64_t)(pow(2, 31))) {
    printf("basis too big\n");
    return 0;
  }

  ctr_lengths   = 1;
  ctr_elements  = (int32_t)nb + 1;

  basis[0]  = (int32_t)nb;
  /* basis[1]  = (int32_t)nb; */
  for (i = 0; i < bload; ++i) {
    if ((long)bs[i] & bred) {
      continue;
    } else {
      /* length of polynomial including this length entry itself */
      basis[ctr_lengths++]  = (bs[i][0]-2)/2 * lterm;
      for (j = 2; j < bs[i][0]; j += 2) {
        basis[ctr_elements++] = bs[i][j+1]; /* coefficient */
        for (k = 0; k < nvars; ++k) {
          basis[ctr_elements++] = (ev + bs[i][j])[k];
        }
      }
    }
  }
  *bp = basis;

  return len;
}
