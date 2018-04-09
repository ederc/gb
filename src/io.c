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
    const md_t *md
    )
{
  /* todo: this needs to be generalized for different monomial orders */
  switch (md->mo) {
    case 0:
      matrix_row_initial_input_cmp  =
        matrix_row_initial_input_cmp_drl;
      monomial_cmp  = monomial_cmp_drl;
      spair_cmp     = spair_cmp_drl;
      hcm_cmp       = hcm_cmp_pivots_drl;
      break;
    case 1:
      matrix_row_initial_input_cmp  =
        matrix_row_initial_input_cmp_lex;
      monomial_cmp  = monomial_cmp_lex;
      spair_cmp     = spair_cmp_deglex;
      hcm_cmp       = hcm_cmp_pivots_lex;
      break;
    default:
      matrix_row_initial_input_cmp  =
        matrix_row_initial_input_cmp_drl;
      monomial_cmp  = monomial_cmp_drl;
      spair_cmp     = spair_cmp_drl;
      hcm_cmp       = hcm_cmp_pivots_drl;
  }

  switch (md->la) {
    case 1:
      linear_algebra  = sparse_linear_algebra;
      break;
    case 42:
      linear_algebra  = probabilistic_sparse_linear_algebra;
      break;
    default:
      linear_algebra  = sparse_linear_algebra;
  }
  
  if (md->fc != 0) {
    if (md->fc < pow(2, 16)) {
      reduce_dense_row_by_known_pivots = reduce_dense_row_by_known_pivots_16_bit;
      import_julia_data = import_julia_data_16;
      export_julia_data = export_julia_data_16;
    } else {
      /* TODO: 32 bit implementation */
      reduce_dense_row_by_known_pivots = reduce_dense_row_by_known_pivots_32_bit;
      import_julia_data = import_julia_data_32;
      export_julia_data = export_julia_data_32;
    }
  } else {
    printf("no implementation for rationals yet!\n");
  }
}

static inline int32_t check_and_set_meta_data(
    md_t *md,
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const int32_t field_char,
    const int32_t mon_order,
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
  md->nv  = nr_vars;
  md->fc  = fc;
  /* note: prime check should be done in julia */
  fc    = field_char;
  /* monomial order */
  if (mon_order != 0 && mon_order != 1) {
    mo  = 0;
  } else {
    mo  = mon_order;
  }
  md->mo  = mo;
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
  md->nt  = nthrds;

  if (max_nr_pairs <= 0) {
    mnsel = 2147483647; /* 2^31-1 */
  } else {
    mnsel = max_nr_pairs;
  }
  md->mp  = mnsel;

  /* set linear algebra option */
  if (la_option <= 0) {
    laopt = 1;
  } else {
    laopt = la_option;
  }
  md->la  = laopt;

  /* reset statistics data */
  md->density           = 0;
  md->select_ctime      = 0;
  md->select_rtime      = 0;
  md->symbol_ctime      = 0;
  md->symbol_rtime      = 0;
  md->update_ctime      = 0;
  md->update_rtime      = 0;
  md->convert_ctime     = 0;
  md->convert_rtime     = 0;
  md->reduce_ctime      = 0;
  md->reduce_rtime      = 0;
  md->la_ctime          = 0;
  md->la_rtime          = 0;
  md->psort_rtime       = 0;
  md->num_pairsred      = 0;
  md->num_gb_crit       = 0;
  md->num_redundant     = 0;
  md->num_rowsred       = 0;
  md->num_zerored       = 0;
  md->num_ht_enlarge    = 0;
  md->num_sdm_found     = 0;
  md->num_not_sdm_found = 0;

  set_function_pointers();

  return 0;
}

/* note that depending on the input data we set the corresponding
 * function pointers for monomial resp. spair comparisons, taking
 * spairs by a given minimal property for symbolic preprocessing, etc. */
static mat_t *import_julia_data_16(
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const md_t *md
    )
{
  int32_t i, j;
  int16_t *cf;
  len_t *ch;
  int32_t off = 0; /* offset in arrays */
  
  mat_t *mat  = (mat_t *)malloc(sizeof(mat_t));
  mat->r      = (row_t **)malloc((unsigned long)md->ng * sizeof(row_t *));
  
  for (i = 0; i < md->ng; ++i) {
    mat->r[i] = (row_t *)malloc(sizeof(row_t));
    mat->r[i]->sz = lens[i];
    mat->r[i]->os = lens[i] % UNROLL;

    mat->r[i]->cf = malloc((unsigned long)mat->r[i]->sz * sizeof(int16_t));
    mat->r[i]->ch = (len_t *)malloc((unsigned long)mat->r[i]->sz * sizeof(len_t));

    cf  = mat->r[i]->cf;
    ch  = mat->r[i]->ch;
    for (j = off; j < off+lens[i]; ++j) {
      cf[j] = (int16_t)cfs[j];
      ch[j] = insert_in_global_hash_table(exps+(md->nv*j));
    }
    /* mark initial generators, they have to be added to the basis first */
    off +=  lens[i];
  }
  mat->nr = mat->np = mat->na = md->ng;

  return mat;
}

static mat_t *import_julia_data_32(
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const md_t *md
    )
{
  int32_t i, j;
  int16_t *cf;
  len_t *ch;
  int32_t off = 0; /* offset in arrays */
  
  mat_t *mat  = (mat_t *)malloc(sizeof(mat_t));
  mat->r      = (row_t **)malloc((unsigned long)md->ng * sizeof(row_t *));
  
  for (i = 0; i < md->ng; ++i) {
    mat->r[i] = (row_t *)malloc(sizeof(row_t));
    mat->r[i]->sz = lens[i];
    mat->r[i]->os = lens[i] % UNROLL;

    mat->r[i]->cf = malloc((unsigned long)mat->r[i]->sz * sizeof(int32_t));
    mat->r[i]->ch = (len_t *)malloc((unsigned long)mat->r[i]->sz * sizeof(len_t));

    cf  = mat->r[i]->cf;
    ch  = mat->r[i]->ch;
    for (j = off; j < off+lens[i]; ++j) {
      cf[j] = cfs[j];
      ch[j] = insert_in_global_hash_table(exps+(md->nv*j));
    }
    /* mark initial generators, they have to be added to the basis first */
    off +=  lens[i];
  }
  mat->nr = mat->np = mat->na = md->ng;

  return mat;
}

static int64_t export_julia_data_16(
    const bs_t *bs,
    const md_t *md,
    int32_t **bp
    )
{
  int32_t i, j, k;
  int16_t *cf;
  len_t *ch;
  int64_t ctr_lengths, ctr_elements;
  int32_t *basis  = *bp;

  int64_t len = 0; /* complete length of exported array */
  int64_t nb  = 0; /* # elemnts in basis */

  const int32_t lterm = 1 + md->nv; /* length of a term */

  /* compute number of terms */
  for (i = 0; i < bs->ld; ++i) {
    if ((long)bs->p[i]->cf & bred) {
      continue;
    } else {
      len +=  (int64_t)(bs->p[i]->sz);
      nb++;
    }
  }

  /* compute the length considering the number of variables per exponent */
  len = len * (int64_t)lterm;
  /* add storage for length of each element */
  len = len + nb;
  /* add storage for number of generators in basis */
  len++;

  basis  = (int32_t *)malloc((unsigned long)len * sizeof(int32_t));

  if (nb > (int64_t)(pow(2, 31))) {
    printf("basis too big\n");
    return 0;
  }

  ctr_lengths   = 1;
  ctr_elements  = (int64_t)nb + 1;

  basis[0]  = (int32_t)nb;
  /* basis[1]  = (int32_t)nb; */
  for (i = 0; i < bload; ++i) {
    if ((long)bs->p[i]->cf & bred) {
      continue;
    } else {
      cf  = bs->p[i]->cf;
      ch  = bs->p[i]->ch;
      /* length of polynomial including this length entry itself */
      basis[ctr_lengths++]  = bs->p[i]->sz * lterm;
      for (j = 0; j < bs->p[i]->sz; ++j) {
        basis[ctr_elements++] = (int32_t)cf[j]; /* coefficient */
        for (k = 0; k < md->nv; ++k) {
          basis[ctr_elements++] = (ev + ch[j])[k];
        }
      }
    }
  }
  *bp = basis;

  return len;
}

static int64_t export_julia_data_32(
    const bs_t *bs,
    const md_t *md,
    int32_t **bp
    )
{
  int32_t i, j, k;
  int32_t *cf;
  len_t *ch;
  int64_t ctr_lengths, ctr_elements;
  int32_t *basis  = *bp;

  int64_t len = 0; /* complete length of exported array */
  int64_t nb  = 0; /* # elemnts in basis */

  const int32_t lterm = 1 + md->nv; /* length of a term */

  /* compute number of terms */
  for (i = 0; i < bs->ld; ++i) {
    if ((long)bs->p[i]->cf & bred) {
      continue;
    } else {
      len +=  (int64_t)(bs->p[i]->sz);
      nb++;
    }
  }

  /* compute the length considering the number of variables per exponent */
  len = len * (int64_t)lterm;
  /* add storage for length of each element */
  len = len + nb;
  /* add storage for number of generators in basis */
  len++;

  basis  = (int32_t *)malloc((unsigned long)len * sizeof(int32_t));

  if (nb > (int64_t)(pow(2, 31))) {
    printf("basis too big\n");
    return 0;
  }

  ctr_lengths   = 1;
  ctr_elements  = (int64_t)nb + 1;

  basis[0]  = (int32_t)nb;
  /* basis[1]  = (int32_t)nb; */
  for (i = 0; i < bload; ++i) {
    if ((long)bs->p[i]->cf & bred) {
      continue;
    } else {
      cf  = bs->p[i]->cf;
      ch  = bs->p[i]->ch;
      /* length of polynomial including this length entry itself */
      basis[ctr_lengths++]  = bs->p[i]->sz * lterm;
      for (j = 0; j < bs->p[i]->sz; ++j) {
        basis[ctr_elements++] = cf[j]; /* coefficient */
        for (k = 0; k < md->nv; ++k) {
          basis[ctr_elements++] = (ev + ch[j])[k];
        }
      }
    }
  }
  *bp = basis;

  return len;
}
