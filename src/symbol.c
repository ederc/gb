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
 * \file symbol.c
 * \brief Symbolic preprocessing routines
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

static inline val_t **generate_matrix_from_spair_selection(
    const len_t *lcms,
    const len_t *gens,
    const int32_t nr
    )
{
  int32_t i;
  len_t m;
  val_t *b;
  
  val_t **mat = (val_t **)malloc(2 * (unsigned long)nr * sizeof(val_t *));
  nrall = 2 * nr;
  nrows = nr;

  for (i = 0; i < nr; ++i) {
    b       = (val_t *)((long)bs[gens[i]] & bmask);
    m       = monomial_division_no_check(lcms[i], b[2]);
    mat[i]  = multiplied_polynomial_to_matrix_row(m, b);
  }

  return mat;
}

static val_t **select_spairs_by_minimal_degree(
    void
    )
{
  int32_t i, j, k, md, n;
  val_t **mat;
  len_t *tmp_lcm, *tmp_gen;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* get minimal degree */
  md  = (int32_t)1 << 30;
  n   = 0;

  for (i = 0; i < pload; ++i) {
    if (ps[i].deg > md) {
      continue;
    } else {
      if (ps[i].deg == md) {
        n++;
      } else {
        md  = ps[i].deg;
        n   = 1;
      }
    }
  }
  GB_DEBUG(SELDBG, "\t %6d/%6d pairs, deg %d", n, pload, md);
  /* statistics */
  num_pairsred  +=  n;

  /* generate matrix out of selected spairs */
  tmp_lcm = (len_t *)malloc(2 * (unsigned long)n * sizeof(len_t));
  tmp_gen = (len_t *)malloc(2 * (unsigned long)n * sizeof(len_t));
  for (i = 0, j = 0; i < pload; ++i) {
    if (ps[i].deg != md) {
      continue;
    }
    tmp_gen[j]    = ps[i].gen1;
    tmp_gen[j+1]  = ps[i].gen2;
    tmp_lcm[j]    = tmp_lcm[j+1]  = ps[i].lcm;
    j +=  2;
  }
  n = j;

  /* remove duplicates */
  for (i = 0, k = 0; i < n; ++i) {
    for (j = 0; j < k; ++j) {
      if (tmp_lcm[i] == tmp_lcm[j]
          && tmp_gen[i] == tmp_gen[j]) {
        break;
      }
    }
    if (j == k) {
      tmp_lcm[k]    = tmp_lcm[i];
      tmp_gen[k++]  = tmp_gen[i];
    }
  }
  num_duplicates  +=  n-k;

  n = k;
  
  /* statistics */
  num_rowsred +=  n;

  /* move lcms to local hash table */
  for (i = 0; i < n; ++i) {
    tmp_lcm[i]  = insert_in_local_hash_table(ev+tmp_lcm[i]);
  }
  
  mat = generate_matrix_from_spair_selection(tmp_lcm, tmp_gen, n);

  free(tmp_lcm);
  free(tmp_gen);

  /* remove selected spairs from pairset */
  for (i = 0, j = 0; i < pload; ++i) {
    if (ps[i].deg == md) {
      continue;
    }
    ps[j++] = ps[i];
  }
  pload = j;

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  select_ctime  +=  ct1 - ct0;
  select_rtime  +=  rt1 - rt0;
  
  return mat;
}

static inline val_t *find_multiplied_reducer(
    len_t m
    )
{
  int32_t i;
  len_t d;
  val_t *b;

  exp_t *e  = evl + m;

  for (i = e[HASH_DIV]; i < bload; ++i) {
    b = (val_t *)((long)bs[i] & bmask);
    if (b != bs[i]) {
      continue;
    }
    d = monomial_division_with_check(m, b[2]);
    if (d == 0) {
      continue;
    }
    e[HASH_DIV] = i;
    return multiplied_polynomial_to_matrix_row(m, b);
  }
  e[HASH_DIV] = i;
  return NULL;
}

static val_t **symbolic_preprocessing(
    val_t **mat
    )
{
  int32_t i, j;
  exp_t *e;
  val_t *red;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* mark leading monomials in HASH_IND entry */
  for (i = 0; i < nrows; ++i) {
    (evl + mat[i][2])[HASH_IND] = 1;
  }

  /* get reducers from basis */
  for (i = 0; i < nrows; ++i) {
    const len_t len = mat[i][0];
    for (j = 4; j < len; j += 2) {
      e = evl + mat[i][j];
      if (e[HASH_IND]) {
        continue;
      }
      e[HASH_IND] = 1;
      red = find_multiplied_reducer(mat[i][j]);
      if (!red) {
        continue;
      }
      if (nrows == nrall) {
        nrall = 2 * nrall;
        mat   = realloc(mat, (unsigned long)nrall * sizeof(val_t *));
      }
      /* add new reducer to matrix */
      mat[nrows++]  = red;
    }
  }

  /* clear markers */
  for (i = 0; i < nrows; ++i) {
    const len_t len = mat[i][0];
    for (j = 2; j < len; j = j+2) {
      (evl + mat[i][j])[HASH_IND] = 0;
    }
  }

  /* realloc to real size */
  mat   = realloc(mat, (unsigned long)nrows * sizeof(val_t *));
  nrall = nrows;

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  symbol_ctime  +=  ct1 - ct0;
  symbol_rtime  +=  rt1 - rt0;

  return mat;
}
