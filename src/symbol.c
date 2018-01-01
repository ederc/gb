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
    m       = monomial_division(lcms[i], b[1]);
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
  DEBUG(SELDBG, "\t %6d/%6d pairs, deg %d", n, pload, md);
  num_reduced +=  n;

  /* generate matrix out of selected spairs */
  tmp_lcm = (len_t *)malloc((unsigned long)n * sizeof(len_t));
  tmp_gen = (len_t *)malloc((unsigned long)n * sizeof(len_t));
  for (i = 0, j = 0; i < pload; ++i) {
    if (ps[i].deg != md) {
      continue;
    }
    tmp_gen[j]    = ps[i].gen1;
    tmp_gen[j+1]  = ps[i].gen1;
    tmp_lcm[j]  = tmp_lcm[j+1]  = ps[i].lcm;
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
      if (j == k) {
        tmp_lcm[k]    = tmp_lcm[i];
        tmp_gen[k++]  = tmp_gen[i];
      }
    }
  }
  num_duplicates  +=  n-k;

  n = k;
  
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
