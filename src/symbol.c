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

static mat_t *select_spairs_by_minimal_degree(
    const md_t * const md
    )
{
  int32_t i, j, k, l, mdeg, npairs;
  val_t *b;
  deg_t d = 0;
  len_t lcm = 0, load = 0, load_old = 0;
  mat_t *mat;
  len_t *gens;

  exp_t *em = (exp_t *)malloc((unsigned long)md->nv * sizeof(exp_t));

  /* timings */
  double ct0, ct1, rt0, rt1, rrt0, rrt1;
  ct0 = cputime();
  rt0 = realtime();

  /* sort pair set */
  rrt0 = realtime();
  qsort(ps, (unsigned long)pload, sizeof(spair_t), spair_cmp);
  rrt1 = realtime();
  pair_sort_rtime +=  rrt1 - rrt0;
  /* get minimal degree */
  mdeg  = ps[0].deg;

  /* select pairs of this degree respecting maximal selection size mnsel */
  for (i = 0; i < pload; ++i) {
    if (ps[i].deg > mdeg || i >= mnsel) {
      break;
    }
  }
  npairs  = i;
  /* if we stopped due to maximal selection size we still get the following
   * pairs of the same lcm in this matrix */
  if (i > mnsel) {
    j = i+1;
    while (ps[j].lcm == ps[i].lcm) {
      ++j;
    }
    npairs = j;
  }
  GB_DEBUG(SELDBG, " %6d/%6d pairs - deg %2d", npairs, pload, mdeg);
  /* statistics */
  num_pairsred  +=  npairs;
  
  gens  = (len_t *)malloc(2 * (unsigned long)npairs * sizeof(len_t));

  /* preset matrix meta data */
  mat   = (val_t **)malloc(2 * (unsigned long)npairs * sizeof(val_t *));
  nrall = 2 * npairs;
  ncols = ncl = ncr = 0;
  nrows = 0;

  i = 0;
  while (i < npairs) {
    /* ncols initially counts number of different lcms */
    ncols++;
    load_old  = load;
    lcm   = ps[i].lcm;
    gens[load++] = ps[i].gen1;
    gens[load++] = ps[i].gen2;
    j = i+1;
    while (j < npairs && ps[j].lcm == lcm) {
      for (k = load_old; k < load; ++k) {
        if (ps[j].gen1 == gens[k]) {
          break;
        }
      }
      if (k == load) {
        gens[load++]  = ps[j].gen1;
      }
      for (k = load_old; k < load; ++k) {
        if (ps[j].gen2 == gens[k]) {
          break;
        }
      }
      if (k == load) {
        gens[load++]  = ps[j].gen2;
      }
      j++;
    }
    for (k = load_old; k < load; ++k) {
      d = 0;
      b = (val_t *)((long)bs[gens[k]] & bmask);
      /* m = monomial_division_no_check(lcm, b[2]); */
      for (l = 0; l < nvars; ++l) {
        em[l] = (ev+lcm)[l] - (ev+b[2])[l];
        d     +=  em[l];
      }
      const val_t h = (ev+lcm)[HASH_VAL] - (ev+b[2])[HASH_VAL];
      mat[nrows]  = multiplied_polynomial_to_matrix_row(h, d, em, b);
      /* mark lcm column as lead term column */
      (ev+mat[nrows][2])[HASH_IND] = 2;
      nrows++;
    }

    i = j;
  }

  num_duplicates  +=  0;
  num_rowsred     +=  load;

  free(gens);
  free(em);

  /* remove selected spairs from pairset */
  for (j = npairs; j < pload; ++j) {
    ps[j-npairs]  = ps[j];
  }
  pload = pload - npairs;

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
  int32_t i, k;
  deg_t d = 0;
  val_t *b;
  const exp_t * const e  = ev+m;
  exp_t *f;
  /* exp_t *r  = (exp_t *)malloc((unsigned long)nvars * sizeof(exp_t)); */

  const int32_t bl  = bload;
  const int32_t os  = nvars & 1 ? 1 : 0;
  i = e[HASH_DIV];

  const sdm_t ns = ~e[HASH_SDM];
start:
  while (i < bl-3) {
    if (lms[i] & ns &&
        lms[i+1] & ns &&
        lms[i+2] & ns &&
        lms[i+3] & ns) {
      num_sdm_found +=  4;
      i +=  4;
      continue;
    }
    while (lms[i] & ns) {
      i++;
    }
    b = (val_t *)((long)bs[i] & bmask);
    f = ev+b[2];
    if ((e[0]-f[0]) < 0) {
      i++;
      goto start;
    }
    for (k = os; k < nvars; k += 2) {
      if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
        i++;
        goto start;
      }
    }
    exp_t *r = (exp_t *)malloc((unsigned long)nvars * sizeof(exp_t));
    for (i = 0; i < nvars; ++i) {
      r[i]  = e[i] - f[i];
    }
    const val_t h = e[HASH_VAL] - f[HASH_VAL];
    for (i = 0; i < nvars; ++i) {
      d += r[i];
    }
    b = multiplied_polynomial_to_matrix_row(h, d, r, b);
    free(r);
    return b;
  }
start2:
  while (i < bl) {
    if (lms[i] & ns) {
      num_sdm_found++;
      i++;
      continue;
    }
    b = (val_t *)((long)bs[i] & bmask);
    f = ev+b[2];
    if ((e[0]-f[0]) < 0) {
      i++;
      goto start2;
    }
    for (k = os; k < nvars; k += 2) {
      if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
        i++;
        goto start2;
      }
    }
    exp_t *r = (exp_t *)malloc((unsigned long)nvars * sizeof(exp_t));
    for (i = 0; i < nvars; ++i) {
      r[i]  = e[i] - f[i];
    }
    const val_t h = e[HASH_VAL] - f[HASH_VAL];
    for (i = 0; i < nvars; ++i) {
      d += r[i];
    }
    b = multiplied_polynomial_to_matrix_row(h, d, r, b);
    free(r);
    return b;
  }
  (ev+m)[HASH_DIV] = i;
  num_not_sdm_found++;
  return NULL;
}

static val_t **symbolic_preprocessing(
    val_t **mat
    )
{
  int32_t i, j;
  val_t *red;
  val_t m;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* note that we have already counted the different lcms, i.e.
   * ncols until this step. moreover, we have also already marked
   * the corresponding hash indices to represent lead terms. so
   * we only have to do the bookkeeping for newly added reducers
   * in the following. */

  /* get reducers from basis */
  for (i = 0; i < nrows; ++i) {
    const len_t len = mat[i][0];
    /* check row reallocation only once per polynomial */
    if ((nrall - nrows) < (mat[i][0]-4)/2) {
      nrall = 2*nrall > mat[i][0] ? 2*nrall : mat[i][0];
      mat   = realloc(mat, (unsigned long)nrall * sizeof(val_t *));
    }
    for (j = 4; j < len; j += 2) {
      m = mat[i][j];
      if (!(ev+m)[HASH_IND]) {
        (ev+m)[HASH_IND] = 1;
        ncols++;
        red = find_multiplied_reducer(m);
        if (red) {
          (ev+m)[HASH_IND] = 2;
          /* add new reducer to matrix */
          mat[nrows++]  = red;
        }
      }
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
