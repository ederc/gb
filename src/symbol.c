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

static void select_spairs_by_minimal_degree(
    mat_t *mat,
    md_t * md,
    const bs_t * const bs
    )
{
  int32_t i, j, k, l, mdeg, npairs;
  row_t *b;
  deg_t d = 0;
  len_t lcm = 0, load = 0, load_old = 0;
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
  md->psort_rtime +=  rrt1 - rrt0;
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
  mat->r  = realloc(mat->r, 2 * (unsigned long)npairs * sizeof(val_t *));
  /* mat   = (val_t **)malloc(2 * (unsigned long)npairs * sizeof(val_t *)); */
  mat->na = 2 * npairs;
  mat->nc = mat->ncl = mat->ncr = 0;
  mat->nr = 0;

  i = 0;
  while (i < npairs) {
    /* ncols initially counts number of different lcms */
    mat->nc++;
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
      b = bs->p[gens[k]];
      /* b = (val_t *)((long)bs[gens[k]] & bmask); */
      /* m = monomial_division_no_check(lcm, b[2]); */
      for (l = 0; l < md->nv; ++l) {
        em[l] = (ev+lcm)[l] - (ev+b->ch[0])[l];
        d     +=  em[l];
      }
      const val_t h = (ev+lcm)[HASH_VAL] - (ev+b->ch[0])[HASH_VAL];
      mat->r[mat->nr]  = multiplied_polynomial_to_matrix_row(h, d, em, b);
      /* mark lcm column as lead term column */
      (ev+mat->r[mat->nr]->ch[0])[HASH_IND] = 2;
      mat->nr++;
    }

    i = j;
  }

  md->num_duplicates  +=  0;
  md->num_rowsred     +=  load;

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
  md->select_ctime  +=  ct1 - ct0;
  md->select_rtime  +=  rt1 - rt0;
}

static inline row_t *find_multiplied_reducer(
    const len_t m,
    const bs_t * const bs,
    const md_t * md
    )
{
  int32_t i, k;
  deg_t d = 0;
  row_t *b;
  const exp_t * const e  = ev+m;
  exp_t *f;
  const len_t * const lms = bs->lm;
  /* exp_t *r  = (exp_t *)malloc((unsigned long)nvars * sizeof(exp_t)); */

  const int32_t bl  = bs->ld;
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
    b = bs->p[i];
    f = ev+b->ch[0];
    if ((e[0]-f[0]) < 0) {
      i++;
      goto start;
    }
    for (k = md->os; k < md->nv; k += 2) {
      if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
        i++;
        goto start;
      }
    }
    exp_t *r = (exp_t *)malloc((unsigned long)nvars * sizeof(exp_t));
    for (i = 0; i < md->nv; ++i) {
      r[i]  =   e[i] - f[i];
      d     +=  r[i];
    }
    const val_t h = e[HASH_VAL] - f[HASH_VAL];
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
    b = bs->p[i];
    f = ev+b->ch[0];
    if ((e[0]-f[0]) < 0) {
      i++;
      goto start2;
    }
    for (k = md->os; k < md->nv; k += 2) {
      if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
        i++;
        goto start2;
      }
    }
    exp_t *r = (exp_t *)malloc((unsigned long)nvars * sizeof(exp_t));
    for (i = 0; i < md->nv; ++i) {
      r[i]  =   e[i] - f[i];
      d     +=  r[i];
    }
    const val_t h = e[HASH_VAL] - f[HASH_VAL];
    b = multiplied_polynomial_to_matrix_row(h, d, r, b);
    free(r);
    return b;
  }
  (ev+m)[HASH_DIV] = i;
  num_not_sdm_found++;
  return NULL;
}

static void symbolic_preprocessing(
    mat_t *mat,
    md_t * md,
    const bs_t * const bs
    )
{
  int32_t i, j;
  row_t *r;
  row_t *red;
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
  for (i = 0; i < mat->nr; ++i) {
    r = mat->r[i];
    const len_t len = r->sz;
    /* printf("%d | %d | %d\n", mat->na, mat->nr, len);  */
    /* check row reallocation only once per polynomial */
    if ((mat->na - mat->nr) < len) {
      mat->na = mat->na > len ?
        2*mat->na : mat->na+len;
      mat->r  = realloc(mat->r, (unsigned long)mat->na * sizeof(row_t *));
    }
    for (j = 1; j < len; ++j) {
      m = r->ch[j];
      if (!(ev+m)[HASH_IND]) {
        (ev+m)[HASH_IND] = 1;
        mat->nc++;
        red = find_multiplied_reducer(m, bs, md);
        if (red) {
          (ev+m)[HASH_IND] = 2;
          /* add new reducer to matrix */
          mat->r[mat->nr]  = red;
          mat->nr++;
        }
      }
    }
  }

  /* realloc to real size */
  mat->r  = realloc(mat->r, (unsigned long)mat->nr * sizeof(row_t *));
  mat->na = mat->nr;

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->symbol_ctime  +=  ct1 - ct0;
  md->symbol_rtime  +=  rt1 - rt0;
}
