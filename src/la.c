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
 * \file la.c
 * \brief Implementation of linear algebra.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "data.h"

static inline void normalize_matrix_row_16(
    row_t *row,
    const md_t * const md
    )
{
  int32_t i;
  int16_t *cf = row->cf;
  int64_t tmp1, tmp2, tmp3, tmp4;

  const int32_t inv = mod_p_inverse_32(cf[0], fc);

  for (i = 1; i < row->os; ++i) {
    tmp1  =   ((int64_t)cf[i] * inv) % md->fc;
    tmp1  +=  (tmp1 >> 63) & md->fc;
    cf[i] =   (int16_t)tmp1;
  }
  /* printf("os %d | sz %d | i %d\n", row->os, row->sz, i); */
  /* UNROLL should be 4 here */
  for (i = row->os; i < row->sz; i += 4) {
    tmp1    =   ((int64_t)cf[i] * inv) % md->fc;
    tmp2    =   ((int64_t)cf[i+1] * inv) % md->fc;
    tmp3    =   ((int64_t)cf[i+2] * inv) % md->fc;
    tmp4    =   ((int64_t)cf[i+3] * inv) % md->fc;
    tmp1    +=  (tmp1 >> 63) & md->fc;
    tmp2    +=  (tmp2 >> 63) & md->fc;
    tmp3    +=  (tmp3 >> 63) & md->fc;
    tmp4    +=  (tmp4 >> 63) & md->fc;
    cf[i]   =   (int16_t)tmp1;
    cf[i+1] =   (int16_t)tmp2;
    cf[i+2] =   (int16_t)tmp3;
    cf[i+3] =   (int16_t)tmp4;
  }
  cf[0] = 1;
}

static inline void normalize_matrix_row_32(
    row_t *row,
    const md_t * const md
    )
{
  int32_t i;
  int32_t *cf = row->cf;
  int64_t tmp1, tmp2, tmp3, tmp4;

  const int32_t inv = mod_p_inverse_32(cf[0], fc);

  for (i = 1; i < row->os; ++i) {
    tmp1  =   ((int64_t)cf[i] * inv) % md->fc;
    tmp1  +=  (tmp1 >> 63) & md->fc;
    cf[i] =   (int32_t)tmp1;
  }
  /* UNROLL should be 4 here */
  /* printf("os %d | sz %d | i %d\n", row->os, row->sz, i); */
  for (i = row->os; i < row->sz; i += 4) {
    tmp1    =   ((int64_t)cf[i] * inv) % md->fc;
    tmp2    =   ((int64_t)cf[i+1] * inv) % md->fc;
    tmp3    =   ((int64_t)cf[i+2] * inv) % md->fc;
    tmp4    =   ((int64_t)cf[i+3] * inv) % md->fc;
    tmp1    +=  (tmp1 >> 63) & md->fc;
    tmp2    +=  (tmp2 >> 63) & md->fc;
    tmp3    +=  (tmp3 >> 63) & md->fc;
    tmp4    +=  (tmp4 >> 63) & md->fc;
    cf[i]   =   (int32_t)tmp1;
    cf[i+1] =   (int32_t)tmp2;
    cf[i+2] =   (int32_t)tmp3;
    cf[i+3] =   (int32_t)tmp4;
  }
  cf[0] = 1;

  /* for (i = 0; i < row->sz; ++i) {
   *   printf("%d ", cf[i]);
   * }
   * printf("\n"); */
}

static row_t *reduce_dense_row_by_known_pivots_16(
    int64_t *dr,
    const len_t dpiv,         /* pivot of dense row at the beginning */
    const len_t nc,
    row_t ** const pivs,
    const md_t * const md
    )
{
  int32_t i, j, k;
  const int64_t mod = (int64_t)md->fc;

  for (k = 0, i = dpiv; i < nc; ++i) {
    if (dr[i] != 0) {
      dr[i] = dr[i] % mod;
    }
    if (dr[i] == 0) {
      continue;
    }
    if (pivs[i] == NULL) {
      k++;
      continue;
    }

    /* found reducer row, get multiplier */
    const int64_t mul = mod - dr[i];

    const int16_t * const rcf = (int16_t *)pivs[i]->cf;
    const len_t * const rch   = pivs[i]->ch;
    const len_t os  = pivs[i]->os;
    const len_t sz  = pivs[i]->sz;
    /* printf("dr: ");
     * for (j = 0; j < nc; ++j) {
     *   printf("%ld ", dr[j]);
     * }
     * printf("\n"); */
    for (j = 0; j < os; ++j) {
      dr[rch[j]]  +=  mul * (int64_t)rcf[j];
    }
    /* UNROLL is set to be 4 in src/data.h */
    for (; j < sz; j += 4) {
      dr[rch[j]]    +=  mul * (int64_t)rcf[j];
      dr[rch[j+1]]  +=  mul * (int64_t)rcf[j+1];
      dr[rch[j+2]]  +=  mul * (int64_t)rcf[j+2];
      dr[rch[j+3]]  +=  mul * (int64_t)rcf[j+3];
    }
    dr[i] = 0;
  }
  if (k == 0) {
    return NULL;
  }
  /* printf("dr done: ");
   * for (j = 0; j < nc; ++j) {
   *   printf("%ld ", dr[j]);
   * }
   * printf("\n"); */

  /* dense row is not reduced to zero, thus generate new sparse
   * pivot row and normalize it */
  row_t *row  = (row_t *)malloc(sizeof(row_t));
  row->sz     = nc-dpiv;
  row->cf     = (int16_t *)malloc((unsigned long)row->sz * sizeof(int16_t));
  row->ch     = (len_t *)malloc((unsigned long)row->sz * sizeof(len_t));

  int16_t *cf = (int16_t *)row->cf;
  len_t *ch   = row->ch;
  j = 0;
  for (i = dpiv; i < nc; ++i) {
    if (dr[i] != 0) {
      dr[i] = dr[i] % mod; // IS THERE A BETTER MODULUS OPERATOR GIVING THE SMALLEST VALUE W.R.T. THE ABSOLUTE VALUE
      if (dr[i] != 0) {
        if (dr[i] > pow(2, 15)) {
          dr[i] = dr[i] - mod;
        }
        cf[j]  = (int16_t)dr[i];
        ch[j]  = (len_t)i;
        ++j;
      }
    }
  }
  row->sz = j;
  row->os = row->sz % UNROLL;
  cf = realloc(cf, (unsigned long)row->sz * sizeof(int16_t));
  ch = realloc(ch, (unsigned long)row->sz * sizeof(len_t));

  row->cf = cf;
  row->ch = ch;

  if (((int16_t *)row->cf)[0] != 1) {
    normalize_matrix_row(row, md);
  }

  return row;
}

static row_t *reduce_dense_row_by_known_pivots_32(
    int64_t *dr,
    const len_t dpiv,         /* pivot of dense row at the beginning */
    const len_t nc,
    row_t ** const pivs,
    const md_t * const md
    )
{
  int32_t i, j, k;
  const int64_t mod   = (int64_t)md->fc;
  const int64_t mod2  = (int64_t)md->fc * md->fc;

  for (k = 0, i = dpiv; i < nc; ++i) {
    if (dr[i] != 0) {
      dr[i] = dr[i] % mod;
    }
    if (dr[i] == 0) {
      continue;
    }
    if (pivs[i] == NULL) {
      k++;
      continue;
    }

    /* printf("dr to be reduced: ");
     * for (int32_t p = 0; p < nc; ++p) {
     *   printf("%ld ", dr[p]);
     * }
     * printf("\n"); */

    /* found reducer row, get multiplier */
    const int64_t mul = (int64_t)dr[i];

    const int32_t * const cf  = (int32_t *)pivs[i]->cf;
    const len_t * const ch    = pivs[i]->ch;
    const len_t os  = pivs[i]->os;
    const len_t sz  = pivs[i]->sz;

    /* printf("reducer: ");
     * for (int32_t p = 0; p < sz; ++p) {
     *   printf("%d | %d  ", cf[p], ch[p]);
     * }
     * printf("\n"); */
    for (j = 0; j < os; ++j) {
      dr[ch[j]]  -=  mul * (int64_t)cf[j];
      dr[ch[j]]  +=  (dr[ch[j]] >> 63) & mod2;
    }
    /* UNROLL is set to be 4 in src/data.h */
    for (; j < sz; j += 4) {
      dr[ch[j]]   -=  mul * (int64_t)cf[j];
      dr[ch[j]]   +=  (dr[ch[j]] >> 63) & mod2;
      dr[ch[j+1]] -=  mul * (int64_t)cf[j+1];
      dr[ch[j+1]] +=  (dr[ch[j+1]] >> 63) & mod2;
      dr[ch[j+2]] -=  mul * (int64_t)cf[j+2];
      dr[ch[j+2]] +=  (dr[ch[j+2]] >> 63) & mod2;
      dr[ch[j+3]] -=  mul * (int64_t)cf[j+3];
      dr[ch[j+3]] +=  (dr[ch[j+3]] >> 63) & mod2;
    }
    dr[i] = 0;
  }
  if (k == 0) {
    return NULL;
  }
    /* printf("dr finally reduced: ");
     * for (i = 0; i < nc; ++i) {
     *   printf("%ld ", dr[i]);
     * }
     * printf("\n"); */

  /* dense row is not reduced to zero, thus generate new sparse
   * pivot row and normalize it */
  row_t *row  = (row_t *)malloc(sizeof(row_t));
  row->sz     = nc-dpiv;
  row->cf     = (int32_t *)malloc((unsigned long)row->sz * sizeof(int32_t));
  row->ch     = (len_t *)malloc((unsigned long)row->sz * sizeof(len_t));

  int32_t *cf = row->cf;
  len_t *ch   = row->ch;
  j = 0;
  for (i = dpiv; i < nc; ++i) {
    if (dr[i] != 0) {
      dr[i] = dr[i] % mod;
      if (dr[i] != 0) {
        cf[j]  = (int32_t)dr[i];
        ch[j]  = (len_t)i;
        /* printf("cf[%d] = %d | ch[%d] = %d\n", j, cf[j], j, ch[j]); */
        ++j;
      }
    }
  }
  row->sz = j;
  row->os = row->sz % UNROLL;
  cf = realloc(cf, (unsigned long)row->sz * sizeof(int32_t));
  ch = realloc(ch, (unsigned long)row->sz * sizeof(len_t));

  row->cf = cf;
  row->ch = ch;

  if (((int32_t *)row->cf)[0] != 1) {
    normalize_matrix_row(row, md);
  }

  return row;
}

static void sparse_linear_algebra_16(
    mat_t *mat,
    md_t *md
    )
{
  int32_t i, j, k;
  len_t sc    = 0;    /* starting column */
  int32_t sz  = 0;
  row_t *npiv = NULL; /* new pivot row */
  row_t *row  = NULL;
  int16_t *cf = NULL;
  len_t *ch   = NULL;

  const len_t nrl = mat->nrl;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* all pivots, first we can only fill in all known lead terms */
  row_t **pivs  = (row_t **)calloc((unsigned long)mat->nc, sizeof(row_t *));
  /* unkown pivot rows we have to reduce with the known pivots first */
  row_t **upivs = (row_t **)malloc((unsigned long)nrl * sizeof(row_t *));

  i = 0;
  j = 1;
  for (i = 0; i < mat->nr; ++i) {
    row = mat->r[i];
    if (!pivs[row->ch[0]]) {
      pivs[row->ch[0]]   = row;
    } else {
      /* shorter rows first */
      upivs[nrl-j]  = row;
      j++;
    }
  }

  int64_t *dr  = (int64_t *)malloc(
      (unsigned long)(md->nt * mat->nc) * sizeof(int64_t));
#pragma omp parallel for num_threads(md->nt) \
  private(i, j, k, npiv, cf, ch, sc, sz) shared(pivs)
  for (i = 0; i < mat->nrl; ++i) {
    int64_t *drl  = dr + (omp_get_thread_num() * mat->nc);
    npiv  = upivs[i];
    k   = 0;
    /* do the reduction */
    do {
      cf  = (int16_t *)npiv->cf;
      ch  = npiv->ch;
      sz  = npiv->sz;
      sc  = ch[0];
      memset(drl, 0, (unsigned long)mat->nc * sizeof(int64_t));
      for (j = 0; j < sz; ++j) {
        drl[ch[j]] = (int64_t)cf[j];
      }
      free(npiv->cf);
      free(npiv->ch);
      free(npiv);
      npiv  = reduce_dense_row_by_known_pivots_16(drl, sc, mat->nc, pivs, md);
      if (!npiv) {
        break;
      }
      k  = __sync_bool_compare_and_swap(&pivs[npiv->ch[0]], NULL, npiv);
    } while (!k);
  }

  free(upivs);
  upivs = NULL;

  /* we do not need the old pivots anymore */
  for (i = 0; i < ncl; ++i) {
    free(pivs[i]->cf);
    free(pivs[i]->ch);
    free(pivs[i]);
    pivs[i] = NULL;
  }

  npivs = 0; /* number of new pivots */

  dr      = realloc(dr, (unsigned long)mat->nc * sizeof(int64_t));
  mat->r  = realloc(mat->r, (unsigned long)(mat->ncr) * sizeof(row_t *));
  /* interreduce new pivots, i.e. pivs[ncl + ...] */
  for (i = (mat->nc-1); i >= mat->nru; --i) {
    if (pivs[i]) {
      memset(dr, 0, (unsigned long)mat->nc * sizeof(int64_t));
      cf  = (int16_t *)pivs[i]->cf;
      ch  = pivs[i]->ch;
      sz  = pivs[i]->sz;
      sc  = ch[0];
      for (j = 0; j < sz; ++j) {
        dr[ch[j]] = (int64_t)cf[j];
      }
      free(pivs[i]->cf);
      free(pivs[i]->ch);
      free(pivs[i]);
      pivs[i] = NULL;
      pivs[i] = mat->r[mat->np++] =
        reduce_dense_row_by_known_pivots_16(dr, sc, mat->nc, pivs, md);
    }
  }
  free(pivs);
  pivs  = NULL;

  free(dr);
  dr  = NULL;
  mat->r  = realloc(mat->r, (unsigned long)mat->np * sizeof(row_t *));
  mat->nr = mat->na = mat->np;

  md->num_zerored += (mat->nrl - mat->np);
  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->la_ctime  +=  ct1 - ct0;
  md->la_rtime  +=  rt1 - rt0;

  GB_DEBUG(LADBG, "%7d new %7d zero - %9.3f sec",
      mat->np, mat->nrl-mat->np, rt1-rt0);
}

static void sparse_linear_algebra_32(
    mat_t *mat,
    md_t *md
    )
{
  int32_t i, j, k;
  len_t sc    = 0;    /* starting column */
  int32_t sz  = 0;
  row_t *npiv = NULL; /* new pivot row */
  row_t *row  = NULL;
  int32_t *cf = NULL;
  len_t *ch   = NULL;

  const len_t nrl = mat->nrl;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* all pivots, first we can only fill in all known lead terms */
  row_t **pivs  = (row_t **)calloc((unsigned long)mat->nc, sizeof(row_t *));
  /* unkown pivot rows we have to reduce with the known pivots first */
  row_t **upivs = (row_t **)malloc((unsigned long)mat->nrl * sizeof(row_t *));

  i = 0;
  j = 1;
  for (i = 0; i < mat->nr; ++i) {
    row = mat->r[i];
    if (!pivs[row->ch[0]]) {
      pivs[row->ch[0]]   = row;
    } else {
      /* shorter rows first */
      upivs[nrl-j]  = row;
      j++;
    }
  }

  mat->np = 0;

  int64_t *dr  = (int64_t *)malloc(
      (unsigned long)(md->nt * mat->nc) * sizeof(int64_t));
#pragma omp parallel for num_threads(md->nt) \
  private(i, j, k, npiv, cf, ch, sc, sz) shared(pivs)
  for (i = 0; i < mat->nrl; ++i) {
    int64_t *drl  = dr + (omp_get_thread_num() * mat->nc);
    npiv  = upivs[i];
    k   = 0;
    /* do the reduction */
    do {
      cf  = (int32_t *)npiv->cf;
      ch  = npiv->ch;
      sz  = npiv->sz;
      sc  = ch[0];
      memset(drl, 0, (unsigned long)mat->nc * sizeof(int64_t));
      for (j = 0; j < sz; ++j) {
        drl[ch[j]] = (int64_t)cf[j];
      }
      free(cf);
      free(ch);
      free(npiv);
      npiv  = reduce_dense_row_by_known_pivots_32(drl, sc, mat->nc, pivs, md);
      if (!npiv) {
        break;
      }
      k  = __sync_bool_compare_and_swap(&pivs[npiv->ch[0]], NULL, npiv);
    } while (!k);
  }

  free(upivs);
  upivs = NULL;

  /* we do not need the old pivots anymore */
  for (i = 0; i < ncl; ++i) {
    free(pivs[i]->cf);
    free(pivs[i]->ch);
    free(pivs[i]);
    pivs[i] = NULL;
  }

  npivs = 0; /* number of new pivots */

  dr      = realloc(dr, (unsigned long)mat->nc * sizeof(int64_t));
  mat->r  = realloc(mat->r, (unsigned long)(mat->ncr) * sizeof(row_t *));
  /* interreduce new pivots, i.e. pivs[ncl + ...] */
  for (i = (mat->nc-1); i >= mat->nru; --i) {
    if (pivs[i]) {
      memset(dr, 0, (unsigned long)mat->nc * sizeof(int64_t));
      cf  = (int32_t *)pivs[i]->cf;
      ch  = pivs[i]->ch;
      sz  = pivs[i]->sz;
      for (j = 0; j < sz; ++j) {
        dr[ch[j]] = (int64_t)cf[j];
      }
      val_t sc  = ch[0];
      free(pivs[i]->cf);
      free(pivs[i]->ch);
      free(pivs[i]);
      pivs[i] = NULL;
      pivs[i] = mat->r[mat->np++] =
        reduce_dense_row_by_known_pivots_32(dr, sc, mat->nc, pivs, md);
    }
  }
  free(pivs);
  pivs  = NULL;

  free(dr);
  dr  = NULL;
  mat->r  = realloc(mat->r, (unsigned long)mat->np * sizeof(row_t *));
  mat->nr = mat->na = mat->np;

  md->num_zerored += (mat->nrl - mat->np);
  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->la_ctime  +=  ct1 - ct0;
  md->la_rtime  +=  rt1 - rt0;

  GB_DEBUG(LADBG, "%7d new %7d zero - %9.3f sec",
      mat->np, mat->nrl-mat->np, rt1-rt0);
}

static void probabilistic_sparse_linear_algebra_16(
    mat_t *mat,
    md_t *md
    )
{
  int32_t i, j, k, l;
  len_t sc    = 0;    /* starting column */
  int32_t sz  = 0;
  row_t *npiv = NULL; /* new pivot row */
  row_t *row  = NULL;
  int16_t *cf = NULL;
  len_t *ch   = NULL;

  const int64_t mod2  = (int64_t)md->fc * md->fc;

  const len_t nrl = mat->nrl;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* all pivots, first we can only fill in all known lead terms */
  row_t **pivs  = (row_t **)calloc((unsigned long)mat->nc, sizeof(row_t *));
  /* unkown pivot rows we have to reduce with the known pivots first */
  row_t **upivs = (row_t **)malloc((unsigned long)nrl * sizeof(row_t *));

  i = 0;
  j = 1;
  for (i = 0; i < mat->nr; ++i) {
    row = mat->r[i];
    if (!pivs[row->ch[0]]) {
      pivs[row->ch[0]]  = row;
    } else {
      /* shorter rows first */
      upivs[nrl-j]  = row;
      j++;
    }
  }

  /* compute rows per block */
  const int32_t nb  = (int32_t)(floor(sqrt(nrl/2))) > 0 ?
    (int32_t)(floor(sqrt(nrl/2))) :
    (int32_t)(floor(sqrt(nrl))) ;
  const int32_t rem = (nrl % nb == 0) ? 0 : 1;
  const int32_t rpb = (nrl / nb) + rem;

  int64_t *dr   = (int64_t *)malloc(
      (unsigned long)(md->nt * mat->nc) * sizeof(int64_t));
  int64_t *mul  = (int64_t *)malloc(
      (unsigned long)(md->nt * rpb) * sizeof(int64_t));

#pragma omp parallel for num_threads(md->nt) \
  private(i, j, k,l, sc, npiv) shared(pivs)
  for (i = 0; i < nb; ++i) {
    int64_t *drl  = dr + (omp_get_thread_num() * mat->nc);
    int64_t *mull = mul + (omp_get_thread_num() * rpb);

    const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
    const int32_t nrbl  = (int32_t) (nbl - i*rpb);

    sc  = 0;

    if (nrbl != 0) {
      int32_t bctr  = 0;
      while (bctr < nrbl) {
        /* fill random value array */
        for (j = 0; j < nrbl; ++j) {
          mull[j] = (int64_t)rand() % md->fc;
        }
        
        /* generate one dense row as random linear combination
         * of the rows of the block */
        memset(drl, 0, (unsigned long)mat->nc * sizeof(int64_t));
        for (l = 0, j = i*rpb; j < nbl; ++j) {
          cf  = (int16_t *)pivs[j]->cf;
          ch  = pivs[j]->ch;
          sz  = pivs[j]->sz;
          sc  = sc < ch[0] ? sc : ch[0];
          for (k = 0; k < sz; ++k) {
            /* we need to do it in this way in order to support 32-bit primes */
            drl[ch[k]]  -=  mull[l] * cf[k];
            drl[ch[k]]  +=  (drl[ch[k]] >> 63) & mod2;
          }
          l++;
        }
        k = 0;
        /* do the reduction */
        do {
          npiv = reduce_dense_row_by_known_pivots_16(drl, sc, mat->nc, pivs, md);
          if (!npiv) {
            bctr  = nrbl;
            break;
          }
          k  = __sync_bool_compare_and_swap(&pivs[npiv->ch[0]], NULL, npiv);
          if (!k) {
            memset(drl, 0, (unsigned long)mat->nc * sizeof(int64_t));
            cf  = (int16_t *)npiv->cf;
            ch  = npiv->ch;
            sz  = npiv->sz;
            sc  = ch[0];
            for (j = 0; j < sz; ++j) {
              drl[ch[j]] = (int64_t)cf[j];
            }
            free(npiv->cf);
            free(npiv->ch);
            free(npiv);
          }
        } while (!k);
        bctr++;
      }
      for (j = i*rpb; j < nbl; ++j) {
        free(upivs[j]->cf);
        free(upivs[j]->ch);
        free(upivs[j]);
        upivs[j]  = NULL;
      }
    }
  }

  free(mul);
  mul = NULL;

  free(upivs);
  upivs = NULL;

  /* we do not need the old pivots anymore */
  for (i = 0; i < mat->ncl; ++i) {
    free(pivs[i]->cf);
    free(pivs[i]->ch);
    free(pivs[i]);
    pivs[i] = NULL;
  }

  dr  = realloc(dr, (unsigned long)mat->nc * sizeof(int64_t));
  npivs = 0; /* number of new pivots */

  mat->r = realloc(mat->r, (unsigned long)(mat->ncr) * sizeof(row_t *));
  /* interreduce new pivots, i.e. pivs[ncl + ...] */
  for (i = (mat->nc-1); i >= mat->nru; --i) {
    if (pivs[i]) {
      memset(dr, 0, (unsigned long)mat->nc * sizeof(int64_t));
      cf  = (int16_t *)pivs[i]->cf;
      ch  = pivs[i]->ch;
      sz  = pivs[i]->sz;
      sc  = ch[0];
      for (j = 0; j < sz; ++j) {
        dr[ch[j]] = (int64_t)cf[j];

      }
      free(pivs[i]->cf);
      free(pivs[i]->ch);
      free(pivs[i]);
      pivs[i] = NULL;
      pivs[i] = mat->r[npivs++] =
        reduce_dense_row_by_known_pivots_16(dr, sc, mat->nc, pivs, md);
    }
  }
  free(dr);
  free(pivs);
  dr    = NULL;
  pivs  = NULL;

  mat->r  = realloc(mat->r, (unsigned long)npivs * sizeof(val_t *));
  mat->np = mat->nr = mat->na = npivs;

  md->num_zerored += (mat->nrl - npivs);
  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->la_ctime  +=  ct1 - ct0;
  md->la_rtime  +=  rt1 - rt0;

  GB_DEBUG(LADBG, "%7d new %7d zero - %9.3f sec",
      mat->np, mat->nrl-mat->np, rt1-rt0);
}

static void probabilistic_sparse_linear_algebra_32(
    mat_t *mat,
    md_t *md
    )
{
  int32_t i, j, k, l;
  len_t sc    = 0;    /* starting column */
  int32_t sz  = 0;
  row_t *npiv = NULL; /* new pivot row */
  row_t *row  = NULL;
  int32_t *cf = NULL;
  len_t *ch   = NULL;

  const int64_t mod2  = (int64_t)md->fc * md->fc;

  const len_t nrl = mat->nrl;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* all pivots, first we can only fill in all known lead terms */
  row_t **pivs  = (row_t **)calloc((unsigned long)mat->nc, sizeof(row_t *));
  /* unkown pivot rows we have to reduce with the known pivots first */
  row_t **upivs = (row_t **)malloc((unsigned long)nrl * sizeof(row_t *));

  i = 0;
  j = 1;
  for (i = 0; i < mat->nr; ++i) {
    row = mat->r[i];
    if (!pivs[row->ch[0]]) {
      pivs[row->ch[0]]  = row;
    } else {
      /* shorter rows first */
      upivs[nrl-j]  = row;
      j++;
    }
  }

  /* compute rows per block */
  const int32_t nb  = (int32_t)(floor(sqrt(nrl/2))) > 0 ?
    (int32_t)(floor(sqrt(nrl/2))) :
    (int32_t)(floor(sqrt(nrl))) ;
  const int32_t rem = (nrl % nb == 0) ? 0 : 1;
  const int32_t rpb = (nrl / nb) + rem;

  int64_t *dr   = (int64_t *)malloc(
      (unsigned long)(md->nt * mat->nc) * sizeof(int64_t));
  int64_t *mul  = (int64_t *)malloc(
      (unsigned long)(md->nt * rpb) * sizeof(int64_t));

#pragma omp parallel for num_threads(md->nt) \
  private(i, j, k,l, sc, npiv) shared(pivs)
  for (i = 0; i < nb; ++i) {
    int64_t *drl  = dr + (omp_get_thread_num() * mat->nc);
    int64_t *mull = mul + (omp_get_thread_num() * rpb);

    const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
    const int32_t nrbl  = (int32_t) (nbl - i*rpb);

    sc  = 0;

    if (nrbl != 0) {
      int32_t bctr  = 0;
      while (bctr < nrbl) {
        /* fill random value array */
        for (j = 0; j < nrbl; ++j) {
          mull[j] = (int64_t)rand() % md->fc;
        }
        
        /* generate one dense row as random linear combination
         * of the rows of the block */
        memset(drl, 0, (unsigned long)mat->nc * sizeof(int64_t));
        for (l = 0, j = i*rpb; j < nbl; ++j) {
          cf  = (int32_t *)pivs[j]->cf;
          ch  = pivs[j]->ch;
          sz  = pivs[j]->sz;
          sc  = sc < ch[0] ? sc : ch[0];
          for (k = 0; k < sz; ++k) {
            /* we need to do it in this way in order to support 32-bit primes */
            drl[ch[k]]  -=  mull[l] * cf[k];
            drl[ch[k]]  +=  (drl[ch[k]] >> 63) & mod2;
          }
          l++;
        }
        k = 0;
        /* do the reduction */
        do {
          npiv = reduce_dense_row_by_known_pivots_32(drl, sc, mat->nc, pivs, md);
          if (!npiv) {
            bctr  = nrbl;
            break;
          }
          k  = __sync_bool_compare_and_swap(&pivs[npiv->ch[0]], NULL, npiv);
          if (!k) {
            memset(drl, 0, (unsigned long)mat->nc * sizeof(int64_t));
            cf  = (int32_t *)npiv->cf;
            ch  = npiv->ch;
            sz  = npiv->sz;
            sc  = ch[0];
            for (j = 0; j < sz; ++j) {
              drl[ch[j]] = (int64_t)cf[j];
            }
            free(npiv->cf);
            free(npiv->ch);
            free(npiv);
          }
        } while (!k);
        bctr++;
      }
      for (j = i*rpb; j < nbl; ++j) {
        free(upivs[j]->cf);
        free(upivs[j]->ch);
        free(upivs[j]);
        upivs[j]  = NULL;
      }
    }
  }

  free(mul);
  mul = NULL;

  free(upivs);
  upivs = NULL;

  /* we do not need the old pivots anymore */
  for (i = 0; i < mat->ncl; ++i) {
    free(pivs[i]->cf);
    free(pivs[i]->ch);
    free(pivs[i]);
    pivs[i] = NULL;
  }

  dr  = realloc(dr, (unsigned long)mat->nc * sizeof(int64_t));
  npivs = 0; /* number of new pivots */

  mat->r = realloc(mat->r, (unsigned long)(mat->ncr) * sizeof(row_t *));
  /* interreduce new pivots, i.e. pivs[ncl + ...] */
  for (i = (mat->nc-1); i >= mat->nru; --i) {
    if (pivs[i]) {
      memset(dr, 0, (unsigned long)mat->nc * sizeof(int64_t));
      cf  = (int32_t *)pivs[i]->cf;
      ch  = pivs[i]->ch;
      sz  = pivs[i]->sz;
      sc  = ch[0];
      for (j = 0; j < sz; ++j) {
        dr[ch[j]] = (int64_t)cf[j];

      }
      free(pivs[i]->cf);
      free(pivs[i]->ch);
      free(pivs[i]);
      pivs[i] = NULL;
      pivs[i] = mat->r[npivs++] =
        reduce_dense_row_by_known_pivots_32(dr, sc, mat->nc, pivs, md);
    }
  }
  free(dr);
  free(pivs);
  dr    = NULL;
  pivs  = NULL;

  mat->r  = realloc(mat->r, (unsigned long)npivs * sizeof(val_t *));
  mat->np = mat->nr = mat->na = npivs;

  md->num_zerored += (mat->nrl - npivs);
  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->la_ctime  +=  ct1 - ct0;
  md->la_rtime  +=  rt1 - rt0;

  GB_DEBUG(LADBG, "%7d new %7d zero - %9.3f sec",
      mat->np, mat->nrl-mat->np, rt1-rt0);
}
