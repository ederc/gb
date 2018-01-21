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

static inline void normalize_matrix_row(
    val_t *row
    )
{
  int32_t i;
  int64_t tmp1, tmp2, tmp3, tmp4;

  const int32_t inv = mod_p_inverse_32(row[3], fc);
  
  for (i = 3; i < row[1]; i += 2) {
    tmp1    =   ((int64_t)row[i] * inv) % fc;
    tmp1    +=  (tmp1 >> 63) & fc;
    row[i]  =   (val_t)tmp1;
  }
  for (; i < row[0]; i += 8) {
    tmp1      =   ((int64_t)row[i] * inv) % fc;
    tmp2      =   ((int64_t)row[i+2] * inv) % fc;
    tmp3      =   ((int64_t)row[i+4] * inv) % fc;
    tmp4      =   ((int64_t)row[i+6] * inv) % fc;
    tmp1      +=  (tmp1 >> 63) & fc;
    tmp2      +=  (tmp2 >> 63) & fc;
    tmp3      +=  (tmp3 >> 63) & fc;
    tmp4      +=  (tmp4 >> 63) & fc;
    row[i]    =   (val_t)tmp1;
    row[i+2]  =   (val_t)tmp2;
    row[i+4]  =   (val_t)tmp3;
    row[i+6]  =   (val_t)tmp4;
  }
  row[3]  = 1;
}

static val_t *reduce_dense_row_by_known_pivots_16_bit(
    int64_t *dr,
    val_t *const *pivs,
    const val_t dpiv    /* pivot of dense row at the beginning */
    )
{
  int32_t i, j, k;
  const int64_t mod = (int64_t)fc;

  for (k = 0, i = dpiv; i < ncols; ++i) {
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

    /* printf("16-bit dr before reduction with piv %d\n", i);
     * for (int32_t l = 0; l < ncols; ++l) {
     *   printf("%ld ", dr[l]);
     * }
     * printf("\n"); */
    /* found reducer row, get multiplier */
    const int64_t mul = mod - dr[i];
    /* printf("mul %ld\n", mul); */

    for (j = 2; j < pivs[i][1]; j += 2) {
      dr[pivs[i][j]]  +=  mul * (int64_t)pivs[i][j+1];
    }
    /* UNROLL is set to be 4 in src/data.h */
    for (; j < pivs[i][0]; j += 8) {
      dr[pivs[i][j]]    +=  mul * (int64_t)pivs[i][j+1];
      dr[pivs[i][j+2]]  +=  mul * (int64_t)pivs[i][j+3];
      dr[pivs[i][j+4]]  +=  mul * (int64_t)pivs[i][j+5];
      dr[pivs[i][j+6]]  +=  mul * (int64_t)pivs[i][j+7];
    }
    dr[i] = 0;
    /* printf("dr after reduction with piv %d\n", i);
     * for (int32_t l = 0; l < ncols; ++l) {
     *   printf("%ld ", dr[l]);
     * }
     * printf("\n"); */
  }
  /* printf("k %d\n", k); */
  if (k == 0) {
    return NULL;
  }

  /* dense row is not reduced to zero, thus generate new sparse
   * pivot row and normalize it */
  val_t *row  = (val_t *)malloc(
      (unsigned long)(2*(ncols-dpiv)+2) * sizeof(val_t));
  j = 2;
  for (i = dpiv; i < ncols; ++i) {
    if (dr[i] != 0) {
      dr[i] = dr[i] % mod;
      if (dr[i] != 0) {
        row[j++]  = (val_t)i;
        row[j++]  = (val_t)dr[i];
      }
    }
  }
  row[0]  = j;
  row[1]  = 2 * (((j-2)/2) % UNROLL) + 2;
  row     = realloc(row, (unsigned long)row[0] * sizeof(val_t));

  /* for (int32_t l = 0; l < row[0]; ++l) {
   *   printf("%2d ", row[l]);
   * }
   * printf("\n"); */

  if (row[3] != 1) {
    normalize_matrix_row(row);
  }

  /* for (int32_t l = 0; l < row[0]; ++l) {
   *   printf("%2d ", row[l]);
   * }
   * printf("\n"); */

  return row;
}

static val_t *reduce_dense_row_by_known_pivots_32_bit(
    int64_t *dr,
    val_t *const *pivs,
    const val_t dpiv    /* pivot of dense row at the beginning */
    )
{
  int32_t i, j, k;
  const int64_t mod2  = (int64_t)fc * fc;

  for (k = 0, i = dpiv; i < ncols; ++i) {
    if (dr[i] != 0) {
      dr[i] = dr[i] % fc;
    }
    if (dr[i] == 0) {
      continue;
    }
    if (pivs[i] == NULL) {
      k++;
      continue;
    }

    /* found reducer row, get multiplier */
    const int64_t mul = (int64_t)dr[i];

    for (j = 2; j < pivs[i][1]; j += 2) {
      dr[pivs[i][j]]    -=  mul * pivs[i][j+1];
      dr[pivs[i][j]]    +=  (dr[pivs[i][j]] >> 63) & mod2;
    }
    for (; j < pivs[i][0]; j += 8) {
      dr[pivs[i][j]]    -=  mul * pivs[i][j+1];
      dr[pivs[i][j]]    +=  (dr[pivs[i][j]] >> 63) & mod2;
      dr[pivs[i][j+2]]  -=  mul * pivs[i][j+3];
      dr[pivs[i][j+2]]  +=  (dr[pivs[i][j+2]] >> 63) & mod2;
      dr[pivs[i][j+4]]  -=  mul * pivs[i][j+5];
      dr[pivs[i][j+4]]  +=  (dr[pivs[i][j+4]] >> 63) & mod2;
      dr[pivs[i][j+6]]  -=  mul * pivs[i][j+7];
      dr[pivs[i][j+6]]  +=  (dr[pivs[i][j+6]] >> 63) & mod2;
    }
    dr[i] = 0;
  }
  if (k == 0) {
    return NULL;
  }

  /* dense row is not reduced to zero, thus generate new sparse
   * pivot row and normalize it */
  val_t *row  = (val_t *)malloc(
      (unsigned long)(2*(ncols-dpiv)+2) * sizeof(val_t));
  j = 2;
  for (i = dpiv; i < ncols; ++i) {
    if (dr[i] != 0) {
      dr[i] = dr[i] % fc;
      if (dr[i] != 0) {
        row[j++]  = (val_t)i;
        row[j++]  = (val_t)dr[i];
      }
    }
  }
  row[0]  = j;
  row[1]  = 2 * (((j-2)/2) % UNROLL) + 2;
  row     = realloc(row, (unsigned long)row[0] * sizeof(val_t));

  if (row[3] != 1) {
    normalize_matrix_row(row);
  }

  return row;
}

static val_t **sparse_linear_algebra(
    val_t **mat
    )
{
  int32_t i, j;
  val_t sc  = 0;  /* starting column */
  val_t *npiv;    /* new pivot row */

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* all pivots, first we can only fill in all known lead terms */
  val_t **pivs  = (val_t **)calloc((unsigned long)ncols, sizeof(val_t *));
  /* unkown pivot rows we have to reduce with the known pivots first */
  val_t **upivs = (val_t **)malloc((unsigned long)nrl * sizeof(val_t *));

  /* for (i = 0; i < ncols; ++i) {
   *   printf("3%d %p\n", i, pivs[i]);
   * } */
  i = 0;
  j = 1;
  for (i = 0; i < nrows; ++i) {
    if (!pivs[mat[i][2]]) {
      pivs[mat[i][2]] = mat[i];
    } else {
      /* shorter rows first */
      upivs[nrl-j]  = mat[i];
      j++;
    }
  }

  free(mat);
  mat = NULL;
  printf("nthrds %d\n", nthrds);
#pragma omp parallel num_threads(nthrds)
  {
  int64_t *dr = (int64_t *)malloc((unsigned long)ncols * sizeof(int64_t));
#pragma omp parallel for num_threads(nthrds)
  for (i = 0; i < nrl; ++i) {
    printf("num of threads running %d\n", omp_get_num_threads());
    /* load next row to dense format for further reduction */
    memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
    for (j = 2; j < upivs[i][0]; j += 2) {
      dr[upivs[i][j]] = (int64_t)upivs[i][j+1];
    }
    sc  = upivs[i][2];
    free(upivs[i]);
    upivs[i]  = NULL;
    /* do the reduction */
    do {
      npiv  = reduce_dense_row_by_known_pivots(dr, pivs, sc);
      if (!npiv) {
        break;
      }
      /* j = compare_and_swap((void *)(&pivs[npiv[2]]), 0, (long)npiv); */
      j  = __sync_bool_compare_and_swap(&pivs[npiv[2]], NULL, npiv);
    } while (j);
  }
  free(dr);
  }
  free(upivs);
  upivs = NULL;

  /* we do not need the old pivots anymore */
  int64_t *dr = (int64_t *)malloc((unsigned long)ncols * sizeof(int64_t));
  for (i = 0; i < nru; ++i) {
    free(pivs[i]);
    pivs[i] = NULL;
  }

  npivs = 0; /* number of new pivots */

  mat = realloc(mat, (unsigned long)(ncr) * sizeof(val_t *));
  /* interreduce new pivots, i.e. pivs[ncl + ...] */
  for (i = (ncols-1); i >= ncl; --i) {
    if (pivs[i]) {
      memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
      for (j = 2; j < pivs[i][0]; j += 2) {
        dr[pivs[i][j]] = (int64_t)pivs[i][j+1];
      }
      sc  = pivs[i][2];
      free(pivs[i]);
      pivs[i] = NULL;
      mat[npivs++] = reduce_dense_row_by_known_pivots(dr, pivs, sc);
    }
  }
  free(dr);
  free(pivs);
  pivs  = NULL;

  mat   = realloc(mat, (unsigned long)npivs * sizeof(val_t *));
  nrows = nrall = npivs;

  num_zerored += (nrl - npivs);
  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  la_ctime  +=  ct1 - ct0;
  la_rtime  +=  rt1 - rt0;

  GB_DEBUG(LADBG, "%7d new %7d zero - %9.3f sec", npivs, nrl, rt1-rt0);

  return mat;
}

static val_t **probabilistic_sparse_linear_algebra(
    val_t **mat
    )
{
  int32_t i, j, k;
  val_t sc  = 0;  /* starting column */
  val_t *npiv;    /* new pivot row */

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* all pivots, first we can only fill in all known lead terms */
  val_t **pivs  = (val_t **)calloc((unsigned long)ncols, sizeof(val_t *));
  /* unkown pivot rows we have to reduce with the known pivots first */
  val_t **upivs = (val_t **)malloc((unsigned long)nrl * sizeof(val_t *));

  i = 0;
  j = 1;
  for (i = 0; i < nrows; ++i) {
    if (!pivs[mat[i][2]]) {
      pivs[mat[i][2]] = mat[i];
    } else {
      /* shorter rows first */
      upivs[nrl-j]  = mat[i];
      j++;
    }
  }

  free(mat);
  mat = NULL;

  /* compute rows per block */
  const int32_t nb  = (int32_t)(floor(sqrt(nrl/2))) > 0 ?
    (int32_t)(floor(sqrt(nrl/2))) :
    (int32_t)(floor(sqrt(nrl))) ;
  const int32_t rem = (nrl % nb == 0) ? 0 : 1;
  const int32_t rpb = (nrl / nb) + rem;

  /* we need to allocate dr statically for openmp */
  int64_t *dr   = (int64_t *)malloc((unsigned long)ncols * sizeof(int64_t));
  int64_t *mul  = (int64_t *)malloc((unsigned long)rpb * sizeof(int64_t));

/* #pragma omp parallel for num_threads(nthrds) shared(drg) \
 *   private(j, k, dr, mul) */
  for (i = 0; i < nb; ++i) {
    const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
    const int32_t nrbl  = (int32_t) (nbl - i*rpb);

    sc  = 0;

    if (nrbl != 0) {
      int32_t bctr  = 0;
      while (bctr < nrbl) {
        /* fill random value array */
        for (j = 0; j < nrbl; ++j) {
          mul[j]  = (int64_t)rand() % fc;
        }

        /* generate one dense row as random linear combination
        * of the rows of the block */
        memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
        int32_t ctr = 0;
        for (j = i*rpb; j < nbl; ++j) {
          sc  = sc < upivs[j][2] ? sc : upivs[j][2];
          for (k = 2; k < upivs[j][0]; k = k+2) {
            dr[upivs[j][k]]  +=  mul[ctr] * upivs[j][k+1];
          }
          ctr++;
        }

        /* do the reduction */
        do {
          npiv  = reduce_dense_row_by_known_pivots(dr, pivs, sc);
          if (!npiv) {
            bctr  = nrbl;
            break;
          }
          j = compare_and_swap((void *)(&pivs[npiv[2]]), 0, (long)npiv);
        } while (j);
        bctr++;
      }
      /* fill global dense row for final check at the end */
      /* for (j = 0; j < nrbl; ++j) {
       *   mul[j]  = (val_t)rand() % fc;
       * }
       * for (j = i*rpb; j < nbl; ++j) {
       *   for (k = 2; k < upivs[j][0]; k = k+2) {
       *     if (upivs[j][k] >= ncols) {
       *       printf("upivs[%d][%d] = %d > %d\n", j, k, upivs[j][k], ncols);
       *     }
       *     drg[upivs[j][k]]  +=  mul[j] * upivs[j][k+1];
       *   }
       * }
       * printf("drg2\n");
       * for (i = 0; i < ncols; ++i) {
       *   printf("%ld ", drg[i]);
       * }
       * printf("\n"); */
    }
  }
  /* do final check, go back if check fails */

  for (i = 0; i < nrl; ++i) {
    free(upivs[i]);
    upivs[i] = NULL;
  }
  free(upivs);
  upivs = NULL;

  /* we do not need the old pivots anymore */
  for (i = 0; i < nru; ++i) {
    free(pivs[i]);
    pivs[i] = NULL;
  }

  npivs = 0; /* number of new pivots */

  mat = realloc(mat, (unsigned long)(ncr) * sizeof(val_t *));
  /* interreduce new pivots, i.e. pivs[ncl + ...] */
  for (i = (ncols-1); i >= ncl; --i) {
    if (pivs[i]) {
      memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
      for (j = 2; j < pivs[i][0]; j += 2) {
        dr[pivs[i][j]] = (int64_t)pivs[i][j+1];
      }
      sc  = pivs[i][2];
      free(pivs[i]);
      pivs[i] = NULL;
      mat[npivs++] = reduce_dense_row_by_known_pivots(dr, pivs, sc);
    }
  }
  free(dr);
  free(mul);
  free(pivs);
  pivs  = NULL;

  mat   = realloc(mat, (unsigned long)npivs * sizeof(val_t *));
  nrows = nrall = npivs;

  num_zerored += (nrl - npivs);
  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  la_ctime  +=  ct1 - ct0;
  la_rtime  +=  rt1 - rt0;

  GB_DEBUG(LADBG, "%7d new %7d zero - %9.3f sec", npivs, nrl, rt1-rt0);

  return mat;
}
