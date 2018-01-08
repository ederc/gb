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

  const len_t len   = row[0];
  const int64_t inv = (int64_t)mod_p_inverse_32(row[3], fc);
  
  i = (len-2)/2;
  i = i & 1 ? 5 : 3;
  for (; i < len; i = i+4) {
    row[i]    = (val_t)(row[i] * inv) % fc;
    row[i+2]  = (val_t)(row[i+2] * inv) % fc;
  }
  /* probably we left out the first coefficient */
  row[3]  = 1;
}

static val_t *reduce_dense_row_by_known_pivots_16_bit(
    int64_t *dr,
    val_t **pivs
    )
{
}

static val_t **sparse_linear_algebra_16_bit(
    val_t **mat
    )
{
  int32_t i, j;
  double ct0, ct1, rt0, rt1;
  val_t *npiv; /* new pivot row */

  /* timings */
  ct0 = cputime();
  rt0 = realtime();

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  la_ctime  +=  ct1 - ct0;
  la_rtime  +=  rt1 - rt0;

  /* all pivots, first we can only fill in all known lead terms */
  val_t **pivs  = (val_t **)calloc((unsigned long)nrows, sizeof(val_t *));
  /* unkown pivot rows we have to reduce with the known pivots first */
  val_t **upivs = (val_t **)malloc((unsigned long)nrl * sizeof(val_t *));

  i = 0;
  j = 1;
  while (i < nrows) {
    if (!pivs[mat[i][2]]) {
      pivs[mat[i][2]] = mat[i];
    } else {
      /* shorter rows first */
      pivs[nru-j] = mat[i];
      j++;
    }
  }

  free(mat);
  mat = NULL;
  const unsigned long nc  = (unsigned long)(ncl + ncr);
  int64_t *dr = (int64_t *)malloc(nc * sizeof(int64_t));
#pragma omp parallel for num_threads(nthrds) private(dr)
  for (i = 0; i < nrl; ++i) {
    /* load next row to dense format for further reduction */
    memset(dr, 0, nc * sizeof(int64_t));
    for (j = 2; j < upivs[i][0]; j += 2) {
      dr[upivs[i][j]] = (int64_t)upivs[i][j+1];
    }
    free(upivs[i]);
    upivs[i]  = NULL;
    /* do the reduction */
    do {
      npiv  = reduce_dense_row_by_known_pivots_16_bit(dr, pivs);
      if (!npiv) {
        break;
      }
      j = __sync_bool_comparse_and_swap(&pivs[npiv[2]], NULL, npiv);
    } while (j);
  }

  /* we do not need the old pivots anymore */
  for (i = 0; i < nru; ++i) {
    free(pivs[i]);
    pivs[i] = NULL;
  }

  npivs = 0; /* number of new pivots */

  /* interreduce new pivots, i.e. pivs[ncl + ...] */
  for (i = (ncl+ncr-1); i > ncl+1; ++i) {
    if (pivs[i]) {
      memset(dr, 0, nc * sizeof(int64_t));
      for (j = 2; j < pivs[i][0]; j += 2) {
        dr[pivs[i][j]] = (int64_t)pivs[i][j+1];
      }
      free(pivs[i]);
      pivs[i] = NULL;
      mat[npivs++] = reduce_dense_row_by_known_pivots_16_bit(dr, pivs);
    }
  }
  mat   = realloc(mat, (unsigned long)npivs * sizeof(val_t *));
  nrows = nrall = npivs;

  free(dr);

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  la_ctime  +=  ct1 - ct0;
  la_rtime  +=  rt1 - rt0;

  return mat;
}
