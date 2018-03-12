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
 * \file convert.c
 * \brief Implementation of hash value / column index conversions
 * for the exponent hashes of polynomials resp. rows.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

/* after calling this procedure we have column indices instead of exponent
 * hashes in the polynomials resp. rows. moreover, we have sorted each row
 * by pivots / non-pivots. thus we get already an A|B splicing of the
 * initial matrix. this is a first step for receiving a full GBLA matrix. */
static len_t *convert_hashes_to_columns(
    val_t **mat
    )
{
  int32_t i, j, k, l;
  val_t *row;
  len_t *hcm; /* hash-to-column map */
  int64_t nterms = 0;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* need to allocate memory for all possible exponents we
   * have in the local hash table since we do not which of
   * them are corresponding to multipliers and which are
   * corresponding to the multiplied terms in reducers. */
  hcm = (len_t *)malloc((unsigned long)ncols * sizeof(len_t));
  /* j counts all columns, k counts known pivots */
  for (j = 0, k = 0, i = 0; i < nrows; ++i) {
    row     =   mat[i];
    nterms  +=  (int64_t)row[0];
    for (l = 2; l < row[0]; l += 2) {
      if ((ev+row[l])[HASH_IND] > 0) {
        hcm[j++]  = row[l];
        if ((ev+row[l])[HASH_IND] == 2) {
          k++;
          (ev+row[l])[HASH_IND]  = -1;
        } else {
          (ev+row[l])[HASH_IND]  = -2;
        }
      }
    }
  }
  /* sort monomials w.r.t known pivots, then w.r.t. to the monomial order */
  qsort(hcm, (unsigned long)j, sizeof(len_t), hcm_cmp);
  
  /* set number of rows and columns in ABCD splicing */
  nru = ncl = k;
  nrl = nrows - nru;
  ncr = j - ncl;

  /* store the other direction (hash -> column) in HASH_IND */
  for (i = 0; i < j; ++i) {
    (ev + hcm[i])[HASH_IND]  = i;
  }


  
  /* map column positions to matrix rows */
#pragma omp parallel for num_threads(nthrds) private(i, j)
  for (i = 0; i < nrows; ++i) {
    row = mat[i];
    for (j = 2; j < row[0]; j += 2) {
      row[j]  = (ev + row[j])[HASH_IND];
    }
  }

  /* next we sort each row by the new colum order due
   * to known / unkown pivots */
#pragma omp parallel for num_threads(nthrds) private(i)
  for (i = 0; i < nrows; ++i) {
    qsort(mat[i]+2, (unsigned long)(mat[i][0]-2)/2, 2 * sizeof(val_t),
        columns_cmp);
  }
  /* compute density of matrix */
  nterms  *=  100; /* for percentage */
  density = (double)nterms / (double)nrows / (double)ncols;

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  convert_ctime +=  ct1 - ct0;
  convert_rtime +=  rt1 - rt0;
  GB_DEBUG(SYMDBG, " %7d x %7d mat - %6.3f%%", nrows, ncols, density);

  return hcm;
}

static val_t **convert_columns_to_hashes(
    val_t **mat,
    const len_t *hcm
    )
{
  int32_t i, j;
  val_t *row;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  for (i = 0; i < ncols; ++i) {
    (ev+hcm[i])[HASH_IND] = 0;
  }

#pragma omp parallel for num_threads(nthrds) private(i, j)
  for (i = 0; i < nrows; ++i) {
    row = mat[i];
    for (j = 2; j < row[1]; j += 2) {
      row[j]  = hcm[row[j]];
    }
    /* loop unrolling, UNROLL = 4 */
    for (; j < row[0]; j += 8) {
      row[j]    = hcm[row[j]];
      row[j+2]  = hcm[row[j+2]];
      row[j+4]  = hcm[row[j+4]];
      row[j+6]  = hcm[row[j+6]];
    }
  }

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  convert_ctime +=  ct1 - ct0;
  convert_rtime +=  rt1 - rt0;

  return mat;
}
