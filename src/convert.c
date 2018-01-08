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
  int32_t i, j, k;
  val_t *row;
  exp_t *e;
  len_t *hcm; /* hash-to-column map */
  uint64_t nterms = 0;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* need to allocate memory for all possible exponents we
   * have in the local hash table since we do not which of
   * them are corresponding to multipliers and which are
   * corresponding to the multiplied terms in reducers. */
  hcm = (len_t *)malloc((unsigned long)(elload/HASH_LEN) * sizeof(len_t));

  for (i = 0, k = 0; i < nrows; ++i) {
    row     =   mat[i];
    nterms  +=  (uint64_t)row[0];
    /* first monomial is special: gives us all known pivots */
    (evl + row[2])[HASH_IND]  = 2;
    for (j = 4; j < row[0]; j += 2) {
      e = evl + row[j];
      if (e[HASH_IND]) {
        continue;
      }
      e[HASH_IND] = 1;
      hcm[k++]    = row[j];
    }
  }
  hcm = realloc(hcm, (unsigned long)k * sizeof(len_t));

  /* sort monomials w.r.t known pivots, then w.r.t. to the monomial order */
  qsort(hcm, (unsigned long)k, sizeof(len_t), hcm_cmp);
  
  /* get number of known pivots, i.e. number of upper rows and number of
   * left columns in ABCD splicing of GBLA matrix. moreover, from this
   * information we can also compute the number of lower rows and the
   * number of right columns */
  j = 0;
  while ((evl+hcm[j])[HASH_IND] == 2) {
    j++;
  }
  /* set number of rows and columns in ABCD splicing */
  nru = ncl = j;
  nrl = nrows - nru;
  ncr = k - ncl;
  nc  = ncl + ncr;

  /* store the other direction (hash -> column) in HASH_IND */
  for (i = 0; i < k; ++i) {
    (evl + hcm[i])[HASH_IND]  = i;
  }
  
  /* map column positions to matrix rows */
  for (i = 0; i < nrows; ++i) {
    row = mat[i];
    for (j = 2; j < row[0]; j += 2) {
      row[j]  = (evl + row[j])[HASH_IND];
    }
  }

  /* next we sort each row by the new colum order due
   * to known / unkown pivots */
#pragma omp parallel for num_threads(nthrds)
  for (int32_t l = 0; l < nrows; ++l) {
    qsort(mat[l]+2, (unsigned long)(mat[l][0]-2)/2, 2 * sizeof(val_t),
        columns_cmp);
  }
  /* compute density of matrix */
  density = (double)nterms / (double)(nrows * k);

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  convert_ctime +=  ct1 - ct0;
  convert_rtime +=  rt1 - rt0;

  return hcm;
}

static val_t **convert_columns_to_hashes(
    val_t **mat,
    len_t *hcm
    )
{
  int32_t i, j, k;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* TODO */

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  convert_ctime +=  ct1 - ct0;
  convert_rtime +=  rt1 - rt0;
}