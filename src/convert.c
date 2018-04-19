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
    mat_t *mat,
    md_t *md
    )
{
  int32_t i, j, k, l;
  row_t *row;
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
  hcm = (len_t *)malloc((unsigned long)mat->nc * sizeof(len_t));
  /* j counts all columns, k counts known pivots */
  for (j = 0, k = 0, i = 0; i < mat->nr; ++i) {
    row     =   mat->r[i];
    nterms  +=  (int64_t)row->sz;
    for (l = 0; l < row->sz; ++l) {
#if ORDER_COLUMNS
      if ((ev+row->ch[l])[HASH_IND] > 0) {
#else
      if ((ev+row->ch[l])[HASH_IND] != 0) {
#endif
        hcm[j++]  = row->ch[l];
        if ((ev+row->ch[l])[HASH_IND] == 2) {
          k++;
#if ORDER_COLUMNS
          (ev+row->ch[l])[HASH_IND]  = -1;
        } else {
          (ev+row->ch[l])[HASH_IND]  = -2;
        }
#else
        }
        (ev+row->ch[l])[HASH_IND] = 0;
#endif
      }
    }
  }
  /* sort monomials w.r.t known pivots, then w.r.t. to the monomial order */
  qsort(hcm, (unsigned long)j, sizeof(len_t), hcm_cmp);
  
  /* set number of rows and columns in ABCD splicing */
  mat->nru  = mat->ncl  = k;
  mat->nrl  = mat->nr - mat->nru;
  mat->ncr  = j - mat->ncl;

  /* store the other direction (hash -> column) in HASH_IND */
  for (i = 0; i < j; ++i) {
    (ev + hcm[i])[HASH_IND]  = i;
  }


  
  /* map column positions to matrix rows */
  printf("mat->nr %d\n", mat->nr);
#pragma omp parallel for num_threads(md->nt) private(i, j)
  for (i = 0; i < mat->nr; ++i) {
    row = mat->r[i];
    printf("i %d -- %p\n", i, row);
    for (j = 0; j < row->os; ++j) {
      row->ch[j]  = (ev + row->ch[j])[HASH_IND];
    }
    for (; j < row->sz; j += 4) {
      row->ch[j]    = (ev + row->ch[j])[HASH_IND];
      row->ch[j+1]  = (ev + row->ch[j+1])[HASH_IND];
      row->ch[j+2]  = (ev + row->ch[j+2])[HASH_IND];
      row->ch[j+3]  = (ev + row->ch[j+3])[HASH_IND];
    }
  }

  /* next we sort each row by the new colum order due
   * to known / unkown pivots */

  /* NOTE: As strange as it may sound, we do not need to sort the rows.
   * When reducing, we copy them to dense rows, there we copy the coefficients
   * at the right place and reduce then. For the reducers itself it is not
   * important in which order the terms are represented as long as the first
   * term is the lead term, which is always true. Once a row is finally reduced
   * it is copied back to a sparse representation, now in the correct term
   * order since it is coming from the correctly sorted dense row. So all newly
   * added elements have all their terms sorted correctly w.r.t. the given
   * monomial order. */

  /* compute density of matrix */
  nterms  *=  100; /* for percentage */
  density = (double)nterms / (double)mat->nr / (double)mat->nc;

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->convert_ctime +=  ct1 - ct0;
  md->convert_rtime +=  rt1 - rt0;
  GB_DEBUG(SYMDBG, " %7d x %7d mat - %6.3f%%", mat->nr, mat->nc, density);

  return hcm;
}

static void convert_columns_to_hashes(
    mat_t *mat,
    md_t *md,
    const len_t *hcm
    )
{
  int32_t i, j;
  row_t *row;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  for (i = 0; i < mat->nc; ++i) {
    (ev+hcm[i])[HASH_IND] = 0;
  }

#pragma omp parallel for num_threads(md->nt) private(i, j)
  for (i = 0; i < mat->nr; ++i) {
    row = mat->r[i];
    for (j = 0; j < row->os; ++j) {
      row->ch[j]  = hcm[row->ch[j]];
    }
    /* loop unrolling, UNROLL = 4 */
    for (; j < row->sz; j += 4) {
      row->ch[j]    = hcm[row->ch[j]];
      row->ch[j+1]  = hcm[row->ch[j+1]];
      row->ch[j+2]  = hcm[row->ch[j+2]];
      row->ch[j+3]  = hcm[row->ch[j+3]];
    }
  }

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->convert_ctime +=  ct1 - ct0;
  md->convert_rtime +=  rt1 - rt0;
}
