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

static val_t **select_spairs_by_deg_lex(
    void
    )
{
  int32_t i, j, k, md, npairs;
  val_t *b;
  len_t m;

  len_t lcm = 0, load = 0, load_old = 0;
  val_t **mat;
  len_t *gens;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* sort pair set */
  qsort(ps, (unsigned long)pload, SP_LEN * sizeof(sp_t), &spair_lex_cmp);
  /* get minimal degree */
  md  = ps[SP_DEG];
  /* printf("\n pairs sorted:\n");
   * for (i = 0; i < pload; ++i) {
   *   printf("p%2d | %d | %d | %d || ", i, ps[i].gen1, ps[i].gen2, ps[i].deg);
   *   for (j = 0; j < nvars; ++j) {
   *     printf("%d ", (ev+ps[i].lcm)[j]);
   *   }
   *   printf("\n\n");
   * } */

  /* select pairs of this degree respecting maximal selection size mnsel */
  for (i = 0; i < pload; i = i + SP_LEN) {
    if ((ps+i)[SP_DEG] > md || i >= mnsel*SP_GEN) {
      break;
    }
  }
  npairs  = i;
  /* if we stopped due to maximal selection size we still get the following
   * pairs of the same lcm in this matrix */
  if (i > mnsel*SP_GEN) {
    j = i + SP_LEN;
    while ((ps+j)[SP_LCM] == (ps+i)[SP_LCM]) {
      j = j + SP_LEN;
    }
    npairs = j;
  }
  npairs  = npairs/SP_LEN;
  GB_DEBUG(SELDBG, " %6d/%6d pairs - deg %2d", npairs, pload/SP_LEN, md);
  /* statistics */
  num_pairsred  +=  npairs;
  
  gens  = (len_t *)malloc(2 * (unsigned long)npairs * sizeof(len_t));

  /* preset matrix meta data */
  mat   = (val_t **)malloc(2 * (unsigned long)npairs * sizeof(val_t *));
  nrall = 2 * npairs;
  ncols = ncl = ncr = 0;
  nrows = 0;

  i = 0;
  npairs  = npairs * SP_LEN;
  while (i < npairs) {
    load_old  = load;
    lcm   = (ps+i)[SP_LCM];
    gens[load++] = (ps+i)[SP_G1];
    gens[load++] = (ps+i)[SP_G2];
    j = i + SP_LEN;
    while (j < npairs && (ps+j)[SP_LCM] == lcm) {
      for (k = load_old; k < load; ++k) {
        if ((ps+j)[SP_G1]== gens[k]) {
          break;
        }
      }
      if (k == load) {
        gens[load++]  = (ps+j)[SP_G1];
      }
      for (k = load_old; k < load; ++k) {
        if ((ps+j)[SP_G2] == gens[k]) {
          break;
        }
      }
      if (k == load) {
        gens[load++]  = (ps+j)[SP_G2];
      }
      j = j + SP_LEN;
    }
    for (k = load_old; k < load; ++k) {
      b = (val_t *)((long)bs[gens[k]] & bmask);
      m = monomial_division_no_check(lcm, b[2]);
      mat[nrows++]  = multiplied_polynomial_to_matrix_row(m, b);
    }

    i = j;
  }

/*   val_t *tmp  = (val_t *)calloc(100000, sizeof(val_t));
 *   int32_t ctr = 0;
 *   for (i=0; i < nrows; ++i) {
 *     for (j = 2; j < mat[i][0]; j += 2) {
 *       k = 0;
 *       while (tmp[k] != 0) {
 *         if (tmp[k] == mat[i][j]) {
 *           break;
 *         }
 *         k++;
 *       }
 *       if (k == ctr) {
 *         tmp[ctr]  = mat[i][j];
 *         ctr++;
 *       }
 *     }
 *   }
 *   printf("\npair poly rows (%d) have %d columns together\n", nrows, ctr);
 *   free(tmp);
 *   tmp = NULL;
 *  */

  num_duplicates  +=  0;
  num_rowsred     +=  load;

  free(gens);

  /* remove selected spairs from pairset */
  for (j = npairs; j < pload; j = j + SP_LEN) {
    memcpy(ps+j-npairs, ps+j, (unsigned long)SP_LEN * sizeof(sp_t));
  }
  pload = pload - npairs;

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  select_ctime  +=  ct1 - ct0;
  select_rtime  +=  rt1 - rt0;
  
  return mat;
}

static val_t **select_spairs_by_minimal_degree(
    void
    )
{
  int32_t i, j, k, md, npairs;
  val_t *b;
  len_t m;

  len_t lcm = 0, load = 0, load_old = 0;
  val_t **mat;
  len_t *gens;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  /* for (i = 0; i < pload; i = i + SP_LEN) {
   *   printf("pair[%d] = %d | %d | %d || %p\n", i/SP_LEN, (ps+i)[SP_G1], (ps+i)[SP_G2], (ps+i)[SP_DEG], ps+i);
   * } */
  /* sort pair set */
  qsort(ps, (unsigned long)pload/SP_LEN, SP_LEN * sizeof(sp_t), &spair_cmp);
  /* for (i = 0; i < pload; i = i + SP_LEN) {
   *   printf("pair[%d] = %d | %d | %d || %p\n", i/SP_LEN, (ps+i)[SP_G1], (ps+i)[SP_G2], (ps+i)[SP_DEG], ps+i);
   * } */
  /* get minimal degree */
  md  = ps[SP_DEG];
  /* printf("\n pairs sorted:\n");
   * for (i = 0; i < pload; ++i) {
   *   printf("p%2d | %d | %d | %d || ", i, ps[i].gen1, ps[i].gen2, ps[i].deg);
   *   for (j = 0; j < nvars; ++j) {
   *     printf("%d ", (ev+ps[i].lcm)[j]);
   *   }
   *   printf("\n\n");
   * } */

  /* select pairs of this degree respecting maximal selection size mnsel */
  for (i = 0; i < pload; i = i + SP_LEN) {
    if ((ps+i)[SP_DEG] > md || i >= mnsel) {
      break;
    }
  }
  npairs  = i;
  /* if we stopped due to maximal selection size we still get the following
   * pairs of the same lcm in this matrix */
  if (i > mnsel) {
    j = i + SP_LEN;
    while ((ps+j)[SP_LCM] == (ps+i)[SP_LCM]) {
      j = j + SP_LEN;
    }
    npairs = j;
  }
  npairs  = npairs/SP_LEN;
  GB_DEBUG(SELDBG, " %6d/%6d pairs - deg %2d", npairs, pload/SP_LEN, md);
  /* statistics */
  num_pairsred  +=  npairs;
  
  gens  = (len_t *)malloc(2 * (unsigned long)npairs * sizeof(len_t));

  /* preset matrix meta data */
  mat   = (val_t **)malloc(2 * (unsigned long)npairs * sizeof(val_t *));
  nrall = 2 * npairs;
  ncols = ncl = ncr = 0;
  nrows = 0;

  i = 0;
  npairs  = npairs * SP_LEN;
  while (i < npairs) {
    load_old  = load;
    lcm   = (ps+i)[SP_LCM];
    gens[load++] = (ps+i)[SP_G1];
    gens[load++] = (ps+i)[SP_G2];
    j = i + SP_LEN;
    while (j < npairs && (ps+j)[SP_LCM] == lcm) {
      for (k = load_old; k < load; ++k) {
        if ((ps+j)[SP_G1]== gens[k]) {
          break;
        }
      }
      if (k == load) {
        gens[load++]  = (ps+j)[SP_G1];
      }
      for (k = load_old; k < load; ++k) {
        if ((ps+j)[SP_G2] == gens[k]) {
          break;
        }
      }
      if (k == load) {
        gens[load++]  = (ps+j)[SP_G2];
      }
      j = j + SP_LEN;
    }
    for (k = load_old; k < load; ++k) {
      b = (val_t *)((long)bs[gens[k]] & bmask);
      m = monomial_division_no_check(lcm, b[2]);
      mat[nrows++]  = multiplied_polynomial_to_matrix_row(m, b);
    }

    i = j;
  }

/*   val_t *tmp  = (val_t *)calloc(100000, sizeof(val_t));
 *   int32_t ctr = 0;
 *   for (i=0; i < nrows; ++i) {
 *     for (j = 2; j < mat[i][0]; j += 2) {
 *       k = 0;
 *       while (tmp[k] != 0) {
 *         if (tmp[k] == mat[i][j]) {
 *           break;
 *         }
 *         k++;
 *       }
 *       if (k == ctr) {
 *         tmp[ctr]  = mat[i][j];
 *         ctr++;
 *       }
 *     }
 *   }
 *   printf("\npair poly rows (%d) have %d columns together\n", nrows, ctr);
 *   free(tmp);
 *   tmp = NULL;
 *  */

  num_duplicates  +=  0;
  num_rowsred     +=  load;

  free(gens);

  /* for (i = 0; i < pload; i = i + SP_LEN) {
   *   printf("3pair[%d] = %d | %d | %d || %p\n", i/SP_LEN, (ps+i)[SP_G1], (ps+i)[SP_G2], (ps+i)[SP_DEG], ps+i);
   * } */
  /* remove selected spairs from pairset */
  for (j = npairs; j < pload; j = j + SP_LEN) {
    memcpy(ps+j-npairs, ps+j, (unsigned long)SP_LEN * sizeof(sp_t));
  }
  pload = pload - npairs;
  /* for (i = 0; i < pload; i = i + SP_LEN) {
   *   printf("4pair[%d] = %d | %d | %d || %p\n", i/SP_LEN, (ps+i)[SP_G1], (ps+i)[SP_G2], (ps+i)[SP_DEG], ps+i);
   * } */

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

  const int32_t bl  = bload;
  /* printf("to be divided: ");
   * for (int32_t j = 0; j < nvars; ++j) {
   *   printf("%d",(ev+m)[j]);
   * }
   * printf("\n"); */

  /* search basis */
  for (i = (ev+m)[HASH_DIV]; i < bl; ++i) {
    b = (val_t *)((long)bs[i] & bmask);
    /* if we only use nonredundant basis elements LEX computation are very slow */
    /* if (b != bs[i]) {
     *   continue;
     * } */
    /* printf("test divisor: ");
     * for (int32_t j = 0; j < nvars; ++j) {
     *   printf("%d", (ev+b[2])[j]);
     * } */
    d = monomial_division_with_check(m, b[2]);
    if (d == 0) {
      continue;
    }
    (ev+m)[HASH_DIV] = i;
    /* printf("\ndivisor found: %d   ", b[2]);
     * for (int32_t j = 0; j < nvars; ++j) {
     *   printf("%d", (ev+b[2])[j]);
     * }
     * printf("\n"); */
    return multiplied_polynomial_to_matrix_row(d, b);
  }
  (ev+m)[HASH_DIV] = i;
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

  /* mark leading monomials in HASH_IND entry */
  for (i = 0; i < nrows; ++i) {
    if (!(ev+mat[i][2])[HASH_IND]) {
      (ev+mat[i][2])[HASH_IND] = 2;
      ncols++;
    }
  }

  /* get reducers from basis */
  for (i = 0; i < nrows; ++i) {
    const len_t len = mat[i][0];
    for (j = 4; j < len; j += 2) {
      m = mat[i][j];
      if ((ev+m)[HASH_IND]) {
        continue;
      }
      ncols++;
      red = find_multiplied_reducer(m);
      if (!red) {
        (ev+m)[HASH_IND] = 1;
        continue;
      }
      (ev+m)[HASH_IND] = 2;
      if (nrows == nrall) {
        nrall = 2 * nrall;
        mat   = realloc(mat, (unsigned long)nrall * sizeof(val_t *));
      }
      /* add new reducer to matrix */
      mat[nrows++]  = red;
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
