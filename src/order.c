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
 * \file order.c
 * \brief Order procedures
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

#if 0
/* sorting stuff for matrices */
static int columns_cmp(
    const void *a,
    const void *b
    )
{
  const val_t ca  = *((val_t *)a);
  const val_t cb  = *((val_t *)b);

  return (int)(ca - cb);
}
#endif

static int matrix_row_initial_input_cmp_lex(
    const void *a,
    const void *b
    )
{
  int32_t i;

  const len_t va  = (*(row_t **)a)->ch[0];
  const len_t vb  = (*(row_t **)b)->ch[0];

  const exp_t * const ea  = ev + va;
  const exp_t * const eb  = ev + vb;

  /* lexicographical */
  for (i = 0; i < nvars; ++i) {
    if (ea[i] < eb[i]) {
      return -1;
    }
    if (ea[i] > eb[i]) {
      return 1;
    }
  }
  return 0;
  /* return memcmp(ea, eb, (unsigned long)nvars * sizeof(exp_t)); */
}

static int matrix_row_initial_input_cmp_drl(
    const void *a,
    const void *b
    )
{
  int32_t i;
  len_t va, vb;

  va  = (*(row_t *)a).ch[0];
  vb  = (*(row_t *)b).ch[0];

  const exp_t * const ea  = ev + va;
  const exp_t * const eb  = ev + vb;

  /* DRL */
  if (ea[HASH_DEG] < eb[HASH_DEG]) {
    return -1;
  } else {
    if (ea[HASH_DEG] != eb[HASH_DEG]) {
      return 1;
    }
  }

  /* note: reverse lexicographical */
  for (i = nvars-1; i >= 0; --i) {
    if (ea[i] < eb[i]) {
      return -1;
    } else {
      if (ea[i] != eb[i]) {
        return 1;
      }
    }
  }
  return 0;
}

static int matrix_row_cmp(
    const void *a,
    const void *b
    )
{
  len_t ca, cb;
  /* compare pivot resp. column index */
  ca  = ((row_t *)a)->ch[0];
  cb  = ((row_t *)b)->ch[0];
  if (ca > cb) {
    return 1;
  }
  if (ca < cb) {
    return -1;
  }
  len_t la, lb;
  /* same column index => compare density of row */
  la  = ((row_t *)a)->sz;
  lb  = ((row_t *)b)->sz;
  if (la > lb) {
    return 1;
  }
  if (la < lb) {
    return -1;
  }
  return 0;
}

static inline void sort_matrix_rows(
    mat_t *mat)
{
  qsort(mat->r, (unsigned long)mat->nr, sizeof(row_t *), &matrix_row_cmp);
}

/* comparison for monomials (in local hash table) */
static int monomial_cmp_pivots_drl(
    const len_t a,
    const len_t b
    )
{
  int32_t i;

  const exp_t * const ea  = ev + a;
  const exp_t * const eb  = ev + b;

#if ORDER_COLUMNS
  /* first known pivots vs. tail terms */
  if (ea[HASH_IND] != eb[HASH_IND]) {
    if (ea[HASH_IND] < eb[HASH_IND]) {
      return 1;
    } else {
      return -1;
    }
  } else {
#endif

    /* then DRL */
    if (ea[HASH_DEG] > eb[HASH_DEG]) {
      return -1;
    } else {
      if (ea[HASH_DEG] != eb[HASH_DEG]) {
        return 1;
      }
    }

    /* note: reverse lexicographical */
    for (i = nvars-1; i >= 0; --i) {
      if (ea[i] > eb[i]) {
        return 1;
      } else {
        if (ea[i] != eb[i]) {
          return -1;
        }
      }
    }
#if ORDER_COLUMNS
  }
#endif

  return 0;
}

static int monomial_cmp_pivots_lex(
    const len_t a,
    const len_t b
    )
{
  int32_t i;
  const exp_t * const ea  = ev + a;
  const exp_t * const eb  = ev + b;

  /* first known pivots vs. tail terms */
  if (ea[HASH_IND] < eb[HASH_IND]) {
    return 1;
  } else {
    if (ea[HASH_IND] != eb[HASH_IND]) {
      return -1;
    }
  }

  /* lexicographical */
  for (i = 0; i < nvars; ++i) {
    if (eb[i] < ea[i]) {
      return -1;
    }
    if (eb[i] > ea[i]) {
      return 1;
    }
  }
  return 0;
  /* return memcmp(eb, ea, (unsigned long)nvars * sizeof(exp_t)); */
}

static inline int monomial_cmp_drl(
    const exp_t * const ea,
    const exp_t * const eb
    )
{
  int32_t i;

  if (ea[HASH_DEG] > eb[HASH_DEG]) {
    return 1;
  } else {
    if (ea[HASH_DEG] != eb[HASH_DEG]) {
      return -1;
    }
  }

  for (i = nvars-1; i >= 0; --i) {
    if (ea[i] < eb[i]) {
      return 1;
    } else {
      if (ea[i] != eb[i]) {
        return -1;
      }
    }
  }
  return 0;
}

static inline int monomial_cmp_lex(
    const exp_t * const ea,
    const exp_t * const eb
    )
{
  int32_t i;

  for (i = 0; i < nvars; ++i) {
    if (ea[i] < eb[i]) {
      return -1;
    }
    if (ea[i] > eb[i]) {
      return 1;
    }
  }
  return 0;
  /* return memcmp(ea, eb, (unsigned long)nvars * sizeof(exp_t)); */
}

/* comparison for hash-column-maps */
static int hcm_cmp_pivots_drl(
    const void *a,
    const void *b
    )
{
  const len_t ma  = ((len_t *)a)[0];
  const len_t mb  = ((len_t *)b)[0];

  return monomial_cmp_pivots_drl(ma, mb);
}

static int hcm_cmp_pivots_lex(
    const void *a,
    const void *b
    )
{
  const len_t ma  = ((len_t *)a)[0];
  const len_t mb  = ((len_t *)b)[0];

  return monomial_cmp_pivots_lex(ma, mb);
}

/* comparison for s-pairs once their lcms are in the global hash table */
static int spair_cmp_deglex(
    const void *a,
    const void *b
    )
{
  spair_t *sa = (spair_t *)a;
  spair_t *sb = (spair_t *)b;
  if (sa->deg != sb->deg) {
    return (sa->deg < sb->deg) ? -1 : 1;
  } else {
    return (int)monomial_cmp(ev+sa->lcm, ev+sb->lcm);
  }
}

static int spair_cmp_drl(
    const void *a,
    const void *b
    )
{
  spair_t *sa = (spair_t *)a;
  spair_t *sb = (spair_t *)b;
  return (int)monomial_cmp(ev+sa->lcm, ev+sb->lcm);
}

/* comparison for s-pairs while their lcms are in the local hash table */
static int spair_local_cmp(
    const void *a,
    const void *b
    )
{
  spair_t sa = *((spair_t *)a);
  spair_t sb = *((spair_t *)b);

  return (int)monomial_cmp(evl+sa.lcm, evl+sb.lcm);
}
