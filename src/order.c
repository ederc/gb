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
  len_t i;
  int64_t va, vb;
  const int64_t hl  = HASH_LEN;

  va  = ((val_t **)a)[0][2] * hl;
  vb  = ((val_t **)b)[0][2] * hl;

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
  len_t i;
  int64_t va, vb;
  const int64_t hl  = HASH_LEN;

  va  = ((val_t **)a)[0][2] * hl;
  vb  = ((val_t **)b)[0][2] * hl;

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
  for (i = nvars; i > 0; --i) {
    if (ea[i-1] < eb[i-1]) {
      return -1;
    } else {
      if (ea[i-1] != eb[i-1]) {
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
  val_t va, vb;
  /* compare pivot resp. column index */
  va  = ((val_t **)a)[0][2];
  vb  = ((val_t **)b)[0][2];
  if (va > vb) {
    return 1;
  }
  if (va < vb) {
    return -1;
  }
  /* same column index => compare density of row */
  va  = ((val_t **)a)[0][0];
  vb  = ((val_t **)b)[0][0];
  if (va > vb) {
    return 1;
  }
  if (va < vb) {
    return -1;
  }
  return 0;
}

static inline val_t **sort_matrix_rows(
    val_t **mat)
{
  qsort(mat, (unsigned long)nrows, sizeof(val_t *), &matrix_row_cmp);
  return mat;
}

/* comparison for monomials (in local hash table) */
static int monomial_cmp_pivots_drl(
    const len_t a,
    const len_t b
    )
{
  len_t i;
  const int64_t hl  = HASH_LEN;

  const exp_t * const ea  = ev + a*hl;
  const exp_t * const eb  = ev + b*hl;

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
    for (i = nvars; i > 0; --i) {
      if (ea[i-1] > eb[i-1]) {
        return 1;
      } else {
        if (ea[i-1] != eb[i-1]) {
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
  len_t i;
  const int64_t hl  = HASH_LEN;

  const exp_t * const ea  = ev + a*hl;
  const exp_t * const eb  = ev + b*hl;

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
  len_t i;

  if (ea[HASH_DEG] > eb[HASH_DEG]) {
    return 1;
  } else {
    if (ea[HASH_DEG] != eb[HASH_DEG]) {
      return -1;
    }
  }

  for (i = nvars; i > 0; --i) {
    if (ea[i-1] < eb[i-1]) {
      return 1;
    } else {
      if (ea[i-1] != eb[i-1]) {
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
  len_t i;

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

  const int64_t hl  = HASH_LEN;
  if ((ev+sa->lcm*hl)[HASH_DEG] != (ev+sb->lcm*hl)[HASH_DEG]) {
    return ((ev+sa->lcm*hl)[HASH_DEG] < (ev+sb->lcm*hl)[HASH_DEG]) ? -1 : 1;
  } else {
    return (int)monomial_cmp(ev+sa->lcm*hl, ev+sb->lcm*hl);
  }
}

static int spair_cmp_drl(
    const void *a,
    const void *b
    )
{
  spair_t *sa = (spair_t *)a;
  spair_t *sb = (spair_t *)b;

  const int64_t hl  = HASH_LEN;

  return (int)monomial_cmp(ev+sa->lcm*hl, ev+sb->lcm*hl);
}

/* comparison for s-pairs while their lcms are in the local hash table */
static int spair_local_cmp(
    const void *a,
    const void *b
    )
{
  spair_t sa = *((spair_t *)a);
  spair_t sb = *((spair_t *)b);

  const int64_t hl  = HASH_LEN;

  return (int)monomial_cmp(evl+sa.lcm*hl, evl+sb.lcm*hl);
}
