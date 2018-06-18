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
  const hl_t ca  = *((hl_t *)a);
  const hl_t cb  = *((hl_t *)b);

  return (int)(ca - cb);
}
#endif

static int matrix_row_initial_input_cmp_lex(
    const void *a,
    const void *b
    )
{
  len_t i;

  const hl_t va  = ((hl_t **)a)[0][2];
  const hl_t vb  = ((hl_t **)b)[0][2];

  const exp_t * const ea  = ev[va];
  const exp_t * const eb  = ev[vb];

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

  const hl_t va  = ((hl_t **)a)[0][2];
  const hl_t vb  = ((hl_t **)b)[0][2];

  const deg_t da = hd[va].deg;
  const deg_t db = hd[vb].deg;

  /* DRL */
  if (da < db) {
    return -1;
  } else {
    if (da != db) {
      return 1;
    }
  }

  const exp_t * const ea  = ev[va];
  const exp_t * const eb  = ev[vb];

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
  hl_t va, vb;
  /* compare pivot resp. column index */
  va  = ((hl_t **)a)[0][2];
  vb  = ((hl_t **)b)[0][2];
  if (va > vb) {
    return 1;
  }
  if (va < vb) {
    return -1;
  }
  /* same column index => compare density of row */
  va  = ((hl_t **)a)[0][0];
  vb  = ((hl_t **)b)[0][0];
  if (va > vb) {
    return 1;
  }
  if (va < vb) {
    return -1;
  }
  return 0;
}

static inline hl_t **sort_matrix_rows(
    hl_t **mat)
{
  qsort(mat, (unsigned long)nrows, sizeof(hl_t *), &matrix_row_cmp);
  return mat;
}

/* comparison for monomials (in local hash table) */
static int monomial_cmp_pivots_drl(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;

  const hd_t ha = hd[a];
  const hd_t hb = hd[b];
#if ORDER_COLUMNS
  /* first known pivots vs. tail terms */
  if (ha.idx != hb.idx) {
    if (ha.idx < hb.idx) {
      return 1;
    } else {
      return -1;
    }
  }
#endif

  /* then DRL */
  if (ha.deg > hb.deg) {
    return -1;
  } else {
    if (ha.deg != hb.deg) {
      return 1;
    }
  }

  const exp_t * const ea  = ev[a];
  const exp_t * const eb  = ev[b];

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

  return 0;
}

static int monomial_cmp_pivots_lex(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;

  const hd_t ha = hd[a];
  const hd_t hb = hd[b];
#if ORDER_COLUMNS
  /* first known pivots vs. tail terms */
  if (ha.idx != hb.idx) {
    if (ha.idx < hb.idx) {
      return 1;
    } else {
      return -1;
    }
  }
#endif

  const exp_t * const ea  = ev[a];
  const exp_t * const eb  = ev[b];

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

static inline int monomial_local_cmp_drl(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;

  const deg_t da = hdl[a].deg;
  const deg_t db = hdl[b].deg;

  /* DRL */
  if (da > db) {
    return 1;
  } else {
    if (da != db) {
      return -1;
    }
  }

  const exp_t * const ea  = evl[a];
  const exp_t * const eb  = evl[b];

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

static inline int monomial_cmp_drl(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;

  const deg_t da = hd[a].deg;
  const deg_t db = hd[b].deg;

  /* DRL */
  if (da > db) {
    return 1;
  } else {
    if (da != db) {
      return -1;
    }
  }

  const exp_t * const ea  = ev[a];
  const exp_t * const eb  = ev[b];

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

static inline int monomial_local_cmp_lex(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;

  const exp_t * const ea  = evl[a];
  const exp_t * const eb  = evl[b];

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

static inline int monomial_cmp_lex(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;

  const exp_t * const ea  = ev[a];
  const exp_t * const eb  = ev[b];

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
  const hl_t ma  = ((hl_t *)a)[0];
  const hl_t mb  = ((hl_t *)b)[0];

  return monomial_cmp_pivots_drl(ma, mb);
}

static int hcm_cmp_pivots_lex(
    const void *a,
    const void *b
    )
{
  const hl_t ma  = ((hl_t *)a)[0];
  const hl_t mb  = ((hl_t *)b)[0];

  return monomial_cmp_pivots_lex(ma, mb);
}

/* comparison for s-pairs once their lcms are in the global hash table */
static int spair_cmp_deglex(
    const void *a,
    const void *b
    )
{
  const hl_t la = ((spair_t *)a)->lcm;
  const hl_t lb = ((spair_t *)b)->lcm;

  if (hd[la].deg != hd[lb].deg) {
    return (hd[la].deg < hd[lb].deg) ? -1 : 1;
  } else {
    return (int)monomial_cmp(la, lb);
  }
}

static int spair_cmp_drl(
    const void *a,
    const void *b
    )
{
  const hl_t la = ((spair_t *)a)->lcm;
  const hl_t lb = ((spair_t *)b)->lcm;

  return (int)monomial_cmp(la, lb);
}

/* comparison for s-pairs while their lcms are in the local hash table */
static int spair_local_cmp(
    const void *a,
    const void *b
    )
{
  const hl_t la = ((spair_t *)a)->lcm;
  const hl_t lb = ((spair_t *)b)->lcm;

  return (int)monomial_local_cmp(la, lb);
}
