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
 * \file hash.c
 * \brief Implementation of the hash tables used in gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

/* One general note:
 * Hashing is not that easy when it comes to "real", i.e. bigger
 * examples. One global hash table is not enough since we would
 * store too many information in the hash table due to intermediate
 * symbolic preprocessing computation, i.e. multiplying polynomials
 * with some other monomials. Freeing the complete hash table after
 * we finished a matrix and reconstruct the hash table with entries
 * only for the elements in the basis is too time consuming. Thus we
 * decided to use two hash tables: One global hash table keeping
 * hashes of exponents of polynomials in the basis. For each new
 * symbolic preprocessing step we construct a new local hash table
 * which stores the intermediate multiplied polynomial exponents.
 * This table is used for the linear algebra. Afterwards only the
 * exponent hashes of the new basis polynomials are added to the
 * global hash table, the local one is freed. Clearly, we might have
 * the same exponents living in both, the global and the local hash
 * table at the same time during some parts of the computation, still
 * it is an OK-ish trade-off between memory and speed. */

/* The idea of the structure of the hash table is taken from an
 * implementation by Roman Pearce and Michael Monagan in Maple. */

static val_t fc     = 0; /* field characteristic */
static len_t nvars  = 0; /* number of variables */
static len_t bpv    = 0; /* bits per variable in divmask */
static len_t ndvars = 0; /* number of variables for divmask */

/* exponent block includes after the exponents also: */
#define HASH_DEG  (nvars+0) /* total degree */
#define HASH_SDM  (nvars+1) /* short divisor check mask */
#define HASH_VAL  (nvars+2) /* hash value */
#define HASH_DIV  (nvars+3) /* last divisor checked */
#define HASH_IND  (nvars+4) /* index into GB matrix */
#define HASH_LEN  (nvars+5) /* length of monomials */

/* global hash table for storing elements in basis */
static exp_t *ev    = NULL;
static len_t esize  = 0;
static len_t eload  = 1;

/* map with index from monomials to exponents */
static len_t *map   = NULL;
static len_t msize  = 0;

/* local hash table during symbolic preprocessing */
static exp_t *evl     = NULL;
static len_t elsize   = 0;
static len_t elload   = 1;

static len_t *mapl    = NULL;
static len_t mlsize   = 0;

/* random values for generating hash values */
static val_t *rv  = NULL;

/* divisor map for short divisibility tests */
static sdm_t *dm  = NULL;

/* pseudo random number generator for hash value
 * generation */
uint32_t rseed  = 2463534242;
static val_t pseudo_random_number_generator(
    void
    )
{
	rseed ^= (rseed << 13);
	rseed ^= (rseed >> 17);
	rseed ^= (rseed << 5);
	return (val_t)rseed;
}

static void initialize_global_hash_table(
    len_t nr_vars,
    len_t ht_size,
    val_t field_char
    )
{
  int32_t i;

  /* generate map */
  nvars = nr_vars;
  fc    = field_char;
  bpv   = (CHAR_BIT * sizeof(sdm_t)) / (unsigned long)nvars;
  if (bpv == 0)
    bpv++;
  ndvars  = (unsigned long)nvars < (CHAR_BIT * sizeof(sdm_t)) ?
    nvars : (CHAR_BIT * sizeof(sdm_t));
  msize = (len_t)pow(2, ht_size);
  map   = calloc((unsigned long)msize, sizeof(len_t));

  /* generate divmask map */
  dm  = calloc((unsigned long)(ndvars * bpv), sizeof(sdm_t));

  /* generate random values */
  rv  = calloc((unsigned long)nvars, sizeof(val_t));
  for (i = nvars-1; i >= 0; --i) {
    /* random values should not be zero */
    rv[i] = pseudo_random_number_generator() | 1;
  }
  /* generate exponent vector */
  esize = HASH_LEN + (msize/2);
  /* keep first entry empty for faster divisibility checks */
  eload = HASH_LEN;
  ev    = calloc((unsigned long)esize, sizeof(exp_t));
}

static void initialize_local_hash_table(
    len_t ht_size
    )
{
  /* generate map */
  mlsize  = (len_t)pow(2, ht_size);
  mapl    = calloc((unsigned long)mlsize, sizeof(len_t));

  /* generate exponent vector */
  elsize  = HASH_LEN + (msize/2);
  /* keep first entry empty for faster divisibility checks */
  elload  = HASH_LEN;
  evl     = calloc((unsigned long)elsize, sizeof(exp_t));
}

static void free_global_hash_table(
    void
    )
{
  if (map) {
    free(map);
    map = NULL;
  }
  if (dm) {
    free(dm);
    dm  = NULL;
  }
  if (rv) {
    free(rv);
    rv  = NULL;
  }
  if (ev) {
    free(ev);
    ev  = NULL;
  }
  fc    = 0;
  nvars = 0;
  esize = 0;
  eload = 0;
  msize = 0;
}

static void free_local_hash_table(
    void
    )
{
  if (mapl) {
    free(mapl);
    mapl  = NULL;
  }
  if (evl) {
    free(evl);
    evl = NULL;
  }
  elsize  = 0;
  elload  = 0;
  mlsize  = 0;
}

/* we just double the hash table size */
static void enlarge_global_hash_table(
    void
    )
{
  int32_t h, i, j, k;
  int32_t *e;

  esize = 2 * esize;
  ev    = realloc(ev, (unsigned long)esize * sizeof(exp_t));

  msize = 2 * msize;
  map   = realloc(map, (unsigned long)msize * sizeof(len_t));

  memset(map, 0, (unsigned long)msize * sizeof(len_t));

  /* reinsert known elements */
  for (i = HASH_LEN; i < eload; i += HASH_LEN) {
    e = ev + i;
    h = e[HASH_VAL];

    /* probing */
    k = h;
    for (j = 0; j < msize; ++j) {
      k = (k+h) & (msize-1);
      if (map[k]) {
        continue;
      }
      map[k]  = i;
      break;
    }
  }
}

static void enlarge_local_hash_table(
    void
    )
{
  int32_t h, i, j, k;
  int32_t *e;

  elsize  = 2 * elsize;
  evl     = realloc(evl, (unsigned long)elsize * sizeof(exp_t));

  mlsize  = 2 * mlsize;
  mapl    = realloc(mapl, (unsigned long)mlsize * sizeof(len_t));

  memset(mapl, 0, (unsigned long)mlsize * sizeof(len_t));

  /* reinsert known elements */
  for (i = HASH_LEN; i < elload; i += HASH_LEN) {
    e = evl + i;
    h = e[HASH_VAL];

    /* probing */
    k = h;
    for (j = 0; j < mlsize; ++j) {
      k = (k+h) & (mlsize-1);
      if (mapl[k]) {
        continue;
      }
      mapl[k]  = i;
      break;
    }
  }
}

static inline sdm_t generate_short_divmask(
    const exp_t *a
    )
{
  int32_t i, j;
  int32_t res = 0;
  int32_t ctr = 0;

  for (i = 0; i < ndvars; ++i) {
    for (j = 0; j < bpv; ++j) {
      if (a[i] >= dm[ctr]) {
        res |= 1 << ctr;
      }
      ctr++;
    }
  }
 
  return res;
}

/* note: we calculate the divmask after reading in the input generators. those
 * are first stored in the local hash table. thus we use the local exponents to
 * generate the divmask */
static inline void calculate_divmask(
    void
    )
{
  int32_t i, j, steps;
  int32_t ctr = 0;
  deg_t *max_exp  = (deg_t *)malloc((unsigned long)ndvars * sizeof(deg_t));
  deg_t *min_exp  = (deg_t *)malloc((unsigned long)ndvars * sizeof(deg_t));

  exp_t *e  = evl + HASH_LEN;

  /* get initial values from first hash table entry */
  for (i = 0; i < ndvars; ++i) {
    max_exp[i]  = min_exp[i]  = e[i];
  }

  /* get maximal and minimal exponent element entries in hash table */
  for (i = 2*HASH_LEN; i < elload; i = i + HASH_LEN) {
    e = evl + i;
    for (j = 0; j < ndvars; ++j) {
      if (e[j] > max_exp[j]) {
        max_exp[j]  = e[j];
        continue;
      }
      if (e[j] < min_exp[j]) {
        min_exp[j]  = e[j];
      }
    }
  }

  /* calculate average values for generating divmasks */
  for (i = 0; i < ndvars; ++i) {
    steps = (max_exp[i] - min_exp[i]) / bpv;
    if (steps == 0)
      steps++;
    for (j = 0; j < bpv; ++j) {
      dm[ctr++] = min_exp[i]++;
    }
  }

  /* initialize divmasks for elements already added to hash table */
  for (i = HASH_LEN; i < elload; i = i + HASH_LEN) {
    e = evl + i;
    e[HASH_SDM] = generate_short_divmask(e);
  }

  free(max_exp);
  free(min_exp);
}

static inline len_t insert_in_global_hash_table(
    const exp_t *a
    )
{
  int32_t i, j, k, pos, deg;
  exp_t *e;
  int32_t h = 0;

  /* generate hash value */
  for (i = 0; i < nvars; ++i) {
    h +=  rv[i] * a[i];
  }

  /* probing */
  k = h;
  for (i = 0; i < msize; ++i) {
    k = (k+i) & (msize-1);
    if (!map[k]) {
      break;
    }
    e = ev + map[k];
    if (e[HASH_VAL] != h) {
      continue;
    }
    for (j = 0; j < nvars; ++j) {
      if (e[i] != a[i]) {
        break;
      }
    }
    if (i == nvars) {
      return map[k];
    }
  }

  /* add element to hash table */
  pos = eload;
  e   = ev + pos;
  deg = 0;
  for (i = 0; i < nvars; ++i) {
    e[i]  =   a[i];
    deg   +=  a[i];
  }
  e[HASH_DEG] = deg;
  e[HASH_SDM] = generate_short_divmask(e);
  e[HASH_VAL] = h;
  e[HASH_DIV] = 0;
  e[HASH_IND] = 0;
  map[k]      = pos;

  eload +=  HASH_LEN;
  if (eload >= esize) {
    enlarge_global_hash_table();
  }

  return pos;
}

static inline len_t insert_in_local_hash_table(
    const exp_t *a
    )
{
  int32_t i, j, k, pos, deg;
  exp_t *e;
  int32_t h = 0;

  /* generate hash value */
  for (i = 0; i < nvars; ++i) {
    h +=  rv[i] * a[i];
  }

  /* probing */
  k = h;
  for (i = 0; i < mlsize; ++i) {
    k = (k+i) & (mlsize-1);
    if (!mapl[k]) {
      break;
    }
    e = evl + mapl[k];
    if (e[HASH_VAL] != h) {
      continue;
    }
    for (j = 0; j < nvars; ++j) {
      if (e[i] != a[i]) {
        break;
      }
    }
    if (i == nvars) {
      return mapl[k];
    }
  }

  /* add element to hash table */
  pos = elload;
  e   = evl + pos;
  deg = 0;
  for (i = 0; i < nvars; ++i) {
    e[i]  =   a[i];
    deg   +=  a[i];
  }
  e[HASH_DEG] = deg;
  e[HASH_VAL] = h;
  e[HASH_DIV] = 0;
  e[HASH_IND] = 0;
  mapl[k]     = pos;

  elload  +=  HASH_LEN;
  if (elload >= elsize) {
    enlarge_local_hash_table();
  }

  return pos;
}

static inline void clear_local_hash_table(
    void
    )
{
  memset(evl, 0, (unsigned long)(elload * HASH_LEN) * sizeof(exp_t));
  memset(mapl, 0, (unsigned long)(mlsize) * sizeof(len_t));
  elload  = 0;
}

/* note that the product insertion, i.e. monomial x polynomial
 * is only needed for the local hash table. in the global one we
 * only add the basis elements, i.e. no multiplication is applied. */
static inline len_t insert_in_local_hash_table_product(
    const exp_t *a1,
    const exp_t *a2
    )
{
  int32_t i, j, k, pos;
  exp_t *e;

  int32_t h = a1[HASH_VAL] + a2[HASH_VAL];

  /* probing */
  k = h;
  for (i = 0; i < mlsize; ++i) {
    k = (k+i) & (mlsize-1);
    if (!mapl[k]) {
      break;
    }
    e = evl + mapl[k];
    if (e[HASH_VAL] != h) {
      continue;
    }
    for (j = 0; j < nvars; ++j) {
      if (e[i] != a1[i] + a2[i]) {
        break;
      }
    }
    if (i == nvars) {
      return mapl[k];
    }
  }

  /* add element to hash table */
  pos = elload;
  e   = evl + pos;
  for (i = 0; i < nvars; ++i) {
    e[i]  =   a1[i] + a2[i];
  }
  e[HASH_DEG] = a1[HASH_DEG] + a2[HASH_DEG];
  e[HASH_VAL] = h;
  e[HASH_DIV] = 0;
  e[HASH_IND] = 0;
  mapl[k]     = pos;

  elload  +=  HASH_LEN;
  if (elload >= elsize) {
    enlarge_local_hash_table();
  }

  return pos;
}

/* we can check equality of lcm and multiplication of two monomials
 * by their hash values. If both hash values are NOT the same, then
 * the corresponding exponent vectors CANNOT be the same. */
static inline int lcm_equals_multiplication(
    const len_t a,
    const len_t b,
    const len_t lcm
    )
{
  const exp_t * const ea  = ev + a;
  const exp_t * const eb  = ev + b;
  const exp_t * const el  = ev + lcm;

  if (el[HASH_DEG] != (ea[HASH_DEG] + eb[HASH_DEG])) {
    return 0;
  }
  if (el[HASH_VAL] != (ea[HASH_VAL] + eb[HASH_VAL])) {
    return 0;
  } else {
    /* both have the same degree and the same hash value, either they
     * are the same or our hashing is broken resp. really bad */
    return 1;
  }
}

static inline len_t get_lcm(
    const len_t a,
    const len_t b
    )
{
  int32_t i;

  /* exponents of basis elements, thus from global hash table */
  const exp_t * const ea = ev + a;
  const exp_t * const eb = ev + b;

  exp_t *e  = alloca((unsigned long)nvars * sizeof(exp_t));
  for (i = 0; i < nvars; ++i) {
    e[i]  = ea[i] < eb[i] ? eb[i] : ea[i];
  }
  /* goes into global hash table for spairs */
  return insert_in_global_hash_table(e);
}

static inline len_t monomial_multiplication(
    const len_t a,
    const len_t b
    )
{
  /* a is the multiplier monomial living the local table,
   * b is a monomial from a polynomial in the basis thus
   * living in the global table */
  const exp_t * const ea = evl  + a;
  const exp_t * const eb = ev   + b;

  return insert_in_local_hash_table_product(ea, eb);
}

/* returns zero if a is not divisible by b, else 1 is returned */
static inline len_t check_monomial_division(
    const len_t a,
    const len_t b
    )
{
  int32_t i;

  const exp_t * const ea  = ev + a;
  const exp_t * const eb  = ev + b;
  /* short divisor mask check */
  if (eb[HASH_SDM] & ~ea[HASH_SDM]) {
    return 0;
  }

  /* degree check */
  if (ea[HASH_DEG] < eb[HASH_DEG]) {
    return 0;
  }

  /* exponent check */
  if (ea[0] < eb[0]) {
    return 0;
  }
  i = nvars & 1 ? 1 : 0;
  for (; i < nvars; i += 2) {
    if (ea[i] < eb[i] || ea[i+1] < eb[i+1]) {
      return 0;
    }
  }
  return 1;
}

/* it is assumed that b divides a, thus no tests for
 * divisibility at all */
static inline len_t monomial_division(
    const len_t a,
    const len_t b
    )
{
  int32_t i;

  const exp_t * const ea  = ev + a;
  const exp_t * const eb  = ev + b;

  exp_t *e = (exp_t *)alloca((unsigned long)nvars * sizeof(exp_t));

  i = nvars & 1 ? 1 : 0;
  for (; i < nvars; i += 2) {
    e[i]    = ea[i]   - eb[i];
    e[i+1]  = ea[i+1] - eb[i+1];
  }
  e[0]  = ea[0] - eb[0];
  return insert_in_local_hash_table(e);
}

static inline val_t *multiplied_polynomial_to_matrix_row(
    const len_t mult,
    const val_t *poly
    )
{
  int32_t i;

  val_t *row  = (val_t *)malloc((unsigned long)poly[0] * sizeof(val_t));
  row[0]  = poly[0];
  row[1]  = 0; /* no new piv */
  for (i = 2; i < poly[0]; i += 2) {
    row[i]    = monomial_multiplication(mult, poly[i]);
    row[i+1]  = poly[i+1];
  }

  return row;
}
