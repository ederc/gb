/* gb: Gröbner Basis
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
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

static len_t nvars  = 0;

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
    len_t ht_size
    )
{
  int32_t i;

  /* generate map */
  nvars = nr_vars;
  msize = (len_t)pow(2, ht_size);
  map   = calloc((unsigned long)msize, sizeof(len_t));

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
  mapl    = calloc((unsigned long)msize, sizeof(len_t));

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
  if (rv) {
    free(rv);
    rv  = NULL;
  }
  if (ev) {
    free(ev);
    ev  = NULL;
  }
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

static inline len_t insert_in_global_hash_table(
    const exp_t *a
    )
{
  int32_t i, j, k, pos, deg;
  int32_t *e;
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
  int32_t *e;
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

/* note that the product insertion, i.e. monomial x polynomial
 * is only needed for the local hash table. in the global one we
 * only add the basis elements, i.e. no multiplication is applied. */
static inline len_t insert_in_local_hash_table_product(
    const exp_t *a1,
    const exp_t *a2
    )
{
  int32_t i, j, k, pos;
  int32_t *e;

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
