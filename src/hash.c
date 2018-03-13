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

static len_t htes   = 0;  /* hash table exponent at start */
static len_t nvars  = 0; /* number of variables */
static len_t bpv    = 0; /* bits per variable in divmask */
static len_t ndvars = 0; /* number of variables for divmask */

/* exponent block includes after the exponents also: */
#define HASH_DEG  (0) /* total degree */
#define HASH_SDM  (1) /* short divisor check mask */
#define HASH_VAL  (2) /* hash value */
#define HASH_DIV  (3) /* last divisor checked */
#define HASH_IND  (4) /* index into GB matrix */
#define HASH_LEN  (5) /* length of monomials */

/* hash tables consist of
 * ---------------------------------------------------
 * - exponent vectors ev
 * - meta data like degrees, hash values, etc. md
 * - the mapping from the hashed position value to the
 *   position in the ev/md array
 * --------------------------------------------------- */

/* global hash table for storing elements in basis */
static exp_t *ev    = NULL;
static len_t *md    = NULL;
static len_t esize  = 0;
static len_t eload  = 0;

/* map with index from monomials to exponents */
static len_t *map   = NULL;
static len_t msize  = 0;

/* local hash table during symbolic preprocessing */
static exp_t *evl     = NULL;
static len_t *mdl     = NULL;
static len_t elsize   = 0;
static len_t elload   = 0;

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
    )
{
  int32_t i, len;

  /* generate map */
  bpv   = (len_t)((CHAR_BIT * sizeof(sdm_t)) / (unsigned long)nvars);
  if (bpv == 0)
    bpv++;
  ndvars  = (unsigned long)nvars < (CHAR_BIT * sizeof(sdm_t)) ?
    nvars : (len_t)((CHAR_BIT * sizeof(sdm_t)));
  msize = (len_t)pow(2, htes);
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
  esize = msize/2;
  len   = HASH_LEN * esize;
  md    = calloc((unsigned long)len, sizeof(len_t));
  /* keep first entry empty for faster divisibility checks */
  len   = nvars * esize;
  ev    = calloc((unsigned long)len, sizeof(exp_t));
  eload = 1;
}

static void initialize_local_hash_table(
    )
{
  int32_t len = 0;

  /* generate map */
  mlsize  = (len_t)pow(2, htes-5);
  mapl    = calloc((unsigned long)mlsize, sizeof(len_t));

  elsize  = mlsize/2;
  len     = HASH_LEN * elsize;
  /* keep first entry empty for faster divisibility checks */
  mdl     = calloc((unsigned long)len, sizeof(len_t));
  len     = nvars * elsize;
  evl     = calloc((unsigned long)len, sizeof(exp_t));
  elload  = 1;
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
  if (md) {
    free(md);
    md  = NULL;
  }
  fc      = 0;
  nvars   = 0;
  esize   = 0;
  eload   = 0;
  msize   = 0;
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
  if (mdl) {
    free(mdl);
    mdl = NULL;
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
  printf("enlarge %d => %d\n", msize, 2*msize);
  exp_t *e;

  int32_t len = 2 * HASH_LEN * eload;
  md  = realloc(md, (unsigned long)len * sizeof(len_t));
  len = len/2;
  memset(md+eload*HASH_LEN, 0, (unsigned long)len * sizeof(len_t));

  len = 2 * nvars * eload;
  ev  = realloc(ev, (unsigned long)len * sizeof(exp_t));
  len = len/2;
  memset(ev+eload*nvars, 0, (unsigned long)len * sizeof(exp_t));

  msize = 2 * msize;
  map   = realloc(map, (unsigned long)msize * sizeof(len_t));
  memset(map, 0, (unsigned long)msize * sizeof(len_t));

  /* reinsert known elements */
  for (i = 1; i < eload; ++i) {
    e = ev + i*nvars;
    h = (md+(i*HASH_LEN))[HASH_VAL];

    /* probing */
    k = h;
    for (j = 0; j < msize; ++j) {
      k = (k+j) & (msize-1);
      if (map[k]) {
        continue;
      }
      map[k]  = i;
      break;
    }
  }
  esize = 2*esize;
}

static inline sdm_t generate_short_divmask(
    const exp_t *a
    )
{
  int32_t i, j;
  sdm_t res = 0;
  sdm_t ctr = 0;

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
  exp_t *max_exp  = (exp_t *)malloc((unsigned long)ndvars * sizeof(exp_t));
  exp_t *min_exp  = (exp_t *)malloc((unsigned long)ndvars * sizeof(exp_t));

  exp_t *e  = ev + 1*nvars;

  len_t *m;

  /* get initial values from first hash table entry */
  for (i = 0; i < ndvars; ++i) {
    max_exp[i]  = min_exp[i]  = e[i];
  }

  /* get maximal and minimal exponent element entries in hash table */
  for (i = 2; i < eload; ++i) {
    e = ev + i*nvars;
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
  for (i = 1; i < eload; ++i) {
    e = ev + i*nvars;
    m = md + i*HASH_LEN;
    m[HASH_SDM] = generate_short_divmask(e);
  }

  free(max_exp);
  free(min_exp);
}

/* returns zero if a is not divisible by b, else 1 is returned */
static inline len_t check_monomial_division_local(
    const len_t a,
    const len_t b
    )
{
  int32_t i;

  const len_t ma = (mdl + a*HASH_LEN)[HASH_SDM];
  const len_t mb = (mdl + b*HASH_LEN)[HASH_SDM];
  /* short divisor mask check */
  if (mb & ~ma) {
    return 0;
  }

  /* degree check */
  /* if (ea[HASH_DEG] < eb[HASH_DEG]) {
   *   return 0;
   * } */

  const exp_t * const ea  = evl + a*nvars;
  const exp_t * const eb  = evl + b*nvars;

  /* exponent check */
  if (ea[0] < eb[0] || ea[nvars-1] < eb[nvars-1]) {
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

static inline len_t check_monomial_division(
    const len_t a,
    const len_t b
    )
{
  int32_t i;

  const len_t ma = (md + a*HASH_LEN)[HASH_SDM];
  const len_t mb = (md + b*HASH_LEN)[HASH_SDM];
  /* short divisor mask check */
  if (mb & ~ma) {
    return 0;
  }

  /* degree check */
  /* if (ea[HASH_DEG] < eb[HASH_DEG]) {
   *   return 0;
   * } */

  const exp_t * const ea  = ev + a*nvars;
  const exp_t * const eb  = ev + b*nvars;

  /* exponent check */
  if (ea[0] < eb[0] || ea[nvars-1] < eb[nvars-1]) {
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

static inline len_t insert_in_global_hash_table(
    const exp_t *a
    )
{
  int32_t i, k, pos, deg;
  exp_t *e;
  len_t *m;
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
    e = ev + map[k]*nvars;
    m = md + map[k]*HASH_LEN;
    if (m[HASH_VAL] != h) {
      continue;
    }
    if (memcmp(e, a, (unsigned long)nvars * sizeof(exp_t)) == 0) {
      return map[k];
    }
  }

  /* add element to hash table */
  map[k]  = pos = eload;
  e   = ev + pos*nvars;
  m   = md + pos*HASH_LEN; 
  deg = 0;
  for (i = 0; i < nvars; ++i) {
    e[i]  =   a[i];
    deg   +=  a[i];
  }
  m[HASH_DEG] = deg;
  m[HASH_SDM] = generate_short_divmask(e);
  m[HASH_VAL] = h;

  eload++;
  if (eload >= esize) {
    enlarge_global_hash_table();
  }

  return pos;
}


static inline len_t insert_in_local_hash_table(
    const exp_t *a
    )
{
  int32_t i, k, pos, deg;
  exp_t *e;
  len_t *m;
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
    e = evl + mapl[k]*nvars;
    m = mdl + mapl[k]*HASH_LEN;
    if (m[HASH_VAL] != h) {
      continue;
    }
    if (memcmp(e, a, (unsigned long)nvars * sizeof(exp_t)) == 0) {
      return mapl[k];
    }
  }

  /* add element to hash table */
  mapl[k] = pos = elload;
  e   = evl + pos*nvars;
  m   = mdl + pos*HASH_LEN;
  deg = 0;
  for (i = 0; i < nvars; ++i) {
    e[i]  =   a[i];
    deg   +=  a[i];
  }
  m[HASH_DEG] = deg;
  m[HASH_SDM] = generate_short_divmask(e);
  m[HASH_VAL] = h;

  elload++;

  return pos;
}

static inline void reset_local_hash_table(
    const len_t size
    )
{
  /* is there still enough space in the local table? */
  if (size >= (elsize-elload)) {
    if (2*size >= mlsize) {
      while (2*size >= mlsize) {
        elsize  = 2 * elsize;
        mlsize  = 2 * mlsize;
      }
      evl   = realloc(evl, (unsigned long)(elsize * nvars) * sizeof(exp_t));
      mdl   = realloc(mdl, (unsigned long)(elsize * HASH_LEN) * sizeof(len_t));
      mapl  = realloc(mapl, (unsigned long)mlsize * sizeof(len_t));
    }
    memset(evl, 0, (unsigned long)(elsize * nvars) * sizeof(exp_t));
    memset(mdl, 0, (unsigned long)(elsize * HASH_LEN) * sizeof(len_t));
    memset(mapl, 0, (unsigned long)mlsize * sizeof(len_t));

    elload  = 1;
  }
}

/* note that the product insertion, i.e. monomial x polynomial
 * is only needed for the local hash table. in the global one we
 * only add the basis elements, i.e. no multiplication is applied. */
static inline len_t insert_in_global_hash_table_product(
    const len_t p1,
    const len_t p2
    )
{
  int32_t i, k, pos;
  exp_t *e, *n;
  len_t *m;

  const len_t * const m1  = md + (p1*HASH_LEN);
  const len_t * const m2  = md + (p2*HASH_LEN);
  const int32_t h = m1[HASH_VAL] + m2[HASH_VAL];

  n = ev + eload*nvars;
  exp_t *a1 = ev + p1*nvars;
  exp_t *a2 = ev + p2*nvars;
  for (i = 0; i < nvars; ++i) {
    n[i]  = a1[i] + a2[i];
  }
  /* probing */
  k = h;
  for (i = 0; i < msize; ++i) {
    k = (k+i) & (msize-1);
    if (!map[k]) {
      break;
    }
    e = ev + map[k]*nvars;
    m = md + map[k]*HASH_LEN;
    if (m[HASH_VAL] != h) {
      continue;
    }
    if (memcmp(n, e, (unsigned long)nvars * sizeof(exp_t)) == 0) {
      return map[k];
    }
  }

  /* add element to hash table */
  map[k]  = pos = eload;
  m = md+eload*HASH_LEN;
  m[HASH_DEG] = m1[HASH_DEG] + m2[HASH_DEG];
  m[HASH_SDM] = generate_short_divmask(n);
  m[HASH_VAL] = h;

  eload++;
  if (eload >= esize) {
    enlarge_global_hash_table();
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
  const len_t * const ma  = md + a*HASH_LEN;
  const len_t * const mb  = md + b*HASH_LEN;
  const len_t * const ml  = mdl + lcm*HASH_LEN;

  if (ml[HASH_DEG] != (ma[HASH_DEG] + mb[HASH_DEG])) {
    return 0;
  }
  if (ml[HASH_VAL] != (ma[HASH_VAL] + mb[HASH_VAL])) {
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
  const exp_t * const ea = ev + a*nvars;
  const exp_t * const eb = ev + b*nvars;

  exp_t *e  = alloca((unsigned long)nvars * sizeof(exp_t));
  for (i = 0; i < nvars; ++i) {
    e[i]  = ea[i] < eb[i] ? eb[i] : ea[i];
  }
  /* goes into local hash table for spairs */
  return insert_in_local_hash_table(e);
}

static inline len_t monomial_multiplication(
    const len_t a,
    const len_t b
    )
{
  return insert_in_global_hash_table_product(a, b);
}

/* we try monomial division including check if divisibility is
 * fulfilled. */
static inline len_t monomial_division_with_check(
    const len_t a,
    const len_t b
    )
{
  int32_t i;

  const len_t ma  = md[a*HASH_LEN + HASH_SDM];
  const len_t mb  = md[b*HASH_LEN + HASH_SDM];
  /* short divisor mask check */
  if (mb & ~ma) {
    return 0;
  }

  const exp_t * const ea  = ev + a*nvars;
  const exp_t * const eb  = ev + b*nvars;
  /* exponent check */
  if (ea[0] < eb[0] || ea[nvars-1] < eb[nvars-1]) {
    return 0;
  }

  exp_t *e = (exp_t *)alloca((unsigned long)nvars * sizeof(exp_t));
  e[0]  = ea[0] - eb[0];

  i = nvars & 1 ? 1 : 0;
  for (; i < nvars; i += 2) {
    if (ea[i] < eb[i] || ea[i+1] < eb[i+1]) {
      return 0;
    } else {
      e[i]    = ea[i] - eb[i];
      e[i+1]  = ea[i+1] - eb[i+1];
    }
  }
  return insert_in_global_hash_table(e);
}

/* it is assumed that b divides a, thus no tests for
 * divisibility at all */
static inline len_t monomial_division_no_check(
    const len_t a,
    const len_t b
    )
{
  int32_t i;

  const exp_t * const ea  = ev + a*nvars;
  const exp_t * const eb  = ev + b*nvars;

  exp_t *e = (exp_t *)alloca((unsigned long)nvars * sizeof(exp_t));

  i = nvars & 1 ? 1 : 0;
  for (; i < nvars; i += 2) {
    e[i]    = ea[i]   - eb[i];
    e[i+1]  = ea[i+1] - eb[i+1];
  }
  e[0]  = ea[0] - eb[0];
  return insert_in_global_hash_table(e);
}

static inline val_t *multiplied_polynomial_to_matrix_row(
    const len_t mult,
    const val_t *poly
    )
{
  int32_t i;

  /* printf("mulitplied row: "); */
  val_t *row  = (val_t *)malloc((unsigned long)poly[0] * sizeof(val_t));
  row[0]  = poly[0]; /* length */
  row[1]  = poly[1]; /* loop unroll offset */
  for (i = 2; i < poly[1]; i += 2) {
    row[i]    = monomial_multiplication(mult, poly[i]);
    row[i+1]  = poly[i+1];
  }
  for (;i < poly[0]; i += 8) {
    row[i]    = monomial_multiplication(mult, poly[i]);
    row[i+1]  = poly[i+1];
    row[i+2]  = monomial_multiplication(mult, poly[i+2]);
    row[i+3]  = poly[i+3];
    row[i+4]  = monomial_multiplication(mult, poly[i+4]);
    row[i+5]  = poly[i+5];
    row[i+6]  = monomial_multiplication(mult, poly[i+6]);
    row[i+7]  = poly[i+7];
  }
  /* printf("multiplied polys added\n");
   * for (int32_t p = 2; p < poly[0]; p += 2) {
   *   for (int32_t o = 0; o < nvars; ++o) {
   *     printf("%d ", (ev+poly[p]*nvars)[o]);
   *   }
   *   printf(" || ");
   * }
   * printf("\n multiplied with  ");
   * for (int32_t o = 0; o < nvars; ++o) {
   *   printf("%d ", (ev+mult*nvars)[o]);
   * }
   * printf(" |&| ");
   * printf("\nmult to\n");
   * for (int32_t p = 2; p < row[0]; p += 2) {
   *   for (int32_t o = 0; o < nvars; ++o) {
   *     printf("%d ", (ev+row[p]*nvars)[o]);
   *   }
   *   printf(" |&| ");
   * }
   * printf("\nttt\n");
   * printf("\n"); */

  return row;
}
