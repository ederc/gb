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
    void
    )
{
  len_t i;
  hl_t j;

  /* generate map */
  bpv   = (len_t)((CHAR_BIT * sizeof(sdm_t)) / (unsigned long)nvars);
  if (bpv == 0) {
    bpv++;
  }
  ndvars  = (unsigned long)nvars < (CHAR_BIT * sizeof(sdm_t)) ?
    nvars : (len_t)((CHAR_BIT * sizeof(sdm_t)));
  hsz   = (hl_t)pow(2, htes);
  hmap  = calloc((unsigned long)hsz, sizeof(hl_t));

  /* generate divmask map */
  dm  = calloc((unsigned long)(ndvars * bpv), sizeof(sdm_t));

  /* generate random values */
  rv  = calloc((unsigned long)nvars, sizeof(val_t));
  for (i = nvars; i > 0; --i) {
    /* random values should not be zero */
    rv[i-1] = pseudo_random_number_generator() | 1;
  }
  /* generate exponent vector */
  esz = hsz/2;
  /* keep first entry empty for faster divisibility checks */
  eld = 1;
  hd  = (hd_t *)calloc((unsigned long)esz, sizeof(hd_t));
  ev  = (exp_t **)malloc((unsigned long)esz * sizeof(exp_t *));
  if (ev == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  exp_t *tmp  = (exp_t *)malloc(
      (unsigned long)nvars * (unsigned long)esz * sizeof(exp_t));
  if (tmp == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  for (j = 0; j < esz; ++j) {
    ev[j]  = tmp + (unsigned long)(j*nvars);
  }

  etmp  = (exp_t *)malloc((unsigned long)esz * sizeof(exp_t));
}

static void initialize_local_hash_table(
    void
    )
{
  hl_t j;

  /* generate map */
  hlsz  = (hl_t)pow(2, htes-5);
  hmapl = calloc((unsigned long)hlsz, sizeof(hl_t));

  /* generate exponent vector */
  elsz  = hlsz/2;
  /* keep first entry empty for faster divisibility checks */
  elld  = 1;
  hdl   = (hd_t *)calloc((unsigned long)elsz, sizeof(hd_t));
  evl   = (exp_t **)malloc((unsigned long)elsz * sizeof(exp_t *));
  if (evl == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  exp_t *tmp  = (exp_t *)malloc(
      (unsigned long)nvars * (unsigned long)elsz * sizeof(exp_t));
  if (tmp == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  for (j = 0; j < elsz; ++j) {
    evl[j]  = tmp + (unsigned long)(j*nvars);
  }
}

static void free_global_hash_table(
    void
    )
{
  if (hmap) {
    free(hmap);
    hmap = NULL;
  }
  if (hd) {
    free(hd);
    hd  = NULL;
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
    /* note: memory is allocated as one big block,
     *       so freeing ev[0] is enough */
    free(ev[0]);
    free(ev);
    ev  = NULL;
  }
  if (etmp) {
    free(etmp);
    etmp  = NULL;
  }
  fc    = 0;
  nvars = 0;
  esz   = 0;
  eld   = 0;
  hsz   = 0;
}

static void free_local_hash_table(
    void
    )
{
  if (hmapl) {
    free(hmapl);
    hmapl  = NULL;
  }
  if (hdl) {
    free(hdl);
    hdl = NULL;
  }
  if (evl) {
    /* note: memory is allocated as one big block,
     *       so freeing evl[0] is enough */
    free(evl[0]);
    free(evl);
    evl = NULL;
  }
  elsz  = 0;
  elld  = 0;
  hlsz  = 0;
}

/* we just double the hash table size */
static void enlarge_global_hash_table(
    void
    )
{
  hl_t i, j;
  val_t h, k;

  j   = esz; /* store old size */
  esz = 2 * esz;
  hd  = realloc(hd, (unsigned long)esz * sizeof(hd_t));
  memset(hd+eld, 0, (unsigned long)(esz-eld) * sizeof(hd_t));
  ev  = realloc(ev, (unsigned long)esz * sizeof(exp_t *));
  if (ev == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  /* note: memory is allocated as one big block, so reallocating
   *       memory from ev[0] is enough    */
  ev[0] = realloc(ev[0],
      (unsigned long)esz * (unsigned long)nvars * sizeof(exp_t));
  if (ev[0] == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  /* due to realloc we have to reset ALL ev entries,
   * memory might have been moved */
  for (i = 1; i < esz; ++i) {
    ev[i] = ev[0] + (unsigned long)(i*nvars);
  }

  hsz   = 2 * hsz;
  hmap  = realloc(hmap, (unsigned long)hsz * sizeof(hl_t));
  memset(hmap, 0, (unsigned long)hsz * sizeof(hl_t));

  /* reinsert known elements */
  for (i = 1; i < eld; ++i) {
    h = hd[i].val;

    /* probing */
    k = h;
    for (j = 0; j < hsz; ++j) {
      k = (k+j) & (hsz-1);
      if (hmap[k]) {
        continue;
      }
      hmap[k] = i;
      break;
    }
  }
  num_htenl++;
}

static inline sdm_t generate_short_divmask(
    const exp_t * const a
    )
{
  len_t i, j;
  int32_t res = 0;
  int32_t ctr = 0;

  for (i = 0; i < ndvars; ++i) {
    for (j = 0; j < bpv; ++j) {
      if ((sdm_t)a[i] >= dm[ctr]) {
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
  hl_t i;
  len_t j, steps;
  int32_t ctr = 0;

  deg_t *max_exp  = (deg_t *)malloc((unsigned long)ndvars * sizeof(deg_t));
  deg_t *min_exp  = (deg_t *)malloc((unsigned long)ndvars * sizeof(deg_t));

  exp_t *e  = ev[1];

  /* get initial values from first hash table entry */
  for (i = 0; i < ndvars; ++i) {
    max_exp[i]  = min_exp[i]  = e[i];
  }

  /* get maximal and minimal exponent element entries in hash table */
  for (i = 2; i < eld; ++i) {
    e = ev[i];
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
      dm[ctr++] = (sdm_t)steps++;
    }
  }

  /* initialize divmasks for elements already added to hash table */
  for (i = 1; i < eld; i++) {
    hd[i].sdm = generate_short_divmask(ev[i]);
  }

  free(max_exp);
  free(min_exp);
}

/* returns zero if a is not divisible by b, else 1 is returned */
static inline hl_t check_monomial_division(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;
  const len_t nv  = nvars;

  /* short divisor mask check */
  if (hd[b].sdm & ~hd[a].sdm) {
    return 0;
  }

  const exp_t *const ea = ev[a];
  const exp_t *const eb = ev[b];
  /* exponent check */
  if (ea[0] < eb[0]) {
    return 0;
  }
  i = nv & 1 ? 1 : 0;
  for (; i < nv; i += 2) {
    if (ea[i] < eb[i] || ea[i+1] < eb[i+1]) {
      return 0;
    }
  }
  return 1;
}

static inline hl_t check_monomial_division_local(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;
  const len_t nv  = nvars;

  /* short divisor mask check */
  if (hdl[b].sdm & ~hdl[a].sdm) {
    return 0;
  }

  const exp_t *const ea = evl[a];
  const exp_t *const eb = evl[b];
  /* exponent check */
  if (ea[0] < eb[0]) {
    return 0;
  }
  i = nv & 1 ? 1 : 0;
  for (; i < nv; i += 2) {
    if (ea[i] < eb[i] || ea[i+1] < eb[i+1]) {
      return 0;
    }
  }
  return 1;
}

static inline hl_t insert_in_global_hash_table(
    const exp_t *a
    )
{
  hl_t i, k, pos;
  len_t j;
  exp_t deg;
  exp_t *e;
  hd_t *d;
  val_t h = 0;

  /* generate hash value */
  for (j = 0; j < nvars; ++j) {
    h +=  rv[j] * a[j];
  }

  /* probing */
  k = h;
  for (i = 0; i < hsz; ++i) {
    k = (k+i) & (hsz-1);
    if (!hmap[k]) {
      break;
    }
    if (hd[hmap[k]].val != h) {
      continue;
    }
    if (memcmp(ev[hmap[k]], a, (unsigned long)nvars * sizeof(exp_t)) == 0) {
      return hmap[k];
    }
  }

  /* add element to hash table */
  hmap[k]  = pos = eld;
  e   = ev[pos];
  d   = hd + pos;
  deg = 0;
  for (j = 0; j < nvars; ++j) {
    e[j]  =   a[j];
    deg   +=  a[j];
  }
  d->deg  = deg;
  d->sdm  = generate_short_divmask(e);
  d->val  = h;

  eld++;
  if (eld >= esz) {
    enlarge_global_hash_table();
  }

  return pos;
}

static inline hl_t insert_in_global_hash_table_no_enlargement_check(
    const exp_t *a
    )
{
  hl_t i, k, pos;
  len_t j;
  exp_t deg;
  exp_t *e;
  hd_t *d;
  val_t h = 0;

  /* generate hash value */
  for (j = 0; j < nvars; ++j) {
    h +=  rv[j] * a[j];
  }

  /* probing */
  k = h;
  for (i = 0; i < hsz; ++i) {
    k = (k+i) & (hsz-1);
    if (!hmap[k]) {
      break;
    }
    if (hd[hmap[k]].val != h) {
      continue;
    }
    if (memcmp(ev[hmap[k]], a, (unsigned long)nvars * sizeof(exp_t)) == 0) {
      return hmap[k];
    }
  }

  /* add element to hash table */
  hmap[k]  = pos = eld;
  e   = ev[pos];
  d   = hd + pos;
  deg = 0;
  for (j = 0; j < nvars; ++j) {
    e[j]  =   a[j];
    deg   +=  a[j];
  }
  d->deg  = deg;
  d->sdm  = generate_short_divmask(e);
  d->val  = h;

  eld++;

  return pos;
}


static inline hl_t insert_in_local_hash_table(
    const exp_t *a
    )
{
  hl_t i, k, pos;
  len_t j;
  exp_t deg;
  exp_t *e;
  hd_t *d;
  val_t h = 0;

  /* generate hash value */
  for (j = 0; j < nvars; ++j) {
    h +=  rv[j] * a[j];
  }

  /* probing */
  k = h;
  for (i = 0; i < hlsz; ++i) {
    k = (k+i) & (hlsz-1);
    if (!hmapl[k]) {
      break;
    }
    if (hd[hmapl[k]].val != h) {
      continue;
    }
    if (memcmp(evl[hmapl[k]], a, (unsigned long)nvars * sizeof(exp_t)) == 0) {
      return hmapl[k];
    }
  }

  /* add element to hash table */
  hmapl[k]  = pos = elld;
  e   = evl[pos];
  d   = hdl + pos;
  deg = 0;
  for (j = 0; j < nvars; ++j) {
    e[j]  =   a[j];
    deg   +=  a[j];
  }
  d->deg  = deg;
  d->sdm  = generate_short_divmask(e);
  d->val  = h;

  elld++;

  return pos;
}

static inline void reset_local_hash_table(
    const len_t size
    )
{
  hl_t i, j;
  /* is there still enough space in the local table? */
  if (size >= (elsz-elld)) {
    j = elsz; /* store ol size */
    if (2*size >= hlsz) {
      while (2*size >= hlsz) {
        elsz  = 2 * elsz;
        hlsz  = 2 * hlsz;
      }
      hdl   = realloc(hdl, (unsigned long)elsz * sizeof(hd_t));
      evl  = realloc(evl, (unsigned long)elsz * sizeof(exp_t *));
      if (evl == NULL) {
        printf("Computation needs too much memory on this machine, \
            segmentation fault will follow.\n");
      }
      /* note: memory is allocated as one big block, so reallocating
      *       memory from evl[0] is enough    */
      evl[0]  = realloc(evl[0],
          (unsigned long)elsz * (unsigned long)nvars * sizeof(exp_t));
      if (evl[0] == NULL) {
        printf("Computation needs too much memory on this machine, \
            segmentation fault will follow.\n");
      }
      /* due to realloc we have to reset ALL evl entries, memory might be moved */
      for (i = 1; i < elsz; ++i) {
        evl[i] = evl[0] + (unsigned long)(i*nvars);
      }
      /* for (i = j; i < elsz; ++i) {
       *   evl[i]  = (exp_t *)malloc((unsigned long)nvars * sizeof(exp_t));
       * } */
      hmapl = realloc(hmapl, (unsigned long)hlsz * sizeof(hl_t));
    }
    memset(hdl, 0, (unsigned long)elsz * sizeof(hd_t));
    memset(hmapl, 0, (unsigned long)hlsz * sizeof(hl_t));

    elld  = 1;
  }
}

static inline hl_t insert_in_global_hash_table_product_special(
    const val_t h1,
    const deg_t deg,
    const exp_t * const ea,
    const hl_t b
    )
{
  hl_t i, k, pos;
  len_t j;
  exp_t *n;
  hd_t *d;

  const val_t h   = h1 + hd[b].val;
  const exp_t * const eb = ev[b];

  n = ev[eld];
  for (j = 0; j < nvars; ++j) {
    n[j]  = ea[j] + eb[j];
  }
  /* probing */
  k = h;
  for (i = 0; i < hsz; ++i) {
    k = (k+i) & (hsz-1);
    if (!hmap[k]) {
      break;
    }
    if (hd[hmap[k]].val != h) {
      continue;
    }
    if (memcmp(n, ev[hmap[k]], (unsigned long)nvars * sizeof(exp_t)) == 0) {
      return hmap[k];
    }
  }

  /* add element to hash table */
  hmap[k] = pos = eld;
  d = hd + eld;
  d->deg  = deg + hd[b].deg;
  d->sdm  = generate_short_divmask(n);
  d->val  = h;

  eld++;
  return pos;
}

/* note that the product insertion, i.e. monomial x polynomial
 * is only needed for the local hash table. in the global one we
 * only add the basis elements, i.e. no multiplication is applied. */
static inline hl_t insert_in_global_hash_table_product(
    const hl_t a,
    const hl_t b
    )
{
  hl_t i, k, pos;
  len_t j;
  exp_t *n;
  hd_t *d;

  const val_t h   = hd[a].val + hd[b].val;

  /* printf("hash %d\n", h); */
  n = ev[eld];
  for (j = 0; j < nvars; ++j) {
    n[j]  = ev[a][j] + ev[b][j];
  }
  k = h;
  for (i = 0; i < hsz; ++i) {
    k = (k+i) & (hsz-1);
    if (!hmap[k]) {
      break;
    }
    if (hd[hmap[k]].val != h) {
      continue;
    }
    if (memcmp(n, ev[hmap[k]], (unsigned long)nvars * sizeof(exp_t)) == 0) {
      return hmap[k];
    }
  }

  /* add element to hash table */
  hmap[k] = pos = eld;
  d = hd + eld;
  d->deg  = hd[a].deg + hd[b].deg;
  d->sdm  = generate_short_divmask(n);
  d->val  = h;

  eld++;
  return pos;
}

static void reset_global_hash_table(
    void
    )
{
  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  len_t i, j;
  hl_t k;
  exp_t *e;
  val_t *b;

  exp_t **oev  = ev;
  ev  = calloc((unsigned long)esz, sizeof(exp_t *));
  if (ev == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  exp_t *tmp  = (exp_t *)malloc(
      (unsigned long)nvars * (unsigned long)esz * sizeof(exp_t));
  if (tmp == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  for (k = 0; k < esz; ++k) {
    ev[k]  = tmp + k*nvars;
  }
  eld = 1;
  memset(hmap, 0, (unsigned long)hsz * sizeof(hl_t));
  memset(hd, 0, (unsigned long)esz * sizeof(hd_t));

  /* reinsert known elements */
  for (i = 0; i < bload; ++i) {
    b = (val_t *)((long)bs[i] & bmask);
    for (j = 2; j < b[0]; j = j+2) {
      e = oev[b[j]];
      b[j]  = insert_in_global_hash_table_no_enlargement_check(e);
    }
  }
  for (i = 0; i < pload; ++i) {
    e = oev[ps[i].lcm];
    ps[i].lcm = insert_in_global_hash_table_no_enlargement_check(e);
  }
  /* note: all memory is allocated as a big block, so it is
   *       enough to free oev[0].       */
  free(oev[0]);
  free(oev);

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  rght_ctime  +=  ct1 - ct0;
  rght_rtime  +=  rt1 - rt0;
}

/* we can check equality of lcm and multiplication of two monomials
 * by their hash values. If both hash values are NOT the same, then
 * the corresponding exponent vectors CANNOT be the same. */
static inline int lcm_equals_multiplication(
    const hl_t a,
    const hl_t b,
    const hl_t lcm
    )
{
  const hd_t ha = hd[a];
  const hd_t hb = hd[b];
  const hd_t hl = hdl[lcm];

  if (hl.deg != ha.deg + hb.deg) {
    return 0;
  }
  if (hl.val != ha.val + hb.val) {
    return 0;
  } else {
    /* both have the same degree and the same hash value, either they
     * are the same or our hashing is broken resp. really bad */
    return 1;
  }
}

static inline hl_t get_lcm(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;

  /* exponents of basis elements, thus from global hash table */
  const exp_t * const ea = ev[a];
  const exp_t * const eb = ev[b];

  for (i = 0; i < nvars; ++i) {
    etmp[i]  = ea[i] < eb[i] ? eb[i] : ea[i];
  }
  /* goes into local hash table for spairs */
  return insert_in_local_hash_table(etmp);
}

static inline hl_t monomial_multiplication(
    const hl_t a,
    const hl_t b
    )
{
  return insert_in_global_hash_table_product(a, b);
}

/* we try monomial division including check if divisibility is
 * fulfilled. */
static inline hl_t monomial_division_with_check(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;

  const hd_t ha = hd[a];
  const hd_t hb = hd[b];
  /* short divisor mask check */
  if (hb.sdm & ~ha.sdm) {
    return 0;
  }

  const exp_t * const ea  = ev[a];
  const exp_t * const eb  = ev[b];
  etmp[0]  = ea[0] - eb[0];

  i = nvars & 1 ? 1 : 0;
  for (; i < nvars; i += 2) {
    if (ea[i] < eb[i] || ea[i+1] < eb[i+1]) {
      return 0;
    } else {
      etmp[i]    = ea[i] - eb[i];
      etmp[i+1]  = ea[i+1] - eb[i+1];
    }
  }
  return insert_in_global_hash_table(etmp);
}

/* it is assumed that b divides a, thus no tests for
 * divisibility at all */
static inline hl_t monomial_division_no_check(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;

  const exp_t * const ea  = ev[a];
  const exp_t * const eb  = ev[b];

  i = nvars & 1 ? 1 : 0;
  for (; i < nvars; i += 2) {
    etmp[i]    = ea[i]   - eb[i];
    etmp[i+1]  = ea[i+1] - eb[i+1];
  }
  etmp[0]  = ea[0] - eb[0];
  return insert_in_global_hash_table(etmp);
}

static inline val_t *multiplied_polynomial_to_matrix_row(
    const val_t hm,
    const deg_t deg,
    const exp_t * const em,
    const val_t *poly
    )
{
  len_t i;

  val_t *row  = (val_t *)malloc((unsigned long)poly[0] * sizeof(val_t));
  memcpy(row, poly, (unsigned long)poly[0] * sizeof(val_t));
  /* hash table product insertions appear only here:
   * we check for hash table enlargements first and then do the insertions
   * without further elargment checks there */
  if (eld+((poly[0]-2)/2) >= esz) {
    enlarge_global_hash_table();
  }
  for (i = 2; i < poly[1]; i += 2) {
    row[i]  = insert_in_global_hash_table_product_special(
                hm, deg, em, poly[i]);
  }
#if 1
  for (;i < poly[0]; i += 8) {
    row[i]    = insert_in_global_hash_table_product_special(
                  hm, deg, em, poly[i]);
    row[i+2]  = insert_in_global_hash_table_product_special(
                  hm, deg, em, poly[i+2]);
    row[i+4]  = insert_in_global_hash_table_product_special(
                  hm, deg, em, poly[i+4]);
    row[i+6]  = insert_in_global_hash_table_product_special(
                  hm, deg, em, poly[i+6]);
  }
#endif

  return row;
}
