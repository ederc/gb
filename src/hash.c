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

/* we have three different hash tables:
 * 1. one hash table for elements in the basis
 * 2. one hash table for the spairs during the update process
 * 3. one hash table for the multiplied elements during symbolic
 *    preprocessing */

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

static void initialize_basis_hash_table(
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

static void initialize_update_hash_table(
    void
    )
{
  hl_t j;

  /* generate map */
  husz  = (hl_t)pow(2, htes-5);
  humap = calloc((unsigned long)husz, sizeof(hl_t));

  /* generate exponent vector */
  eusz  = husz/2;
  /* keep first entry empty for faster divisibility checks */
  euld  = 1;
  hdu   = (hd_t *)calloc((unsigned long)eusz, sizeof(hd_t));
  evu   = (exp_t **)malloc((unsigned long)eusz * sizeof(exp_t *));
  if (evu == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  exp_t *tmp  = (exp_t *)malloc(
      (unsigned long)nvars * (unsigned long)eusz * sizeof(exp_t));
  if (tmp == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  for (j = 0; j < eusz; ++j) {
    evu[j]  = tmp + (unsigned long)(j*nvars);
  }
}

static void initialize_symbolic_hash_table(
    void
    )
{
  hl_t j;

  /* generate map */
  hssz  = (hl_t)pow(2, htes-2);
  hmaps = calloc((unsigned long)hssz, sizeof(hl_t));

  /* generate exponent vector */
  essz  = hssz/2;
  /* keep first entry empty for faster divisibility checks */
  esld  = 1;
  hds   = (hd_t *)calloc((unsigned long)essz, sizeof(hd_t));
  evs   = (exp_t **)malloc((unsigned long)essz * sizeof(exp_t *));
  if (evs == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  exp_t *tmp  = (exp_t *)malloc(
      (unsigned long)nvars * (unsigned long)essz * sizeof(exp_t));
  if (tmp == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  for (j = 0; j < essz; ++j) {
    evs[j]  = tmp + (unsigned long)(j*nvars);
  }
}

static void free_basis_hash_table(
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

static void free_update_hash_table(
    void
    )
{
  if (humap) {
    free(humap);
    humap  = NULL;
  }
  if (hdu) {
    free(hdu);
    hdu = NULL;
  }
  if (evu) {
    /* note: memory is allocated as one big block,
     *       so freeing evu[0] is enough */
    free(evu[0]);
    free(evu);
    evu = NULL;
  }
  eusz  = 0;
  euld  = 0;
  husz  = 0;
}

static void free_symbolic_hash_table(
    void
    )
{
  if (hmaps) {
    free(hmaps);
    hmaps  = NULL;
  }
  if (hds) {
    free(hds);
    hds = NULL;
  }
  if (evs) {
    /* note: memory is allocated as one big block,
     *       so freeing evl[0] is enough */
    free(evs[0]);
    free(evs);
    evs = NULL;
  }
  essz  = 0;
  esld  = 0;
  hssz  = 0;
}

/* we just double the hash table size */
static void enlarge_basis_hash_table(
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
}

static void enlarge_symbolic_hash_table(
    void
    )
{
  hl_t i, j;
  val_t h, k;

  j   = essz; /* store old size */
  essz = 2 * essz;
  hds  = realloc(hds, (unsigned long)essz * sizeof(hd_t));
  memset(hds+esld, 0, (unsigned long)(essz-esld) * sizeof(hd_t));
  evs  = realloc(evs, (unsigned long)essz * sizeof(exp_t *));
  if (evs == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  /* note: memory is allocated as one big block, so reallocating
   *       memory from ev[0] is enough    */
  evs[0] = realloc(evs[0],
      (unsigned long)essz * (unsigned long)nvars * sizeof(exp_t));
  if (evs[0] == NULL) {
    printf("Computation needs too much memory on this machine, \
        segmentation fault will follow.\n");
  }
  /* due to realloc we have to reset ALL ev entries,
   * memory might have been moved */
  for (i = 1; i < essz; ++i) {
    evs[i] = evs[0] + (unsigned long)(i*nvars);
  }

  hssz   = 2 * hssz;
  hmaps  = realloc(hmaps, (unsigned long)hssz * sizeof(hl_t));
  memset(hmaps, 0, (unsigned long)hssz * sizeof(hl_t));

  /* reinsert known elements */
  for (i = 1; i < esld; ++i) {
    h = hds[i].val;

    /* probing */
    k = h;
    for (j = 0; j < hssz; ++j) {
      k = (k+j) & (hssz-1);
      if (hmaps[k]) {
        continue;
      }
      hmaps[k] = i;
      break;
    }
  }
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

static inline hl_t check_monomial_division_update(
    const hl_t a,
    const hl_t b
    )
{
  len_t i;
  const len_t nv  = nvars;

  /* short divisor mask check */
  if (hdu[b].sdm & ~hdu[a].sdm) {
    return 0;
  }

  const exp_t *const ea = evu[a];
  const exp_t *const eb = evu[b];
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

static inline hl_t insert_in_basis_hash_table(
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
  i = 0;
restart:
  for (; i < hsz; ++i) {
    k = (k+i) & (hsz-1);
    const hl_t hm = hmap[k];
    if (!hm) {
      break;
    }
    if (hd[hm].val != h) {
      continue;
    }
    for (j = 0; j < nvars; ++j) {
        const exp_t * const ehm = ev[hm];
        if (a[j] != ehm[j]) {
            i++;
            goto restart;
        }
    }
    return hm;
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
    enlarge_basis_hash_table();
  }

  return pos;
}

static inline hl_t insert_in_basis_hash_table_no_enlargement_check(
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
  i = 0;
restart:
  for (; i < hsz; ++i) {
    k = (k+i) & (hsz-1);
    const hl_t hm = hmap[k];
    if (!hm) {
      break;
    }
    if (hd[hm].val != h) {
      continue;
    }
    for (j = 0; j < nvars; ++j) {
        const exp_t * const ehm = ev[hm];
        if (a[j] != ehm[j]) {
            i++;
            goto restart;
        }
    }
    return hm;
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


static inline hl_t insert_in_update_hash_table(
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
  i = 0;
restart:
  for (; i < husz; ++i) {
    k = (k+i) & (husz-1);
    const hl_t hm = humap[k];
    if (!hm) {
      break;
    }
    if (hdu[hm].val != h) {
      continue;
    }
    for (j = 0; j < nvars; ++j) {
        const exp_t * const ehm = evu[hm];
        if (a[j] != ehm[j]) {
            i++;
            goto restart;
        }
    }
    return hm;
  }

  /* add element to hash table */
  humap[k]  = pos = euld;
  e   = evu[pos];
  d   = hdu + pos;
  deg = 0;
  for (j = 0; j < nvars; ++j) {
    e[j]  =   a[j];
    deg   +=  a[j];
  }
  d->deg  = deg;
  d->sdm  = generate_short_divmask(e);
  d->val  = h;

  euld++;

  return pos;
}

static inline void reset_update_hash_table(
    const len_t size
    )
{
  hl_t i, j;
  /* is there still enough space in the local table? */
  if (size >= (eusz-euld)) {
    j = eusz; /* store ol size */
    if (2*size >= husz) {
      while (2*size >= husz) {
        eusz  = 2 * eusz;
        husz  = 2 * husz;
      }
      hdu   = realloc(hdu, (unsigned long)eusz * sizeof(hd_t));
      evu  = realloc(evu, (unsigned long)eusz * sizeof(exp_t *));
      if (evu == NULL) {
        printf("Computation needs too much memory on this machine, \
            segmentation fault will follow.\n");
      }
      /* note: memory is allocated as one big block, so reallocating
      *       memory from evl[0] is enough    */
      evu[0]  = realloc(evu[0],
          (unsigned long)eusz * (unsigned long)nvars * sizeof(exp_t));
      if (evu[0] == NULL) {
        printf("Computation needs too much memory on this machine, \
            segmentation fault will follow.\n");
      }
      /* due to realloc we have to reset ALL evl entries, memory might be moved */
      for (i = 1; i < eusz; ++i) {
        evu[i] = evu[0] + (unsigned long)(i*nvars);
      }
      /* for (i = j; i < elsz; ++i) {
       *   evl[i]  = (exp_t *)malloc((unsigned long)nvars * sizeof(exp_t));
       * } */
      humap = realloc(humap, (unsigned long)husz * sizeof(hl_t));
    }
    memset(hdu, 0, (unsigned long)eusz * sizeof(hd_t));
    memset(humap, 0, (unsigned long)husz * sizeof(hl_t));

    euld  = 1;
  }
}

static inline void reset_symbolic_hash_table(
        void
    )
{
    /* is there still enough space in the local table? */
    memset(hds, 0, (unsigned long)essz * sizeof(hd_t));
    memset(hmaps, 0, (unsigned long)hssz * sizeof(hl_t));

    esld  = 1;
}

static inline hl_t insert_in_symbolic_hash_table_product_special(
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

  /* printf("b %d | bload %d\n", b, bload); */
  const val_t h   = h1 + hd[b].val;
  const exp_t * const eb = ev[b];

  /* printf("esld %d / %d essz\n", esld, essz); */
  n = evs[esld];
  for (j = 0; j < nvars; ++j) {
    n[j]  = ea[j] + eb[j];
  }
  k = h;
  i = 0;
restart:
  for (; i < hssz; ++i) {
    k = (k+i) & (hssz-1);
    const hl_t hm  = hmaps[k];
    if (!hm) {
      break;
    }
    if (hds[hm].val != h) {
      continue;
    }
    for (j = 0; j < nvars; ++j) {
        const exp_t * const ehm = evs[hm];
        if (n[j] != ehm[j]) {
            i++;
            goto restart;
        }
    }
    return hm;
  }

  /* add element to hash table */
  hmaps[k] = pos = esld;
  d = hds + esld;
  d->deg  = deg + hd[b].deg;
  d->sdm  = generate_short_divmask(n);
  d->val  = h;

  esld++;
  return pos;
}

static inline hl_t insert_in_basis_hash_table_product_special(
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

  /* printf("b %d | bload %d\n", b, bload); */
  const val_t h   = h1 + hd[b].val;
  const exp_t * const eb = ev[b];

  n = ev[eld];
  for (j = 0; j < nvars; ++j) {
    n[j]  = ea[j] + eb[j];
  }
  /* probing */
  k = h;
  i = 0;
restart:
  for (; i < hsz; ++i) {
    k = (k+i) & (hsz-1);
    const hl_t hm = hmap[k];
    if (!hm) {
      break;
    }
    if (hd[hm].val != h) {
      continue;
    }
    for (j = 0; j < nvars; ++j) {
        const exp_t * const ehm = ev[hm];
        if (n[j] != ehm[j]) {
            i++;
            goto restart;
        }
    }
    return hm;
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
 * is only needed for the update hash table. in the basis one we
 * only add the basis elements, i.e. no multiplication is applied. */
static inline hl_t insert_in_basis_hash_table_product(
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
  i = 0;
restart:
  for (; i < hsz; ++i) {
    k = (k+i) & (hsz-1);
    const hl_t hm = hmap[k];
    if (!hmap[k]) {
      break;
    }
    if (hd[hmap[k]].val != h) {
      continue;
    }
    for (j = 0; j < nvars; ++j) {
        const exp_t * const ehm = ev[hm];
        if (n[j] != ehm[j]) {
            i++;
            goto restart;
        }
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

/* this function is deprecated, we no longer need to
 * reset the basis hash table */
static void reset_basis_hash_table(
    ps_t *psl,
    stat_t *st
    )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j;
    hl_t k;
    exp_t *e;
    dt_t *b;

    spair_t *ps = psl->p;

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
        const len_t os  = gbdt[i][1];
        const len_t len = gbdt[i][2];
        b = gbdt[i] + 3;
        for (j = 0; j < os; ++j) {
            e = oev[b[j]];
            b[j]  = insert_in_basis_hash_table_no_enlargement_check(e);
        }
        for (; j < len; j += 4) {
            e       = oev[b[j]];
            b[j]    = insert_in_basis_hash_table_no_enlargement_check(e);
            e       = oev[b[j+1]];
            b[j+1]  = insert_in_basis_hash_table_no_enlargement_check(e);
            e       = oev[b[j+2]];
            b[j+2]  = insert_in_basis_hash_table_no_enlargement_check(e);
            e       = oev[b[j+3]];
            b[j+3]  = insert_in_basis_hash_table_no_enlargement_check(e);
        }
    }
    const len_t pld = psl->ld;
    for (i = 0; i < pld; ++i) {
        e = oev[ps[i].lcm];
        ps[i].lcm = insert_in_basis_hash_table_no_enlargement_check(e);
    }
    /* note: all memory is allocated as a big block, so it is
     *       enough to free oev[0].       */
    free(oev[0]);
    free(oev);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->rht_ctime  +=  ct1 - ct0;
    st->rht_rtime  +=  rt1 - rt0;
}

static inline int prime_monomials(
    const hl_t a,
    const hl_t b
    )
{
    len_t i;

    const exp_t * const ea = ev[a];
    const exp_t * const eb = ev[b];

    const len_t nv  = nvars;
    const len_t os  = nvars % 2;
    if (ea[0] != 0 && eb[0] != 0) {
        return 0;
    }
    for (i = os; i < nv; ++i) {
        if ((ea[i] != 0 && eb[i] != 0) || (ea[i+1] != 0 && eb[i+1] != 0)) {
            return 0;
        }
    }
    return 1;
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
  const hd_t hl = hdu[lcm];

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

  /* exponents of basis elements, thus from basis hash table */
  const exp_t * const ea = ev[a];
  const exp_t * const eb = ev[b];

  for (i = 0; i < nvars; ++i) {
    etmp[i]  = ea[i] < eb[i] ? eb[i] : ea[i];
  }
  return insert_in_update_hash_table(etmp);
}

static inline hl_t monomial_multiplication(
    const hl_t a,
    const hl_t b
    )
{
  return insert_in_basis_hash_table_product(a, b);
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
  return insert_in_basis_hash_table(etmp);
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
  return insert_in_basis_hash_table(etmp);
}

static inline dt_t *multiplied_polynomial_to_matrix_row(
    const val_t hm,
    const deg_t deg,
    const exp_t * const em,
    const dt_t *poly
    )
{
  len_t i;

  dt_t *row = (dt_t *)malloc((unsigned long)(poly[2]+3) * sizeof(dt_t));
  row[0]    = poly[0];
  row[1]    = poly[1];
  row[2]    = poly[2];
  /* hash table product insertions appear only here:
   * we check for hash table enlargements first and then do the insertions
   * without further elargment checks there */
  while (esld+poly[2] >= essz) {
    enlarge_symbolic_hash_table();
  }
  const len_t os  = row[1];
  const len_t len = row[2];
  dt_t *rs  = row + 3;
  const dt_t * const ps  = poly + 3;
  /* printf("poly[1] %d | poly[2] %d\n", poly[1], poly[2]); */
  for (i = 0; i < os; ++i) {
    rs[i] = insert_in_symbolic_hash_table_product_special(
                hm, deg, em, ps[i]);
  }
  for (;i < len; i += 4) {
    rs[i]   = insert_in_symbolic_hash_table_product_special(
                  hm, deg, em, ps[i]);
    rs[i+1] = insert_in_symbolic_hash_table_product_special(
                  hm, deg, em, ps[i+1]);
    rs[i+2] = insert_in_symbolic_hash_table_product_special(
                  hm, deg, em, ps[i+2]);
    rs[i+3] = insert_in_symbolic_hash_table_product_special(
                  hm, deg, em, ps[i+3]);
  }

  return row;
}

/* deprecated once we use an own hash table for symbolic preprocessing data */
#if 0
static inline hl_t *reset_idx_in_basis_hash_table_and_free_hcm(
        hl_t *hcm
        )
{
    len_t i;

    for (i = 0; i < ncols; ++i) {
        hd[hcm[i]].idx  = 0;
    }
    free(hcm);

    return NULL;
}
#endif
