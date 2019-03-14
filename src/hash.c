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
 * 1. one hash table for elements in the basis (bht)
 * 2. one hash table for the spairs during the update process (uht)
 * 3. one hash table for the multiplied elements during symbolic
 *    preprocessing (sht) */

/* The idea of the structure of the hash table is taken from an
 * implementation by Roman Pearce and Michael Monagan in Maple. */

static val_t pseudo_random_number_generator(
    uint32_t *seed
    )
{
    uint32_t rseed  = *seed;
    rseed ^=  (rseed << 13);
    rseed ^=  (rseed >> 17);
    rseed ^=  (rseed << 5);
    *seed =   rseed;
    return (val_t)rseed;
}

static ht_t *initialize_basis_hash_table(
    const stat_t *st
    )
{
    len_t i;
    hl_t j;

    const len_t nv  = st->nvars;

    ht_t *ht  = (ht_t *)malloc(sizeof(ht_t));
    ht->nv    = nv;
    /* generate map */
    ht->bpv = (len_t)((CHAR_BIT * sizeof(sdm_t)) / (unsigned long)nv);
    if (ht->bpv == 0) {
        ht->bpv++;
    }
    ht->ndv = (unsigned long)nv < (CHAR_BIT * sizeof(sdm_t)) ?
        nv : (len_t)((CHAR_BIT * sizeof(sdm_t)));
    ht->hsz   = (hl_t)pow(2, st->init_hts);
    ht->hmap  = calloc((unsigned long)ht->hsz, sizeof(hl_t));

    /* generate divmask map */
    ht->dm  = (sdm_t *)calloc(
            (unsigned long)(ht->ndv * ht->bpv), sizeof(sdm_t));

    /* generate random values */
    ht->rsd = 2463534242;
    ht->rn  = calloc((unsigned long)nv, sizeof(val_t));
    for (i = nv; i > 0; --i) {
        /* random values should not be zero */
        ht->rn[i-1] = pseudo_random_number_generator(&(ht->rsd)) | 1;
    }
    /* generate exponent vector */
    ht->esz = ht->hsz/2;
    /* keep first entry empty for faster divisibility checks */
    ht->eld = 1;
    ht->hd  = (hd_t *)calloc((unsigned long)ht->esz, sizeof(hd_t));
    ht->ev  = (exp_t **)malloc((unsigned long)ht->esz * sizeof(exp_t *));
    if (ht->ev == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    exp_t *tmp  = (exp_t *)malloc(
            (unsigned long)ht->nv * (unsigned long)ht->esz * sizeof(exp_t));
    if (tmp == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    const hl_t esz  = ht->esz;
    for (j = 0; j < esz; ++j) {
        ev[j]  = tmp + (unsigned long)(j*nv);
    }
    return ht;
}

static ht_t *initialize_secondary_hash_table(
    const ht_t *bht,
    const stat_t *st
    )
{
    hl_t j;
    const hl_t nv = bht->nv;

    ht_t *ht  = (ht_t *)malloc(sizeof(ht_t)); 
    ht->nv    = nv;

    /* generate map */
    ht->hsz   = (hl_t)pow(2, st->init_hts-5);
    ht->hmap  = calloc((unsigned long)ht->hsz, sizeof(hl_t));

    /* divisor mask and hash value seeds from basis hash table */
    ht->dm  = bht->dm;
    ht->rn  = bht->rn;

    /* generate exponent vector */
    ht->esz = ht->hsz/2;
    /* keep first entry empty for faster divisibility checks */
    ht->eld = 1;
    ht->hd  = (hd_t *)calloc((unsigned long)ht->esz, sizeof(hd_t));
    ht->ev  = (exp_t **)malloc((unsigned long)ht->esz * sizeof(exp_t *));
    if (ht->ev == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    exp_t *tmp  = (exp_t *)malloc(
            (unsigned long)nv * (unsigned long)eusz * sizeof(exp_t));
    if (tmp == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    const hl_t esz  = ht->esz;
    for (j = 0; j < esz; ++j) {
        evu[j]  = tmp + (unsigned long)(j*nv);
    }
    return ht;
}

static void free_divmask(
        ht_t *ht
        )
{
    if (ht->dm) {
        free(ht->dm);
        ht->dm  = NULL;
    }
}

static void free_random_numbers(
        ht_t *ht
        )
{
    if (ht->rn) {
        free(ht->rn);
        ht->rn  = NULL;
    }
}

static void free_hash_table(
        ht_t **htp
    )
{
    ht_t *ht  = *htp;
  if (ht->hmap) {
    free(ht->hmap);
    ht->hmap = NULL;
  }
  if (ht->hd) {
    free(ht->hd);
    ht->hd  = NULL;
  }
  if (ht->ev) {
    /* note: memory is allocated as one big block,
     *       so freeing ev[0] is enough */
    free(ht->ev[0]);
    free(ht->ev);
    ht->ev  = NULL;
  }
  free(ht);
  ht    = NULL;
  *htp  = ht;
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
    const exp_t * const a,
    const ht_t *ht
    )
{
  len_t i, j;
  int32_t res = 0;
  int32_t ctr = 0;
  const len_t ndv = ht->ndv;
  const len_t bpv = ht->bpv;

  for (i = 0; i < ndv; ++i) {
    for (j = 0; j < bpv; ++j) {
      if ((sdm_t)a[i] >= ht->dm[ctr]) {
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
    ht_t *ht
    )
{
  hl_t i;
  len_t j, steps;
  int32_t ctr = 0;
  exp_t **ev  = ht->ev;

  deg_t *max_exp  = (deg_t *)malloc((unsigned long)ndvars * sizeof(deg_t));
  deg_t *min_exp  = (deg_t *)malloc((unsigned long)ndvars * sizeof(deg_t));

  exp_t *e  = ev[1];

  /* get initial values from first hash table entry */
  for (i = 0; i < ht->ndv; ++i) {
    max_exp[i]  = min_exp[i]  = e[i];
  }

  /* get maximal and minimal exponent element entries in hash table */
  for (i = 2; i < ht->eld; ++i) {
    e = ev[i];
    for (j = 0; j < ht->ndv; ++j) {
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
  for (i = 0; i < ht->ndv; ++i) {
    steps = (max_exp[i] - min_exp[i]) / ht->bpv;
    if (steps == 0)
      steps++;
    for (j = 0; j < ht->bpv; ++j) {
      ht->dm[ctr++] = (sdm_t)steps++;
    }
  }

  /* initialize divmasks for elements already added to hash table */
  for (i = 1; i < ht->eld; i++) {
    ht->hd[i].sdm = generate_short_divmask(ev[i], ht);
  }

  free(max_exp);
  free(min_exp);
}

/* returns zero if a is not divisible by b, else 1 is returned */
static inline hl_t check_monomial_division(
    const hl_t a,
    const hl_t b,
    const ht_t *ht
    )
{
  len_t i;
  const len_t nv  = ht->nv;

  /* short divisor mask check */
  if (ht->hd[b].sdm & ~ht->hd[a].sdm) {
    return 0;
  }

  const exp_t *const ea = ht->ev[a];
  const exp_t *const eb = ht->ev[b];
  /* exponent check */
  for (i = nv-1; i > 0; i -= 2) {
    if (ea[i] < eb[i] || ea[i-1] < eb[i-1]) {
      return 0;
    }
  }
  if (ea[0] < eb[0]) {
    return 0;
  }
  return 1;
}

static inline void check_monomial_division_in_update(
    hl_t *a,
    const len_t start,
    const len_t end,
    const hl_t b,
    const ht_t *ht
    )
{
    len_t i, j;
    const len_t nv  = ht->nv;

    const sdm_t sb        = ht->hd[b].sdm;
    const exp_t *const eb = ht->ev[b];
    /* pairs are sorted, we only have to search entries
     * above the starting point */
        j = start+1;
restart:
    for (; j < end; ++j) {
        if (a[j] < 0) {
            continue;
        }
        /* short divisor mask check */
        if (~ht->hd[a[j]].sdm & sb) {
            continue;
        }
        const exp_t *const ea = ht->ev[a[j]];
        /* exponent check */
        for (i = nv-1; i > 0; i -= 2) {
            if (ea[i] < eb[i] || ea[i-1] < eb[i-1]) {
                j++;
                goto restart;
            }
        }
        if (ea[0] < eb[0]) {
            continue;
        }
        a[j]  = -1;
    }
}

static inline hl_t insert_in_basis_hash_table(
    const exp_t *a
    )
{
  hl_t i, k, pos;
  len_t j;
  deg_t deg;
  exp_t *e;
  hd_t *d;
  val_t h = 0;

  const len_t nv  = nvars;

  /* generate hash value */
  for (j = 0; j < nv; ++j) {
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
    const exp_t * const ehm = ev[hm];
    for (j = nv-1; j > 0; j -= 2) {
        if (a[j] != ehm[j] || a[j-1] != ehm[j-1]) {
            i++;
            goto restart;
        }
    }
    if (a[0] != ehm[0]) {
        i++;
        goto restart;
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

static inline hl_t insert_in_hash_table(
    const exp_t *a,
    ht_t *ht
    )
{
  hl_t i, k, pos;
  len_t j;
  deg_t deg;
  exp_t *e;
  hd_t *d;
  val_t h = 0;
  const len_t nv  = ht->nv;
  const hl_t hsz  = ht->hsz;

  /* generate hash value */
  for (j = 0; j < nv; ++j) {
    h +=  rv[j] * a[j];
  }

  /* probing */
  k = h;
  i = 0;
restart:
  for (; i < hsz; ++i) {
    k = (k+i) & (hsz-1);
    const hl_t hm = ht->hmap[k];
    if (!hm) {
      break;
    }
    if (ht->hd[hm].val != h) {
      continue;
    }
    const exp_t * const ehm = ht->ev[hm];
    for (j = nv-1; j > 0; j -= 2) {
        if (a[j] != ehm[j] || a[j-1] != ehm[j-1]) {
            i++;
            goto restart;
        }
    }
    if (a[0] != ehm[0]) {
        i++;
        goto restart;
    }
    return hm;
  }

  /* add element to hash table */
  ht->hmap[k]  = pos = eld;
  e   = ht->ev[pos];
  d   = ht->hd + pos;
  deg = 0;
  for (j = 0; j < nv; ++j) {
    e[j]  =   a[j];
    deg   +=  a[j];
  }
  d->deg  = deg;
  d->sdm  = generate_short_divmask(e, ht);
  d->val  = h;

  ht->eld++;

  return pos;
}


static inline hl_t insert_in_update_hash_table(
    const exp_t *a
    )
{
  hl_t i, k, pos;
  len_t j;
  deg_t deg;
  exp_t *e;
  hd_t *d;
  val_t h = 0;
  const len_t nv  = nvars;

  /* generate hash value */
  for (j = 0; j < nv; ++j) {
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
    const exp_t * const ehm = evu[hm];
    for (j = nv-1; j > 0; j -= 2) {
        if (a[j] != ehm[j] || a[j-1] != ehm[j-1]) {
            i++;
            goto restart;
        }
    }
    if (a[0] != ehm[0]) {
        i++;
        goto restart;
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

static inline void reinitialize_hash_table(
    ht_t *ht,
    const len_t size
    )
{
    hl_t i;
    /* is there still enough space in the local table? */
    if (size >= (ht->esz)) {
        while (size >= ht->esz) {
            ht->esz = 2 * ht->esz;
            ht->hsz = 2 * ht->hsz;
        }
        const hl_t esz  = ht->esz;
        const hl_t hsz  = ht->hsz;
        const len_t nv  = ht->nv;
        ht->hd  = realloc(ht->hd, (unsigned long)esz * sizeof(hd_t));
        ht->ev  = realloc(ht->ev, (unsigned long)esz * sizeof(exp_t *));
        if (ht->ev == NULL) {
            printf("Computation needs too much memory on this machine, \
                    segmentation fault will follow.\n");
        }
        /* note: memory is allocated as one big block, so reallocating
         *       memory from evl[0] is enough    */
        ht->ev[0]  = realloc(ht->ev[0],
                (unsigned long)esz * (unsigned long)nv * sizeof(exp_t));
        if (evu[0] == NULL) {
            printf("Computation needs too much memory on this machine, \
                    segmentation fault will follow.\n");
        }
        /* due to realloc we have to reset ALL evl entries, memory might be moved */
        for (i = 1; i < esz; ++i) {
            ht->ev[i] = ht->ev[0] + (unsigned long)(i*nv);
        }
        ht->hmap  = realloc(ht->hmap, (unsigned long)hsz * sizeof(hl_t));
    }
    memset(ht->hd, 0, (unsigned long)esz * sizeof(hd_t));
    memset(ht->hmap, 0, (unsigned long)hsz * sizeof(hl_t));

    ht->eld  = 1;
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

static inline int prime_monomials(
    const hl_t a,
    const hl_t b,
    const ht_t *ht
    )
{
    len_t i;

    const exp_t * const ea = ht->ev[a];
    const exp_t * const eb = ht->ev[b];

    const len_t nv  = ht->nv;
    for (i = nv-1; i > 0; i -= 2) {
        if ((ea[i] != 0 && eb[i] != 0) || (ea[i-1] != 0 && eb[i-1] != 0)) {
            return 0;
        }
    }
    if (ea[0] != 0 && eb[0] != 0) {
        return 0;
    }
    return 1;
}

static inline void insert_plcms_in_basis_hash_table(
    ps_t *psl,
    spair_t *pp,
    ht_t *bht,
    const ht_t *uht,
    const hl_t * const lcms,
    const len_t start,
    const len_t end
    )
{
    hl_t i, k, pos;
    len_t j, l, m;
    hd_t *d;

    spair_t *ps     = psl->p;
    const len_t nv  = bht->nv;
    const hl_t hsz  = bht->hsz;
    m = start;
    l = 0;
letsgo:
    for (; l < end; ++l) {
        if (lcms[l] < 0) {
            continue;
        }
        if (prime_monomials(gbdt[pp[l].gen1][3], gbdt[pp[0].gen2][3], bht)) {
            continue;
        }
        ps[m] = pp[l];
        const val_t h = uht->hd[lcms[l]].val;
        memcpy(bht->ev[bht->eld], uht->ev[lcms[l]],
                (unsigned long)nv * sizeof(exp_t));
        const exp_t * const n = bht->ev[eld];
        k = h;
        i = 0;
restart:
        for (; i < hsz; ++i) {
            k = (k+i) & (hsz-1);
            const hl_t hm = bht->hmap[k];
            if (!hm) {
                break;
            }
            if (bht->hd[hm].val != h) {
                continue;
            }
            const exp_t * const ehm = bht->ev[hm];
            for (j = nv-1; j > 0; j -= 2) {
                if (n[j] != ehm[j] || n[j-1] != ehm[j-1]) {
                    i++;
                    goto restart;
                }
            }
            if (n[0] != ehm[0]) {
                i++;
                goto restart;
            }
            ps[m++].lcm = hm;
            l++;
            goto letsgo;
        }

        /* add element to hash table */
        bht->hmap[k] = pos = bht->eld;
        d = bht->hd + bht->eld;
        d->deg  = uht->hd[lcms[l]].deg;
        d->sdm  = uht->hd[lcms[l]].sdm;
        d->val  = h;

        bht->eld++;
        ps[m++].lcm =  pos;
    }
    psl->ld = m;
}

static inline void insert_in_basis_hash_table_pivots(
    hm_t *row,
    const hl_t * const hcm
    )
{
    hl_t i, k, pos;
    len_t j, l;
    hd_t *d;

    const len_t len = row[2]+3;
    const len_t nv  = nvars;
    l = 3;
letsgo:
    for (; l < len; ++l) {
        const val_t h = hds[hcm[row[l]]].val;
        memcpy(ev[eld], evs[hcm[row[l]]], (unsigned long)nv * sizeof(exp_t));
        const exp_t * const n = ev[eld];
        k = h;
        i = 0;
restart:
        for (; i < hsz; ++i) {
            k = (k+i) & (hsz-1);
            const hl_t hm  = hmap[k];
            if (!hm) {
                break;
            }
            if (hd[hm].val != h) {
                continue;
            }
            const exp_t * const ehm = ev[hm];
            for (j = nv-1; j > 0; j -= 2) {
                if (n[j] != ehm[j] || n[j-1] != ehm[j-1]) {
                    i++;
                    goto restart;
                }
            }
            if (n[0] != ehm[0]) {
                i++;
                goto restart;
            }
            row[l] = hm;
            l++;
            goto letsgo;
        }

        /* add element to hash table */
        hmap[k] = pos = eld;
        d = hd + eld;
        d->deg  = hds[hcm[row[l]]].deg;
        d->sdm  = hds[hcm[row[l]]].sdm;
        d->val  = h;

        eld++;
        row[l] =  pos;
    }
}

static inline void insert_multiplied_poly_in_hash_table(
    hm_t *row,
    const val_t h1,
    const deg_t deg,
    const exp_t * const ea,
    const hm_t * const b,
    const ht_t * const ht1,
    ht_t *ht2
    )
{
    hl_t i, k, pos;
    len_t j, l;
    exp_t *n;
    hd_t *d;

    const len_t len = b[2]+3;
    const len_t nv  = ht1->nv;

    exp_t * const *ev1      = ht1->ev;
    const hd_t * const hd1  = ht1->hd;
    
    exp_t **ev2     = ht2->ev;
    hd_t *hd2       = ht2->hd;
    const hl_t hsz2 = ht2->hsz;

    l = 3;
letsgo:
    for (; l < len; ++l) {
        /* printf("b %d | bload %d\n", b, bload); */
        const val_t h   = h1 + hd1[b[l]].val;
        const exp_t * const eb = ev1[b[l]];

        /* printf("esld %d / %d essz\n", esld, essz); */
        n = ev2[ht2->eld];
        for (j = 0; j < nv; ++j) {
            n[j]  = (exp_t)(ea[j] + eb[j]);
        }
        k = h;
        i = 0;
restart:
        for (; i < hsz2; ++i) {
            k = (k+i) & (hsz2-1);
            const hl_t hm  = ht2->hmap[k];
            if (!hm) {
                break;
            }
            if (hd2[hm].val != h) {
                continue;
            }
            const exp_t * const ehm = ev2[hm];
            for (j = nv-1; j > 0; j -= 2) {
                if (n[j] != ehm[j] || n[j-1] != ehm[j-1]) {
                    i++;
                    goto restart;
                }
            }
            if (n[0] != ehm[0]) {
                i++;
                goto restart;
            }
            row[l] = hm;
            l++;
            goto letsgo;
        }

        /* add element to hash table */
        ht2->hmap[k]  = pos = ht2->eld;
        d = hd2 + ht2->eld;
        d->deg  = deg + hd1[b[l]].deg;
        d->sdm  = generate_short_divmask(n, ht2);
        d->val  = h;

        ht2->eld++;
        row[l] =  pos;
    }
}

static inline void reinsert_in_hash_table(
    hm_t *row,
    exp_t * const *oev,
    ht_t *ht
    )
{
    hl_t i, k, pos;
    len_t j, l;
    exp_t *e;
    hd_t *d;
    val_t h;

    const len_t len = row[2]+3;
    const len_t nv  = ht->nv;
    const hl_t hsz  = ht->hsz;
    l = 3;
letsgo:
    for (; l < len; ++l) {
        const exp_t * const n = oev[row[l]];
        /* generate hash value */
        h = 0;
        for (j = 0; j < nv; ++j) {
            h +=  rv[j] * n[j];
        }
        k = h;
        i = 0;
restart:
        for (; i < hsz; ++i) {
            k = (k+i) & (hsz-1);
            const hl_t hm  = ht->hmap[k];
            if (!hm) {
                break;
            }
            if (ht->hd[hm].val != h) {
                continue;
            }
            const exp_t * const ehm = ht->ev[hm];
            for (j = nv-1; j > 0; j -= 2) {
                if (n[j] != ehm[j] || n[j-1] != ehm[j-1]) {
                    i++;
                    goto restart;
                }
            }
            if (n[0] != ehm[0]) {
                i++;
                goto restart;
            }
            row[l] = hm;
            l++;
            goto letsgo;
        }

        /* add element to hash table */
        ht->hmap[k] = pos = ht->eld;
        e = ht->ev[ht->eld];
        d = ht->hd + ht->eld;
        for (j = 0; j < nv; ++j) {
            e[j]    =   n[j];
            d->deg  +=  n[j];
        }
        d->sdm  = generate_short_divmask(e, ht);
        d->val  = h;

        ht->eld++;
        row[l] =  pos;
    }
}

static void reset_hash_table(
    ht_t *ht,
    bs_t *bs,
    ps_t *psl,
    stat_t *st
    )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i;
    hl_t k;
    exp_t *e;

    spair_t *ps = psl->p;
    exp_t **oev  = ht->ev;

    const len_t nv  = ht->nv;
    const hl_t esz  = ht->esz;
    const bl_t bld  = bs->ld;
    const len_t pld = psl->ld;

    ht->ev  = calloc((unsigned long)esz, sizeof(exp_t *));
    if (ht->ev == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    exp_t *tmp  = (exp_t *)malloc(
            (unsigned long)nv * (unsigned long)esz * sizeof(exp_t));
    if (tmp == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    for (k = 0; k < esz; ++k) {
        ht->ev[k]  = tmp + k*nv;
    }
    ht->eld = 1;
    memset(ht->hmap, 0, (unsigned long)hsz * sizeof(hl_t));
    memset(ht->hd, 0, (unsigned long)esz * sizeof(hd_t));

    /* reinsert known elements */
    for (i = 0; i < bld; ++i) {
        reinsert_in_hash_table(gbdt[i], oev, ht);
    }
    for (i = 0; i < pld; ++i) {
        e = oev[ps[i].lcm];
        ps[i].lcm = insert_in_hash_table(e, ht);
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

/* computes lcm of a and b from ht1 and inserts it in ht2 */
static inline hl_t get_lcm(
    const hl_t a,
    const hl_t b,
    const ht_t *ht1,
    ht_t *ht2
    )
{
    len_t i;

    /* exponents of basis elements, thus from basis hash table */
    const exp_t * const ea = ht1->ev[a];
    const exp_t * const eb = ht1->ev[b];
    exp_t *etmp = ht1->ev[0];
    const len_t nv  = ht1->nv;

    for (i = 0; i < nv; ++i) {
        etmp[i]  = ea[i] < eb[i] ? eb[i] : ea[i];
    }
    return insert_in_hash_table(etmp, ht2);
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
  etmp[0]  = (exp_t)(ea[0] - eb[0]);

  i = nvars & 1 ? 1 : 0;
  for (; i < nvars; i += 2) {
    if (ea[i] < eb[i] || ea[i+1] < eb[i+1]) {
      return 0;
    } else {
      etmp[i]    = (exp_t)(ea[i] - eb[i]);
      etmp[i+1]  = (exp_t)(ea[i+1] - eb[i+1]);
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
    etmp[i]    = (exp_t)(ea[i] - eb[i]);
    etmp[i+1]  = (exp_t)(ea[i+1] - eb[i+1]);
  }
  etmp[0]  = (exp_t)(ea[0] - eb[0]);
  return insert_in_basis_hash_table(etmp);
}

static inline hm_t *multiplied_poly_to_matrix_row(
    ht_t *sht,
    const ht_t *bht,
    const val_t hm,
    const deg_t deg,
    const exp_t * const em,
    const hm_t *poly
    )
{
  hm_t *row = (hm_t *)malloc((unsigned long)(poly[2]+3) * sizeof(hm_t));
  row[0]    = poly[0];
  row[1]    = poly[1];
  row[2]    = poly[2];
  /* hash table product insertions appear only here:
   * we check for hash table enlargements first and then do the insertions
   * without further elargment checks there */
  while (sht->eld+poly[2] >= sht->esz) {
    enlarge_symbolic_hash_table();
  }
  insert_multiplied_poly_in_hash_table(row, hm, deg, em, poly, bht, sht);

  return row;
}
