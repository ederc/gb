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
val_t pseudo_random_number_generator(
    ht_t *ht
    )
{
    ht->rseed ^= (ht->rseed << 13);
    ht->rseed ^= (ht->rseed >> 17);
    ht->rseed ^= (ht->rseed << 5);

    return (val_t)ht->rseed;
}

ht_t *initialize_global_hash_table(
    const stat_t *st
    )
{
    len_t i;
    hl_t j;

    ht_t *ht  = (ht_t *)calloc(1, sizeof(ht_t));
    ht->nv    = st->nr_vars;
    ht->rseed = 2463534242;

    /* generate map */
    ht->bpv = (len_t)((CHAR_BIT * sizeof(sdm_t)) / (unsigned long)nvars);
    if (ht->bpv == 0) {
        ht->bpv++;
    }
    ht->ndv = (unsigned long)ht->nv < (CHAR_BIT * sizeof(sdm_t)) ?
        ht->nv : (len_t)((CHAR_BIT * sizeof(sdm_t)));
    ht->hsz = (hl_t)pow(2, st->init_ht_sz);
    ht->map = calloc((unsigned long)ht->hsz, sizeof(hl_t));

    /* generate divmask map */
    ht->dm  = calloc((unsigned long)(ht->ndv * ht->bpv), sizeof(sdm_t));

    /* generate random values */
    ht->rv  = calloc((unsigned long)ht->nv, sizeof(val_t));
    for (i = 0 ; i < ht->nv; ++i) {
        /* random values should not be zero */
        rv[i] = pseudo_random_number_generator() | 1;
    }
    /* generate exponent vector */
    ht->esz = ht->hsz / 2;
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
    for (j = 0; j < ht->esz; ++j) {
        ht->ev[j]  = tmp + (unsigned long)(j * ht->nv);
    }

    return ht;
}

ht_t *initialize_local_hash_table(
    const stat_t *st,
    const ht_t *ght
    )
{
    hl_t j;

    ht_t *ht  = (ht_t *)calloc(1, sizeof(ht_t));
    ht->nv    = st->nr_vars;

    ht->rv  = ght->rv;

    /* generate map */
    ht->sz  = (hl_t)pow(2, st->init_ht_sz - 5);
    ht->map = calloc((unsigned long)ht->sz, sizeof(hl_t));

    /* generate exponent vector */
    ht->esz = ht->sz / 2;
    /* keep first entry empty for faster divisibility checks */
    ht->eld = 1;
    ht->hd  = (hd_t *)calloc((unsigned long)ht->esz, sizeof(hd_t));
    ht->ev  = (exp_t **)malloc((unsigned long)ht->esz * sizeof(exp_t *));
    if (ht->ev == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    exp_t *tmp  = (exp_t *)malloc(
            (unsigned long)ht->nv * (unsigned long)elsz * sizeof(exp_t));
    if (tmp == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    for (j = 0; j < ht->esz; ++j) {
        ht->ev[j]  = tmp + (unsigned long)(j*ht->nv);
    }

    return ht;
}

void free_hash_table(
    ht_t **htp
    )
{
    ht_t *ht  = *htp;
    if (ht->map) {
        free(ht->map);
    }
    if (ht->hd) {
        free(ht->hd);
    }
    if (ht->dm) {
        free(ht->dm);
    }
    if (ht->rv) {
        free(ht->rv);
    }
    if (ht->ev) {
        /* note: memory is allocated as one big block,
         *       so freeing ev[0] is enough */
        free(ht->ev[0]);
        free(ht->ev);
    }
    free(ht);
    ht    = NULL;
    *htp  = ht;
}

/* we just double the hash table size */
void enlarge_global_hash_table(
    ht_t *ht
    )
{
    hl_t i, j;
    val_t h, k;

    ht->esz = 2 * ht->esz;
    ht->hd  = realloc(ht->hd, (unsigned long)ht->esz * sizeof(hd_t));
    memset(ht->hd+ht->eld, 0, (unsigned long)(ht->esz-ht->eld) * sizeof(hd_t));
    ht->ev  = realloc(ht->ev, (unsigned long)ht->esz * sizeof(exp_t *));
    if (ht->ev == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    /* note: memory is allocated as one big block, so reallocating
     *       memory from ev[0] is enough    */
    ht->ev[0] = realloc(ht->ev[0],
            (unsigned long)ht->esz * (unsigned long)ht->nv * sizeof(exp_t));
    if (ht->ev[0] == NULL) {
        printf("Computation needs too much memory on this machine, \
                segmentation fault will follow.\n");
    }
    /* due to realloc we have to reset ALL ev entries,
     * memory might have been moved */
    for (i = 1; i < ht->esz; ++i) {
        ht->ev[i] = ht->ev[0] + (unsigned long)(i*nvars);
    }

    ht->hsz = 2 * ht->hsz;
    ht->map = realloc(ht->map, (unsigned long)ht->hsz * sizeof(hl_t));
    memset(ht->map, 0, (unsigned long)ht->hsz * sizeof(hl_t));

    const hl_t hsz  = ht->hsz;
    const hk_t eld  = ht->eld;
    /* reinsert known elements */
    for (i = 1; i < eld; ++i) {
        h = ht->hd[i].val;

        /* probing */
        k = h;
        for (j = 0; j < hsz; ++j) {
            k = (k+j) & (hsz-1);
            if (ht->map[k]) {
                continue;
            }
            ht->map[k]  = i;
            break;
        }
    }
}

inline sdm_t generate_short_divmask(
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
inline void calculate_divmask(
    ht_t *ht
    )
{
    hl_t i;
    len_t j, steps;
    int32_t ctr = 0;

    deg_t *max_exp  = (deg_t *)malloc((unsigned long)ht->ndv * sizeof(deg_t));
    deg_t *min_exp  = (deg_t *)malloc((unsigned long)ht->ndv * sizeof(deg_t));

    exp_t *e  = ht->ev[1];

    /* get initial values from first hash table entry */
    for (i = 0; i < ht->ndv; ++i) {
        max_exp[i]  = min_exp[i]  = e[i];
    }

    /* get maximal and minimal exponent element entries in hash table */
    for (i = 2; i < ht->eld; ++i) {
        e = ht->ev[i];
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
        ht->hd[i].sdm = generate_short_divmask(ht->ev[i], ht);
    }

    free(max_exp);
    free(min_exp);
}

/* returns zero if a is not divisible by b, else 1 is returned */
inline hl_t check_monomial_division(
    const hl_t a,
    const hl_t b,
    const ht_t *ht
    )
{
    len_t i;
    const len_t nv  = ht->nv;

    /* short divisor mask check */
    if (hd[b].sdm & ~hd[a].sdm) {
        return 0;
    }

    const exp_t *const ea = ht->ev[a];
    const exp_t *const eb = ht->ev[b];
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

inline hl_t insert_in_hash_table_no_enlargement_check(
    const exp_t *a,
    ht_t *ht
    )
{
    hl_t i, k, pos;
    len_t j;
    exp_t deg;
    exp_t *e;
    hd_t *d;
    val_t h = 0;

    /* generate hash value */
    for (j = 0; j < ht->nv; ++j) {
        h +=  ht->rv[j] * a[j];
    }

    const hl_t hsz  = ht->hsz;
    /* probing */
    k = h;
    for (i = 0; i < hsz; ++i) {
        k = (k+i) & (hsz-1);
        if (!ht->map[k]) {
            break;
        }
        if (ht->hd[ht->map[k]].val != h) {
            continue;
        }
        if (memcmp(ht->ev + ht->map[k], a,
                    (unsigned long)ht->nv * sizeof(exp_t)) == 0) {
            return ht->map[k];
        }
    }

    /* add element to hash table */
    ht->map[k]  = pos = ht->eld;
    e   = ht->ev + pos;
    d   = ht->hd + pos;
    deg = 0;
    for (j = 0; j < ht->nv; ++j) {
        e[j]  =   a[j];
        deg   +=  a[j];
    }
    d->deg  = deg;
    d->sdm  = generate_short_divmask(e, ht);
    d->val  = h;

    return pos;
}

inline hl_t insert_in_hash_table(
    const exp_t *a,
    ht_t *ht
    )
{
    hl_t pos  = insert_in_hash_table_no_enlargement_check(a, ht);

    ht->eld++;
    if (ht->eld >= ht->esz) {
        enlarge_global_hash_table(ht);
    }

    return pos;
}

inline void reset_local_hash_table(
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

inline hl_t insert_in_hash_table_product_special(
    const val_t h1,
    const deg_t deg,
    const exp_t * const ea,
    const hl_t b,
    ht_t *ht
    )
{
    hl_t i, k, pos;
    len_t j;
    exp_t *n;
    hd_t *d;

    /* printf("b %d | bload %d\n", b, bload); */
    const val_t h   = h1 + ht->hd[b].val;
    const exp_t * const eb = ht->ev[b];

    n = ht->ev[ht->eld];
    const hl_t hsz  = ht->hsz;
    for (j = 0; j < ht->nv; ++j) {
        n[j]  = ea[j] + eb[j];
    }
    /* probing */
    k = h;
    for (i = 0; i < hsz; ++i) {
        k = (k+i) & (hsz-1);
        if (!ht->map[k]) {
            break;
        }
        if (ht->hd[ht->map[k]].val != h) {
            continue;
        }
        if (memcmp(n, ht->ev[ht->map[k]],
                    (unsigned long)ht->nv * sizeof(exp_t)) == 0) {
            return ht->map[k];
        }
    }

    /* add element to hash table */
    ht->map[k]  = pos = ht->eld;
    d           = ht->hd + ht->eld;
    d->deg      = deg + ht->hd[b].deg;
    d->sdm      = generate_short_divmask(n, ht);
    d->val      = h;

    ht->eld++;
    return pos;
}

void reset_global_hash_table(
    ht_t *ht,
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

    exp_t **oev  = ht->ev;
    ht->ev  = calloc((unsigned long)ht->esz, sizeof(exp_t *));
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
    for (k = 0; k < ht->esz; ++k) {
        ht->ev[k]  = tmp + k*ht-<nv;
    }
    ht->eld = 1;
    memset(ht->map, 0, (unsigned long)ht->hsz * sizeof(hl_t));
    memset(ht->hd, 0, (unsigned long)ht->esz * sizeof(hd_t));

    /* reinsert known elements */
    for (i = 0; i < bs->ld; ++i) {
        b = bs->hd[i];
        for (j = 3; j < b[1]; ++j) {
            e = oev[b[j]];
            b[j]  = insert_in_hash_table_no_enlargement_check(e, ht);
        }
        for (; j < b[2]; j += 4) {
            e       = oev[b[j]];
            b[j]    = insert_in_hash_table_no_enlargement_check(e, ht);
            e       = oev[b[j+1]];
            b[j+1]  = insert_in_hash_table_no_enlargement_check(e, ht);
            e       = oev[b[j+2]];
            b[j+2]  = insert_in_hash_table_no_enlargement_check(e, ht);
            e       = oev[b[j+3]];
            b[j+3]  = insert_in_hash_table_no_enlargement_check(e, ht);
        }
    }

    const len_t pld = psl->ld;
    spair_t *ps = psl->p;
    for (i = 0; i < pld; ++i) {
        e = oev[ps[i].lcm];
        ps[i].lcm = insert_in_hash_table_no_enlargement_check(e, ht);
    }
    /* note: all memory is allocated as a big block, so it is
     *       enough to free oev[0].       */
    free(oev[0]);
    free(oev);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->rght_ctime  +=  ct1 - ct0;
    st->rght_rtime  +=  rt1 - rt0;
}

/* we can check equality of lcm and multiplication of two monomials
 * by their hash values. If both hash values are NOT the same, then
 * the corresponding exponent vectors CANNOT be the same. */
static inline int lcm_equals_multiplication(
    const hl_t a,
    const hl_t b,
    const ht_t * const ght,
    const hl_t lcm,
    const ht_t *const lht
    )
{
    const hd_t ha = ght->hd[a];
    const hd_t hb = ght->hd[b];
    const hd_t hl = lht->hd[lcm];

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
    const hl_t b,
    const ht_t * const ght,
    ht_t *lht
    )
{
    len_t i;

    /* exponents of basis elements, thus from global hash table */
    const exp_t * const ea = ght->ev + a;
    const exp_t * const eb = ght->ev + b;

    for (i = 0; i < ght->nv; ++i) {
        etmp[i]  = ea[i] < eb[i] ? eb[i] : ea[i];
    }
    /* goes into local hash table for spairs */
    return insert_in_hash_table(etmp, lht);
}

/* we try monomial division including check if divisibility is
 * fulfilled. */
static inline hl_t monomial_division_with_check(
    const hl_t a,
    const hl_t b,
    ht_t * ht
    )
{
    len_t i;

    const hd_t ha = ht->hd[a];
    const hd_t hb = ht->hd[b];
    /* short divisor mask check */
    if (hb.sdm & ~ha.sdm) {
        return 0;
    }

    const len_t nv  = ht->nv;
    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];
    etmp[0]  = ea[0] - eb[0];

    i = nv & 1 ? 1 : 0;
    for (; i < nv; i += 2) {
        if (ea[i] < eb[i] || ea[i+1] < eb[i+1]) {
            return 0;
        } else {
            etmp[i]    = ea[i] - eb[i];
            etmp[i+1]  = ea[i+1] - eb[i+1];
        }
    }
    return insert_in_hash_table(etmp, ht);
}

/* it is assumed that b divides a, thus no tests for
 * divisibility at all */
static inline hl_t monomial_division_no_check(
    const hl_t a,
    const hl_t b,
    ht_t *ht
    )
{
    len_t i;

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    const len_t nv  = ht->nv;
    i = nv & 1 ? 1 : 0;
    for (; i < nv; i += 2) {
        etmp[i]    = ea[i]   - eb[i];
        etmp[i+1]  = ea[i+1] - eb[i+1];
    }
    etmp[0]  = ea[0] - eb[0];
    return insert_in_hash_table(etmp, ht);
}

static inline dt_t *multiplied_polynomial_to_matrix_row(
    const val_t hm,
    const deg_t deg,
    const exp_t * const em,
    const dt_t *poly,
    ht_t *ht
    )
{
  len_t i;

  dt_t *row = (dt_t *)malloc((unsigned long)poly[2] * sizeof(dt_t));
  row[0]    = poly[0];
  row[1]    = poly[1];
  row[2]    = poly[2];
  /* hash table product insertions appear only here:
   * we check for hash table enlargements first and then do the insertions
   * without further elargment checks there */
  while (ht->eld+poly[2]-3 >= ht->esz) {
    enlarge_global_hash_table(ht);
  }
  /* printf("poly[1] %d | poly[2] %d\n", poly[1], poly[2]); */
  for (i = 3; i < poly[1]; ++i) {
    row[i]  = insert_in_hash_table_product_special(
                hm, deg, em, poly[i], ht);
  }
  for (;i < poly[2]; i += 4) {
    row[i]    = insert_in_hash_table_product_special(
                  hm, deg, em, poly[i], ht);
    row[i+1]  = insert_in_hash_table_product_special(
                  hm, deg, em, poly[i+1], ht);
    row[i+2]  = insert_in_hash_table_product_special(
                  hm, deg, em, poly[i+2], ht);
    row[i+3]  = insert_in_hash_table_product_special(
                  hm, deg, em, poly[i+3], ht);
  }

  return row;
}

static inline hl_t *reset_idx_in_global_hash_table_and_free_hcm(
        hl_t *hcm,
        ht_t *ht,
        const len_t nc
        )
{
    len_t i;

    hd_t *hd  = ht->hd;
    for (i = 0; i < nc; ++i) {
        hd[hcm[i]].idx  = 0;
    }
    free(hcm);

    return NULL;
}
