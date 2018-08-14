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

#include "hash.h"

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
void enlarge_hash_table(
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

void regenerate_hash_table(
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
