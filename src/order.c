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
static int initial_input_cmp_lex(
        const void *a,
        const void *b,
        void *htp
        )
{
    len_t i;
    ht_t *ht  = htp;

    const hm_t ha  = ((hm_t **)a)[0][3];
    const hm_t hb  = ((hm_t **)b)[0][3];

    const exp_t * const ea  = ht->ev[ha];
    const exp_t * const eb  = ht->ev[hb];

    /* lexicographical */
    const len_t nv  = ht->nv;

    i = 0;
    while(i < nv-1 && ea[i] == eb[i]) {
        ++i;
    }
    return ea[i] - eb[i];
}

static int initial_input_cmp_drl(
        const void *a,
        const void *b,
        void *htp
        )
{
    len_t i;
    ht_t *ht  = htp;

    const hm_t ha  = ((hm_t **)a)[0][3];
    const hm_t hb  = ((hm_t **)b)[0][3];

    const deg_t da = ht->hd[ha].deg;
    const deg_t db = ht->hd[hb].deg;

    /* DRL */
    if (da < db) {
        return -1;
    } else {
        if (da != db) {
            return 1;
        }
    }

    const exp_t * const ea  = ht->ev[ha];
    const exp_t * const eb  = ht->ev[hb];

    /* note: reverse lexicographical */
    i = ht->nv - 1;
    while (i > 0 && ea[i] == eb[i]) {
        --i;
    }
    return ea[i] - eb[i];
}

static int matrix_row_cmp(
        const void *a,
        const void *b
        )
{
    hm_t va, vb;
    /* compare pivot resp. column index */
    va  = ((hm_t **)a)[0][3];
    vb  = ((hm_t **)b)[0][3];
    if (va > vb) {
        return 1;
    }
    if (va < vb) {
        return -1;
    }
    /* same column index => compare density of row */
    va  = ((hm_t **)a)[0][2];
    vb  = ((hm_t **)b)[0][2];
    if (va > vb) {
        return 1;
    }
    if (va < vb) {
        return -1;
    }
    return 0;
}

static inline hm_t **sort_matrix_rows(
        hm_t **matdt)
{
    qsort(matdt, (unsigned long)nrows, sizeof(hm_t *), &matrix_row_cmp);
    return matdt;
}

static int dense_matrix_row_cmp(
        const void *a,
        const void *b
        )
{
    const cf32_t pa = ((cf32_t **)a)[0][0];
    const cf32_t pb = ((cf32_t **)b)[0][0];

    return pa-pb;
}

static inline cf32_t **sort_dense_matrix_rows(
        cf32_t **dm,
        const len_t nr
        )
{
    qsort(dm, (unsigned long)nr, sizeof(cf32_t *), &dense_matrix_row_cmp);
    return dm;
}

static int monomial_cmp_pivots_drl(
        const hl_t a,
        const hl_t b
        )
{
    len_t i;

    const hd_t ha = hds[a];
    const hd_t hb = hds[b];
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

    const exp_t * const ea  = evs[a];
    const exp_t * const eb  = evs[b];

    /* printf("ea ");
     * for (int ii=0; ii < nvars; ++ii) {
     *     printf("%d ", ea[ii]);
     * }
     * printf("\n");
     * printf("eb ");
     * for (int ii=0; ii < nvars; ++ii) {
     *     printf("%d ", eb[ii]);
     * }
     * printf("\n"); */

    /* note: reverse lexicographical */
    i = nvars - 1;
    while (i > 0 && ea[i] == eb[i]) {
        --i;
    }
    return ea[i] - eb[i];
}

static int monomial_cmp_pivots_lex(
        const hl_t a,
        const hl_t b
        )
{
    len_t i;

    const hd_t ha = hds[a];
    const hd_t hb = hds[b];
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

    const exp_t * const ea  = evs[a];
    const exp_t * const eb  = evs[b];

    /* lexicographical */
    const len_t nv  = nvars;

    i = 0;
    while(i < nv-1 && ea[i] == eb[i]) {
        ++i;
    }
    return eb[i] - ea[i];
}

static inline int monomial_cmp_drl(
        const hl_t a,
        const hl_t b,
        const ht_t *ht
        )
{
    len_t i;

    const deg_t da = ht->hd[a].deg;
    const deg_t db = ht->hd[b].deg;

    /* DRL */
    if (da > db) {
        return 1;
    } else {
        if (da != db) {
            return -1;
        }
    }

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    i = ht->nv - 1;
    while (i > 0 && ea[i] == eb[i]) {
        --i;
    }
    return eb[i] - ea[i];
}

static inline int monomial_cmp_lex(
        const hl_t a,
        const hl_t b,
        const ht_t *ht
        )
{
    len_t i;

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];
    const len_t nv  = ht->nv;

    i = 0;
    while(i < nv-1 && ea[i] == eb[i]) {
        ++i;
    }
    return ea[i] - eb[i];
}

/* comparison for spair generators */
static int gens_cmp(
        const void *a,
        const void *b
        )
{
    const len_t ga = *((len_t *)a);
    const len_t gb = *((len_t *)b);

    return (ga - gb);
}

/* comparison for sparse rows in preparation for generation of pbm files */
static int pbm_cmp(
        const void *a,
        const void *b
        )
{
    const hm_t ca = *((hm_t *)a);
    const hm_t cb = *((hm_t *)b);

    return (ca - cb);
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
        const void *b,
        void *htp
        )
{
    const hl_t la   = ((spair_t *)a)->lcm;
    const hl_t lb   = ((spair_t *)b)->lcm;
    const ht_t *ht  = (ht_t *)htp;

    if (ht->hd[la].deg != ht->hd[lb].deg) {
        return (ht->hd[la].deg < ht->hd[lb].deg) ? -1 : 1;
    } else {
        return (int)monomial_cmp(la, lb, ht);
    }
}

static int spair_cmp_drl(
        const void *a,
        const void *b,
        void *htp
        )
{
    const hl_t la   = ((spair_t *)a)->lcm;
    const hl_t lb   = ((spair_t *)b)->lcm;
    const ht_t *ht  = (ht_t *)htp;

    return (int)monomial_cmp(la, lb, ht);
}

static int spair_degree_cmp(
        const void *a,
        const void *b
        )
{
    const deg_t da  = hd[((spair_t *)a)->lcm].deg;
    const deg_t db  = hd[((spair_t *)b)->lcm].deg;

    return (da-db);
}
