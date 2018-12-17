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

#include "order.h"

int initial_basis_cmp_lex(
        const void *a,
        const void *b
        )
{
    len_t i;

    const hd_t * const ha = ((mon_t *)a)->h[0];
    const hd_t * const hb = ((mon_t *)b)->h[0];

    const exp_t * const ea  = ha->exp;
    const exp_t * const eb  = hb->exp;

    /* lexicographical */
    for (i = 1; i < gbnv; ++i) {
        if (ea[i] < eb[i]) {
            return -1;
        }
        if (ea[i] > eb[i]) {
            return 1;
        }
    }
    return 0;
}

int initial_basis_cmp_drl(
        const void *a,
        const void *b
        )
{
    len_t i;

    const hd_t * const ha = ((mon_t *)a)->h[0];
    const hd_t * const hb = ((mon_t *)b)->h[0];

    const deg_t da = ha->deg;
    const deg_t db = hb->deg;

    /* DRL */
    if (da < db) {
        return -1;
    } else {
        if (da != db) {
            return 1;
        }
    }

    const exp_t * const ea  = ha->exp;
    const exp_t * const eb  = hb->exp;

    /* note: reverse lexicographical */
    for (i = gbnv; i > 0; --i) {
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

int matrix_row_cmp(
        const void *a,
        const void *b
        )
{
    dt_t va, vb;
    /* compare pivot resp. column index */
    va  = ((dt_t **)a)[0][3];
    vb  = ((dt_t **)b)[0][3];
    if (va > vb) {
        return 1;
    }
    if (va < vb) {
        return -1;
    }
    /* same column index => compare density of row */
    va  = ((dt_t **)a)[0][2];
    vb  = ((dt_t **)b)[0][2];
    if (va > vb) {
        return 1;
    }
    if (va < vb) {
        return -1;
    }
    return 0;
}

int dense_matrix_row_cmp(
        const void *a,
        const void *b
        )
{
    const cf_t pa = ((cf_t **)a)[0][0];
    const cf_t pb = ((cf_t **)b)[0][0];

    return pa-pb;
}

/* comparison for monomials (in local hash table) */
int monomial_cmp_pivots_drl(
        const hl_t a,
        const hl_t b,
        const ht_t *const ht
        )
{
    len_t i;

    const hd_t ha = ht->hd[a];
    const hd_t hb = ht->hd[b];
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

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    /* note: reverse lexicographical */
    for (i = ht->nv; i > 0; --i) {
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

int monomial_cmp_pivots_lex(
        const hl_t a,
        const hl_t b,
        const ht_t *const ht
        )
{
    len_t i;

    const hd_t ha = ht->hd[a];
    const hd_t hb = ht->hd[b];
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

    const exp_t * const ea  = ht->ev[a];
    const exp_t * const eb  = ht->ev[b];

    /* lexicographical */
    for (i = 0; i < ht->nv; ++i) {
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

/* comparison for hash-column-maps */
int hcm_cmp_pivots_drl(
        const void *a,
        const void *b,
        void *ht
        )
{
    const hl_t ma  = ((hl_t *)a)[0];
    const hl_t mb  = ((hl_t *)b)[0];

    return monomial_cmp_pivots_drl(ma, mb, ht);
}

int hcm_cmp_pivots_lex(
        const void *a,
        const void *b,
        void *ht
        )
{
    const hl_t ma  = ((hl_t *)a)[0];
    const hl_t mb  = ((hl_t *)b)[0];

    return monomial_cmp_pivots_lex(ma, mb, ht);
}

/* comparison for s-pairs lcms in the update hash table */
int spair_cmp_deglex(
        const void *a,
        const void *b
        )
{
    const hd_t * const la = ((spair_t *)a)->lcm;
    const hd_t * const lb = ((spair_t *)b)->lcm;

    if (la->deg != lb->deg) {
        return (la->deg < lb->deg) ? -1 : 1;
    } else {
        return (int)monomial_cmp(la, lb);
    }
}

int spair_cmp_lht_drl(
        const void *a,
        const void *b
        )
{
    const hl_t la = ((spair_t *)a)->lcm;
    const hl_t lb = ((spair_t *)b)->lcm;

    return (int)monomial_cmp(la, lb, lht);
}

/* comparison for s-pairs once their lcms are in the global hash table */
int spair_cmp_ght_deglex(
        const void *a,
        const void *b
        )
{
    const hl_t la = ((spair_t *)a)->lcm;
    const hl_t lb = ((spair_t *)b)->lcm;

    if (hd[la].deg != hd[lb].deg) {
        return (hd[la].deg < hd[lb].deg) ? -1 : 1;
    } else {
        return (int)monomial_cmp(la, lb, ght);
    }
}

int spair_cmp_ght_drl(
        const void *a,
        const void *b,
        void *htl
        )
{
    const hl_t la = ((spair_t *)a)->lcm;
    const hl_t lb = ((spair_t *)b)->lcm;

    return (int)monomial_cmp(la, lb, ght);
}
