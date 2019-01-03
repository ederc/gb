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
 * \file order.h
 * \brief Order procedures
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_ORDER_H
#define GB_ORDER_H

#include "data.h"

int initial_basis_cmp_lex(
        const void *a,
        const void *b
        );

int initial_basis_cmp_drl(
        const void *a,
        const void *b
        );

int matrix_row_cmp(
        const void *a,
        const void *b
        );

int dense_matrix_row_cmp(
        const void *a,
        const void *b
        );

int monomial_cmp_pivots_drl(
        const hd_t *a,
        const hd_t *b
        );

int monomial_cmp_pivots_lex(
        const hd_t *a,
        const hd_t *b
        );

int hcm_cmp_pivots_drl(
        const void *a,
        const void *b
        );

int hcm_cmp_pivots_lex(
        const void *a,
        const void *b
        );

int spair_cmp_deglex(
        const void *a,
        const void *b
        );

int spair_cmp_drl(
        const void *a,
        const void *b
        );

inline mat_t *sort_matrix_rows(
        mat_t *mat
        )
{
    qsort(mat->pv, (unsigned long)mat->ncl, sizeof(row_t *), &matrix_row_cmp);
    return mat;
}

static inline int monomial_cmp_drl(
        const hd_t *a,
        const hd_t *b
        )
{
    len_t i;

    const deg_t da = a->deg;
    const deg_t db = b->deg;

    /* DRL */
    if (da > db) {
        return 1;
    } else {
        if (da != db) {
            return -1;
        }
    }

    const exp_t * const ea  = a->exp;
    const exp_t * const eb  = b->exp;

    for (i = gbnv; i > 0; --i) {
        if (ea[i-1] < eb[i-1]) {
            return 1;
        } else {
            if (ea[i-1] != eb[i-1]) {
                return -1;
            }
        }
    }
    return 0;
}

static inline int monomial_cmp_lex(
        const hd_t *a,
        const hd_t *b
        )
{
    len_t i;

    const exp_t * const ea  = a->exp;
    const exp_t * const eb  = b->exp;

    const len_t nv  = gbnv;

    for (i = 0; i < nv; ++i) {
        if (ea[i] < eb[i]) {
            return -1;
        }
        if (ea[i] > eb[i]) {
            return 1;
        }
    }
    return 0;
}
#endif
