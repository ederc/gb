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
 * \file symbol.h
 * \brief Symbolic preprocessing routines
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_SYMBOL_H
#define GB_SYMBOL_H

#include "data.h"
#include "hash.h"
#include "time.h"

mat_t *symbolic_preprocessing(
        mat_t *mat,
        ht_t *ht,
        const bs_t *const bs,
        stat_t *st
        );

mat_t *select_spairs_by_minimal_degree(
        mat_t *mat,
        ps_t *psl,
        ht_t *ht,
        stat_t *st
        );

inline dt_t *find_multiplied_reducer(
        const dt_t m,
        ht_t *ht,
        const bs_t *const bs
        )
{
    len_t i, k;
    deg_t d = 0;
    dt_t *b;
    const exp_t * const e  = ht->ev[m];
    exp_t *f;

    const sdm_t *const lms  = bs->lm;

    const len_t bl  = bs->ld;
    const len_t os  = ht->nv & 1 ? 1 : 0;
    i = ht->hd[m].div;

    const sdm_t ns  = ~ht->hd[m].sdm;

    exp_t etmp[ht->nv];

start:
    while (i < bl-3) {
        if (lms[i] & ns &&
                lms[i+1] & ns &&
                lms[i+2] & ns &&
                lms[i+3] & ns) {
            i +=  4;
            continue;
        }
        while (lms[i] & ns) {
            i++;
        }
        b = bs->hd[i];
        f = ht->ev[b[3]];
        if ((e[0]-f[0]) < 0) {
            i++;
            goto start;
        }
        for (k = os; k < ht->nv; k += 2) {
            if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
                i++;
                goto start;
            }
        }
        for (k = 0; k < ht->nv; ++k) {
            etmp[k] = e[k] - f[k];
        }
        const hl_t h  = ht->hd[m].val - ht->hd[b[3]].val;
        for (k = 0; k < ht->nv; ++k) {
            d += etmp[k];
        }
        b = multiplied_polynomial_to_matrix_row(h, d, etmp, b, ht);
        ht->hd[m].div = i;
        return b;
    }
start2:
    while (i < bl) {
        if (lms[i] & ns) {
            i++;
            continue;
        }
        b = bs->hd[i];
        f = ht->ev[b[3]];
        if ((e[0]-f[0]) < 0) {
            i++;
            goto start2;
        }
        for (k = os; k < ht->nv; k += 2) {
            if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
                i++;
                goto start2;
            }
        }
        for (k = 0; k < ht->nv; ++k) {
            etmp[k] = e[k] - f[k];
        }
        const hl_t h  = ht->hd[m].val - ht->hd[b[3]].val;
        for (k = 0; k < ht->nv; ++k) {
            d += etmp[k];
        }
        b = multiplied_polynomial_to_matrix_row(h, d, etmp, b, ht);
        ht->hd[m].div = i;
        return b;
    }
    ht->hd[m].div = i;
    return NULL;
}
#endif
