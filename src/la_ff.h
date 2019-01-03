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
 * \file la_ff.h
 * \brief Implementation of the linear algebra parts over finite fields.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_LA_FF_H
#define GB_LA_FF_H

#include "data.h"
#include "hash.h"
#include "time.h"
#include "tools.h"

inline mat_t *initialize_matrix(
        )
{
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    return mat;
}

inline void free_matrix(
        mat_t **matp
        )
{
    mat_t *mat  = *matp;
    free(mat->pv);
    free(mat->npv);
    free(mat->tbr);
    free(mat);
    *matp = mat;
}

static inline row_t normalize_sparse_matrix_row_16(
        row_t row
        )
{
    len_t i;

    int64_t tmp1, tmp2, tmp3, tmp4;

    cf16_t *cf        = (cf16_t *)row.cl;
    const int32_t fc  = gbfc;
    const int32_t inv = mod_p_inverse_32((int32_t)cf[0], (int32_t)fc);

    for (i = 0; i < row.of; ++i) {
        tmp1  =   ((int64_t)cf[i] * inv) % fc;
        tmp1  +=  (tmp1 >> 63) & fc;
        cf[i] =   (cf16_t)tmp1;
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = row.of; i < row.sz; i += 4) {
        tmp1    =   ((int64_t)cf[i] * inv) % fc;
        tmp2    =   ((int64_t)cf[i+1] * inv) % fc;
        tmp3    =   ((int64_t)cf[i+2] * inv) % fc;
        tmp4    =   ((int64_t)cf[i+3] * inv) % fc;
        tmp1    +=  (tmp1 >> 63) & fc;
        tmp2    +=  (tmp2 >> 63) & fc;
        tmp3    +=  (tmp3 >> 63) & fc;
        tmp4    +=  (tmp4 >> 63) & fc;
        cf[i]   =   (cf16_t)tmp1;
        cf[i+1] =   (cf16_t)tmp2;
        cf[i+2] =   (cf16_t)tmp3;
        cf[i+3] =   (cf16_t)tmp4;
    }

    return row;
}

static inline row_t normalize_sparse_matrix_row_32(
        row_t row
        )
{
    len_t i;

    int64_t tmp1, tmp2, tmp3, tmp4;

    cf32_t *cf        = (cf32_t *)row.cl;
    const int32_t fc  = gbfc;
    const int32_t inv = mod_p_inverse_32((int32_t)cf[0], (int32_t)fc);

    for (i = 0; i < row.of; ++i) {
        tmp1  =   ((int64_t)cf[i] * inv) % fc;
        tmp1  +=  (tmp1 >> 63) & fc;
        cf[i] =   (cf32_t)tmp1;
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = row.of; i < row.sz; i += 4) {
        tmp1    =   ((int64_t)cf[i] * inv) % fc;
        tmp2    =   ((int64_t)cf[i+1] * inv) % fc;
        tmp3    =   ((int64_t)cf[i+2] * inv) % fc;
        tmp4    =   ((int64_t)cf[i+3] * inv) % fc;
        tmp1    +=  (tmp1 >> 63) & fc;
        tmp2    +=  (tmp2 >> 63) & fc;
        tmp3    +=  (tmp3 >> 63) & fc;
        tmp4    +=  (tmp4 >> 63) & fc;
        cf[i]   =   (cf32_t)tmp1;
        cf[i+1] =   (cf32_t)tmp2;
        cf[i+2] =   (cf32_t)tmp3;
        cf[i+3] =   (cf32_t)tmp4;
    }

    return row;
}

row_t reduce_dense_row_by_sparse_pivots_16(
        int64_t *dr,
        mat_t *mat,
        const ci_t sc
        );

mat_t *exact_sparse_reduced_echelon_form_16(
        mat_t *mat,
        const stat_t *st
        );

mat_t *exact_sparse_linear_algebra_16(
        mat_t *mat,
        stat_t *st
        );
#endif
