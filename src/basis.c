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
 * \file basis.c
 * \brief General and global data
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "basis.h"

void normalize_initial_basis_16(
        bs_t *bs,
        const len_t ne
        )
{
    len_t i, j;
    int64_t tmp1, tmp2, tmp3, tmp4;

    const int32_t fc  = gbfc;

    for (i = 0; i < ne; ++i) {
        cf16_t *p = (cf16_t *)bs->m[i].cl;

        const int32_t inv = mod_p_inverse_32((int32_t)p[0], fc);

        for (j = 0; j < bs->m[i].of; ++j) {
            tmp1  =   ((int64_t)p[j] * inv) % fc;
            tmp1  +=  (tmp1 >> 63) & fc;
            p[j]  =   (cf16_t)tmp1;
        }
        for (; j < bs->m[i].sz; j += 4) {
            tmp1    =   ((int64_t)p[j] * inv) % fc;
            tmp2    =   ((int64_t)p[j+1] * inv) % fc;
            tmp3    =   ((int64_t)p[j+2] * inv) % fc;
            tmp4    =   ((int64_t)p[j+3] * inv) % fc;
            tmp1    +=  (tmp1 >> 63) & fc;
            tmp2    +=  (tmp2 >> 63) & fc;
            tmp3    +=  (tmp3 >> 63) & fc;
            tmp4    +=  (tmp4 >> 63) & fc;
            p[j]    =   (cf16_t)tmp1;
            p[j+1]  =   (cf16_t)tmp2;
            p[j+2]  =   (cf16_t)tmp3;
            p[j+3]  =   (cf16_t)tmp4;
        }
    }
}

void normalize_initial_basis_32(
        bs_t *bs,
        const len_t ne
        )
{
    len_t i, j;
    int64_t tmp1, tmp2, tmp3, tmp4;

    const int32_t fc  = gbfc;

    for (i = 0; i < ne; ++i) {
        cf32_t *p = (cf32_t *)bs->m[i].cl;

        const int32_t inv = mod_p_inverse_32((int32_t)p[0], fc);

        for (j = 0; j < bs->m[i].of; ++j) {
            tmp1  =   ((int64_t)p[j] * inv) % fc;
            tmp1  +=  (tmp1 >> 63) & fc;
            p[j]  =   (cf32_t)tmp1;
        }
        for (; j < bs->m[i].sz; j += 4) {
            tmp1    =   ((int64_t)p[j] * inv) % fc;
            tmp2    =   ((int64_t)p[j+1] * inv) % fc;
            tmp3    =   ((int64_t)p[j+2] * inv) % fc;
            tmp4    =   ((int64_t)p[j+3] * inv) % fc;
            tmp1    +=  (tmp1 >> 63) & fc;
            tmp2    +=  (tmp2 >> 63) & fc;
            tmp3    +=  (tmp3 >> 63) & fc;
            tmp4    +=  (tmp4 >> 63) & fc;
            p[j]    =   (cf32_t)tmp1;
            p[j+1]  =   (cf32_t)tmp2;
            p[j+2]  =   (cf32_t)tmp3;
            p[j+3]  =   (cf32_t)tmp4;
        }
    }
}
