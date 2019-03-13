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

#include "data.h"

static void free_basis(
        bs_t **bsp
        )
{
    len_t i;
    bs_t *bs  = *bsp;
    if (bs->cf_ff) {
        for (i = 0; i < bs->ld; ++i) {
            free(bs->cf_ff[i]);
            free(bs->hm[i]);
        }
        free(bs->cf_ff);
        bs->cf_ff = NULL;
        free(bs->hm);
        bs->hm  = NULL;
    }
    if (bs->cf_q) {
        for (i = 0; i < bs->ld; ++i) {
            free(bs->cf_q[i]);
            free(bs->hm[i]);
        }
        free(bs->cf_q);
        bs->cf_q  = NULL;
        free(bs->hm);
        bs->hm  = NULL;
    }
    free(bs->lm);
    bs->lm  = NULL;
    free(bs->red);
    bs->red = NULL;
    free(bs);
    bs  = NULL;
    *bsp  = bs;
}

/* finite field stuff */
static bs_t *initialize_basis_ff(
        const int32_t ngens
        )
{
    bs_t *bs  = (bs_t *)malloc(sizeof(bs_t));
    bs->ld  = 0;
    bs->sz  = 2*ngens;

    bs->cf_ff = (cf32_t **)malloc((unsigned long)bs->sz * sizeof(cf32_t *));
    bs->hm    = (dt_t **)malloc((unsigned long)bs->sz * sizeof(dt_t *));
    bs->lm    = (sdm_t *)malloc((unsigned long)bs->sz * sizeof(sdm_t));
    bs->red   = (int8_t *)calloc((unsigned long)bs->sz, sizeof(int8_t));

    return bs;
}

static inline void check_enlarge_basis_ff(
        bs_t *bs,
        const len_t added
        )
{
    if (bs->ld + added >= bs->sz) {
        bs->sz    = bs->sz * 2 > bs->ld + added ? bs->sz * 2 : bs->ld + added;
        bs->cf_ff = realloc(bs->cf_ff,
                (unsigned long)bs->sz * sizeof(cf32_t *));
        bs->hm    = realloc(bs->hm, (unsigned long)bs->sz * sizeof(dt_t *));
        bs->lm    = realloc(bs->lm, (unsigned long)bs->sz * sizeof(sdm_t));
        bs->red   = realloc(bs->red, (unsigned long)bs->sz * sizeof(int8_t));
        memset(bs->red+bs->ld, 0,
                (unsigned long)(bs->sz-bs->ld) * sizeof(int8_t));
    }
}

static inline void normalize_initial_basis_ff(
        void
        )
{
    len_t i, j;
    int64_t tmp1, tmp2, tmp3, tmp4;
    dt_t **dt = gbdt;

    for (i = 0; i < nrows; ++i) {
        cf32_t *row = gbcf_ff[dt[i][0]];

        const int32_t inv = mod_p_inverse_32((int32_t)row[0], (int32_t)fc);
        const len_t os    = dt[i][1]; 
        const len_t len   = dt[i][2]; 

        for (j = 0; j < os; ++j) {
            tmp1    =   ((int64_t)row[j] * inv) % fc;
            tmp1    +=  (tmp1 >> 63) & fc;
            row[j]  =   (cf32_t)tmp1;
        }
        for (j = os; j < len; j += 4) {
            tmp1      =   ((int64_t)row[j] * inv) % fc;
            tmp2      =   ((int64_t)row[j+1] * inv) % fc;
            tmp3      =   ((int64_t)row[j+2] * inv) % fc;
            tmp4      =   ((int64_t)row[j+3] * inv) % fc;
            tmp1      +=  (tmp1 >> 63) & fc;
            tmp2      +=  (tmp2 >> 63) & fc;
            tmp3      +=  (tmp3 >> 63) & fc;
            tmp4      +=  (tmp4 >> 63) & fc;
            row[j]    =   (hl_t)tmp1;
            row[j+1]  =   (hl_t)tmp2;
            row[j+2]  =   (hl_t)tmp3;
            row[j+3]  =   (hl_t)tmp4;
        }
    }
}

/* characteristic zero stuff */
static bs_t *initialize_basis_q(
        const int32_t ngens
        )
{
    bs_t *bs  = (bs_t *)malloc(sizeof(bs_t));
    bs->ld  = 0;
    bs->sz  = 2*ngens;

    bs->cf_q  = (mpz_t **)malloc((unsigned long)bs->sz * sizeof(mpz_t *));
    bs->hm    = (dt_t **)malloc((unsigned long)bs->sz * sizeof(dt_t *));
    bs->lm    = (sdm_t *)malloc((unsigned long)bs->sz * sizeof(sdm_t));
    bs->red   = (int8_t *)calloc((unsigned long)bs->sz, sizeof(int8_t));

    return bs;
}

static inline void check_enlarge_basis_q(
        bs_t *bs,
        const len_t added
        )
{
    if (bs->ld + added >= bs->sz) {
        bs->sz    = bs->sz * 2 > bs->ld + added ? bs->sz * 2 : bs->ld + added;
        bs->cf_q  = realloc(bs->cf_q,
                (unsigned long)bs->sz * sizeof(mpz_t *));
        bs->hm    = realloc(bs->hm, (unsigned long)bs->sz * sizeof(dt_t *));
        bs->lm    = realloc(bs->lm, (unsigned long)bs->sz * sizeof(sdm_t));
        bs->red   = realloc(bs->red, (unsigned long)bs->sz * sizeof(int8_t));
        memset(bs->red+bs->ld, 0,
                (unsigned long)(bs->sz-bs->ld) * sizeof(int8_t));
    }
}

static inline void normalize_initial_basis_q(
        void
        )
{
}
