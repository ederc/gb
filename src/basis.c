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

/* finite field stuff */
static void initialize_basis_ff(
        int32_t ngens
        )
{
    bload = 0;
    bsize = 2*ngens;

    gbcf_ff = (cf32_t **)malloc((unsigned long)bsize * sizeof(cf32_t *));
    gbdt    = (dt_t **)malloc((unsigned long)bsize * sizeof(dt_t *));
    lms     = (sdm_t *)malloc((unsigned long)bsize * sizeof(sdm_t));
    red     = (int8_t *)malloc((unsigned long)bsize * sizeof(int8_t));
    memset(red, 0, (unsigned long)bsize * sizeof(int8_t));
}

static inline void check_enlarge_basis_ff(
        len_t added
        )
{
    if (bload+added >= bsize) {
        bsize = bsize*2 > bload+added ? bsize*2 : bload+added;
        gbcf_ff = realloc(gbcf_ff, (unsigned long)bsize * sizeof(cf32_t *));
        gbdt    = realloc(gbdt, (unsigned long)bsize * sizeof(dt_t *));
        lms     = realloc(lms, (unsigned long)bsize * sizeof(sdm_t));
        red     = realloc(red, (unsigned long)bsize * sizeof(int8_t));
        memset(red+bload, 0, (unsigned long)(bsize-bload) * sizeof(int8_t));
    }
}

static void free_basis_ff(
        void
        )
{
    len_t i;
    if (gbcf_ff) {
        for (i = 0; i < bload; ++i) {
            free(gbcf_ff[i]);
            free(gbdt[i]);
        }
        free(gbcf_ff);
        gbcf_ff = NULL;
        free(gbdt);
        gbdt  = NULL;
        free(lms);
        lms = NULL;
        free(red);
        red = NULL;
        blold = 0;
        bload = 0;
        bsize = 0;
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
static void initialize_basis_q(
        int32_t ngens
        )
{
    bload = 0;
    bsize = 2*ngens;

    gbcf_q  = (mpz_t **)malloc((unsigned long)bsize * sizeof(mpz_t *));
    gbdt    = (dt_t **)malloc((unsigned long)bsize * sizeof(dt_t *));
    lms     = (sdm_t *)malloc((unsigned long)bsize * sizeof(sdm_t));
    red     = (int8_t *)malloc((unsigned long)bsize * sizeof(int8_t));
    memset(red, 0, (unsigned long)bsize * sizeof(int8_t));
}

static inline void check_enlarge_basis_q(
        len_t added
        )
{
    if (bload+added >= bsize) {
        bsize = bsize*2 > bload+added ? bsize*2 : bload+added;
        gbcf_q  = realloc(gbcf_q, (unsigned long)bsize * sizeof(mpz_t *));
        gbdt    = realloc(gbdt, (unsigned long)bsize * sizeof(dt_t *));
        lms     = realloc(lms, (unsigned long)bsize * sizeof(sdm_t));
        red     = realloc(red, (unsigned long)bsize * sizeof(int8_t));
        memset(red+bload, 0, (unsigned long)(bsize-bload) * sizeof(int8_t));
    }
}

static void free_basis_q(
        void
        )
{
    len_t i;
    if (gbcf_q) {
        for (i = 0; i < bload; ++i) {
            free(gbcf_q[i]);
            free(gbdt[i]);
        }
        free(gbcf_q);
        gbcf_q  = NULL;
        free(gbdt);
        gbdt  = NULL;
        free(lms);
        lms = NULL;
        free(red);
        red = NULL;
        blold = 0;
        bload = 0;
        bsize = 0;
    }
}

static inline void normalize_initial_basis_q(
        void
        )
{
}
