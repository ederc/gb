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
 * \file la.c
 * \brief Implementation of linear algebra.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "data.h"

static inline cf_t **free_dense_matrix(
        cf_t **dm
        )
{
    len_t i;

    for (i = 0; i < ncr; ++i) {
        free(dm[i]);
    }
    free(dm);
    dm  = NULL;

    return dm;
}

static inline void normalize_matrix_rows(
        cf_t **mat
        )
{
    len_t i, j;
    int64_t tmp1, tmp2, tmp3, tmp4;

    for (i = 0; i < nrows; ++i) {
        cf_t *row = mat[i];

        const int32_t inv = mod_p_inverse_32((int32_t)row[3], (int32_t)fc);

        for (j = 3; j < row[1]; ++j) {
            tmp1    =   ((int64_t)row[j] * inv) % fc;
            tmp1    +=  (tmp1 >> 63) & fc;
            row[j]  =   (cf_t)tmp1;
        }
        for (j = row[1]; j < row[2]; j += 4) {
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

static inline cf_t *normalize_dense_matrix_row(
        cf_t *row,
        const dt_t pc
        )
{
    len_t i;
    /* printf("normalize: %d != 1\n", row[0]); */

    const dt_t len  = ncr - pc;
    const dt_t os   = len % 4;
    /* printf("ncr %d | pc %d | os %d\n", ncr, pc, os); */
    int64_t tmp1, tmp2, tmp3, tmp4;

    const int32_t inv = mod_p_inverse_32((int32_t)row[0], (int32_t)fc);

    for (i = 1; i < os; ++i) {
        tmp1    =   ((int64_t)row[i] * inv) % fc;
        tmp1    +=  (tmp1 >> 63) & fc;
        row[i]  =   (cf_t)tmp1;
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = os; i < len; i += 4) {
        tmp1      =   ((int64_t)row[i] * inv) % fc;
        tmp2      =   ((int64_t)row[i+1] * inv) % fc;
        tmp3      =   ((int64_t)row[i+2] * inv) % fc;
        tmp4      =   ((int64_t)row[i+3] * inv) % fc;
        tmp1      +=  (tmp1 >> 63) & fc;
        tmp2      +=  (tmp2 >> 63) & fc;
        tmp3      +=  (tmp3 >> 63) & fc;
        tmp4      +=  (tmp4 >> 63) & fc;
        row[i]    =   (cf_t)tmp1;
        row[i+1]  =   (cf_t)tmp2;
        row[i+2]  =   (cf_t)tmp3;
        row[i+3]  =   (cf_t)tmp4;
    }
    row[0]  = 1;

    return row;
}

static inline cf_t *normalize_sparse_matrix_row(
        cf_t *row
        )
{
    len_t i;

    int64_t tmp1, tmp2, tmp3, tmp4;

    const int32_t inv = mod_p_inverse_32((int32_t)row[3], (int32_t)fc);

    for (i = 3; i < row[1]; ++i) {
        tmp1    =   ((int64_t)row[i] * inv) % fc;
        tmp1    +=  (tmp1 >> 63) & fc;
        row[i]  =   (cf_t)tmp1;
    }
    /* we need to set i to os since os < 1 is possible */
    for (i = row[1]; i < row[2]; i += 4) {
        tmp1      =   ((int64_t)row[i] * inv) % fc;
        tmp2      =   ((int64_t)row[i+1] * inv) % fc;
        tmp3      =   ((int64_t)row[i+2] * inv) % fc;
        tmp4      =   ((int64_t)row[i+3] * inv) % fc;
        tmp1      +=  (tmp1 >> 63) & fc;
        tmp2      +=  (tmp2 >> 63) & fc;
        tmp3      +=  (tmp3 >> 63) & fc;
        tmp4      +=  (tmp4 >> 63) & fc;
        row[i]    =   (cf_t)tmp1;
        row[i+1]  =   (cf_t)tmp2;
        row[i+2]  =   (cf_t)tmp3;
        row[i+3]  =   (cf_t)tmp4;
    }
    row[3]  = 1;

    return row;
}

static dt_t *reduce_dense_row_by_known_pivots_sparse_17_bit(
        int64_t *dr,
        dt_t *const *pivs,
        const hl_t dpiv,    /* pivot of dense row at the beginning */
        const dt_t tmp_pos  /* position of new coeffs array in tmpcf */
        )
{
    hl_t i, j;
    dt_t *dts;
    cf_t *cfs;
    const int64_t mod = (int64_t)fc;

    for (i = dpiv; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        /* printf("i %d\n", i); */
        if (pivs[i] == NULL) {
            continue;
        }

        /* printf("17-bit dr before reduction with piv %d\n", i);
         * for (int32_t l = 0; l < ncols; ++l) {
         *   printf("%ld ", dr[l]);
         * }
         * printf("\n"); */
        /* found reducer row, get multiplier */
        const int64_t mul = mod - dr[i];
        dts   = pivs[i];
        if (i < ncl) {
            cfs   = gbcf[dts[0]];
        } else {
            cfs   = tmpcf[dts[0]];
        }
        /* if (i >= ncl) {
         *     printf("cfs[%d]: ", dts[0]);
         *     for (j = 3; j < dts[2]; ++j) {
         *         printf("%d ", cfs[j]);
         *     }
         *     printf("\n");
         * } */
         /* for (j = 0; j < ncols; ++j) {
         *     printf("%ld ", dr[j]);
         * }
         * printf("\n"); */
        for (j = 3; j < dts[1]; ++j) {
            dr[dts[j]]  +=  mul * cfs[j];
        }
        for (; j < dts[2]; j += 4) {
            dr[dts[j]]    +=  mul * cfs[j];
            dr[dts[j+1]]  +=  mul * cfs[j+1];
            dr[dts[j+2]]  +=  mul * cfs[j+2];
            dr[dts[j+3]]  +=  mul * cfs[j+3];
        }
        dr[i] = 0;
        /* for (j = 0; j < ncols; ++j) {
         *     printf("%ld ", dr[j]);
         * }
         * printf("\n----------\n"); */
    }
    /* printf("reduction step done\n"); */

    dt_t *dt  = (dt_t *)malloc((unsigned long)(ncr+3) * sizeof(dt_t));
    cf_t *cf  = (cf_t *)malloc((unsigned long)(ncr+3) * sizeof(cf_t));
    j = 3;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                dt[j] = (dt_t)i;
                cf[j] = (cf_t)dr[i];
                j++;
            }
        }
    }
    if (j == 3) {
        free(dt);
        dt = NULL;
    } else {
        dt  = realloc(dt, (unsigned long)j * sizeof(dt_t));
        cf  = realloc(cf, (unsigned long)j * sizeof(cf_t));
        dt[0] = tmp_pos;
        cf[0] = 0;
        dt[1] = (j-3) % 4 + 3;
        cf[1] = dt[1];
        dt[2] = cf[2] = j;
        tmpcf[tmp_pos]  = cf;
    }
    return dt;
}

static dt_t *reduce_dense_row_by_known_pivots_sparse_31_bit(
        int64_t *dr,
        dt_t *const *pivs,
        const hl_t dpiv,    /* pivot of dense row at the beginning */
        const dt_t tmp_pos  /* position of new coeffs array in tmpcf */
        )
{
    hl_t i, j;
    cf_t *cfs;
    dt_t *dts;
    const int64_t mod   = (int64_t)fc;
    const int64_t mod2  = (int64_t)fc * fc;

    for (i = dpiv; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = (int64_t)dr[i];
        dts   = pivs[i];
        if (i < ncl) {
            cfs   = gbcf[dts[0]];
        } else {
            cfs   = tmpcf[dts[0]];
        }
        for (j = 3; j < dts[1]; ++j) {
            dr[dts[j]]  -=  mul * cfs[j];
            dr[dts[j]]  +=  (dr[dts[j]] >> 63) & mod2;
        }
        for (; j < dts[2]; j += 4) {
            dr[dts[j]]    -=  mul * cfs[j];
            dr[dts[j+1]]  -=  mul * cfs[j+1];
            dr[dts[j+2]]  -=  mul * cfs[j+2];
            dr[dts[j+3]]  -=  mul * cfs[j+3];
            dr[dts[j]]    +=  (dr[dts[j]] >> 63) & mod2;
            dr[dts[j+1]]  +=  (dr[dts[j+1]] >> 63) & mod2;
            dr[dts[j+2]]  +=  (dr[dts[j+2]] >> 63) & mod2;
            dr[dts[j+3]]  +=  (dr[dts[j+3]] >> 63) & mod2;
        }
        dr[i] = 0;
    }

    dt_t *dt  = (dt_t *)malloc((unsigned long)(ncr+3) * sizeof(dt_t));
    cf_t *cf  = (cf_t *)malloc((unsigned long)(ncr+3) * sizeof(cf_t));
    j = 3;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                dt[j] = (dt_t)i;
                cf[j] = (cf_t)dr[i];
                j++;
            }
        }
    }
    if (j == 3) {
        free(dt);
        dt = NULL;
    } else {
        dt  = realloc(dt, (unsigned long)j * sizeof(dt_t));
        cf  = realloc(cf, (unsigned long)j * sizeof(cf_t));
        dt[0] = tmp_pos;
        cf[0] = 0;
        dt[1] = (j-3) % 4 + 3;
        cf[1] = dt[1];
        dt[2] = cf[2] = j;
        tmpcf[tmp_pos]  = cf;
    }
    return dt;
}

static len_t reduce_dense_row_by_all_pivots_17_bit(
        int64_t *dr,
        cf_t **tbr,
        dt_t *const *pivs,
        cf_t *const *dpivs,
        const len_t pc
        )
{
    hl_t i, j, k, l;
    const int64_t mod = (int64_t)fc;
    len_t np  = -1;
    cf_t *red;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();
    /* step 1: reduce by sparse known pivots */
    for (i = pc; i < ncl; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        /* printf("i %d\n", i); */
        if (pivs[i] == NULL) {
            continue;
        }

        /* printf("17-bit dr before reduction with piv %d\n", i);
         * for (int32_t l = 0; l < ncols; ++l) {
         *   printf("%ld ", dr[l]);
         * }
         * printf("\n"); */
        /* found reducer row, get multiplier */
        const int64_t mul = mod - dr[i];
        const dt_t *dts   = pivs[i];
        const cf_t *cfs   = gbcf[dts[0]];
        /* printf("cfs[%d]: ", dts[0]);
         * for (j = 3; j < dts[2]; ++j) {
         *     printf("%d ", cfs[j]);
         * }
         * printf("\n"); */
        for (j = 3; j < dts[1]; ++j) {
            dr[dts[j]]  +=  mul * cfs[j];
        }
        for (; j < dts[2]; j += 4) {
            dr[dts[j]]    +=  mul * cfs[j];
            dr[dts[j+1]]  +=  mul * cfs[j+1];
            dr[dts[j+2]]  +=  mul * cfs[j+2];
            dr[dts[j+3]]  +=  mul * cfs[j+3];
        }
        dr[i] = 0;
        /* for (j = 0; j < ncols; ++j) {
         *     printf("%ld ", dr[j]);
         * }
         * printf("\n---------------------\n"); */
    }
    k = 0;
    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    step1_ctime  +=  ct1 - ct0;
    step1_rtime  +=  rt1 - rt0;
    /* step 2: reduce by new dense pivots */
    for (i = ncl; i < ncols; ++i) {
        /* printf("dr[%d] = %ld\n", i, dr[i]); */
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (dpivs[i-ncl] == NULL) {
            if (np == -1) {
                np  = i;
            }
            k++;
            continue;
        }

        /* printf("reduce dr[%d] = %ld by dpivs[%d] where ncl = %d\n", i, dr[i], i-ncl, ncl); */
        red = dpivs[i-ncl];
        const int64_t mul = mod - dr[i];
        const len_t os    = (ncols - i) % 4;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] +=  mul * red[l];
        }
        for (; j < ncols; l += 4, j += 4) {
            dr[j]   +=  mul * red[l];
            dr[j+1] +=  mul * red[l+1];
            dr[j+2] +=  mul * red[l+2];
            dr[j+3] +=  mul * red[l+3];
        }
        /* printf("red[%d]: ", i-ncl);
         * for (j = 0; j < ncols-i ; ++j) {
         *     printf("%d ", red[j]);
         * }
         * printf("\n");
         * for (j = 0; j < ncols; ++j) {
         *     printf("%ld ", dr[j]);
         * }
         * printf("\n"); */
    }
    if (k == 0) {
        return -1;
    }

    cf_t *row = (cf_t *)calloc((unsigned long)(ncols-np), sizeof(cf_t));
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row(row, np-ncl);
    }
    *tbr  = row;

    return np-ncl;
}

static len_t reduce_dense_row_by_all_pivots_31_bit(
        int64_t *dr,
        cf_t **tbr,
        dt_t *const *pivs,
        cf_t *const *dpivs,
        const len_t pc
        )
{
    hl_t i, j, k, l;
    const int64_t mod   = (int64_t)fc;
    const int64_t mod2  = (int64_t)fc * fc;
    len_t np  = -1;
    cf_t *red;

    /* step 1: reduce by sparse known pivots */
    for (i = pc; i < ncl; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        /* printf("i %d\n", i); */
        if (pivs[i] == NULL) {
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = (int64_t)dr[i];
        const dt_t *dts   = pivs[i];
        const cf_t *cfs   = gbcf[dts[0]];
        for (j = 3; j < dts[1]; ++j) {
            dr[dts[j]]  -=  mul * cfs[j];
            dr[dts[j]]  +=  (dr[dts[j]] >> 63) & mod2;
        }
        for (; j < dts[2]; j += 4) {
            dr[dts[j]]    -=  mul * cfs[j];
            dr[dts[j+1]]  -=  mul * cfs[j+1];
            dr[dts[j+2]]  -=  mul * cfs[j+2];
            dr[dts[j+3]]  -=  mul * cfs[j+3];
            dr[dts[j]]    +=  (dr[dts[j]] >> 63) & mod2;
            dr[dts[j+1]]  +=  (dr[dts[j+1]] >> 63) & mod2;
            dr[dts[j+2]]  +=  (dr[dts[j+2]] >> 63) & mod2;
            dr[dts[j+3]]  +=  (dr[dts[j+3]] >> 63) & mod2;
        }
        dr[i] = 0;
    }
    k = 0;
    /* step 2: reduce by new dense pivots */
    for (i = ncl; i < ncols; ++i) {
        /* printf("dr[%d] = %ld\n", i, dr[i]); */
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (dpivs[i-ncl] == NULL) {
            if (np == -1) {
                np  = i;
            }
            k++;
            continue;
        }

        red = dpivs[i-ncl];
        const int64_t mul = (int64_t) dr[i];
        const len_t os    = (ncols - i) % 4;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] -=  mul * red[l];
            dr[j] +=  (dr[j] >> 63) & mod2;
        }
        for (; j < ncols; l+=4, j += 4) {
            dr[j]   -=  mul * red[l];
            dr[j+1] -=  mul * red[l+1];
            dr[j+2] -=  mul * red[l+2];
            dr[j+3] -=  mul * red[l+3];
            dr[j]   +=  (dr[j] >> 63) & mod2;
            dr[j+1] +=  (dr[j+1] >> 63) & mod2;
            dr[j+2] +=  (dr[j+2] >> 63) & mod2;
            dr[j+3] +=  (dr[j+3] >> 63) & mod2;
        }
    }
    if (k == 0) {
        return -1;
    }

    cf_t *row = (cf_t *)calloc((unsigned long)(ncols-np), sizeof(cf_t));
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row(row, np-ncl);
    }
    *tbr  = row;

    return np-ncl;
}

static cf_t *reduce_dense_row_by_known_pivots_17_bit(
        int64_t *dr,
        dt_t *const *pivs,
        const hl_t dpiv    /* pivot of dense row at the beginning */
        )
{
    hl_t i, j;
    const int64_t mod = (int64_t)fc;

    for (i = dpiv; i < ncl; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        /* printf("i %d\n", i); */
        if (pivs[i] == NULL) {
            continue;
        }

        /* printf("17-bit dr before reduction with piv %d\n", i);
         * for (int32_t l = 0; l < ncols; ++l) {
         *   printf("%ld ", dr[l]);
         * }
         * printf("\n"); */
        /* found reducer row, get multiplier */
        const int64_t mul = mod - dr[i];
        const dt_t *dts   = pivs[i];
        const cf_t *cfs   = gbcf[dts[0]];
        /* printf("cfs[%d]: ", dts[0]);
         * for (j = 3; j < dts[2]; ++j) {
         *     printf("%d ", cfs[j]);
         * }
         * printf("\n");
         * for (j = 0; j < ncols; ++j) {
         *     printf("%ld ", dr[j]);
         * }
         * printf("\n"); */
        for (j = 3; j < dts[1]; ++j) {
            dr[dts[j]]  +=  mul * cfs[j];
        }
        for (; j < dts[2]; j += 4) {
            dr[dts[j]]    +=  mul * cfs[j];
            dr[dts[j+1]]  +=  mul * cfs[j+1];
            dr[dts[j+2]]  +=  mul * cfs[j+2];
            dr[dts[j+3]]  +=  mul * cfs[j+3];
        }
        dr[i] = 0;
        /* for (j = 0; j < ncols; ++j) {
         *     printf("%ld ", dr[j]);
         * }
         * printf("\n----------\n"); */
    }
    /* printf("reduction step done\n"); */

    /* store a dense row for further dense gaussian elimination */
    cf_t *row  = (cf_t *)calloc(
            (unsigned long)(ncr), sizeof(cf_t));
    j = 0;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                j++;
                row[i-ncl]  = (cf_t)dr[i];
            }
        }
    }
    if (j == 0) {
        free(row);
        row = NULL;
    } else {
        /* printf("returned row: ");
         * for (i= 0; i < ncr; ++i) {
         *     printf("%d ", row[i]);
         * }
         * printf("\n"); */
    }
    return row;
}

static cf_t *reduce_dense_row_by_known_pivots_31_bit(
        int64_t *dr,
        dt_t *const *pivs,
        const hl_t dpiv    /* pivot of dense row at the beginning */
        )
{
    hl_t i, j;
    const int64_t mod   = (int64_t)fc;
    const int64_t mod2  = (int64_t)fc * fc;

    for (i = dpiv; i < ncl; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = (int64_t)dr[i];
        const dt_t *dts   = pivs[i];
        const cf_t *cfs   = gbcf[dts[0]];
        for (j = 3; j < dts[1]; ++j) {
            dr[dts[j]]  -=  mul * cfs[j];
            dr[dts[j]]  +=  (dr[dts[j]] >> 63) & mod2;
        }
        for (; j < dts[2]; j += 4) {
            dr[dts[j]]    -=  mul * cfs[j];
            dr[dts[j+1]]  -=  mul * cfs[j+1];
            dr[dts[j+2]]  -=  mul * cfs[j+2];
            dr[dts[j+3]]  -=  mul * cfs[j+3];
            dr[dts[j]]    +=  (dr[dts[j]] >> 63) & mod2;
            dr[dts[j+1]]  +=  (dr[dts[j+1]] >> 63) & mod2;
            dr[dts[j+2]]  +=  (dr[dts[j+2]] >> 63) & mod2;
            dr[dts[j+3]]  +=  (dr[dts[j+3]] >> 63) & mod2;
        }
        dr[i] = 0;
    }

    /* store a dense row for further dense gaussian elimination */
    cf_t *row  = (cf_t *)calloc(
            (unsigned long)(ncr), sizeof(cf_t));

    j = 0;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                j++;
                row[i-ncl]  = (cf_t)dr[i];
            }
        }
    }
    if (j == 0) {
        free(row);
        row = NULL;
    }

    return row;
}

static cf_t *reduce_dense_row_by_dense_new_pivots_17_bit(
        int64_t *dr,
        len_t *pc,
        cf_t *const *pivs
        )
{
    hl_t i, j, k, l;
    len_t np  = -1;
    const int64_t mod = (int64_t)fc;

    /* printf("START DENSE ROW REDUCTION\n"); */
    for (k = 0, i = *pc; i < ncr; ++i) {
        /* printf("dr[%d] = %ld\n", i, dr[i]); */
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == -1) {
                np  = i;
            }
            k++;
            continue;
        }

/*         printf("cfs ");
 *         for (int32_t p = 0; p < ncr-i; ++p) {
 *             printf("%d ", pivs[i][p]);
 *         }
 *         printf("\n");
 *
 *         for (int32_t p = 0; p < ncr; ++p) {
 *             printf("%ld ", dr[p]);
 *         }
 *         printf("\n"); */

        const int64_t mul = mod - dr[i];
        const len_t os    = (ncr - i) % 4;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] +=  mul * pivs[i][l];
        }
        for (; j < ncr; l += 4, j += 4) {
            dr[j]   +=  mul * pivs[i][l];
            dr[j+1] +=  mul * pivs[i][l+1];
            dr[j+2] +=  mul * pivs[i][l+2];
            dr[j+3] +=  mul * pivs[i][l+3];
        }
        /* for (int32_t p = 0; p < ncr; ++p) {
         *     printf("%ld ", dr[p]);
         * }
         * printf("\n--------------------------------\n"); */
    }
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    /* printf("dense reduction step done\n"); */

    cf_t *row = (cf_t *)calloc((unsigned long)(ncr-np), sizeof(cf_t));
    for (i = np; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row(row, np);
    }
    *pc = np;
    return row;
}

static cf_t *reduce_dense_row_by_dense_new_pivots_31_bit(
        int64_t *dr,
        len_t *pc,
        cf_t *const *pivs
        )
{
    hl_t i, j, k, l;
    len_t np  = -1;
    const int64_t mod = (int64_t)fc;
    const int64_t mod2  = (int64_t)fc * fc;

    for (k = 0, i = *pc; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == -1) {
                np  = i;
            }
            k++;
            continue;
        }

        const int64_t mul = (int64_t)dr[i];
        const len_t os    = (ncr - i) % 4;
        for (l = 0, j = i; l < os; ++l, ++j) {
            dr[j] -=  mul * pivs[i][l];
            dr[j] +=  (dr[j] >> 63) & mod2;
        }
        for (; j < ncr; l+=4, j += 4) {
            dr[j]   -=  mul * pivs[i][l];
            dr[j+1] -=  mul * pivs[i][l+1];
            dr[j+2] -=  mul * pivs[i][l+2];
            dr[j+3] -=  mul * pivs[i][l+3];
            dr[j]   +=  (dr[j] >> 63) & mod2;
            dr[j+1] +=  (dr[j+1] >> 63) & mod2;
            dr[j+2] +=  (dr[j+2] >> 63) & mod2;
            dr[j+3] +=  (dr[j+3] >> 63) & mod2;
        }
    }
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    cf_t *row = (cf_t *)calloc((unsigned long)(ncr-np), sizeof(cf_t));
    for (i = np; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row(row, np);
    }
    /* printf("tbr[%d] = %p\n", idx, tbr[idx]); */
    *pc = np;

    return row;
}

static cf_t **sparse_reduced_echelon_form(
        dt_t **mat
        )
{
    len_t i, j, k;
    hl_t sc    = 0;    /* starting column */

    /* we fill in all known lead terms in pivs */
    dt_t **pivs   = (dt_t **)calloc((unsigned long)ncols, sizeof(dt_t *));
    /* unkown pivot rows we have to reduce with the known pivots first */
    dt_t **upivs  = (dt_t **)malloc((unsigned long)nrl * sizeof(dt_t *));

    i = 0;
    j = 1;
    for (i = 0; i < nrows; ++i) {
        if (!pivs[mat[i][3]]) {
            pivs[mat[i][3]]   = mat[i];
            /* printf("pivs %d is set by row %d\n", mat[i][3], i); */
        } else {
            /* shorter rows first */
            upivs[nrl-j]  = mat[i];
            j++;
        }
    }

    /* for (int32_t o = 0; o < nru; ++o) {
     *     printf("%d | %d | %d | %d || %d\n", o, pivs[o][0], pivs[0][1], pivs[o][2], pivs[o][3]);
     * }
     * printf("nru %d | nrl %d | nrows %d || ncl %d | ncr %d | ncols %d\n",
     *         nru, nrl, nrows, ncl, ncr, ncols); */

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)(nthrds * ncols) * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(nthrds) \
    private(i, j, k, sc)
    for (i = 0; i < nrl; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        dt_t *npiv    = upivs[i];
        cf_t *cfs     = gbcf[npiv[0]];
        k = 0;
        do {
            memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
            for (j = 3; j < npiv[1]; ++j) {
                drl[npiv[j]] = (int64_t)cfs[j];
            }
            for (; j < npiv[2]; j += 4) {
                drl[npiv[j]]    = (int64_t)cfs[j];
                drl[npiv[j+1]]  = (int64_t)cfs[j+1];
                drl[npiv[j+2]]  = (int64_t)cfs[j+2];
                drl[npiv[j+3]]  = (int64_t)cfs[j+3];
            }
            sc  = npiv[3];
            free(npiv);
            npiv  = reduce_dense_row_by_known_pivots_sparse(drl, pivs, sc, i);
            if (!npiv) {
                break;
            }
            /* normalize coefficient array
             * NOTE: this has to be done here, otherwise the reduction may
             * lead to wrong results in a parallel computation since other
             * threads might directly use the new pivot once it is synced. */
            if (tmpcf[npiv[0]][3] != 1) {
                normalize_sparse_matrix_row(tmpcf[npiv[0]]);
            }
            k   = __sync_bool_compare_and_swap(&pivs[npiv[3]], NULL, npiv);
            cfs = tmpcf[i];
        } while (!k);
    }
    free(upivs);
    upivs = NULL;

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    npivs = 0; /* number of new pivots */

    dr  = realloc(dr, (unsigned long)ncols * sizeof(uint64_t));
    mat = realloc(mat, (unsigned long)ncr * sizeof(dt_t *));

    /* interreduce new pivots */
    cf_t *cfs;
    dt_t cf_array_pos;
    for (i = (ncols-1); i >= nru; --i) {
        if (pivs[i]) {
            memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
            cfs = tmpcf[pivs[i][0]];
            cf_array_pos  = pivs[i][0];
            sc  = pivs[i][3];
            for (j = 3; j < pivs[i][1]; ++j) {
                dr[pivs[i][j]] = (int64_t)cfs[j];
            }
            for (; j < pivs[i][2]; j += 4) {
                dr[pivs[i][j]]    = (int64_t)cfs[j];
                dr[pivs[i][j+1]]  = (int64_t)cfs[j+1];
                dr[pivs[i][j+2]]  = (int64_t)cfs[j+2];
                dr[pivs[i][j+3]]  = (int64_t)cfs[j+3];
            }
            free(pivs[i]);
            free(cfs);
            pivs[i] = NULL;
            pivs[i] = mat[npivs++] =
                reduce_dense_row_by_known_pivots_sparse(dr, pivs, sc, cf_array_pos);
        }
    }
    free(pivs);
    pivs  = NULL;

    free(dr);
    dr  = NULL;
    mat   = realloc(mat, (unsigned long)npivs * sizeof(hl_t *));
    nrows = nrall = npivs;

    return mat;
}

static cf_t **sparse_AB_CD_linear_algebra(
        dt_t ***matp
        )
{
    dt_t **mat  = *matp;
    len_t i, j;
    hl_t sc    = 0;    /* starting column */

    /* we fill in all known lead terms in pivs */
    dt_t **pivs   = (dt_t **)calloc((unsigned long)nru, sizeof(dt_t *));
    /* unkown pivot rows we have to reduce with the known pivots first */
    dt_t **upivs  = (dt_t **)malloc((unsigned long)nrl * sizeof(dt_t *));
    /* dense rows representing updated D part;
     * after reducing CD part with AB */
    cf_t **drs    = (cf_t **)calloc((unsigned long)nrl, sizeof(cf_t *));

    i = 0;
    j = 1;
    for (i = 0; i < nrows; ++i) {
        if (!pivs[mat[i][3]]) {
            pivs[mat[i][3]]   = mat[i];
            /* printf("pivs %d is set by row %d\n", mat[i][3], i); */
        } else {
            /* shorter rows first */
            upivs[nrl-j]  = mat[i];
            j++;
        }
    }
    free(mat);
    mat   = NULL;
    *matp = mat;

    /* for (int32_t o = 0; o < nru; ++o) {
     *     printf("%d | %d | %d | %d || %d\n", o, pivs[o][0], pivs[0][1], pivs[o][2], pivs[o][3]);
     * }
     * printf("nru %d | nrl %d | nrows %d || ncl %d | ncr %d | ncols %d\n",
     *         nru, nrl, nrows, ncl, ncr, ncols); */

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)(nthrds * ncols) * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(nthrds) \
    private(i, j, sc)
    for (i = 0; i < nrl; ++i) {
        cf_t *cfs = NULL;
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        dt_t *npiv  = upivs[i];
        /* do the reduction */
        memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
        cfs = gbcf[npiv[0]];
        for (j = 3; j < npiv[1]; ++j) {
            drl[npiv[j]] = (int64_t)cfs[j];
        }
        for (; j < npiv[2]; j += 4) {
            drl[npiv[j]]    = (int64_t)cfs[j];
            drl[npiv[j+1]]  = (int64_t)cfs[j+1];
            drl[npiv[j+2]]  = (int64_t)cfs[j+2];
            drl[npiv[j+3]]  = (int64_t)cfs[j+3];
        }
        sc  = npiv[3];
        free(npiv);
        drs[i]  = reduce_dense_row_by_known_pivots(drl, pivs, sc);
    }
    free(dr);
    dr  = NULL;
    free(upivs);
    upivs = NULL;

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }
    free(pivs);
    pivs  = NULL;

    /* remove NULL dense rows */
    npivs = 0; /* number of new pivots */
    for (i = 0; i < nrl; ++i) {
        if (drs[i] != NULL) {
            drs[npivs++]  = drs[i];
        }
    }
    /* for (i = 0; i < npivs; ++i) {
     *     printf("drs[%d] = %p\n", i, drs[i]);
     * } */

    if (npivs == 0) {
        free(drs);
        drs = NULL;
    }

    return drs;
}

static cf_t **interreduce_dense_matrix(
    cf_t **dm
    )
{
    len_t i, j, l;
    dt_t npc, os;
    int64_t *dr = malloc((unsigned long)ncr * sizeof(int64_t));

    for (i = ncr-1; i > -1; --i) {
        if (dm[i]) {
            /* printf("interreduce dm[%d]\n", i); */
            memset(dr, 0, (unsigned long)ncr * sizeof(int64_t));
            npc  = ncr - i;
            os   = npc % 4;
            for (j = i, l = 0; l < os; ++j, ++l) {
                dr[j] = (int64_t)dm[i][l];
            }
            for (; l < npc; j += 4, l += 4) {
                dr[j]   = (int64_t)dm[i][l];
                dr[j+1] = (int64_t)dm[i][l+1];
                dr[j+2] = (int64_t)dm[i][l+2];
                dr[j+3] = (int64_t)dm[i][l+3];
            }
            /* for (j = 0; j < ncr; ++j) {
             *     printf("%ld ", dr[j]);
             * }
             * printf("\n"); */
            free(dm[i]);
            dm[i] = NULL;
            /* start with previous pivot the reduction process, so keep the
             * pivot element as it is */
            dm[i] = reduce_dense_row_by_dense_new_pivots(dr, &i, dm);
            /* printf("tbr to dm[%d] = %p\n", i, dm[i]);
             * for (j = i; j < ncr; ++j) {
             *     printf("%d ", dm[i][j-i]);
             * }
             * printf("\n"); */
        }
    }
    free(dr);
    return dm;
}

static cf_t **probabilistic_sparse_dense_echelon_form(
        dt_t ***matp
        )
{
    len_t i, j, k, l, m;
    dt_t **mat  = *matp;
    /* we fill in all known lead terms in pivs */
    dt_t **pivs   = (dt_t **)calloc((unsigned long)nru, sizeof(dt_t *));
    /* unkown pivot rows we have to reduce with the known pivots first */
    dt_t **upivs  = (dt_t **)malloc((unsigned long)nrl * sizeof(dt_t *));

    j = 1;
    for (i = 0; i < nrows; ++i) {
        if (!pivs[mat[i][3]]) {
            pivs[mat[i][3]]   = mat[i];
            /* printf("pivs %d is set by row %d\n", mat[i][3], i); */
        } else {
            /* shorter rows first */
            upivs[nrl-j]  = mat[i];
            j++;
        }
    }
    free(mat);
    mat   = NULL;
    *matp = mat;
    /* rows already representing new pivots */
    cf_t **nps  = (cf_t **)calloc((unsigned long)ncr, sizeof(cf_t *));

    const int64_t mod2  = (int64_t)fc * fc;

    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(nrl/3)))+1;
    const len_t rem = (nrl % nb == 0) ? 0 : 1;
    const len_t rpb = (nrl / nb) + rem;

    int64_t *dr   = (int64_t *)malloc(
        (unsigned long)(nthrds * ncols) * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (unsigned long)(nthrds * rpb) * sizeof(int64_t));

    /* reduction process to get all possible pivots, no interreduction here */
#pragma omp parallel for num_threads(nthrds) \
    private(i, j, k, l, m) shared(nps, npivs)
    for (i = 0; i < nb; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        int64_t *mull = mul + (omp_get_thread_num() * rpb);
        const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);
        if (nrbl != 0) {
            dt_t *npiv;
            cf_t *cfs;
            cf_t *tmp;
            dt_t npc;
            dt_t os;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                npc = 0;
                os  = ncr % 4;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    mull[j] = (int64_t)rand() % fc;
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    npiv  = upivs[m];
                    cfs   = gbcf[npiv[0]];
                    for (l = 3; l < npiv[1]; ++l) {
                        drl[npiv[l]]  -=  mull[k] * cfs[l];
                        drl[npiv[l]]  +=  (drl[npiv[l]] >> 63) & mod2;
                    }
                    for (; l < npiv[2]; l += 4) {
                        drl[npiv[l]]    -=  mull[k] * cfs[l];
                        drl[npiv[l]]    +=  (drl[npiv[l]] >> 63) & mod2;
                        drl[npiv[l+1]]  -=  mull[k] * cfs[l+1];
                        drl[npiv[l+1]]  +=  (drl[npiv[l+1]] >> 63) & mod2;
                        drl[npiv[l+2]]  -=  mull[k] * cfs[l+2];
                        drl[npiv[l+2]]  +=  (drl[npiv[l+2]] >> 63) & mod2;
                        drl[npiv[l+3]]  -=  mull[k] * cfs[l+3];
                        drl[npiv[l+3]]  +=  (drl[npiv[l+3]] >> 63) & mod2;
                    }
                }
                k   = 0;
                npc = 0;
                /* do the reduction */
                do {
                    npc = reduce_dense_row_by_all_pivots(drl, &tmp, pivs, nps, npc);
                    if (npc == -1) {
                        bctr  = nrbl;
                        break;
                    }
                    k = __sync_bool_compare_and_swap(&nps[npc], NULL, tmp);
                    /* some other thread has already added a pivot so we have to
                    * recall the dense reduction process */
                    if (!k) {
                        memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
                        os  = (ncr-npc) % 4;
                        for (l = 0, j = npc+ncl; l < os; ++l, ++j) {
                            drl[j]  = tmp[l];
                        }
                        for (; l < ncr-npc; l += 4, j += 4) {
                            drl[j]    = tmp[l];
                            drl[j+1]  = tmp[l+1];
                            drl[j+2]  = tmp[l+2];
                            drl[j+3]  = tmp[l+3];
                        }
                    }
                } while (!k);
                bctr++;
            }
            for (j = i*rpb; j < nbl; ++j) {
                free(upivs[j]);
                upivs[j]  = NULL;
            }
        }
    }
    /* count number of pivots */
    const len_t os  = ncr % 4;
    for (i = 0; i < os; ++i) {
        if (nps[i] != NULL) {
            npivs++;
        }
    }
    for (; i < ncr; i += 4) {
        if (nps[i] != NULL) {
            npivs++;
        }
        if (nps[i+1] != NULL) {
            npivs++;
        }
        if (nps[i+2] != NULL) {
            npivs++;
        }
        if (nps[i+3] != NULL) {
            npivs++;
        }
    }


    for (i = 0; i < nru; ++i) {
        free(pivs[i]);
    }
    free(pivs);
    pivs  = NULL;
    free(mul);
    mul = NULL;
    free(upivs);
    upivs = NULL;
    free(dr);
    dr  = NULL;

    return nps;
}

static cf_t **exact_dense_linear_algebra(
        cf_t **dm,
        const len_t nr
        )
{
    len_t i, j, k, l;
    /* rows already representing new pivots */
    cf_t **nps  = (cf_t **)calloc((unsigned long)ncr, sizeof(cf_t *));
    /* rows to be further reduced */
    cf_t **tbr  = (cf_t **)calloc((unsigned long)nr, sizeof(cf_t *));
    int64_t *dr   = (int64_t *)malloc(
            (unsigned long)(nthrds * ncr) * sizeof(int64_t));

    /* separate rows already representing new pivots and rows to
     * be further reduced by these new pivots */
    j     = 0;
    npivs = 0;
    for (i = 0; i < nr; ++i) {
        /* printf("dense row[%d] ", i);
         * for (int32_t p = 0; p < ncr; ++p) {
         *     printf("%d ", dm[i][p]);
         * }
         * printf("\n"); */
        if (dm[i] != NULL) {
            k = 0;
            while (dm[i][k] == 0) {
                ++k;
            }
            /* printf("pc %d -> %p\n", k, nps[k]); */
            if (nps[k] == NULL) {
                /* we have a pivot, cut the dense row down to start
                 * at the first nonzero entry */
                /* printf("ncr - k = %d - %d\n", ncr, k); */
                memmove(dm[i], dm[i]+k, (unsigned long)(ncr-k) * sizeof(cf_t));
                dm[i] = realloc(dm[i], (unsigned long)(ncr-k) * sizeof(cf_t));
                nps[k] = dm[i];
                /* printf("nps[%d] = ", k); */
                if (nps[k][0] != 1) {
                    nps[k]  = normalize_dense_matrix_row(nps[k], k);
                }
                /* npivs++; */
            } else {
                /* printf("set tbr[%d]\n", j); */
                tbr[j++]  = dm[i];
            }
        }
    }
    free(dm);
    dm  = NULL;

    const len_t ntr = j;
    tbr = realloc(tbr, (unsigned long)ntr * sizeof(cf_t *));
    /* offset modulo 4 for loop unrolling, +1 due to storing the first
     * nonzero entry at the first position */

    /* reduction process to get all possible pivots, no interreduction here */
#pragma omp parallel for num_threads(nthrds) \
    private(i, j, k, l) shared(nps, tbr, npivs)
    for (i = 0; i < ntr; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncr);
        memset(drl, 0, (unsigned long)ncr * sizeof(int64_t));
        /* printf("tbr[%d] = %p\n", i, tbr[i]); */
        dt_t npc  = 0;
        dt_t os   = 0;
        cf_t *npiv  = tbr[i];;
        do {
            os   = (ncr-npc) % 4;
            memset(drl, 0, (unsigned long)ncr * sizeof(int64_t));
            for (l = 0, j = npc; l < os; ++l, ++j) {
                drl[j]  = (int64_t)npiv[l];
            }
            for (; j < ncr; l += 4, j += 4) {
                drl[j]    = (int64_t)npiv[l];
                drl[j+1]  = (int64_t)npiv[l+1];
                drl[j+2]  = (int64_t)npiv[l+2];
                drl[j+3]  = (int64_t)npiv[l+3];
            }
            free(npiv);
            npiv = NULL;
            npiv = reduce_dense_row_by_dense_new_pivots(drl, &npc, nps);
            if (npc == -1) {
                break;
            }
            k = __sync_bool_compare_and_swap(&nps[npc], NULL, npiv);
            /* some other thread has already added a pivot so we have to
             * recall the dense reduction process */
        } while (!k);
    }
    /* count number of pivots */
    const len_t os  = ncr % 4;
    for (i = 0; i < os; ++i) {
        if (nps[i] != NULL) {
            npivs++;
        }
    }
    for (; i < ncr; i += 4) {
        if (nps[i] != NULL) {
            npivs++;
        }
        if (nps[i+1] != NULL) {
            npivs++;
        }
        if (nps[i+2] != NULL) {
            npivs++;
        }
        if (nps[i+3] != NULL) {
            npivs++;
        }
    }
    free(tbr);
    free(dr);

    return nps;
}

static dt_t **exact_sparse_linear_algebra(
        dt_t **mat
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* allocate temporary storage space for sparse
     * coefficients of new pivot rows */
    tmpcf = realloc(tmpcf, (unsigned long)nrl * sizeof(cf_t *));
    mat   = sparse_reduced_echelon_form(mat);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    la_ctime  +=  ct1 - ct0;
    la_rtime  +=  rt1 - rt0;

    num_zerored += (nrl - npivs);
    GB_DEBUG(LADBG, "%7d new %7d zero - %9.3f sec", npivs, nrl-npivs, rt1-rt0);
    
    return mat;
}

static dt_t **convert_to_sparse_matrix_rows(
        cf_t * const *dm,
        dt_t **mat
        )
{
    if (npivs == 0) {
        return NULL;
    }

    len_t i, j, k, l;
    cf_t *cfs;
    dt_t *dts;

    mat   = realloc(mat, (unsigned long)npivs * sizeof(dt_t *));
    tmpcf = realloc(tmpcf, (unsigned long)npivs * sizeof(cf_t *));

    l = 0;
    for (i = ncr-1; i > -1; --i) {
        if (dm[i] != NULL) {
            cfs = malloc((unsigned long)(ncr-i+3) * sizeof(cf_t));
            dts = malloc((unsigned long)(ncr-i+3) * sizeof(dt_t));
            const dt_t len    = ncr-i;
            const dt_t os     = len % 4;
            const dt_t shift  = ncl+i;

            for (k = 3, j = 0; j < os; ++j) {
                if (dm[i][j] != 0) {
                    cfs[k]    = dm[i][j];
                    dts[k++]  = j+shift;
                }
            }
            for (; j < len; j += 4) {
                if (dm[i][j] != 0) {
                    cfs[k]    = dm[i][j];
                    dts[k++]  = j+shift;
                }
                if (dm[i][j+1] != 0) {
                    cfs[k]    = dm[i][j+1];
                    dts[k++]  = j+1+shift;
                }
                if (dm[i][j+2] != 0) {
                    cfs[k]    = dm[i][j+2];
                    dts[k++]  = j+2+shift;
                }
                if (dm[i][j+3] != 0) {
                    cfs[k]    = dm[i][j+3];
                    dts[k++]  = j+3+shift;
                }
            }

            /* store meta data in first entries */
            dts[0]  = l; /* position of coefficient array in tmpcf */
            dts[1]  = ((k-3) % 4) + 3;
            dts[2]  = k;
            cfs[0]  = 0;
            cfs[1]  = dts[1];
            cfs[2]  = dts[2];

            /* adjust memory usage */
            dts = realloc(dts, (unsigned long)k * sizeof(dt_t));
            cfs = realloc(cfs, (unsigned long)k * sizeof(cf_t));

            /* link to basis */
            mat[l]    = dts;
            tmpcf[l]  = cfs;
            l++;
        }
    }

    return mat;
}

/* NOTE: this note is about the different linear algebra implementations:
 * exact and probabilistic linear algebra differ only in the last,
 * dense reduction step: the reduction of CD via AB is sparse and
 * the same for both. this generates a dense D' part which is then
 * either reduced via exact linear algebra or via probabilistic
 * linear algebra */
static dt_t **exact_sparse_dense_linear_algebra(
        dt_t **mat
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* generate updated dense D part via reduction of CD with AB */
    cf_t **dm;
    dm  = sparse_AB_CD_linear_algebra(&mat);
    if (npivs > 0) {      
        dm  = exact_dense_linear_algebra(dm, npivs);
        dm  = interreduce_dense_matrix(dm);
    }

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    mat = convert_to_sparse_matrix_rows(dm, mat);

    /* free dm */
    if (dm) {
        for (i = 0; i < ncr; ++i) {
            free(dm[i]);
        }
        free(dm);
        dm  = NULL;
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    la_ctime  +=  ct1 - ct0;
    la_rtime  +=  rt1 - rt0;

    num_zerored += (nrl - npivs);
    GB_DEBUG(LADBG, "%7d new %7d zero - %9.3f sec", npivs, nrl-npivs, rt1-rt0);
    
    return mat;
}

static dt_t **probabilistic_sparse_dense_linear_algebra(
        dt_t **mat
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* generate updated dense D part via reduction of CD with AB */
    cf_t **dm;
    npivs = 0;
    dm    = probabilistic_sparse_dense_echelon_form(&mat);
    dm    = interreduce_dense_matrix(dm);

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    mat = convert_to_sparse_matrix_rows(dm, mat);

    /* free dm */
    if (dm) {
        for (i = 0; i < ncr; ++i) {
            free(dm[i]);
        }
        free(dm);
        dm  = NULL;
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    la_ctime  +=  ct1 - ct0;
    la_rtime  +=  rt1 - rt0;

    num_zerored += (nrl - npivs);
    GB_DEBUG(LADBG, "%7d new %7d zero - %9.3f sec", npivs, nrl-npivs, rt1-rt0);
    
    return mat;
}
