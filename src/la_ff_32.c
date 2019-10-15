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
 * \file la_ff.c
 * \brief Implementation of finite field linear algebra.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "data.h"

static void interreduce_matrix_rows_ff_32(
        mat_t *mat,
        bs_t *bs,
        const stat_t *st
        )
{
    len_t i, j, k;

    const len_t nrows = mat->nr;
    const len_t ncols = mat->nc;

    hm_t **pivs = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    for (i = 0; i < nrows; ++i) {
        pivs[mat->r[i][3]]  = mat->r[i];
    }
    mat->cf_32  = realloc(mat->cf_32,
            (unsigned long)nrows * sizeof(cf32_t *));
    int64_t *dr = (int64_t *)malloc((unsigned long)ncols * sizeof(int64_t));
    /* interreduce new pivots */
    cf32_t *cfs;
    /* starting column, coefficient array position in tmpcf */
    hm_t sc, cfp;
    k = nrows - 1;
    for (i = (ncols-1); i >= 0; --i) {
        if (pivs[i] != NULL) {
            memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
            cfs = bs->cf_32[pivs[i][0]];
            cfp = pivs[i][0]; // this is the position in bs
            const len_t os  = pivs[i][1];
            const len_t len = pivs[i][2];
            const hm_t * const ds = pivs[i] + 3;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                dr[ds[j]] = (int64_t)cfs[j];
            }
            for (; j < len; j += 4) {
                dr[ds[j]]   = (int64_t)cfs[j];
                dr[ds[j+1]] = (int64_t)cfs[j+1];
                dr[ds[j+2]] = (int64_t)cfs[j+2];
                dr[ds[j+3]] = (int64_t)cfs[j+3];
            }
            free(pivs[i]);
            free(cfs);
            pivs[i] = NULL;
            pivs[i] = mat->r[k--] =
                reduce_dense_row_by_known_pivots_sparse_ff_32(
                        dr, mat, bs, pivs, sc, cfp, st->fc);
        }
    }
    mat->np = nrows;
    free(pivs);
    free(dr);
}

static inline cf32_t *normalize_dense_matrix_row_ff_32(
        cf32_t *row,
        const hm_t len,
        const int32_t fc
        )
{
    len_t i;

    const hm_t os   = len % 4;
    int64_t tmp1, tmp2, tmp3, tmp4;

    const int32_t inv = mod_p_inverse_32((int32_t)row[0], (int32_t)fc);

    for (i = 1; i < os; ++i) {
        tmp1    =   ((int64_t)row[i] * inv) % fc;
        tmp1    +=  (tmp1 >> 63) & fc;
        row[i]  =   (cf32_t)tmp1;
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
        row[i]    =   (cf32_t)tmp1;
        row[i+1]  =   (cf32_t)tmp2;
        row[i+2]  =   (cf32_t)tmp3;
        row[i+3]  =   (cf32_t)tmp4;
    }
    row[0]  = 1;

    return row;
}

static inline cf32_t *normalize_sparse_matrix_row_ff_32(
        cf32_t *row,
        const len_t os,
        const len_t len,
        const int32_t fc
        )
{
    len_t i;

    int64_t tmp1, tmp2, tmp3, tmp4;

    const int32_t inv = mod_p_inverse_32((int32_t)row[0], (int32_t)fc);

    for (i = 0; i < os; ++i) {
        tmp1    =   ((int64_t)row[i] * inv) % fc;
        tmp1    +=  (tmp1 >> 63) & fc;
        row[i]  =   (cf32_t)tmp1;
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
        row[i]    =   (cf32_t)tmp1;
        row[i+1]  =   (cf32_t)tmp2;
        row[i+2]  =   (cf32_t)tmp3;
        row[i+3]  =   (cf32_t)tmp4;
    }
    row[0]  = 1;

    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_17_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hl_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const int32_t fc
        )
{
    hl_t i, j, k;
    hm_t *dts;
    cf32_t *cfs;
    len_t np  = 0;
    const int64_t mod           = (int64_t)fc;
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf32_t * const * const mcf  = mat->cf_32;

    k = 0;
    for (i = dpiv; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == 0) {
                np  = i;
            }
            k++;
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = mod - dr[i];
        dts   = pivs[i];
        if (i < ncl) {
            cfs   = bs->cf_32[dts[0]];
        } else {
            cfs   = mcf[dts[0]];
        }
        const len_t os  = dts[1];
        const len_t len = dts[2];
        const hm_t * const ds  = dts + 3;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] +=  mul * cfs[j];
        }
        for (; j < len; j += 4) {
            dr[ds[j]]   +=  mul * cfs[j];
            dr[ds[j+1]] +=  mul * cfs[j+1];
            dr[ds[j+2]] +=  mul * cfs[j+2];
            dr[ds[j+3]] +=  mul * cfs[j+3];
        }
        dr[i] = 0;
    }
    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((unsigned long)(ncols-np+3) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(ncols-np) * sizeof(cf32_t));
    j = 0;
    hm_t *rs = row + 3;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                rs[j] = (hm_t)i;
                cf[j] = (cf32_t)dr[i];
                j++;
            }
        }
    }
    if (j == 0) {
        free(row);
        row = NULL;
        free(cf);
        cf  = NULL;
    } else {
        row     = realloc(row, (unsigned long)(j+3) * sizeof(hm_t));
        cf      = realloc(cf, (unsigned long)j * sizeof(cf32_t));
        row[0]  = tmp_pos;
        row[1]  = j % 4;
        row[2]  = j;
        mat->cf_32[tmp_pos]  = cf;
    }
    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_31_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hl_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos, /* position of new coeffs array in tmpcf */
        const int32_t fc
        )
{
    hl_t i, j, k;
    cf32_t *cfs;
    hm_t *dts;
    len_t np  = 0;
    const int64_t mod           = (int64_t)fc;
    const int64_t mod2          = (int64_t)fc * fc;
    const len_t ncols           = mat->nc;
    const len_t ncl             = mat->ncl;
    cf32_t * const * const mcf  = mat->cf_32;

    k = 0;
    for (i = dpiv; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        if (dr[i] == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == 0) {
                np  = i;
            }
            k++;
            continue;
        }

        /* found reducer row, get multiplier */
        const int64_t mul = (int64_t)dr[i];
        dts   = pivs[i];
        if (i < ncl) {
            cfs   = bs->cf_32[dts[0]];
        } else {
            cfs   = mcf[dts[0]];
        }

        const len_t os  = dts[1];
        const len_t len = dts[2];
        const hm_t * const ds = dts + 3;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 4) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j+1]] -=  mul * cfs[j+1];
            dr[ds[j+2]] -=  mul * cfs[j+2];
            dr[ds[j+3]] -=  mul * cfs[j+3];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
            dr[ds[j+1]] +=  (dr[ds[j+1]] >> 63) & mod2;
            dr[ds[j+2]] +=  (dr[ds[j+2]] >> 63) & mod2;
            dr[ds[j+3]] +=  (dr[ds[j+3]] >> 63) & mod2;
        }
        dr[i] = 0;
    }

    if (k == 0) {
        return NULL;
    }

    hm_t *row   = (hm_t *)malloc((unsigned long)(ncols-np+3) * sizeof(hm_t));
    cf32_t *cf  = (cf32_t *)malloc((unsigned long)(ncols-np) * sizeof(cf32_t));
    j = 0;
    hm_t *rs  = row + 3;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                rs[j] = (hm_t)i;
                cf[j] = (cf32_t)dr[i];
                j++;
            }
        }
    }
    if (j == 0) {
        free(row);
        row = NULL;
        free(cf);
        cf  = NULL;
    } else {
        row     = realloc(row, (unsigned long)(j+3) * sizeof(hm_t));
        cf      = realloc(cf, (unsigned long)j * sizeof(cf32_t));
        row[0]  = tmp_pos;
        row[1]  = j % 4;
        row[2]  = j;
        mat->cf_32[tmp_pos]  = cf;
    }
    return row;
}

static cf32_t *reduce_dense_row_by_all_pivots_17_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t *const * const pivs,
        cf32_t *const * const dpivs,
        const int32_t fc
        )
{
    hl_t i, j, k, l;
    const int64_t mod = (int64_t)fc;
    len_t np  = -1;
    cf32_t *red;

    const len_t ncl   = mat->ncl;
    const len_t ncols = mat->nc;

    /* step 1: reduce by sparse known pivots */
    for (i = *pc; i < ncl; ++i) {
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
        const int64_t mul = mod - dr[i];
        const cf32_t *cfs = bs->cf_32[pivs[i][0]];
        const len_t os    = pivs[i][1];
        const len_t len   = pivs[i][2];
        const hm_t * const ds  = pivs[i] + 3;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] +=  mul * cfs[j];
        }
        for (; j < len; j += 4) {
            dr[ds[j]]   +=  mul * cfs[j];
            dr[ds[j+1]] +=  mul * cfs[j+1];
            dr[ds[j+2]] +=  mul * cfs[j+2];
            dr[ds[j+3]] +=  mul * cfs[j+3];
        }
        dr[i] = 0;
    }
    k = 0;
    /* step 2: reduce by new dense pivots */
    for (i = ncl; i < ncols; ++i) {
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
    }
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    cf32_t *row = (cf32_t *)calloc((unsigned long)(ncols-np), sizeof(cf32_t));
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncols-np, fc);
    }

    *pc = np - ncl;
    return row;
}

static cf32_t *reduce_dense_row_by_all_pivots_31_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t * const * const pivs,
        cf32_t * const * const dpivs,
        const int32_t fc
        )
{
    hl_t i, j, k, l;
    const int64_t mod   = (int64_t)fc;
    const int64_t mod2  = (int64_t)fc * fc;
    len_t np  = -1;
    cf32_t *red;

    const len_t ncl   = mat->ncl;
    const len_t ncols = mat->nc;

    /* step 1: reduce by sparse known pivots */
    for (i = *pc; i < ncl; ++i) {
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
        const cf32_t *cfs = bs->cf_32[pivs[i][0]];
        const len_t os    = pivs[i][1];
        const len_t len   = pivs[i][2];
        const hm_t * const ds = pivs[i] + 3;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 4) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j+1]] -=  mul * cfs[j+1];
            dr[ds[j+2]] -=  mul * cfs[j+2];
            dr[ds[j+3]] -=  mul * cfs[j+3];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
            dr[ds[j+1]] +=  (dr[ds[j+1]] >> 63) & mod2;
            dr[ds[j+2]] +=  (dr[ds[j+2]] >> 63) & mod2;
            dr[ds[j+3]] +=  (dr[ds[j+3]] >> 63) & mod2;
        }
        dr[i] = 0;
    }
    k = 0;
    /* step 2: reduce by new dense pivots */
    for (i = ncl; i < ncols; ++i) {
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
        *pc = -1;
        return NULL;
    }

    cf32_t *row = (cf32_t *)calloc((unsigned long)(ncols-np), sizeof(cf32_t));
    for (i = np; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncols-np, fc);
    }

    *pc = np - ncl;

    return row;
}

static cf32_t *reduce_dense_row_by_old_pivots_17_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const * const pivs,
        const hl_t dpiv,    /* pivot of dense row at the beginning */
        const int32_t fc
        )
{
    hl_t i, j;
    const int64_t mod = (int64_t)fc;
    const len_t ncols = mat->nc;
    const len_t ncl   = mat->ncl;
    const len_t ncr   = mat->ncr;

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
        const int64_t mul = mod - dr[i];
        const cf32_t *cfs = bs->cf_32[pivs[i][0]];
        const len_t os    = pivs[i][1];
        const len_t len   = pivs[i][2];
        const hm_t * const ds = pivs[i] + 3;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] +=  mul * cfs[j];
        }
        for (; j < len; j += 4) {
            dr[ds[j]]   +=  mul * cfs[j];
            dr[ds[j+1]] +=  mul * cfs[j+1];
            dr[ds[j+2]] +=  mul * cfs[j+2];
            dr[ds[j+3]] +=  mul * cfs[j+3];
        }
        dr[i] = 0;
    }

    /* store a dense row for further dense gaussian elimination */
    cf32_t *row  = (cf32_t *)calloc(
            (unsigned long)(ncr), sizeof(cf32_t));
    j = 0;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                j++;
                row[i-ncl]  = (cf32_t)dr[i];
            }
        }
    }
    if (j == 0) {
        free(row);
        row = NULL;
    }
    return row;
}

static cf32_t *reduce_dense_row_by_old_pivots_31_bit(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hl_t dpiv,  /* pivot of dense row at the beginning */
        const int32_t fc
        )
{
    hl_t i, j;
    const int64_t mod   = (int64_t)fc;
    const int64_t mod2  = (int64_t)fc * fc;
    const len_t ncols   = mat->nc;
    const len_t ncl     = mat->ncl;
    const len_t ncr     = mat->ncr;

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
        const cf32_t *cfs = bs->cf_32[pivs[i][0]];
        const len_t os    = pivs[i][1];
        const len_t len   = pivs[i][2];
        const hm_t * const ds = pivs[i] + 3;
        for (j = 0; j < os; ++j) {
            dr[ds[j]] -=  mul * cfs[j];
            dr[ds[j]] +=  (dr[ds[j]] >> 63) & mod2;
        }
        for (; j < len; j += 4) {
            dr[ds[j]]   -=  mul * cfs[j];
            dr[ds[j+1]] -=  mul * cfs[j+1];
            dr[ds[j+2]] -=  mul * cfs[j+2];
            dr[ds[j+3]] -=  mul * cfs[j+3];
            dr[ds[j]]   +=  (dr[ds[j]] >> 63) & mod2;
            dr[ds[j+1]] +=  (dr[ds[j+1]] >> 63) & mod2;
            dr[ds[j+2]] +=  (dr[ds[j+2]] >> 63) & mod2;
            dr[ds[j+3]] +=  (dr[ds[j+3]] >> 63) & mod2;
        }
        dr[i] = 0;
    }

    /* store a dense row for further dense gaussian elimination */
    cf32_t *row  = (cf32_t *)calloc(
            (unsigned long)(ncr), sizeof(cf32_t));

    j = 0;
    for (i = ncl; i < ncols; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
            if (dr[i] != 0) {
                j++;
                row[i-ncl]  = (cf32_t)dr[i];
            }
        }
    }
    if (j == 0) {
        free(row);
        row = NULL;
    }

    return row;
}

static cf32_t *reduce_dense_row_by_dense_new_pivots_17_bit(
        int64_t *dr,
        len_t *pc,
        cf32_t * const * const pivs,
        const len_t ncr,
        const int32_t fc
        )
{
    hl_t i, j, k, l;
    len_t np  = -1;
    const int64_t mod = (int64_t)fc;

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
    }
    if (k == 0) {
        *pc = -1;
        return NULL;
    }

    cf32_t *row = (cf32_t *)calloc((unsigned long)(ncr-np), sizeof(cf32_t));
    for (i = np; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncr-np, fc);
    }
    *pc = np;
    return row;
}

static cf32_t *reduce_dense_row_by_dense_new_pivots_31_bit(
        int64_t *dr,
        len_t *pc,
        cf32_t * const * const pivs,
        const len_t ncr,
        const int32_t fc
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

    cf32_t *row = (cf32_t *)calloc((unsigned long)(ncr-np), sizeof(cf32_t));
    for (i = np; i < ncr; ++i) {
        if (dr[i] != 0) {
            dr[i] = dr[i] % mod;
        }
        row[i-np]  = (cf32_t)dr[i];
    }
    if (row[0] != 1) {
        row = normalize_dense_matrix_row_ff_32(row, ncr-np, fc);
    }
    *pc = np;

    return row;
}

static void probabilistic_sparse_reduced_echelon_form_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        const stat_t * const st
        )
{
    len_t i, j, k, l, m;

    const len_t nrows = mat->nr;
    const len_t ncols = mat->nc;
    const len_t nru   = mat->nru;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    hm_t **rows = mat->r;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = (hm_t **)malloc((unsigned long)nrl * sizeof(hm_t *));

    i = 0;
    j = 1;
    for (i = 0; i < nrows; ++i) {
        if (!pivs[rows[i][3]]) {
            pivs[rows[i][3]]   = rows[i];
        } else {
            /* shorter rows first */
            upivs[nrl-j]  = rows[i];
            j++;
        }
    }

    const int32_t fc    = st->fc;
    const int64_t mod2  = (int64_t)fc * fc;

    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(nrl/3)))+1;
    const len_t rem = (nrl % nb == 0) ? 0 : 1;
    const len_t rpb = (nrl / nb) + rem;

    int64_t *dr   = (int64_t *)malloc(
        (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (unsigned long)(st->nthrds * rpb) * sizeof(int64_t));

    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
    const int distribution = nrl/st->nthrds/100 == 0 ? 1 : nrl/st->nthrds/100;
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l, m) \
    schedule(static, distribution)
    for (i = 0; i < nb; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        int64_t *mull = mul + (omp_get_thread_num() * rpb);
        const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);
        if (nrbl != 0) {
            hm_t *npiv  = NULL;
            cf32_t *cfs;
            /* starting column, offset, coefficient array position in tmpcf */
            hm_t sc, cfp;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                cfp = bctr + i*rpb;
                sc  = 0;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    mull[j] = (int64_t)rand() % fc;
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    npiv  = upivs[m];
                    cfs   = bs->cf_32[npiv[0]];
                    const len_t os  = npiv[1];
                    const len_t len = npiv[2];
                    const hm_t * const ds = npiv + 3;
                    sc    = sc < ds[0] ? sc : ds[0];
                    for (l = 0; l < os; ++l) {
                        drl[ds[l]]  -=  mull[k] * cfs[l];
                        drl[ds[l]]  +=  (drl[ds[l]] >> 63) & mod2;
                    }
                    for (; l < len; l += 4) {
                        drl[ds[l]]    -=  mull[k] * cfs[l];
                        drl[ds[l]]    +=  (drl[ds[l]] >> 63) & mod2;
                        drl[ds[l+1]]  -=  mull[k] * cfs[l+1];
                        drl[ds[l+1]]  +=  (drl[ds[l+1]] >> 63) & mod2;
                        drl[ds[l+2]]  -=  mull[k] * cfs[l+2];
                        drl[ds[l+2]]  +=  (drl[ds[l+2]] >> 63) & mod2;
                        drl[ds[l+3]]  -=  mull[k] * cfs[l+3];
                        drl[ds[l+3]]  +=  (drl[ds[l+3]] >> 63) & mod2;
                    }
                }
                k     = 0;
                cfs   = NULL;
                npiv  = NULL;
                /* do the reduction */
                do {
                    free(cfs);
                    cfs = NULL;
                    free(npiv);
                    npiv  = NULL;
                    npiv  = reduce_dense_row_by_known_pivots_sparse_ff_32(
                            drl, mat, bs, pivs, sc, cfp, st->fc);
                    if (!npiv) {
                        bctr  = nrbl;
                        break;
                    }
                    /* normalize coefficient array
                    * NOTE: this has to be done here, otherwise the reduction may
                    * lead to wrong results in a parallel computation since other
                    * threads might directly use the new pivot once it is synced. */
                    if (mat->cf_32[npiv[0]][0] != 1) {
                        normalize_sparse_matrix_row_ff_32(
                                mat->cf_32[npiv[0]], npiv[1], npiv[2], st->fc);
                    }
                    cfs = mat->cf_32[npiv[0]];
                    sc  = npiv[3];
                    k   = __sync_bool_compare_and_swap(&pivs[npiv[3]], NULL, npiv);
                } while (!k);
                bctr++;
            }
            for (j = i*rpb; j < nbl; ++j) {
                free(upivs[j]);
                upivs[j]  = NULL;
            }
        }
    }
    free(upivs);
    upivs = NULL;
    free(mul);
    mul   = NULL;

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    dr    = realloc(dr, (unsigned long)ncols * sizeof(uint64_t));
    rows  = realloc(rows, (unsigned long)ncr * sizeof(hm_t *));

    /* interreduce new pivots */
    cf32_t *cfs;
    /* starting column, coefficient array position in tmpcf */
    hm_t sc, cfp;
    for (i = (ncols-1); i >= nru; --i) {
        if (pivs[i]) {
            memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
            cfs = mat->cf_32[pivs[i][0]];
            cfp = pivs[i][0];
            const len_t os  = pivs[i][1];
            const len_t len = pivs[i][2];
            const hm_t * const ds = pivs[i] + 3;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                dr[ds[j]] = (int64_t)cfs[j];
            }
            for (; j < len; j += 4) {
                dr[ds[j]]   = (int64_t)cfs[j];
                dr[ds[j+1]] = (int64_t)cfs[j+1];
                dr[ds[j+2]] = (int64_t)cfs[j+2];
                dr[ds[j+3]] = (int64_t)cfs[j+3];
            }
            free(pivs[i]);
            free(cfs);
            pivs[i] = NULL;
            pivs[i] = rows[npivs++] =
                reduce_dense_row_by_known_pivots_sparse_ff_32(
                        dr, mat, bs, pivs, sc, cfp, st->fc);
        }
    }
    free(pivs);
    pivs  = NULL;

    free(dr);
    dr  = NULL;

    rows    = realloc(rows, (unsigned long)npivs * sizeof(hl_t *));
    mat->r  = rows;
    mat->np = mat->nr = mat->sz = npivs;
}

static void exact_sparse_reduced_echelon_form_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        const stat_t * const st
        )
{
    len_t i, j, k;
    hl_t sc    = 0;    /* starting column */

    const len_t nrows = mat->nr;
    const len_t ncols = mat->nc;
    const len_t nru   = mat->nru;
    const len_t nrl   = mat->nrl;
    const len_t ncr   = mat->ncr;
    const len_t ncl   = mat->ncl;

    hm_t **rows = mat->r;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)ncols, sizeof(hm_t *));
    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = (hm_t **)malloc((unsigned long)nrl * sizeof(hm_t *));

    i = 0;
    j = 1;
    for (i = 0; i < nrows; ++i) {
        if (!pivs[rows[i][3]]) {
            pivs[rows[i][3]]   = rows[i];
        } else {
            /* shorter rows first */
            upivs[nrl-j]  = rows[i];
            j++;
        }
    }

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
    const int distribution = nrl/st->nthrds/100 == 0 ? 1 : nrl/st->nthrds/100;
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, sc) \
    schedule(static, distribution)
    for (i = 0; i < nrl; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        hm_t *npiv      = upivs[i];
        cf32_t *cfs       = bs->cf_32[npiv[0]];
        const len_t os  = npiv[1];
        const len_t len = npiv[2];
        const hm_t * const ds = npiv + 3;
        k = 0;
        memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
        for (j = 0; j < os; ++j) {
            drl[ds[j]]  = (int64_t)cfs[j];
        }
        for (; j < len; j += 4) {
            drl[ds[j]]    = (int64_t)cfs[j];
            drl[ds[j+1]]  = (int64_t)cfs[j+1];
            drl[ds[j+2]]  = (int64_t)cfs[j+2];
            drl[ds[j+3]]  = (int64_t)cfs[j+3];
        }
        cfs = NULL;
        do {
            sc  = npiv[3];
            free(npiv);
            free(cfs);
            npiv  = reduce_dense_row_by_known_pivots_sparse_ff_32(
                    drl, mat, bs, pivs, sc, i, st->fc);
            if (!npiv) {
                break;
            }
            /* normalize coefficient array
             * NOTE: this has to be done here, otherwise the reduction may
             * lead to wrong results in a parallel computation since other
             * threads might directly use the new pivot once it is synced. */
            if (mat->cf_32[npiv[0]][0] != 1) {
                normalize_sparse_matrix_row_ff_32(
                        mat->cf_32[npiv[0]], npiv[1], npiv[2], st->fc);
            }
            k   = __sync_bool_compare_and_swap(&pivs[npiv[3]], NULL, npiv);
            cfs = mat->cf_32[npiv[0]];
        } while (!k);
    }
    free(upivs);
    upivs = NULL;

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    dr      = realloc(dr, (unsigned long)ncols * sizeof(uint64_t));
    mat->r  = realloc(mat->r, (unsigned long)ncr * sizeof(hm_t *));
    rows    = mat->r;

    /* interreduce new pivots */
    cf32_t *cfs;
    hm_t cf_array_pos;
    for (i = (ncols-1); i >= nru; --i) {
        if (pivs[i]) {
            memset(dr, 0, (unsigned long)ncols * sizeof(int64_t));
            cfs = mat->cf_32[pivs[i][0]];
            cf_array_pos    = pivs[i][0];
            const len_t os  = pivs[i][1];
            const len_t len = pivs[i][2];
            const hm_t * const ds = pivs[i] + 3;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                dr[ds[j]] = (int64_t)cfs[j];
            }
            for (; j < len; j += 4) {
                dr[ds[j]]    = (int64_t)cfs[j];
                dr[ds[j+1]]  = (int64_t)cfs[j+1];
                dr[ds[j+2]]  = (int64_t)cfs[j+2];
                dr[ds[j+3]]  = (int64_t)cfs[j+3];
            }
            free(pivs[i]);
            free(cfs);
            pivs[i] = NULL;
            pivs[i] = rows[npivs++] =
                reduce_dense_row_by_known_pivots_sparse_ff_32(
                        dr, mat, bs, pivs, sc, cf_array_pos, st->fc);
        }
    }
    free(pivs);
    pivs  = NULL;
    free(dr);
    dr  = NULL;

    mat->r  = realloc(mat->r, (unsigned long)npivs * sizeof(hl_t *));
    mat->np = mat->nr = mat->sz = npivs;
}

static cf32_t **sparse_AB_CD_linear_algebra_ff_32(
        mat_t *mat,
        const bs_t * bs,
        const stat_t * const st
        )
{
    len_t i, j;
    hl_t sc = 0;    /* starting column */

    const len_t nrows = mat->nr;
    const len_t ncols = mat->nc;
    const len_t nrl   = mat->nrl;
    const len_t nru   = mat->nru;
    const len_t ncl   = mat->ncl;

    hm_t **rows = mat->r;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)nru, sizeof(hm_t *));
    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = (hm_t **)malloc((unsigned long)nrl * sizeof(hm_t *));
    /* dense rows representing updated D part;
     * after reducing CD part with AB */
    cf32_t **drs    = (cf32_t **)calloc((unsigned long)nrl, sizeof(cf32_t *));

    i = 0;
    j = 1;
    for (i = 0; i < nrows; ++i) {
        if (!pivs[rows[i][3]]) {
            pivs[rows[i][3]]   = rows[i];
        } else {
            /* shorter rows first */
            upivs[nrl-j]  = rows[i];
            j++;
        }
    }
    free(mat->r);
    mat->r  = NULL;

    int64_t *dr  = (int64_t *)malloc(
            (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
    const int distribution = nrl/st->nthrds/100 == 0 ? 1 : nrl/st->nthrds/100;
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, sc) \
    schedule(static, distribution)
    for (i = 0; i < nrl; ++i) {
        cf32_t *cfs = NULL;
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        hm_t *npiv  = upivs[i];
        /* do the reduction */
        memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));
        cfs = bs->cf_32[npiv[0]];
        const len_t os  = npiv[1];
        const len_t len = npiv[2];
        const hm_t * const ds = npiv +3;
        for (j = 0; j < os; ++j) {
            drl[ds[j]]  = (int64_t)cfs[j];
        }
        for (; j < len; j += 4) {
            drl[ds[j]]    = (int64_t)cfs[j];
            drl[ds[j+1]]  = (int64_t)cfs[j+1];
            drl[ds[j+2]]  = (int64_t)cfs[j+2];
            drl[ds[j+3]]  = (int64_t)cfs[j+3];
        }
        sc  = ds[0];
        free(npiv);
        drs[i]  = reduce_dense_row_by_old_pivots_ff_32(
                drl, mat, bs, pivs, sc, st->fc);
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
    len_t npivs = 0; /* number of new pivots */
    for (i = 0; i < nrl; ++i) {
        if (drs[i] != NULL) {
            drs[npivs++]  = drs[i];
        }
    }
    if (npivs == 0) {
        free(drs);
        drs = NULL;
    }
    mat->np = npivs;

    return drs;
}

static cf32_t **interreduce_dense_matrix_ff_32(
    cf32_t **dm,
    const len_t ncr,
    const int32_t fc
    )
{
    len_t i, j, l;
    int64_t *dr = malloc((unsigned long)ncr * sizeof(int64_t));

    for (i = ncr-1; i > -1; --i) {
        if (dm[i]) {
            /* printf("interreduce dm[%d]\n", i); */
            memset(dr, 0, (unsigned long)ncr * sizeof(int64_t));
            const len_t npc  = ncr - i;
            const len_t os   = npc % 4;
            for (j = i, l = 0; l < os; ++j, ++l) {
                dr[j] = (int64_t)dm[i][l];
            }
            for (; l < npc; j += 4, l += 4) {
                dr[j]   = (int64_t)dm[i][l];
                dr[j+1] = (int64_t)dm[i][l+1];
                dr[j+2] = (int64_t)dm[i][l+2];
                dr[j+3] = (int64_t)dm[i][l+3];
            }
            free(dm[i]);
            dm[i] = NULL;
            /* start with previous pivot the reduction process, so keep the
             * pivot element as it is */
            dm[i] = reduce_dense_row_by_dense_new_pivots_ff_32(
                    dr, &i, dm, ncr, fc);
        }
    }
    free(dr);
    return dm;
}

static cf32_t **exact_dense_linear_algebra_ff_32(
        cf32_t **dm,
        mat_t *mat,
        const stat_t * const st
        )
{
    len_t i, j, k, l, npivs;

    const len_t nrows = mat->np; /* we need the pivots until now here */
    const len_t ncr   = mat->ncr;

    /* rows already representing new pivots */
    cf32_t **nps  = (cf32_t **)calloc((unsigned long)ncr, sizeof(cf32_t *));
    /* rows to be further reduced */
    cf32_t **tbr  = (cf32_t **)calloc((unsigned long)nrows, sizeof(cf32_t *));
    int64_t *dr   = (int64_t *)malloc(
            (unsigned long)(st->nthrds * ncr) * sizeof(int64_t));

    /* separate rows already representing new pivots and rows to
     * be further reduced by these new pivots */
    j     = 0;
    npivs = 0;
    for (i = 0; i < nrows; ++i) {
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
                memmove(dm[i], dm[i]+k, (unsigned long)(ncr-k) * sizeof(cf32_t));
                dm[i] = realloc(dm[i], (unsigned long)(ncr-k) * sizeof(cf32_t));
                nps[k] = dm[i];
                /* printf("nps[%d] = ", k); */
                if (nps[k][0] != 1) {
                    nps[k]  = normalize_dense_matrix_row_ff_32(nps[k], ncr-k, st->fc);
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
    tbr = realloc(tbr, (unsigned long)ntr * sizeof(cf32_t *));
    /* offset modulo 4 for loop unrolling, +1 due to storing the first
     * nonzero entry at the first position */

    /* reduction process to get all possible pivots, no interreduction here */
    const int distribution = ntr/st->nthrds/100 == 0 ? 1 : ntr/st->nthrds/100;
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l) shared(nps, tbr) \
    schedule(static, distribution)
    for (i = 0; i < ntr; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncr);
        memset(drl, 0, (unsigned long)ncr * sizeof(int64_t));
        hm_t npc  = 0;
        hm_t os   = 0;
        cf32_t *npiv  = tbr[i];
        os   = (ncr-npc) % 4;
        for (l = 0, j = npc; l < os; ++l, ++j) {
            drl[j]  = (int64_t)npiv[l];
        }
        for (; j < ncr; l += 4, j += 4) {
            drl[j]    = (int64_t)npiv[l];
            drl[j+1]  = (int64_t)npiv[l+1];
            drl[j+2]  = (int64_t)npiv[l+2];
            drl[j+3]  = (int64_t)npiv[l+3];
        }
        do {
            free(npiv);
            npiv = NULL;
            npiv = reduce_dense_row_by_dense_new_pivots_ff_32(
                    drl, &npc, nps, mat->ncr, st->fc);
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
    mat->np = npivs;

    free(tbr);
    free(dr);

    return nps;
}

static cf32_t **probabilistic_dense_linear_algebra_ff_32(
        cf32_t **dm,
        mat_t *mat,
        const stat_t * const st
        )
{
    len_t i, j, k, l, m, npivs;

    const int32_t fc  = st->fc;
    const len_t nrows = mat->np; /* we need the pivots until now here */
    const len_t ncols = mat->nc;
    const len_t ncr   = mat->ncr;

    /* rows already representing new pivots */
    cf32_t **nps  = (cf32_t **)calloc((unsigned long)ncr, sizeof(cf32_t *));
    /* rows to be further reduced */
    cf32_t **tbr  = (cf32_t **)calloc((unsigned long)nrows, sizeof(cf32_t *));

    /* separate rows already representing new pivots and rows to
     * be further reduced by these new pivots */
    j     = 0;
    npivs = 0;
    for (i = 0; i < nrows; ++i) {
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
                memmove(dm[i], dm[i]+k, (unsigned long)(ncr-k) * sizeof(cf32_t));
                dm[i] = realloc(dm[i], (unsigned long)(ncr-k) * sizeof(cf32_t));
                nps[k] = dm[i];
                /* printf("nps[%d] = ", k); */
                if (nps[k][0] != 1) {
                    nps[k]  = normalize_dense_matrix_row_ff_32(nps[k], ncr-k, st->fc);
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
    tbr = realloc(tbr, (unsigned long)ntr * sizeof(cf32_t *));

    const int64_t mod2  = (int64_t)st->fc * st->fc;

    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(ntr/3)))+1;
    const len_t rem = (ntr % nb == 0) ? 0 : 1;
    const len_t rpb = (ntr / nb) + rem;

    int64_t *dr   = (int64_t *)malloc(
        (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (unsigned long)(st->nthrds * rpb) * sizeof(int64_t));

    /* reduction process to get all possible pivots, no interreduction here */
    const int distribution = ntr/st->nthrds/100 == 0 ? 1 : ntr/st->nthrds/100;
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l) shared(nps, tbr) \
    schedule(static, distribution)
    for (i = 0; i < ntr; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncr);
        int64_t *mull = mul + (omp_get_thread_num() * rpb);
        const int32_t nbl   = (int32_t) (ntr > (i+1)*rpb ? (i+1)*rpb : ntr);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);

        if (nrbl != 0) {
            hm_t npc;
            hm_t os;
            cf32_t *tmp;
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
                memset(drl, 0, (unsigned long)ncr * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    tmp  = tbr[m];
                    for (l = 0, j = npc; l < os; ++l, ++j) {
                        drl[j]  -=  mull[k] * tmp[l];
                        drl[j]  +=  (drl[j] >> 63) & mod2;
                    }
                    for (; j < ncr; l += 4, j += 4) {
                        drl[j]    -=  mull[k] * tmp[l];
                        drl[j]    +=  (drl[j] >> 63) & mod2;
                        drl[j+1]  -=  mull[k] * tmp[l+1];
                        drl[j+1]  +=  (drl[j+1] >> 63) & mod2;
                        drl[j+2]  -=  mull[k] * tmp[l+2];
                        drl[j+2]  +=  (drl[j+2] >> 63) & mod2;
                        drl[j+3]  -=  mull[k] * tmp[l+3];
                        drl[j+3]  +=  (drl[j+3] >> 63) & mod2;
                    }
                }
                k   = 0;
                npc = 0;
                /* do the reduction */
                tmp = NULL;
                do {
                    free(tmp);
                    tmp = reduce_dense_row_by_dense_new_pivots_ff_32(
                            drl, &npc, nps, mat->ncr, st->fc);
                    if (npc == -1) {
                        bctr  = nrbl;
                        break;
                    }
                    k = __sync_bool_compare_and_swap(&nps[npc], NULL, tmp);
                    /* some other thread has already added a pivot so we have to
                    * recall the dense reduction process */
                } while (!k);
                bctr++;
            }
            for (j = i*rpb; j < nbl; ++j) {
                free(tbr[j]);
                tbr[j]  = NULL;
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
    mat->np = npivs;

    free(mul);
    free(tbr);
    free(dr);

    return nps;
}

static cf32_t **probabilistic_sparse_dense_echelon_form_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        const stat_t * const st
        )
{
    len_t i, j, k, l, m, npivs;

    const len_t nru   = mat->nru;
    const len_t nrl   = mat->nrl;
    const len_t nrows = mat->nr;
    const len_t ncr   = mat->ncr;
    const len_t ncols = mat->nc;

    hm_t **rows = mat->r;

    /* we fill in all known lead terms in pivs */
    hm_t **pivs   = (hm_t **)calloc((unsigned long)nru, sizeof(hm_t *));
    /* unkown pivot rows we have to reduce with the known pivots first */
    hm_t **upivs  = (hm_t **)malloc((unsigned long)nrl * sizeof(hm_t *));

    j = 1;
    for (i = 0; i < nrows; ++i) {
        if (!pivs[rows[i][3]]) {
            pivs[rows[i][3]]   = rows[i];
            /* printf("pivs %d is set by row %d\n", mat[i][3], i); */
        } else {
            /* shorter rows first */
            upivs[nrl-j]  = rows[i];
            j++;
        }
    }
    free(mat->r);
    mat->r  = NULL;

    /* rows already representing new pivots */
    cf32_t **nps  = (cf32_t **)calloc((unsigned long)ncr, sizeof(cf32_t *));

    const int32_t fc    = st->fc;
    const int64_t mod2  = (int64_t)fc * fc;

    /* compute rows per block */
    const len_t nb  = (len_t)(floor(sqrt(nrl/3)))+1;
    const len_t rem = (nrl % nb == 0) ? 0 : 1;
    const len_t rpb = (nrl / nb) + rem;

    int64_t *dr   = (int64_t *)malloc(
        (unsigned long)(st->nthrds * ncols) * sizeof(int64_t));
    int64_t *mul  = (int64_t *)malloc(
        (unsigned long)(st->nthrds * rpb) * sizeof(int64_t));

    /* reduction process to get all possible pivots, no interreduction here */
    const int distribution = nrl/st->nthrds/100 == 0 ? 1 : nrl/st->nthrds/100;
#pragma omp parallel for num_threads(st->nthrds) \
    private(i, j, k, l, m) shared(nps) \
    schedule(static, distribution)
    for (i = 0; i < nb; ++i) {
        int64_t *drl  = dr + (omp_get_thread_num() * ncols);
        int64_t *mull = mul + (omp_get_thread_num() * rpb);
        const int32_t nbl   = (int32_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
        const int32_t nrbl  = (int32_t) (nbl - i*rpb);
        if (nrbl != 0) {
            hm_t *npiv;
            cf32_t *tmp;
            hm_t npc;
            len_t bctr  = 0;
            while (bctr < nrbl) {
                npc = 0;

                /* fill random value array */
                for (j = 0; j < nrbl; ++j) {
                    mull[j] = (int64_t)rand() % fc;
                }
                /* generate one dense row as random linear combination
                 * of the rows of the block */
                memset(drl, 0, (unsigned long)ncols * sizeof(int64_t));

                for (k = 0, m = i*rpb; m < nbl; ++k, ++m) {
                    npiv  = upivs[m];
                    tmp   = bs->cf_32[npiv[0]];
                    const len_t os  = npiv[1];
                    const len_t len = npiv[2];
                    const hm_t * const ds = npiv + 3;
                    for (l = 0; l < os; ++l) {
                        drl[ds[l]]  -=  mull[k] * tmp[l];
                        drl[ds[l]]  +=  (drl[ds[l]] >> 63) & mod2;
                    }
                    for (; l < len; l += 4) {
                        drl[ds[l]]    -=  mull[k] * tmp[l];
                        drl[ds[l]]    +=  (drl[ds[l]] >> 63) & mod2;
                        drl[ds[l+1]]  -=  mull[k] * tmp[l+1];
                        drl[ds[l+1]]  +=  (drl[ds[l+1]] >> 63) & mod2;
                        drl[ds[l+2]]  -=  mull[k] * tmp[l+2];
                        drl[ds[l+2]]  +=  (drl[ds[l+2]] >> 63) & mod2;
                        drl[ds[l+3]]  -=  mull[k] * tmp[l+3];
                        drl[ds[l+3]]  +=  (drl[ds[l+3]] >> 63) & mod2;
                    }
                }
                k   = 0;
                npc = 0;
                /* do the reduction */
                tmp = NULL;
                do {
                    free(tmp);
                    tmp = reduce_dense_row_by_all_pivots_ff_32(
                            drl, mat, bs, &npc, pivs, nps, st->fc);
                    if (npc == -1) {
                        bctr  = nrbl;
                        break;
                    }
                    k = __sync_bool_compare_and_swap(&nps[npc], NULL, tmp);
                    /* some other thread has already added a pivot so we have to
                    * recall the dense reduction process */
                } while (!k);
                bctr++;
            }
            for (j = i*rpb; j < nbl; ++j) {
                free(upivs[j]);
                upivs[j]  = NULL;
            }
        }
    }
    npivs = 0;
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
    mat->np = npivs;


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

static void convert_to_sparse_matrix_rows_ff_32(
        mat_t *mat,
        cf32_t * const * const dm
        )
{
    if (mat->np == 0) {
        return;
    }

    len_t i, j, k, l;
    cf32_t *cfs;
    hm_t *dts, *dss;

    const len_t ncr = mat->ncr;
    const len_t ncl = mat->ncl;

    mat->r      = realloc(mat->r, (unsigned long)mat->np * sizeof(hm_t *));
    mat->cf_32  = realloc(mat->cf_32,
            (unsigned long)mat->np * sizeof(cf32_t *));

    l = 0;
    for (i = ncr-1; i > -1; --i) {
        if (dm[i] != NULL) {
            cfs = malloc((unsigned long)(ncr-i) * sizeof(cf32_t));
            dts = malloc((unsigned long)(ncr-i+3) * sizeof(hm_t));
            const hm_t len    = ncr-i;
            const hm_t os     = len % 4;
            const hm_t shift  = ncl+i;
            dss = dts + 3;

            for (k = 0, j = 0; j < os; ++j) {
                if (dm[i][j] != 0) {
                    cfs[k]    = dm[i][j];
                    dss[k++]  = j+shift;
                }
            }
            for (; j < len; j += 4) {
                if (dm[i][j] != 0) {
                    cfs[k]    = dm[i][j];
                    dss[k++]  = j+shift;
                }
                if (dm[i][j+1] != 0) {
                    cfs[k]    = dm[i][j+1];
                    dss[k++]  = j+1+shift;
                }
                if (dm[i][j+2] != 0) {
                    cfs[k]    = dm[i][j+2];
                    dss[k++]  = j+2+shift;
                }
                if (dm[i][j+3] != 0) {
                    cfs[k]    = dm[i][j+3];
                    dss[k++]  = j+3+shift;
                }
            }

            /* store meta data in first entries */
            dts[0]  = l; /* position of coefficient array in tmpcf */
            dts[1]  = k % 4;
            dts[2]  = k;

            /* adjust memory usage */
            dts = realloc(dts, (unsigned long)(k+3) * sizeof(hm_t));
            cfs = realloc(cfs, (unsigned long)k * sizeof(cf32_t));

            /* link to basis */
            mat->r[l]     = dts;
            mat->cf_32[l] = cfs;
            l++;
        }
    }
}

/* NOTE: this note is about the different linear algebra implementations:
 * exact and probabilistic linear algebra differ only in the last,
 * dense reduction step: the reduction of CD via AB is sparse and
 * the same for both. this generates a dense D' part which is then
 * either reduced via exact linear algebra or via probabilistic
 * linear algebra */
static void probabilistic_sparse_linear_algebra_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* allocate temporary storage space for sparse
     * coefficients of new pivot rows */
    mat->cf_32 = realloc(mat->cf_32,
            (unsigned long)mat->nrl * sizeof(cf32_t *));
    probabilistic_sparse_reduced_echelon_form_ff_32(mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%7d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void exact_sparse_linear_algebra_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* allocate temporary storage space for sparse
     * coefficients of new pivot rows */
    mat->cf_32  = realloc(mat->cf_32,
            (unsigned long)mat->nrl * sizeof(cf32_t *));
    exact_sparse_reduced_echelon_form_ff_32(mat, bs, st);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%7d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void exact_sparse_dense_linear_algebra_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    const len_t ncr = mat->ncr;

    /* generate updated dense D part via reduction of CD with AB */
    cf32_t **dm;
    dm  = sparse_AB_CD_linear_algebra_ff_32(mat, bs, st);
    if (mat->np > 0) {      
        dm  = exact_dense_linear_algebra_ff_32(dm, mat, st);
        dm  = interreduce_dense_matrix_ff_32(dm, ncr, st->fc);
    }

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    convert_to_sparse_matrix_rows_ff_32(mat, dm);

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
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%7d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void probabilistic_sparse_dense_linear_algebra_ff_32_2(
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    const len_t ncr = mat->ncr;

    /* generate updated dense D part via reduction of CD with AB */
    cf32_t **dm;
    dm  = sparse_AB_CD_linear_algebra_ff_32(mat, bs, st);
    if (mat->np > 0) {      
        dm  = probabilistic_dense_linear_algebra_ff_32(dm, mat, st);
        dm  = interreduce_dense_matrix_ff_32(dm, mat->ncr, st->fc);
    }

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    convert_to_sparse_matrix_rows_ff_32(mat, dm);

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
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%7d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}

static void probabilistic_sparse_dense_linear_algebra_ff_32(
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    const len_t ncr = mat->ncr;

    /* generate updated dense D part via reduction of CD with AB */
    cf32_t **dm = NULL;
    mat->np = 0;
    dm      = probabilistic_sparse_dense_echelon_form_ff_32(mat, bs, st);
    dm      = interreduce_dense_matrix_ff_32(dm, mat->ncr, st->fc);

    /* convert dense matrix back to sparse matrix representation,
     * use tmpcf for storing the coefficient arrays */
    convert_to_sparse_matrix_rows_ff_32(mat, dm);

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
    st->la_ctime  +=  ct1 - ct0;
    st->la_rtime  +=  rt1 - rt0;

    st->num_zerored += (mat->nrl - mat->np);
    if (st->info_level > 1) {
        printf("%7d new %7d zero", mat->np, mat->nrl - mat->np);
        fflush(stdout);
    }
}
