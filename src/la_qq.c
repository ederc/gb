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
 * \file la_qq.c
 * \brief Implementation of rational linear algebra.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "data.h"

static inline mpz_t *remove_content_of_sparse_matrix_row_qq(
        mpz_t *row,
        const len_t os,
        const len_t len
        )
{
    len_t i;
    long remove = 1;

    mpz_t content;
    mpz_init(content);
    /* compute content, i.e. gcd of all coefficients */
    mpz_set(content, row[0]);
    /* printf("before removing content: ");
     * for (i = 0; i < len; ++i) {
     *     gmp_printf("%Zd ", row[i]);
     * }
     * printf("\n"); */
    for (i = 1; i < len; ++i) {
        /* gmp_printf("gcd of %Zd and %Zd is ", content, row[i]); */
        mpz_gcd(content, content, row[i]);
        /* gmp_printf("-> %Zd\n", content); */
        if (mpz_cmp_si(content, 1) == 0) {
            remove = 0;
            break;
        }
    }
    if (remove == 1) {
        /* remove content */
        for (i = 0; i < os; ++i) {
            mpz_divexact(row[i], row[i], content);
        }
        for (; i < len; i += 4) {
            mpz_divexact(row[i], row[i], content);
            mpz_divexact(row[i+1], row[i+1], content);
            mpz_divexact(row[i+2], row[i+2], content);
            mpz_divexact(row[i+3], row[i+3], content);
        }
    }
    mpz_clear(content);

    /* make lead coefficient positive */
    if (mpz_sgn(row[0]) < 0) {
        for (i = 0; i < os; ++i) {
            mpz_neg(row[i], row[i]);
        }
        for (; i < len; i += 4) {
            mpz_neg(row[i], row[i]);
            mpz_neg(row[i+1], row[i+1]);
            mpz_neg(row[i+2], row[i+2]);
            mpz_neg(row[i+3], row[i+3]);
        }
    }

    /* printf("after removing content: ");
     * for (i = 0; i < len; ++i) {
     *     gmp_printf("%Zd ", row[i]);
     * }
     * printf("\n"); */
    return row;
}

static hm_t *reduce_dense_row_by_known_pivots_sparse_qq(
        mpz_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hl_t dpiv,    /* pivot of dense row at the beginning */
        const hm_t tmp_pos  /* position of new coeffs array in tmpcf */
        )
{
    hl_t i, j;
    hm_t *dts;
    mpz_t *cfs;
    len_t np  = -1;
    const len_t ncols         = mat->nc;
    const len_t ncl           = mat->ncl;
    mpz_t * const * const mcf = mat->cf_qq;

    hm_t *row = NULL;
    mpz_t *cf = NULL;
    len_t rlen  = 0;

    mpz_t mul1, mul2;
    mpz_inits(mul1, mul2, NULL);
    for (i = dpiv; i < ncols; ++i) {
        /* uses mpz_sgn for checking if dr[i] = 0 */
        if (mpz_sgn(dr[i]) == 0) {
            continue;
        }
        if (pivs[i] == NULL) {
            if (np == -1) {
                row = (hm_t *)malloc(
                        (unsigned long)(ncols-np+3) * sizeof(hm_t));
                cf  = (mpz_t *)malloc(
                        (unsigned long)(ncols-np) * sizeof(mpz_t));
                np  = i;
            }
            mpz_init(cf[rlen]);
            mpz_swap(cf[rlen], dr[i]);
            row[rlen+3] = i;
            rlen++;
            continue;
        }
        /* found reducer row, get multiplier */
        dts = pivs[i];
        if (i < ncl) {
            cfs   = bs->cf_qq[dts[0]];
        } else {
            cfs   = mcf[dts[0]];
        }
        const len_t os  = dts[1];
        const len_t len = dts[2];
        const hm_t * const ds  = dts + 3;

        /* check if lead coefficient of dr is multiple of lead coefficient
         * of cfs, generate corresponding multipliers respectively */
        if (mpz_divisible_p(dr[i], cfs[0]) != 0) {
            mpz_divexact(mul2, dr[i], cfs[0]);
        } else {
            mpz_lcm(mul1, dr[i], cfs[0]);
            mpz_divexact(mul2, mul1, cfs[0]);
            mpz_divexact(mul1, mul1, dr[i]);
            for (j = 0; j < rlen; ++j) {
                mpz_mul(cf[j], cf[j], mul1);
            }
            for (j = i+1; j < ncols; ++j) {
                if (mpz_sgn(dr[j]) != 0) {
                    mpz_mul(dr[j], dr[j], mul1);
                }
            }
        }
        for (j = 0; j < os; ++j) {
            mpz_submul(dr[ds[j]], mul2, cfs[j]);
        }
        for (; j < len; j += 4) {
            mpz_submul(dr[ds[j]], mul2, cfs[j]);
            mpz_submul(dr[ds[j+1]], mul2, cfs[j+1]);
            mpz_submul(dr[ds[j+2]], mul2, cfs[j+2]);
            mpz_submul(dr[ds[j+3]], mul2, cfs[j+3]);
        }
    }
    if (rlen != 0) {
        row     = realloc(row, (unsigned long)(rlen+3) * sizeof(hm_t));
        cf      = realloc(cf, (unsigned long)rlen * sizeof(mpz_t));
        row[0]  = tmp_pos;
        row[1]  = rlen % 4;
        row[2]  = rlen;
        mat->cf_qq[tmp_pos]  = cf;
    }
    mpz_clears(mul1, mul2, NULL);
    return row;
}

static void exact_sparse_reduced_echelon_form_qq(
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

    const len_t drlen = st->nthrds * ncols;
    mpz_t *dr  = (mpz_t *)malloc(
            (unsigned long)drlen * sizeof(mpz_t));
    for (i = 0; i < drlen; ++i) {
        mpz_init(dr[i]);
    }
    /* mo need to have any sharing dependencies on parallel computation,
     * no data to be synchronized at this step of the linear algebra */
#pragma omp parallel for num_threads(st->nthrds) private(i, j, k, sc)
    for (i = 0; i < nrl; ++i) {
        mpz_t *drl  = dr + (omp_get_thread_num() * ncols);
        hm_t *npiv  = upivs[i];
        mpz_t *cfs  = bs->cf_qq[npiv[0]];
        len_t os    = npiv[1];
        len_t len   = npiv[2];
        hm_t * ds   = npiv + 3;
        k = 0;
        /* reset entries to zero */
        for (j = 0; j < ncols; ++j) {
            mpz_set_si(drl[j], 0);
        }
        for (j = 0; j < os; ++j) {
            mpz_set(drl[ds[j]], cfs[j]);
        }
        for (; j < len; j += 4) {
            mpz_set(drl[ds[j]], cfs[j]);
            mpz_set(drl[ds[j+1]], cfs[j+1]);
            mpz_set(drl[ds[j+2]], cfs[j+2]);
            mpz_set(drl[ds[j+3]], cfs[j+3]);
        }
        cfs = NULL;
        k   = 1;
        do {
            os   = npiv[1];
            len  = npiv[2];
            sc  = npiv[3];
            ds  = npiv + 3;
            if (k == 0) {
                /* if we redo the row than we only handle elements
                 * with column index >= sc, so we do not need to
                 * reset the other entries in the row to zero */
                for (j = sc; j < ncols; ++j) {
                    mpz_set_si(drl[j], 0);
                }
                for (j = 0; j < os; ++j) {
                    mpz_set(drl[ds[j]], cfs[j]);
                }
                for (; j < len; j += 4) {
                    mpz_set(drl[ds[j]], cfs[j]);
                    mpz_set(drl[ds[j+1]], cfs[j+1]);
                    mpz_set(drl[ds[j+2]], cfs[j+2]);
                    mpz_set(drl[ds[j+3]], cfs[j+3]);
                }
                for (j = 0; j < len; ++j) {
                    mpz_clear(cfs[j]);
                }
            }
            free(cfs);
            free(npiv);
            npiv  = reduce_dense_row_by_known_pivots_sparse_qq(
                    drl, mat, bs, pivs, sc, i);
            if (!npiv) {
                break;
            }
            /* remove content of coefficient array for better usage as reducer
             * later on.
             * NOTE: this has to be done here, otherwise the reduction may
             * lead to wrong results in a parallel computation since other
             * threads might directly use the new pivot once it is synced. */
            if (mpz_cmp_si(mat->cf_qq[npiv[0]][0], 1) != 0) {
                remove_content_of_sparse_matrix_row_qq(
                        mat->cf_qq[npiv[0]], npiv[1], npiv[2]);
            }
            k   = __sync_bool_compare_and_swap(&pivs[npiv[3]], NULL, npiv);
            cfs = mat->cf_qq[npiv[0]];
        } while (k == 0);
        cfs = NULL;
    }
    free(upivs);
    upivs = NULL;

    /* we do not need the old pivots anymore */
    for (i = 0; i < ncl; ++i) {
        free(pivs[i]);
        pivs[i] = NULL;
    }

    len_t npivs = 0; /* number of new pivots */

    for (i = ncols; i < st->nthrds*ncols; ++i) {
        mpz_clear(dr[i]);
    }
    dr      = realloc(dr, (unsigned long)ncols * sizeof(mpz_t));
    /* mpz_t *dr  = (mpz_t *)malloc(
     *         (unsigned long)(ncols) * sizeof(mpz_t)); */
    /* for (i = 0; i < ncols; ++i) {
     *     mpz_set_si(dr[i], 0);
     * } */
    mat->r  = realloc(mat->r, (unsigned long)ncr * sizeof(hm_t *));
    rows    = mat->r;

    /* interreduce new pivots */
    mpz_t *cfs;
    hm_t cf_array_pos;
    for (i = (ncols-1); i >= nru; --i) {
        if (pivs[i]) {
            /* reset entries to zero */
    /*         rows[npivs++] = pivs[i];
     *     }
     * } */
            for (j = 0; j < ncols; ++j) {
                mpz_set_si(dr[j], 0);
            }
            cfs = mat->cf_qq[pivs[i][0]];
            cf_array_pos    = pivs[i][0];
            const len_t os  = pivs[i][1];
            const len_t len = pivs[i][2];
            const hm_t * const ds = pivs[i] + 3;
            sc  = ds[0];
            for (j = 0; j < os; ++j) {
                mpz_set(dr[ds[j]], cfs[j]);
            }
            for (; j < len; j += 4) {
                mpz_set(dr[ds[j]], cfs[j]);
                mpz_set(dr[ds[j+1]], cfs[j+1]);
                mpz_set(dr[ds[j+2]], cfs[j+2]);
                mpz_set(dr[ds[j+3]], cfs[j+3]);
            }
            free(pivs[i]);
            free(cfs);
            pivs[i] = NULL;
            pivs[i] = rows[npivs] =
                reduce_dense_row_by_known_pivots_sparse_qq(
                        dr, mat, bs, pivs, sc, cf_array_pos);
            remove_content_of_sparse_matrix_row_qq(
                    mat->cf_qq[rows[npivs][0]], rows[npivs][1], rows[npivs][2]);
            npivs++;
        }
    }
    free(pivs);
    pivs  = NULL;
    for (i = 0; i < ncols; ++i) {
        mpz_clear(dr[i]);
    }
    free(dr);
    dr  = NULL;

    mat->r  = realloc(mat->r, (unsigned long)npivs * sizeof(hl_t *));
    mat->np = mat->nr = mat->sz = npivs;
}

static void exact_sparse_linear_algebra_qq(
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
    mat->cf_qq  = realloc(mat->cf_qq,
            (unsigned long)mat->nrl * sizeof(mpz_t *));
    exact_sparse_reduced_echelon_form_qq(mat, bs, st);

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
