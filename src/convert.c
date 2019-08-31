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
 * \file convert.c
 * \brief Implementation of hash value / column index conversions
 * for the exponent hashes of polynomials resp. rows.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

/* after calling this procedure we have column indices instead of exponent
 * hashes in the polynomials resp. rows. moreover, we have sorted each row
 * by pivots / non-pivots. thus we get already an A|B splicing of the
 * initial matrix. this is a first step for receiving a full GBLA matrix. */
static void convert_hashes_to_columns(
        hl_t **hcmp,
        mat_t *mat,
        stat_t *st,
        ht_t *sht
        )
{
    len_t i, j, k;
    hm_t *row;
    int64_t nterms = 0;

    hl_t *hcm = *hcmp;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t hi;

    const len_t mnr = mat->nr;
    const hl_t esld = sht->eld;
    hd_t *hds       = sht->hd;
    hm_t **rows     = mat->r;

    /* all elements in the sht hash table represent
     * exactly one column of the matrix */
    hcm = realloc(hcm, (unsigned long)(esld-1) * sizeof(hl_t));
    for (k = 0, j = 0, i = 1; i < esld; ++i) {
        hi  = hds[i].idx;

        hcm[j++]  = i;
        if (hi == 2) {
            k++;
        }
    }
    sort_r(hcm, (unsigned long)j, sizeof(hl_t), hcm_cmp, sht);

    /* printf("hcm\n");
     * for (int ii=0; ii<j; ++ii) {
     *     printf("hcm[%d] = %d | ", ii, hcm[ii]);
     *     for (int jj = 0; jj < sht->nv; ++jj) {
     *         printf("%d ", sht->ev[hcm[ii]][jj]);
     *     }
     *     printf("\n");
     * } */

    mat->nru  = mat->ncl = k;
    mat->nrl  = mat->nr - mat->nru;
    mat->ncr  = esld - 1 - mat->ncl;

    st->num_rowsred +=  mat->nrl;

    /* store the other direction (hash -> column) */
    for (k = 0; k < esld-1; ++k) {
        hds[hcm[k]].idx  = k;
    }

    /* map column positions to matrix rows */
#pragma omp parallel for num_threads(st->nthrds) private(i, j)
    for (i = 0; i < mnr; ++i) {
        const len_t os  = rows[i][1];
        const len_t len = rows[i][2];
        row = rows[i] + 3;
        for (j = 0; j < os; ++j) {
            row[j]  = hds[row[j]].idx;
        }
        for (; j < len; j += 4) {
            row[j]    = hds[row[j]].idx;
            row[j+1]  = hds[row[j+1]].idx;
            row[j+2]  = hds[row[j+2]].idx;
            row[j+3]  = hds[row[j+3]].idx;
        }
    }
    for (i = 0; i < mnr; ++i) {
        nterms  +=  rows[i][2];
    }

    /* next we sort each row by the new colum order due
     * to known / unkown pivots */

    /* NOTE: As strange as it may sound, we do not need to sort the rows.
     * When reducing, we copy them to dense rows, there we copy the coefficients
     * at the right place and reduce then. For the reducers itself it is not
     * important in which order the terms are represented as long as the first
     * term is the lead term, which is always true. Once a row is finally reduced
     * it is copied back to a sparse representation, now in the correct term
     * order since it is coming from the correctly sorted dense row. So all newly
     * added elements have all their terms sorted correctly w.r.t. the given
     * monomial order. */

    /* compute density of matrix */
    nterms  *=  100; /* for percentage */
    double density = (double)nterms / (double)mnr / (double)mat->nc;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
    if (st->info_level > 1) {
        printf(" %7d x %-7d %8.3f%%", mat->nr, mat->nc, density);
        fflush(stdout);
    }
    *hcmp = hcm;
}

static void convert_sparse_matrix_rows_to_basis_elements(
        mat_t *mat,
        bs_t *bs,
        ht_t *bht,
        const ht_t * const sht,
        const hl_t * const hcm,
        stat_t *st
        )
{
    len_t i, ctr;

    const len_t bl  = bs->ld;
    const len_t np  = mat->np;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(bs, mat->np);

    hm_t **rows = mat->r;

    ctr = 0;
    for (i = 0; i < np; ++i) {
        ctr +=  rows[i][2];
    }
    while (bht->esz - bht->eld < ctr) {
        enlarge_hash_table(bht);
    }
    if (st->fc > 0) {
        for (i = 0; i < np; ++i) {
            insert_in_basis_hash_table_pivots(rows[i], bht, sht, hcm);
            bs->cf_ff[bl+i] = mat->cf_ff[rows[i][0]];
            rows[i][0]      = bl+i;
            bs->hm[bl+i]    = rows[i];
            /* printf("\n");
             * for (int jj = 0; jj < bht->nv; ++jj) {
             *     printf("%d ", bht->ev[bs->hm[bl+i][3]][jj]);
             * }
             * printf(" + %d terms\n", rows[i][2]-1); */
            /* for (int ii=0; ii < rows[i][2]; ++ii) {
             *     printf("%d | ", bs->cf_ff[bl+i][ii]);
             *     for (int jj = 0; jj < bht->nv; ++jj) {
             *         printf("%d ", bht->ev[bs->hm[bl+i][ii]][jj]);
             *     }
             *     printf(" || ");
             * }
             * printf("\n"); */
        }
    } else {
        if (st->fc == 0) {
            for (i = 0; i < np; ++i) {
                insert_in_basis_hash_table_pivots(rows[i], bht, sht, hcm);
                bs->cf_qq[bl+i] = mat->cf_qq[rows[i][0]];
                rows[i][0]      = bl+i;
                bs->hm[bl+i]    = rows[i];
                /* printf("\n bl+i = %d = %d + %d\n", bl+i, bl, i);
                 * for (int jj = 0; jj < bht->nv; ++jj) {
                 *     printf("%d ", bht->ev[bs->hm[bl+i][3]][jj]);
                 * }
                 * printf(" + %d terms\n", rows[i][2]-1);
                 * for (int ii=0; ii < rows[i][2]; ++ii) {
                 *     gmp_printf("%Zd | ", bs->cf_qq[bl+i][ii]);
                 *     for (int jj = 0; jj < bht->nv; ++jj) {
                 *         printf("%d ", bht->ev[bs->hm[bl+i][ii+3]][jj]);
                 *     }
                 *     printf(" || ");
                 * }
                 * printf("\n"); */
            }
        }
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}
