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
static hl_t *convert_hashes_to_columns(
        hl_t *hcm,
        mat_t *mat,
        stat_t *st,
        ht_t *sht
        )
{
    len_t i, j, k;
    hm_t *row;
    int64_t nterms = 0;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t hi;

    const hl_t esld = sht->eld;
    hd_t *hds       = sht->hd;
    hm_t **rows     = mat->r;

    /* need to allocate memory for all possible exponents we
     * have in the local hash table since we do not which of
     * them are corresponding to multipliers and which are
     * corresponding to the multiplied terms in reducers. */
    hcm = (hl_t *)malloc((unsigned long)(esld-1) * sizeof(hl_t));
    for (k = 0, j = 0, i = 1; i < esld; ++i) {
        hi  = hds[i].idx;

        hcm[j++]  = i;
        if (hi == 2) {
            k++;
        }
    }
    sort_r(hcm, (unsigned long)j, sizeof(hl_t), hcm_cmp, sht);

    mat->nru  = mat->ncl = k;
    mat->nrl  = mat->nr - mat->nru;
    mat->ncr  = j - mat->ncl;

    st->num_rowsred     +=  nrl;

    const len_t hld = j;
    /* store the other direction (hash -> column) */
    for (k = 0; k < hld; ++k) {
        hds[hcm[k]].idx  = k;
    }

    const len_t mnr = mat->nr;
    /* map column positions to matrix rows */
#pragma omp parallel for num_threads(nthrds) private(i, j)
    for (i = 0; i < mnr; ++i) {
        const len_t os  = rows[i][1];
        const len_t len = rows[i][2];
        row = rows[i] + 3;
        nterms  +=  len;
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
    double density = (double)nterms / (double)mat->nr / (double)mat->nc;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
    if (st->info_level > 1) {
        printf(" %7d x %-7d %8.3f%%", mat->nr, mat->nc, density);
        fflush(stdout);
    }
}

static void convert_sparse_matrix_rows_to_basis_elements(
        hm_t **mat,
        const hl_t *hcm,
        stat_t *st
        )
{
    len_t i, ctr;

    len_t bl  = bload;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(npivs);

/* #pragma omp parallel for num_threads(nthrds) private(i, j) */
    ctr = 0;
    for (i = 0; i < npivs; ++i) {
        ctr +=  mat[i][2];
    }
    while (esz - eld < ctr) {
        enlarge_basis_hash_table();
    }
    for (i = 0; i < npivs; ++i) {
        insert_in_basis_hash_table_pivots(mat[i], hcm);
        gbcf_ff[bl+i] = tmpcf_ff[mat[i][0]];
        mat[i][0]     = bl+i;
        gbdt[bl+i]    = mat[i];
    }

    free(tmpcf_ff);
    tmpcf_ff = NULL;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}
