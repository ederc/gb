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

#include "convert.h"

/* after calling this procedure we have column indices instead of exponent
 * hashes in the polynomials resp. rows. moreover, we have sorted each row
 * by pivots / non-pivots. thus we get already an A|B splicing of the
 * initial matrix. this is a first step for receiving a full GBLA matrix. */
hl_t *convert_hashes_to_columns(
        mat_t *mat,
        ht_t *ht,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j, k, l, hi;
    hl_t *row;

    hd_t *hd  = ht->hd;

    int64_t nterms = 0;

    const len_t nr  = mat->nr;

    /* need to allocate memory for all possible exponents we
     * have in the local hash table since we do not which of
     * them are corresponding to multipliers and which are
     * corresponding to the multiplied terms in reducers. */
    hd_t **hcm = (hd_t **)malloc((unsigned long)(ht->eld-1) * sizeof(hd_t *));
    /* j counts all columns, k counts known pivots */
    for (j = 0, k = 0, i = 1; i < ht->eld; ++i) {
        hcm[j++]  = hd[i];
        if (hd[i].idx == 2) {
            k++;
        }
    }
    /* sort monomials w.r.t known pivots, then w.r.t. to the monomial order */
    qsort(hcm, (unsigned long)j, sizeof(hl_t), hcm_cmp);

    /* set number of rows and columns in ABCD splicing */
    mat->nru  = mat->ncl  = k;
    mat->nrl  = mat->nr - mat->nru;
    mat->ncr  = j - mat->ncl;

    st->num_rowsred     +=  mat->nrl;

    /* store the other direction (hash -> column) in HASH_IND */
    for (i = 0; i < j; ++i) {
        hd[hcm[i]].idx  = i;
    }



    /* map column positions to matrix rows */
#pragma omp parallel for num_threads(nthrds) private(i, j)
    for (i = 0; i < nr; ++i) {
        row = mat->r[i];
        for (j = 3; j < row[1]; ++j) {
            row[j]  = hd[row[j]].idx;
        }
        for (; j < row[2]; j += 4) {
            row[j]    = hd[row[j]].idx;
            row[j+1]  = hd[row[j+1]].idx;
            row[j+2]  = hd[row[j+2]].idx;
            row[j+3]  = hd[row[j+3]].idx;
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
    double density = (double)nterms / (double)mat->nr/ (double)mat->nc;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
    if (il > 1) {
        printf(" %7d x %-7d %8.3f%%", mat->nr, mat->nc, density);
        fflush(stdout);
    }

    return hcm;
}

void convert_sparse_matrix_rows_to_basis_elements(
        mat_t *mat,
        bs_t *bs,
        ht_t *ht,
        const hl_t *hcm,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j;

    hl_t *row = NULL;
    len_t bl  = bs->ld;

    const np  = mat->np;

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(bs, mat->np);


#pragma omp parallel for num_threads(nthrds) private(i, j)
    for (i = 0; i < np; ++i) {
        row = mat->r[i];
        for (j = 3; j < row[1]; ++j) {
            row[j] = hcm[row[j]];
        }
        for (; j < row[2]; j += 4) {
            row[j]   = hcm[row[j]];
            row[j+1] = hcm[row[j+1]];
            row[j+2] = hcm[row[j+2]];
            row[j+3] = hcm[row[j+3]];
        }
        bs->cf[bl+i]  = bs->tcf[row[0]];
        row[0]        = bl+i;
        bs->hd[bl+i]  = row;

        bs->lm[bl+i]  = ht->hd[row[3]].sdm;

        /* printf("new element [%d]\n", bl);
         * for (int32_t p = 3; p < gbcf[bl][2]; ++p) {
         *     printf("%d | ", gbcf[bl][p]);
         *     for (int32_t q = 0; q < nvars; ++q) {
         *         printf("%d", ev[gbdt[bl][p]][q]);
         *     }
         *     printf(" || ");
         * }
         * printf("\n"); */
    }

    /* printf("thread %d frees tmpcf %p\n", omp_get_thread_num(), tmpcf); */
    free(bs->tcf);
    bs->tcf = NULL;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}
