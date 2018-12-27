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
        dt_t **matdt,
        stat_t *st
        )
{
    len_t i, j, k, l;
    dt_t *row;
    hl_t *hcm; /* hash-to-column map */
    int64_t nterms = 0;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t hi;

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
    /* for (int ii=0; ii < j; ++ii) {
     *     printf("hcm[%d] = %d | %d | ", ii, hcm[ii], hds[hcm[ii]].idx);
     *     for (int jj=0; jj < nvars; ++jj) {
     *         printf("%d ", evs[hcm[ii]][jj]);
     *     }
     *     printf("\n");
     * } */
    qsort(hcm, (unsigned long)j, sizeof(hl_t), hcm_cmp);
    /* for (int ii=0; ii < j; ++ii) {
     *     printf("sorted hcm[%d] = %d | ", ii, hcm[ii]);
     *     for (int jj=0; jj < nvars; ++jj) {
     *         printf("%d ", evs[hcm[ii]][jj]);
     *     }
     *     printf("\n");
     * } */

    nru = ncl = k;
    nrl = nrows - nru;
    ncr = j - ncl;

    st->num_rowsred     +=  nrl;

    /* store the other direction (hash -> column) in HASH_IND */
    for (k = 0; k < j; ++k) {
        hds[hcm[k]].idx  = k;
    }
#if 0
    /* j counts all columns, k counts known pivots */
    for (j = 0, k = 0, i = 0; i < nrows; ++i) {
        row     =   matdt[i];
        nterms  +=  (int64_t)row[2];
        for (l = 3; l < row[2]; ++l) {
            /* printf("row[%d][%d] = %d | ", i, l, row[l]);
             * for (int32_t p = 0; p < nvars; ++p) {
             *     printf("%d ", ev[row[l]][p]);
             * }
             * printf("\n"); */
            hi  = hd[row[l]].idx;
            /* printf("hi %d\n", hi); */
#if ORDER_COLUMNS
            if (hi > 0) {
#else
                if (hi != 0) {
#endif
                    /* printf("j %d / %d ncols\n", j, ncols); */
                    hcm[j++]  = row[l];
                    if (hi == 2) {
                        k++;
#if ORDER_COLUMNS
                        hd[row[l]].idx  = -1;
                    } else {
                        hd[row[l]].idx  = -2;
                    }
#else
                }
                hd[row[l]].idx  = 0;
#endif
            }
        }
    }
    /* for (i = 0; i < j; ++i) {
     *   printf("hcm[%d] = %d | %d || ", i, hcm[i], (ev+hcm[i])[HASH_DEG]);
     *   for (l = 0; l < nvars; ++l) {
     *     printf("%d ", (ev+hcm[i])[l]);
     *   }
     *   printf("\n");
     * } */
    /* sort monomials w.r.t known pivots, then w.r.t. to the monomial order */
    qsort(hcm, (unsigned long)j, sizeof(hl_t), hcm_cmp);
    /* printf("ncl %d\n", k);
     * for (i = 0; i < j; ++i) {
     *   printf("hcm[%d] = ", i);
     *   for (l = 0; l < nvars; ++l) {
     *     printf("%d ", ev[hcm[i]][l]);
     *   }
     *   printf("\n");
     * } */

    /* set number of rows and columns in ABCD splicing */
    nru = ncl = k;
    nrl = nrows - nru;
    ncr = j - ncl;

    st->num_rowsred     +=  nrl;

    /* store the other direction (hash -> column) in HASH_IND */
    for (i = 0; i < j; ++i) {
        hd[hcm[i]].idx  = i;
    }
#endif


    /* map column positions to matrix rows */
#pragma omp parallel for num_threads(nthrds) private(i, j)
    for (i = 0; i < nrows; ++i) {
        row = matdt[i];
        nterms  +=  row[2];
        for (j = 3; j < row[1]; ++j) {
            row[j]  = hds[row[j]].idx;
        }
        for (; j < row[2]; j += 4) {
            row[j]    = hds[row[j]].idx;
            row[j+1]  = hds[row[j+1]].idx;
            row[j+2]  = hds[row[j+2]].idx;
            row[j+3]  = hds[row[j+3]].idx;
        }
        /* for (j = 3; j < row[2]; ++j) {
         *     printf("%d ", row[j]);
         * }
         * printf("\n"); */
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
    double density = (double)nterms / (double)nrows / (double)ncols;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
    if (il > 1) {
        printf(" %7d x %-7d %8.3f%%", nrows, ncols, density);
        fflush(stdout);
    }

    return hcm;
}

static void convert_sparse_matrix_rows_to_basis_elements(
        dt_t **mat,
        const hl_t *hcm,
        stat_t *st
        )
{
    len_t i, j;

    len_t bl  = bload;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(npivs);

/* #pragma omp parallel for num_threads(nthrds) private(i, j) */
    for (i = 0; i < npivs; ++i) {
        for (j = 3; j < mat[i][1]; ++j) {
            mat[i][j] = insert_in_global_hash_table(evs[hcm[mat[i][j]]]);
        }
        for (; j < mat[i][2]; j += 4) {
            mat[i][j]   = insert_in_global_hash_table(evs[hcm[mat[i][j]]]);
            mat[i][j+1] = insert_in_global_hash_table(evs[hcm[mat[i][j+1]]]);
            mat[i][j+2] = insert_in_global_hash_table(evs[hcm[mat[i][j+2]]]);
            mat[i][j+3] = insert_in_global_hash_table(evs[hcm[mat[i][j+3]]]);
        }
        gbcf[bl+i]  = tmpcf[mat[i][0]];
        mat[i][0]   = bl+i;
        gbdt[bl+i]  = mat[i];
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
    free(tmpcf);
    tmpcf = NULL;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}
