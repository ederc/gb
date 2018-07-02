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
        dt_t **matdt
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
    hcm = (hl_t *)malloc((unsigned long)ncols * sizeof(hl_t));
    /* j counts all columns, k counts known pivots */
    for (j = 0, k = 0, i = 0; i < nrows; ++i) {
        row     =   matdt[i];
        nterms  +=  (int64_t)row[2];
        for (l = 3; l < row[2]; ++l) {
            hi  = hd[row[l]].idx;
#if ORDER_COLUMNS
            if (hi > 0) {
#else
                if (hi != 0) {
#endif
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
    /* for (i = 0; i < j; ++i) {
     *   printf("hcm[%d] = ", i);
     *   for (l = 0; l < nvars; ++l) {
     *     printf("%d ", (ev+hcm[i]*hl)[l]);
     *   }
     *   printf("\n");
     * } */

    /* set number of rows and columns in ABCD splicing */
    nru = ncl = k;
    nrl = nrows - nru;
    ncr = j - ncl;

    /* store the other direction (hash -> column) in HASH_IND */
    for (i = 0; i < j; ++i) {
        hd[hcm[i]].idx  = i;
    }



    /* map column positions to matrix rows */
#pragma omp parallel for num_threads(nthrds) private(i, j)
    for (i = 0; i < nrows; ++i) {
        row = matdt[i];
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
    density = (double)nterms / (double)nrows / (double)ncols;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    convert_ctime +=  ct1 - ct0;
    convert_rtime +=  rt1 - rt0;
    GB_DEBUG(SYMDBG, " %7d x %7d mat - %6.3f%%", nrows, ncols, density);

    return hcm;
}

static void convert_dense_matrix_to_basis_elements(
        cf_t * const *dm,
        const dt_t *hcm
        )
{
    len_t i, j, k;
    cf_t *cfs, *dr;
    dt_t *dts;

    len_t bl  = bload;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* fix size of basis for entering new elements directly */
    check_enlarge_basis(npivs);

#pragma omp parallel for num_threads(nthrds) private(i, j)
    for (i = 0; i < ncr; ++i) {
        if (dm[i] != NULL) {
            dr  = dm[i];
            cfs = malloc((unsigned long)(ncr-i+3) * sizeof(cf_t));
            dts = malloc((unsigned long)(ncr-i+3) * sizeof(dt_t));
            const dt_t os = (ncr-i) % 4;

            for (k = 3, j = 0; j < os; ++j) {
                if (dr[j] != 0) {
                    cfs[k]    = dr[j];
                    dts[k++]  = hcm[j+ncl];
                }
            }
            /* store meta data in first entries */
            dts[0]  = bl;
            dts[1]  = k + 3;
            dts[2]  = (k % 4) + 3;
            cfs[0]  = 0;
            cfs[1]  = k + 3;
            cfs[2]  = (k % 4) + 3;

            /* adjust memory usage */
            dts = realloc(dts, (unsigned long)(k+3) * sizeof(dt_t));
            cfs = realloc(cfs, (unsigned long)(k+3) * sizeof(cf_t));

            /* link to basis */
            gbdt[bl]  = dts;
            gbcf[bl]  = cfs;
            bl++;
        }
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    convert_ctime +=  ct1 - ct0;
    convert_rtime +=  rt1 - rt0;
}
