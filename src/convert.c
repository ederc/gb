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
hd_t **convert_hashes_to_columns(
        mat_t *mat,
        ht_t *ht,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j, k;

    hd_t *hd  = ht->hd;
    hd_t **h  = NULL;

    int64_t nterms = 0;

    const len_t nru = mat->nru;
    const len_t nrl = mat->nrl;
    const int32_t nthrds  = st->nthrds;

    mon_t p;
    ci_t *c;

    printf("nru %d || nrl %d\n", nru, nrl);
    mat->pv   = realloc(mat->pv, (unsigned long)nru * sizeof(row_t));
    mat->tbr  = realloc(mat->tbr, (unsigned long)nrl * sizeof(row_t));
    /* need to allocate memory for all possible exponents we
     * have in the local hash table since we do not which of
     * them are corresponding to multipliers and which are
     * corresponding to the multiplied terms in reducers. */
    hd_t **hcm = (hd_t **)malloc((unsigned long)(ht->eld-1) * sizeof(hd_t *));
    /* j counts all columns, k counts known pivots */
    for (j = 0, k = 0, i = 1; i < ht->eld; ++i) {
        hcm[j++]  = hd+i;
        printf("hcm[%d] = %p | [%d]idx %d\n", j-1, hcm[j-1], i, hd[i].idx);
        if (hd[i].idx == 2) {
            k++;
        }
    }
    printf("j %d | ht->eld %d\n", j, ht->eld);
    /* sort monomials w.r.t known pivots, then w.r.t. to the monomial order */
    qsort(hcm, (unsigned long)j, sizeof(hd_t *), hcm_cmp);

    /* set number of rows and columns in ABCD splicing */
    mat->ncl  = k;
    mat->ncr  = j - mat->ncl;

    st->num_rowsred +=  mat->nrl;

    /* store the other direction (hash -> column) in HASH_IND */
#pragma omp parallel for num_threads(nthrds) private(i)
    for (i = 0; i < j; ++i) {
        hcm[i]->idx  = i;
    }

    /* map column positions to matrix rows */
#pragma omp parallel for num_threads(nthrds) private(i, j)
    for (i = 0; i < nru; ++i) {
        p = mat->pmp[i];
        mat->pv[i].sz  = p.sz;
        mat->pv[i].of  = p.of;
        mat->pv[i].cl  = p.cl;
        c = (ci_t *)malloc((unsigned long)(p.sz) * sizeof(ci_t));
        h = p.h;
        for (j = 0; j < p.of; ++j) {
            c[j]  = h[j]->idx;
        }
        for (; j < p.sz; j += 4) {
            c[j]    = h[j]->idx;
            c[j+1]  = h[j+1]->idx;
            c[j+2]  = h[j+2]->idx;
            c[j+3]  = h[j+3]->idx;
        }
        mat->pv[i].ci  = c;
        free(mat->pmp[i].h);
        mat->pmp[i].h  = NULL;
    }
    free(mat->pmp);
    mat->pmp = NULL;

#pragma omp parallel for num_threads(nthrds) private(i, j)
    for (i = 0; i < nrl; ++i) {
        p = mat->npmp[i];
        mat->tbr[i].sz = p.sz;
        mat->tbr[i].of = p.of;
        mat->tbr[i].cl = p.cl;
        c = (ci_t *)malloc((unsigned long)(p.sz) * sizeof(ci_t));
        h = p.h;
        for (j = 0; j < p.of; ++j) {
            c[j]  = h[j]->idx;
        }
        for (; j < p.sz; j += 4) {
            c[j]    = h[j]->idx;
            c[j+1]  = h[j+1]->idx;
            c[j+2]  = h[j+2]->idx;
            c[j+3]  = h[j+3]->idx;
        }
        mat->tbr[i].ci = c;
        free(mat->npmp[i].h);
        mat->npmp[i].h  = NULL;
    }
    free(mat->npmp);
    mat->npmp = NULL;

    printf("mat->ncr %d\n", mat->ncr);
    /* allocate space for coefficient arrays of possible new pivots */
    printf("mat->npv before %p\n", mat->npv);
    mat->npv  = realloc(mat->npv, (unsigned long)mat->ncr * sizeof(row_t));
    memset(mat->npv, 0, (unsigned long)mat->ncr * sizeof(row_t));

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
    if (st->info_level > 1) {
        printf(" %7d x %-7d %8.3f%%", mat->nr, mat->nc, density);
        fflush(stdout);
    }

    return hcm;
}

void convert_sparse_matrix_rows_to_basis_elements(
        mat_t *mat,
        bs_t *bs,
        ht_t *ht,
        hd_t **hcm,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j, k;

    row_t row;
    len_t bl  = bs->ld;

    const len_t np  = mat->np;
    const len_t ncr = mat->ncr;
    const int32_t nthrds  = st->nthrds;

    /* fix size of basis for entering new elements directly */
    printf("np %d\n", np);
    check_enlarge_basis(bs, np);

    for (k = 0, i = 0; i < ncr; ++i) {
        row = mat->npv[i];
        if (row.sz > 0) {
            hd_t **h  =
                (hd_t **)malloc((unsigned long)row.sz * sizeof(hd_t *));
            for (j = 0; j < row.of; ++j) {
                h[j] = insert_in_hash_table(hcm[row.ci[j]]->exp, ht);
            }
            for (; j < row.sz; j += 4) {
                h[j]    = insert_in_hash_table(hcm[row.ci[j]]->exp, ht);
                h[j+1]  = insert_in_hash_table(hcm[row.ci[j+1]]->exp, ht);
                h[j+2]  = insert_in_hash_table(hcm[row.ci[j+2]]->exp, ht);
                h[j+3]  = insert_in_hash_table(hcm[row.ci[j+3]]->exp, ht);
            }

            bs->cf[bl+k]    = row.cl;
            bs->m[bl+k].h   = h;
            bs->m[bl+k].of  = row.of;
            bs->m[bl+k].sz  = row.sz;
            bs->m[bl+k].cl  = row.cl;
            bs->red[bl+k]   = 0;
            bs->lm[bl+k]    = h[0]->sdm;
            free(mat->npv[i].ci);
        printf("new element [%d]\n", bl+k);
        for (int32_t p = 0; p < bs->m[bl+k].sz; ++p) {
            cf16_t *cf  = (cf16_t *)bs->m[bl+k].cl;
            printf("%d | ", cf[p]);
            for (int32_t q = 0; q < gbnv; ++q) {
                printf("%d", bs->m[bl+k].h[p]->exp[q]);
            }
            printf(" || ");
        }
        printf("\n");
        k++;
        }
    }

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->convert_ctime +=  ct1 - ct0;
    st->convert_rtime +=  rt1 - rt0;
}
