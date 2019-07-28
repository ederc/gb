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
 * \file io.c
 * \brief Input and output handling
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

/* note that depending on the input data we set the corresponding
 * function pointers for monomial resp. spair comparisons, taking
 * spairs by a given minimal property for symbolic preprocessing, etc. */
static void import_julia_data_ff(
        bs_t *bs,
        ht_t *ht,
        stat_t *st,
        const int32_t *lens,
        const void *vcfs,
        const int32_t *exps
        )
{
    int32_t i, j;
    len_t k;
    cf32_t * cf;
    hm_t *hm;

    int32_t *cfs  = (int32_t *)vcfs;

    int32_t off       = 0; /* offset in arrays */
    const len_t nv    = st->nvars;
    const len_t ngens = st->ngens;

    int32_t nterms  = 0;
    for (i = 0; i < st->ngens; ++i) {
        nterms  +=  lens[i];
    }
    while (nterms >= ht->esz) {
        enlarge_hash_table(ht);
    }
    exp_t *e  = ht->ev[0]; /* use as temporary storage */
    for (i = 0; i < ngens; ++i) {
        hm  = (hm_t *)malloc(((unsigned long)lens[i]+3) * sizeof(hm_t));
        cf  = (cf32_t *)malloc((unsigned long)(lens[i]) * sizeof(cf32_t));
        bs->hm[i]     = hm;
        bs->cf_ff[i]  = cf;

        hm[0]  = i; /* link to matcf entry */
        hm[1]  = (lens[i] % UNROLL); /* offset */
        hm[2]  = lens[i]; /* length */

        bs->red[i] = 0;

        for (j = off; j < off+lens[i]; ++j) {
            for (k = 0; k < nv; ++k) {
                e[k]  = (exp_t)(exps+(nv*j))[k];
            }
            hm[j-off+3] = insert_in_hash_table(e, ht);
            cf[j-off]   = (cf32_t)cfs[j];
        }
        /* mark initial generators, they have to be added to the basis first */
        off +=  lens[i];
    }
    deg_t deg = 0;
    for (i = 0; i < ngens; ++i) {
        hm  = bs->hm[i];
        deg = ht->hd[hm[3]].deg;
        k   = hm[2] + 3;
        for (j = 4; j < k; ++j) {
            if (deg != ht->hd[hm[j]].deg) {
                st->homogeneous = 0;
                goto done;
            }
        }
    }
    st->homogeneous = 1;
done:

    /* we have to reset the ld value once we have normalized the initial
     * elements in order to start update correctly */
    bs->ld  = st->ngens;
}

static int64_t export_julia_data_ff(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        const bs_t * const bs,
        const ht_t * const ht
        )
{
    len_t i, j, k;

    const len_t nv  = ht->nv;
    const bl_t bld  = bs->ld;

    hm_t *dt;

    int64_t nterms  = 0; /* # of terms in basis */
    int64_t nelts  = 0; /* # elemnts in basis */

    /* compute number of terms */
    for (i = 0; i < bld; ++i) {
        if (bs->red[i] != 0) {
            continue;
        } else {
            nterms +=  (int64_t)bs->hm[i][2];
            nelts++;
        }
    }

    if (nelts > (int64_t)(pow(2, 31))) {
        printf("Basis has more than 2^31 elements, cannot store it.\n");
        return 0;
    }

    int32_t *len  = (int32_t *)malloc(
            (unsigned long)(nelts) * sizeof(int32_t));
    int32_t *exp  = (int32_t *)malloc(
            (unsigned long)(nterms) * (unsigned long)(nv) * sizeof(int32_t));
    cf32_t *cf   = (cf32_t *)malloc(
            (unsigned long)(nterms) * sizeof(cf32_t));

    /* counters for lengths, exponents and coefficients */
    int32_t cl = 0, ce = 0, cc = 0;
    for (i = 0; i < bld; ++i) {
        if (bs->red[i] != 0) {
            continue;
        } else {
            len[cl] = bs->hm[i][2];
            memcpy(cf+cc, bs->cf_ff[bs->hm[i][0]],
                    (unsigned long)(len[cl]) * sizeof(cf32_t));

            dt  = bs->hm[i] + 3;
            for (j = 0; j < len[cl]; ++j) {
                for (k = 0; k < nv; ++k) {
                    exp[ce++] = (int32_t)ht->ev[dt[j]][k]; 
                }
            }
            cc  +=  len[cl];
            cl++;
        }
    }

    *bload  = (int32_t)nelts;
    *blen   = len;
    *bexp   = exp;
    *bcf    = (void *)cf;

    return nterms;
}

static inline void set_function_pointers(
        const stat_t *st
        )
{
    /* todo: this needs to be generalized for different monomial orders */
    switch (st->mo) {
        case 0:
            initial_input_cmp   = initial_input_cmp_drl;
            monomial_cmp        = monomial_cmp_drl;
            spair_cmp           = spair_cmp_drl;
            hcm_cmp             = hcm_cmp_pivots_drl;
            break;
        case 1:
            initial_input_cmp   = initial_input_cmp_lex;
            monomial_cmp        = monomial_cmp_lex;
            spair_cmp           = spair_cmp_deglex;
            hcm_cmp             = hcm_cmp_pivots_lex;
            break;
        default:
            initial_input_cmp   = initial_input_cmp_drl;
            monomial_cmp        = monomial_cmp_drl;
            spair_cmp           = spair_cmp_drl;
            hcm_cmp             = hcm_cmp_pivots_drl;
    }

    switch (st->laopt) {
        case 1:
            linear_algebra  = exact_sparse_dense_linear_algebra_ff;
            break;
        case 2:
            linear_algebra  = exact_sparse_linear_algebra_ff;
            break;
        case 42:
            linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff;
            break;
        case 43:
            linear_algebra  = probabilistic_sparse_dense_linear_algebra_ff_2;
            break;
        case 44:
            linear_algebra  = probabilistic_sparse_linear_algebra_ff;
            break;
        default:
            linear_algebra  = exact_sparse_linear_algebra_ff;
    }

    /* up to 17 bits we can use one modular operation for reducing a row. this works
     * for matrices with #rows <= 54 million */
    if (st->fc == 0) {
        initialize_basis        = initialize_basis_q;
        import_julia_data       = import_julia_data_ff;
        export_julia_data       = export_julia_data_ff;
        check_enlarge_basis     = check_enlarge_basis_q;
        normalize_initial_basis = normalize_initial_basis_q;
    } else {
        initialize_basis        = initialize_basis_ff;
        import_julia_data       = import_julia_data_ff;
        export_julia_data       = export_julia_data_ff;
        check_enlarge_basis     = check_enlarge_basis_ff;
        normalize_initial_basis = normalize_initial_basis_ff;
    }
    if (st->fc < pow(2, 17)) {
        reduce_dense_row_by_all_pivots =
            reduce_dense_row_by_all_pivots_17_bit;
        reduce_dense_row_by_old_pivots =
            reduce_dense_row_by_old_pivots_17_bit;
        reduce_dense_row_by_known_pivots_sparse =
            reduce_dense_row_by_known_pivots_sparse_17_bit;
        reduce_dense_row_by_dense_new_pivots  =
            reduce_dense_row_by_dense_new_pivots_17_bit;
    } else {
        reduce_dense_row_by_all_pivots =
            reduce_dense_row_by_all_pivots_31_bit;
        reduce_dense_row_by_old_pivots =
            reduce_dense_row_by_old_pivots_31_bit;
        reduce_dense_row_by_known_pivots_sparse =
            reduce_dense_row_by_known_pivots_sparse_31_bit;
        reduce_dense_row_by_dense_new_pivots  =
            reduce_dense_row_by_dense_new_pivots_31_bit;
    }
}

static inline int32_t check_and_set_meta_data(
        ps_t *ps,
        stat_t *st,
        const int32_t *lens,
        const int32_t *cfs,
        const int32_t *exps,
        const int32_t field_char,
        const int32_t mon_order,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_hash_table,
        const int32_t la_option,
        const int32_t pbm_file,
        const int32_t info_level
        )
{
    if (nr_gens <= 0
            || nr_vars <= 0
            || field_char <= 0
            || lens == NULL
            || cfs == NULL
            || exps == NULL) {
        return 1;
    }

    st->ngens = nr_gens;
    st->nvars = nr_vars;
    /* note: prime check should be done in julia */
    st->fc    = field_char;
    /* monomial order */
    if (mon_order != 0 && mon_order != 1) {
        st->mo  = 0;
    } else {
        st->mo  = mon_order;
    }
    /* set hash table size */
    st->init_hts  = ht_size;
    if (st->init_hts <= 0) {
        st->init_hts  = 12;
    }
    /* info level */
    st->info_level  = info_level >= 0 ? info_level : 0;
    if (st->info_level > 2) {
        st->info_level = 2;
    }

    /* generation of pbm files on the fly? */
    st->gen_pbm_file  = pbm_file > 0 ? 1 : 0;

    /* resetting basis hash table */
    st->reset_ht = reset_hash_table > 0 ? reset_hash_table : 2147483647; /* 2^31-1 */;

    /* set number of threads */
    if (nr_threads <= 0) {
        st->nthrds  = 1;
    } else {
        st->nthrds  = nr_threads;
    }

    if (max_nr_pairs <= 0) {
        ps->mnsel = 2147483647; /* 2^31-1 */
    } else {
        ps->mnsel = max_nr_pairs;
    }
    st->mnsel = ps->mnsel;

    /* set linear algebra option */
    if (la_option <= 0) {
        st->laopt = 1;
    } else {
        st->laopt = la_option;
    }

    set_function_pointers(st);

    return 0;
}

static void write_pbm_file(
    mat_t *mat,
    const stat_t * const st
    )
{
    len_t i, j, k;
    unsigned char b = 0;
    char buffer[512];
    char fn[200];

    hm_t **rows = mat->r;
    const len_t ncols = mat->nc;
    const len_t nrows = mat->nr;

    sprintf(fn, "%d-%d-%d-%d.pbm", st->current_rd, nrows, ncols, st->current_deg);
    FILE *fh  = fopen(fn, "wb");

    /* magic header */
    sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", nrows, ncols, ncols, nrows);

    fwrite(buffer, sizeof(char), strlen(buffer), fh);


    for (i = 0; i < nrows; ++i) {
        const len_t len = rows[i][2];
        hm_t row[len];
        memcpy(row, rows[i]+3, (unsigned long)len * sizeof(hm_t));
        qsort(row, (unsigned long)len, sizeof(hm_t), pbm_cmp);
        /* the rows may not be sorted by column indices, thus we
         * have to go over them again and again and again */
        k = 0;
        for (j = 0; j < ncols; ++j) {
            if (k < len && row[k] == j) {
                b |=  (unsigned char)(1 << (7 - (j % 8)));
                k++;
            } else {
                b &= (unsigned char) ~(1 << (7 - (j % 8)));
            }
            /* write each byte */
            if (j % 8 == 7) {
                fwrite(&b, sizeof(unsigned char), 1, fh);
                b = 0;
            }
        }
        /* write leftover bits */
        if (j % 8 != 0) {
            fwrite(&b, sizeof(unsigned char), 1, fh);
        }
        fflush(fh);
    }
    fclose(fh);
}
