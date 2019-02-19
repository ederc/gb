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

static inline void set_function_pointers(
        void
        )
{
    /* todo: this needs to be generalized for different monomial orders */
    switch (mo) {
        case 0:
            matrix_row_initial_input_cmp  =
                matrix_row_initial_input_cmp_drl;
            monomial_cmp        = monomial_cmp_drl;
            monomial_update_cmp = monomial_update_cmp_drl;
            spair_cmp           = spair_cmp_drl;
            hcm_cmp             = hcm_cmp_pivots_drl;
            break;
        case 1:
            matrix_row_initial_input_cmp  =
                matrix_row_initial_input_cmp_lex;
            monomial_cmp        = monomial_cmp_lex;
            monomial_update_cmp = monomial_update_cmp_lex;
            spair_cmp           = spair_cmp_deglex;
            hcm_cmp             = hcm_cmp_pivots_lex;
            break;
        default:
            matrix_row_initial_input_cmp  =
                matrix_row_initial_input_cmp_drl;
            monomial_cmp        = monomial_cmp_drl;
            monomial_update_cmp = monomial_update_cmp_drl;
            spair_cmp           = spair_cmp_drl;
            hcm_cmp             = hcm_cmp_pivots_drl;
    }

    switch (laopt) {
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
            linear_algebra  = exact_sparse_dense_linear_algebra_ff;
    }

    /* up to 17 bits we can use one modular operation for reducing a row. this works
     * for matrices with #rows <= 54 million */
    if (fc == 0) {
        initialize_basis        = initialize_basis_q;
        check_enlarge_basis     = check_enlarge_basis_q;
        free_basis              = free_basis_q;
        normalize_initial_basis = normalize_initial_basis_q;
    } else {
        initialize_basis        = initialize_basis_ff;
        check_enlarge_basis     = check_enlarge_basis_ff;
        free_basis              = free_basis_ff;
        normalize_initial_basis = normalize_initial_basis_ff;
    }
    if (fc < pow(2, 17)) {
        reduce_dense_row_by_all_pivots =
            reduce_dense_row_by_all_pivots_17_bit;
        reduce_dense_row_by_known_pivots =
            reduce_dense_row_by_known_pivots_17_bit;
        reduce_dense_row_by_known_pivots_sparse =
            reduce_dense_row_by_known_pivots_sparse_17_bit;
        reduce_dense_row_by_dense_new_pivots  =
            reduce_dense_row_by_dense_new_pivots_17_bit;
    } else {
        reduce_dense_row_by_all_pivots =
            reduce_dense_row_by_all_pivots_31_bit;
        reduce_dense_row_by_known_pivots =
            reduce_dense_row_by_known_pivots_31_bit;
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

    ngens = nr_gens;
    nvars = nr_vars;
    /* note: prime check should be done in julia */
    fc    = field_char;
    /* monomial order */
    if (mon_order != 0 && mon_order != 1) {
        mo  = 0;
    } else {
        mo  = mon_order;
    }
    /* set hash table size */
    htes  = ht_size;
    if (htes <= 0) {
        htes  = 12;
    }
    /* info level */
    st->info_level  = info_level >= 0 ? info_level : 0;
    if (st->info_level > 2) {
        st->info_level = 2;
    }

    /* generation of pbm files on the fly? */
    st->gen_pbm_file  = pbm_file > 0 ? 1 : 0;

    /* resetting basis hash table */
    rht = reset_hash_table > 0 ? reset_hash_table : 2147483647; /* 2^31-1 */;

    /* set number of threads */
    if (nr_threads <= 0) {
        nthrds  = 1;
    } else {
        nthrds  = nr_threads;
    }

    if (max_nr_pairs <= 0) {
        ps->mnsel = 2147483647; /* 2^31-1 */
    } else {
        ps->mnsel = max_nr_pairs;
    }
    mnsel = ps->mnsel;

    /* set linear algebra option */
    if (la_option <= 0) {
        laopt = 1;
    } else {
        laopt = la_option;
    }

    set_function_pointers();

    return 0;
}

/* note that depending on the input data we set the corresponding
 * function pointers for monomial resp. spair comparisons, taking
 * spairs by a given minimal property for symbolic preprocessing, etc. */
static void import_julia_data_ff(
        const int32_t *lens,
        const int32_t *cfs,
        const int32_t *exps,
        const int32_t nr_gens
        )
{
    int32_t i, j;
    len_t k;
    cf32_t * cf;
    dt_t *dt;

    int32_t off = 0; /* offset in arrays */
    exp_t *e  = (exp_t *)malloc((unsigned long)nvars * sizeof(exp_t));

    for (i = 0; i < nr_gens; ++i) {
        /* gbdt rows store three more elemets at the beginning:
         * gbdt[0] is the index j of the corresponding coefficient array gbcf[j]
         * gbdt[1] is the offset of the length of the array for loop unrolling
         * gbdt[2] is the real length of the array for looping */
        gbdt[i]     = (dt_t *)malloc(((unsigned long)lens[i]+3) * sizeof(dt_t));
        gbcf_ff[i]  = (cf32_t *)malloc((unsigned long)(lens[i]) * sizeof(cf32_t));
        gbdt[i][0]  = i; /* link to matcf entry */
        red[i]      = 0;
        gbdt[i][1]  = (lens[i] % UNROLL); /* offset */
        gbdt[i][2]  = lens[i]; /* length */

        cf  = gbcf_ff[i];
        dt  = gbdt[i] + 3;
        for (j = off; j < off+lens[i]; ++j) {
            for (k = 0; k < nvars; ++k) {
                e[k]  = (exp_t)(exps+(nvars*j))[k];
            }
            dt[j-off]  = insert_in_basis_hash_table(e);
            cf[j-off]  = (cf32_t)cfs[j];
        }
        /* mark initial generators, they have to be added to the basis first */
        off +=  lens[i];
    }
    npivs = nrows = nrall = nr_gens;
    free(e);
}

static int64_t export_julia_data_ff(
        int32_t **bp
        )
{
    len_t i, j, k;
    int64_t ctr_lengths, ctr_elements;
    int32_t *basis  = *bp;

    cf32_t *cf;
    dt_t *dt;

    int64_t len = 0; /* complete length of exported array */
    int64_t nb  = 0; /* # elemnts in basis */

    const len_t lterm = 1 + nvars; /* length of a term */

    /* compute number of terms */
    for (i = 0; i < bload; ++i) {
        if (red[i]) {
            continue;
        } else {
            len +=  (int64_t)gbdt[i][2];
            nb++;
        }
    }

    /* compute the length considering the number of variables per exponent */
    len = len * (int64_t)lterm;
    /* add storage for length of each element */
    len = len + nb;
    /* add storage for number of generators in basis */
    len++;

    basis  = (int32_t *)malloc((unsigned long)len * sizeof(int32_t));

    if (nb > (int64_t)(pow(2, 31))) {
        printf("basis too big\n");
        return 0;
    }

    ctr_lengths   = 1;
    ctr_elements  = (int64_t)nb + 1;

    basis[0]  = (int32_t)nb;
    for (i = 0; i < bload; ++i) {
        if (red[i]) {
            continue;
        } else {
            /* length of polynomial including this length entry itself */
            basis[ctr_lengths++]  = (int32_t)((gbdt[i][2]) * lterm);
            cf  = gbcf_ff[gbdt[i][0]];
            dt  = gbdt[i] + 3;
            for (j = 0; j < gbdt[i][2]; ++j) {
                basis[ctr_elements++] = (int32_t)cf[j]; /* coefficient */
                for (k = 0; k < nvars; ++k) {
                    basis[ctr_elements++] = (int32_t)ev[dt[j]][k];
                }
            }
        }
    }
    *bp = basis;

    return len;
}

static void write_pbm_file(
    dt_t **mat,
    stat_t *st
    )
{
    len_t i, j, k;
    unsigned char b = 0;
    char buffer[512];
    char fn[200];
    sprintf(fn, "%d-%d-%d-%d.pbm", st->current_rd, nrows, ncols, st->current_deg);
    FILE *fh  = fopen(fn, "wb");

    /* magic header */
    sprintf(buffer, "P4\n# matrix size(%u, %u)\n%u %u\n", nrows, ncols, ncols, nrows);

    fwrite(buffer, sizeof(char), strlen(buffer), fh);


    for (i = 0; i < nrows; ++i) {
        const len_t len = mat[i][2];
        dt_t row[len];
        memcpy(row, mat[i]+3, (unsigned long)len * sizeof(dt_t));
        qsort(row, (unsigned long)len, sizeof(dt_t), pbm_cmp);
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
