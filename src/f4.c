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
 * \file f4.c
 * \brief Implementation of the F4 Algorithm
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "f4.h"

/* we get from julia the generators as three arrays:
 * 0.  a pointer to an int32_t array for returning the basis to julia
 * 1.  an array of the lengths of each generator
 * 2.  an array of all coefficients of all generators in the order:
 *     first all coefficients of generator 1, then all of generator 2, ...
 * 3.  an array of all exponents of all generators in the order:
 *     first all exponents of generator 1, then all of generator 2, ...
 *
 *  RETURNs the length of the jl_basis array */
int64_t f4_julia_ff(
        int32_t **jl_basis,
        const int32_t *lens,
        const void *cfs,
        const int32_t *exps,
        const int32_t field_char,
        const int32_t mon_order,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t regenerate_ht,
        const int32_t la_option,
        const int32_t info_level
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    double rct0, rct1, rrt0, rrt1; /* for one round only */
    ct0 = cputime();
    rt0 = realtime();

    int32_t rd; /* number of F4 round */
    int32_t lr; /* last regeneration of global hash table */
    hd_t **hcm; /* hash-column-map */

    /* initialize stuff */
    stat_t *st  = initialize_statistics();
    /* checks and set all meta data. if a nonzero value is returned then
     * some of the input data is corrupted. */
    gbnv  = 2;
    printf("gbnv global %d\n", gbnv);
    printf("gbfc global %d\n", gbfc);
    if (check_and_set_meta_data(st, lens, cfs, exps, field_char, mon_order,
                nr_vars, nr_gens, ht_size, nr_threads, max_nr_pairs,
                regenerate_ht, la_option, info_level)) {
        return 0;
    }

    printf("gbnv global %d\n", gbnv);
    printf("gbfc global %d\n", gbfc);
    if (st->info_level > 0) {
        print_initial_statistics(st);
    }

    /* global hash table */
    ht_t *ght = initialize_global_hash_table(st);
    /* upate hash table for spairs */
    ht_t *uht = initialize_local_hash_table(st, ght);
    /* hash table for symbolic preprocessing */
    ht_t *sht = initialize_local_hash_table(st, ght);
    etmp  = (exp_t *)malloc((unsigned long)gbnv * sizeof(exp_t));

    bs_t *bs    = initialize_basis(st);
    ps_t *ps    = initialize_pairset(st);
    mat_t *mat  = initialize_matrix();

    import_julia_data(bs, ght, lens, cfs, exps, nr_gens);

    /* for faster divisibility checks, needs to be done after we have
     * read some input data for applying heuristics */
    calculate_divmask(ght);
    /* set short divisor mask also for local hash table */
    sht->dm = uht->dm = ght->dm;

    st->num_matrices++;

    for (int ii=0; ii < st->nr_gens; ++ii) {
        printf("gen[%d] = ", ii);
        cf16_t *cf = (cf16_t *)bs->m[ii].cl;
        for (int jj=0; jj < bs->m[ii].sz; ++jj) {
            printf("%d | ", cf[jj]);
            for (int kk=0; kk < gbnv; ++kk) {
                printf("%d", bs->m[ii].h[jj]->exp[kk]);
            }
            printf(" || ");
        }
        printf("\n");
    }
    /* sort initial elements, smallest lead term first */
    qsort(bs->m, (unsigned long)st->nr_gens,
            sizeof(mon_t), initial_basis_cmp);
    /* normalize input generators */
    normalize_initial_basis(bs, st->nr_gens);
    /* normalize_matrix_rows(bs->cf, bs->ld, st->field_char); */

    /* move input generators to basis and generate first spairs */
    check_enlarge_pairset(ps, bs->ld, st->nr_gens);
    /* generate initial spairs */
    update_basis(ps, bs, ght, uht, st, st->nr_gens);

    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    lr  = 0;
    if (st->info_level > 1) {
        print_round_statistics_header();
    }
    for (rd = 0; ps->ld > 0; ++rd) {
        st->num_matrices++;
        rct0 = cputime();
        rrt0 = realtime();

        lr  = check_regenerate_hash_table(ght, ps, bs, st, lr, rd);

        /* preprocess data for next reduction round */
        mat = select_spairs_by_minimal_degree(mat, ps, sht, bs, st);
        mat = symbolic_preprocessing(mat, sht, bs, st);
        hcm = convert_hashes_to_columns(mat, sht, st);
        printf("pv\n");
        for (int ii=0; ii<mat->nru; ++ii) {
            for (int jj=0; jj<mat->pv[ii].sz; ++jj) {
                printf("%d ", mat->pv[ii].ci[jj]);
            }
            printf("\n");
        }
        printf("tbr\n");
        for (int ii=0; ii<mat->nrl; ++ii) {
            for (int jj=0; jj<mat->tbr[ii].sz; ++jj) {
                printf("%d ", mat->tbr[ii].ci[jj]);
            }
            printf("\n");
        }
        mat = sort_matrix_rows(mat);
        printf("pv\n");
        for (int ii=0; ii<mat->nru; ++ii) {
            for (int jj=0; jj<mat->pv[ii].sz; ++jj) {
                printf("%d ", mat->pv[ii].ci[jj]);
            }
            printf("\n");
        }
        printf("tbr\n");
        for (int ii=0; ii<mat->nrl; ++ii) {
            for (int jj=0; jj<mat->tbr[ii].sz; ++jj) {
                printf("%d ", mat->tbr[ii].ci[jj]);
            }
            printf("\n");
        }
        /* linear algebra, depending on choice, see set_function_pointers() */
        mat = linear_algebra(mat, st);
        /* columns indices are mapped back to exponent hashes */
        if (mat->np > 0) {
            convert_sparse_matrix_rows_to_basis_elements(
                    mat, bs, ght, hcm, st);
        }
        free(hcm);
        hcm = NULL;
        reset_symbolic_hash_table(sht);


        check_enlarge_pairset(ps, bs->ld, mat->np);
        update_basis(ps, bs, ght, uht, st, mat->np);

        rct1 = cputime();
        rrt1 = realtime();
        if (st->info_level > 1) {
            printf("%13.3f sec\n", rrt1-rrt0);
        }
    }
    if (st->info_level > 1) {
        print_round_statistics_footer();
    }

    st->len_output  = export_julia_data(bs, jl_basis);
    st->size_basis  = (*jl_basis)[0];

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    if (st->info_level > 0) {
        print_final_statistics(st, ght->eld, uht->eld);
    }

    /* free and clean up */
    free_hash_table(&ght);
    /* rv and dm is shared with ght, thus already removed. we have
     * to set them to NULL in order to prevent a second free'ing. */
    uht->rv = NULL;
    uht->dm = NULL;
    sht->rv = NULL;
    sht->dm = NULL;
    free_hash_table(&uht);
    free_hash_table(&sht);
    free_pairset(&ps);
    free_matrix(&mat);
    free_basis(&bs);
    free(etmp);

    int64_t output  = st->len_output;
    free(st);

    return output;
}
