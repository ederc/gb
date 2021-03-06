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

#include "data.h"

void free_julia_data(
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,      /* coefficients of basis elements */
        const int64_t ngens,
        const int64_t field_char
        )
{
    int64_t i;
    int64_t len = 0;

    /* lengths resp. nterms */
    int32_t *lens  = *blen;
    for (i = 0; i < ngens; ++i) {
        len += (int64_t)lens[i];
    }

    free(lens);
    lens = NULL;
    *blen = lens;

    /* exponent vectors */
    int32_t *exps = *bexp;
    free(exps);
    exps  = NULL;
    *bexp = exps;
    bexp  = NULL;

    /* coefficients */
    if (field_char == 0) {
        mpz_t **cfs = (mpz_t **)bcf;
        for (i = 0; i < len; ++i) {
            mpz_clear((*cfs)[i]);
        }
        free(*cfs);
        free(cfs);
        cfs = NULL;
    } else {
        if (field_char > 0) {
            int32_t *cfs  = *((int32_t **)bcf);
            free(cfs);
            cfs = NULL;
        }
    }
    bcf = NULL;
}

static void clear_matrix(
        mat_t *mat
        )
{
    free(mat->r);
    mat->r  = NULL;
    free(mat->cf_8);
    mat->cf_8 = NULL;
    free(mat->cf_16);
    mat->cf_16  = NULL;
    free(mat->cf_32);
    mat->cf_32  = NULL;
    free(mat->cf_qq);
    mat->cf_qq  = NULL;
    free(mat->cf_ab_qq);
    mat->cf_ab_qq  = NULL;
}

static void reduce_basis(
        bs_t *bs,
        mat_t *mat,
        hl_t **hcmp,
        ht_t *bht,
        ht_t *sht,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j, k;

    hl_t *hcm   = *hcmp;
    exp_t *etmp = bht->ev[0];
    memset(etmp, 0, (unsigned long)bht->nv * sizeof(exp_t));

    mat->r  = (hm_t **)malloc((unsigned long)bs->lml * 2 * sizeof(hm_t *));
    mat->nr = 0;
    mat->sz = 2 * bs->lml;

    /* add all non-redundant basis elements as matrix rows */
    for (i = 0; i < bs->lml; ++i) {
        mat->r[mat->nr] = multiplied_poly_to_matrix_row(
                sht, bht, 0, 0, etmp, bs->hm[bs->lmps[i]]);
        sht->hd[mat->r[mat->nr][3]].idx  = 1;
        mat->nr++;
    }
    symbolic_preprocessing(mat, bs, st, sht, bht);
    /* no known pivots, we need mat->ncl = 0, so set all indices to 1 */
    for (i = 0; i < sht->eld; ++i) {
        sht->hd[i].idx = 1;
    }

    /* generate hash <-> column mapping */
    if (st->info_level > 1) {
        printf("reduce final basis ");
        fflush(stdout);
    }
    convert_hashes_to_columns(&hcm, mat, st, sht);
    mat->nc = mat->ncl + mat->ncr;
    /* sort rows */
    sort_matrix_rows(mat);
    /* do the linear algebra reduction */
    interreduce_matrix_rows(mat, bs, st);
    /* remap rows to basis elements (keeping their position in bs) */
    /* convert_sparse_matrix_rows_to_final_basis(mat, bs, bht, sht, hcm, st); */
    convert_sparse_matrix_rows_to_basis_elements(
        mat, bs, bht, sht, hcm, st);

    bs->ld  = mat->np;

    clean_hash_table(sht);
    clear_matrix(mat);

    /* we may have added some multiples of reduced basis polynomials
     * from the matrix, so we get rid of them. */
    k = 0;
    i = bs->ld-1;
start:
    for (; i >= 0; --i) {
        for (j = 0; j < k; ++j) {
            if (check_monomial_division(
                        bs->hm[i][3],
                        bs->hm[bs->lmps[j]][3], bht)) {
                --i;
                goto start;
            }
        }
        bs->lmps[k++] = i;
    }
    *hcmp = hcm;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->reduce_gb_ctime = ct1 - ct0;
    st->reduce_gb_rtime = rt1 - rt0;
    if (st->info_level > 1) {
        printf("%13.3f sec\n", rt1-rt0);
    }

    if (st->info_level > 1) {
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
}

/* we get from julia the generators as three arrays:
 * 0.  a pointer to an int32_t array for returning the basis to julia
 * 1.  an array of the lengths of each generator
 * 2.  an array of all coefficients of all generators in the order:
 *     first all coefficients of generator 1, then all of generator 2, ...
 * 3.  an array of all exponents of all generators in the order:
 *     first all exponents of generator 1, then all of generator 2, ...
 *
 *  RETURNs the length of the jl_basis array */
int64_t f4_julia(
        /* return values */
        int32_t *bld,   /* basis load */
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,     /* coefficients of basis elements */
        /* input values */
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const int32_t field_char,
        const int32_t mon_order,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_ht,
        const int32_t la_option,
        const int32_t reduce_gb,
        const int32_t pbm_file,
        const int32_t info_level
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    double rrt0, rrt1; /* for one round only */
    ct0 = cputime();
    rt0 = realtime();

    int32_t round, i, j;
    /* hashes-to-columns map, initialized with length 1, is reallocated
     * in each call when generating matrices for linear algebra */
    hl_t *hcm = (hl_t *)malloc(sizeof(hl_t));
    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    ps_t * ps   = initialize_pairset();

    /* initialize stuff */
    stat_t *st  = initialize_statistics();

    /* checks and set all meta data. if a nonzero value is returned then
     * some of the input data is corrupted. */
    if (check_and_set_meta_data(ps, st, lens, exps, cfs, field_char, mon_order,
                nr_vars, nr_gens, ht_size, nr_threads, max_nr_pairs,
                reset_ht, la_option, reduce_gb, pbm_file, info_level)) {
        return 0;
    }

    /* initialize basis */
    bs_t *bs  = initialize_basis(st->ngens);
    /* initialize basis hash table, update hash table, symbolic hash table */
    ht_t *bht = initialize_basis_hash_table(st);
    ht_t *uht = initialize_secondary_hash_table(bht, st);
    ht_t *sht = initialize_secondary_hash_table(bht, st);

    import_julia_data(bs, bht, st, lens, exps, cfs);

    if (st->info_level > 0) {
        print_initial_statistics(st);
    }

    /* for faster divisibility checks, needs to be done after we have
     * read some input data for applying heuristics */
    calculate_divmask(bht);

    /* sort initial elements, smallest lead term first */
    sort_r(bs->hm, (unsigned long)bs->ld, sizeof(hm_t *),
            initial_input_cmp, bht);
    /* normalize input generators */
    if (st->fc > 0) {
        normalize_initial_basis(bs, st->fc);
    } else {
        if (st->fc == 0) {
            remove_content_of_initial_basis(bs);
        }
    }

    /* reset bs->ld for first update process */
    bs->ld  = 0;

    /* move input generators to basis and generate first spairs.
     * always check redundancy since input generators may be redundant
     * even so they are homogeneous. */
    update_basis(ps, bs, bht, uht, st, st->ngens, 1);

    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    if (st->info_level > 1) {
        printf("\ndeg     sel   pairs        mat          density \
          new data             time(rd)\n");
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    for (round = 1; ps->ld > 0; ++round) {
      if (round % st->reset_ht == 0) {
        reset_hash_table(bht, bs, ps, st);
        st->num_rht++;
      }
      rrt0  = realtime();
      st->max_bht_size  = st->max_bht_size > bht->hsz ?
        st->max_bht_size : bht->hsz;
      st->current_rd  = round;

      /* preprocess data for next reduction round */
      select_spairs_by_minimal_degree(mat, bs, ps, st, sht, bht);
      symbolic_preprocessing(mat, bs, st, sht, bht);
      convert_hashes_to_columns(&hcm, mat, st, sht);
      sort_matrix_rows(mat);
      /* print pbm files of the matrices */
      if (st->gen_pbm_file != 0) {
        write_pbm_file(mat, st);
      }
      /* linear algebra, depending on choice, see set_function_pointers() */
      linear_algebra(mat, bs, st);
      /* columns indices are mapped back to exponent hashes */
      if (mat->np > 0) {
        convert_sparse_matrix_rows_to_basis_elements(
            mat, bs, bht, sht, hcm, st);
      }
      clean_hash_table(sht);
      /* all rows in mat are now polynomials in the basis,
       * so we do not need the rows anymore */
      clear_matrix(mat);

      /* check redundancy only if input is not homogeneous */
      update_basis(ps, bs, bht, uht, st, mat->np, 1-st->homogeneous);

      rrt1 = realtime();
      if (st->info_level > 1) {
        printf("%13.3f sec\n", rrt1-rrt0);
      }
    }
    if (st->info_level > 1) {
        printf("-------------------------------------------------\
----------------------------------------\n");
    }
    /* remove possible redudant elements */
    j = 0;
    for (i = 0; i < bs->lml; ++i) {
        if (bs->red[bs->lmps[i]] == 0) {
            bs->lm[j]   = bs->lm[i];
            bs->lmps[j] = bs->lmps[i];
            ++j;
        }
    }
    bs->lml = j;

    /* reduce final basis? */
    if (st->reduce_gb == 1) {
        reduce_basis(bs, mat, &hcm, bht, sht, st);
    }

    st->nterms_basis  = export_julia_data(bld, blen, bexp, bcf, bs, bht);
    st->size_basis    = *bld;

    const int64_t nterms  = st->nterms_basis;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->overall_ctime = ct1 - ct0;
    st->overall_rtime = rt1 - rt0;

    if (st->info_level > 0) {
        print_final_statistics(st);
    }

    /* free and clean up */
    free(hcm);
    free_shared_hash_data(bht);
    free_hash_table(&sht);
    free_hash_table(&uht);
    free_hash_table(&bht);
    free_pairset(&ps);
    free_basis(&bs);
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);

    free(st);

    return nterms;
}
