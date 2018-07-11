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
        int32_t **jl_basis,
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
        const int32_t la_option
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    int32_t round, last_reset;
    hl_t *hcm; /* hash-column-map */
    /* dense matrices, results from linear algebra, before converting
     * the rows back to polynomial representations in the basis */
    cf_t **dm;
    /* matrix holding sparse information generated
     * during symbolic preprocessing */
    dt_t **mat;

    /* checks and set all meta data. if a nonzero value is returned then
     * some of the input data is corrupted. */
    if (check_and_set_meta_data(lens, cfs, exps, field_char, mon_order,
                nr_vars, nr_gens, ht_size, nr_threads, max_nr_pairs, reset_hash_table,
                la_option)) {
        return 0;
    }

    /* initialize stuff */
    initialize_statistics();
    GB_DEBUG(GBDBG, "-------------------------------------------------\n");
    GB_DEBUG(GBDBG, "#variables             %15d\n", nvars);
    GB_DEBUG(GBDBG, "#equations             %15d\n", nr_gens);
    GB_DEBUG(GBDBG, "field characteristic   %15d\n", fc);
    if (mo == 0) {
        GB_DEBUG(GBDBG, "monomial order                     DRL\n");
    }
    if (mo == 1) {
        GB_DEBUG(GBDBG, "monomial order                     LEX\n");
    }
    if ((mo != 0) && (mo != 1)) {
        GB_DEBUG(GBDBG, "monomial order               DONT KNOW\n");
    }
    GB_DEBUG(GBDBG, "linear algebra option  %15d\n", laopt);
    GB_DEBUG(GBDBG, "intial hash table size %15d (2^%d)\n",
            (int32_t)pow(2,htes), htes);
    GB_DEBUG(GBDBG, "reset global hash table %14d\n", rght);
    GB_DEBUG(GBDBG, "maximal pair selection %15d\n", mnsel);
    GB_DEBUG(GBDBG, "#threads               %15d\n", nthrds);
    GB_DEBUG(GBDBG, "-------------------------------------------------\n");

    initialize_basis(nr_gens);
    initialize_pairset();
    initialize_global_hash_table();
    initialize_local_hash_table();

    import_julia_data(lens, cfs, exps, nr_gens);

    /* for faster divisibility checks, needs to be done after we have
     * read some input data for applying heuristics */
    calculate_divmask();

    /* sort initial elements, smallest lead term first */
    qsort(gbdt, (unsigned long)nrows, sizeof(dt_t *),
            matrix_row_initial_input_cmp);
    /* normalize input generators */
    normalize_matrix_rows(gbcf);
    for (int32_t v = 0; v < nr_gens; ++v) {
        printf("new basis element [%d] (%d terms) ", v, gbdt[v][2]-3);
        for (int32_t p = 3; p < gbdt[v][2]; ++p) {
            printf("%d | ", gbcf[v][p]);
            for (int32_t o = 0; o < nvars; ++o) {
                printf("%d ", ev[gbdt[v][p]][o]);
            }
            printf(" || ");
        }
        printf("\n");
    }

    /* move input generators to basis and generate first spairs */
    update_basis();

    /* let's start the f4 rounds,  we are done when no more spairs
     * are left in the pairset */
    last_reset  = 0;
    for (round = 0; pload > 0; ++round) {
        if (round - last_reset == rght) {
            last_reset  = round;
            reset_global_hash_table(); 
        }
        GB_DEBUG(GBDBG, "%3d", round);

        /* preprocess data for next reduction round */
        mat = select_spairs_by_minimal_degree(mat);
        printf("coeffs of selected pairs\n");
        for (int32_t i= 0 ;i < nrows; ++i) {
            for (int32_t j = 3; j < mat[i][2]; ++j) {
                printf("%d ", gbcf[mat[i][0]][j]);
            }
            printf("\n");
        }
        mat = symbolic_preprocessing(mat);
        /* exponent hashes mapped to column indices for linear algebra */
        printf("first column entries in matrix rows\n");
        for (int32_t i= 0 ;i < nrows; ++i) {
            printf("%d -> ", i);
            for (int32_t j = 0; j < nvars; ++j) {
                printf("%d ", ev[mat[i][3]][j]);
            }
            printf("\n");
        }
        hcm = convert_hashes_to_columns(mat);
        printf("first column entries in matrix rows\n");
        for (int32_t i= 0 ;i < nrows; ++i) {
            printf("%d -> %d\n", i, mat[i][3]);
        }
        /* sort matrix rows by decreasing pivots */
        /* for (int32_t o = 0; o < nrows; ++o) {
         *   printf("%d | %d | %d\n", o, mat[o][0], mat[o][2]);
         * } */
        mat = sort_matrix_rows(mat);
        for (int32_t o = 0; o < nrows; ++o) {
          printf("%d | %d | %d | %d\n", o, mat[o][0], mat[0][1], mat[o][2]);
        }
        /* linear algebra, depending on choice, see set_function_pointers() */
        dm = linear_algebra(mat);
        free(mat);
        mat = NULL;
        /* columns indices are mapped back to exponent hashes */
        if (npivs > 0) {
            convert_dense_matrix_to_basis_elements(dm, hcm);
            free_dense_matrix(dm);
        }
        hcm = reset_idx_in_global_hash_table_and_free_hcm(hcm);
        update_basis();

        GB_DEBUG(GBDBG, "\n");
    }

    int64_t len = export_julia_data(jl_basis);

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    GB_DEBUG(GBDBG, "-------------------------------------------------\n");
    GB_DEBUG(GBDBG, "overall                %15.3f sec\n", rt1-rt0);
    GB_DEBUG(GBDBG, "overall(cpu)           %15.3f sec\n", ct1-ct0);
    GB_DEBUG(GBDBG, "select                 %15.3f sec\n", select_rtime);
    GB_DEBUG(GBDBG, "pair sort              %15.3f sec\n", pair_sort_rtime);
    GB_DEBUG(GBDBG, "symbol                 %15.3f sec\n", symbol_rtime);
    GB_DEBUG(GBDBG, "update                 %15.3f sec\n", update_rtime);
    GB_DEBUG(GBDBG, "update1                %15.3f sec\n", update1_rtime);
    GB_DEBUG(GBDBG, "convert                %15.3f sec\n", convert_rtime);
    if (rght != 0) {
        GB_DEBUG(GBDBG, "rght                   %15.3f sec\n", rght_rtime);
    }
    GB_DEBUG(GBDBG, "la                     %15.3f sec\n", la_rtime);
    GB_DEBUG(GBDBG, "-------------------------------------------------\n");
    GB_DEBUG(GBDBG, "size of basis          %15d\n", (*jl_basis[0]));
    GB_DEBUG(GBDBG, "#terms in basis        %15ld\n",
            (len-(*jl_basis)[0]-1)/(1+nvars));
    GB_DEBUG(GBDBG, "#pairs reduced         %15ld\n", num_pairsred);
    GB_DEBUG(GBDBG, "#GM criterion          %15ld\n", num_gb_crit);
    GB_DEBUG(GBDBG, "#redundant             %15ld\n", num_redundant);
    GB_DEBUG(GBDBG, "#rows reduced          %15ld\n", num_rowsred);
    GB_DEBUG(GBDBG, "#zero reductions       %15ld\n", num_zerored);
    GB_DEBUG(GBDBG, "global hash table load %15d <= 2^%d\n",
            eld, (int32_t)((ceil(log(eld)/log(2)))));
    GB_DEBUG(GBDBG, "local hash table load  %15d <= 2^%d\n",
            elld, (int32_t)(ceil(log(elld)/log(2))));
    GB_DEBUG(GBDBG, "#ht enlargements       %15ld\n", num_htenl);
    GB_DEBUG(GBDBG, "#no reducer found      %15ld\n", num_sdm_found + num_not_sdm_found);
    GB_DEBUG(GBDBG, "sdm findings           %14.3f%% \n",
            (double)100*(double)num_sdm_found/(double)(num_sdm_found + num_not_sdm_found));
    GB_DEBUG(GBDBG, "-------------------------------------------------\n");

    /* free and clean up */
    free_local_hash_table();
    free_global_hash_table();
    free_pairset();
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);
    free_basis();

    return len;
}
