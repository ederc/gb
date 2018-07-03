#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    int32_t i, j, k;

    const int32_t lens[]  = {2,2,2}; 
    const int32_t cfs[]   = {1, 1, 1, 3, 5, 3};
    const int32_t exps[]  = {1, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0};

    const int32_t nr_vars           = 2;
    const int32_t nr_gens           = 3;
    const int32_t ht_size           = 12;
    const int32_t field_char        = 101;
    const int32_t mon_order         = 0;
    const int32_t nr_threads        = 1;
    const int32_t la_option         = 1;
    const int32_t max_nr_pairs      = 100;
    const int32_t reset_hash_table  = 0;

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
        mat = symbolic_preprocessing(mat);
        /* for (int32_t o = 0; o < nrows; ++o) {
         *   printf("%d | %d | %d (%d) || ", o, mat[o][0], mat[o][2], (ev+mat[o][2])[HASH_IND]);
         *   for (int32_t p = 0; p < nvars; ++p) {
         *     printf("%d ", (ev+mat[o][2])[p]);
         *   }
         *   printf("\n");
         * } */
        /* exponent hashes mapped to column indices for linear algebra */
        hcm = convert_hashes_to_columns(mat);
        /* sort matrix rows by decreasing pivots */
        /* for (int32_t o = 0; o < nrows; ++o) {
         *   printf("%d | %d | %d\n", o, mat[o][0], mat[o][2]);
         * } */
        mat = sort_matrix_rows(mat);
        /* for (int32_t o = 0; o < nrows; ++o) {
         *   printf("%d | %d | %d\n", o, mat[o][0], mat[o][2]);
         * } */
        /* linear algebra, depending on choice, see set_function_pointers() */
        dm = linear_algebra(mat);
        /* columns indices are mapped back to exponent hashes */
        convert_dense_matrix_to_basis_elements(dm, hcm);
        /* mat = convert_columns_to_hashes(mat, hcm); */

        free(hcm);
        hcm = NULL;

        update_basis();

        GB_DEBUG(GBDBG, "\n");
    }


    /* free and clean up */
    free_local_hash_table();
    free_global_hash_table();
    free_pairset();
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free(mat);
    free_basis();

    int32_t **basis = (int32_t **)malloc(sizeof(int32_t *));
    int64_t len = export_julia_data(basis);

    int32_t val[12] = {2, 6, 3, 1, 1, 0, 1, 0, 1, 1, 0, 1};

    int32_t failure = 0;
    for (i = 0; i < 12; ++i) {
        if (val[i] != (*basis)[i]) {
            failure = 1;
            break;
        }
    }
    free(*basis);
    free(basis);

    basis = NULL;

    /* free and clean up */
    free_local_hash_table();
    free_global_hash_table();
    free_pairset();
    /* note that all rows kept from mat during the overall computation are
     * basis elements and thus we do not need to free the rows itself, but
     * just the matrix structure */
    free_basis();

    return failure;
}
