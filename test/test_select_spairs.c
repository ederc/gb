#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    int32_t i, j;

    const int32_t lens[]  = {2,2,2}; 
    const int32_t cfs[]   = {1, 1, 1, 3, 5, 3};
    const int32_t exps[]  = {1, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0};

    const int32_t nr_vars           = 2;
    const int32_t nr_gens           = 3;
    const int32_t ht_size           = 12;
    const int32_t field_char        = 101;
    const int32_t mon_order         = 0;
    const int32_t nr_threads        = 2;
    const int32_t info_level				=	2;
		const int32_t la_option         = 1;
    const int32_t max_nr_pairs      = 100;
    const int32_t reset_hash_table  = 0;

    ps_t *ps  = initialize_pairset();
    if (check_and_set_meta_data(ps, lens, cfs, exps, field_char, mon_order,
                nr_vars, nr_gens, ht_size, nr_threads, max_nr_pairs, reset_hash_table,
                la_option, info_level)) {
        return 1;
    }

    dt_t **mat;

    /* initialize stuff */
    initialize_basis(nr_gens);
    initialize_global_hash_table();
    initialize_local_hash_table();

    stat_t *st  = initialize_statistics();
    if (hsz/hlsz != 32) {
        return 1;
    }

    import_julia_data(lens, cfs, exps, nr_gens);

    /* for faster divisibility checks, needs to be done after we have
     * read some input data for applying heuristics */
    calculate_divmask();

    /* sort initial elements, smallest lead term first */
    qsort(gbdt, (unsigned long)nrows, sizeof(dt_t *),
            matrix_row_initial_input_cmp);
    /* normalize input generators */
    normalize_matrix_rows(gbcf, gbdt);

    /* move input generators to basis and generate first spairs */
    update_basis(ps, st);

    mat = select_spairs_by_minimal_degree(ps, mat, st);

    /* free and clean up */
    free_local_hash_table();
    free_global_hash_table();
    free_basis();
    free_pairset(&ps);
    for (i = 0; i < nrows; ++i) {
        free(mat[i]);
    }
    free(mat);
    return 0;
}
