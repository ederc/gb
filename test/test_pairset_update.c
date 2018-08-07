#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    int32_t i;

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
    const int32_t max_nr_pairs      = 10;
    const int32_t reset_hash_table  = 0;

    if (check_and_set_meta_data(lens, cfs, exps, field_char, mon_order,
                nr_vars, nr_gens, ht_size, nr_threads, max_nr_pairs, reset_hash_table,
                la_option, info_level)) {
        return 1;
    }

    /* initialize stuff */
    initialize_basis(nr_gens);
    initialize_pairset();
    initialize_global_hash_table();
    initialize_local_hash_table();

    import_julia_data(lens, cfs, exps, nr_gens);
    calculate_divmask();
    /* sort initial elements, smallest lead term first */
    qsort(gbdt, (unsigned long)nrows, sizeof(dt_t *),
            matrix_row_initial_input_cmp);
    /* normalize input generators */
    normalize_matrix_rows(gbcf);

    stat_t *st  = initialize_statistics();
    /* move input generators to basis and generate first spairs */
    update_basis(st);

    if (pload != 2) {
        printf("pload wrong - %d != 2\n", pload);
        return 1;
    }

    if (ps[0].gen1 != 0) {
        printf("ps[0].gen1 wrong - %d != 1\n", ps[0].gen1);
        return 1;
    }
    if (ps[0].gen2 != 1) {
        printf("ps[0].gen2 wrong - %d != 2\n", ps[0].gen2);
        return 1;
    }
    if (hd[ps[0].lcm].deg != 2) {
        printf("ps[0].lcm.deg wrong - %d != 2\n", hd[ps[0].lcm].deg);
        return 1;
    }

    if (ps[1].gen1 != 0) {
        printf("ps[1].gen1 wrong - %d != 0\n", ps[1].gen1);
        return 1;
    }
    if (ps[1].gen2 != 2) {
        printf("ps[1].gen2 wrong - %d != 2\n", ps[1].gen2);
        return 1;
    }
    if (hd[ps[1].lcm].deg != 3) {
        printf("ps[1].lcm.deg wrong - %d != 3\n", hd[ps[1].lcm].deg);
        return 1;
    }

    /* free and clean up */
    free_local_hash_table();
    free_global_hash_table();
    free_pairset();
    free_basis();
    return 0;
}
