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
    const int32_t mon_order         = 1;
    const int32_t nr_threads        = 2;
    const int32_t info_level				=	2;
		const int32_t la_option         = 1;
    const int32_t pbm_file          = 0;
    const int32_t max_nr_pairs      = 100;
    const int32_t reset_hash_table  = 0;

    ps_t *ps    = initialize_pairset();
    stat_t *st  = initialize_statistics();
    if (check_and_set_meta_data(ps, st, lens, cfs, exps, field_char, mon_order,
                nr_vars, nr_gens, ht_size, nr_threads, max_nr_pairs,
                reset_hash_table, la_option, pbm_file, info_level)) {
        return 1;
    }

    /* initialize stuff */
    initialize_basis(nr_gens);
    initialize_basis_hash_table();
    if (fc != field_char) {
        return 1;
    }
    if (nvars != nr_vars) {
        return 1;
    }
    initialize_update_hash_table();
    if (hsz/husz != 32) {
        return 1;
    }

    import_julia_data(lens, cfs, exps, nr_gens);

    /* free and clean up */
    free_update_hash_table();
    free_basis_hash_table();
    free_basis();
    free_pairset(&ps);
    return 0;
}
