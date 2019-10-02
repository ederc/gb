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
    const int32_t pbm_file          = 0;
    const int32_t max_nr_pairs      = 10;
    const int32_t reset_hash_table  = 0;

    ps_t *ps    = initialize_pairset();
    stat_t *st  = initialize_statistics();
    if (check_and_set_meta_data(ps, st, lens, exps, cfs, field_char, mon_order,
                nr_vars, nr_gens, ht_size, nr_threads, max_nr_pairs,
                reset_hash_table, la_option, pbm_file, info_level)) {
        return 1;
    }

    /* initialize stuff */
    bs_t *bs  = initialize_basis_ff(st->ngens);
    ht_t *bht = initialize_basis_hash_table(st);
    ht_t *uht = initialize_secondary_hash_table(bht, st);

    import_julia_data_ff(bs, bht, st, lens, exps, cfs);
    calculate_divmask(bht);
    /* sort initial elements, smallest lead term first */
    sort_r(bs->hm, (unsigned long)bs->ld, sizeof(hm_t *),
            initial_input_cmp, bht);
    /* normalize input generators */
    normalize_initial_basis_ff(bs, st->fc);
    /* reset bs->ld for first update process */
    bs->ld  = 0;
    /* move input generators to basis and generate first spairs */
    update_basis(ps, bs, bht, uht, st, st->ngens, 1);

    if (ps->ld != 2) {
        printf("pload wrong - %d != 2\n", ps->ld);
        return 1;
    }

    if (ps->p[0].gen1 != 0) {
        printf("ps[0].gen1 wrong - %d != 1\n", ps->p[0].gen1);
        return 1;
    }
    if (ps->p[0].gen2 != 1) {
        printf("ps[0].gen2 wrong - %d != 2\n", ps->p[0].gen2);
        return 1;
    }
    if (bht->hd[ps->p[0].lcm].deg != 2) {
        printf("ps[0].lcm.deg wrong - %d != 2\n", bht->hd[ps->p[0].lcm].deg);
        return 1;
    }

    if (ps->p[1].gen1 != 0) {
        printf("ps[1].gen1 wrong - %d != 0\n", ps->p[1].gen1);
        return 1;
    }
    if (ps->p[1].gen2 != 2) {
        printf("ps[1].gen2 wrong - %d != 2\n", ps->p[1].gen2);
        return 1;
    }
    if (bht->hd[ps->p[1].lcm].deg != 3) {
        printf("ps[1].lcm.deg wrong - %d != 3\n", bht->hd[ps->p[1].lcm].deg);
        return 1;
    }

    /* free and clean up */
    free_shared_hash_data(bht);
    free_hash_table(&uht);
    free_hash_table(&bht);
    free_pairset(&ps);
    free_basis(&bs);
    free(st);
    return 0;
}
