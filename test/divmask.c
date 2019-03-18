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

    const int32_t nr_gens     = 3;

    /* initialize stuff */
    stat_t st;
    st.ngens    = 3;
    st.nvars    = 2;
    st.fc       = 101;
    st.init_hts = 12;

    bs_t * bs = initialize_basis_ff(st.ngens);
    ht_t *bht = initialize_basis_hash_table(&st);

    if (bht->ndv != bht->nv) {
        return 1;
    }
    if (bht->bpv != 16) {
        return 1;
    }

    import_julia_data_ff(
            bs, bht, &st, lens, cfs, exps);

    calculate_divmask(bht);
    for (i = 0; i < 16; ++i) {
        if (dm[i] != i+1) {
            return 1;
        }
    }
    for (i = 16; i < 32; ++i) {
        if (dm[i] + 16 != i+1) {
            return 1;
        }
    }

    free_basis(&bs);
    free_hash_table(&bht);
    return 0;
}
