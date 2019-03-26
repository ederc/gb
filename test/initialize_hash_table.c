#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    /* initialize stuff */
    stat_t st;
    st.nvars    = 34;
    st.fc       = 101;
    st.init_hts = 10;

    ht_t *bht = initialize_basis_hash_table(&st);
    if (bht->nv != 34) {
        return 1;
    }
    if (bht->ndv != 32) {
        return 1;
    }
    if (bht->bpv != 1) {
        return 1;
    }
    if (bht->eld != 1) {
        return 1;
    }
    ht_t *sht = initialize_secondary_hash_table(bht, &st);
    if (bht->hsz/sht->hsz != 32) {
        return 1;
    }
    if (sht->eld != 1) {
        return 1;
    }
    free_shared_hash_data(bht);
    free_hash_table(&sht);
    free_hash_table(&bht);

    return 0;
}
