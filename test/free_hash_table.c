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
    ht_t *sht = initialize_secondary_hash_table(bht, &st);
    free_hash_table(&sht);
    if (sht != NULL) {
        return 1;
    }
    free_hash_table(&bht);
    if (bht != NULL) {
        return 1;
    }
    return 0;
}
