#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    /* initialize stuff */
    bs_t *bs  = initialize_basis_ff(4);
    if (bs->ld != 0) {
        return 1;
    }
    if (bs->sz != 8) {
        return 1;
    }
    if (bs->hm == NULL) {
        return 1;
    }
    if (bs->cf_ff == NULL) {
        return 1;
    }
    if (bs->cf_qq != NULL) {
        return 1;
    }
    free_basis(&bs);
    bs  = initialize_basis_qq(4);
    if (bs->ld != 0) {
        return 1;
    }
    if (bs->sz != 8) {
        return 1;
    }
    if (bs->hm == NULL) {
        return 1;
    }
    if (bs->cf_qq == NULL) {
        return 1;
    }
    if (bs->cf_ff != NULL) {
        return 1;
    }
    free_basis(&bs);
    return 0;
}
