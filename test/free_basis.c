#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    bs_t * bs = initialize_basis_ff_32(4);
    free_basis(&bs);
    if (bs != NULL) {
        return 1;
    }
    bs  = initialize_basis_qq(4);
    free_basis(&bs);
    if (bs != NULL) {
        return 1;
    }
    return 0;
}
