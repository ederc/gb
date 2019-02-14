#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    /* initialize stuff */
    initialize_basis_ff(4);
    if (bload != 0) {
        return 1;
    }
    if (bsize != 8) {
        return 1;
    }
    free_basis_ff();
    initialize_basis_q(4);
    if (bload != 0) {
        return 1;
    }
    if (bsize != 8) {
        return 1;
    }
    free_basis_q();
    return 0;
}
