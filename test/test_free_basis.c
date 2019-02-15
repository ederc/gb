#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    initialize_basis_ff(4);
    free_basis_ff();
    if (gbcf_ff != NULL) {
        return 1;
    }
    if (gbdt != NULL) {
        return 1;
    }
    if (lms!= NULL) {
        return 1;
    }
    if (bload != 0) {
        return 1;
    }
    if (bsize != 0) {
        return 1;
    }
    initialize_basis_q(4);
    free_basis_q();
    if (gbcf_q != NULL) {
        return 1;
    }
    if (gbdt != NULL) {
        return 1;
    }
    if (lms!= NULL) {
        return 1;
    }
    if (bload != 0) {
        return 1;
    }
    if (bsize != 0) {
        return 1;
    }
    return 0;
}
