#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    /* initialize stuff */
    initialize_basis(4);
    free_basis();
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
    return 0;
}
