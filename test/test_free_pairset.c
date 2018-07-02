#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    /* initialize stuff */
    initialize_pairset();
    free_pairset();
    if (ps != NULL) {
        return 1;
    }
    if (pload != 0) {
        return 1;
    }
    if (psize != 0) {
        return 1;
    }
    return 0;
}
