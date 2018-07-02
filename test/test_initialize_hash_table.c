#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    /* initialize stuff */
    nvars = 34;
    fc    = 101;
    htes  = 10;

    initialize_global_hash_table();
    if (fc != 101) {
        return 1;
    }
    if (nvars != 34) {
        return 1;
    }
    if (ndvars != 32) {
        return 1;
    }
    if (bpv != 1) {
        return 1;
    }
    if (eld != 1) {
        return 1;
    }
    initialize_local_hash_table();
    if (hsz/hlsz != hlsz) {
        return 1;
    }
    if (elld != 1) {
        return 1;
    }
    free_local_hash_table();
    free_global_hash_table();

    return 0;
}
