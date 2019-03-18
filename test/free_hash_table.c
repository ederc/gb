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
    htes  = 4;

    initialize_basis_hash_table();
    initialize_update_hash_table();
    free_update_hash_table();
    if (husz != 0) {
        return 1;
    }
    if (euld != 0) {
        return 1;
    }
    if (eusz != 0) {
        return 1;
    }
    if (humap != NULL) {
        return 1;
    }
    if (evu != NULL) {
        return 1;
    }
    free_basis_hash_table();
    if (fc != 0) {
        return 1;
    }
    if (hsz != 0) {
        return 1;
    }
    if (eld != 0) {
        return 1;
    }
    if (esz != 0) {
        return 1;
    }
    if (hmap != NULL) {
        return 1;
    }
    if (ev != NULL) {
        return 1;
    }
    if (rv != NULL) {
        return 1;
    }

    return 0;
}
