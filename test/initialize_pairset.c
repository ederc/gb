#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    /* initialize stuff */
    ps_t *ps  = initialize_pairset();
    if (ps == NULL) {
        return 1;
    }
    if (ps->sz == 0) {
        return 1;
    }
    if (ps->ld != 0) {
        return 1;
    }
    free_pairset(&ps);
    return 0;
}
