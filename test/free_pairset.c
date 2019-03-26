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
    free_pairset(&ps);
    if (ps != NULL) {
        return 1;
    }
    return 0;
}
