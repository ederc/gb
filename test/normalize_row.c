#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    int32_t i;

    /* set some field characteristic */
    int32_t fc  = 101;

    /* initialize row with 14 entries plus two meta data */
    int32_t *row  = (int32_t *)malloc(16 * sizeof(int32_t));

    row[0]  = 16;
    row[1]  = 8;
    for (i = 1; i < 16; ++i) {
        row[i]  = 300 - i;
    }
    len_t len = 16;
    normalize_dense_matrix_row(row, len, fc);
    for (i=1; i < 16; ++i) {
        printf("%d | %d\n", i, row[i]);
    }
    printf("\n");

    if (row[1] != 25) {
        return 1;
    }
    if (row[2] != 6) {
        return 1;
    }
    if (row[7] != 12) {
        return 1;
    }
    if (row[11] != 37) {
        return 1;
    }
    if (row[15] != 62) {
        return 1;
    }
    free(row);

    return 0;
}
