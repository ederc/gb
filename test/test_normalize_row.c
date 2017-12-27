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
  fc  = 101;

  /* initialize row with 14 entries plus two meta data */
  int32_t *row  = (int32_t *)malloc(16 * sizeof(int32_t));

  row[0]  = 16;
  for (i = 1; i < 16; ++i) {
    row[i]  = 300 - i;
  }
  normalize_matrix_row(row);
  if (row[2] != 1) {
    return 1;
  }
  if (row[4] != 62) {
    return 1;
  }
  if (row[6] != 22) {
    return 1;
  }
  if (row[8] != 83) {
    return 1;
  }
  if (row[10] != 43) {
    return 1;
  }
  if (row[12] != 3) {
    return 1;
  }
  if (row[14] != 64) {
    return 1;
  }
  return 0;
}
