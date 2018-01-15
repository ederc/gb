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
  row[1]  = 8;
  for (i = 1; i < 16; ++i) {
    row[i]  = 300 - i;
  }
  for (i=1; i < 16; ++i) {
    printf("%d ", row[i]);
  }
  printf("\n");
  normalize_matrix_row(row);
  for (i=1; i < 16; ++i) {
    printf("%d ", row[i]);
  }
  printf("\n");

  if (row[3] != 1) {
    return 1;
  }
  if (row[5] != 35) {
    return 1;
  }
  if (row[7] != 69) {
    return 1;
  }
  if (row[9] != 2) {
    return 1;
  }
  if (row[11] != 36) {
    return 1;
  }
  if (row[13] != 70) {
    return 1;
  }
  if (row[15] != 3) {
    return 1;
  }
  free(row);

  return 0;
}
