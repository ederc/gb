#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
    void
    )
{
  /* initialize stuff */
  initialize_global_hash_table(34, 4, 101);
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
  if (eload != HASH_LEN) {
    return 1;
  }
  initialize_local_hash_table(4);
  if (mlsize != msize) {
    return 1;
  }
  if (elload != HASH_LEN) {
    return 1;
  }
  free_local_hash_table();
  free_global_hash_table();

  return 0;
}
