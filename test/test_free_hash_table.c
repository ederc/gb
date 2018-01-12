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

  initialize_global_hash_table();
  initialize_local_hash_table();
  free_local_hash_table();
  if (mlsize != 0) {
    return 1;
  }
  if (elload != 0) {
    return 1;
  }
  if (elsize != 0) {
    return 1;
  }
  if (mapl != NULL) {
    return 1;
  }
  if (evl != NULL) {
    return 1;
  }
  free_global_hash_table();
  if (fc != 0) {
    return 1;
  }
  if (msize != 0) {
    return 1;
  }
  if (eload != 0) {
    return 1;
  }
  if (esize != 0) {
    return 1;
  }
  if (map != NULL) {
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
