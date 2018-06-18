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
  if (hlsz != 0) {
    return 1;
  }
  if (elld != 0) {
    return 1;
  }
  if (elsz != 0) {
    return 1;
  }
  if (hmapl != NULL) {
    return 1;
  }
  if (evl != NULL) {
    return 1;
  }
  free_global_hash_table();
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
