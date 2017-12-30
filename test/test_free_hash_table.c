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
  initialize_local_hash_table(4);
  psize = 120;
  initialize_lcm_hash_table();
  free_lcm_hash_table();
  if (mssize != 0) {
    return 1;
  }
  if (esload != 0) {
    return 1;
  }
  if (essize != 0) {
    return 1;
  }
  if (maps != NULL) {
    return 1;
  }
  if (evs != NULL) {
    return 1;
  }
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
