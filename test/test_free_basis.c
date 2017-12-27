#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
    void
    )
{
  /* initialize stuff */
  initialize_basis(4);
  free_basis();
  if (bs != NULL) {
    return 1;
  }
  if (bload != 0) {
    return 1;
  }
  if (bsize != 0) {
    return 1;
  }
  return 0;
}
