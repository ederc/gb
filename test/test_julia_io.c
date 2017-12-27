#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
    void
    )
{
  const int32_t lens[]  = {2,2,2}; 
  const int32_t cfs[]   = {1, 1, 1, 3, 5, 3};
  const int32_t exps[]  = {1, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0};

  const int32_t nr_vars = 2;
  const int32_t nr_gens = 3;
  const int32_t ht_size = 12;

  val_t **mat;
  
  /* initialize stuff */
  initialize_basis(nr_gens);
  initialize_global_hash_table(nr_vars, ht_size);
  if (nvars != nr_vars) {
    return 1;
  }
  initialize_local_hash_table(ht_size);
  if (mlsize != msize) {
    return 1;
  }

  mat = import_julia_data(lens, cfs, exps, nr_gens);

  /* free and clean up */
  free_local_hash_table();
  free_global_hash_table();
  free_basis();
  return 0;
}
