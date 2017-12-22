#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/data.h"

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

  f4_julia(lens, cfs, exps, nr_vars, nr_gens, ht_size);
  return 0;
}
