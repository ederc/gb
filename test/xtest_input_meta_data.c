#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
    void
    )
{
  int32_t i;

  const int32_t lens[]  = {2,2,2}; 
  const int32_t cfs[]   = {1, 1, 1, 3, 5, 3};
  const int32_t exps[]  = {1, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0};

  const int32_t nr_vars     = -2;
  const int32_t nr_gens     = -3;
  const int32_t ht_size     = 12;
  const int32_t field_char  = 101;
  const int32_t nr_threads  = 2;
  const int32_t la_option   = 1;

  if (check_and_set_meta_data(lens, cfs, exps, field_char, nr_vars,
      nr_gens, ht_size, nr_threads, la_option)) {
    return 1;
  }

  return 0;
}
