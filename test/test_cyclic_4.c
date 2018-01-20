#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
    void
    )
{
  int32_t i, j, k;
  len_t *hcm;

  int32_t round = 0;

  const int32_t lens[]  = {4, 4, 4, 2}; 
  const int32_t cfs[]   = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1};
  const int32_t exps[]  = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
    1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1,
    1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0};

  const int32_t nr_vars     = 4;
  const int32_t nr_gens     = 4;
  const int32_t ht_size     = 12;
  const int32_t field_char  = 65521;
  const int32_t nr_threads  = 1;
  const int32_t la_option   = 42;

  int32_t *basis = f4_julia(
      lens, cfs, exps, field_char, nr_vars, nr_gens, ht_size,
      nr_threads, la_option);

  int32_t val[159]  = {159, 7, 20, 15, 20, 30, 20, 20, 25, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 2, 0, 0, 2, 0, 1, 0, 1, 1, 0, 0, 0, 2, 1, 0, 1, 2, 0, 1, 0, 0, 2, 1, 65520, 0, 1, 0, 2, 65520, 0, 0, 0, 3, 1, 0, 1, 1, 2, 1, 0, 0, 2, 2, 65520, 0, 1, 0, 3, 1, 0, 0, 1, 3, 65520, 0, 0, 0, 4, 65520, 0, 0, 0, 0, 1, 0, 0, 3, 2, 1, 0, 0, 2, 3, 65520, 0, 0, 1, 0, 65520, 0, 0, 0, 1, 1, 0, 1, 0, 4, 1, 0, 0, 0, 5, 65520, 0, 1, 0, 0, 65520, 0, 0, 0, 1, 1, 0, 0, 2, 4, 1, 0, 1, 1, 0, 65520, 0, 1, 0, 1, 1, 0, 0, 1, 1, 65519, 0, 0, 0, 2};

  int32_t failure = 0;
  for (i = 0; i < val[0]; ++i) {
    if (val[i] != basis[i]) {
      failure = 1;
      break;
    }
  }
  free(basis);
  basis = NULL;

  return failure;
}
