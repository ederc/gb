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
  const int32_t exps[]  = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0};

  const int32_t nr_vars       = 4;
  const int32_t nr_gens       = 4;
  const int32_t ht_size       = 12;
  const int32_t field_char    = 65521;
  const int32_t mon_order     = 1;
  const int32_t nr_threads    = 1;
  const int32_t la_option     = 1;
  const int32_t max_nr_pairs  = 0;

  int32_t failure = 0;

  int32_t **basis = (int32_t **)malloc(sizeof(int32_t *));
  int64_t len     = f4_julia(
      basis, lens, cfs, exps, field_char, mon_order, nr_vars,
      nr_gens, ht_size, nr_threads, max_nr_pairs, la_option);

  if (len != 142) {
    failure = 1;
    return failure;
  }

  int32_t val[142]  = {6, 20, 15, 20, 30, 30, 20, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0,
    1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 2, 0, 0, 2, 0, 1, 0, 1, 1, 0, 0, 0, 2,
    1, 0, 0, 3, 2, 1, 0, 0, 2, 3, 65520, 0, 0, 1, 0, 65520, 0, 0, 0, 1, 1, 0,
    1, 1, 0, 65519, 0, 1, 0, 5, 1, 0, 1, 0, 1, 1, 0, 0, 2, 4, 1, 0, 0, 1, 1,
    65519, 0, 0, 0, 6, 1, 0, 1, 0, 4, 65520, 0, 1, 0, 0, 32760, 0, 0, 2, 7,
    32761, 0, 0, 2, 3, 32762, 0, 0, 0, 5, 32759, 0, 0, 0, 1, 1, 0, 0, 2, 6,
    65520, 0, 0, 2, 2, 65520, 0, 0, 0, 4, 1, 0, 0, 0, 0};

  for (i = 0; i < len; ++i) {
    if (val[i] != (*basis)[i]) {
      failure = 1;
      break;
    }
  }
  free(*basis);
  free(basis);
  basis = NULL;

  return failure;
}
