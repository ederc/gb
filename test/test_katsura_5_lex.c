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

  const int32_t lens[]  = {6, 6, 5, 5, 4}; 
  const int32_t cfs[]   = {1, 2, 2, 2, 2, -1, 1, -1, 2, 2, 2, 2, 2, 2, -1, 2, 2, 2, 1, 2, 2, -1, 2, 2, 2, -1};
  const int32_t exps[]  = {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0};

  const int32_t nr_vars       = 5;
  const int32_t nr_gens       = 5;
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

  free(*basis);
  free(basis);
  basis = NULL;

  return failure;
}
