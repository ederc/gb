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

  const int32_t nr_vars           = 2;
  const int32_t nr_gens           = 3;
  const int32_t ht_size           = 12;
  const int32_t field_char        = 101;
  const int32_t mon_order         = 0;
  const int32_t nr_threads        = 2;
  const int32_t la_option         = 1;
  const int32_t max_nr_pairs      = 10;
  const int32_t reset_hash_table  = 0;

  if (check_and_set_meta_data(lens, cfs, exps, field_char, mon_order,
      nr_vars, nr_gens, ht_size, nr_threads, max_nr_pairs, reset_hash_table,
      la_option)) {
    return 1;
  }

  val_t **mat;
  
  /* initialize stuff */
  initialize_basis(nr_gens);
  initialize_pairset();
  initialize_global_hash_table();
  initialize_local_hash_table();

  mat = import_julia_data(lens, cfs, exps, nr_gens);
  calculate_divmask();
  /* sort initial elements, smallest lead term first */
  qsort(mat, (unsigned long)nrows, sizeof(val_t *),
      matrix_row_initial_input_cmp);

  for (i = 0; i < nrows; ++i) {
    normalize_matrix_row(mat[i]);
  }

  update_basis(mat);

  if (pload != 2) {
    return 1;
  }

  if (ps[0].gen1 != 1) {
    return 1;
  }
  if (ps[0].gen2 != 2) {
    return 1;
  }
  if (ps[0].deg != 2) {
    return 1;
  }
  
  if (ps[1].gen1 != 0) {
    return 1;
  }
  if (ps[1].gen2 != 2) {
    return 1;
  }
  if (ps[1].deg != 3) {
    return 1;
  }

  /* free and clean up */
  free_local_hash_table();
  free_global_hash_table();
  free_pairset();
  free(mat);
  free_basis();
  return 0;
}
