#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
    void
    )
{
  int32_t i, j;

  const int32_t lens[]  = {2,2,2}; 
  const int32_t cfs[]   = {1, 1, 1, 3, 5, 3};
  const int32_t exps[]  = {1, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0};

  const int32_t nr_vars       = 2;
  const int32_t nr_gens       = 3;
  const int32_t ht_size       = 12;
  const int32_t field_char    = 101;
  const int32_t nr_threads    = 2;
  const int32_t la_option     = 1;
  const int32_t max_nr_pairs  = 100;

  if (check_and_set_meta_data(lens, cfs, exps, field_char, nr_vars,
      nr_gens, ht_size, nr_threads, max_nr_pairs, la_option)) {
    return 1;
  }

  val_t **mat;
  
  /* initialize stuff */
  initialize_basis(nr_gens);
  initialize_global_hash_table();
  initialize_local_hash_table();
  if (msize/mlsize != 32) {
    return 1;
  }

  mat = import_julia_data(lens, cfs, exps, nr_gens);

  /* for faster divisibility checks, needs to be done after we have
   * read some input data for applying heuristics */
  calculate_divmask();

  /* normalize input generators */
  for (i = 0; i < nrows; ++i) {
    normalize_matrix_row(mat[i]);
  }

  /* move input generators to basis and generate first spairs */
  update_basis(mat);

  free(mat);
  mat = NULL;

  mat = select_spairs();

  /* free and clean up */
  free_local_hash_table();
  free_global_hash_table();
  /* since we have not moved data from matrix to basis we have to
   * free the matrix rows in this test case */
  for (i = 0; i < nrall; ++i) {
    free(mat[i]);
    mat[i]  = NULL;
  }
  free(mat);
  free_basis();
  return 0;
}
