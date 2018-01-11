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

  const int32_t lens[]  = {2,2,2}; 
  const int32_t cfs[]   = {1, 1, 1, 3, 5, 3};
  const int32_t exps[]  = {1, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0};

  const int32_t nr_vars     = 2;
  const int32_t nr_gens     = 3;
  const int32_t ht_size     = 12;
  const int32_t field_char  = 101;

  val_t **mat;

  nthrds  = 1;
  laopt   = 1;
  
  /* initialize stuff */
  initialize_basis(nr_gens);
  initialize_global_hash_table(nr_vars, ht_size, field_char);
  if (fc != field_char) {
    return 1;
  }
  if (nvars != nr_vars) {
    return 1;
  }
  initialize_local_hash_table(ht_size);
  if (mlsize != msize) {
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

    printf("new elements (round %d):\n", round);
    for (i = 0; i < npivs; ++i) {
      for (j = 2; j < mat[i][0]; j += 2) {
        printf("%3d | ", mat[i][j+1]);
        for (k = 0; k < nvars; ++k) {
          printf("%d", (evl+mat[i][j])[k]);
        }
        printf(" || ");
      }
      printf("\n");
    }

  /* move input generators to basis and generate first spairs */
  update_basis(mat);

  free(mat);
  mat = NULL;

  /* let's start the f4 rounds,  we are done when no more spairs
   * are left in the pairset */
  for (round = 1; pload > 0; ++round) {

    /* preprocess data for next reduction round */
    mat = select_spairs();
    mat = symbolic_preprocessing(mat);
    /* exponent hashes mapped to column indices for linear algebra */
    hcm = convert_hashes_to_columns(mat);
    /* sort matrix rows by decreasing pivots */
    mat = sort_matrix_rows(mat);

    /* here starts the linear algebra part depending on
     * the chosen options */
    switch (laopt) {
      case 1:
        mat = sparse_linear_algebra(mat);
        break;
      default:
        mat = sparse_linear_algebra(mat);
    }

    mat = convert_columns_to_hashes(mat, hcm);
    
    free(hcm);
    hcm = NULL;

    printf("new elements (round %d):\n", round);
    for (i = 0; i < npivs; ++i) {
      for (j = 2; j < mat[i][0]; j += 2) {
        printf("%3d | ", mat[i][j+1]);
        for (k = 0; k < nvars; ++k) {
          printf("%d", (evl+mat[i][j])[k]);
        }
        printf(" || ");
      }
      printf("\n");
    }
    update_basis(mat);
  }

  /* free and clean up */
  free_local_hash_table();
  free_global_hash_table();
  free_pairset();
  /* note that all rows kept from mat during the overall computation are
   * basis elements and thus we do not need to free the rows itself, but
   * just the matrix structure */
  free(mat);
  free_basis();

  return 0;
}
