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
  const int32_t la_option   = 1;

  if (check_and_set_meta_data(lens, cfs, exps, field_char, nr_vars,
      nr_gens, ht_size, nr_threads, la_option)) {
    return 1;
  }

  val_t **mat;
  
  /* initialize stuff */
  initialize_basis(nr_gens);
  initialize_global_hash_table();
  /* initialize_local_hash_table(); */

  mat = import_julia_data(lens, cfs, exps, nr_gens);

  /* for faster divisibility checks, needs to be done after we have
   * read some input data for applying heuristics */
  calculate_divmask();

  /* sort initial elements, smallest lead term first */
  qsort(mat, (unsigned long)nrows, sizeof(val_t *),
      matrix_row_initial_input_cmp);
  /* normalize input generators */
  for (i = 0; i < nrows; ++i) {
    normalize_matrix_row(mat[i]);
  }

    printf("new elements (round %d):\n", round);
    for (i = 0; i < npivs; ++i) {
      for (j = 2; j < mat[i][0]; j += 2) {
        printf("%3d | ", mat[i][j+1]);
        for (k = 0; k < nvars; ++k) {
          printf("%d", (ev+mat[i][j])[k]);
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
    GB_DEBUG(GBDBG, "%3d", round);

    for (int32_t o = 0; o < eload; o += HASH_LEN) {
      if ((ev+o)[HASH_IND] != 0) {
        printf("1 HASH_IND = %d at position %d\n", (ev+o)[HASH_IND], o);
      }
    }
    /* preprocess data for next reduction round */
    mat = select_spairs();
    mat = symbolic_preprocessing(mat);
    printf("before conversion\n");
    for (int32_t l = 0; l < nrows; ++l) {
      for (int32_t m = 0; m < mat[l][0]; ++m) {
        printf("%d ", mat[l][m]);
        if (m > 0 && m % 2 == 0) {
          printf("(");
          for (int32_t o = 0; o < nvars; ++o) {
            printf("%d", (ev+mat[l][m])[o]);
          }
          printf(",%d)  ",(ev+mat[l][m])[HASH_IND]);
        }
      }
      printf("\n");
    }
    /* exponent hashes mapped to column indices for linear algebra */
    hcm = convert_hashes_to_columns(mat);
    printf("after conversion\n");
    for (int32_t l = 0; l < nrows; ++l) {
      for (int32_t m = 0; m < mat[l][0]; ++m) {
        printf("%d ", mat[l][m]);
      }
      printf("\n");
    }
    /* sort matrix rows by decreasing pivots */
    mat = sort_matrix_rows(mat);
    printf("after sorting\n");
    for (int32_t l = 0; l < nrows; ++l) {
      for (int32_t m = 0; m < mat[l][0]; ++m) {
        printf("%d ", mat[l][m]);
      }
      printf("\n");
    }

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
    for (int32_t o = 0; o < eload; o += HASH_LEN) {
      if ((ev+o)[HASH_IND] != 0) {
        printf("2 HASH_IND = %d at position %d\n", (ev+o)[HASH_IND], o);
      }
    }

    printf("new elements (round %d):\n", round);
    for (i = 0; i < npivs; ++i) {
      for (j = 2; j < mat[i][0]; j += 2) {
        printf("%3d | ", mat[i][j+1]);
        for (k = 0; k < nvars; ++k) {
          printf("%d", (ev+mat[i][j])[k]);
        }
        printf(" || ");
      }
      printf("\n");
    }
    update_basis(mat);

    free(mat);
    mat = NULL;

    GB_DEBUG(GBDBG, "\n");
  }

  int32_t *test = export_julia_data();

  int32_t val[159]  = {159, 7, 20, 15, 20, 30, 20, 20, 25, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 2, 0, 0, 2, 0, 1, 0, 1, 1, 0, 0, 0, 2, 1, 0, 1, 2, 0, 1, 0, 0, 2, 1, 65520, 0, 1, 0, 2, 65520, 0, 0, 0, 3, 1, 0, 1, 1, 2, 1, 0, 0, 2, 2, 65520, 0, 1, 0, 3, 1, 0, 0, 1, 3, 65520, 0, 0, 0, 4, 65520, 0, 0, 0, 0, 1, 0, 0, 3, 2, 1, 0, 0, 2, 3, 65520, 0, 0, 1, 0, 65520, 0, 0, 0, 1, 1, 0, 1, 0, 4, 1, 0, 0, 0, 5, 65520, 0, 1, 0, 0, 65520, 0, 0, 0, 1, 1, 0, 0, 2, 4, 1, 0, 1, 1, 0, 65520, 0, 1, 0, 1, 1, 0, 0, 1, 1, 65519, 0, 0, 0, 2};

  int32_t failure = 0;
  for (i = 0; i < val[0]; ++i) {
    if (val[i] != test[i]) {
      failure = 1;
      break;
    }
  }
  free(test);
  test = NULL;

  /* free and clean up */
  /* free_local_hash_table(); */
  free_global_hash_table();
  free_pairset();
  /* note that all rows kept from mat during the overall computation are
   * basis elements and thus we do not need to free the rows itself, but
   * just the matrix structure */
  free_basis();

  return failure;
}
