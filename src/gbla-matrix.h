/* gbla: Gröbner Basis Linear Algebra
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
 * This file is part of gbla.
 * gbla is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * gbla is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with gbla . If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * \file matrix.h
 * \brief Interfaces for matrix structures
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GBLA_MATRIX_H
#define GBLA_MATRIX_H

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "../config.h"
#include "types.h"
/* #if GBLA_WITH_FFLAS
 * #include "../draft/matrix.h"
 * #endif */

/**
 * \brief returns number of row blocks for intermediate block submatrix A
 *
 * \param intermediate block submatrix A
 *
 * \return number of row blocks for A
 */
static inline ri_t get_number_intermediate_row_blocks(const ibm_fl_t *A)
{
  return (ri_t) ceil((float) A->nrows / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of column blocks for intermediate block submatrix A
 *
 * \param intermediate block submatrix A
 *
 * \return number of column blocks for A
 */
static inline ci_t get_number_intermediate_col_blocks(const ibm_fl_t *A)
{
  return (ci_t) ceil((float) A->ncols / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of row blocks for hybrid block submatrix A
 *
 * \param hybrid block submatrix A
 *
 * \return number of row blocks for A
 */
static inline ri_t get_number_hybrid_row_blocks(const hbm_fl_t *A)
{
  return (ri_t) ceil((float) A->nrows / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of column blocks for hybrid block submatrix A
 *
 * \param hybrid block submatrix A
 *
 * \return number of column blocks for A
 */
static inline ci_t get_number_hybrid_col_blocks(const hbm_fl_t *A)
{
  return (ci_t) ceil((float) A->ncols / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of row blocks for dense block submatrix A
 *
 * \param dense block submatrix A
 *
 * \return number of row blocks for A
 */
static inline ri_t get_number_dense_row_blocks(const dbm_fl_t *A)
{
  return (ri_t) ceil((float) A->nrows / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of column blocks for dense block submatrix A
 *
 * \param dense block submatrix A
 *
 * \return number of column blocks for A
 */
static inline ci_t get_number_dense_col_blocks(const dbm_fl_t *A)
{
  return (ci_t) ceil((float) A->ncols / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of row blocks for sparse block submatrix A
 *
 * \param sparse block submatrix A
 *
 * \return number of row blocks for A
 */
static inline ri_t get_number_sparse_row_blocks(const sb_fl_t *A)
{
  return (ri_t) ceil((float) A->nrows / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief returns number of column blocks for sparse block submatrix A
 *
 * \param sparse block submatrix A
 *
 * \return number of column blocks for A
 */
static inline ci_t get_number_sparse_col_blocks(const sb_fl_t *A)
{
  return (ci_t) ceil((float) A->ncols / __GBLA_SIMD_BLOCK_SIZE);
}

/**
 * \brief Initializes sparse submatrices without preallocating memory for rows
 *
 * \param sparse submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
static inline void init_sm_no_row_allocation(sm_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  ri_t i;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for rows
  A->row  = (re_t **)malloc(nrows * sizeof(re_t *));
  A->pos  = (ci_t **)malloc(nrows * sizeof(ci_t *));
  A->sz   = (ci_t *)malloc(nrows * sizeof(ci_t));
  A->buf  = (ci_t *)malloc(nrows * sizeof(ci_t));
  memset(A->sz, 0, nrows * sizeof(ci_t));
  // allocate memory already, we will have at least one entry in each row!
  for (i=0; i<nrows; ++i) {
    A->row[i] = NULL;
    A->pos[i] = NULL;
    A->buf[i] = 0;
  }
}

/**
 * \brief Initializes sparse submatrices
 *
 * \param sparse submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
static inline void init_sm(sm_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  ri_t i;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for rows
  A->row  = (re_t **)malloc(nrows * sizeof(re_t *));
  A->pos  = (ci_t **)malloc(nrows * sizeof(ci_t *));
  A->sz   = (ci_t *)malloc(nrows * sizeof(ci_t));
  A->buf  = (ci_t *)malloc(nrows * sizeof(ci_t));
  memset(A->sz, 0, nrows * sizeof(ci_t));
  // allocate memory already, we will have at least one entry in each row!
  for (i=0; i<nrows; ++i) {
    A->row[i] = (re_t *)malloc(2 * __GBLA_SIMD_BLOCK_SIZE * sizeof(re_t));
    A->pos[i] = (ci_t *)malloc(2 * __GBLA_SIMD_BLOCK_SIZE * sizeof(ci_t));
    A->buf[i] = 2 * __GBLA_SIMD_BLOCK_SIZE;
  }
}

/**
 * \brief Initializes dense row submatrices
 *
 * \param dense row submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
static inline void init_dm(dm_t *A, const ri_t nrows, const ri_t ncols)
{
  ri_t i;
  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->rank   = nrows;  // rank assumption
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for matrix
  A->row  = (dr_t **)malloc(nrows * sizeof(dr_t *));
  for (i=0; i<nrows; ++i) {
    A->row[i]           = (dr_t *)malloc(sizeof(dr_t));
    A->row[i]->init_val = (re_t *)malloc(ncols * sizeof(re_t));
    A->row[i]->val      = NULL;
    A->row[i]->piv_val  = NULL;
  }
}

/**
 * \brief Initializes sparse block submatrices
 *
 * \param sparse block submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
static inline void init_sb(sb_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  ri_t i, j;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for blocks

  // row and column loops
  const ci_t clA  = get_number_sparse_col_blocks(A);
  const ri_t rlA  = get_number_sparse_row_blocks(A);
  // we know that no block of A will be empty, thus we can already allocate
  // corresponding memory
  A->blocks = (sbl_t **)malloc(rlA * sizeof(sbl_t *));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (sbl_t *)malloc(clA * sizeof(sbl_t));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j].val = NULL;
      A->blocks[i][j].pos = NULL;
      A->blocks[i][j].sz  = NULL;
      A->blocks[i][j].buf = NULL;
      /*
      A->blocks[i][j].val = (re_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
      A->blocks[i][j].pos = (bi_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
      A->blocks[i][j].sz  = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
      A->blocks[i][j].buf = (bi_t *)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
        A->blocks[i][j].row[k]  = NULL;
        A->blocks[i][j].pos[k]  = NULL;
        A->blocks[i][j].sz[k]   = 0;
        A->blocks[i][j].buf[k]  = 0;
      }
      */
    }
  }
}

/**
 * \brief Initializes hybrid block submatrices
 *
 * \param hybrid block submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
static inline void init_hbm(hbm_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  ri_t i, j;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for blocks

  // row and column loops
  const ci_t clA  = get_number_hybrid_col_blocks(A);
  const ri_t rlA  = get_number_hybrid_row_blocks(A);

  A->blocks = (dbl_t ****)malloc(rlA * sizeof(dbl_t ***));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (dbl_t ***)malloc(clA * sizeof(dbl_t **));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j] = NULL;
    }
  }
}

/**
 * \brief Initializes intermediate block submatrices
 *
 * \param intermediate block submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
static inline void init_ibm(ibm_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  ri_t i, j;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for blocks

  // row and column loops
  const ci_t clA  = get_number_intermediate_col_blocks(A);
  const ri_t rlA  = get_number_intermediate_row_blocks(A);

  A->blocks = (ibl_t **)malloc(rlA * sizeof(ibl_t *));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (ibl_t *)malloc(clA * sizeof(ibl_t));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j].val = NULL;
      A->blocks[i][j].row = NULL;
      A->blocks[i][j].pos = NULL;
      A->blocks[i][j].sz  = NULL;
    }
  }
}

/**
 * \brief Initializes dense block submatrices
 *
 * \param dense block submatrix A
 *
 * \param number of rows nrows
 *
 * \param number of columns ncols
 */
static inline void init_dbm(dbm_fl_t *A, const ri_t nrows, const ri_t ncols)
{
  ri_t i, j;

  // initialize meta data for block submatrices
  A->nrows  = nrows;  // row dimension
  A->ncols  = ncols;  // col dimension
  A->nnz    = 0;      // number nonzero elements

  // allocate memory for blocks

  // row and column loops
  const ci_t clA  = get_number_dense_col_blocks(A);
  const ri_t rlA  = get_number_dense_row_blocks(A);

  A->blocks = (dbl_t **)malloc(rlA * sizeof(dbl_t *));
  for (i=0; i<rlA; ++i) {
    A->blocks[i]  = (dbl_t *)malloc(clA * sizeof(dbl_t));
    for (j=0; j<clA; ++j) {
      A->blocks[i][j].val = NULL;
    }
  }
}

/**
 * \brief Frees given sparse submatrix A
 *
 * \param sparse submatrix A
 *
 * \param number of threads for parallel computation
 */
static inline void free_sparse_matrix(sm_fl_t **A_in, int nthrds)
{
  sm_fl_t *A      = *A_in;
  ri_t i;
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i)
    for (i=0; i<A->nrows; ++i) {
      free(A->row[i]);
      free(A->pos[i]);
    }
  }
  free(A->row);
  free(A->pos);
  free(A->sz);
  free(A->buf);
  free(A);
  A = NULL;
  *A_in  = A;
}

/**
 * \brief Frees given sparse block submatrix A
 *
 * \param sparse block submatrix A
 *
 * \param number of threads for parallel computation
 */
static inline void free_sparse_submatrix(sb_fl_t **A_in, int nthrds)
{
  sb_fl_t *A      = *A_in;
  const ci_t clA  = get_number_sparse_col_blocks(A);
  const ri_t rlA  = get_number_sparse_row_blocks(A);
  ri_t i, j, k;
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i, j, k)
    for (i=0; i<rlA; ++i) {
      for (j=0; j<clA; ++j) {
        if (A->blocks[i][j].val != NULL) {
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
            free(A->blocks[i][j].val[k]);
            free(A->blocks[i][j].pos[k]);
          }
          free(A->blocks[i][j].val);
          free(A->blocks[i][j].pos);
          free(A->blocks[i][j].sz);
          free(A->blocks[i][j].buf);
        }
      }
      free(A->blocks[i]);
    }
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
}

/**
 * \brief Frees given hybrid block submatrix A
 *
 * \param hybrid block submatrix A
 *
 * \param number of threads for parallel computation
 */
static inline void free_hybrid_submatrix(hbm_fl_t **A_in, int nthrds)
{
  hbm_fl_t *A     = *A_in;
  const ci_t clA  = get_number_hybrid_col_blocks(A);
  const ri_t rlA  = get_number_hybrid_row_blocks(A);
  ri_t i, j, k, l;
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i, j, k, l)
    for (i=0; i<rlA; ++i) {
      for (j=0; j<clA; ++j) {
        if (A->blocks[i][j] != NULL) {
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
            if (A->blocks[i][j][k] != NULL) {
              for (l=0; l<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++l) {
                if (A->blocks[i][j][k][l].val != NULL) {
                  free(A->blocks[i][j][k][l].val);
                  A->blocks[i][j][k][l].val  = NULL;
                }
              }
              free(A->blocks[i][j][k]);
              A->blocks[i][j][k]  = NULL;
            }
          }
          free(A->blocks[i][j]);
          A->blocks[i][j] = NULL;
        }
      }
      free(A->blocks[i]);
      A->blocks[i]  = NULL;
    }
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
}

void realloc_block_rows_B(sbm_fl_t *A, const ri_t rbi, const ci_t bir,
    const bi_t lib, const bi_t init_buffer_A, bi_t *buffer_A);

/**
 * \brief Frees given intermediate block submatrix A
 *
 * \param intermediate block submatrix A
 */
static inline void free_intermediate_submatrix(ibm_fl_t **A_in)
{
  ibm_fl_t *A     = *A_in;
  const ci_t clA  = get_number_intermediate_col_blocks(A);
  const ri_t rlA  = get_number_intermediate_row_blocks(A);
  ri_t j;
  ci_t i;
  bi_t k;
  // free A
  for (j=0; j<rlA; ++j) {
    for (i=0; i<clA; ++i) {
      if (A->blocks[j][i].val != NULL) {
        free(A->blocks[j][i].val);
        A->blocks[j][i].val  = NULL;
      }
      if (A->blocks[j][i].row != NULL) {
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
          free(A->blocks[j][i].row[k]);
          free(A->blocks[j][i].pos[k]);
        }
        free(A->blocks[j][i].row);
        free(A->blocks[j][i].pos);
        free(A->blocks[j][i].sz);
      }
    }
    free(A->blocks[j]);
    A->blocks[j]  = NULL;
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
}

/**
 * \brief Frees given dense block submatrix A
 *
 * \param dense block submatrix A
 *
 * \param number of threads for parallel computation
 */
static inline void free_dense_submatrix(dbm_fl_t **A_in, int nthrds)
{
  dbm_fl_t *A     = *A_in;
  const ci_t clA  = get_number_dense_col_blocks(A);
  const ri_t rlA  = get_number_dense_row_blocks(A);
  ri_t j;
  ci_t i;
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(i,j)
    for (j=0; j<rlA; ++j) {
      for (i=0; i<clA; ++i) {
        if (A->blocks[j][i].val != NULL) {
          free(A->blocks[j][i].val);
          A->blocks[j][i].val  = NULL;
        }
      }
      free(A->blocks[j]);
      A->blocks[j]  = NULL;
    }
  }
  free(A->blocks);
  A->blocks = NULL;
  free(A);
  A = NULL;
  *A_in  = A;
}

/**
 * \brief Frees given dense row submatrix A
 *
 * \param dense row submatrix A
 *
 * \param number of threads for parallel computation
 */
static inline void free_dense_row_submatrix(dm_t **A_in, int nthrds)
{
  dm_t *A   = *A_in;
  ri_t j;
  // free A
#pragma omp parallel num_threads(nthrds)
  {
#pragma omp for private(j)
    for (j=0; j<A->nrows; ++j) {
      free(A->row[j]->piv_val);
      free(A->row[j]->init_val);
      free(A->row[j]->val);
      free(A->row[j]);
    }
  }
  free(A->row);
  free(A);
  A = NULL;
  *A_in  = A;
}

/**
 * \brief Comparison function for qsort for dense pivot rows
 *
 * \param dense row a
 *
 * \param dense row b
 *
 * \return b.piv_lead - a.piv_lead
 */
static inline int cmp_dr_pivs(const void *a, const void *b)
{
  if ((dr_t *)a == NULL)
    return -1;
  if ((dr_t *)b == NULL)
    return 1;
  return (int)((*((dr_t **)(a)))->piv_lead) - (int)((*((dr_t **)(b)))->piv_lead);
}

/**
 * \brief Sorts dense row matrix A by its pivot leads: A->piv_lead[i] <= A->piv_lead[j] for
 * i <= j.
 *
 * \param dense row submatrix A
 *
 * \param upto which row to be sorted nrows
 */
static inline void sort_partly_reduced_matrix_by_pivots(dm_t *A, ri_t nrows)
{
  qsort(A->row, nrows, sizeof(dr_t **), cmp_dr_pivs);
  //return A->row[0]->piv_lead;
}

/**
 * \brief Comparison function for qsort for dense rows
 *
 * \param dense row a
 *
 * \param dense row b
 *
 * \return b.lead - a.lead
 */
static inline int cmp_dr(const void *a, const void *b)
{
  if ((dr_t *)a == NULL)
    return -1;
  if ((dr_t *)b == NULL)
    return 1;
  return (int)((*((dr_t **)(a)))->lead) - (int)((*((dr_t **)(b)))->lead);
}

/**
 * \brief Sorts dense row matrix A by its pivots: A->lead[i] <= A->lead[j] for
 * i <= j.
 *
 * \param dense row submatrix A
 *
 * \return column index of first row in A after sorting
 */
static inline ci_t sort_dense_matrix_by_pivots(dm_t *A)
{
  qsort(A->row, A->rank, sizeof(dr_t **), cmp_dr);
  return A->row[0]->lead;
}

/**
 * \brief Takes a captured zero row with index ridx from A and moves it to the
 * first zero row position. If there are nonzero rows below row ridx then it
 * swaps the bottom one of these with row ridx.
 *
 * \param dense row matrix A
 *
 * \param row index ridx
 *
 * \return index of row it swapped row ridx with (this is possibly ridx itself).
 */
static inline ri_t swap_zero_row(dm_t *A, const ri_t ridx)
{
  const ri_t swp_idx  = A->rank - 1;
  if (ridx < swp_idx) {
    dr_t *tmp = A->row[swp_idx];
    A->row[swp_idx] = A->row[ridx];
    A->row[ridx]    = tmp;
  }
  A->rank--;

  return A->rank;
}

/**
 * \brief Takes a captured row with index ridx from A and moves it to the
 * first zero row position. If there are nonzero rows below row ridx then it
 * swaps the bottom one of these with row ridx.
 *
 * \param dense row matrix A
 *
 * \param row index ridx
 *
 * \return index of row it swapped row ridx with (this is possibly ridx itself).
 */
static inline ri_t swap_row(dm_t *A, const ri_t ridx)
{
  const ri_t swp_idx  = A->rank - 1;
  if (ridx < swp_idx) {
    dr_t *tmp = A->row[swp_idx];
    A->row[swp_idx] = A->row[ridx];
    A->row[ridx]    = tmp;
  }

  return A->rank;
}

/**
 * \brief Copies block matrix D to a dense row matrix
 *
 * \note Zero rows are directly removed, rows are sorted by known pivots if
 * parameter sort is 1.
 *
 * \param block submatrix A
 *
 * \param number of threads for parallel computation
 *
 * \param sorting of dense matrix by pivots? sort
 *
 * \return dense row matrix
 */
static inline dm_t *copy_block_to_dense_matrix(dbm_fl_t **A,
    const int nthrds, const int sort)
{
  ri_t i, j;
  dbm_fl_t *in  = *A;

  dm_t *out  = (dm_t *)malloc(sizeof(dm_t));
  init_dm(out, in->nrows, in->ncols);

  const ri_t blocks_per_col = get_number_dense_row_blocks(in);
  const ci_t blocks_per_row = get_number_dense_col_blocks(in);

  // set entries in out to zero
  for (i=0; i<out->nrows; ++i)
    memset(out->row[i]->init_val, 0, out->ncols * sizeof(re_t));

#pragma omp parallel num_threads(nthrds)
  {
    ri_t i, min_k;
    ci_t j, min_l;
    bi_t k;
#pragma omp for
    for (i=0; i<blocks_per_col; ++i) {
      min_k = __GBLA_SIMD_BLOCK_SIZE < (out->nrows-i*__GBLA_SIMD_BLOCK_SIZE) ?
        __GBLA_SIMD_BLOCK_SIZE : (out->nrows-i*__GBLA_SIMD_BLOCK_SIZE);
      for (j=0; j<blocks_per_row; ++j) {
        if (in->blocks[i][j].val != NULL) {
          min_l = __GBLA_SIMD_BLOCK_SIZE < (out->ncols-j*__GBLA_SIMD_BLOCK_SIZE) ?
            __GBLA_SIMD_BLOCK_SIZE : (out->ncols-j*__GBLA_SIMD_BLOCK_SIZE);
          for (k=0; k<min_k; ++k) {
              memcpy(
                  out->row[i*__GBLA_SIMD_BLOCK_SIZE+k]->init_val+(j*__GBLA_SIMD_BLOCK_SIZE),
                in->blocks[i][j].val+(k*__GBLA_SIMD_BLOCK_SIZE),
                min_l * sizeof(re_t));
          }
          free(in->blocks[i][j].val);
        }
      }
      free(in->blocks[i]);
    }
  }
  free(in->blocks);
  free(in);
  in  = NULL;
  *A  = in;

  // get first nonzero entry in each row
  i = 0;
next_round:
  for (; i<out->nrows; ++i) {
    for (j=0; j<out->ncols; ++j) {
      if (out->row[i]->init_val[j] != 0) {
        out->row[i]->lead  = j;
        i++;
        goto next_round;
      }
    }
    // if we found a zero row, i.e. we have not broken the j-loop beforehand
    out->row[i]->lead = out->ncols;
    free(out->row[i]->init_val);
    out->row[i]->init_val = NULL;
  }

  // now sort out and get rank
  if (sort)
    sort_dense_matrix_by_pivots(out);
  i = out->nrows-1;
  while (i < (re_t)(-1) && out->row[i]->lead == out->ncols) {
    out->rank--;
    i--;
  }
  /*
     for (int ii=0; ii<out->rank; ++ii)
     printf("%u --------------------> %u\n",ii, out->row[ii]->lead);
     */
  return out;
}

/**
 * \brief Updates the lead index of the dense row 
 *
 * \param dense row matrix A
 *
 * \param index of row ridx
 */
static inline void update_lead_of_row(dm_t *A, const ri_t ridx)
{
  ci_t i;
  register re_l_t *values = A->row[ridx]->val;
  for (i=A->row[ridx]->lead; i<A->ncols; ++i) {
    if (values[i] != 0) {
      values[i] = MODP(values[i], A->mod);
      if (values[i] != 0) {
        A->row[ridx]->lead  = i;
        return;
      }
    }
  }
  A->row[ridx]->lead  = A->ncols;
}

/**
 * \brief Copies sparse matrix A to a sparse block matrix
 *
 * \param sparse submatrix A
 *
 * \param number of threads for parallel computation
 *
 * \return sparse block matrix
 */
static inline sb_fl_t *copy_sparse_to_block_matrix(const sm_fl_t *A,
    const int nthrds)
{
  sb_fl_t *out  = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  init_sb(out, A->nrows, A->ncols);

  const ri_t blocks_per_col = get_number_sparse_row_blocks(out);
  const ci_t blocks_per_row = get_number_sparse_col_blocks(out);

#pragma omp parallel num_threads(nthrds)
  {
    bi_t *block_length  = (bi_t *)calloc(blocks_per_row, sizeof(bi_t));
    bi_t col_idx;
    ri_t i, j;
    ci_t k, l, m;
    ci_t ctr;
    ri_t max;
#pragma omp for 
    for (i=0; i<blocks_per_col; ++i) {
      max = A->nrows < (i+1)*__GBLA_SIMD_BLOCK_SIZE ? A->nrows : (i+1)*__GBLA_SIMD_BLOCK_SIZE;
      for (j=i*__GBLA_SIMD_BLOCK_SIZE; j<max; ++j) {
        memset(block_length, 0, blocks_per_row * sizeof(bi_t));
        // get lengths for the sparse block rows
        for (k=0; k<A->sz[j]; ++k) {
          block_length[(A->ncols - 1 - A->pos[j][k]) / __GBLA_SIMD_BLOCK_SIZE]++;
        }

        col_idx = j % __GBLA_SIMD_BLOCK_SIZE;
        // allocate memory
        for (l=0; l<blocks_per_row; ++l) {
          if (block_length[l] != 0) {
            if (out->blocks[i][l].val == NULL) {
              out->blocks[i][l].val = (re_t **)malloc(
                  __GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
              out->blocks[i][l].pos = (bi_t **)malloc(
                  __GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
              out->blocks[i][l].sz = (bi_t *)malloc(
                  __GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
              for (m=0; m<__GBLA_SIMD_BLOCK_SIZE; ++m) {
                out->blocks[i][l].val[m]  = NULL;
                out->blocks[i][l].pos[m]  = NULL;
                out->blocks[i][l].sz[m]   = 0;
              }
            }
          // memory for row j % __GBLA_SIMD_BLOCK_SIZE
          out->blocks[i][l].val[col_idx]  = (re_t *)malloc(
              block_length[l] * sizeof(re_t));
          out->blocks[i][l].pos[col_idx]  = (bi_t *)malloc(
              block_length[l] * sizeof(bi_t));
          }
        }

        // copy elements
        // we loop at "+1" in order to not get into trouble with comparing
        // unsigned ints with 0 at the very end. thus we decrement "m" to "m-1"
        // in each loop.
        ctr = A->sz[j];
        for (l=0; l<blocks_per_row; ++l) {
          for (m=ctr; m>(ctr-block_length[l]); --m) {
            out->blocks[i][l].val[col_idx][out->blocks[i][l].sz[col_idx]] =
              A->row[j][m-1];
            out->blocks[i][l].pos[col_idx][out->blocks[i][l].sz[col_idx]] =
              (bi_t)((A->ncols - 1 -A->pos[j][m-1]) % __GBLA_SIMD_BLOCK_SIZE);
            out->blocks[i][l].sz[col_idx]++;
          }
          ctr -=  block_length[l];
        }
      }
    }    
    free(block_length);
  }
  return out;
}

/**
 * \brief Copies sparse matrix A to a intermediate block matrix
 *
 * \param sparse submatrix A
 *
 * \param number of threads for parallel computation
 *
 * \return intermediate block matrix
 */
static inline ibm_fl_t *copy_sparse_to_intermediate_block_matrix(const sm_fl_t *A,
    const int nthrds)
{
  ibm_fl_t *out  = (ibm_fl_t *)malloc(sizeof(ibm_fl_t));
  init_ibm(out, A->nrows, A->ncols);

  const ri_t blocks_per_col = get_number_intermediate_row_blocks(out);
  const ci_t blocks_per_row = get_number_intermediate_col_blocks(out);

#pragma omp parallel num_threads(nthrds)
  {
    bi_t *block_length  = (bi_t *)calloc(blocks_per_row, sizeof(bi_t));
    bi_t col_idx;
    ri_t i, j;
    ci_t k, l, m;
    ci_t ctr;
    ri_t max;
#pragma omp for   
    for (i=0; i<blocks_per_col; ++i) {
      max = A->nrows < (i+1)*__GBLA_SIMD_BLOCK_SIZE ? A->nrows : (i+1)*__GBLA_SIMD_BLOCK_SIZE;
      for (j=i*__GBLA_SIMD_BLOCK_SIZE; j<max; ++j) {
        memset(block_length, 0, blocks_per_row * sizeof(bi_t));
        // get lengths for the sparse block rows
        /*
        tmp2  = (A->ncols - 1 - A->pos[j][0]) / __GBLA_SIMD_BLOCK_SIZE;
        for (k=1; k<A->sz[j]; ++k) {
          tmp1  = tmp2;
          tmp2  = (A->ncols - 1 - A->pos[j][k]) / __GBLA_SIMD_BLOCK_SIZE;
          printf("tmp1 %d -- tmp2 %d\n",tmp1,tmp2);
          if (tmp1 == tmp2)
            block_length[tmp1]++;
        }
        if (tmp1 == tmp2)
          block_length[tmp1]++;
        else
          block_length[tmp2]++;
        */
        for (k=0; k<A->sz[j]; ++k) {
          block_length[(A->ncols - 1 - A->pos[j][k]) / __GBLA_SIMD_BLOCK_SIZE]++;
        }

        col_idx = j % __GBLA_SIMD_BLOCK_SIZE;
        // allocate memory
        for (l=0; l<blocks_per_row; ++l) {
          if (block_length[l] != 0) {
            if (out->blocks[i][l].row == NULL) {
              out->blocks[i][l].row = (re_t **)malloc(
                  __GBLA_SIMD_BLOCK_SIZE * sizeof(re_t *));
              out->blocks[i][l].pos = (bi_t **)malloc(
                  __GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t *));
              out->blocks[i][l].sz = (bi_t *)malloc(
                  __GBLA_SIMD_BLOCK_SIZE * sizeof(bi_t));
              for (m=0; m<__GBLA_SIMD_BLOCK_SIZE; ++m) {
                out->blocks[i][l].row[m]  = NULL;
                out->blocks[i][l].pos[m]  = NULL;
                out->blocks[i][l].sz[m]   = 0;
              }
            }
          // memory for row j % __GBLA_SIMD_BLOCK_SIZE
          out->blocks[i][l].row[col_idx]  = (re_t *)malloc(
              block_length[l] * sizeof(re_t));
          out->blocks[i][l].pos[col_idx]  = (bi_t *)malloc(
              block_length[l] * sizeof(bi_t));
          }
        }

        // copy elements
        // we loop at "+1" in order to not get into trouble with comparing
        // unsigned ints with 0 at the very end. thus we decrement "m" to "m-1"
        // in each loop.
        ctr = A->sz[j];
        for (l=0; l<blocks_per_row; ++l) {
          for (m=ctr; m>(ctr-block_length[l]); --m) {
            out->blocks[i][l].row[col_idx][out->blocks[i][l].sz[col_idx]] =
              A->row[j][m-1];
            out->blocks[i][l].pos[col_idx][out->blocks[i][l].sz[col_idx]] =
              (bi_t)((A->ncols - 1 -A->pos[j][m-1]) % __GBLA_SIMD_BLOCK_SIZE);
            out->blocks[i][l].sz[col_idx]++;
          }
          ctr -=  block_length[l];
        }
      }
    }    
    free(block_length);
  }
#pragma omp parallel num_threads(nthrds)
  {
    ri_t i, j;
    ci_t k, l;
    unsigned long nnz;
#pragma omp for   
  // now switch blocks that are dense to dense representation
  for (i=0; i<blocks_per_col; ++i) {
    for (j=0; j<blocks_per_row; ++j) {
      if (out->blocks[i][j].sz != NULL) {
        nnz = 0;
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
          nnz +=  out->blocks[i][j].sz[k];
        }
        if (nnz > __GBLA_SIMD_BLOCK_SIZE_RECT/2) {
          out->blocks[i][j].val = (re_t *)calloc(
              __GBLA_SIMD_BLOCK_SIZE_RECT, sizeof(re_t));
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
            //printf("sz[%u][%u][%u]\n",i,j,k);
            //printf("= %p\n",out->blocks[i][j].sz);
            for (l=0; l<out->blocks[i][j].sz[k]; ++l) {
              out->blocks[i][j].val[k * __GBLA_SIMD_BLOCK_SIZE + out->blocks[i][j].pos[k][l]] =
                out->blocks[i][j].row[k][l];
            }
            free(out->blocks[i][j].row[k]);
            free(out->blocks[i][j].pos[k]);
            //free(out->blocks[i][j].sz);
            //out->blocks[i][j].row[k]  = NULL;
            //out->blocks[i][j].pos[k]  = NULL;
          }
          free(out->blocks[i][j].row);
          free(out->blocks[i][j].pos);
          free(out->blocks[i][j].sz);
          out->blocks[i][j].row = NULL;
          out->blocks[i][j].pos = NULL;
          out->blocks[i][j].sz  = NULL;
        }
      }
    }
  }
}
  return out;
}

/**
 * \brief Copies sparse matrix A to a dense block matrix
 *
 * \param sparse submatrix A
 *
 * \param number of threads for parallel computation
 *
 * \return dense block matrix
 */
static inline dbm_fl_t *copy_sparse_to_dense_block_matrix(const sm_fl_t *A,
    const int nthrds)
{
  dbm_fl_t *out  = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  init_dbm(out, A->nrows, A->ncols);

  const ri_t blocks_per_col = get_number_dense_row_blocks(out);
  const ci_t blocks_per_row = get_number_dense_col_blocks(out);

#pragma omp parallel num_threads(nthrds)
  {
    bi_t *block_length  = (bi_t *)calloc(blocks_per_row, sizeof(bi_t));
    bi_t col_idx;
    ri_t i, j;
    ci_t k, l, m;
    ci_t ctr;
    ri_t max;
#pragma omp for   
    for (i=0; i<blocks_per_col; ++i) {
      max = A->nrows < (i+1)*__GBLA_SIMD_BLOCK_SIZE ? A->nrows : (i+1)*__GBLA_SIMD_BLOCK_SIZE;
      for (j=i*__GBLA_SIMD_BLOCK_SIZE; j<max; ++j) {
        memset(block_length, 0, blocks_per_row * sizeof(bi_t));
        // get lengths for the sparse block rows
        /*
        tmp2  = (A->ncols - 1 - A->pos[j][0]) / __GBLA_SIMD_BLOCK_SIZE;
        for (k=1; k<A->sz[j]; ++k) {
          tmp1  = tmp2;
          tmp2  = (A->ncols - 1 - A->pos[j][k]) / __GBLA_SIMD_BLOCK_SIZE;
          printf("tmp1 %d -- tmp2 %d\n",tmp1,tmp2);
          if (tmp1 == tmp2)
            block_length[tmp1]++;
        }
        if (tmp1 == tmp2)
          block_length[tmp1]++;
        else
          block_length[tmp2]++;
        */
        for (k=0; k<A->sz[j]; ++k) {
          block_length[(A->ncols - 1 - A->pos[j][k]) / __GBLA_SIMD_BLOCK_SIZE]++;
        }

        col_idx = j % __GBLA_SIMD_BLOCK_SIZE;
        // allocate memory
        for (l=0; l<blocks_per_row; ++l) {
          if (block_length[l] != 0) {
            if (out->blocks[i][l].val == NULL) {
              out->blocks[i][l].val = (re_t *)calloc(
                  __GBLA_SIMD_BLOCK_SIZE_RECT, sizeof(re_t));
            }
          }
        }

        // copy elements
        // we loop at "+1" in order to not get into trouble with comparing
        // unsigned ints with 0 at the very end. thus we decrement "m" to "m-1"
        // in each loop.
        ctr = A->sz[j];
        for (l=0; l<blocks_per_row; ++l) {
          for (m=ctr; m>(ctr-block_length[l]); --m) {
            out->blocks[i][l].val[(col_idx * __GBLA_SIMD_BLOCK_SIZE)+
            ((A->ncols - 1 -A->pos[j][m-1]) % __GBLA_SIMD_BLOCK_SIZE)] =
              A->row[j][m-1];
          }
          ctr -=  block_length[l];
        }
      }
    }    
  }
  /*
  printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  for (int ii=0; ii<blocks_per_col; ++ii) {
    for (int jj=0; jj<blocks_per_row; ++jj) {
      if (out->blocks[ii][jj].row != NULL) {
        printf("%d .. %d\n", ii, jj);
        for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
          for (int ll=0; ll<out->blocks[ii][jj].sz[kk]; ++ll) {
            printf("%d | %d || ", out->blocks[ii][jj].row[kk][ll], out->blocks[ii][jj].pos[kk][ll]);
          }
          printf("\n");
        }
      }
    }
  }
  printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  */
  return out;
}

/**
 * \brief Copies data from block matrix in to input sparse matrix format. If deleteIn
 * is set the input block matrix in is deleted.
 *
 * \note This procedure is only used for complete reduced row echelon forms
 *
 * \param pointer to block matrix input_bl
 *
 * \param pointer to multiline matrix input_ml
 *
 * \param rank of multiline matrix rank_input_ml
 *
 * \param sparse matrix output
 *
 * \param deleteIn, if 1 the input matrix in is deleted, if 0 in is not deleted.
 *
 * \param number of threads to copy data in parallel nthrds
 */
void copy_block_ml_matrices_to_sparse_matrix(sbm_fl_t **input_bl,
    sm_fl_ml_t **input_ml, ri_t rank_input_ml, sm_t **output,
    int deleteIn, int nthrds);

#if GBLA_WITH_FFLAS
/**
 * \brief Copies data from block matrix into DNS matrix format.
 *
 * \param pointer to block matrix source
 *
 * \param pointer to DNS matrix destination
 */
void copy_block_ml_matrix_to_dns_matrix(sbm_fl_t **source, DNS **destination);
#endif

/**
 * \brief Copies data from block matrix in to multiline matrix out. If deleteIn
 * is set the input block matrix in is deleted.
 *
 * \param pointer to block matrix in, input matrix
 *
 * \param multiline matrix out, output matrix
 *
 * \param deleteIn, if 1 the input matrix in is deleted, if 0 in is not deleted.
 *
 * \param number of threads to copy data in parallel nthrds
 *
 * \return multiline matrix generated out of block input matrix
 */
sm_fl_ml_t *copy_block_matrix_to_multiline_matrix(sbm_fl_t **input,
     sm_fl_ml_t *out, int deleteIn, int nthrds);

/**
 * \brief Copies multiline format matrix A_in to block format matrix B
 *
 * \param pointer to pointer to multiline matrix A_in
 *
 * \param block height bheight
 *
 * \param block width bwidth
 -*
 * \param freeing memory? if 1 then memory is freed, otherwise not, free_memory
 *
 * \param number of threads for parallel computations nthrds
 *
 * \return matrix A in block format
 */
sbm_fl_t *copy_multiline_to_block_matrix_rl(sm_fl_ml_t **A_in,
    ri_t bheight, ci_t bwidth, int free_memory, int nthrds);
#endif

double compute_density(nnz_t nnz, ri_t nrows, ri_t ncols) ;

sm_t *sort_schreyer_matrix(sm_t *M);


/* vim:sts=2:sw=2:ts=2:
 */
