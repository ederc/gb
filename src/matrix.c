/* gb: Gröbner Basis
 * Copyright (C) 2015 Christian Eder <ederc@mathematik.uni-kl.de>
 * This file is part of gb.
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
 * \file matrix.c
 * \brief Implementation of the construction and conversion from and to groebner
 * basis matrices.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "matrix.h"

int reduce_gbla_matrix(mat_t * mat, int verbose, int nthreads)
{
  /*  timing structs */
  struct timeval t_load_start;
  struct timeval t_complete;
  if (verbose > 2)
    gettimeofday(&t_complete, NULL);
  // A^-1 * B
  if (verbose > 2) {
    printf("---------------------------------------------------------------------------\n");
    printf("GBLA Matrix Reduction\n");
    printf("---------------------------------------------------------------------------\n");
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing A ...");
    fflush(stdout);
  }
  if (mat->A->blocks != NULL) {
    if (elim_fl_A_sparse_dense_block(&(mat->A), mat->B, mat->mod, nthreads)) {
      printf("Error while reducing A.\n");
      return -1;
    }
  }
  if (verbose > 2) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  // reducing submatrix C to zero using methods of Faugère & Lachartre
  if (verbose > 2) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing C ...");
    fflush(stdout);
  }
  if (mat->C->blocks != NULL) {
    if (elim_fl_C_sparse_dense_block(mat->B, &(mat->C), mat->D, 1, mat->mod, nthreads)) {
      printf("Error while reducing C.\n");
      return -1;
    }
  }
  if (verbose > 2) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  // copy block D to dense wide (re_l_t) representation
  mat->DR = copy_block_to_dense_matrix(&(mat->D), nthreads, 1);
  mat->DR->mod  = mat->mod;

  // eliminate mat->DR using a structured Gaussian Elimination process on the rows
  ri_t rank_D = 0;
  // echelonizing D to zero using methods of Faugère & Lachartre
  if (verbose > 2) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing D ...");
    fflush(stdout);
  }
  if (mat->DR->nrows > 0)
    //rank_D = elim_fl_dense_D(mat->DR, nthreads);
    rank_D = elim_fl_dense_D_completely(mat->DR, nthreads);
  if (verbose > 2) {
    printf("%9.3f sec %5d %5d %5d\n",
        walltime(t_load_start) / (1000000), rank_D, mat->DR->nrows - rank_D, mat->DR->nrows);
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  if (verbose > 2) {
    printf("---------------------------------------------------------------------------\n");
    printf("%-38s","Reduction completed ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        walltime(t_complete) / (1000000));
    if (verbose > 3)
      print_mem_usage();
  }

  return rank_D;
}

int reduce_gbla_matrix_keep_A(mat_t *mat, int verbose, int nthreads)
{
  /*  timing structs */
  struct timeval t_load_start;
  struct timeval t_complete;
  if (verbose > 2)
    gettimeofday(&t_complete, NULL);
  // A^-1 * B
  if (verbose > 2) {
    printf("---------------------------------------------------------------------------\n");
    printf("GBLA Matrix Reduction\n");
    printf("---------------------------------------------------------------------------\n");
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Storing A in C ...");
    fflush(stdout);
  }
  if (mat->AR != NULL) {
    if (elim_fl_C_sparse_dense_keep_A(mat->CR, &(mat->AR), mat->mod, nthreads)) {
      printf("Error while reducing A.\n");
      return -1;
    }
  }
  if (verbose > 2) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  // reducing submatrix C to zero using methods of Faugère & Lachartre
  if (verbose > 2) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Copying C to sparse block representation ...");
    fflush(stdout);
  }
  if (mat->CR != NULL) {
    mat->C  = copy_sparse_to_block_matrix(mat->CR, nthreads);
    free_sparse_matrix(&(mat->CR), nthreads);
  }
  if (verbose > 2) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  if (verbose > 2) {
    printf("%-38s","Reducing C to zero ...");
    fflush(stdout);
  }
  if (mat->AR != NULL) {
    if (elim_fl_C_sparse_dense_block(mat->B, &(mat->C), mat->D, 0, mat->mod, nthreads)) {
      printf("Error while reducing A.\n");
      return -1;
    }
  }
  if (verbose > 2) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  // copy block D to dense wide (re_l_t) representation
  mat->DR = copy_block_to_dense_matrix(&(mat->D), nthreads, 1);
  mat->DR->mod  = mat->mod;

  // eliminate mat->DR using a structured Gaussian Elimination process on the rows
  ri_t rank_D = 0;
  // echelonizing D to zero using methods of Faugère & Lachartre
  if (verbose > 2) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing D ...");
    fflush(stdout);
  }
  if (mat->DR->nrows > 0)
    //rank_D = elim_fl_dense_D(mat->DR, nthreads);
    rank_D = elim_fl_dense_D_completely(mat->DR, nthreads);
  if (verbose > 2) {
    printf("%9.3f sec %5d %5d %5d\n",
        walltime(t_load_start) / (1000000), rank_D, mat->DR->nrows - rank_D, mat->DR->nrows);
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  if (verbose > 2) {
    printf("---------------------------------------------------------------------------\n");
    printf("%-38s","Reduction completed ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        walltime(t_complete) / (1000000));
    if (verbose > 3)
      print_mem_usage();
  }

  return rank_D;
}
