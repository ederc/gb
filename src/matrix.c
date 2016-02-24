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
  if (verbose > 1)
    gettimeofday(&t_complete, NULL);
  // A^-1 * B
  if (verbose > 1) {
    printf("---------------------------------------------------------------------\n");
    printf("GBLA Matrix Reduction\n");
    printf("---------------------------------------------------------------------\n");
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing A ...");
    fflush(stdout);
  }
  /*
  printf("---- A1 -----\n");
  if (mat->A != NULL && mat->A->blocks != NULL) {
  for (int ii=0; ii<mat->A->nrows; ++ii) {
    printf("%u || ", ii);
    for (int jj=0; jj<mat->A->blocks[0][0].sz[ii]; ++jj) {
      printf("%u - %u | ",mat->A->blocks[0][0].val[ii][jj],mat->A->blocks[0][0].pos[ii][jj]);
    }
    printf("\n");
  }
  }
  if (mat->B != NULL && mat->B->blocks != NULL) {
  printf("---- B1 -----\n");
  for (int ii=0; ii<mat->bs; ++ii) {
    printf("%u || ", ii);
    for (int jj=0; jj<mat->bs; ++jj) {
      printf("%u ",mat->B->blocks[0][0].val[mat->bs*ii+jj]);
    }
    printf("\n");
  }
  }
  */
  if (mat->A->blocks != NULL) {
    if (elim_fl_A_sparse_dense_block(&(mat->A), mat->B, mat->mod, nthreads)) {
      printf("Error while reducing A.\n");
      return -1;
    }
  }
  if (verbose > 1) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 2) {
    print_mem_usage();
  }
  /*
  if (mat->B != NULL && mat->B->blocks != NULL) {
  printf("---- B2 -----\n");
  for (int ii=0; ii<mat->bs; ++ii) {
    printf("%u || ", ii);
    for (int jj=0; jj<mat->bs; ++jj) {
      printf("%u ",mat->B->blocks[0][0].val[mat->bs*ii+jj]);
    }
    printf("\n");
  }
  }
  printf("\n");
  */
  // reducing submatrix C to zero using methods of Faugère & Lachartre
  if (verbose > 1) {
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
  if (verbose > 1) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 2) {
    print_mem_usage();
  }
  // copy block D to dense wide (re_l_t) representation
  dm_t *D_red = copy_block_to_dense_matrix(&(mat->D), nthreads);
  D_red->mod  = mat->mod;

  // eliminate D_red using a structured Gaussian Elimination process on the rows
  ri_t rank_D = 0;
  // echelonizing D to zero using methods of Faugère & Lachartre
  if (verbose > 1) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing D ...");
    fflush(stdout);
  }
  if (D_red->nrows > 0)
    //rank_D = elim_fl_dense_D(D_red, nthreads);
    rank_D = elim_fl_dense_D_completely(D_red, nthreads);
  if (verbose > 1) {
    printf("%9.3f sec (rank D: %u)\n",
        walltime(t_load_start) / (1000000), rank_D);
  }
  if (verbose > 2) {
    print_mem_usage();
  }
  if (verbose > 1) {
    printf("---------------------------------------------------------------------\n");
    printf("%-38s","Reduction completed ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        walltime(t_complete) / (1000000));
    if (verbose > 1) 
      print_mem_usage();
  }

  mat->DR = D_red;

  return rank_D;
}
