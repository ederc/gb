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

ri_t reduce_gbla_matrix(mat_t * mat, int verbose, int nthreads)
{
  /*  timing structs */
  struct timeval t_load_start;
  struct timeval t_complete;
  if (verbose > 2)
    gettimeofday(&t_complete, NULL);
  /* A^-1 * B */
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
      return 1;
    }
  }
  if (verbose > 2) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  /* reducing submatrix C to zero using methods of Faugère & Lachartre */
  if (verbose > 2) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing C ...");
    fflush(stdout);
  }
  if (mat->C->blocks != NULL) {
    if (elim_fl_C_sparse_dense_block(mat->B, &(mat->C), mat->D, mat->mod, nthreads)) {
      printf("Error while reducing C.\n");
      return 1;
    }
  }
  if (verbose > 2) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  /* copy block D to dense wide (re_l_t) representation */
  mat->DR = copy_block_to_dense_matrix(&(mat->D), nthreads, 1);
  mat->DR->mod  = mat->mod;
#if 0
  printf("number of rows of DR %u\n", mat->DR->nrows);
  for (int ii=0; ii<mat->DR->nrows; ++ii) {
    printf("ROW %d\n",ii);
    if (mat->DR->row[ii]->init_val == NULL)
      printf("NULL!");
    else {
      printf("%u || ", mat->DR->row[ii]->lead);
      for (int jj=0; jj<mat->DR->ncols; ++jj)
#if defined(GBLA_USE_UINT16) || defined(GBLA_USE_UINT32)
        printf("%u (%u)  ", mat->DR->row[ii]->init_val[jj], jj+mat->ncl);
#else
        printf("%.0f  ", mat->DR->row[ii]->init_val[jj]);
#endif
    }
    printf("\n");
  }
#endif

  /* eliminate mat->DR using a structured Gaussian Elimination process on the rows */
  nelts_t rank_D = 0;
  /* echelonizing D to zero using methods of Faugère & Lachartre */
  if (verbose > 2) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing D ...");
    fflush(stdout);
  }
  if (mat->DR->nrows > 0) {
    if (nthreads == 1) {
      rank_D = elim_fl_dense_D_completely(mat->DR, nthreads);
    } else {
      rank_D = elim_fl_dense_D(mat->DR, nthreads);
      nelts_t l;
      for (l=1; l<mat->DR->rank; ++l) {
      /* for (l=(int)(mat->DR->rank-1); l>0; --l) { */
        copy_piv_to_val(mat->DR, mat->DR->rank-l-1);
        completely_reduce_D(mat->DR, mat->DR->rank-l-1);
      }
    }
  }
  if (verbose > 2) {
    printf("%9.3f sec %5d %5d %5d\n",
        walltime(t_load_start) / (1000000), rank_D, mat->DR->nrows - rank_D, mat->DR->nrows);
  }
  /* if we simplify, then copy B to dense row representation */
  if (mat->sl > 0 && mat->B->blocks != NULL) {
    /* first copy B to BR (dense row format) */
    mat->BR = copy_block_to_dense_matrix(&(mat->B), nthreads, 0);
    mat->BR->mod  = mat->mod;
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

ri_t reduce_gbla_matrix_keep_A(mat_t *mat, int verbose, int nthreads)
{
  /*  timing structs */
  struct timeval t_load_start;
  struct timeval t_complete;

  if (verbose > 2)
    gettimeofday(&t_complete, NULL);
  /* A^-1 * B */
  if (verbose > 2) {
    printf("---------------------------------------------------------------------------\n");
    printf("GBLA Matrix Reduction\n");
    printf("---------------------------------------------------------------------------\n");
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Storing A in C ...");
    fflush(stdout);
  }

  if (mat->AR->row != NULL) {
    if (elim_fl_C_sparse_dense_keep_A(mat->CR, &(mat->AR), mat->mod, nthreads)) {
      printf("Error while reducing A.\n");
      return 1;
    }
  }
  if (verbose > 2) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  /* reducing submatrix C to zero using methods of Faugère & Lachartre */
  if (verbose > 2) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Copying C to sparse block representation ...");
    fflush(stdout);
  }
  if (mat->CR->row != NULL) {
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
  if (mat->C != NULL) {
    if (elim_fl_C_sparse_dense_block(mat->B, &(mat->C), mat->D, mat->mod, nthreads)) {
      printf("Error while reducing A.\n");
      return 1;
    }
  }
  if (verbose > 2) {
    printf("%9.3f sec\n",
        walltime(t_load_start) / (1000000));
  }
  if (verbose > 3) {
    print_mem_usage();
  }
  /* copy block D to dense wide (re_l_t) representation */
  mat->DR = copy_block_to_dense_matrix(&(mat->D), nthreads, 1);
  mat->DR->mod  = mat->mod;

  /* eliminate mat->DR using a structured Gaussian Elimination process on the rows */
  nelts_t rank_D = 0;
  /* echelonizing D to zero using methods of Faugère & Lachartre */
  if (verbose > 2) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Reducing D ...");
    fflush(stdout);
  }
  if (mat->DR->nrows > 0)
    /* rank_D = elim_fl_dense_D(mat->DR, nthreads); */
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

static inline void write_sparse_compact_row(src_t **rows,
    const mat_gb_block_t *om, const nelts_t idx,
    const mat_gb_meta_data_t *meta)
{
  const nelts_t max = meta->bs < meta->nr_CD - idx*meta->bs ?
    meta->bs : meta->nr_CD - idx*meta->bs;

  nelts_t i, j, k;
  nelts_t ctr;

  for (i=0; i<max; ++i) {
    src_t *row  = (src_t *)malloc((2*meta->nc_BD+1) * sizeof(src_t));
    /* go only over D part, C is already zero */
    ctr = 1;
    for (j=meta->ncb_AC; j<meta->ncb; ++j) {
      printf("j %u\n", j);
      if (om[j].len != NULL) {
        printf("drin?");
        for (k=om[j].len[i]; k<om[j].len[i+1]; ++k) {
          row[ctr]  = (src_t)om[j].pos[k] + (j-meta->ncb_AC)*meta->bs + meta->nc_AC;
          row[ctr+1]  = om[j].val[k];
          printf("%u | %u || ", row[ctr+1], row[ctr]);
          ctr = ctr+2;
        }
      }
      printf("\n");
    }
    if (ctr>1) {
      row[0]  = ctr;
      row = realloc(row, ctr * sizeof(src_t));
      rows[i + idx*meta->bs] = row;
      if (row[2] != 1) {
        normalize_row_c(row, meta->mod);
      }  
    } else {
      free(row);
      row = NULL;
      rows[i + idx*meta->bs] = row;
    }
  }
}

smc_t *convert_mat_gb_to_smc_format(const mat_gb_block_t *om,
    const mat_gb_meta_data_t *meta, const int t)
{
  nelts_t i, ctr;

  smc_t *nm = initialize_sparse_compact_matrix(meta->nr_CD, meta->nc_AC,
      meta->nc_BD, meta->mod);

  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (i=0; i<nm->nr; i=i+meta->bs) {
      #pragma omp task
      {
        /* printf("constructs row %u\n",i); */
        write_sparse_compact_row(nm->row, om, i, meta);
      }
    }
    }
  }
  ctr = 0;
  for (i=0; i<nm->nr; ++i)
    if (nm->row[i] != NULL)
      ctr++;
  nm->rk  = ctr;

  return nm;
}
