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

inline mat_t *initialize_gbla_matrix(const spd_t *spd, const gb_t *basis)
{
  mat_t *mat  = (mat_t *)malloc(sizeof(mat_t));

  mat->A  = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  mat->B  = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  mat->C  = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  mat->D  = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  mat->DR = NULL;

  mat->mod  = basis->mod;
  mat->bs   = __GBLA_SIMD_BLOCK_SIZE;
  mat->ncl  = spd->col->nlm;
  mat->ncr  = spd->col->load - spd->col->nlm;
  mat->nru  = spd->selu->load;
  mat->nrl  = spd->sell->load;
  mat->rbu  = get_number_of_row_blocks(spd->selu, mat->bs);
  mat->rbl  = get_number_of_row_blocks(spd->sell, mat->bs);
  mat->cbl  = get_number_of_left_column_blocks(spd->col, mat->bs);
  mat->cbr  = get_number_of_right_column_blocks(spd->col, mat->bs);
#if MATRIX_DEBUG
  printf("cbl %u | cbr %u\n",mat->cbl, mat->cbr);
#endif

  // initialize parts of gbla matrix with known dimensions

  // D exists always
  init_dbm(mat->D, spd->sell->load, spd->col->load - spd->col->nlm);
  
  // if no upper part A & B are NULL and so is C
  if (spd->selu->load == 0 || spd->col->nlm == 0) {
    mat->A->blocks  = NULL;
    mat->A->nrows   = spd->selu->load;
    mat->A->ncols   = spd->col->nlm;;
    mat->B->blocks  = NULL;
    mat->A->nrows   = spd->selu->load;
    mat->B->ncols   = spd->col->load - spd->col->nlm;;
    mat->C->blocks  = NULL;
    mat->C->nrows   = spd->sell->load;
    mat->C->ncols   = spd->col->nlm;;
  } else {
    init_sb(mat->A, spd->selu->load, spd->col->nlm);
    init_sb(mat->C, spd->sell->load, spd->col->nlm);
    init_dbm(mat->B, spd->selu->load, spd->col->load - spd->col->nlm);
  }

#if MATRIX_DEBUG
    printf("A (%u x %u)\n", mat->A->nrows, mat->A->ncols);
    printf("B (%u x %u)\n", mat->B->nrows, mat->B->ncols);
    printf("C (%u x %u)\n", mat->C->nrows, mat->C->ncols);
    printf("D (%u x %u)\n", mat->D->nrows, mat->D->ncols);
#endif

  return mat;
}

void free_gbla_matrix(mat_t *mat)
{
  nelts_t i, j;

  // A, C and D are already freed, just check again
  
  // B is dense block matrix
  for (i=0; i<mat->rbu; ++i) {
    for (j=0; j<mat->cbr; ++j) {
      if (mat->B->blocks[i][j].val != NULL)
        free(mat->B->blocks[i][j].val);
    }
    free(mat->B->blocks[i]);
  }
  free(mat->B->blocks);
  free(mat->B);

  // DR is a dense row matrix
  for (i=0; i<mat->DR->nrows; ++i) {
    free(mat->DR->row[i]->piv_val);
    free(mat->DR->row[i]->init_val);
    free(mat->DR->row[i]->val);
  }
  free(mat->DR->row);
  free(mat->DR);
  mat->DR = NULL;

  free(mat);
  mat = NULL;
}

inline dbr_t *initialize_dense_block_row(const nelts_t nb, const bi_t bs)
{
  nelts_t i;

  dbr_t *dbr  = (dbr_t *)malloc(sizeof(dbr_t));
  dbr->ctr    = (nelts_t *)malloc(nb * sizeof(nelts_t));
  dbr->cf     = (coeff_t **)malloc(nb * sizeof(coeff_t *));
  for (i=0; i<nb; ++i)
    dbr->cf[i]  = (coeff_t *)malloc(bs * sizeof(coeff_t));

  return dbr;
}

inline void free_dense_block_row(dbr_t *dbr, const nelts_t nb)
{
  nelts_t i;

  free(dbr->ctr);
  for (i=0; i<nb; ++i)
    free(dbr->cf[i]);
  free(dbr->cf);
  free(dbr);
  dbr = NULL;
}

inline ci_t get_number_of_left_column_blocks(const pre_t *col, const nelts_t bs)
{
  return (ci_t) ceil((float) col->nlm/ bs);
}

inline ci_t get_number_of_right_column_blocks(const pre_t *col, const nelts_t bs)
{
  return (ci_t) ceil((float) (col->load - col->nlm)/ bs);
}

inline ri_t get_number_of_row_blocks(const sel_t *sel, const nelts_t bs)
{
  return (ri_t) ceil((float) sel->load / bs);
}

inline void reset_buffer(dbr_t *dbr, const nelts_t ncb, const bi_t bs)
{
  nelts_t i;

  memset(dbr->ctr, 0, ncb * sizeof(nelts_t));
  for (i=0; i<ncb; ++i)
    memset(dbr->cf[i], 0, bs * sizeof(coeff_t));
}

inline void allocate_dense_block(dbm_fl_t *A, const nelts_t rbi, const nelts_t bir,
    const bi_t bs)
{
  A->blocks[rbi][bir].val = (re_t *)calloc(bs * bs, sizeof(re_t));
}

inline void allocate_sparse_block(sb_fl_t *A, const nelts_t rbi, const nelts_t bir,
    const bi_t bs)
{
  bi_t i;

  A->blocks[rbi][bir].val = (re_t **)malloc(bs * sizeof(re_t *));
  A->blocks[rbi][bir].pos = (bi_t **)malloc(bs * sizeof(bi_t *));
  A->blocks[rbi][bir].sz  = (bi_t *)malloc(bs * sizeof(bi_t));
  for (i=0; i<bs; ++i) {
    A->blocks[rbi][bir].val[i]  = NULL;
    A->blocks[rbi][bir].pos[i]  = NULL;
    A->blocks[rbi][bir].sz[i]   = 0;
  }
}

inline void allocate_sparse_row_in_block(sb_fl_t *A, const nelts_t rbi,
    const nelts_t bir, const bi_t rib, const nelts_t sz, const bi_t bs)
{
  //printf("rbi %u | bir %u | rib %u | sz %u\n",rbi, bir, rib, sz);
  A->blocks[rbi][bir].val[rib]  = (re_t *)malloc(sz * sizeof(re_t));
  A->blocks[rbi][bir].pos[rib]  = (bi_t *)malloc(sz * sizeof(bi_t));
}

inline void write_to_sparse_row(sb_fl_t *A, const coeff_t *cf, const nelts_t rbi,
    const nelts_t bir, const nelts_t rib, const nelts_t sz,
    const bi_t bs, const coeff_t mod)
{
  bi_t i;
  for (i=bs; i>0; --i) {
  //for (i=0; i<bs; ++i) {
    if (cf[i-1] != 0) {
      A->blocks[rbi][bir].sz[rib]++;
      A->blocks[rbi][bir].val[rib][sz - A->blocks[rbi][bir].sz[rib]]  =
        (re_t)((re_m_t)mod - cf[i-1]);
      A->blocks[rbi][bir].pos[rib][sz - A->blocks[rbi][bir].sz[rib]]  = i-1;
      //printf("i %u | A->blocks[%u][%u].pos[%u][%u - %u] = %u\n",i, rbi,bir,rib, sz,A->blocks[rbi][bir].sz[rib],A->blocks[rbi][bir].pos[rib][sz - A->blocks[rbi][bir].sz[rib]]);
    }
  }
}

inline void write_to_dense_row(dbm_fl_t *A, const coeff_t *cf, const nelts_t rbi,
    const nelts_t bir, const nelts_t rib, const bi_t bs)
{
  bi_t i;

  for (i=0; i<bs; ++i) {
    A->blocks[rbi][bir].val[(rib*bs)+i] = cf[i];
  }
}

inline void store_in_matrix(sb_fl_t *A, dbm_fl_t *B, const dbr_t *dbr,
    const nelts_t rbi, const nelts_t rib, const nelts_t ncb, const nelts_t fr,
    const bi_t bs, const coeff_t mod)
{
  nelts_t i;

  // calculate index of last block on left side
  // if there is nothing on the lefthand side what can happen when interreducing
  // the initial input elements then we have to adjust fbr to 0
  const nelts_t fbr = fr == 0 ? 0 : (fr-1)/bs + 1;

  // do sparse left side A
  for (i=0; i<fbr; ++i) {
    //printf("dbr->ctr[%u] = %u\n",i,dbr->ctr[i]);
    if (dbr->ctr[i] > 0) {
      if (A->blocks[rbi][i].val == NULL)
        allocate_sparse_block(A, rbi, i, bs);
      allocate_sparse_row_in_block(A, rbi, i, rib, dbr->ctr[i], bs);
      write_to_sparse_row(A, dbr->cf[i], rbi, i, rib, dbr->ctr[i], bs, mod);
    }
  }

  // do dense right side B
  for (i=0; i<ncb-fbr; ++i) {
    if (dbr->ctr[fbr+i] > 0) {
      if (B->blocks[rbi][i].val == NULL)
        allocate_dense_block(B, rbi, i, bs);
      write_to_dense_row(B, dbr->cf[fbr+i], rbi, i, rib, bs);
    }
  }
}


inline mat_t *generate_gbla_matrix(const gb_t *basis, const spd_t *spd, const int nthreads)
{
  mat_t *mat  = initialize_gbla_matrix(spd, basis);
  #pragma omp parallel num_threads(nthreads)
  {
    #pragma omp single nowait
    {
    // fill the upper part AB
    for (ri_t i=0; i<mat->rbu; ++i) {
      #pragma omp task
      {
        generate_row_blocks(mat->A, mat->B, i, spd->selu->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, spd->selu, spd->col);
      }
    }
    // fill the lower part CD
    for (ri_t i=0; i<mat->rbl; ++i) {
      #pragma omp task
      {
        generate_row_blocks(mat->C, mat->D, i, spd->sell->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, spd->sell, spd->col);
      }
    }
    }
    #pragma omp taskwait
  }

  return mat;
}

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
