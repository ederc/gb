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
 * \file matrix.h
 * \brief Implementation of the construction and conversion from and to groebner
 * basis matrices.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_MATRIX_H
#define GB_MATRIX_H

#include "gb_config.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <gbla/elimination.h>
#include <cli/io.h>
#include "types.h"
#include "hash.h"

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#ifndef MATRIX_DEBUG
#define MATRIX_DEBUG 0
#endif

/**
 * \brief Computes the number of column blocks on the left side in the gbla
 * matrix for a given set of monomials and a given block size.
 *
 * \param symbolic preprocessing monomials col
 *
 * \param block size bs
 *
 * \return number of left hand side column blocks
 */
static inline ci_t get_number_of_left_column_blocks(const pre_t *col, const nelts_t bs)
{
  return (ci_t) ceil((float) col->nlm/ bs);
}

/**
 * \brief Computes the number of column blocks on the right side in the gbla
 * matrix for a given set of monomials and a given block size.
 *
 * \param symbolic preprocessing monomials col
 *
 * \param block size bs
 *
 * \return number of right hand side column blocks
 */
static inline ci_t get_number_of_right_column_blocks(const pre_t *col, const nelts_t bs)
{
  return (ci_t) ceil((float) (col->load - col->nlm)/ bs);
}

/**
 * \brief Computes the number of row blocks in the gbla matrix for a given
 * selection (could be upper or lower part) and a given block size.
 *
 * \param symbolic preprocessing selection sel
 *
 * \param block size bs
 *
 * \return number of row blocks
 */
static inline ri_t get_number_of_row_blocks(const sel_t *sel, const nelts_t bs)
{
  return (ri_t) ceil((float) sel->load / bs);
}

/**
 * \brief Initilializes the gbla matrix corresponding to the selection done
 * during symbolic preprocessing. It also converts the polynomial representation
 * to rows resp. blocks in the matrix. This version has non-block submatrices A
 * and C resp. mat->AR and mat->CR
 *
 * \param symbolic preprocessing data spd
 *
 * \param intermediate groebner basis basis
 *
 * \return gbla matrix
 */
static inline mat_t *initialize_gbla_matrix_keep_A(const spd_t *spd, const gb_t *basis)
{
  mat_t *mat  = (mat_t *)malloc(sizeof(mat_t));

  mat->AR = (sm_fl_t *)malloc(sizeof(sm_fl_t));
  mat->B  = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  mat->CR = (sm_fl_t *)malloc(sizeof(sm_fl_t));
  mat->D  = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  mat->A  = NULL;
  mat->C  = NULL;
  mat->BR = NULL;
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
    mat->AR->row    = NULL;
    mat->AR->nrows  = spd->selu->load;
    mat->AR->ncols  = spd->col->nlm;;
    mat->B->blocks  = NULL;
    mat->B->nrows   = spd->selu->load;
    mat->B->ncols   = spd->col->load - spd->col->nlm;;
    mat->CR->row    = NULL;
    mat->CR->nrows  = spd->sell->load;
    mat->CR->ncols  = spd->col->nlm;;
  } else {
    init_sm_no_row_allocation(mat->AR, spd->selu->load, spd->col->nlm);
    init_sm_no_row_allocation(mat->CR, spd->sell->load, spd->col->nlm);
    init_dbm(mat->B, spd->selu->load, spd->col->load - spd->col->nlm);
  }

#if MATRIX_DEBUG
    printf("A (%u x %u)\n", mat->AR->nrows, mat->AR->ncols);
    printf("B (%u x %u)\n", mat->B->nrows, mat->B->ncols);
    printf("C (%u x %u)\n", mat->CR->nrows, mat->CR->ncols);
    printf("D (%u x %u)\n", mat->D->nrows, mat->D->ncols);
#endif

  return mat;
}

/**
 * \brief Initilializes the gbla matrix corresponding to the selection done
 * during symbolic preprocessing. It also converts the polynomial representation
 * to rows resp. blocks in the matrix.
 *
 * \param symbolic preprocessing data spd
 *
 * \param intermediate groebner basis basis
 *
 * \return gbla matrix
 */
static inline mat_t *initialize_gbla_matrix(const spd_t *spd, const gb_t *basis)
{
  mat_t *mat  = (mat_t *)malloc(sizeof(mat_t));

  mat->A  = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  mat->B  = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  mat->C  = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  mat->D  = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  mat->AR = NULL;
  mat->BR = NULL;
  mat->CR = NULL;
  mat->DR = NULL;

  mat->sl   = basis->sl;
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

/**
 * \brief Frees gbla matrix after reduction.
 *
 * \param gbla matrix mat
 */
static inline void free_gbla_matrix(mat_t **mat_in)
{
  mat_t *mat  = *mat_in;

  // A, C and D are already freed, just check again
  if (mat->AR != NULL) {
    if (mat->AR->row != NULL) {
      free_sparse_matrix(&(mat->AR), 1);
    } else {
      free(mat->AR);
      mat->AR = NULL;
    }
  }
  if (mat->A != NULL) {
    if (mat->A->blocks != NULL) {
      free_sparse_submatrix(&(mat->A), 1);
    } else {
      free(mat->A);
      mat->A  = NULL;
    }
  }
  // mat->CR already freed after copying it to sparse block matrix mat->C
  if (mat->C != NULL) {
    if (mat->C->blocks != NULL) {
      free_sparse_submatrix(&(mat->C), 1);
    } else {
      free(mat->C);
      mat->C  = NULL;
    }
  }

  if (mat->D != NULL) {
    if (mat->D->blocks != NULL) {
      free_dense_submatrix(&(mat->D), 1);
    } else {
      free(mat->D);
      mat->D  = NULL;
    }
  }
  
  // B is dense block matrix
  if (mat->B != NULL) {
    if (mat->B->blocks != NULL) {
      free_dense_submatrix(&(mat->B), 1);
    } else {
      free(mat->B);
      mat->B  = NULL;
    }
  }
  if (mat->BR != NULL) {
    // DR is a dense row matrix
    free_dense_row_submatrix(&(mat->BR), 1);
  }

  if (mat->DR != NULL) {
    // DR is a dense row matrix
    free_dense_row_submatrix(&(mat->DR), 1);
  }

  free(mat);
  mat = NULL;

  *mat_in = mat;
}

/**
 * \brief Allocates memory for a sparse block in gbla matrix.
 *
 * \note Does not allocate memory for buffer in sb_fl_t data structure from
 * gbla.
 *
 * \param sparse block matrix A
 *
 * \param row block index rbi
 *
 * \param block index in row bir
 *
 * \param block size bs
 */
static inline void allocate_sparse_block(sb_fl_t *A, const nelts_t rbi, const nelts_t bir,
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

/**
 * \brief Allocates memory for a dense block in gbla matrix.
 *
 * \param dense block matrix A
 *
 * \param row block index rbi
 *
 * \param block index in row bir
 *
 * \param block size bs
 */
static inline void allocate_dense_block(dbm_fl_t *A, const nelts_t rbi, const nelts_t bir,
    const bi_t bs)
{
  A->blocks[rbi][bir].val = (re_t *)calloc(bs * bs, sizeof(re_t));
}

/**
 * \brief Initializes a dense block row for buffering values when converting a
 * polynomial to a row in the gbla matrix.
 *
 * \param number of blocks nb
 *
 * \param number of the first block on the righthandside of the gbla matrix fbr
 *
 * \param block size bs
 *
 * \return dense block row
 */
static inline dbr_t *initialize_dense_block_row(const nelts_t nb, const bi_t bs)
{
  nelts_t i;

  dbr_t *dbr  = (dbr_t *)malloc(sizeof(dbr_t));
  dbr->ctr    = (nelts_t *)calloc(nb, sizeof(nelts_t));
  dbr->cf     = (coeff_t **)malloc(nb * sizeof(coeff_t *));
  // allocate buffer for sparse rows
  for (i=0; i<nb; ++i)
    dbr->cf[i]  = (coeff_t *)calloc(bs, sizeof(coeff_t));
  return dbr;
}

/**
 * \brief Initializes a dense block row for buffering values when converting a
 * polynomial to a row in the gbla matrix. The sparse part (lefthand side) is
 * stored row-wise, the righthand side (dense part) is stored already as dense
 * block and thus these blocks are just linked into the matrix at the very end.
 *
 * \param number of blocks nb
 *
 * \param number of the first block on the righthandside of the gbla matrix fbr
 *
 * \param block size bs
 *
 * \return dense block row
 */
static inline dbr_t *initialize_dense_block_row_new(const nelts_t nb, const nelts_t fbr, const bi_t bs)
{
  nelts_t i;

  dbr_t *dbr  = (dbr_t *)malloc(sizeof(dbr_t));
  dbr->ctr    = (nelts_t *)calloc(nb, sizeof(nelts_t));
  if (fbr>0) {
    dbr->cf     = (coeff_t **)malloc(fbr * sizeof(coeff_t *));
    // allocate buffer for sparse rows
    for (i=0; i<fbr; ++i)
      dbr->cf[i]  = (coeff_t *)calloc(bs, sizeof(coeff_t));
  } else {
    dbr->cf = NULL;
  }
  dbr->bl     = (coeff_t **)malloc((nb-fbr) * sizeof(coeff_t *));
  // allocate buffer for dense blocks, those will be the final blocks in the
  // matrix
  for (i=0; i<nb-fbr; ++i)
    dbr->bl[i]  = (coeff_t *)calloc(bs*bs, sizeof(coeff_t));

  return dbr;
}

/**
 * \brief Frees dense block row.
 *
 * \param dense block row dbr
 *
 * \param number of blocks nb
 */
static inline void free_dense_block_row(dbr_t *dbr, const nelts_t nb)
{
  nelts_t i;

  free(dbr->ctr);
  for (i=0; i<nb; ++i)
    free(dbr->cf[i]);
  free(dbr->cf);
  free(dbr);
  dbr = NULL;
}


/**
 * \brief Frees dense block row.
 *
 * \param dense block row dbr
 *
 * \param index of first block on the righthand side of the gbla matrix fbr
 */
static inline void free_dense_block_row_new(dbr_t *dbr, const nelts_t fbr)
{
  nelts_t i;

  free(dbr->ctr);
  if (fbr>0) {
    for (i=0; i<fbr; ++i)
      free(dbr->cf[i]);
    free(dbr->cf);
  }
  free(dbr->bl);
  free(dbr);
  dbr = NULL;
}

/**
 * \brief Writes buffered data to sparse matrix row in block
 *
 * \param sparse block matrix A
 *
 * \param coefficients to be written cf
 *
 * \param row block index rbi
 *
 * \param block index in row bir
 *
 * \param row in block rib
 *
 * \param size of row sz
 *
 * \param block size bs
 *
 * \param field characteristic mod
 */
static inline void write_to_sparse_row(sm_fl_t *A, const coeff_t *cf, const nelts_t ri,
    const nelts_t bir, const bi_t bs)
{
  bi_t i;
  //for (i=bs; i>0; --i) {
  for (i=0; i<bs; ++i) {
    //if (cf[i-1] != 0) {
    if (cf[i] != 0) {
      A->row[ri][A->sz[ri]] = cf[i];
      A->pos[ri][A->sz[ri]] = i+ bir*bs;
      A->sz[ri]++;
    }
  }
}

/**
 * \brief Writes buffered data to sparse matrix row in block
 *
 * \param sparse block matrix A
 *
 * \param coefficients to be written cf
 *
 * \param row block index rbi
 *
 * \param block index in row bir
 *
 * \param row in block rib
 *
 * \param size of row sz
 *
 * \param block size bs
 *
 * \param field characteristic mod
 */
static inline void write_to_sparse_row_in_block(sb_fl_t *A, const coeff_t *cf, const nelts_t rbi,
    const nelts_t bir, const nelts_t rib, const nelts_t sz,
    const bi_t bs, const coeff_t mod)
{
#if 1
  bi_t i;
  int ctr  = -1;
  for (i=0; i<bs; i=i+8) {
    if (cf[i] != 0) {
      A->blocks[rbi][bir].val[rib][++ctr]  =
        (re_t)((re_m_t)mod - cf[i]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i;
    }
    if (cf[i+1] != 0) {
      A->blocks[rbi][bir].val[rib][++ctr]  =
        (re_t)((re_m_t)mod - cf[i+1]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i+1;
    }
    if (cf[i+2] != 0) {
      A->blocks[rbi][bir].val[rib][++ctr]  =
        (re_t)((re_m_t)mod - cf[i+2]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i+2;
    }
    if (cf[i+3] != 0) {
      A->blocks[rbi][bir].val[rib][++ctr]  =
        (re_t)((re_m_t)mod - cf[i+3]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i+3;
    }
    if (cf[i+4] != 0) {
      A->blocks[rbi][bir].val[rib][++ctr]  =
        (re_t)((re_m_t)mod - cf[i+4]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i+4;
    }
    if (cf[i+5] != 0) {
      A->blocks[rbi][bir].val[rib][++ctr]  =
        (re_t)((re_m_t)mod - cf[i+5]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i+5;
    }
    if (cf[i+6] != 0) {
      A->blocks[rbi][bir].val[rib][++ctr]  =
        (re_t)((re_m_t)mod - cf[i+6]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i+6;
    }
    if (cf[i+7] != 0) {
      A->blocks[rbi][bir].val[rib][++ctr]  =
        (re_t)((re_m_t)mod - cf[i+7]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i+7;
    }
  }
  A->blocks[rbi][bir].sz[rib] = sz;
#else
  bi_t i;
  bi_t ctr  = sz;
  for (i=bs; i>0; i=i-8) {
    if (cf[i-1] != 0) {
      A->blocks[rbi][bir].val[rib][--ctr]  =
        (re_t)((re_m_t)mod - cf[i-1]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i-1;
    }
    if (cf[i-2] != 0) {
      A->blocks[rbi][bir].val[rib][--ctr]  =
        (re_t)((re_m_t)mod - cf[i-2]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i-2;
    }
    if (cf[i-3] != 0) {
      A->blocks[rbi][bir].val[rib][--ctr]  =
        (re_t)((re_m_t)mod - cf[i-3]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i-3;
    }
    if (cf[i-4] != 0) {
      A->blocks[rbi][bir].val[rib][--ctr]  =
        (re_t)((re_m_t)mod - cf[i-4]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i-4;
    }
    if (cf[i-5] != 0) {
      A->blocks[rbi][bir].val[rib][--ctr]  =
        (re_t)((re_m_t)mod - cf[i-5]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i-5;
    }
    if (cf[i-6] != 0) {
      A->blocks[rbi][bir].val[rib][--ctr]  =
        (re_t)((re_m_t)mod - cf[i-6]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i-6;
    }
    if (cf[i-7] != 0) {
      A->blocks[rbi][bir].val[rib][--ctr]  =
        (re_t)((re_m_t)mod - cf[i-7]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i-7;
    }
    if (cf[i-8] != 0) {
      A->blocks[rbi][bir].val[rib][--ctr]  =
        (re_t)((re_m_t)mod - cf[i-8]);
      A->blocks[rbi][bir].pos[rib][ctr]  = i-8;
    }
  }
  A->blocks[rbi][bir].sz[rib] = sz;
#endif
}

/**
 * \brief Writes buffered data to dense matrix row in block
 *
 * \param sparse block matrix A
 *
 * \param coefficients to be written cf
 *
 * \param row block index rbi
 *
 * \param block index in row bir
 *
 * \param row in block rib
 *
 * \param block size bs
 */
static inline void write_to_dense_row(dbm_fl_t *A, const coeff_t *cf, const nelts_t rbi,
    const nelts_t bir, const nelts_t rib, const bi_t bs)
{
#if 1
    memcpy(A->blocks[rbi][bir].val+(rib*bs), cf, bs*sizeof(coeff_t));
#else
  bi_t i;

  for (i=0; i<bs; i=i+4) {
    A->blocks[rbi][bir].val[(rib*bs)+i]   = cf[i];
    A->blocks[rbi][bir].val[(rib*bs)+i+1] = cf[i+1];
    A->blocks[rbi][bir].val[(rib*bs)+i+2] = cf[i+2];
    A->blocks[rbi][bir].val[(rib*bs)+i+3] = cf[i+3];
  }
#endif
}


/**
 * \brief Allocates memory for sparse row in gbla matrix block.
 *
 * \param sparse block matrix A
 *
 * \param row block index rbi
 *
 * \param block index in row bir
 *
 * \param row in block rib
 *
 * \param size of row sz
 *
 * \param block size bs
 */ 
static inline void allocate_sparse_row_in_block(sb_fl_t *A, const nelts_t rbi,
    const nelts_t bir, const bi_t rib, const nelts_t sz, const bi_t bs)
{
  A->blocks[rbi][bir].val[rib]  = (re_t *)malloc(sz * sizeof(re_t));
  A->blocks[rbi][bir].pos[rib]  = (bi_t *)malloc(sz * sizeof(bi_t));
#if MARTIX_DEBUG
  printf("rbi %u | bir %u | rib %u | sz %u || %p | %p\n",rbi, bir, rib, sz, A->blocks[rbi][bir].val[rib],A->blocks[rbi][bir].pos[rib]);
#endif
}


/**
 * \brief Allocates memory for sparse row in gbla matrix.
 *
 * \param sparse block matrix A
 *
 * \param row index ri
 *
 * \param size of row sz
 */ 
static inline void allocate_sparse_row(sm_fl_t *A, const nelts_t ri, const nelts_t sz)
{
  A->row[ri]  = (re_t *)malloc(sz * sizeof(re_t));
  A->pos[ri]  = (ci_t *)malloc(sz * sizeof(ci_t));
}

/**
 * \brief After we have stored the data from one polynomial in the dense buffer,
 * we now write from the buffer to the matrix.
 *
 * \param sparse row matrix A
 *
 * \param dense block matrix B
 *
 * \param dense block row buffer dbr
 *
 * \param row block index rbi
 *
 * \param row index n block rib
 *
 * \param number of column blocks ncb
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param field characteristic mod
 */
static inline void store_in_matrix_keep_A(sm_fl_t *A, dbm_fl_t *B, const dbr_t *dbr,
    const nelts_t rbi, const nelts_t rib, const nelts_t ncb, const nelts_t fr,
    const bi_t bs, const coeff_t mod)
{
  nelts_t i;

  // calculate index of last block on left side
  // if there is nothing on the lefthand side what can happen when interreducing
  // the initial input elements then we have to adjust fbr to 0
  const nelts_t fbr = fr == 0 ? 0 : (fr-1)/bs + 1;

  nelts_t nl  = 0;
  // get number of elements in sparse row in A
  for (i=0; i<fbr; ++i)
    nl  +=  dbr->ctr[i];

  if (nl > 0) {
    // allocate memory for rows in A
    allocate_sparse_row(A, rbi*bs + rib, nl);
    // do sparse left side A
    for (i=0; i<fbr; ++i) {
      //printf("dbr->ctr[%u] = %u\n",i,dbr->ctr[i]);
      if (dbr->ctr[i] > 0) {
        write_to_sparse_row(A, dbr->cf[i], rbi*bs + rib, i, bs);
      }
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

/**
 * \brief After we have stored the data from one polynomial in the dense buffer,
 * we now write from the buffer to the matrix.
 *
 * \param sparse block matrix A
 *
 * \param dense block matrix B
 *
 * \param dense block row buffer dbr
 *
 * \param row block index rbi
 *
 * \param row index n block rib
 *
 * \param number of column blocks ncb
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param field characteristic mod
 */
static inline void store_in_matrix(sb_fl_t *A, dbm_fl_t *B, const dbr_t *dbr,
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
      write_to_sparse_row_in_block(A, dbr->cf[i], rbi, i, rib, dbr->ctr[i], bs, mod);
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

/**
 * \brief After we have stored the data from one polynomial in the dense buffer,
 * we now write from the buffer to the sparse part of the matrix.
 *
 * \param sparse block matrix A
 *
 * \param dense block row buffer dbr
 *
 * \param row block index rbi
 *
 * \param row index n block rib
 *
 * \param number of column blocks ncb
 *
 * \param index of first block on the righthand side in matrix fbr
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param field characteristic mod
 */
static inline void store_in_matrix_new_sparse(sb_fl_t *A, const dbr_t *dbr,
    const nelts_t rbi, const nelts_t rib, const nelts_t ncb, const nelts_t fbr,
    const bi_t bs, const coeff_t mod)
{
  nelts_t i;

  // do sparse left side A
  for (i=0; i<fbr; ++i) {
    if (dbr->ctr[i] > 0) {
      if (A->blocks[rbi][i].val == NULL)
        allocate_sparse_block(A, rbi, i, bs);
      allocate_sparse_row_in_block(A, rbi, i, rib, dbr->ctr[i], bs);
      write_to_sparse_row_in_block(A, dbr->cf[i], rbi, i, rib, dbr->ctr[i], bs, mod);
    }
  }
}

/**
 * \brief After we have stored the data from one polynomial in the dense buffer,
 * we now write from the buffer to the sparse part of the matrix.
 *
 * \param dense block matrix B
 *
 * \param dense block row buffer dbr
 *
 * \param number of column blocks ncb
 *
 * \param index of first block on the righthand side in matrix fbr
 *
 */
static inline void store_in_matrix_new_dense(dbm_fl_t *B, const dbr_t *dbr,
    const nelts_t rbi, const nelts_t ncb, const nelts_t fbr)
{
  nelts_t i;
  // do dense right side B
  for (i=0; i<ncb-fbr; ++i) {
    if (dbr->ctr[fbr+i] > 0) {
      // link buffered nonzero block into gbla matrix
      B->blocks[rbi][i].val = dbr->bl[i];
    } else {
      free(dbr->bl[i]);
    }
  }
}

/**
 * \brief Resets buffer to all entries zero once the a row is done.
 *
 * \param dense block row buffer dbr
 *
 * \param number of column blocks ncb
 *
 * \param block size bs
 */
static inline void reset_buffer(dbr_t *dbr, const nelts_t ncb, const bi_t bs)
{
  nelts_t i;

  memset(dbr->ctr, 0, ncb * sizeof(nelts_t));
  for (i=0; i<ncb; ++i)
    memset(dbr->cf[i], 0, bs * sizeof(coeff_t));
}

/**
 * \brief Stores polynomial data in dense block row which is a buffer for the
 * gbla matrix.
 *
 * \note Buffering is used in order to not allocate too much memory at the same
 * time for the gbla matrix. Also it improves the writing to the matrix since we
 * can write a full buffer at once which increases locality of data.
 *
 * \note In order to find the right blocks when buffering the coefficients we
 * have to use fr as marking point where the righthand side of the gbla matrix
 * blocks start.
 *
 * \param dense block row dbr
 *
 * \param polynomial index in basis bi
 *
 * \param polynomial index in simplifier list si
 *
 * \param hash position of multiplier mul
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param intermediate groebner basis basis
 *
 * \param simplifier list sf
 *
 * \param hash table ht
 */
static inline void store_in_buffer(dbr_t *dbr, const hash_t mul, const nelts_t nt,
    const hash_t *eh, const coeff_t *cf, const nelts_t fr, const bi_t bs,
    const gb_t *basis, const gb_t *sf, const mp_cf4_ht_t *ht)
//static inline void store_in_buffer(dbr_t *dbr, const nelts_t bi, const nelts_t si,
//    const hash_t mul, const nelts_t fr, const bi_t bs, const gb_t *basis,
//    const gb_t *sf, const mp_cf4_ht_t *ht)
//
//static inline void store_in_buffer(dbr_t *dbr, const nelts_t pi, const hash_t mul,
//    const nelts_t fr, const bi_t bs, const gb_t *basis, const mp_cf4_ht_t *ht)
{
  nelts_t j, tmp;
  // hash position and column position
  hash_t hp, cp;

  // calculate index of last block on left side
  // if there is nothing on the lefthand side what can happen when interreducing
  // the initial input elements then we have to adjust fbr to 0
  const nelts_t fbr = fr == 0 ? 0 : (fr-1)/bs + 1;

  // do some loop unrollinga
  j = 0;
  if (nt > 3) {
    for (j=0; j<nt-3; j=j+4) {
      hp  = find_in_hash_table_product(mul, eh[j], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j, basis->eh[pi][j], cp, basis->cf[pi][j]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = cf[j];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = cf[j];
        dbr->ctr[fbr+cp/bs]++;
      }

      hp  = find_in_hash_table_product(mul, eh[j+1], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+1, basis->eh[pi][j+1], cp, basis->cf[pi][j+1]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = cf[j+1];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = cf[j+1];
        dbr->ctr[fbr+cp/bs]++;
      }

      hp  = find_in_hash_table_product(mul, eh[j+2], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+2, basis->eh[pi][j+2], cp, basis->cf[pi][j+2]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = cf[j+2];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = cf[j+2];
        dbr->ctr[fbr+cp/bs]++;
      }

      hp  = find_in_hash_table_product(mul, eh[j+3], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+3, basis->eh[pi][j+3], cp, basis->cf[pi][j+3]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = cf[j+3];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = cf[j+3];
        dbr->ctr[fbr+cp/bs]++;
      }
    }
  }
  tmp = j;
  for (j=tmp; j<nt; ++j) {
    hp  = find_in_hash_table_product(mul, eh[j], ht);
    //hp  = find_in_hash_table_product(mul, basis->eh[pi][j], ht);
#if MATRIX_DEBUG
    for (int ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii]);
    printf(" ||| ");
    for (int ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[eh[j]][ii]);
    printf(" ||| ");
    for (int ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii] + ht->exp[eh[j]][ii]);
    printf(" ------------> %lu\n", hp);
    printf("fr %u | hp %u | eh[%u] = %u | cp %u | cf %u\n", fr, hp, j, eh[j], cp, cf[j]);
#endif
    cp  = ht->idx[hp];
    if (cp<fr) {
      dbr->cf[cp/bs][cp%bs]  = cf[j];
      dbr->ctr[cp/bs]++;
    } else {
      cp = cp - fr;
      dbr->cf[fbr+cp/bs][cp%bs]  = cf[j];
      dbr->ctr[fbr+cp/bs]++;
    }
  }
}

static inline void store_in_matrix_direct(sb_fl_t *A, dbm_fl_t *B, const hash_t mul, const nelts_t nt,
    const hash_t *eh, const coeff_t *cf, const nelts_t fbr, const nelts_t fr, const nelts_t rbi,  const bi_t rib, const bi_t bs,
    const coeff_t mod, const gb_t *basis, const gb_t *sf, const mp_cf4_ht_t *ht, const nelts_t i)
//static inline void store_in_buffer(dbr_t *dbr, const nelts_t bi, const nelts_t si,
//    const hash_t mul, const nelts_t fr, const bi_t bs, const gb_t *basis,
//    const gb_t *sf, const mp_cf4_ht_t *ht)
//
//static inline void store_in_buffer(dbr_t *dbr, const nelts_t pi, const hash_t mul,
//    const nelts_t fr, const bi_t bs, const gb_t *basis, const mp_cf4_ht_t *ht)
{
  int j, tmp;
  // hash position and column position
  hash_t hp, cp;

  for (j=nt-1; j>-1; --j) {
    hp  = find_in_hash_table_product(mul, eh[j], ht);
    //hp  = find_in_hash_table_product(mul, basis->eh[pi][j], ht);
#if MATRIX_DEBUG
    for (int ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii]);
    printf(" ||| ");
    for (int ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[eh[j]][ii]);
    printf(" ||| ");
    for (int ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii] + ht->exp[eh[j]][ii]);
    printf(" ------------> %lu\n", hp);
    printf("fr %u | hp %u | eh[%u] = %u | cp %u | cf %u\n", fr, hp, j, eh[j], cp, cf[j]);
#endif
    cp  = ht->idx[hp];
    if (cp<fr) {
      //printf("A thd %d writes to UPPER %p ||| %u | %u || %u %u value %6u from term %5u of poly %3u | co %u --> %lu * %lu = %lu\n", omp_get_thread_num(), A->blocks[rbi][cp/bs].val[rib], rbi, cp/bs, rib, A->blocks[rbi][cp/bs].sz[rib],cf[j],j,i, cp, mul, eh[j], hp);
      //printf("cp %u | rib %u | sz %u\n", cp, rib, A->blocks[rbi][cp/bs].sz[rib]);
      A->blocks[rbi][cp/bs].val[rib][A->blocks[rbi][cp/bs].sz[rib]] = (re_t)((re_m_t)mod - cf[j]);
      A->blocks[rbi][cp/bs].pos[rib][A->blocks[rbi][cp/bs].sz[rib]] = cp%bs;
      A->blocks[rbi][cp/bs].sz[rib]++;
    } else {
      //printf("B thd %d writes to LOWER %p ||| %u | %u || %u %u value %6u from term %5u of poly %3u | co %u --> %lu * %lu = %lu\n", omp_get_thread_num(), B->blocks[rbi][cp/bs].val+rib*bs, rbi, cp/bs, rib, rib*bs+cp%bs,cf[j],j,i, cp, mul, eh[j], hp);
      cp = cp - fr;
      B->blocks[rbi][cp/bs].val[rib*bs+cp%bs] = cf[j];
    }
  }
}

/**
 * \brief Stores polynomial data in dense block row which is a buffer for the
 * gbla matrix.
 *
 * \note Buffering is used in order to not allocate too much memory at the same
 * time for the gbla matrix. Also it improves the writing to the matrix since we
 * can write a full buffer at once which increases locality of data.
 *
 * \note In order to find the right blocks when buffering the coefficients we
 * have to use fr as marking point where the righthand side of the gbla matrix
 * blocks start.
 *
 * \param dense block row dbr
 *
 * \param polynomial index in basis bi
 *
 * \param polynomial index in simplifier list si
 *
 * \param hash position of multiplier mul
 *
 * \param index of first block on the righthand side in matrix fbr
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param intermediate groebner basis basis
 *
 * \param simplifier list sf
 *
 * \param hash table ht
 */
static inline void store_in_buffer_new(dbr_t *dbr, const bi_t rib,  const hash_t mul, const nelts_t nt,
    const hash_t *eh, const coeff_t *cf, const nelts_t fbr, const nelts_t fr, const bi_t bs,
    const gb_t *basis, const gb_t *sf, const mp_cf4_ht_t *ht)
{
  nelts_t j, tmp;
  // hash position and column position
  hash_t hp, cp;

  // do some loop unrolling
  j = 0;
  if (nt > 3) {
    for (j=0; j<nt-3; j=j+4) {
      hp  = find_in_hash_table_product(mul, eh[j], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j, basis->eh[pi][j], cp, basis->cf[pi][j]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = cf[j];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->bl[cp/bs][rib*bs+cp%bs]  = cf[j];
        dbr->ctr[fbr+cp/bs] = 1;
      }

      hp  = find_in_hash_table_product(mul, eh[j+1], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+1, basis->eh[pi][j+1], cp, basis->cf[pi][j+1]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = cf[j+1];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->bl[cp/bs][rib*bs+cp%bs]  = cf[j+1];
        dbr->ctr[fbr+cp/bs] = 1;
      }

      hp  = find_in_hash_table_product(mul, eh[j+2], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+2, basis->eh[pi][j+2], cp, basis->cf[pi][j+2]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = cf[j+2];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->bl[cp/bs][rib*bs+cp%bs]  = cf[j+2];
        dbr->ctr[fbr+cp/bs] = 1;
      }

      hp  = find_in_hash_table_product(mul, eh[j+3], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+3, basis->eh[pi][j+3], cp, basis->cf[pi][j+3]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = cf[j+3];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->bl[cp/bs][rib*bs+cp%bs]  = cf[j+3];
        dbr->ctr[fbr+cp/bs] = 1;
      }
    }
  }
  tmp = j;
  for (j=tmp; j<nt; ++j) {
    hp  = find_in_hash_table_product(mul, eh[j], ht);
    //hp  = find_in_hash_table_product(mul, basis->eh[pi][j], ht);
#if MATRIX_DEBUG
    for (int ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii]);
    printf(" ||| ");
    for (int ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[eh[j]][ii]);
    printf(" ||| ");
    for (int ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii] + ht->exp[eh[j]][ii]);
    printf(" ------------> %lu\n", hp);
    printf("fr %u | hp %u | eh[%u] = %u | cp %u | cf %u\n", fr, hp, j, eh[j], cp, cf[j]);
#endif
    cp  = ht->idx[hp];
    if (cp<fr) {
      dbr->cf[cp/bs][cp%bs]  = cf[j];
      dbr->ctr[cp/bs]++;
      //printf("%7u in %p at %5u %5u by thread %d\n", cf[j], dbr->cf[cp/bs], cp/bs, cp%bs, omp_get_thread_num());
    } else {
      cp = cp - fr;
      dbr->bl[cp/bs][rib*bs+cp%bs]  = cf[j];
      dbr->ctr[fbr+cp/bs] = 1;
      //printf("%7u in %p at %5u %5u by thread %d\n", cf[j], dbr->bl[cp/bs], cp/bs, rib*bs+cp%bs, omp_get_thread_num());
    }
  }
}

/**
 * \brief Generates one row block of gbla matrix.
 *
 * \note In order to find the right blocks when buffering the coefficients we
 * have to use fr as marking point where the righthand side of the gbla matrix
 * blocks start.
 *
 * \param sparse block matrix A
 *
 * \param dense block matrix B
 *
 * \param row block index rbi
 *
 * \param number of rows nr
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param number of column blocks (i.e. how many blocks are needed over the full
 * column range) ncb
 *
 * \param intermediate groebner basis gb
 *
 * \param simplifier list sf
 *
 * \param symbolic preprocessing selection sel
 *
 * \param symbolic preprocessing monomials col
 */
static inline void generate_row_blocks(sb_fl_t * A, dbm_fl_t *B, const nelts_t rbi,
    const nelts_t nr, const nelts_t fr, const bi_t bs, const nelts_t ncb,
    const gb_t *basis, const gb_t *sf, const sel_t *sel, const pre_t *col)
{
  nelts_t i;
  // get new row index in block rib
  bi_t rib;
  // multiplier
  hash_t mul;
  // polynomial exponent array
  hash_t *eh;
  // polynomial coefficient array
  coeff_t *cf;
  // polynomial number of terms
  nelts_t nt; // preallocate buffer to store row in dense format
  const nelts_t min = (rbi+1)*bs > nr ? nr : (rbi+1)*bs;
  // for each row we allocate memory in the sparse, left side and go through the
  // polynomials and add corresponding entries in the matrix

  dbr_t *dbr  = initialize_dense_block_row(ncb, bs);

  for (i=rbi*bs; i<min; ++i) {
    // zero out buffer data
    reset_buffer(dbr, ncb, bs);

    rib = i % bs;
    mul = sel->mpp[i].mul;
    eh  = sel->mpp[i].eh;
    cf  = sel->mpp[i].cf;
    nt  = sel->mpp[i].nt;

    store_in_buffer(dbr, mul, nt, eh, cf, fr, bs, basis, sf, ht);
#if MATRIX_DEBUG
    printf("ROW %u\n",i);
    for (int ii=0; ii<ncb; ++ii)
      for (int jj=0; jj<bs; ++jj)
        printf("%u ",dbr->cf[ii][jj]);
    printf("\n");
#endif

    store_in_matrix(A, B, dbr, rbi, rib, ncb, fr, bs, basis->mod);
  }
  free_dense_block_row(dbr, ncb);
}

/**
 * \brief Generates one row block of gbla matrix.
 *
 * \note In order to find the right blocks when buffering the coefficients we
 * have to use fr as marking point where the righthand side of the gbla matrix
 * blocks start.
 *
 * \param sparse block matrix A
 *
 * \param dense block matrix B
 *
 * \param row block index rbi
 *
 * \param number of rows nr
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param number of column blocks (i.e. how many blocks are needed over the full
 * column range) ncb
 *
 * \param intermediate groebner basis gb
 *
 * \param simplifier list sf
 *
 * \param symbolic preprocessing selection sel
 *
 * \param symbolic preprocessing monomials col
 */
static inline void generate_row_blocks_new(sb_fl_t * A, dbm_fl_t *B, const nelts_t rbi,
    const nelts_t nr, const nelts_t fr, const bi_t bs, const nelts_t ncb,
    const gb_t *basis, const gb_t *sf, const sel_t *sel, const pre_t *col)
{
  nelts_t i;
  // get new row index in block rib
  bi_t rib;
  // multiplier
  hash_t mul;
  // polynomial exponent array
  hash_t *eh;
  // polynomial coefficient array
  coeff_t *cf;
  // polynomial number of terms
  nelts_t nt; // preallocate buffer to store row in dense format
  const nelts_t min = (rbi+1)*bs > nr ? nr : (rbi+1)*bs;
  // calculate index of last block on left side
  // if there is nothing on the lefthand side what can happen when interreducing
  // the initial input elements then we have to adjust fbr to 0
  const nelts_t fbr = fr == 0 ? 0 : (fr-1)/bs + 1;
  // for each row we allocate memory in the sparse, left side and go through the
  // polynomials and add corresponding entries in the matrix

  dbr_t *dbr  = initialize_dense_block_row_new(ncb, fbr, bs);

  for (i=rbi*bs; i<min; ++i) {
    // zero out buffer data
    reset_buffer(dbr, fbr, bs);

    rib = i % bs;
    mul = sel->mpp[i].mul;
    eh  = sel->mpp[i].eh;
    cf  = sel->mpp[i].cf;
    nt  = sel->mpp[i].nt;

    store_in_buffer_new(dbr, rib, mul, nt, eh, cf, fbr, fr, bs, basis, sf, ht);
#if MATRIX_DEBUG
    printf("ROW %u\n",i);
    for (int ii=0; ii<ncb; ++ii)
      for (int jj=0; jj<bs; ++jj)
        printf("%u ",dbr->cf[ii][jj]);
    printf("\n");
#endif
    // we store the sparse part (lefthand side of the matrix per row and write
    // them row-wise
    store_in_matrix_new_sparse(A, dbr, rbi, rib, ncb, fbr, bs, basis->mod);
  }
  // for the dense matrix part we have allocated full dense blocks already, thus
  // we only link the nonzero blocks into B
  store_in_matrix_new_dense(B, dbr, rbi, ncb, fbr);
  free_dense_block_row_new(dbr, fbr);
}

static inline void generate_row_blocks_no_buffer(sb_fl_t * A, dbm_fl_t *B, const nelts_t rbi,
    const nelts_t nr, const nelts_t fr, const bi_t bs, const nelts_t ncb,
    const gb_t *basis, const gb_t *sf, const sel_t *sel, const pre_t *col)
{
  nelts_t i;
  // get new row index in block rib
  bi_t rib;
  // multiplier
  hash_t mul;
  // polynomial exponent array
  hash_t *eh;
  // polynomial coefficient array
  coeff_t *cf;
  // polynomial number of terms
  nelts_t nt; // preallocate buffer to store row in dense format
  const nelts_t min = (rbi+1)*bs > nr ? nr : (rbi+1)*bs;
  // calculate index of last block on left side
  // if there is nothing on the lefthand side what can happen when interreducing
  // the initial input elements then we have to adjust fbr to 0
  const nelts_t fbr = fr == 0 ? 0 : (fr-1)/bs + 1;
  // for each row we allocate memory in the sparse, left side and go through the
  // polynomials and add corresponding entries in the matrix

  // allocate all possible memory in matrix for this block row
  /*
  for (int k=0; k<fbr; ++k) {
    A->blocks[rbi][k].val = (re_t **)malloc(bs * sizeof(re_t *));
    A->blocks[rbi][k].pos = (bi_t **)malloc(bs * sizeof(bi_t *));
    A->blocks[rbi][k].sz  = (bi_t *)malloc(bs * sizeof(bi_t));
    for (int l=0; l<bs; ++l) {
      A->blocks[rbi][k].val[l]  = (re_t *)malloc(bs * sizeof(re_t));
      A->blocks[rbi][k].pos[l]  = (bi_t *)malloc(bs * sizeof(bi_t));
      A->blocks[rbi][k].sz[l]   = 0;
    }
  }
  for (int k=0; k<(ncb-fbr); ++k) {
    B->blocks[rbi][k].val = (re_t *)calloc(bs * bs, sizeof(re_t));
  }
  */
  for (i=rbi*bs; i<min; ++i) {

    rib = i % bs;
    mul = sel->mpp[i].mul;
    eh  = sel->mpp[i].eh;
    cf  = sel->mpp[i].cf;
    nt  = sel->mpp[i].nt;

    store_in_matrix_direct(A, B, mul, nt, eh, cf, fbr, fr, rbi, rib, bs, basis->mod, basis, sf, ht, i);
  }

  // free useless allocated memory in A and B
  /*
  int cz  = 0;
  for (int l=0; l<fbr; ++l) {
    cz  = 0;
    for (int k=0; k<bs; ++k) {
      if (A->blocks[rbi][l].sz[k] == 0) {
        cz++;
        free(A->blocks[rbi][l].val[k]);
        free(A->blocks[rbi][l].pos[k]);
        A->blocks[rbi][l].val[k]  = NULL;
        A->blocks[rbi][l].pos[k]  = NULL;
      }
    }
    if (cz == bs) {
      free(A->blocks[rbi][l].val);
      free(A->blocks[rbi][l].pos);
      free(A->blocks[rbi][l].sz);
      A->blocks[rbi][l].val = NULL;
      A->blocks[rbi][l].pos = NULL;
      A->blocks[rbi][l].sz  = NULL;
    }
  }
  coeff_t zb[bs*bs];
  memset(zb, 0, bs*bs*sizeof(coeff_t));
  for (int l=0; l<(ncb-fbr); ++l) {
    if (memcmp(B->blocks[rbi][l].val, zb, bs*bs*sizeof(coeff_t)) == 0) {
      free(B->blocks[rbi][l].val);
      B->blocks[rbi][l].val = NULL;
    }
  }
  */
}

/**
 * \brief Generates one row of gbla matrix.
 *
 * \note A is a sparse row matrix that has no block structure. In the reduction
 * process A and B will not be updated at all, thus it is not useful when also
 * applying simplify.
 *
 * \note In order to find the right blocks when buffering the coefficients we
 * have to use fr as marking point where the righthand side of the gbla matrix
 * blocks start.
 *
 * \param sparse row matrix A
 *
 * \param dense block matrix B
 *
 * \param row block index rbi
 *
 * \param number of rows nr
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param number of column blocks (i.e. how many blocks are needed over the full
 * column range) ncb
 *
 * \param intermediate groebner basis gb
 *
 * \param simplifier list sf
 *
 * \param symbolic preprocessing selection sel
 *
 * \param symbolic preprocessing monomials col
 */
static inline void generate_row_blocks_keep_A(sm_fl_t *A, dbm_fl_t *B, const nelts_t rbi,
    const nelts_t nr, const nelts_t fr, const bi_t bs, const nelts_t ncb,
    const gb_t *basis, const gb_t *sf, const sel_t *sel, const pre_t *col)
{
  nelts_t i;
  // get new row index in block rib
  bi_t rib;
  // multiplier
  hash_t mul;
  // polynomial exponent array
  hash_t *eh;
  // polynomial coefficient array
  coeff_t *cf;
  // polynomial number of terms
  nelts_t nt;
  // preallocate buffer to store row in dense format
  dbr_t *dbr  = initialize_dense_block_row(ncb, bs);
  nelts_t min = (rbi+1)*bs > nr ? nr : (rbi+1)*bs;
  //printf("rbi %u | nr %u | min %u\n",rbi,nr,min);
  // for each row we allocate memory in the sparse, left side and go through the
  // polynomials and add corresponding entries in the matrix
  for (i=rbi*bs; i<min; ++i) {
    // zero out buffer data
    reset_buffer(dbr, ncb, bs);

    rib = i % bs;
    mul = sel->mpp[i].mul;
    eh  = sel->mpp[i].eh;
    cf  = sel->mpp[i].cf;
    nt  = sel->mpp[i].nt;

    store_in_buffer(dbr, mul, nt, eh, cf, fr, bs, basis, sf, ht);
#if MATRIX_DEBUG
    printf("ROW %u\n",i);
    for (int ii=0; ii<ncb; ++ii)
      for (int jj=0; jj<bs; ++jj)
        printf("%u ",dbr->cf[ii][jj]);
    printf("\n");
#endif

    store_in_matrix_keep_A(A, B, dbr, rbi, rib, ncb, fr, bs, basis->mod);
  }
  free_dense_block_row(dbr, ncb);
}

/**
 * \brief Generates gbla matrix: enters all coefficients from the essential data
 * given by symbolic preprocessing
 *
 * \note It splices the matrix in block rows of mat->bs row size and handles
 * each block row separately. Thus, if available, the generation of the matrix
 * can be done in parallel on different threads using OpenMP.
 *
 * \param intermediate groebner basis basis
 *
 * \param simplifier list sf
 *
 * \param symbplic preprocessing data spd
 *
 * \param number of threads nthreads
 *
 * \return gbla matrix mat
 */
static inline mat_t *generate_gbla_matrix(const gb_t *basis, const gb_t *sf,
    const spd_t *spd, const int nthreads)
{
  // constructing gbla matrices is not threadsafe at the moment
  const int t = nthreads;
  mat_t *mat  = initialize_gbla_matrix(spd, basis);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    // fill the upper part AB
    for (int i=0; i<mat->rbu; ++i) {
      #pragma omp task
      {
        generate_row_blocks_new(mat->A, mat->B, i, spd->selu->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, sf, spd->selu, spd->col);
      }
    }
    // fill the lower part CD
    for (int i=0; i<mat->rbl; ++i) {
      #pragma omp task
      {
        generate_row_blocks_new(mat->C, mat->D, i, spd->sell->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, sf, spd->sell, spd->col);
      }
    }
    #pragma omp taskwait
    }
  }
  /*
  if (mat->A != NULL && mat->A->blocks != NULL) {
    for (int ii=0; ii<mat->rbu; ++ii) {
      for (int jj=0; jj<mat->cbl; ++jj) {
        if (mat->A->blocks[ii][jj].val != NULL) {
          printf("%d .. %d\n", ii, jj);
          for (int kk=0; kk<mat->bs; ++kk) {
            for (int ll=0; ll<mat->A->blocks[ii][jj].sz[kk]; ++ll) {
              printf("%d | %d || ", mat->A->blocks[ii][jj].val[kk][ll], mat->A->blocks[ii][jj].pos[kk][ll]);
            }
            printf("\n");
          }
        }
      }
    }
  }
  */
  return mat;
}

static inline mat_t *generate_gbla_matrix_test(const gb_t *basis, const gb_t *sf,
    const spd_t *spd, const int nthreads)
{
  // constructing gbla matrices is not threadsafe at the moment
  const int t = nthreads;
  mat_t *mat  = initialize_gbla_matrix(spd, basis);
  const nelts_t fbr = spd->col->nlm == 0 ? 0 : (spd->col->nlm-1)/mat->bs + 1;
  const bi_t bs = mat->bs;
  const nelts_t ncb = mat->cbl + mat->cbr;
  // allocate all possible memory in matrix for this block row
  for (int j=0; j<mat->rbu; ++j) {
    for (int k=0; k<fbr; ++k) {
      mat->A->blocks[j][k].val = (re_t **)malloc(bs * sizeof(re_t *));
      mat->A->blocks[j][k].pos = (bi_t **)malloc(bs * sizeof(bi_t *));
      mat->A->blocks[j][k].sz  = (bi_t *)malloc(bs * sizeof(bi_t));
      for (int l=0; l<bs; ++l) {
        mat->A->blocks[j][k].val[l]  = (re_t *)malloc(bs * sizeof(re_t));
        mat->A->blocks[j][k].pos[l]  = (bi_t *)malloc(bs * sizeof(bi_t));
        mat->A->blocks[j][k].sz[l]   = 0;
      }
    }
    for (int k=0; k<(ncb-fbr); ++k) {
      mat->B->blocks[j][k].val = (re_t *)calloc(bs * bs, sizeof(re_t));
    }
  }
  for (int j=0; j<mat->rbl; ++j) {
    for (int k=0; k<fbr; ++k) {
      mat->C->blocks[j][k].val = (re_t **)malloc(bs * sizeof(re_t *));
      mat->C->blocks[j][k].pos = (bi_t **)malloc(bs * sizeof(bi_t *));
      mat->C->blocks[j][k].sz  = (bi_t *)malloc(bs * sizeof(bi_t));
      for (int l=0; l<bs; ++l) {
        mat->C->blocks[j][k].val[l]  = (re_t *)malloc(bs * sizeof(re_t));
        mat->C->blocks[j][k].pos[l]  = (bi_t *)malloc(bs * sizeof(bi_t));
        mat->C->blocks[j][k].sz[l]   = 0;
      }
    }
    for (int k=0; k<(ncb-fbr); ++k) {
      mat->D->blocks[j][k].val = (re_t *)calloc(bs * bs, sizeof(re_t));
    }
  }
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    // fill the upper part AB
    for (int i=0; i<mat->rbu; ++i) {
      #pragma omp task
      {
        generate_row_blocks_no_buffer(mat->A, mat->B, i, spd->selu->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, sf, spd->selu, spd->col);
      }
    }
    // fill the lower part CD
    for (int i=0; i<mat->rbl; ++i) {
      #pragma omp task
      {
        generate_row_blocks_no_buffer(mat->C, mat->D, i, spd->sell->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, sf, spd->sell, spd->col);
      }
    }
    #pragma omp taskwait
    }
  }
  // free useless allocated memory in A and B
  for (int j=0; j<mat->rbu; ++j) {
    int cz  = 0;
    for (int l=0; l<fbr; ++l) {
      cz  = 0;
      for (int k=0; k<bs; ++k) {
        if (mat->A->blocks[j][l].sz[k] == 0) {
          cz++;
          free(mat->A->blocks[j][l].val[k]);
          free(mat->A->blocks[j][l].pos[k]);
          mat->A->blocks[j][l].val[k]  = NULL;
          mat->A->blocks[j][l].pos[k]  = NULL;
        }
      }
      if (cz == bs) {
        free(mat->A->blocks[j][l].val);
        free(mat->A->blocks[j][l].pos);
        free(mat->A->blocks[j][l].sz);
        mat->A->blocks[j][l].val = NULL;
        mat->A->blocks[j][l].pos = NULL;
        mat->A->blocks[j][l].sz  = NULL;
      }
    }
    coeff_t zb[bs*bs];
    memset(zb, 0, bs*bs*sizeof(coeff_t));
    for (int l=0; l<(ncb-fbr); ++l) {
      if (memcmp(mat->B->blocks[j][l].val, zb, bs*bs*sizeof(coeff_t)) == 0) {
        free(mat->B->blocks[j][l].val);
        mat->B->blocks[j][l].val = NULL;
      }
    }
  }
  for (int j=0; j<mat->rbl; ++j) {
    int cz  = 0;
    for (int l=0; l<fbr; ++l) {
      cz  = 0;
      for (int k=0; k<bs; ++k) {
        if (mat->C->blocks[j][l].sz[k] == 0) {
          cz++;
          free(mat->C->blocks[j][l].val[k]);
          free(mat->C->blocks[j][l].pos[k]);
          mat->C->blocks[j][l].val[k]  = NULL;
          mat->C->blocks[j][l].pos[k]  = NULL;
        }
      }
      if (cz == bs) {
        free(mat->C->blocks[j][l].val);
        free(mat->C->blocks[j][l].pos);
        free(mat->C->blocks[j][l].sz);
        mat->C->blocks[j][l].val = NULL;
        mat->C->blocks[j][l].pos = NULL;
        mat->C->blocks[j][l].sz  = NULL;
      }
    }
    coeff_t zb[bs*bs];
    memset(zb, 0, bs*bs*sizeof(coeff_t));
    for (int l=0; l<(ncb-fbr); ++l) {
      if (memcmp(mat->D->blocks[j][l].val, zb, bs*bs*sizeof(coeff_t)) == 0) {
        free(mat->D->blocks[j][l].val);
        mat->D->blocks[j][l].val = NULL;
      }
    }
  }
  /*
  if (mat->A != NULL && mat->A->blocks != NULL) {
    for (int ii=0; ii<mat->rbu; ++ii) {
      for (int jj=0; jj<mat->cbl; ++jj) {
        if (mat->A->blocks[ii][jj].val != NULL) {
          printf("%d .. %d\n", ii, jj);
          for (int kk=0; kk<mat->bs; ++kk) {
            for (int ll=0; ll<mat->A->blocks[ii][jj].sz[kk]; ++ll) {
              printf("%d | %d || ", mat->A->blocks[ii][jj].val[kk][ll], mat->A->blocks[ii][jj].pos[kk][ll]);
            }
            printf("\n");
          }
        }
      }
    }
  }
  */
  return mat;
}

/**
 * \brief Generates gbla matrix: enters all coefficients from the essential data
 * given by symbolic preprocessing. This version keeps A and is thus not useful
 * when applying simplify since we do not update A and B.
 *
 * \param intermediate groebner basis basis
 *
 * \param simplifier list sf
 *
 * \param symbplic preprocessing data spd
 *
 * \param number of threads nthreads
 *
 * \return gbla matrix mat
 */
static inline mat_t *generate_gbla_matrix_keep_A(const gb_t *basis, const gb_t *sf,
    const spd_t *spd, const int nthreads)
{
  // constructing gbla matrices is not threadsafe at the moment
  const int t = 1;
  mat_t *mat  = initialize_gbla_matrix_keep_A(spd, basis);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    // fill the upper part AB
    for (int i=0; i<mat->rbu; ++i) {
      #pragma omp task
      {
        generate_row_blocks_keep_A(mat->AR, mat->B, i, spd->selu->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, sf, spd->selu, spd->col);
      }
    }
    // fill the lower part CD
    for (int i=0; i<mat->rbl; ++i) {
      #pragma omp task
      {
        generate_row_blocks_keep_A(mat->CR, mat->D, i, spd->sell->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, sf, spd->sell, spd->col);
      }
    }
    }
    #pragma omp taskwait
  }

  return mat;
}

/**
 * \brief Reduces gbla matrix generated by symbolic preprocessing data. Stores
 * reduced D in mat->DR in dense row format. This computation is completely done
 * in gbla.
 *
 * \param gbla matrix mat
 *
 * \param verbosity level verbose
 *
 * \param number of threads nthreads
 *
 * \return rank of D, negative value if failure happens
 */
int reduce_gbla_matrix(mat_t * mat, int verbose, int nthreads);

/**
 * \brief Reduces gbla matrix generated by symbolic preprocessing data. Stores
 * reduced D in mat->DR in dense row format. This computation is completely done
 * in gbla. This version keeps A and is not useful when applying simplify
 *
 * \param gbla matrix mat
 *
 * \param verbosity level verbose
 *
 * \param number of threads nthreads
 *
 * \return rank of D, negative value if failure happens
 */
int reduce_gbla_matrix_keep_A(mat_t * mat, int verbose, int nthreads);
#endif
