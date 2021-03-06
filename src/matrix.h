/* gb: Gröbner Basis
 * Copyright (C) 2018 Christian Eder <ederc@mathematik.uni-kl.de>
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

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include "gbla-elimination.h"
#include <cli/io.h>
#include <config.h>
#include "types.h"
#include "hash.h"

#include "immintrin.h"

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#ifndef MATRIX_DEBUG
#define MATRIX_DEBUG 0
#endif

#define newred 0
#define BLOCK 10000

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
  return (ci_t) ceil((float)col->nlm/(float)bs);
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
  return (ci_t) ceil((float)(col->load - col->nlm)/(float)bs);
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
  return (ri_t) ceil((float)sel->load/(float)bs);
}

/**
 * \brief Computes inverse value of x modulo y:
 * We compute the inverse using the extended GCD. So we are only interested in x.
 * Note that internally we need signed types, but we return only unsigned type
 * re_t for x.
 *
 * \param x
 *
 * \param modulur
 */
static inline void inverse_val_new(cf_t *x, const cf_t modulus) {
  assert(*x);
  if ( *x == 1 ) return ;
  assert((int32_t)modulus > 0);
  int32_t u1 = 1, u2 = 0;
  int32_t v1 = 0, v3 = (int32_t)modulus;
  int32_t u3 = (int32_t)*x, v2 = 1;
  while (v3 != 0) {
    int32_t q  = u3 / v3;
    int32_t t1 = u1 - v1 * q;
    u1  = v1; v1  = t1;

    int32_t t3 = u3 - v3 * q;
    u3  = v3; v3  = t3;

    int32_t t2 = u2 - v2 * q;
    u2  = v2; v2  = t2;
  }
  if (u1 < 0) {
    u1  +=  (int32_t)modulus;
    /* check_inverse(*x,u1,modulus); */
    *x  =   (cf_t)u1;
    return;
  }
  if (u1 > (int32_t)modulus) {
    u1  -=  (int32_t)modulus;
    /* check_inverse(*x,u1,modulus); */
    *x  =   (cf_t) u1;
    return;
  }
  /* check_inverse(*x,u1,modulus); */
  *x  = (cf_t)u1;
  return;
}

static inline cf_t inverse_coeff(int32_t cf, const int32_t mod) {
  int32_t a, b, x, y, q, t;
  a = mod;
  b = cf % mod;
  b +=  (b >> 31) & mod;
  x =   1;
  y =   0;

  while (b != 0) {
    t = b;
    q = a / t;
    b = a - q * t;
    a = t;
    t = x;
    x = y - q * t;
    y = t;
  }
  y +=  (y >> 31) & mod;
  while (y < 0)
    y +=  mod;
  return (cf_t)y;
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

  /* initialize parts of gbla matrix with known dimensions */

  /* D exists always */
  init_dbm(mat->D, spd->sell->load, spd->col->load - spd->col->nlm);
  
  /* if no upper part A & B are NULL and so is C */
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

static inline smat_t *initialize_sparse_matrix(const nelts_t nr, const nelts_t ncl, const nelts_t ncr, cf_t mod)
{
  smat_t *mat  = (smat_t *)malloc(sizeof(smat_t));

  mat->mod  = mod;
  mat->ncl  = ncl;
  mat->ncr  = ncr;
  mat->nr   = nr;
  mat->rk   = nr;

  mat->row  = (sr_t **)malloc(mat->nr * sizeof(sr_t *));

  return mat;
}

static inline smc_t *initialize_sparse_compact_matrix(const nelts_t nr, const nelts_t ncl, const nelts_t ncr, cf_t mod)
{
  smc_t *mat  = (smc_t *)malloc(sizeof(smc_t));

  mat->mod  = mod;
  mat->ncl  = ncl;
  mat->ncr  = ncr;
  mat->nr   = nr;
  mat->rk   = nr;

  mat->row  = (src_t **)malloc(mat->nr * sizeof(src_t *));

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
static inline mat_t *initialize_gbla_matrix(
    const smc_t *AB, const smc_t *CD)
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

  mat->mod  = AB->mod;
  mat->bs   = __GBLA_SIMD_BLOCK_SIZE;
  mat->ncl  = AB->ncl;
  mat->ncr  = AB->ncr;
  mat->nru  = AB->nr;
  mat->nrl  = CD->nr;
  mat->rbu  = (ri_t) ceil((float)mat->nru / (float)mat->bs);
  mat->rbl  = (ri_t) ceil((float)mat->nrl / (float)mat->bs);
  mat->cbl  = (ci_t) ceil((float)mat->ncl / (float)mat->bs);
  mat->cbr  = (ci_t) ceil((float)mat->ncr / (float)mat->bs);
#if MATRIX_DEBUG
  printf("cbl %u | cbr %u\n",mat->cbl, mat->cbr);
#endif

  /* initialize parts of gbla matrix with known dimensions */

  /* D exists always */
  init_dbm(mat->D, mat->nrl, mat->ncr);
  
  /* if no upper part A & B are NULL and so is C */
  if (AB->rk == 0) {
    mat->A->blocks  = NULL;
    mat->A->nrows   = mat->nru;
    mat->A->ncols   = mat->ncl;
    mat->B->blocks  = NULL;
    mat->B->nrows   = mat->nru;
    mat->B->ncols   = mat->ncr;
    mat->C->blocks  = NULL;
    mat->C->nrows   = mat->nrl;
    mat->C->ncols   = mat->ncl;
  } else {
    init_sb(mat->A, mat->nru, mat->ncl);
    init_sb(mat->C, mat->nrl, mat->ncl);
    init_dbm(mat->B, mat->nru, mat->ncr);
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

  /* A, C and D are already freed, just check again */
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
  /* mat->CR already freed after copying it to sparse block matrix mat->C */
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
  
  /* B is dense block matrix */
  if (mat->B != NULL) {
    if (mat->B->blocks != NULL) {
      free_dense_submatrix(&(mat->B), 1);
    } else {
      free(mat->B);
      mat->B  = NULL;
    }
  }
  if (mat->BR != NULL) {
    /* DR is a dense row matrix */
    free_dense_row_submatrix(&(mat->BR), 1);
  }

  if (mat->DR != NULL) {
    /* DR is a dense row matrix */
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
  dbr->cf     = (re_t **)malloc(nb * sizeof(re_t *));
  /* allocate buffer for sparse rows */
  for (i=0; i<nb; ++i)
    dbr->cf[i]  = (re_t *)calloc(bs, sizeof(re_t));
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
    dbr->cf     = (re_t **)malloc(fbr * sizeof(re_t *));
    /* allocate buffer for sparse rows */
    for (i=0; i<fbr; ++i)
      dbr->cf[i]  = (re_t *)calloc(bs, sizeof(re_t));
  } else {
    dbr->cf = NULL;
  }
  dbr->bl     = (re_t **)malloc((nb-fbr) * sizeof(re_t *));
  /* allocate buffer for dense blocks, those will be the final blocks in the
   * matrix */
  for (i=0; i<nb-fbr; ++i)
    dbr->bl[i]  = (re_t *)calloc(bs*bs, sizeof(re_t));

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
static inline void write_to_sparse_row(sm_fl_t *A, const re_t *cf, const nelts_t ri,
    const nelts_t bir, const bi_t bs)
{
  bi_t i;
  /* for (i=bs; i>0; --i) { */
  for (i=0; i<bs; ++i) {
    /* if (cf[i-1] != 0) { */
    if (cf[i] != 0) {
      A->row[ri][A->sz[ri]] = (re_t)cf[i];
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
static inline void write_to_sparse_row_in_block(sb_fl_t *A, const re_t *cf, const nelts_t rbi,
    const nelts_t bir, const nelts_t rib, const nelts_t sz,
    const bi_t bs, const cf_t mod)
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
static inline void write_to_dense_row(dbm_fl_t *A, const re_t *cf, const nelts_t rbi,
    const nelts_t bir, const nelts_t rib, const bi_t bs)
{
#if 1
    memcpy(A->blocks[rbi][bir].val+(rib*bs), cf, bs*sizeof(cf_t));
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
 */ 
static inline void allocate_sparse_row_in_block(sb_fl_t *A, const nelts_t rbi,
    const nelts_t bir, const bi_t rib, const nelts_t sz)
{
  A->blocks[rbi][bir].val[rib]  = (re_t *)malloc(sz * sizeof(re_t));
  A->blocks[rbi][bir].pos[rib]  = (bi_t *)malloc(sz * sizeof(bi_t));
#if MATRIX_DEBUG
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
 */
static inline void store_in_matrix_keep_A(sm_fl_t *A, dbm_fl_t *B, const dbr_t *dbr,
    const nelts_t rbi, const nelts_t rib, const nelts_t ncb, const nelts_t fr,
    const bi_t bs)
{
  nelts_t i;

  /* calculate index of last block on left side
   * if there is nothing on the lefthand side what can happen when interreducing
   * the initial input elements then we have to adjust fbr to 0 */
  const nelts_t fbr = fr == 0 ? 0 : (fr-1)/bs + 1;

  nelts_t nl  = 0;
  /* get number of elements in sparse row in A */
  for (i=0; i<fbr; ++i)
    nl  +=  dbr->ctr[i];

  if (nl > 0) {
    /* allocate memory for rows in A */
    allocate_sparse_row(A, rbi*bs + rib, nl);
    /* do sparse left side A */
    for (i=0; i<fbr; ++i) {
      /* printf("dbr->ctr[%u] = %u\n",i,dbr->ctr[i]); */
      if (dbr->ctr[i] > 0) {
        write_to_sparse_row(A, dbr->cf[i], rbi*bs + rib, i, bs);
      }
    }
  }

  /* do dense right side B */
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
    const bi_t bs, const cf_t mod)
{
  nelts_t i;

  /* calculate index of last block on left side
   * if there is nothing on the lefthand side what can happen when interreducing
   * the initial input elements then we have to adjust fbr to 0 */
  const nelts_t fbr = fr == 0 ? 0 : (fr-1)/bs + 1;

  /* do sparse left side A */
  for (i=0; i<fbr; ++i) {
    /* printf("dbr->ctr[%u] = %u\n",i,dbr->ctr[i]); */
    if (dbr->ctr[i] > 0) {
      if (A->blocks[rbi][i].val == NULL)
        allocate_sparse_block(A, rbi, i, bs);
      allocate_sparse_row_in_block(A, rbi, i, rib, dbr->ctr[i]);
      write_to_sparse_row_in_block(A, dbr->cf[i], rbi, i, rib, dbr->ctr[i], bs, mod);
    }
  }

  /* do dense right side B */
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
 * \param index of first block on the righthand side in matrix fbr
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param field characteristic mod
 */
static inline void store_in_matrix_new_sparse(sb_fl_t *A, const dbr_t *dbr,
    const nelts_t rbi, const nelts_t rib, const nelts_t fbr, const bi_t bs,
    const cf_t mod)
{
  nelts_t i;

  /* do sparse left side A */
  for (i=0; i<fbr; ++i) {
    if (dbr->ctr[i] > 0) {
      if (A->blocks[rbi][i].val == NULL)
        allocate_sparse_block(A, rbi, i, bs);
      allocate_sparse_row_in_block(A, rbi, i, rib, dbr->ctr[i]);
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
  /* do dense right side B */
  for (i=0; i<ncb-fbr; ++i) {
    if (dbr->ctr[fbr+i] > 0) {
      /* link buffered nonzero block into gbla matrix */
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
    memset(dbr->cf[i], 0, bs * sizeof(re_t));
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
 * \param hash table ht
 */
static inline void store_in_buffer(dbr_t *dbr, const hash_t mul, const nelts_t nt,
    const hash_t *eh, const cf_t *cf, const nelts_t fr, const bi_t bs,
    const ht_t *ht)
{
  nelts_t j, tmp;
  /* hash position and column position */
  hash_t hp, cp;

  /* calculate index of last block on left side
   * if there is nothing on the lefthand side what can happen when interreducing
   * the initial input elements then we have to adjust fbr to 0 */
  const nelts_t fbr = fr == 0 ? 0 : (fr-1)/bs + 1;

  /* do some loop unrollinga */
  j = 0;
  if (nt > 3) {
    for (j=0; j<nt-3; j=j+4) {
      hp  = find_in_hash_table_product(mul, eh[j], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j, basis->eh[pi][j], cp, basis->cf[pi][j]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = (re_t)cf[j];
        dbr->ctr[fbr+cp/bs]++;
      }

      hp  = find_in_hash_table_product(mul, eh[j+1], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+1, basis->eh[pi][j+1], cp, basis->cf[pi][j+1]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j+1];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = (re_t)cf[j+1];
        dbr->ctr[fbr+cp/bs]++;
      }

      hp  = find_in_hash_table_product(mul, eh[j+2], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+2, basis->eh[pi][j+2], cp, basis->cf[pi][j+2]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j+2];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = (re_t)cf[j+2];
        dbr->ctr[fbr+cp/bs]++;
      }

      hp  = find_in_hash_table_product(mul, eh[j+3], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+3, basis->eh[pi][j+3], cp, basis->cf[pi][j+3]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j+3];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = (re_t)cf[j+3];
        dbr->ctr[fbr+cp/bs]++;
      }
    }
  }
  tmp = j;
  for (j=tmp; j<nt; ++j) {
    hp  = find_in_hash_table_product(mul, eh[j], ht);
    /* hp  = find_in_hash_table_product(mul, basis->eh[pi][j], ht); */
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
      dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j];
      dbr->ctr[cp/bs]++;
    } else {
      cp = cp - fr;
      dbr->cf[fbr+cp/bs][cp%bs]  = (re_t)cf[j];
      dbr->ctr[fbr+cp/bs]++;
    }
  }
}

static inline void store_in_matrix_direct(sb_fl_t *A, dbm_fl_t *B, const hash_t mul, const nelts_t nt,
    const hash_t *eh, const cf_t *cf, const nelts_t fr, const nelts_t rbi,  const bi_t rib, const bi_t bs,
    const cf_t mod, const ht_t *ht)
{
  int j = (int)(nt-1);
  /* hash position and column position */
  hash_t hp, cp;

  for (; j>-1; --j) {
    hp  = find_in_hash_table_product(mul, eh[j], ht);
    /* hp  = find_in_hash_table_product(mul, basis->eh[pi][j], ht); */
#if MATRIX_DEBUG
    for (nelts_t ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii]);
    printf(" ||| ");
    for (nelts_t ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[eh[j]][ii]);
    printf(" ||| ");
    for (nelts_t ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii] + ht->exp[eh[j]][ii]);
    printf(" ------------> %lu\n", hp);
    printf("fr %u | hp %u | eh[%u] = %u | cp %u | cf %u\n", fr, hp, j, eh[j], cp, cf[j]);
#endif
    cp  = ht->idx[hp];
    if (cp<fr) {
      /* printf("A thd %d writes to UPPER %p ||| %u | %u || %u %u value %6u from term %5u of poly %3u | co %u --> %lu * %lu = %lu\n", omp_get_thread_num(), A->blocks[rbi][cp/bs].val[rib], rbi, cp/bs, rib, A->blocks[rbi][cp/bs].sz[rib],cf[j],j,i, cp, mul, eh[j], hp); */
      /* printf("cp %u | rib %u | sz %u\n", cp, rib, A->blocks[rbi][cp/bs].sz[rib]); */
      A->blocks[rbi][cp/bs].val[rib][A->blocks[rbi][cp/bs].sz[rib]] = (re_t)((re_m_t)mod - cf[j]);
      A->blocks[rbi][cp/bs].pos[rib][A->blocks[rbi][cp/bs].sz[rib]] = (bi_t)(cp%bs);
      A->blocks[rbi][cp/bs].sz[rib]++;
    } else {
      /* printf("B thd %d writes to LOWER %p ||| %u | %u || %u %u value %6u from term %5u of poly %3u | co %u --> %lu * %lu = %lu\n", omp_get_thread_num(), B->blocks[rbi][cp/bs].val+rib*bs, rbi, cp/bs, rib, rib*bs+cp%bs,cf[j],j,i, cp, mul, eh[j], hp); */
      cp = cp - fr;
      B->blocks[rbi][cp/bs].val[rib*bs+cp%bs] = (re_t)cf[j];
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
 * \param hash table ht
 */
static inline void store_in_buffer_new(dbr_t *dbr, const bi_t rib,  const hash_t mul, const nelts_t nt,
    const hash_t *eh, const cf_t *cf, const nelts_t fbr, const nelts_t fr, const bi_t bs,
    const ht_t *ht)
{
  nelts_t j, tmp;
  /* hash position and column position */
  hash_t hp, cp;

  /* do some loop unrolling */
  j = 0;
  if (nt > 3) {
    for (j=0; j<nt-3; j=j+4) {
      hp  = find_in_hash_table_product(mul, eh[j], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j, basis->eh[pi][j], cp, basis->cf[pi][j]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->bl[cp/bs][rib*bs+cp%bs]  = (re_t)cf[j];
        dbr->ctr[fbr+cp/bs] = 1;
      }

      hp  = find_in_hash_table_product(mul, eh[j+1], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+1, basis->eh[pi][j+1], cp, basis->cf[pi][j+1]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j+1];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->bl[cp/bs][rib*bs+cp%bs]  = (re_t)cf[j+1];
        dbr->ctr[fbr+cp/bs] = 1;
      }

      hp  = find_in_hash_table_product(mul, eh[j+2], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+2, basis->eh[pi][j+2], cp, basis->cf[pi][j+2]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j+2];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->bl[cp/bs][rib*bs+cp%bs]  = (re_t)cf[j+2];
        dbr->ctr[fbr+cp/bs] = 1;
      }

      hp  = find_in_hash_table_product(mul, eh[j+3], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+3, basis->eh[pi][j+3], cp, basis->cf[pi][j+3]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j+3];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->bl[cp/bs][rib*bs+cp%bs]  = (re_t)cf[j+3];
        dbr->ctr[fbr+cp/bs] = 1;
      }
    }
  }
  tmp = j;
  for (j=tmp; j<nt; ++j) {
    hp  = find_in_hash_table_product(mul, eh[j], ht);
    /* hp  = find_in_hash_table_product(mul, basis->eh[pi][j], ht); */
#if MATRIX_DEBUG
    for (nelts_t ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii]);
    printf(" ||| ");
    for (nelts_t ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[eh[j]][ii]);
    printf(" ||| ");
    for (nelts_t ii=0; ii<basis->nv; ++ii)
      printf("%u ",ht->exp[mul][ii] + ht->exp[eh[j]][ii]);
    printf(" ------------> %lu\n", hp);
    printf("fr %u | hp %u | eh[%u] = %u | cp %u | cf %u\n", fr, hp, j, eh[j], cp, cf[j]);
#endif
    cp  = ht->idx[hp];
    if (cp<fr) {
      dbr->cf[cp/bs][cp%bs]  = (re_t)cf[j];
      dbr->ctr[cp/bs]++;
      /* printf("%7u in %p at %5u %5u by thread %d\n", cf[j], dbr->cf[cp/bs], cp/bs, cp%bs, omp_get_thread_num()); */
    } else {
      cp = cp - fr;
      dbr->bl[cp/bs][rib*bs+cp%bs]  = (re_t)cf[j];
      dbr->ctr[fbr+cp/bs] = 1;
      /* printf("%7u in %p at %5u %5u by thread %d\n", cf[j], dbr->bl[cp/bs], cp/bs, rib*bs+cp%bs, omp_get_thread_num()); */
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
* \param symbolic preprocessing selection sel
*/
static inline void generate_gbla_row_blocks(sb_fl_t *L,
    dbm_fl_t *R, const mat_t *GM, const nelts_t rbi, const smc_t *M)
{
  const nelts_t min = (rbi+1) * GM->bs > M->rk ? M->rk : (rbi+1) * GM->bs;
  const nelts_t bs2 = GM->bs * GM->bs;
  size_t j    = 0;
  nelts_t ctr = 0;

  /* allocate all possible dense blocks on the righthand side */
  for (size_t i = 0; i < GM->cbr; ++i) {
    R->blocks[rbi][i].val = (re_t *)malloc(bs2 * sizeof(re_t));
  }

  /* loop over the rows in this block */
  for (size_t i = rbi * GM->bs; i < min; ++i) {
    /* allocate all possible memory for the sparse rows */
    for (size_t j = 0; j < GM->cbl; ++j) {
      L->blocks[rbi][j].val[i]  = (re_t *)malloc(GM->bs * sizeof(re_t));
      L->blocks[rbi][j].pos[i]  = (bi_t *)malloc(GM->bs * sizeof(bi_t));
    }

    ctr = 0;
    j   = 2;

    /* fill sparse lefthand side */
    while (j < M->row[i][0] && M->row[i][j] < GM->ncl) {
      ctr = L->blocks[rbi][j/GM->bs].sz[i];
      L->blocks[rbi][j/GM->bs].pos[i][ctr]  = M->row[i][j] % GM->bs;
      L->blocks[rbi][j/GM->bs].val[i][ctr]  = (re_t)M->row[i][j+1];
      L->blocks[rbi][j/GM->bs].sz[i]++;
      j += 2;
    }
    
    /* fill dense righthand side */
    while (j < M->row[i][0]) {
      R->blocks[rbi][j/GM->bs - GM->cbl].val[GM->bs * i + j % GM->bs] =
        (re_t)M->row[i][j+1];
    }

    /* free useless sparse memory */
    for (j = 0; j < GM->cbl; ++j) {
      if (L->blocks[rbi][j].sz[i] == 0) {
        free(L->blocks[rbi][j].val[i]);
        L->blocks[rbi][j].val[i]  = NULL;
        free(L->blocks[rbi][j].pos[i]); 
        L->blocks[rbi][j].pos[i]  = NULL; 
      }
    }
  }

  /* free useless dense memory */
  for (size_t i = 0; i < GM->cbr; ++i) {
    for (j = 0; j < bs2; ++j) {
      if (R->blocks[rbi][i].val[j] != 0) {
        break;
      }
    }
    if (j == bs2) {
      free(R->blocks[rbi][i].val);
      R->blocks[rbi][i].val = NULL;
    }
  }
}

static inline int cmp_src_tmp(const void *a, const void *b)
{
  const src_tmp_t ra  = *((src_tmp_t *)a);
  const src_tmp_t rb  = *((src_tmp_t *)b);

  return (int)ra.pos - (int)rb.pos;
}

static inline void set_column_indices(smc_t *M, const ht_t *ht)
{
  for (size_t i = 0; i < M->nr; ++i) {
    for (size_t j = 2; j < M->row[i][1]; j = j+2) {
      M->row[i][j]  = ht->idx[M->row[i][j]];
    }
  }
}

static inline int cmp_by_column_indices(const void *a, const void *b)
{
  const src_t pa = *((src_t *)a);
  const src_t pb = *((src_t *)b);

  return (int)(pa - pb);
}

static inline void sort_rows_by_column_indices(smc_t *M)
{
  for (size_t i = 0; i < M->rk; ++i)
    qsort(M->row[i]+2, (M->row[i][1]-2)/2, 2*sizeof(src_t),
        cmp_by_column_indices);
}



static inline int cmp_rows_by_decreasing_lm_drl(const void *a, const void *b)
{
  const src_t *ra = *((src_t **)a);
  const src_t *rb = *((src_t **)b);

  if (ht->deg[ra[2]] != ht->deg[rb[2]])
    return (int)(ht->deg[ra[2]]-ht->deg[rb[2]]);

  return memcmp(ht->exp[rb[2]], ht->exp[ra[2]], sizeof(exp_t) * ht->nv);
}

static inline void sort_rows_by_decreasing_lm_drl(smc_t *M)
{
  qsort(M->row, M->rk, sizeof(src_t *), cmp_rows_by_decreasing_lm_drl);
}



static inline int cmp_rows_by_decreasing_lm(const void *a, const void *b)
{
  const src_t *ra = *((src_t **)a);
  const src_t *rb = *((src_t **)b);

  return (int)(ra[2] - rb[2]);
}

static inline int cmp_rows_by_increasing_lm(const void *a, const void *b)
{
  const src_t *ra = *((src_t **)a);
  const src_t *rb = *((src_t **)b);

  return (int)(rb[2] - ra[2]);
}

static inline int cmp_rows_by_increasing_lm_drl(const void *a, const void *b)
{
  const src_t *ra = *((src_t **)a);
  const src_t *rb = *((src_t **)b);

  if (ht->deg[rb[2]] != ht->deg[ra[2]])
    return (int)(ht->deg[rb[2]]-ht->deg[ra[2]]);

  return memcmp(ht->exp[ra[2]], ht->exp[rb[2]], sizeof(exp_t) * ht->nv);
}

static inline void sort_rows_by_increasing_lm_drl(smc_t *M)
{
  qsort(M->row, M->rk, sizeof(src_t *), cmp_rows_by_increasing_lm_drl);
}



static inline int cmp_rows_by_decreasing_lm_lex(const void *a, const void *b)
{
  const src_t *ra = *((src_t **)a);
  const src_t *rb = *((src_t **)b);
  return memcmp(ht->exp[rb[2]], ht->exp[ra[2]], sizeof(exp_t) * ht->nv);
}

static inline void sort_rows_by_decreasing_lm_lex(smc_t *M)
{
  qsort(M->row, M->rk, sizeof(src_t *), cmp_rows_by_decreasing_lm_lex);
}



static inline int cmp_rows_by_increasing_lm_lex(const void *a, const void *b)
{
  const src_t *ra = *((src_t **)a);
  const src_t *rb = *((src_t **)b);
  return memcmp(ht->exp[ra[2]], ht->exp[rb[2]], sizeof(exp_t) * ht->nv);
}

static inline void sort_rows_by_increasing_lm_lex(smc_t *AB)
{
  qsort(AB->row, AB->rk, sizeof(src_t *), cmp_rows_by_increasing_lm_lex);
}

#define OLD_POLY_REPRESENTATION 0
#if OLD_POLY_REPRESENTATION
static inline void poly_to_sparse_compact_matrix_row_true_columns(const mpp_t *mpp,
    const nelts_t idx, const gb_t *basis, src_t **rows)
{
  const nelts_t nt  = basis->nt[mpp->bi];
  const cf_t *cf    = basis->cf[mpp->bi];
  const hash_t *eh  = basis->eh[mpp->bi];
  rows[idx]         = (src_t *)malloc((2*nt+1) * sizeof(src_t));
  rows[idx][0]      = 2*nt+1;

#if newred
  printf("corresponding polynomial in basis = %u\n", rows[idx][0]);
#endif
  nelts_t ctr = 1;
  for (nelts_t i=0; i<nt; ++i) {
    rows[idx][ctr++]  = ht->idx[find_in_hash_table_product(mpp->mul, eh[i], ht)];
    rows[idx][ctr++]  = cf[i];
  }
}

static inline void poly_to_sparse_compact_matrix_row_offset_true_columns(const mpp_t *mpp,
    const gb_t *basis, src_t **rows)
{
  const nelts_t nt  = basis->nt[mpp->bi];
  const hash_t *eh  = basis->eh[mpp->bi];
  const nelts_t idx = ht->idx[find_in_hash_table_product(mpp->mul, eh[0], ht)];
  rows[idx]         = (src_t *)malloc((nt+2) * sizeof(src_t));
  rows[idx][0]      = mpp->bi;
  rows[idx][1]      = nt % 2 + 2; /* stores offset mod 4 */

#if newred
  printf("corresponding polynomial in basis = %u\n", row[0]);
#endif
  for (nelts_t i=0; i<nt; ++i) {
    rows[idx][i+2]  = ht->idx[find_in_hash_table_product(mpp->mul, eh[i], ht)];
  }
}

static inline void poly_to_sparse_compact_matrix_row_offset_test(const mpp_t *mpp_start,
   const nelts_t bs, const gb_t *basis, src_t **rows)
{
  for (nelts_t k=0; k<bs; ++k) {
    const mpp_t *mpp  = mpp_start+k;
    const nelts_t nt  = basis->nt[mpp->bi];
    const cf_t *cf    = basis->cf[mpp->bi];
    const hash_t *eh  = basis->eh[mpp->bi];
    src_tmp_t *tmp    = (src_tmp_t *)malloc(nt * sizeof(src_tmp_t));
    rows[k]           = (src_t *)malloc((2*nt+2) * sizeof(src_t));
    rows[k][0]        = 2*nt+2;
    rows[k][1]        = nt % 2 + 1; /* stores offset mod 4 */

    nelts_t cp; /* column position */

    for (nelts_t i=0; i<nt; ++i) {
      cp  = ht->idx[find_in_hash_table_product(mpp->mul, eh[i], ht)];
      tmp[i].pos  = cp;
      tmp[i].cf   = cf[i];
    }
    /* sort temporary src_t holder */
    qsort(tmp, nt, sizeof(src_tmp_t *), cmp_src_tmp);
    nelts_t ctr = 2;
    for (nelts_t i=0; i<nt; ++i) {
      rows[k][ctr++]  = tmp[i].pos;
      rows[k][ctr++]  = tmp[i].cf;
    }
    free(tmp);
  }
}

static inline void poly_to_sparse_compact_matrix_row(const mpp_t *mpp_start,
    const nelts_t ncl, const nelts_t nc, const nelts_t bs, const gb_t *basis,
    src_t **rows)
{
  src_t *bf = (src_t *)malloc(nc * sizeof(src_t));
  nelts_t ner;

  for (nelts_t k=0; k<bs; ++k) {
    ner = 0;
    memset(bf, 0, nc*sizeof(src_t));
    const mpp_t *mpp  = mpp_start+k;
    const nelts_t nt  = basis->nt[mpp->bi];
    const cf_t *cf    = basis->cf[mpp->bi];
    const hash_t *eh  = basis->eh[mpp->bi];
    /* printf("nt[%u] = %u\n", k, mpp->nt); */
    rows[k]           = (src_t *)malloc((2*nt+1) * sizeof(src_t));
    rows[k][0]        = 2*nt+1;

    nelts_t cp; /* column position */

#if newred
    printf("nc %u\n", nc);
#endif

#if newred
    printf("row size = %u\n", row[0]);
#endif
    for (nelts_t i=0; i<nt; ++i) {
      cp  = ht->idx[find_in_hash_table_product(mpp->mul, eh[i], ht)];
      if (cp>ncl)
        ner++;
      bf[cp]  = cf[i];
      /* printf("cp %u | bf[%u] = %u for term %u | hpos %lu\n", cp, cp, bf[cp], i,find_in_hash_table_product(mpp->mul, mpp->eh[i], ht)); */
    }
    nelts_t ctr = 1;
    for (nelts_t i=0; i<nc; ++i) {
      if (bf[i] != 0) {
        /* printf("write bf[%u] = %u to position %u || cp0 %u\n", i, bf[i], ctr, cp0); */
        rows[k][ctr++]  = i;
        rows[k][ctr++]  = bf[i];
      }
    }
    /* printf("ctr %u == %u rows[%u][0] ?\n", ctr, rows[k][0], k); */
#if newred
    printf("\n");
    printf("row[0] ===== %u\n", rows[k][0]);
    for (int ii=1; ii<rows[k][0]; ii += 2)
      printf("%u ! %u -- ", rows[k][ii+1], rows[k][ii]);
    printf("\n");
#endif
  }
  free(bf);
}

static inline sr_t *poly_to_sparse_matrix_row(const mpp_t *mpp, const nelts_t nc,
    const gb_t *basis)
{
  sr_t *row = (sr_t *)malloc(sizeof(sr_t));
  const nelts_t nt  = basis->nt[mpp->bi];
  const cf_t *cf    = basis->cf[mpp->bi];
  const hash_t *eh  = basis->eh[mpp->bi];
  row->sz   = nt;
  row->val  = (cf_t *)malloc(row->sz * sizeof(cf_t));
  row->pos  = (nelts_t *)malloc(row->sz * sizeof(nelts_t));

  nelts_t cp; /* column position */

#if newred
  printf("nc %u\n", nc);
#endif
  cf_t *bf = (cf_t *)calloc(nc, sizeof(cf_t));
#if newred
  printf("row size = %u\n", row->sz);
#endif
  for (nelts_t i=0; i<nt; ++i) {
    cp  = ht->idx[find_in_hash_table_product(mpp->mul, eh[i], ht)];
    bf[cp]  = cf[i];
  }
  nelts_t ctr = 0;
  for (nelts_t i=0; i<nc; ++i) {
    if (bf[i] != 0) {
      row->pos[ctr]   = i;
      row->val[ctr] = bf[i];
      ctr++;
    }
  }

  free(bf);
#if newred
  printf("\n");
  printf("length ===== %u\n", row->sz);
#endif
  return row;
}

static inline void generate_row_blocks_no_buffer(sb_fl_t * A, dbm_fl_t *B,
    const nelts_t rbi, const nelts_t nr, const nelts_t fr, const bi_t bs,
    const gb_t *basis, const gb_t *sf,const sel_t *sel)
{
  nelts_t i;
  /* get new row index in block rib */
  bi_t rib;
  /* multiplier */
  hash_t mul;
  /* polynomial exponent array */
  hash_t *eh;
  /* polynomial coefficient array */
  cf_t *cf;
  /* polynomial number of terms */
  nelts_t nt; /* preallocate buffer to store row in dense format */
  const nelts_t min = (rbi+1)*bs > nr ? nr : (rbi+1)*bs;
  /* for each row we allocate memory in the sparse, left side and go through the
   * polynomials and add corresponding entries in the matrix */

  /* allocate all possible memory in matrix for this block row */
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
    if (sel->mpp[i].sf == 0) {
      eh  = basis->eh[sel->mpp[i].bi];
      cf  = basis->cf[sel->mpp[i].bi];
      nt  = basis->nt[sel->mpp[i].bi];
    } else {
      eh  = sf->eh[sel->mpp[i].sf];
      cf  = sf->cf[sel->mpp[i].sf];
      nt  = sf->nt[sel->mpp[i].sf];
    }

    store_in_matrix_direct(A, B, mul, nt, eh, cf, fr, rbi, rib, bs, basis->mod, ht);
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
  cf_t zb[bs*bs];
  memset(zb, 0, bs*bs*sizeof(cf_t));
  for (int l=0; l<(ncb-fbr); ++l) {
    if (memcmp(B->blocks[rbi][l].val, zb, bs*bs*sizeof(cf_t)) == 0) {
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
    const gb_t *basis, const gb_t *sf, const sel_t *sel)
{
  nelts_t i;
  /* get new row index in block rib */
  bi_t rib;
  /* multiplier */
  hash_t mul;
  /* polynomial exponent array */
  hash_t *eh;
  /* polynomial coefficient array */
  cf_t *cf;
  /* polynomial number of terms */
  nelts_t nt;
  /* preallocate buffer to store row in dense format */
  dbr_t *dbr  = initialize_dense_block_row(ncb, bs);
  nelts_t min = (rbi+1)*bs > nr ? nr : (rbi+1)*bs;
  /* printf("rbi %u | nr %u | min %u\n",rbi,nr,min); */
  /* for each row we allocate memory in the sparse, left side and go through the
   * polynomials and add corresponding entries in the matrix */
  for (i=rbi*bs; i<min; ++i) {
    /* zero out buffer data */
    reset_buffer(dbr, ncb, bs);

    rib = i % bs;
    mul = sel->mpp[i].mul;
    if (sel->mpp[i].sf == 0) {
      eh  = basis->eh[sel->mpp[i].bi];
      cf  = basis->cf[sel->mpp[i].bi];
      nt  = basis->nt[sel->mpp[i].bi];
    } else {
      eh  = sf->eh[sel->mpp[i].sf];
      cf  = sf->cf[sel->mpp[i].sf];
      nt  = sf->nt[sel->mpp[i].sf];
    }

    store_in_buffer(dbr, mul, nt, eh, cf, fr, bs, ht);
#if MATRIX_DEBUG
    printf("ROW %u\n",i);
    for (int ii=0; ii<ncb; ++ii)
      for (int jj=0; jj<bs; ++jj)
        printf("%u ",dbr->cf[ii][jj]);
    printf("\n");
#endif

    store_in_matrix_keep_A(A, B, dbr, rbi, rib, ncb, fr, bs);
  }
  free_dense_block_row(dbr, ncb);
}
#endif

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
static inline mat_t *generate_gbla_matrix(const smc_t *AB,
    const smc_t *CD, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;
  mat_t *mat  = initialize_gbla_matrix(AB, CD);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->rbu; ++i) {
      #pragma omp task
      {
        generate_gbla_row_blocks(mat->A, mat->B, mat, i, AB);
      }
    }
    /* fill the lower part CD */
    for (nelts_t i=0; i<mat->rbl; ++i) {
      #pragma omp task
      {
        generate_gbla_row_blocks(mat->C, mat->D, mat, i, CD);
      }
    }
    #pragma omp taskwait
    }
  }
   return mat;
}

#if OLD_POLY_REPRESENTATION
static inline smc_t *generate_sparse_compact_multiline_matrix(const gb_t *basis,
    const sel_t *sel, const nelts_t nr, const nelts_t ncl,
    const nelts_t ncr, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smc_t *mat  = initialize_sparse_compact_matrix(nr/2+nr%2, ncl, ncr, basis->mod);
  /* printf("nr %u --> nr/2 %u == mat->nr %u\n", nr, nr/2+nr%2, mat->nr); */
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<nr; i=i+BLOCK) {
      #pragma omp task
      {
        const nelts_t bs  = BLOCK > nr-i ? nr-i : BLOCK;
        /* printf("mat->nr %u | i %u | bs %u\n", mat->nr, i, bs); */
        /* printf("constructs row %u\n",i); */
        poly_to_sparse_compact_matrix_multiline_row(sel->mpp+i, mat->ncl + mat->ncr,
            bs, basis, mat->row+i/2);
      }
    }
    }
  }
  /* mat->nr = mat->nr/2 + mat->nr%2; */
  return mat;
}

static inline smc_t *generate_sparse_compact_matrix_pos_val(const gb_t *basis,
    const sel_t *sel, const nelts_t nr, const nelts_t ncl,
    const nelts_t ncr, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smc_t *mat  = initialize_sparse_compact_matrix(nr, ncl, ncr, basis->mod);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->nr; i=i+BLOCK) {
      #pragma omp task
      {
        const nelts_t bs  = BLOCK > mat->nr-i ? mat->nr-i : BLOCK;
        /* printf("mat->nr %u | i %u | bs %u\n", mat->nr, i, bs); */
        /* printf("constructs row %u\n",i); */
        poly_to_sparse_compact_matrix_row_pos_val(sel->mpp+i, mat->ncl + mat->ncr,
            bs, basis, mat->row+i);
      }
    }
    }
  }
  return mat;
}

static inline void add_poly_to_block(src_t **pivs, const size_t ri,
    const nelts_t idx, const gb_t *basis, const sel_t* sel)
{
  const mpp_t *mpp  = sel->mpp+ri;
  const nelts_t nt  = basis->nt[mpp->bi];
  const hash_t *eh  = basis->eh[mpp->bi];
  const cf_t *cf    = basis->cf[mpp->bi];

  pivs[idx]         = (src_t *)malloc((2*nt+2) * sizeof(src_t));

  pivs[idx][0]      = mpp->bi;
  pivs[idx][1]      = 2 * nt + 2; /* bi + length + pos/val */

  for (nelts_t i = 0; i < nt; ++i) {
    pivs[idx][2*i+2]  = ht->idx[find_in_hash_table_product(mpp->mul, eh[i], ht)];
    pivs[idx][2*i+3]  = cf[i];
  }
}

static inline void add_poly_to_pivs(src_t **pivs, const size_t ri,
    const gb_t *basis, const sel_t* sel)
{
  const mpp_t *mpp  = sel->mpp+ri;
  const nelts_t nt  = basis->nt[mpp->bi];
  const hash_t *eh  = basis->eh[mpp->bi];
  const cf_t *cf    = basis->cf[mpp->bi];
  const nelts_t idx = ht->idx[find_in_hash_table_product(mpp->mul, eh[0], ht)];

  pivs[idx]         = (src_t *)malloc((2*nt+2) * sizeof(src_t));

  pivs[idx][0]      = mpp->bi;
  pivs[idx][1]      = 2 * nt + 2; /* bi + length + pos/val */

  for (nelts_t i = 0; i < nt; ++i) {
    pivs[idx][2*i+2]  = ht->idx[find_in_hash_table_product(mpp->mul, eh[i], ht)];
    pivs[idx][2*i+3]  = cf[i];
  }
}

static inline smc_t *generate_sparse_compact_matrix_new(const gb_t *basis,
    const sel_t *sel, const nelts_t nr, const nelts_t ncl,
    const nelts_t ncr, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smc_t *mat  = initialize_sparse_compact_matrix(nr, ncl, ncr, basis->mod);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->nr; i=i+BLOCK) {
      #pragma omp task
      {
        const nelts_t bs  = BLOCK > mat->nr-i ? mat->nr-i : BLOCK;
        /* printf("mat->nr %u | i %u | bs %u\n", mat->nr, i, bs); */
        /* printf("constructs row %u\n",i); */
        poly_to_sparse_compact_matrix_row_new(sel->mpp+i, mat->ncl + mat->ncr,
            bs, basis, mat->row+i);
      }
    }
    }
  }
  return mat;
}

static inline smc_t *generate_sparse_compact_matrix_offset_block(const gb_t *basis,
    const sel_t *sel, const nelts_t shift, const nelts_t nr, const nelts_t ncl,
    const nelts_t ncr)
{
  smc_t *mat  = initialize_sparse_compact_matrix(nr, ncl, ncr, basis->mod);

  poly_to_sparse_compact_matrix_row_offset_test(sel->mpp+shift, nr, basis, mat->row);
  return mat;
}

static inline smc_t *generate_sparse_compact_matrix_offset_test(const gb_t *basis,
    const sel_t *sel, const nelts_t nr, const nelts_t ncl,
    const nelts_t ncr, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smc_t *mat  = initialize_sparse_compact_matrix(nr, ncl, ncr, basis->mod);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->nr; i=i+BLOCK) {
      #pragma omp task
      {
        const nelts_t bs  = BLOCK > mat->nr-i ? mat->nr-i : BLOCK;
        /* printf("mat->nr %u | i %u | bs %u\n", mat->nr, i, bs); */
        /* printf("constructs row %u\n",i); */
        poly_to_sparse_compact_matrix_row_offset_test(
            sel->mpp+i, bs, basis, mat->row+i);
      }
    }
    }
  }
  return mat;
}

static inline smc_t *generate_sparse_compact_matrix_test(const gb_t *basis,
    const sel_t *sel, const nelts_t nr, const nelts_t ncl,
    const nelts_t ncr, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smc_t *mat  = initialize_sparse_compact_matrix(nr, ncl, ncr, basis->mod);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->nr; i=i+BLOCK) {
      #pragma omp task
      {
        const nelts_t bs  = BLOCK > mat->nr-i ? mat->nr-i : BLOCK;
        /* printf("mat->nr %u | i %u | bs %u\n", mat->nr, i, bs); */
        /* printf("constructs row %u\n",i); */
        poly_to_sparse_compact_matrix_row_test(sel->mpp+i,
            bs, basis, mat->row+i);
      }
    }
    }
  }
  return mat;
}

static inline smc_t *generate_sparse_compact_matrix_offset(const gb_t *basis,
    const sel_t *sel, const nelts_t nr, const nelts_t ncl,
    const nelts_t ncr, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smc_t *mat  = initialize_sparse_compact_matrix(nr, ncl, ncr, basis->mod);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->nr; i=i+BLOCK) {
      #pragma omp task
      {
        const nelts_t bs  = BLOCK > mat->nr-i ? mat->nr-i : BLOCK;
        /* printf("mat->nr %u | i %u | bs %u\n", mat->nr, i, bs); */
        /* printf("constructs row %u\n",i); */
        poly_to_sparse_compact_matrix_row_offset(sel->mpp+i, mat->ncl + mat->ncr,
            bs, basis, mat->row+i);
      }
    }
    }
  }
  return mat;
}

static inline smc_t *generate_sparse_compact_matrix_offset_true_columns(const gb_t *basis,
    const sel_t *sel, const nelts_t dim, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smc_t *mat  = initialize_sparse_compact_matrix(dim, dim, dim, basis->mod);
  for (nelts_t i=0; i<dim; ++i)
    mat->row[i] = NULL;
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<sel->load; ++i) {
      #pragma omp task
      {
        poly_to_sparse_compact_matrix_row_offset_true_columns(sel->mpp+i, basis, mat->row);
      }
    }
    }
  }
  return mat;
}

static inline smc_t *generate_sparse_compact_matrix_true_columns(const gb_t *basis,
    const sel_t *sel, const nelts_t dim, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smc_t *mat  = initialize_sparse_compact_matrix(sel->load, sel->load, dim, basis->mod);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->nr; ++i) {
      #pragma omp task
      {
        poly_to_sparse_compact_matrix_row_true_columns(sel->mpp+i, i, basis, mat->row);
      }
    }
    }
  }
  return mat;
}

static inline smc_t *generate_sparse_compact_matrix(const gb_t *basis,
    const sel_t *sel, const nelts_t nr, const nelts_t ncl,
    const nelts_t ncr, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smc_t *mat  = initialize_sparse_compact_matrix(nr, ncl, ncr, basis->mod);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->nr; i=i+BLOCK) {
      #pragma omp task
      {
        const nelts_t bs  = BLOCK > mat->nr-i ? mat->nr-i : BLOCK;
        /* printf("mat->nr %u | i %u | bs %u\n", mat->nr, i, bs); */
        /* printf("constructs row %u\n",i); */
        poly_to_sparse_compact_matrix_row(sel->mpp+i, mat->ncl, mat->ncl + mat->ncr,
            bs, basis, mat->row+i);
      }
    }
    }
  }
  return mat;
}

static inline smat_t *generate_sparse_matrix(const gb_t *basis,
    const sel_t *sel, const nelts_t nr, const nelts_t ncl,
    const nelts_t ncr, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;

  smat_t *mat  = initialize_sparse_matrix(nr, ncl, ncr, basis->mod);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->nr; ++i) {
      #pragma omp task
      {
        mat->row[i] = poly_to_sparse_matrix_row(sel->mpp+i, mat->ncl + mat->ncr, basis);
      }
    }
    }
  }
  return mat;
}


static inline mat_t *generate_gbla_matrix_test(const gb_t *basis,
    const gb_t *sf, const spd_t *spd, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;
  mat_t *mat  = initialize_gbla_matrix(spd, basis);
  const nelts_t fbr = spd->col->nlm == 0 ? 0 : (spd->col->nlm-1)/mat->bs + 1;
  const bi_t bs = mat->bs;
  const nelts_t ncb = mat->cbl + mat->cbr;
  /* allocate all possible memory in matrix for this block row */
  for (nelts_t j=0; j<mat->rbu; ++j) {
    for (nelts_t k=0; k<fbr; ++k) {
      mat->A->blocks[j][k].val = (re_t **)malloc(bs * sizeof(re_t *));
      mat->A->blocks[j][k].pos = (bi_t **)malloc(bs * sizeof(bi_t *));
      mat->A->blocks[j][k].sz  = (bi_t *)malloc(bs * sizeof(bi_t));
      for (nelts_t l=0; l<bs; ++l) {
        mat->A->blocks[j][k].val[l]  = (re_t *)malloc(bs * sizeof(re_t));
        mat->A->blocks[j][k].pos[l]  = (bi_t *)malloc(bs * sizeof(bi_t));
        mat->A->blocks[j][k].sz[l]   = 0;
      }
    }
    for (nelts_t k=0; k<(ncb-fbr); ++k) {
      mat->B->blocks[j][k].val = (re_t *)calloc(bs * bs, sizeof(re_t));
    }
  }
  for (nelts_t j=0; j<mat->rbl; ++j) {
    for (nelts_t k=0; k<fbr; ++k) {
      mat->C->blocks[j][k].val = (re_t **)malloc(bs * sizeof(re_t *));
      mat->C->blocks[j][k].pos = (bi_t **)malloc(bs * sizeof(bi_t *));
      mat->C->blocks[j][k].sz  = (bi_t *)malloc(bs * sizeof(bi_t));
      for (nelts_t l=0; l<bs; ++l) {
        mat->C->blocks[j][k].val[l]  = (re_t *)malloc(bs * sizeof(re_t));
        mat->C->blocks[j][k].pos[l]  = (bi_t *)malloc(bs * sizeof(bi_t));
        mat->C->blocks[j][k].sz[l]   = 0;
      }
    }
    for (nelts_t k=0; k<(ncb-fbr); ++k) {
      mat->D->blocks[j][k].val = (re_t *)calloc(bs * bs, sizeof(re_t));
    }
  }
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->rbu; ++i) {
      #pragma omp task
      {
        generate_row_blocks_no_buffer(mat->A, mat->B, i, spd->selu->load, spd->col->nlm,
            mat->bs, basis, sf, spd->selu);
      }
    }
    /* fill the lower part CD */
    for (nelts_t i=0; i<mat->rbl; ++i) {
      #pragma omp task
      {
        generate_row_blocks_no_buffer(mat->C, mat->D, i, spd->sell->load, spd->col->nlm,
            mat->bs, basis, sf, spd->sell);
      }
    }
    #pragma omp taskwait
    }
  }
  /* free useless allocated memory in A and B */
  for (nelts_t j=0; j<mat->rbu; ++j) {
    nelts_t cz  = 0;
    for (nelts_t l=0; l<fbr; ++l) {
      cz  = 0;
      for (nelts_t k=0; k<bs; ++k) {
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
    cf_t zb[bs*bs];
    memset(zb, 0, bs*bs*sizeof(cf_t));
    for (nelts_t l=0; l<(ncb-fbr); ++l) {
      if (memcmp(mat->B->blocks[j][l].val, zb, bs*bs*sizeof(cf_t)) == 0) {
        free(mat->B->blocks[j][l].val);
        mat->B->blocks[j][l].val = NULL;
      }
    }
  }
  for (nelts_t j=0; j<mat->rbl; ++j) {
    nelts_t cz  = 0;
    for (nelts_t l=0; l<fbr; ++l) {
      cz  = 0;
      for (nelts_t k=0; k<bs; ++k) {
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
    cf_t zb[bs*bs];
    memset(zb, 0, bs*bs*sizeof(cf_t));
    for (nelts_t l=0; l<(ncb-fbr); ++l) {
      if (memcmp(mat->D->blocks[j][l].val, zb, bs*bs*sizeof(cf_t)) == 0) {
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
static inline mat_t *generate_gbla_matrix_keep_A(const gb_t *basis,
    const gb_t *sf, const spd_t *spd, const int nthreads)
{
  /* constructing gbla matrices is not threadsafe at the moment */
  const int t = nthreads;
  mat_t *mat  = initialize_gbla_matrix_keep_A(spd, basis);
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
    /* fill the upper part AB */
    for (nelts_t i=0; i<mat->rbu; ++i) {
      #pragma omp task
      {
        generate_row_blocks_keep_A(mat->AR, mat->B, i, spd->selu->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, sf, spd->selu);
      }
    }
    /* fill the lower part CD */
    for (nelts_t i=0; i<mat->rbl; ++i) {
      #pragma omp task
      {
        generate_row_blocks_keep_A(mat->CR, mat->D, i, spd->sell->load, spd->col->nlm,
            mat->bs, mat->cbl+mat->cbr, basis, sf, spd->sell);
      }
    }
    }
    #pragma omp taskwait
  }

  return mat;
}
#endif

static inline sr_t *normalize_row(sr_t *row, const cf_t mod)
{
  nelts_t i;

  cf_t inv  = row->val[0];
  inverse_val_new(&inv, mod);
  const cf_t cinv = inv;

  i = row->sz % 2;
  i = i & 1 ? 1 : 0;
  while (i<row->sz) {
    row->val[i]   = (cf_t)MODP((bf_t)(row->val[i]) * cinv, mod);
    row->val[i+1] = (cf_t)MODP((bf_t)(row->val[i+1]) * cinv, mod);
    i +=  2;
  }
  /* possibly not set in the unrolled loop above */
  row->val[0] = 1;

  return row;
}

static inline int cmp_sparse_rows_by_lead_column(const void *a,
    const void *b)
{
  const sr_t *ra  = *((sr_t **)a);
  const sr_t *rb  = *((sr_t **)b);

  return (int)ra->pos[0] - (int)rb->pos[0];
}

static inline void interreduce_pivots(sr_t **pivs, const nelts_t rk,
    const nelts_t ncl, const nelts_t ncr, const cf_t mod)
{

#if newred
  printf("PIVS BEGINNING LOWER REDUCTION PROCESS\n");
  for (int i=0; i<ncr; ++i)
    printf("%u %p\n", i, pivs[i]);
#endif
  const nelts_t nc  = ncl + ncr;
  nelts_t j;

  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf  = NULL;
  for (nelts_t k=0; k<rk; ++k) {
    const nelts_t shift = pivs[rk-k-1]->pos[0];
    bf = realloc(bf, (nc-shift) * sizeof(bf_t));
    memset(bf, 0, (nc-shift)*sizeof(bf_t));
    bf_t *bft = bf-shift;
    for (nelts_t j=0; j<pivs[rk-k-1]->sz; ++j) {
      bft[pivs[rk-k-1]->pos[j]]  = (bf_t)pivs[rk-k-1]->val[j];
    }
    free(pivs[rk-k-1]->pos);
    free(pivs[rk-k-1]->val);
    
    /* loop over all other pivots and reduce if possible */
    for (nelts_t l=0; l<k; ++l) {
      const nelts_t plc = pivs[rk-l-1]->pos[0];
      if (bf[plc-shift] != 0) {
        const sr_t *red = pivs[rk-l-1];
        const bf_t mul = (bf_t)(mod) - bft[plc]; /* it is already reduced modulo mod */
        /* loop unrolling by 2, test length in advance and probably add the useless
         * reduction on the pivot entry in order to get a even loop length */
        j = red->sz % 2;
        j = j & 1 ? 1 : 0;
        for (; j<red->sz; j+=2) {
#if newred
          printf("pos[%u] = %u\n",j, red->pos[j]-shift);
          printf("%lu + %lu * %u = ",bf[red->pos[j]-shift], mul, red->val[j]);
#endif
          bft[red->pos[j]]   +=  (bf_t)(mul) * red->val[j];
#if newred
          printf(" %lu\n", bf[red->pos[j]-shift]);
          printf("pos[%u] = %u\n",j+1, red->pos[j+1]-shift);
          printf("%lu - %lu * %u = ",bf[red->pos[j+1]-shift], mul, red->val[j+1]);
#endif
          bft[red->pos[j+1]] +=  (bf_t)(mul) * red->val[j+1];
#if newred
          printf(" %lu\n", bf[red->pos[j+1]-shift]);
#endif
        }
        bft[plc] = 0;
      }
    }
    /* write pivot back to sparse representation */
    pivs[rk-k-1]->sz   = nc-shift;
    pivs[rk-k-1]->val  = (cf_t *)malloc(pivs[rk-k-1]->sz * sizeof(cf_t));
    pivs[rk-k-1]->pos  = (nelts_t *)malloc(pivs[rk-k-1]->sz * sizeof(nelts_t));
    nelts_t ctr = 0;
    for (j=shift; j<nc; ++j) {
      bft[j] = bft[j] % mod;
      if (bft[j] != 0) {
        pivs[rk-k-1]->pos[ctr] = j;
        pivs[rk-k-1]->val[ctr] = (src_t)bft[j];
        ctr++;
      }
    } 
    /* fix memory */
    pivs[rk-k-1]->sz = ctr;
    pivs[rk-k-1]->val  = realloc(pivs[rk-k-1]->val, pivs[rk-k-1]->sz * sizeof(cf_t));
    pivs[rk-k-1]->pos  = realloc(pivs[rk-k-1]->pos, pivs[rk-k-1]->sz * sizeof(nelts_t));
  }
  free(bf);
}

static inline sr_t *reduce_lower_rows_by_pivots(sr_t *row, sr_t **pivs,
    const nelts_t start, const nelts_t ncl, const nelts_t ncr,
    const cf_t mod)
{

#if newred
  printf("PIVS BEGINNING LOWER REDUCTION PROCESS\n");
  for (int i=0; i<ncr; ++i)
    printf("%u %p\n", i, pivs[i]);
#endif
  const nelts_t nc  = ncl + ncr;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row->pos[0];
  nelts_t i;
  nelts_t j;
  nelts_t gate; /* used to check if we might have found a new pivot to be considered */
  const nelts_t shift = lc;

  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
#if newred
  printf("row->sz = %u\n", row->sz);
  printf("sz %u | shift %u | ncl %u | nc %u\n", row->sz, shift, ncl, nc);
#endif
  for (nelts_t i=0; i<row->sz; ++i) {
#if newred
    printf("i %u ||| ncr %u ||| ncl %u ||| shift %u\n",i, ncr, ncl, shift);
#endif
    bft[row->pos[i]]  = (bf_t)row->val[i];
  }

#if newred
  printf("buffer: ");
  for (nelts_t i=shift; i<nc; ++i)
    printf("%lu ",bf[i-shift]);
  printf("\n");
#endif
  /* free old lr data */
  free(row->pos);
  free(row->val);
  free(row);
  row = NULL;

  bf_t mul;

#if newred
  printf("nr %u\n",nr);
#endif
red_with_piv:
  i     = lc;
  lc    = UINT32_MAX;
  gate  = 0;
  /* while (i<nc) { */
  for (; i<nc; ++i) {
    if (bft[i])
      bft[i] = MODP(bft[i], mod);
    if (!bft[i])
      continue;
#if newred
    printf("RUN %u of %u\n", i, nc);
    printf("lc %u\n", lc);
    printf("i %u | gate %u\n", i, gate);
#endif
    const sr_t *red = pivs[i-start];
#if newred
    printf("i %u | start %u\n", i, start);
    printf("red %p\n", red);
#endif
    if (!red) {
      if (!gate)
        lc  = i;
      ++gate;
      continue;
    }
#if newred
    printf("we reduce with:\n");
    for (int ii=0; ii<red->sz; ++ii) {
      printf("%u | %u -- ", red->val[ii], red->pos[ii]);
    }
    printf("\n");
#endif
    mul = (bf_t)(mod) - bft[i]; /*it is already reduced modulo mod */
#if newred
    printf("get multiple out of %lu --> %lu\n", bf[i-shift],mul);
#endif
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    j = red->sz % 2;
    j = j & 1 ? 1 : 0;
    for (; j<red->sz; j+=2) {
#if newred
      printf("pos[%u] = %u\n",j, red->pos[j]-shift);
      printf("%lu + %lu * %u = ",bf[red->pos[j]-shift], mul, red->val[j]);
#endif
      bft[red->pos[j]]   +=  (bf_t)(mul) * red->val[j];
#if newred
      printf(" %lu\n", bf[red->pos[j]-shift]);
      printf("pos[%u] = %u\n",j+1, red->pos[j+1]-shift);
      printf("%lu - %lu * %u = ",bf[red->pos[j+1]-shift], mul, red->val[j+1]);
#endif
      bft[red->pos[j+1]] +=  (bf_t)(mul) * red->val[j+1];
#if newred
      printf(" %lu\n", bf[red->pos[j+1]-shift]);
#endif
    }
    /* due to above loop unrolling we cannot be sure if we have set the old
     * pivot entry to zero */
    bft[i] = 0;
#if newred
    printf("buffer:  ");
    for (int ii=0; ii<nc-shift; ++ii)
      printf("%lu ", bf[ii]);
    printf("\n");
#endif
    /* next pivot */
    /* lc  = j; */
#if newred
    printf("i %u | j %u | lc %u\n", i, j, lc);
#endif
    if (lc == nc) {
      free(bf);
      return row;
    }
    /* probably new pivots are added by other threads */
    if (gate && pivs[lc-start])
      goto red_with_piv;
  }
  /* row is completely reduced */
  if (gate==0) {
    free(bf);
#if newred
    printf("REDUCED TO NULL!\n");
#endif
    return row; /* i.e. row = NULL */
  }
  /* probably new pivots are added by other threads */
  if (gate && pivs[lc-start])
    goto red_with_piv;

#if newred
  printf("gate %u\n", gate);
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  row = (sr_t *)malloc(sizeof(sr_t));
  /* allocate maximal possibly needed memory */
  row->sz   = nc-shift;
  row->val  = (cf_t *)malloc(row->sz * sizeof(cf_t));
  row->pos  = (nelts_t *)malloc(row->sz * sizeof(nelts_t));
  nelts_t ctr = 0;
  for (j=shift; j<nc; ++j) {
    bft[j] = bft[j] % mod;
    if (bft[j] != 0) {
      row->pos[ctr] = j;
      row->val[ctr] = (src_t)bft[j];
      ctr++;
    }
  } 
  /* fix memory */
  row->sz = ctr;
#if newred
  printf("row->sz %u\n", row->sz);
#endif
  row->val  = realloc(row->val, row->sz * sizeof(cf_t));
  row->pos  = realloc(row->pos, row->sz * sizeof(nelts_t));

#if newred
  printf("2 resulting row: ");
  for (i=0; i<row->sz; ++i)
    printf("%u : %u | ", row->val[i], row->pos[i]);
  printf("\n");
#endif

  if (row->val[0] != 1)
    row = normalize_row(row, mod);

#if newred
  printf("2 resulting row normalized: ");
  for (i=0; i<row->sz; ++i)
    printf("%u : %u | ", row->val[i], row->pos[i]);
  printf("\n");
#endif

  free(bf);

  return row;
}

static inline void compute_new_pivots(sr_t *row, sr_t **pivs,
    const nelts_t shift, const nelts_t ncl, const nelts_t ncr, const cf_t mod)
{
  int done;

  do {
    row = reduce_lower_rows_by_pivots(row, pivs, shift, ncl, ncr, mod);
    if (!row) 
      return;
    done  = __sync_bool_compare_and_swap(&pivs[row->pos[0]-shift], NULL, row);
  } while (!done);
}

static inline void reduce_lower_rows(smat_t *mat, const nelts_t shift, const int nthreads)
{
#if newred
  printf("mat %p, mat->row %p, mat->row[0] %p\n", mat, mat->row, mat->row[0]);
  printf("mat->nr %u\n", mat->nr);
#endif
  /* sort rows by lead columns */
  qsort(mat->row, mat->nr, sizeof(sr_t **), cmp_sparse_rows_by_lead_column);
  /* initialize holder for pivots that we find during the reduction */
  /* sr_t **pivs = (sr_t **)malloc(mat->nr * sizeof(sr_t *)); */
  sr_t **pivs = (sr_t **)malloc((mat->ncr+mat->ncl-shift) * sizeof(sr_t *));
  for (nelts_t i=0; i<(mat->ncr+mat->ncl-shift); ++i) {
    pivs[i] = NULL;
  }
  nelts_t j = 0;
  for (nelts_t i=0; i<mat->nr; ++i) {
#if newred
    printf("i %u | nrl %u | pos[0] %u | ncl %u\n", i, mat->nr, mat->row[i]->pos[0], mat->ncl);
#endif
    if (!pivs[mat->row[i]->pos[0]-shift])
      pivs[mat->row[i]->pos[0]-shift] = mat->row[i];
    else
      mat->row[j++] = mat->row[i];
  }
  /* printf("j = %u\n", j); */
#pragma omp parallel for num_threads(nthreads)
  for (nelts_t i=0; i<j; ++i) {
    /* printf("row[%u] = %p for rank %u\n", i, mat->row[i], mat->rk); */
    compute_new_pivots(mat->row[i], pivs, shift, mat->ncl, mat->ncr, mat->mod);
  }
#if newred
  printf("FINAL PIVOTS:\n");
  for (nelts_t i=0; i<(mat->ncr+mat->ncl-shift); ++i) {
    printf("piv[%u] (%p)",i, pivs[i]);
    if (pivs[i] == NULL)
      printf("NULL\n");
    else {
      for (int jj=0; jj<pivs[i]->sz; ++jj) {
        printf("%u at %u | ",pivs[i]->val[jj],pivs[i]->pos[jj]);
      }
      printf("\n");
    }
  }
#endif
  nelts_t ctr=0;
  for (nelts_t i=0; i<(mat->ncr+mat->ncl-shift); ++i) {
    if (pivs[i] != NULL) {
      mat->row[ctr++] = pivs[i];
    }
  }
  mat->rk = ctr;
  for (nelts_t i=ctr; i<mat->nr; ++i)
    mat->row[i] = NULL;
  interreduce_pivots(mat->row, mat->rk, mat->ncl, mat->ncr, mat->mod);
  free(pivs);
}

static inline sr_t *reduce_lower_by_one_upper_row(sr_t *row, const sr_t *piv, const cf_t mod, const nelts_t nc)
{
  if (row->pos[0] > piv->pos[0]) {
    return row;
  }

  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row->pos[0];
  nelts_t i, j;
  const nelts_t shift = lc;
#if newred
  printf("ROW TO BE REDUCED: ");
  for (i=0; i<row->sz; ++i) {
    printf("%u | %u || ", row->val[i], row->pos[i]);
  }
  printf("\n");

  printf("lc %u | ncl %u\n", lc, ncl);
#endif
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc((nc-shift), sizeof(bf_t));
  bf_t *bft = bf-shift;
  i = 0;
  while (i<row->sz) {
    bft[row->pos[i]]  = row->val[i];
    ++i;
  }
  free(row->pos);
  free(row->val);
  free(row);
  row = NULL;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul;
#if newred
    printf("lc %u == %u ?\n", lc, pivs->row[lc]->pos[0]);
#endif
    const sr_t *red = piv;
#if newred
    printf("we reduce with:\n");
    for (int ii=0; ii<red->sz; ++ii) {
      printf("%u | %u -- ", red->val[ii], red->pos[ii]);
    }
    printf("\n");
#endif
    mul = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    j = red->sz % 2;
    j = j & 1 ? 1 : 0;
    for (; j<red->sz; j+=2) {
#if newred
      printf("pos[%u] = %u\n",j, red->pos[j]-shift);
      printf("%lu + %lu * %u = ",bf[red->pos[j]-shift], mul, red->val[j]);
#endif
      bft[red->pos[j]]   +=  (bf_t)(mul) * red->val[j];
#if newred
      printf(" %lu\n", bf[red->pos[j]-shift]);
      printf("pos[%u] = %u\n",j+1, red->pos[j+1]-shift);
      printf("%lu - %lu * %u = ",bf[red->pos[j+1]-shift], mul, red->val[j+1]);
#endif
      bft[red->pos[j+1]] +=  (bf_t)(mul) * red->val[j+1];
#if newred
      printf(" %lu\n", bf[red->pos[j+1]-shift]);
#endif
    }
#if newred
    printf("buffer-int:  ");
    for (i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* due to above loop unrolling we cannot be sure if we have set the old
     * pivot entry to zero */
    bft[lc] = 0;
    /* get new pivot element in row */
    for (j=lc+1; j<nc; ++j) {
      if (bft[j] != 0) {
        bft[j] = MODP(bft[j], mod);
      }
      if (bft[j] != 0) {
        break;
      }
    }
#if newred
    printf("buffer:  ");
    for (i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* next pivot */
    lc  = j;

  if (lc == nc) {
#if newred
    printf("NULL RETURNED!\n");
#endif
    free(bf);
    return row; /* i.e. row = NULL */
  }
  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  row = (sr_t *)malloc(sizeof(sr_t));
  /* allocate maximal possibly needed memory */
  row->sz   = nc-lc;
  row->val  = (cf_t *)malloc(row->sz * sizeof(cf_t));
  row->pos  = (nelts_t *)malloc(row->sz * sizeof(nelts_t));
  nelts_t ctr = 0;
  for (j=lc; j<nc; ++j) {
    bft[j] = bft[j] % mod;
    if (bft[j] != 0) {
      row->pos[ctr] = j;
      row->val[ctr] = (src_t)bft[j];
      ctr++;
    }
  }
  /* fix memory */
  row->sz = ctr;
#if newred
  printf("row->sz %u\n", row->sz);
#endif
  row->val  = realloc(row->val, row->sz * sizeof(cf_t));
  row->pos  = realloc(row->pos, row->sz * sizeof(nelts_t));

#if newred
  printf("1 resulting row: ");
  for (i=0; i<row->sz; ++i)
    printf("%u : %u | ", row->val[i], row->pos[i]);
  printf("\n");
#endif

  free(bf);
/*
  if (row->val[0] != 1)
    row = normalize_row(row, mod);
*/
#if newred
  printf("1 resulting row normalized: ");
  for (i=0; i<row->sz; ++i)
    printf("%u : %u | ", row->val[i], row->pos[i]);
  printf("\n");
#endif

  return row;
}
static inline sr_t *reduce_lower_by_upper_rows(sr_t *row, const smat_t *pivs)
{
  if (row->pos[0] >= pivs->ncl)
    return row;

  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  const nelts_t ncl = pivs->ncl;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row->pos[0];
  nelts_t i, j;
  const nelts_t shift = lc;
#if newred
  printf("ROW TO BE REDUCED: ");
  for (i=0; i<row->sz; ++i) {
    printf("%u | %u || ", row->val[i], row->pos[i]);
  }
  printf("\n");

  printf("lc %u | ncl %u\n", lc, pivs->ncl);
#endif
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc((nc-shift), sizeof(bf_t));
  bf_t *bft = bf-shift;
  i = 0;
  while (i<row->sz) {
    bft[row->pos[i]]  = row->val[i];
    ++i;
  }
  free(row->pos);
  free(row->val);
  free(row);
  row = NULL;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul;
  while (lc < ncl) {
#if newred
    printf("lc %u == %u ?\n", lc, pivs->row[lc]->pos[0]);
#endif
    const sr_t *red = pivs->row[lc];
#if newred
#endif

    mul = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    j = red->sz % 2;
    j = j & 1 ? 1 : 0;
    for (; j<red->sz; j+=2) {
#if newred
      printf("pos[%u] = %u\n",j, red->pos[j]-shift);
      printf("%lu + %lu * %u = ",bf[red->pos[j]-shift], mul, red->val[j]);
#endif
      bft[red->pos[j]]   +=  (bf_t)(mul) * red->val[j];
#if newred
      printf(" %lu\n", bf[red->pos[j]-shift]);
      printf("pos[%u] = %u\n",j+1, red->pos[j+1]-shift);
      printf("%lu - %lu * %u = ",bf[red->pos[j+1]-shift], mul, red->val[j+1]);
#endif
      bft[red->pos[j+1]] +=  (bf_t)(mul) * red->val[j+1];
#if newred
      printf(" %lu\n", bf[red->pos[j+1]-shift]);
#endif
    }
#if newred
    printf("buffer-int:  ");
    for (i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* due to above loop unrolling we cannot be sure if we have set the old
     * pivot entry to zero */
    bft[lc] = 0;
    /* get new pivot element in row */
    for (j=lc+1; j<nc; ++j) {
      if (bft[j] != 0) {
        bft[j] = MODP(bft[j], mod);
      }
      if (bft[j] != 0) {
        break;
      }
    }
#if newred
    printf("buffer:  ");
    for (i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* next pivot */
    lc  = j;
  }

  if (lc == nc) {
#if newred
    printf("NULL RETURNED!\n");
#endif
    free(bf);
    return row; /* i.e. row = NULL */
  }
  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  row = (sr_t *)malloc(sizeof(sr_t));
  /* allocate maximal possibly needed memory */
  row->sz   = nc-lc;
  row->val  = (cf_t *)malloc(row->sz * sizeof(cf_t));
  row->pos  = (nelts_t *)malloc(row->sz * sizeof(nelts_t));
  nelts_t ctr = 0;
  for (j=lc; j<nc; ++j) {
    bft[j] = bft[j] % mod;
    if (bft[j] != 0) {
      row->pos[ctr] = j;
      row->val[ctr] = (src_t)bft[j];
      ctr++;
    }
  } 
  /* fix memory */
  row->sz = ctr;
#if newred
  printf("row->sz %u\n", row->sz);
#endif
  row->val  = realloc(row->val, row->sz * sizeof(cf_t));
  row->pos  = realloc(row->pos, row->sz * sizeof(nelts_t));

#if newred
  printf("1 resulting row: ");
  for (i=0; i<row->sz; ++i)
    printf("%u : %u | ", row->val[i], row->pos[i]);
  printf("\n");
#endif

  free(bf);

  if (row->val[0] != 1)
    row = normalize_row(row, mod);

#if newred
  printf("1 resulting row normalized: ");
  for (i=0; i<row->sz; ++i)
    printf("%u : %u | ", row->val[i], row->pos[i]);
  printf("\n");
#endif

  return row;
}
static inline src_t *normalize(src_t *row, const cf_t mod)
{
  if (row[3] != 1) {
    nelts_t i;

    /* cf_t inv  = row[3];
     * inverse_val_new(&inv, mod);
     * const bf_t cinv = (bf_t)inv; */
    const cf_t cinv = (bf_t)inverse_coeff((int32_t)row[3], (int32_t)mod);

    i = (row[1]-2)/2;
    i = i & 1 ? 5 : 3;
    while (i < row[1]) {
      row[i]   = (src_t)(row[i] * cinv) % mod;
      row[i+2]   = (src_t)(row[i+2] * cinv) % mod;
      i +=  4;
    }
    /* possibly not set in the unrolled loop above */
    row[3] = 1;
  }

  return row;
}

static inline src_t *normalize_new_pivot(src_t *row, const cf_t mod)
{
  nelts_t i;

  /* cf_t inv  = row[3];
   * inverse_val_new(&inv, mod);
   * const bf_t cinv = (bf_t)inv; */
  const cf_t cinv = (bf_t)inverse_coeff((int32_t)row[3], (int32_t)mod);

  i = (row[1]-2)/2;
  i = i & 1 ? 5 : 3;
  while (i<row[1]) {
    row[i]   = (src_t)(row[i] * cinv) % mod;
    row[i+2]   = (src_t)(row[i+2] * cinv) % mod;
    i +=  4;
  }
  /* possibly not set in the unrolled loop above */
  row[3] = 1;

  return row;
}

static inline src_t *normalize_row_c(src_t *row, const cf_t mod)
{
  nelts_t i;

  cf_t inv  = row[2];
  inverse_val_new(&inv, mod);
  const cf_t cinv = inv;

  i = (row[0]-1)/2;
  i = i & 1 ? 4 : 2;
  while (i<row[0]) {
    row[i]   = (src_t)MODP((bf_t)(row[i]) * cinv, mod);
    row[i+2] = (src_t)MODP((bf_t)(row[i+2]) * cinv, mod);
    i +=  4;
  }
  /* possibly not set in the unrolled loop above */
  row[2] = 1;

  return row;
}

static inline int cmp_sparse_rows_by_lead_column_offset_c(const void *a,
    const void *b)
{
  const src_t *ra  = *((src_t **)a);
  const src_t *rb  = *((src_t **)b);

  return (int)ra[2] - (int)rb[2];
}

static inline int cmp_sparse_rows_by_lead_column_c(const void *a,
    const void *b)
{
  const src_t *ra  = *((src_t **)a);
  const src_t *rb  = *((src_t **)b);

  return (int)ra[1] - (int)rb[1];
}

static inline void interreduce_pivots_offset_c(src_t **pivs, const nelts_t rk,
    const nelts_t ncl, const nelts_t ncr, const cf_t mod)
{

  /* printf("rank %u\n", rk); */
  const nelts_t nc  = ncl + ncr;
  nelts_t j;

  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf  = NULL;
  for (nelts_t k=0; k<rk; ++k) {
    const nelts_t idx   = rk-k-1;
    const nelts_t shift = pivs[idx][2];
    /* printf("k %u -- shift %u\n", k, shift); */
    bf = realloc(bf, (nc-shift) * sizeof(bf_t));
    bf_t *bft = bf-shift;
    memset(bf, 0, (nc-shift)*sizeof(bf_t));
    for (nelts_t j=2; j<pivs[idx][0]; j=j+2) {
      /*
      printf("pivs[%u] = %p\n", idx, pivs[idx]);
      printf("shift %u | nc-shift %u\n", shift, nc-shift);
      printf("rk %u | k %u\n", rk, k);
      printf("nc %u | j %u\n", nc, j);
      printf("sz %u\n", pivs[idx][0]);
      printf("%u\n", pivs[idx][j]);
      printf("%u\n", pivs[idx][j+1]);
      */
      bft[pivs[idx][j]]  = (bf_t)pivs[idx][j+1];
    }
    free(pivs[idx]);
    
    /* loop over all other pivots and reduce if possible */
    for (nelts_t l=0; l<k; ++l) {
      const nelts_t plc = pivs[rk-l-1][2];
      /* printf("rk-l-1 %u | plc %u\n", rk-l-1,plc); */
      if (bft[plc] != 0) {
        const src_t *red = pivs[rk-l-1];
        const bf_t mul = (bf_t)(mod) - bft[plc]; /* it is already reduced modulo mod */
        /* loop unrolling by 2, test length in advance and probably add the useless
         * reduction on the pivot entry in order to get a even loop length */
#if 0
        j = j & 1 ? 3 : 1;
        for (; j<red[0]; j+=4) {
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
          bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
        }
        /* due to above loop unrolling we cannot be sure if we have set the old
         * pivot entry to zero */
#else
        /*
        const nelts_t mod8  = j % 8;
        for (j=1; j<2*mod8; j=j+2) {
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
        }
        for (; j<red[0]; j+=16) {
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
          bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
          bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
          bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
          bft[red[j+8]] +=  (bf_t)(mul) * red[j+9];
          bft[red[j+10]] +=  (bf_t)(mul) * red[j+11];
          bft[red[j+12]] +=  (bf_t)(mul) * red[j+13];
          bft[red[j+14]] +=  (bf_t)(mul) * red[j+15];
        }
        */
        for (j=2; j<2*red[1]; j=j+2)
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
        /* for (; j<red[0]; j+=8) { */
        for (; j<red[0]; j+=4) {
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
          bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
          /* bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
           * bft[red[j+6]] +=  (bf_t)(mul) * red[j+7]; */
        }
#endif
        bft[plc] = 0;
      }
    }
    /* write pivot back to sparse representation */
    pivs[idx]    = (src_t *)malloc((2*(nc-shift)+2)*sizeof(src_t));
    nelts_t ctr = 2;
    for (j=0; j<nc-shift; ++j) {
      bf[j] = bf[j] % mod;
      if (bf[j] != 0) {
        pivs[idx][ctr] = j+shift;
        pivs[idx][ctr+1] = (src_t)bf[j];
        ctr +=  2;
      }
    } 
    /* fix memory */
    pivs[idx][0]  = ctr;
    pivs[idx][1]  = ((ctr-2)/2) % 2 +1;
    /*
    printf("rk %u -- k %u \n", rk, k);
    printf("sz %u \n",pivs[idx][0]);
    printf("%p\n",pivs[idx]);
    */
    pivs[idx]    = realloc(pivs[idx], pivs[idx][0]*sizeof(src_t));
  }
  free(bf);
}


static inline void interreduce_pivots_c(src_t **pivs, const nelts_t rk,
    const nelts_t ncl, const nelts_t ncr, const cf_t mod)
{

#if newred
  printf("PIVS BEGINNING LOWER REDUCTION PROCESS\n");
  for (int i=0; i<ncr; ++i)
    printf("%u %p\n", i, pivs[i]);
#endif
  const nelts_t nc  = ncl + ncr;
  nelts_t j;

  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf  = NULL;
  for (nelts_t k=0; k<rk; ++k) {
    const nelts_t idx   = rk-k-1;
    const nelts_t shift = pivs[idx][1];
    bf = realloc(bf, (nc-shift) * sizeof(bf_t));
    bf_t *bft = bf-shift;
    memset(bf, 0, (nc-shift)*sizeof(bf_t));
    for (nelts_t j=1; j<pivs[idx][0]; j=j+2) {
      /*
      printf("pivs[%u] = %p\n", idx, pivs[idx]);
      printf("shift %u | nc-shift %u\n", shift, nc-shift);
      printf("rk %u | k %u\n", rk, k);
      printf("nc %u | j %u\n", nc, j);
      printf("sz %u\n", pivs[idx][0]);
      printf("%u\n", pivs[idx][j]);
      printf("%u\n", pivs[idx][j+1]);
      */
      bft[pivs[idx][j]]  = (bf_t)pivs[idx][j+1];
    }
    free(pivs[idx]);
    
    /* loop over all other pivots and reduce if possible */
    for (nelts_t l=0; l<k; ++l) {
      const nelts_t plc = pivs[rk-l-1][1];
      if (bft[plc] != 0) {
        const src_t *red = pivs[rk-l-1];
        const bf_t mul = (bf_t)(mod) - bft[plc]; /* it is already reduced modulo mod */
        /* loop unrolling by 2, test length in advance and probably add the useless
         * reduction on the pivot entry in order to get a even loop length */
        j = (red[0]-1)/2;
#if 0
        j = j & 1 ? 3 : 1;
        for (; j<red[0]; j+=4) {
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
          bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
        }
        /* due to above loop unrolling we cannot be sure if we have set the old
         * pivot entry to zero */
#else
        /*
        const nelts_t mod8  = j % 8;
        for (j=1; j<2*mod8; j=j+2) {
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
        }
        for (; j<red[0]; j+=16) {
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
          bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
          bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
          bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
          bft[red[j+8]] +=  (bf_t)(mul) * red[j+9];
          bft[red[j+10]] +=  (bf_t)(mul) * red[j+11];
          bft[red[j+12]] +=  (bf_t)(mul) * red[j+13];
          bft[red[j+14]] +=  (bf_t)(mul) * red[j+15];
        }
        */
        const nelts_t mod4  = j % 4;
        for (j=1; j<2*mod4; j=j+2)
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
        for (; j<red[0]; j+=8) {
          bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
          bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
          bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
          bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
        }
#endif
        bft[plc] = 0;
      }
    }
    /* write pivot back to sparse representation */
    pivs[idx]    = (src_t *)malloc((2*(nc-shift)+1)*sizeof(src_t));
    nelts_t ctr = 1;
    for (j=0; j<nc-shift; ++j) {
      bf[j] = bf[j] % mod;
      if (bf[j] != 0) {
        pivs[idx][ctr] = j+shift;
        pivs[idx][ctr+1] = (src_t)bf[j];
        ctr +=  2;
      }
    } 
    /* fix memory */
    pivs[idx][0] = ctr;
    /*
    printf("rk %u -- k %u \n", rk, k);
    printf("sz %u \n",pivs[idx][0]);
    printf("%p\n",pivs[idx]);
    */
    pivs[idx]    = realloc(pivs[idx], pivs[idx][0]*sizeof(src_t));
  }
  free(bf);
}

static inline src_t *reduce_lower_rows_by_pivots_c(src_t *row, src_t **pivs,
    const nelts_t ncl, const nelts_t ncr, const cf_t mod)
{

#if newred
  printf("PIVS BEGINNING LOWER REDUCTION PROCESS\n");
  for (int i=0; i<ncr; ++i)
    printf("%u %p\n", i, pivs[i]);
#endif
  const nelts_t nc  = ncl + ncr;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row[1];
  nelts_t i;
  nelts_t j;
  nelts_t gate; /* used to check if we might have found a new pivot to be considered */
  const nelts_t shift = lc;

  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
  for (nelts_t j=1; j<row[0]; j=j+2) {
    bft[row[j]]  = (bf_t)row[j+1];
  }
  /* free old lr data */
  free(row);
  row = NULL;

#if newred
  printf("buffer: ");
  for (nelts_t i=shift; i<nc; ++i)
    printf("%lu ",bf[i-shift]);
  printf("\n");
#endif

  bf_t mul;
  bf_t t0;

#if newred
  printf("nr %u\n",nr);
#endif
red_with_piv:
  i     = lc;
  lc    = UINT32_MAX;
  gate  = 0;
  /* while (i<nc) { */
  for (; i<nc; ++i) {
    t0  = bf[i-shift];
    if (t0)
      t0 = MODP(t0, mod);
    if (!t0)
      continue;
#if newred
    printf("RUN %u of %u\n", i, nc);
    printf("lc %u\n", lc);
    printf("i %u | gate %u\n", i, gate);
#endif
    const src_t *red = pivs[i-ncl];
#if newred
    printf("red %p\n", red);
#endif
    if (!red) {
      if (!gate)
        lc  = i;
      ++gate;
      continue;
    }
#if newred
    printf("we reduce with:\n");
    for (int ii=0; ii<red[0]; ++ii) {
      printf("%u | %u -- ", red[ii+2], red[ii+1]);
    }
    printf("\n");
#endif
    /* mul = (bf_t)(mod) - bf[i-shift]; it is already reduced modulo mod */
    mul = (bf_t)(mod) - t0; /* it is already reduced modulo mod */
#if newred
    printf("get multiple out of %lu --> %lu\n", bf[i-shift],mul);
#endif
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    if (red[0] == 0)
      j = nc-red[1];
    else 
      j = (red[0]-1)/2;
#if 0
    j = j & 1 ? 3 : 1;
    for (; j<red[0]; j+=4) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
    }
    /* due to above loop unrolling we cannot be sure if we have set the old
     * pivot entry to zero */
#else
    /*
    const nelts_t mod8  = j % 8;
    for (j=1; j<2*mod8; j=j+2) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
    }
    for (; j<red[0]; j+=16) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
      bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
      bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
      bft[red[j+8]] +=  (bf_t)(mul) * red[j+9];
      bft[red[j+10]] +=  (bf_t)(mul) * red[j+11];
      bft[red[j+12]] +=  (bf_t)(mul) * red[j+13];
      bft[red[j+14]] +=  (bf_t)(mul) * red[j+15];
    }
    */
    const nelts_t mod4  = j % 4;
    for (j=1; j<2*mod4; j=j+2)
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
    for (; j<red[0]; j+=8) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
      bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
      bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
    }
#endif
    bft[i] = 0;
#if newred
    printf("buffer:  ");
    for (int ii=0; ii<nc-shift; ++ii)
      printf("%lu ", bf[ii]);
    printf("\n");
#endif
    /* next pivot */
    /* lc  = j; */
#if newred
    printf("i %u | j %u | lc %u\n", i, j, lc);
#endif
    if (lc == nc) {
      free(bf);
      return row;
    }
    /* probably new pivots are added by other threads */
    if (gate && pivs[lc-ncl])
      goto red_with_piv;
  }
  /* row is completely reduced */
  if (gate==0) {
    free(bf);
#if newred
    printf("REDUCED TO NULL!\n");
#endif
    return row; /* i.e. row = NULL */
  }
  /* probably new pivots are added by other threads */
  if (gate && pivs[lc-ncl])
    goto red_with_piv;

#if newred
  printf("gate %u\n", gate);
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  /* allocate maximal possibly needed memory */
  row = (src_t *)malloc((2*gate+1)*sizeof(src_t));
  /* row = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t)); */
  nelts_t ctr = 1;
  for (j=0; j<nc-shift; ++j) {
    bf[j] = bf[j] % mod;
    if (bf[j] != 0) {
      row[ctr++] = j+shift;
      row[ctr++] = (src_t)bf[j];
    }
  } 
  /* fix memory */
  row[0] = ctr;
  /* row    = realloc(row, row[0]*sizeof(src_t)); */
#if newred
  printf("row->sz %u\n", row[0]);
#endif

  if (row[2] != 1)
    row = normalize_row_c(row, mod);

  free(bf);
  return row;
}

static inline void compute_new_pivots_c(src_t *row, src_t **pivs,
    const nelts_t ncl, const nelts_t ncr, const cf_t mod)
{
  int done;

  do {
    row = reduce_lower_rows_by_pivots_c(row, pivs, ncl, ncr, mod);
    if (!row) 
      return;
    done  = __sync_bool_compare_and_swap(&pivs[row[1]-ncl], NULL, row);
  } while (!done);
}

static inline void interreduce_upper_rows_offset_c(
    smc_t *mat, const nelts_t shift)
{
#if newred
  printf("mat %p, mat->row %p, mat->row[0] %p\n", mat, mat->row, mat->row[0]);
  printf("mat->nr %u\n", mat->nr);
#endif
  /* sort rows by lead columns */
  qsort(mat->row, mat->nr, sizeof(src_t **), cmp_sparse_rows_by_lead_column_offset_c);
  /* initialize holder for pivots that we find during the reduction */
  src_t **pivs = (src_t **)malloc(mat->rk * sizeof(src_t *));
  for (nelts_t i=0; i<mat->rk; ++i) {
    pivs[i] = NULL;
  }
  nelts_t j = 0;
  for (nelts_t i=0; i<mat->rk; ++i) {
    /* printf("row[%u][2] = %u | shift = %u\n", i, mat->row[i][2], shift); */
    if (!pivs[mat->row[i][2]-shift]) {
      pivs[mat->row[i][2]-shift] = mat->row[i];
    } else {
      mat->row[j++] = mat->row[i];
    }
  }
#if newred
  printf("FINAL PIVOTS:\n");
  for (nelts_t i=0; i<mat->rk; ++i) {
    printf("piv[%u] (%p)",i, pivs[i]);
    if (pivs[i] == NULL)
      printf("NULL\n");
    else {
      for (int jj=2; jj<pivs[i][0]; jj = jj+2) {
        printf("%u at %u | ",pivs[i][jj+1],pivs[i][jj]);
      }
      printf("\n");
    }
  }
#endif
  interreduce_pivots_offset_c(mat->row, mat->rk, mat->ncl, mat->ncr, mat->mod);
  free(pivs);
}

static inline void reduce_lower_rows_c(smc_t *mat, const nelts_t shift, const int nthreads)
{
#if newred
  printf("mat %p, mat->row %p, mat->row[0] %p\n", mat, mat->row, mat->row[0]);
  printf("mat->nr %u\n", mat->nr);
#endif
  /* sort rows by lead columns */
  qsort(mat->row, mat->nr, sizeof(src_t **), cmp_sparse_rows_by_lead_column_c);
  /* initialize holder for pivots that we find during the reduction */
  src_t **pivs = (src_t **)malloc(mat->ncr * sizeof(src_t *));
  for (nelts_t i=0; i<mat->ncr; ++i) {
    pivs[i] = NULL;
  }
  nelts_t j = 0;
  for (nelts_t i=0; i<mat->nr; ++i) {
#if newred
    printf("nrl %u | pos[0] %u | shift %u\n", mat->nr, mat->row[i][1], shift);
#endif
    if (!pivs[mat->row[i][1]-shift])
      pivs[mat->row[i][1]-shift] = mat->row[i];
    else
      mat->row[j++] = mat->row[i];
  }
#pragma omp parallel for num_threads(nthreads)
  for (nelts_t i=0; i<j; ++i)
    compute_new_pivots_c(mat->row[i], pivs, shift, mat->ncr, mat->mod);
#if newred
  printf("FINAL PIVOTS:\n");
  for (nelts_t i=0; i<mat->ncr; ++i) {
    printf("piv[%u] (%p)",i, pivs[i]);
    if (pivs[i] == NULL)
      printf("NULL\n");
    else {
      for (int jj=1; jj<pivs[i]->sz; ++jj) {
        printf("%u at %u | ",pivs[i][2*jj+1],pivs[i][2*jj]);
      }
      printf("\n");
    }
  }
#endif
  nelts_t ctr=0;
  for (nelts_t i=0; i<mat->ncr; ++i) {
    if (pivs[i] != NULL) {
      mat->row[ctr++] = pivs[i];
    }
  }
  mat->rk = ctr;
  for (nelts_t i=ctr; i<mat->nr; ++i)
    mat->row[i] = NULL;
  interreduce_pivots_c(mat->row, mat->rk, shift, mat->ncr, mat->mod);
  free(pivs);
}

static inline src_t *reduce_lower_by_upper_multiline_rows_c(src_t *row, const smc_t *pivs)
{
  if (row[1] >= pivs->ncl)
    return row;

  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  const nelts_t ncl = pivs->ncl;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row[1];
  nelts_t j;
  const nelts_t shift = lc;
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
  for (nelts_t j=1; j<row[0]; j=j+2) {
    bft[row[j]]  = (bf_t)row[j+1];
  }
  free(row);
  row = NULL;

  bf_t mul1, mul2;
  while (lc < ncl) {
    const src_t *red = pivs->row[lc/2];
    if (lc % 2) {
      mul1  = 0;
      mul2  = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
      /* loop unrolling by 2, test length in advance and probably add the useless
       * reduction on the pivot entry in order to get a even loop length */
      j = (red[0]-1)/3;
      const nelts_t mod4  = j % 4;
      /* printf("mod4 %u\n", mod4); */
      for (j=1; j<3*mod4; j=j+3) {
        /* printf("red[%u] %u == %u lc || red[0] %u ?\n", j, red[j], lc, red[0]); */
        bft[red[j]]   +=  (bf_t)(mul2) * red[j+2];
      }
      /* bft[lc] = 0; */
      for (; j<red[0]; j+=12) {
        bft[red[j]]   +=  (bf_t)(mul2) * red[j+2];
        bft[red[j+3]] +=  (bf_t)(mul2) * red[j+5];
        bft[red[j+6]] +=  (bf_t)(mul2) * red[j+8];
        bft[red[j+9]] +=  (bf_t)(mul2) * red[j+11];
      }
    } else {
      mul1  = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
      mul2  = (bf_t)(mod) - ((bft[lc+1] + mul1 * red[5])%mod); /* it is already reduced modulo mod */
      /* loop unrolling by 2, test length in advance and probably add the useless
       * reduction on the pivot entry in order to get a even loop length */
      j = (red[0]-1)/3;
      const nelts_t mod4  = j % 4;
      /* printf("mod4 %u\n", mod4); */
      for (j=1; j<3*mod4; j=j+3) {
        /* printf("red[%u] %u == %u lc || red[0] %u ?\n", j, red[j], lc, red[0]); */
        bft[red[j]]   +=  (bf_t)(mul1) * red[j+1] + (bf_t)(mul2) * red[j+2];
      }
      /* bft[lc] = 0; */
      for (; j<red[0]; j+=12) {
        bft[red[j]]   +=  (bf_t)(mul1) * red[j+1] + (bf_t)(mul2) * red[j+2];
        bft[red[j+3]] +=  (bf_t)(mul1) * red[j+4] + (bf_t)(mul2) * red[j+5];
        bft[red[j+6]] +=  (bf_t)(mul1) * red[j+7] + (bf_t)(mul2) * red[j+8];
        bft[red[j+9]] +=  (bf_t)(mul1) * red[j+10] + (bf_t)(mul2) * red[j+11];
      }
    }
    /* get new pivot element in row */
    for (j=lc+1-shift; j<nc-shift; ++j) {
      if (bf[j] != 0) {
        bf[j] = MODP(bf[j], mod);
      }
      if (bf[j] != 0) {
        break;
      }
    }
    /* next pivot */
    lc  = j+shift;
    if (lc >= nc) {
      free(bf);
      return row; /* i.e. row = NULL */
    }
  }

  /* move data from buffer back to sparse row
   * normalize row if needed
   * allocate maximal possibly needed memory */
  row = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t));
  nelts_t ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf[j] = bf[j] % mod;
    if (bf[j] != 0) {
      row[ctr++] = j+shift;
      row[ctr++] = (src_t)bf[j];
    }
  } 
  /* fix memory */
  row[0] = ctr;
  row    = realloc(row, row[0]*sizeof(src_t));

  if (row[2] != 1)
    row = normalize_row_c(row, mod);

  free(bf);
  return row;
}
static inline src_t *reduce_lower_by_upper_rows_new_c(src_t *row, const smc_t *pivs)
{
  if (row[1] >= pivs->ncl)
    return row;

  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  const nelts_t ncl = pivs->ncl;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row[1];
  nelts_t j;
  const nelts_t shift = lc;
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
  for (nelts_t j=1; j<row[0]; j=j+2) {
    bft[row[j]]  = (bf_t)row[j+1];
  }
  free(row);
  row = NULL;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul;
  while (lc < ncl) {
    const src_t *red = pivs->row[lc];
#if newred
    printf("lc %u\n", lc);
    printf("red %p --> red[0] %u\n", red, red[0]);
    printf("we reduce with:\n");
    for (int ii=1; ii<red[0]; ++ii) {
      printf("%u | %u -- ", red[2*ii+1], red[2*ii]);
    }
    printf("\n");
#endif
    mul = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    j = 1;
    nelts_t len;
    while (j<red[0]) {
      bft = bf - shift + red[j++];
      len = red[j++];
      /* printf("len %u | start %u\n", len, start); */
      /* bft[0]  +=  (bf_t)(mul) * red[j]; */
      for (nelts_t k=0; k<len; k=k+1) {
      /* for (nelts_t k=start; k<len; k=k+1) { */
        bft[k]  +=  (bf_t)(mul) * red[j+k];
        /* bft[pos+k+1]  +=  (bf_t)(mul) * red[j+k+1]; */
      }
      j +=  len;
    }
    bft = bf - shift;
    bft[lc]  = 0;
    /* get new pivot element in row */
    for (j=lc+1-shift; j<nc-shift; ++j) {
      if (bf[j] != 0) {
        bf[j] = MODP(bf[j], mod);
      }
      if (bf[j] != 0) {
        break;
      }
    }
#if newred
    printf("buffer:  ");
    for (i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* next pivot */
    lc  = j+shift;
    if (lc >= nc) {
#if newred
      printf("NULL RETURNED!\n");
#endif
      free(bf);
      return row; /* i.e. row = NULL */
    }
  }

  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  /* allocate maximal possibly needed memory */
  row = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t));
  nelts_t ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf[j] = bf[j] % mod;
    if (bf[j] != 0) {
      row[ctr++] = j+shift;
      row[ctr++] = (src_t)bf[j];
    }
  } 
  /* fix memory */
  row[0] = ctr;
  row    = realloc(row, row[0]*sizeof(src_t));
#if newred
  printf("row->sz %u\n", row[0]);
#endif

  if (row[2] != 1)
    row = normalize_row_c(row, mod);

  free(bf);
  return row;
}

static inline src_t *reduce_lower_by_upper_rows_pos_val_2_c(src_t *row, const smc_t *pivs)
{
  if (row[1] >= pivs->ncl)
    return row;

  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  const nelts_t ncl = pivs->ncl;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row[1];
  nelts_t j, k, len;
  const nelts_t shift = lc;
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
  for (nelts_t j=1; j<row[0]; j=j+2) {
    bft[row[j]]  = (bf_t)row[j+1];
  }
  free(row);
  row = NULL;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul;
  while (lc < ncl) {
    const src_t *red = pivs->row[lc];
#if newred
    printf("lc %u\n", lc);
    printf("red %p --> red[0] %u\n", red, red[0]);
    printf("we reduce with:\n");
    for (int ii=1; ii<red[0]; ++ii) {
      printf("%u | %u -- ", red[2*ii+1], red[2*ii]);
    }
    printf("\n");
#endif
    mul = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    j = len = (red[0]-1)/2+1;
    const nelts_t mod8  = (j-1) % 8;
    /* printf("mod4 %u\n", mod4); */
    k = len;
    for (j=1; j<mod8+1; j++) {
      /* printf("red[%u] %u == %u lc || red[0] %u ?\n", j, red[j], lc, red[0]); */
      bft[red[j]]   +=  (bf_t)(mul) * red[k++];
    }
    /* bft[lc] = 0; */
    for (; j<len; j+=8) {
      /* printf("begin k %u j %u | len %u | red[0] %u\n", k, j, len, red[0]); */
      bft[red[j]]   +=  (bf_t)(mul) * red[k++];
      bft[red[j+1]] +=  (bf_t)(mul) * red[k++];
      bft[red[j+2]] +=  (bf_t)(mul) * red[k++];
      bft[red[j+3]] +=  (bf_t)(mul) * red[k++];
      bft[red[j+4]] +=  (bf_t)(mul) * red[k++];
      bft[red[j+5]] +=  (bf_t)(mul) * red[k++];
      bft[red[j+6]] +=  (bf_t)(mul) * red[k++];
      bft[red[j+7]] +=  (bf_t)(mul) * red[k++];
    }
    /*
    const nelts_t mod4  = j % 4;
    for (j=1; j<2*mod4; j=j+2) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
    }
    for (; j<red[0]; j+=8) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
      bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
      bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
    }
    */
    /* get new pivot element in row */
    for (j=lc+1-shift; j<nc-shift; ++j) {
      if (bf[j] != 0) {
        bf[j] = MODP(bf[j], mod);
      }
      if (bf[j] != 0) {
        break;
      }
    }
#if newred
    printf("buffer:  ");
    for (i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* next pivot */
    lc  = j+shift;
    if (lc >= nc) {
#if newred
      printf("NULL RETURNED!\n");
#endif
      free(bf);
      return row; /* i.e. row = NULL */
    }
  }

  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  /* allocate maximal possibly needed memory */
  row = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t));
  nelts_t ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf[j] = bf[j] % mod;
    if (bf[j] != 0) {
      row[ctr++] = j+shift;
      row[ctr++] = (src_t)bf[j];
    }
  } 
  /* fix memory */
  row[0] = ctr;
  row    = realloc(row, row[0]*sizeof(src_t));
#if newred
  printf("row->sz %u\n", row[0]);
#endif

  if (row[2] != 1)
    row = normalize_row_c(row, mod);

  free(bf);
  return row;
}
static inline src_t *reduce_lower_by_upper_rows_pos_val_c(src_t *row, const smc_t *pivs)
{
  if (row[1] >= pivs->ncl)
    return row;

  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  const nelts_t ncl = pivs->ncl;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row[1];
  nelts_t j;
  const nelts_t shift = lc;
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bf2 = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
  bf_t *bft2;
  for (nelts_t j=1; j<row[0]; j=j+2) {
    bft[row[j]]  = (bf_t)row[j+1];
  }
  free(row);
  row = NULL;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul;
  while (lc < ncl) {
    const src_t *red = pivs->row[lc];
#if newred
    printf("lc %u\n", lc);
    printf("red %p --> red[0] %u\n", red, red[0]);
    printf("we reduce with:\n");
    for (int ii=1; ii<red[0]; ++ii) {
      printf("%u | %u -- ", red[2*ii+1], red[2*ii]);
    }
    printf("\n");
#endif
    mul = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    j = (red[0]-1)/2+1;
    /* read values */
    bft2 = bf2 - 1;
    for (nelts_t l=1; l<j; ++l) {
      bft2[l]  = bft[red[l]];
    }
    /* compute new values */
    bft2 = bf2 - j;
    for (nelts_t l=j; l<red[0]; ++l) {
      bft2[l]  +=  (bf_t)(mul) * red[l];
    }
    /* write values */
    bft2 = bf2 - 1;
    for (nelts_t l=1; l<j; ++l) {
      bft[red[l]] = bft2[l];
    }
    /* get new pivot element in row */
    for (j=lc+1-shift; j<nc-shift; ++j) {
      if (bf[j] != 0) {
        bf[j] = MODP(bf[j], mod);
      }
      if (bf[j] != 0) {
        break;
      }
    }
#if newred
    printf("buffer:  ");
    for (i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* next pivot */
    lc  = j+shift;
    if (lc >= nc) {
#if newred
      printf("NULL RETURNED!\n");
#endif
      free(bf);
      free(bf2);
      return row; /* i.e. row = NULL */
    }
  }

  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  /* allocate maximal possibly needed memory */
  row = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t));
  nelts_t ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf[j] = bf[j] % mod;
    if (bf[j] != 0) {
      row[ctr++] = j+shift;
      row[ctr++] = (src_t)bf[j];
    }
  } 
  /* fix memory */
  row[0] = ctr;
  row    = realloc(row, row[0]*sizeof(src_t));
#if newred
  printf("row->sz %u\n", row[0]);
#endif

  if (row[2] != 1)
    row = normalize_row_c(row, mod);

  free(bf);
  free(bf2);
  return row;
}
static inline void reduce_lower_by_upper_rows_double_c(src_t **row1p, src_t **row2p, const smc_t *pivs)
{
  src_t *row1 = *row1p;
  src_t *row2 = *row2p;
  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  const nelts_t ncl = pivs->ncl;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row1[1]<row2[1] ? row1[1] : row2[1];
  nelts_t j, j1, j2;
  const nelts_t shift = lc;
  /* printf("shift %u | lc %u | row1[1] %u | row2[1] %u\n", shift, lc, row1[1], row2[1]);   */
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf1 = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft1 = bf1-shift;
  for (nelts_t j=1; j<row1[0]; j=j+2) {
    bft1[row1[j]]  = (bf_t)row1[j+1];
  }
  free(row1);
  row1 = NULL;
  bf_t *bf2 = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft2 = bf2-shift;
  for (nelts_t j=1; j<row2[0]; j=j+2) {
    bft2[row2[j]]  = (bf_t)row2[j+1];
  }
  free(row2);
  row2 = NULL;
  *row1p  = row1;
  *row2p  = row2;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul1, mul2;
  while (lc < ncl) {
    const src_t *red = pivs->row[lc];
#if newred
    printf("lc %u\n", lc);
    printf("red %p --> red[0] %u\n", red, red[0]);
    printf("we reduce with:\n");
    for (int ii=1; ii<red[0]; ++ii) {
      printf("%u | %u -- ", red[2*ii+1], red[2*ii]);
    }
    printf("\n");
#endif
    mul1 = (bf_t)(mod) - bft1[lc]; /* it is already reduced modulo mod */
    mul2 = (bf_t)(mod) - bft2[lc]; /* it is already reduced modulo mod */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    j = (red[0]-1)/2;
#if 0
    j = j & 1 ? 3 : 1;
    for (; j<red[0]; j+=4) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
    }
#else
    const nelts_t mod8  = j % 8;
    /* printf("mod4 %u\n", mod4); */
    for (j=1; j<2*mod8; j=j+2) {
      /* printf("red[%u] %u == %u lc || red[0] %u ?\n", j, red[j], lc, red[0]); */
      bft1[red[j]]   +=  (bf_t)(mul1) * red[j+1];
      bft2[red[j]]   +=  (bf_t)(mul2) * red[j+1];
    }
    /* bft[lc] = 0; */
    for (; j<red[0]; j+=16) {
      bft1[red[j]]   +=  (bf_t)(mul1) * red[j+1];
      bft1[red[j+2]] +=  (bf_t)(mul1) * red[j+3];
      bft1[red[j+4]] +=  (bf_t)(mul1) * red[j+5];
      bft1[red[j+6]] +=  (bf_t)(mul1) * red[j+7];
      bft1[red[j+8]] +=  (bf_t)(mul1) * red[j+9];
      bft1[red[j+10]] +=  (bf_t)(mul1) * red[j+11];
      bft1[red[j+12]] +=  (bf_t)(mul1) * red[j+13];
      bft1[red[j+14]] +=  (bf_t)(mul1) * red[j+15];
      bft2[red[j]]   +=  (bf_t)(mul2) * red[j+1];
      bft2[red[j+2]] +=  (bf_t)(mul2) * red[j+3];
      bft2[red[j+4]] +=  (bf_t)(mul2) * red[j+5];
      bft2[red[j+6]] +=  (bf_t)(mul2) * red[j+7];
      bft2[red[j+8]] +=  (bf_t)(mul2) * red[j+9];
      bft2[red[j+10]] +=  (bf_t)(mul2) * red[j+11];
      bft2[red[j+12]] +=  (bf_t)(mul2) * red[j+13];
      bft2[red[j+14]] +=  (bf_t)(mul2) * red[j+15];
    }
    /*
    const nelts_t mod4  = j % 4;
    for (j=1; j<2*mod4; j=j+2) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
    }
    for (; j<red[0]; j+=8) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
      bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
      bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
    }
    */
#endif
    /* get new pivot element in row */
    for (j1=lc+1-shift; j1<nc-shift; ++j1) {
      if (bf1[j1] != 0) {
        bf1[j1] = MODP(bf1[j1], mod);
      }
      if (bf1[j1] != 0) {
        break;
      }
    }
    for (j2=lc+1-shift; j2<nc-shift; ++j2) {
      if (bf2[j2] != 0) {
        bf2[j2] = MODP(bf2[j2], mod);
      }
      if (bf2[j2] != 0) {
        break;
      }
    }
#if newred
    printf("buffer:  ");
    for (i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* next pivot */
    lc  =   j1<j2 ? j1 : j2;
    lc  +=  shift;
    if (lc >= nc) {
#if newred
      printf("NULL RETURNED!\n");
#endif
      free(bf1);
      free(bf2);
      return; /* i.e. row = NULL */
    }
  }

  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  /* allocate maximal possibly needed memory */
  row1 = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t));
  row2 = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t));
  nelts_t ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf1[j] = bf1[j] % mod;
    if (bf1[j] != 0) {
      row1[ctr++] = j+shift;
      row1[ctr++] = (src_t)bf1[j];
    }
  } 
  /* fix memory */
  row1[0] = ctr;
  if (ctr>1) {
    row1    = realloc(row1, row1[0]*sizeof(src_t));
    if (row1[2] != 1)
      row1 = normalize_row_c(row1, mod);
  } else {
    free(row1);
    row1  = NULL;
  }
#if newred
  printf("row->sz %u\n", row1[0]);
#endif
  ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf2[j] = bf2[j] % mod;
    if (bf2[j] != 0) {
      row2[ctr++] = j+shift;
      row2[ctr++] = (src_t)bf2[j];
    }
  } 
  /* fix memory */
  row2[0] = ctr;
  if (ctr>1) {
    row2    = realloc(row2, row2[0]*sizeof(src_t));
    if (row2[2] != 1)
      row2 = normalize_row_c(row2, mod);
  } else {
    free(row2);
    row2  = NULL;
  }
#if newred
  printf("row->sz %u\n", row2[0]);
#endif


  free(bf1);
  free(bf2);

  *row1p  = row1;
  *row2p  = row2;
  return;
}
static inline smc_t *reduce_upper_rows_c(smc_t *pivs)
{

  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  const nelts_t ncl = pivs->ncl;
  
  nelts_t lc, j;
  bf_t mul;
  for (int rn=(int)(pivs->nr-1); rn>-1; --rn) {
    src_t *row  = pivs->row[rn];
    if (row[0] == 1 || row[4] >= ncl)
      continue;
    /* printf("ROW IN\n");
     * for (int ii=2; ii<row[0]; ii=ii+2)
     *   printf("%u at pos %u\n", row[ii+1], row[ii]); */
    /* go to second term, lead term will be kept */
    lc  = row[2];
    const nelts_t shift = lc;
    /* store row in dense format using bf_t data type for delayed modulus
    * computations */
    bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
    bf_t *bft = bf-shift;
    for (nelts_t j=2; j<row[0]; j=j+2) {
      bft[row[j]]  = (bf_t)row[j+1];
    }
    lc  = row[4];
    free(row);
    row = NULL;
    while (lc < ncl) {
      const src_t *red = pivs->row[lc];
      /* printf("lc %u\n", lc);
      * printf("red %p --> red[0] %u -- red[1] %u\n", red, red[0], red[1]); */
#if newred
      printf("lc %u\n", lc);
      printf("red %p --> red[0] %u\n", red, red[0]);
      printf("we reduce with:\n");
      for (int ii=1; ii<red[0]; ++ii) {
        printf("%u | %u -- ", red[2*ii+1], red[2*ii]);
      }
      printf("\n");
#endif
      mul = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
      /* if (mul == 0)
      *   printf("mul is zero!\n"); */
      /* loop unrolling by 2, test length in advance and probably add the useless
      * reduction on the pivot entry in order to get a even loop length */
      /* printf("mod4 %u\n", mod4); */
      /* printf("red[0] %u | red[1] %u\n", red[0], red[1]); */
      for (j=2; j<2*red[1]; j=j+2) {
        /* printf("red[%u] %u == %u lc || red[0] %u ?\n", j, red[j], lc, red[0]); */
        bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      }
      /* bft[lc] = 0; */
      for (; j<red[0]; j+=4) {
        bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
        bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
        /* bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
        * bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
        * bft[red[j+8]] +=  (bf_t)(mul) * red[j+9];
        * bft[red[j+10]] +=  (bf_t)(mul) * red[j+11];
        * bft[red[j+12]] +=  (bf_t)(mul) * red[j+13];
        * bft[red[j+14]] +=  (bf_t)(mul) * red[j+15]; */
      }
      /* get new pivot element in row */
      for (j=lc+1-shift; j<nc-shift; ++j) {
        if (bf[j] != 0) {
          bf[j] = MODP(bf[j], mod);
        }
        if (bf[j] != 0) {
          break;
        }
      }
#if newred
      printf("buffer:  ");
      for (nelts_t i=0; i<nc-shift; ++i)
        printf("%lu ", bf[i]);
      printf("\n");
#endif
      /* next pivot */
      lc  = j+shift;
    }
    /* move data from buffer back to sparse row
    * normalize row if needed */
#if newred
    printf("nc - lc = %u - %u\n", nc, lc);
#endif
    /* allocate maximal possibly needed memory */
    row = (src_t *)malloc((2*(nc-shift))*sizeof(src_t));
    nelts_t ctr = 2;
    for (j=0; j<nc-shift; ++j) {
      bf[j] = bf[j] % mod;
      if (bf[j] != 0) {
        row[ctr++] = j+shift;
        row[ctr++] = (src_t)bf[j];
      }
    } 
    /* fix memory */
    row[0]  = ctr;
    row[1]  = ((ctr-2)/2) % 2 +1;
    row    = realloc(row, row[0]*sizeof(src_t));
     /* printf("row %p\n", row);
     * for (int ii=2; ii<row[0]; ii=ii+2)
     *   printf("%u at pos %u\n", row[ii+1], row[ii]); */
#if newred
    printf("row->sz %u\n", row[0]);
#endif
    free(bf);
    pivs->row[rn] = row;
  }
  return pivs;
}

static inline src_t *reduce_dense_row_by_known_pivots(bf_t *dr,
    src_t **pivs, const nelts_t nc, const mod_t mod)
{
  size_t i, j;
  nelts_t k = 0;
  for (i = 0; i < nc; ++i) {
    if (dr[i] != 0)
      dr[i]  = dr[i] % mod;
    if (dr[i] == 0)
      continue;
    if (pivs[i] == NULL) {
      ++k;
      continue;
    }
    
    /* reduce dense row with found pivot row */
    const bf_t mul = (bf_t)(mod) - dr[i]; /* it is already reduced modulo mod */
    
    j = (pivs[i][1] - 2) / 2;
    j = (j & 1) ? 4 : 2;
    for (; j < pivs[i][1]; j = j+4) {
      dr[pivs[i][j]]    +=  (bf_t)mul * pivs[i][j+1];
      dr[pivs[i][j+2]]  +=  (bf_t)mul * pivs[i][j+3];
    }
    dr[i]  = 0;

    /* get new pivot element in dense row */
    /* for (j = i+1; j < nc; ++j) {
     *   if (dr[j] != 0)
     *     dr[j] = dr[j] % mod;
     *   if (dr[j] != 0)
     *     break;
     * } */
    /* zero reduction of dense row */
    /* if (j == nc)
     *   return NULL; */
  }
  if (k == 0)
    return NULL;

  /* dense row is not reduced to zero, thus we have found a new pivot row */
  src_t *row = (src_t *)malloc((2*nc+2)*sizeof(src_t));
  nelts_t ctr = 2;
  for (j = 0; j < nc; ++j) {
    if (dr[j] != 0) {
      dr[j] = dr[j] % mod;
      if (dr[j] != 0) {
        row[ctr++] = (src_t)j;
        row[ctr++] = (src_t)dr[j];
      }
    }
  } 
  /* fix memory */
  row[0]  = 0;
  row[1]  = ctr;
  row     = realloc(row, row[1]*sizeof(src_t));

  if (row[3] != 1)
    row = normalize_new_pivot(row, mod);

  return row;
}

#if OLD_POLY_REPRESENTATION
static inline src_t *reduce_lower_by_upper_rows_offset_true_columns(src_t *row, const smc_t * pivs, const gb_t *basis)
{
  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  /* const nelts_t nc  = pivs->ncl + pivs->ncr; */
  const nelts_t nc  = pivs->ncr;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row[1];
  nelts_t i, j, k;
  const nelts_t shift = lc;
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
  for (nelts_t j=1; j<row[0]; j=j+2) {
    bft[row[j]]  = (bf_t)row[j+1];
  }
  free(row);
  row = NULL;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul;
  /* printf("go in\n"); */
  reduce_again:
  for (k=0, i=lc, lc=0; i<nc; ++i) {
    if (bft[i] != 0)
      bft[i]  = bft[i] % mod;
    if (bft[i] == 0)
      continue;

    const src_t *red  = pivs->row[i];
    /* printf("pivs %p, pivs->row %p, pivs->row[%u] %p\n", pivs, pivs->row, i,pivs->row[i]); */
    if (red == NULL) {
      if (!k)
        lc  = i;
      k++;
      continue;
    }
    const cf_t *cf    = basis->cf[red[0]]; 
    const nelts_t sz  = basis->nt[red[0]]+2;
    /* printf("sz %u | red[1] %u\n", sz, red[1]);
     * printf("red has values at the following column positions:\n"); */
    mul = (bf_t)(mod) - bft[i]; /* it is already reduced modulo mod */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    /* printf("mod4 %u\n", mod4); */
    /* printf("red[0] %u | red[1] %u\n", red[0], red[1]); */
    for (j=2; j<red[1]; ++j) {
      /* printf("red[%u] %u == %u lc || red[0] %u ?\n", j, red[j], lc, red[0]); */
      bft[red[j]]   +=  (bf_t)(mul) * cf[j-2];
    }
    /* bft[lc] = 0; */
    for (; j<sz; j+=2) {
      /* printf("%u * %u = %lu\n", mul, cf[j-2],(bf_t)(mul) * cf[j-2]);
       * printf("%u * %u = %lu\n", mul, cf[j-1],(bf_t)(mul) * cf[j-1]); */
      bft[red[j]]   +=  (bf_t)(mul) * cf[j-2];
      bft[red[j+1]] +=  (bf_t)(mul) * cf[j-1];
    }
    if (k && pivs->row[lc])
      goto reduce_again;
  }
  if (k == 0) {
    free(bf);
    return row;
  }

  if (k && pivs->row[lc])
    goto reduce_again;

  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  /* allocate maximal possibly needed memory */
  row = (src_t *)malloc((2*k+1)*sizeof(src_t));
  nelts_t ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf[j] = bf[j] % mod;
    if (bf[j] != 0) {
      row[ctr++] = j+shift;
      row[ctr++] = (src_t)bf[j];
    }
  } 
  /* fix memory */
  row[0] = ctr;
  row    = realloc(row, row[0]*sizeof(src_t));
#if newred
  printf("row->sz %u\n", row[0]);
#endif

  if (row[2] != 1)
    row = normalize_row_c(row, mod);

  free(bf);
  return row;
}
#endif

static inline src_t *reduce_lower_by_upper_rows_offset_c_block(src_t *row, const smc_t * pivs,
    const nelts_t offset)
{
  if (row[1] >= pivs->ncl) {
    return row;
  }

  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row[1];
  nelts_t j;
  const nelts_t shift = lc;
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
  for (nelts_t j=1; j<row[0]; j=j+2) {
    bft[row[j]]  = (bf_t)row[j+1];
  }
  free(row);
  row = NULL;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul;
  /* printf("go in\n"); */
  /* we do not have all pivots but only one block of them, so we cannot
   * test against the number of left columns! */
  /* printf("new round %u + %u\n", pivs->row[0][2], pivs->nr); */
  while (lc < pivs->row[0][2] + pivs->nr) {
    /* printf("lc %u | up to %u\n", lc, pivs->row[0][1] + pivs->nr); */
    const src_t *red = pivs->row[lc-offset];
    /* printf("lc %u\n", lc);
     * printf("red %p --> red[0] %u -- red[1] %u\n", red, red[0], red[1]); */
#if newred
    printf("lc %u\n", lc);
    printf("red %p --> red[0] %u\n", red, red[0]);
    printf("we reduce with:\n");
    for (int ii=1; ii<red[0]; ++ii) {
      printf("%u | %u -- ", red[2*ii+1], red[2*ii]);
    }
    printf("\n");
#endif
    mul = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
    /* if (mul == 0)
     *   printf("mul is zero!\n"); */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    /* printf("mod4 %u\n", mod4); */
    /* printf("red[0] %u | red[1] %u\n", red[0], red[1]); */
    for (j=2; j<2*red[1]; j=j+2) {
      /* printf("red[%u] %u == %u lc || red[0] %u ?\n", j, red[j], lc, red[0]); */
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
    }
    /* bft[lc] = 0; */
    for (; j<red[0]; j+=4) {
      /* printf("%u * %u = %lu\n", mul, red[j+1],(bf_t)(mul) * red[j+1]);
       * printf("%u * %u = %lu\n", mul, red[j+3],(bf_t)(mul) * red[j+3]); */
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
      /* bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
       * bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
       * bft[red[j+8]] +=  (bf_t)(mul) * red[j+9];
       * bft[red[j+10]] +=  (bf_t)(mul) * red[j+11];
       * bft[red[j+12]] +=  (bf_t)(mul) * red[j+13];
       * bft[red[j+14]] +=  (bf_t)(mul) * red[j+15]; */
    }
    /* get new pivot element in row */
    for (j=lc+1-shift; j<nc-shift; ++j) {
      if (bf[j] != 0) {
        bf[j] = MODP(bf[j], mod);
      }
      if (bf[j] != 0) {
        break;
      }
    }
#if newred
    printf("buffer:  ");
    for (nelts_t i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* next pivot */
    lc  = j+shift;
    if (lc >= nc) {
#if newred
      printf("NULL RETURNED!\n");
#endif
      free(bf);
      return row; /* i.e. row = NULL */
    }
  }

  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  /* allocate maximal possibly needed memory */
  row = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t));
  nelts_t ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf[j] = bf[j] % mod;
    if (bf[j] != 0) {
      row[ctr++] = j+shift;
      row[ctr++] = (src_t)bf[j];
    }
  } 
  /* fix memory */
  row[0] = ctr;
  row    = realloc(row, row[0]*sizeof(src_t));
#if newred
  printf("row->sz %u\n", row[0]);
#endif

  if (row[2] != 1)
    row = normalize_row_c(row, mod);

  free(bf);
  return row;
}

static inline src_t *reduce_lower_by_upper_rows_offset_c(src_t *row, const smc_t * pivs)
{
  if (row[1] >= pivs->ncl) {
    return row;
  }

  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  const nelts_t ncl = pivs->ncl;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row[1];
  nelts_t j;
  const nelts_t shift = lc;
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
  for (nelts_t j=1; j<row[0]; j=j+2) {
    bft[row[j]]  = (bf_t)row[j+1];
  }
  free(row);
  row = NULL;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul;
  /* printf("go in\n"); */
  while (lc < ncl) {
    const src_t *red = pivs->row[lc];
    /* printf("lc %u\n", lc);
     * printf("red %p --> red[0] %u -- red[1] %u\n", red, red[0], red[1]); */
#if newred
    printf("lc %u\n", lc);
    printf("red %p --> red[0] %u\n", red, red[0]);
    printf("we reduce with:\n");
    for (int ii=1; ii<red[0]; ++ii) {
      printf("%u | %u -- ", red[2*ii+1], red[2*ii]);
    }
    printf("\n");
#endif
    mul = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
    /* if (mul == 0)
     *   printf("mul is zero!\n"); */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    /* printf("mod4 %u\n", mod4); */
    /* printf("red[0] %u | red[1] %u\n", red[0], red[1]); */
    for (j=2; j<2*red[1]; j=j+2) {
      /* printf("red[%u] %u == %u lc || red[0] %u ?\n", j, red[j], lc, red[0]); */
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
    }
    /* bft[lc] = 0; */
    for (; j<red[0]; j+=4) {
      /* printf("%u * %u = %lu\n", mul, red[j+1],(bf_t)(mul) * red[j+1]);
       * printf("%u * %u = %lu\n", mul, red[j+3],(bf_t)(mul) * red[j+3]); */
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
      /* bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
       * bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
       * bft[red[j+8]] +=  (bf_t)(mul) * red[j+9];
       * bft[red[j+10]] +=  (bf_t)(mul) * red[j+11];
       * bft[red[j+12]] +=  (bf_t)(mul) * red[j+13];
       * bft[red[j+14]] +=  (bf_t)(mul) * red[j+15]; */
    }
    /* get new pivot element in row */
    for (j=lc+1-shift; j<nc-shift; ++j) {
      if (bf[j] != 0) {
        bf[j] = MODP(bf[j], mod);
      }
      if (bf[j] != 0) {
        break;
      }
    }
#if newred
    printf("buffer:  ");
    for (nelts_t i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* next pivot */
    lc  = j+shift;
    if (lc >= nc) {
#if newred
      printf("NULL RETURNED!\n");
#endif
      free(bf);
      return row; /* i.e. row = NULL */
    }
  }

  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  /* allocate maximal possibly needed memory */
  row = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t));
  nelts_t ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf[j] = bf[j] % mod;
    if (bf[j] != 0) {
      row[ctr++] = j+shift;
      row[ctr++] = (src_t)bf[j];
    }
  } 
  /* fix memory */
  row[0] = ctr;
  row    = realloc(row, row[0]*sizeof(src_t));
#if newred
  printf("row->sz %u\n", row[0]);
#endif

  if (row[2] != 1)
    row = normalize_row_c(row, mod);

  free(bf);
  return row;
}

static inline void reduce_many_lower_by_upper_rows_offset_c(src_t **rows, 
    const nelts_t nr, const smc_t *pivs)
{
  for (nelts_t i=0; i<nr; ++i)
    rows[i] = reduce_lower_by_upper_rows_offset_c(rows[i], pivs);
}

static inline src_t *reduce_lower_by_upper_rows_c(src_t *row, const smc_t *pivs)
{
  if (row[1] >= pivs->ncl)
    return row;

  /* printf("row %p | lc %u\n", row, row[1]); */
  const cf_t mod    = pivs->mod;
  const nelts_t nc  = pivs->ncl + pivs->ncr;
  const nelts_t ncl = pivs->ncl;
  /* lc is the lead column of the row to be reduced */
  nelts_t lc  = row[1];
  nelts_t j;
  const nelts_t shift = lc;
  /* store row in dense format using bf_t data type for delayed modulus
   * computations */
  bf_t *bf = (bf_t *)calloc(nc-shift, sizeof(bf_t));
  bf_t *bft = bf-shift;
  for (nelts_t j=1; j<row[0]; j=j+2) {
    bft[row[j]]  = (bf_t)row[j+1];
  }
  free(row);
  row = NULL;

#if newred
  printf("buffer:  ");
  for (i=0; i<nc-shift; ++i)
    printf("%lu ", bf[i]);
  printf("\n");
#endif

  bf_t mul;
  while (lc < ncl) {
    const src_t *red = pivs->row[lc];
#if newred
    printf("lc %u\n", lc);
    printf("red %p --> red[0] %u\n", red, red[0]);
    printf("we reduce with:\n");
    for (int ii=1; ii<red[0]; ++ii) {
      printf("%u | %u -- ", red[2*ii+1], red[2*ii]);
    }
    printf("\n");
#endif
    mul = (bf_t)(mod) - bft[lc]; /* it is already reduced modulo mod */
    /* loop unrolling by 2, test length in advance and probably add the useless
     * reduction on the pivot entry in order to get a even loop length */
    j = (red[0]-1)/2;
#if 0
    j = j & 1 ? 3 : 1;
    for (; j<red[0]; j+=4) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
    }
#else
    /*
    const nelts_t mod8  = j % 8;
    for (j=1; j<2*mod8; j=j+2) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
    }
    for (; j<red[0]; j+=16) {
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
      bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
      bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
      bft[red[j+8]] +=  (bf_t)(mul) * red[j+9];
      bft[red[j+10]] +=  (bf_t)(mul) * red[j+11];
      bft[red[j+12]] +=  (bf_t)(mul) * red[j+13];
      bft[red[j+14]] +=  (bf_t)(mul) * red[j+15];
    }
    */
    /*
     *
     *
     *
     *
     * TODO: For each upper row we can store the offset when constructing the
     * row.
     *
     *
     *
     */
    const nelts_t mod4  = j % 4;
    /* printf("mod4 %u\n", mod4); */
    for (j=1; j<2*mod4; j=j+2) {
      /* printf("red[%u] %u == %u lc || red[0] %u ?\n", j, red[j], lc, red[0]); */
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
    }
    /* bft[lc] = 0; */
#define AVX_TEST 0
#if AVX_TEST
    bf_t *tmp = (bf_t *)malloc(4 * sizeof(bf_t));
    __m256i mulv  = _mm256_set1_epi64x(mul);
    for (; j<red[0]; j+=8) {
      /* printf("---in ? %u j | %u red[0]\n", j, red[0]); */
      __m256i bftv  = _mm256_set_epi64x(
          bft[red[j]],
          bft[red[j+2]],
          bft[red[j+4]],
          bft[red[j+6]]);
      __m256i redv  = _mm256_set_epi64x(
          red[j+1],
          red[j+3],
          red[j+5],
          red[j+7]);
      bftv  = _mm256_add_epi64(bftv,_mm256_mul_epu32(mulv,redv));
      _mm256_storeu_si256((bf_t *)tmp, bftv);
      bft[red[j]]   = tmp[3];
    /* printf("bft[%u] %u = %lu\n", red[j], j, bft[red[j]]); */
      bft[red[j+2]] = tmp[2];
      bft[red[j+4]] = tmp[1];
      bft[red[j+6]] = tmp[0];
    }
    free(tmp);
#else
    for (; j<red[0]; j+=8) {
      /* printf("in ? %u j | %u red[0]\n", j, red[0]); */
      bft[red[j]]   +=  (bf_t)(mul) * red[j+1];
    /* printf("bft[%u] %u = %lu\n", red[j], j, bft[red[j]]); */
      bft[red[j+2]] +=  (bf_t)(mul) * red[j+3];
      bft[red[j+4]] +=  (bf_t)(mul) * red[j+5];
      bft[red[j+6]] +=  (bf_t)(mul) * red[j+7];
    }
#endif
#endif
    /* get new pivot element in row */
    for (j=lc+1-shift; j<nc-shift; ++j) {
      if (bf[j] != 0) {
        bf[j] = MODP(bf[j], mod);
      }
      if (bf[j] != 0) {
        break;
      }
    }
#if newred
    printf("buffer:  ");
    for (nelts_t i=0; i<nc-shift; ++i)
      printf("%lu ", bf[i]);
    printf("\n");
#endif
    /* next pivot */
    lc  = j+shift;
    if (lc >= nc) {
#if newred
      printf("NULL RETURNED!\n");
#endif
      free(bf);
      return row; /* i.e. row = NULL */
    }
  }

  /* move data from buffer back to sparse row
   * normalize row if needed */
#if newred
  printf("nc - lc = %u - %u\n", nc, lc);
#endif
  /* allocate maximal possibly needed memory */
  row = (src_t *)malloc((2*(nc-shift)-1)*sizeof(src_t));
  nelts_t ctr = 1;
  for (j=lc-shift; j<nc-shift; ++j) {
    bf[j] = bf[j] % mod;
    if (bf[j] != 0) {
      row[ctr++] = j+shift;
      row[ctr++] = (src_t)bf[j];
    }
  } 
  /* fix memory */
  row[0] = ctr;
  row    = realloc(row, row[0]*sizeof(src_t));
#if newred
  printf("row->sz %u\n", row[0]);
#endif

  /* no normalization in block reduction style,
   * we normalize first when all upper block reductions are done */
  /* if (row[2] != 1)
   *   row = normalize_row_c(row, mod); */

  free(bf);
  return row;
}

#if OLD_POLY_REPRESENTATION
static inline mat_gb_meta_data_t *generate_matrix_meta_data(
    const smc_t *AB, const smc_t *CD, const nelts_t bs)
{
  mat_gb_meta_data_t *mat =
    (mat_gb_meta_data_t *)malloc(sizeof(mat_gb_meta_data_t));

  mat->mod    = AB->mod;
  mat->bs     = (nelts_t)bs;

#define ADAPTIVE_BLOCK_SIZE 1

#if ADAPTIVE_BLOCK_SIZE
  if (AB->nr < mat->bs) {
    mat->bs = AB->nr > 0 ? AB->nr : (AB->ncl + AB->ncr);
  } else {
    nelts_t factor  = AB->nr / mat->bs + 1;
    mat->bs = AB->nr / factor;
  }
  /* make block size a multiple of 8 */
  mat->bs +=  (8 - mat->bs % 8);
  /* printf("%u\n",mat->bs); */

#endif
  mat->nc     = AB->ncl + AB->ncr;
  mat->nr     = AB->nr + CD->nr;
  mat->nc_AC  = mat->nr_AB = AB->nr;
  mat->nc_BD  = AB->ncr;
  mat->nr_CD  = CD->nr;
  mat->nrb_AB = (nelts_t) ceil((float)mat->nr_AB/(float)mat->bs);
  mat->nrb_CD = (nelts_t) ceil((float)mat->nr_CD/(float)mat->bs);
  mat->ncb_AC = (nelts_t) ceil((float)mat->nc_AC/(float)mat->bs);
  mat->ncb_BD = (nelts_t) ceil((float)mat->nc_BD/(float)mat->bs);
  mat->ncb    = mat->ncb_AC + mat->ncb_BD;
  mat->rk     = mat->nr;

  return mat;
}


/* Paramter "start" has one single usage: If generating AB we can 
 * leave out the lead term of the poly since it need not be stored
 * in the A blocks. This also allows us to use the very same reduction
 * resp. update functions for interreducing AB and reducing CD with AB.
 * Thus generation of AB always calls this procedure with "start=1" and
 * the generation of CD is done via "start=0". */
static inline void write_poly_to_matrix(mat_gb_block_t *mat,
    const mat_gb_meta_data_t *meta, const nelts_t idx, const nelts_t start,
    const mpp_t *mpp, const gb_t *basis, const ht_t *ht)
{
  nelts_t i;
  nelts_t cp;
  mat_gb_block_t *bl;

  /* update length entries first */
  for (i=0; i<meta->ncb; ++i) {

    mat[i].len[idx+1] = mat[i].len[idx];
  }

  /* printf("mpp->mul %lu | bi %u | nt %u\n", mpp->mul, mpp->bi, basis->nt[mpp->bi]); */
  for (i=start; i<basis->nt[mpp->bi]; ++i) {
    /* printf("i %u\n", i); */
    cp  = ht->idx[find_in_hash_table_product(mpp->mul,
        basis->eh[mpp->bi][i], ht)];

    if (cp < meta->nc_AC) {
      /* printf("C ");
       * printf("block %u\n", cp /meta->bs); */
      bl  = mat + (cp / meta->bs);
    } else {
      cp  = cp - meta->nc_AC;
      bl  = mat + (meta->ncb_AC + cp / meta->bs);
      /* printf("D ");
       * printf("block %u\n", meta->ncb_AC + cp /meta->bs); */
    }
    /* printf("coeff %u\n", basis->cf[mpp->bi][i]); */
    bl->val[bl->len[idx+1]] = (cf_t)basis->cf[mpp->bi][i];
    bl->pos[bl->len[idx+1]] = (bs_t)(cp % meta->bs);
    /* printf("idx+1 %u bl->val[%u] = %u | bl->pos[%u] = %u\n", idx+1,
     *  bl->len[idx+1],bl->val[bl->len[idx+1]],bl->len[idx+1],bl->pos[bl->len[idx+1]]); */
    bl->len[idx+1]++;
  }
}


static inline void adjust_block_row_types_including_dense_righthand(
    mat_gb_block_t *mat, const mat_gb_meta_data_t *meta)
{
  nelts_t i,j,k;

  const nelts_t bs_square = (nelts_t)meta->bs * meta->bs;
  
  /* lefthand side, always kept sparse */
  for (i=0; i<meta->ncb_AC; ++i) {
    /* printf("i %u\n", i); */
    /* sparse to dense ? */
    if (mat[i].len != NULL) {
      /* printf("%u len %u\n", i, mat[i].len[mat[i].nr]); */
      if (mat[i].len[mat[i].nr] > 0) {
        /* if (mat[i].len[mat[i].nr] > bs_square/2)
         *   printf("DENSE\n"); */
        /* printf("len %u\n", mat[i].len[mat[i].nr]); */
        mat[i].pos =
          realloc(mat[i].pos, mat[i].len[mat[i].nr] * sizeof(bs_t));
        mat[i].val =
          realloc(mat[i].val, mat[i].len[mat[i].nr] * sizeof(cf_t));
      } else {
        free(mat[i].len);
        mat[i].len  = NULL;
        free(mat[i].pos);
        mat[i].pos  = NULL;
        free(mat[i].val);
        mat[i].val  = NULL;
      }
    }
  }
  /* righthand side, probably to be shifted to dense format */
  for (i=meta->ncb_AC; i<meta->ncb; ++i) {
    /* printf("i %u\n", i); */
    /* sparse to dense ? */
    if (mat[i].len != NULL) {
      /* printf("%u len %u\n", i, mat[i].len[mat[i].nr]); */
      if (mat[i].len[mat[i].nr] > 0) {
        if (mat[i].len[mat[i].nr] > bs_square/2) {
          /* printf("make dense\n"); */
          cf_t *val = (cf_t *)calloc(bs_square, sizeof(cf_t));
          for (j=0; j<mat[i].nr; ++j) {
            for (k=mat[i].len[j]; k<mat[i].len[j+1]; ++k) {
              val[j*meta->bs+mat[i].pos[k]]  = mat[i].val[k];
            }
          }
          free(mat[i].len);
          mat[i].len  = NULL;
          free(mat[i].pos);
          mat[i].pos  = NULL;
          free(mat[i].val);
          mat[i].val  = val;
          /* for (nelts_t ii=0; ii<meta->bs; ++ii) {
           *   for (nelts_t jj=0; jj<meta->bs; ++jj) {
           *     printf("%u ",mat[i].val[ii*meta->bs+jj]);
           *   }
           *   printf("\n");
           * } */
          continue;
        } else {
          /* printf("len %u\n", mat[i].len[mat[i].nr]); */
          mat[i].pos =
            realloc(mat[i].pos, mat[i].len[mat[i].nr] * sizeof(bs_t));
          mat[i].val =
            realloc(mat[i].val, mat[i].len[mat[i].nr] * sizeof(cf_t));
          continue;
        }
      } else {
        free(mat[i].len);
        mat[i].len  = NULL;
        free(mat[i].pos);
        mat[i].pos  = NULL;
        free(mat[i].val);
        mat[i].val  = NULL;
        continue;
      }
    }
  }
}

static inline void adjust_block_row_types_including_dense(
    mat_gb_block_t *mat, const mat_gb_meta_data_t *meta)
{
  nelts_t i,j,k;

  const nelts_t bs_square = (nelts_t)meta->bs * meta->bs;
  
  for (i=0; i<meta->ncb; ++i) {
    /* printf("i %u\n", i); */
    /* sparse to dense ? */
    if (mat[i].len != NULL) {
      /* printf("%u len %u\n", i, mat[i].len[mat[i].nr]); */
      if (mat[i].len[mat[i].nr] > 0) {
        if (mat[i].len[mat[i].nr] > bs_square/2) {
          /* printf("make dense\n"); */
          cf_t *val = (cf_t *)calloc(bs_square, sizeof(cf_t));
          for (j=0; j<mat[i].nr; ++j) {
            for (k=mat[i].len[j]; k<mat[i].len[j+1]; ++k) {
              val[j*meta->bs+mat[i].pos[k]]  = mat[i].val[k];
            }
          }
          free(mat[i].len);
          mat[i].len  = NULL;
          free(mat[i].pos);
          mat[i].pos  = NULL;
          free(mat[i].val);
          mat[i].val  = val;
          /* for (nelts_t ii=0; ii<meta->bs; ++ii) {
           *   for (nelts_t jj=0; jj<meta->bs; ++jj) {
           *     printf("%u ",mat[i].val[ii*meta->bs+jj]);
           *   }
           *   printf("\n");
           * } */
          continue;
        } else {
          /* printf("len %u\n", mat[i].len[mat[i].nr]); */
          mat[i].pos =
            realloc(mat[i].pos, mat[i].len[mat[i].nr] * sizeof(bs_t));
          mat[i].val =
            realloc(mat[i].val, mat[i].len[mat[i].nr] * sizeof(cf_t));
          continue;
        }
      } else {
        free(mat[i].len);
        mat[i].len  = NULL;
        free(mat[i].pos);
        mat[i].pos  = NULL;
        free(mat[i].val);
        mat[i].val  = NULL;
        continue;
      }
    }
  }
}

static inline void adjust_block_row_types(mat_gb_block_t *mat,
    const mat_gb_meta_data_t *meta)
{
  for (size_t i = 0; i < meta->ncb; ++i) {
    /* sparse to dense ? */
    if (mat[i].len != NULL) {
      if (mat[i].len[mat[i].nr] > 0) {
        mat[i].pos =
          realloc(mat[i].pos, mat[i].len[mat[i].nr] * sizeof(bs_t));
        mat[i].val =
          realloc(mat[i].val, mat[i].len[mat[i].nr] * sizeof(cf_t));
      } else {
        free(mat[i].len);
        mat[i].len  = NULL;
        free(mat[i].pos);
        mat[i].pos  = NULL;
        free(mat[i].val);
        mat[i].val  = NULL;
      }
    }
  }
}

static inline mat_gb_block_t *generate_mat_gb_block(
    const mat_gb_meta_data_t *meta, const nelts_t nr)
{
  mat_gb_block_t *bl  = (mat_gb_block_t *)malloc(sizeof(mat_gb_block_t));
  const nelts_t bs_square = (nelts_t)meta->bs * meta->bs;

  bl->len = (nelts_t *)calloc((meta->bs + 1), sizeof(nelts_t));
  bl->pos = (bs_t *)malloc(bs_square * sizeof(bs_t));
  bl->val = (cf_t *)malloc(bs_square * sizeof(cf_t));
  bl->nr  = nr;

  return bl;
}

static inline void initialize_mat_gb_block(mat_gb_block_t *bl,
    const mat_gb_meta_data_t *meta, const nelts_t nr)
{
  const nelts_t bs_square = (nelts_t)meta->bs * meta->bs;

  bl->len = (nelts_t *)calloc((meta->bs + 1), sizeof(nelts_t));
  bl->pos = (bs_t *)malloc(bs_square * sizeof(bs_t));
  bl->val = (cf_t *)malloc(bs_square * sizeof(cf_t));
  bl->nr  = nr;
}

static inline void write_to_mat_gb_row_block_inverted_order(
    mat_gb_block_t *mat, const mat_gb_meta_data_t *meta, const nelts_t idx,
    const sel_t *sel, const gb_t *basis, const ht_t *ht)
{
  nelts_t i, j;
  
  /* mat_gb_block_t *start = mat + (idx * meta->ncb); */
  mat_gb_block_t *start = mat;
  /* printf("idx %u | ncb %u\n", idx, meta->ncb); */

  const nelts_t offset  = idx*meta->bs;
  const nelts_t max     =
    (idx+1)*meta->bs < sel->load ? meta->bs : sel->load - offset;

  for (i=0; i<meta->ncb; ++i) {
    /* printf("block allocation %u / %u\n",i, meta->ncb); */
    initialize_mat_gb_block(start+i, meta, max);
  }

  for (i=max; i>0; --i) {
    /* printf("max %u | i %u | offset %u\n", max, i, offset); */
    /* leaves out lead term 1, thus starts at position "1" (3rd parameter) */
    write_poly_to_matrix(start, meta, max-i, 1, sel->mpp+(i-1+offset), basis, ht);
    /* printf("---\n"); */
  }

  /* check len entries */
  if (max < meta->bs) {
    for (i=0; i<meta->ncb; ++i) {
      for (j=max+1; j<meta->bs+1; ++j) {
        start[i].len[j] = start[i].len[max];
      }
    }
  }

  /* check density of blocks */
  /* initially we keep them sparsely represented */
  adjust_block_row_types(start, meta);
}

static inline void write_to_src_mat_row_block(src_t **mat,
    const mat_gb_meta_data_t *meta, const nelts_t idx,
    const sel_t *sel, const gb_t *basis)
{
  const nelts_t offset  = idx*meta->bs;
  const nelts_t max     =
    (idx+1)*meta->bs < sel->load ? meta->bs : sel->load - offset;
  src_t **row   = mat + offset;

  poly_to_sparse_compact_matrix_row_test(sel->mpp+(idx+offset),
      max, basis, row);
}

static inline void write_to_mat_gb_row_block(mat_gb_block_t *mat,
    const mat_gb_meta_data_t *meta, const nelts_t idx, const nelts_t ci,
    const sel_t *sel, const gb_t *basis, const ht_t *ht)
{
  nelts_t i;
  
  mat_gb_block_t *start   = mat + (ci * meta->ncb);

  const nelts_t offset  = idx*meta->bs;
  const nelts_t max     =
    (idx+1)*meta->bs < sel->load ? meta->bs : sel->load - offset;

  for (i=0; i<meta->ncb; ++i) {
    initialize_mat_gb_block(start+i, meta, max);
  }


  for (i=0; i<max; ++i) {
    write_poly_to_matrix(start, meta, i, 0, sel->mpp+(i+offset), basis, ht);
  }

  /* check len entries */
  if (max < meta->bs) {
    for (i=0; i<meta->ncb; ++i) {
      for (nelts_t j=max+1; j<meta->bs+1; ++j) {
        start[i].len[j] = start[i].len[max];
      }
    }
  }

  /* initially we keep them sparsely represented */
  adjust_block_row_types(start, meta);
  /* adjust_block_row_types_including_dense(start, meta); */
}

static inline void invert_first_block(mat_gb_block_t *mat,
    const mat_gb_meta_data_t *meta)
{
  /* sparse */
  if (mat[0].len != NULL) {
    for (size_t i = 0; i < meta->bs; ++i) {
      if (mat[0].len[i]<mat[0].len[i+1]) {
        mat[0].val[mat[0].len[i]] = meta->mod - mat[0].val[mat[0].len[i]];
      }
    }
  } else { /* dense */
    for (size_t i = 0; i < meta->bs; ++i) {
      mat[0].val[i+i*meta->bs] = meta->mod - mat[0].val[i+i*meta->bs];
    }
  }
}

static inline mat_gb_block_t *generate_mat_gb_upper_row_block(
    const nelts_t idx, const mat_gb_meta_data_t *meta, const gb_t *basis,
    const spd_t *spd, const ht_t *ht)
{
  mat_gb_block_t *mat  = (mat_gb_block_t *)malloc(
      meta->ncb * sizeof(mat_gb_block_t));
  /* printf("ncb %u\n", meta->ncb);
   * printf("idx %u\n", idx);
   * printf("selu load %u\n", spd->selu->load);
   * printf("sell load %u\n", spd->sell->load); */

  /* note the inverted row order in each block: we make the first block the
   * unit matrix block, thus we have to reduce rows with smaller column
   * pivot entries by rows with higher column pivot. */
  write_to_mat_gb_row_block_inverted_order(mat, meta, idx, spd->selu,
      basis, ht);

  /* printf("----------------------------------------------------\n");
   * for (int ii=0; ii<meta->ncb; ++ii) {
   *   if (mat[ii].len != NULL) {
   *     for (int jj=0; jj<mat[ii].nr; ++jj)
   *       printf("%u ", mat[ii].len[jj]);
   *     printf("\n");
   *     printf("matlen %p\n", mat[ii].len);
   *     printf("matpos %p\n", mat[ii].pos);
   *     printf("matval %p\n", mat[ii].val);
   *     printf("nc_AC %u\n", meta->nc_AC);
   *     printf("mat[%u].val[0] = %u\n", ii, mat[ii].val[0]);
   *   }
   * }
   * printf("----------------------------------------------------\n"); */
  /* invert_first_block(mat, meta);
   * printf("----------------------------------------------------\n");
   * for (int ii=0; ii<meta->ncb; ++ii) {
   *   if (mat[ii].len != NULL) {
   *     for (int jj=0; jj<mat[ii].nr; ++jj)
   *       printf("%u ", mat[ii].len[jj]);
   *     printf("\n");
   *     printf("matlen %p\n", mat[ii].len);
   *     printf("matpos %p\n", mat[ii].pos);
   *     printf("matval %p\n", mat[ii].val);
   *     printf("nc_AC %u\n", meta->nc_AC);
   *     printf("mat[%u].val[0] = %u\n", ii, mat[ii].val[0]);
   *   }
   * }
   * printf("----------------------------------------------------\n"); */

  return mat;
}

static inline mat_gb_block_t *generate_mat_gb_row_block(
    const nelts_t idx, const mat_gb_meta_data_t *meta, const gb_t *basis,
    const spd_t *spd, const ht_t *ht)
{
  mat_gb_block_t *mat  = (mat_gb_block_t *)malloc(
      meta->ncb * sizeof(mat_gb_block_t));

  write_to_mat_gb_row_block(mat, meta, idx, 0, spd->selu, basis, ht);

  return mat;
}

static inline mat_gb_block_t *generate_mat_gb_upper(
    const mat_gb_meta_data_t *meta, const gb_t *basis, const spd_t *spd, 
    const ht_t *ht, const int t)
{
  mat_gb_block_t *mat  = (mat_gb_block_t *)malloc(
      (meta->nrb_AB * meta->ncb * sizeof(mat_gb_block_t)));
  /* printf("ncb %u | nrb_CD %u\n", meta->ncb, meta->nrb_CD); */

  /* printf("meta %u\n", meta->nrb_CD); */
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
      for (nelts_t i=0; i<meta->nrb_AB; ++i) {
        #pragma omp task
        write_to_mat_gb_row_block(mat, meta, i, i, spd->selu, basis, ht);
      }
    }
  }
  return mat;
}

static inline mat_gb_block_t *generate_mat_gb_lower(
    const mat_gb_meta_data_t *meta, const gb_t *basis, const spd_t *spd, 
    const ht_t *ht, const int t)
{
  mat_gb_block_t *mat  = (mat_gb_block_t *)malloc(
      (meta->nrb_CD * meta->ncb * sizeof(mat_gb_block_t)));
  /* printf("ncb %u | nrb_CD %u\n", meta->ncb, meta->nrb_CD); */

  /* printf("meta %u\n", meta->nrb_CD); */
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
      for (nelts_t i=0; i<meta->nrb_CD; ++i) {
        #pragma omp task
        write_to_mat_gb_row_block(mat, meta, i, i, spd->sell, basis, ht);
      }
    }
  }

  /* printf("%p\n", mat[0].len);
   * printf("%p\n", mat[0].pos);
   * printf("%p\n", mat[0].val);
   * printf("nr %u\n", mat[0].nr); */
  return mat;
}

static inline void free_mat_gb_block(mat_gb_block_t *bl)
{
  free(bl->len);
  bl->len = NULL;
  free(bl->pos);
  bl->pos = NULL;
  free(bl->val);
  bl->val = NULL;
  bl->nr  = 0;
}

smc_t *convert_mat_gb_to_smc_offset_format(const mat_gb_block_t **om,
    const mat_gb_meta_data_t *meta, const int t);

smc_t *convert_mat_gb_to_smc_format(const mat_gb_block_t *om,
    const mat_gb_meta_data_t *meta, const int t);

#endif


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
ri_t reduce_gbla_matrix(mat_t * mat, int verbose, int nthreads);

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
ri_t reduce_gbla_matrix_keep_A(mat_t * mat, int verbose, int nthreads);
#endif
