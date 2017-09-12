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
 * \file elimination.h
 * \brief Different Gaussian Elimination methods
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GBLA_ELIMINATION_H
#define GBLA_ELIMINATION_H

#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <x86intrin.h>
#include "../config.h"
/* #include "mapping.h" */
#include "types.h"
#include "gbla-matrix.h"
#if GBLA_WITH_FFLAS
#include "../draft/echelon.h"
#endif

// opencl-stuff
#if __GBLA_HAVE_OPENCL
#include "cl-timing.h"
#include "cl-helper.h"
#endif 

#undef AVX
#undef SSE
//#define NOSSETEST
#define NOSSE
//#define AVX
//#define NOSSE2

#define DEBUG_NEW_ELIM  0
#define DEBUG_ECHELONIZE  0
#define COUNT_REDS  0

#if COUNT_REDS
static unsigned long nreductions;
#endif

// global variable storing the number of pre_eliminated rows of D done in
// sequential.
static ri_t global_pre_elim;

/**
 * \brief Structure of a waiting list element
 */
typedef struct wle_t {
	ri_t idx; /*!<  row index */
	ri_t lp;  /*!<  row index of last pivot that idx was reduced with */
} wle_t;

/**
 * \brief Structure defining the waiting list or queue for parallel dense
 * echelonization of the D part.
 *
 * \note We use a mutiline concept for the row index lists.
 */
typedef struct wl_t {
	wle_t *list;  /*!<  row indices of rows in waiting list */
	ri_t sidx;    /*!<  smallest row index in rows */
	ri_t slp;     /*!<  last pivot row index by which sidx is
									reduced by already */
	ri_t sz;      /*!<  size of waiting list */
} wl_t;

/**
 * \brief Comparison function for qsort for waiting list
 *
 * \param waiting list element a
 *
 * \param waiting list element b
 *
 * \return *b.idx - *a.idx
 */
int cmp_wle(const void *a, const void *b) ;


/**
 * \brief Initializes wide (=uint64_t) blocks for cumulative multiplications and
 * additions.
 *
 * \param wide block wide_block
 */
static inline void init_wide_blocks(re_l_t ***wide_block)
{
  int ret;
  re_l_t **wb = (re_l_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(re_l_t *));
  uint64_t size = __GBLA_SIMD_BLOCK_SIZE * sizeof(re_l_t);
  bi_t i;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    do {
      ret =posix_memalign((void **)&wb[i], 16, size);
    } while (ret != 0);
  }
  *wide_block  = wb;
}

/**
 * \brief Initializes wide (=uint64_t) row for cumulative multiplications and
 * additions.
 *
 * \param wide row wide_row
 *
 * \param lebngth of row
 */
static inline void init_wide_rows(re_l_t **wide_row, ci_t length)
{
  int ret;
  re_l_t *wr = *wide_row;
  do {
    ret = posix_memalign((void **)&wr, 16, length * sizeof(re_l_t));
  } while (ret != 0);
  *wide_row  = wr;
}

/**
 * \brief Set wide block array to zero.
 *
 * \param wide block array wide_blocks
 *
 * \param length of wide_blocks array length
 */
static inline void set_wide_block_to_zero(re_l_t **wide_block, const bi_t length)
{
    bi_t k;
    for (k=0; k<length; ++k)
      memset(wide_block[k], 0, length * sizeof(re_l_t));
}

/**
 * \brief Frees memory allocated for a wide block.
 *
 * \param wide block to be freed wide_block
 */
static inline void free_wide_block(re_l_t ***wide_block)
{
  re_l_t **wb = *wide_block;
  bi_t i;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i)
    free(wb[i]);
  free(wb);
  wb  = NULL;
  *wide_block = wb;
}

/**
 * \brief Copies entries from hybrid block hybrid_block to dense wide block
 * wide_block for elimination purposes.
 *
 * \note hybrid_block is already checked to be != NULL
 *
 * \param hybrid block hybrid_block
 *
 * \param wide block wide_block
 */
static inline void copy_hybrid_to_wide_block(dbl_t ***hybrid_block_in,
    re_l_t **wide_block)
{
  bi_t i, j, k;
  dbl_t ** hybrid_block = *hybrid_block_in;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    if (hybrid_block[i] != NULL) {
      for (j=0; j<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++j) {
        if (hybrid_block[i][j].val != NULL) {
          for (k=0; k<__GBLA_SIMD_INNER_SIZE; ++k) {
            wide_block[i][j*__GBLA_SIMD_INNER_SIZE + k] =
              hybrid_block[i][j].val[k];
          }
          free(hybrid_block[i][j].val);
        }
      }
      free(hybrid_block[i]);
    }
  }
  free(hybrid_block);
  hybrid_block      = NULL;
  *hybrid_block_in  = hybrid_block;
}

/**
 * \brief Copies entries from dense wide block wide_block back to the hybrid block
 * hybrid_block after elimination purposes.
 *
 * \param wide block wide_block
 *
 * \param hybrid block hybrid_block
 */
static inline void copy_wide_to_hybrid_block(re_l_t **wide_block,
    dbl_t ***hybrid_block)
{
  bi_t i, j, k;
  dbl_t **hb  = (dbl_t **)malloc(__GBLA_SIMD_BLOCK_SIZE * sizeof(dbl_t *));
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
      //printf("wb[%d][%d] = %lu\n",i,j,wide_block[i][j]);
      if (wide_block[i][j] != 0) {
        goto not_zero;
      }
    }
  }
  free(hb);
  return;

not_zero:
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    hb[i] = NULL;
  }
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
      if (wide_block[i][j] != 0) {
        if (hb[i] == NULL) {
          hb[i] = (dbl_t *)malloc(__GBLA_SIMD_INNER_BLOCKS_PER_ROW * sizeof(dbl_t));
          for (k=0; k<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++k) {
            hb[i][k].val  = NULL;
          }
          hb[i][j/__GBLA_SIMD_INNER_SIZE].val =  (re_t *)calloc(
              __GBLA_SIMD_INNER_SIZE, sizeof(re_t));

        } else {
          if (hb[i][j/__GBLA_SIMD_INNER_SIZE].val == NULL) {
            hb[i][j/__GBLA_SIMD_INNER_SIZE].val =  (re_t *)calloc(
                __GBLA_SIMD_INNER_SIZE, sizeof(re_t));
          }
        }
        hb[i][j/__GBLA_SIMD_INNER_SIZE].val[j%__GBLA_SIMD_INNER_SIZE] =
          (re_t) wide_block[i][j];
      }
    }
  }
  *hybrid_block = hb;
  printf("--------------------------------------------------------------\n");
}

/**
 * \brief Copies entries from dense block dense_block to dense wide block
 * wide_block for elimination purposes.
 *
 * \param dense block dense_block
 *
 * \param wide block wide_block
 */
static inline void copy_dense_to_wide_block(re_t *dense_block,
    re_l_t **wide_block)
{
  bi_t i, j;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i)
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j)
      wide_block[i][j]  = dense_block[i*__GBLA_SIMD_BLOCK_SIZE+j];
}

/**
 * \brief Copies entries from sparse matrix A to dense wide row
 * wide_row for elimination purposes.
 *
 * \param sparse matrix A
 *
 * \param row index in sparse matrix idx
 *
 * \param wide row wide_row
 */
static inline void copy_sparse_to_wide_row(re_l_t *wide_row, const sm_fl_t *A,
    const ri_t idx)
{
  ci_t i = 0;
  if (A->sz[idx] > 3) {
    for (; i<A->sz[idx]-3; i = i+4) {
      wide_row[A->pos[idx][i]]    = (re_l_t)A->row[idx][i];
      wide_row[A->pos[idx][i+1]]  = (re_l_t)A->row[idx][i+1];
      wide_row[A->pos[idx][i+2]]  = (re_l_t)A->row[idx][i+2];
      wide_row[A->pos[idx][i+3]]  = (re_l_t)A->row[idx][i+3];
    }
  }
  for (; i<A->sz[idx]; ++i)
    wide_row[A->pos[idx][i]]    = (re_l_t)A->row[idx][i];
}

/**   
 * \brief Copies entries from wide dense row to a sparse matrix row
 *
 * \param sparse matrix C
 *
 * \param wide row wide_row
 *
 * \param row index in sparse matrix idx
 *
 * \param counter of nonzero entries in wide_row nze_ctr
 */
static inline void copy_wide_to_sparse_row(sm_fl_t **C_in, const re_l_t *wide_row,
    const ri_t idx, const ci_t nze_ctr)
{
  sm_fl_t *C  = *C_in;
  ci_t i;
  ci_t c_ncols = C->ncols;
  //printf("realloc %u | %u\n",idx,nze_ctr);
  //printf("%u\n",idx);
  //printf("%u | %p | %p\n",idx,C->row[idx],C->pos[idx]);
  C->row[idx] = realloc(C->row[idx], nze_ctr * sizeof(re_t));
  //printf("--\n");
  C->pos[idx] = realloc(C->pos[idx], nze_ctr * sizeof(ci_t));
  //printf("---\n");
  C->sz[idx]  = 0;
  C->buf[idx] = nze_ctr;
  /*
  for (i=0; i<c_ncols-3; i=i+4) {
    if (wide_row[i] != 0) {
      C->row[idx][C->sz[idx]] = (re_t)wide_row[i];
      C->pos[idx][C->sz[idx]] = i;
      C->sz[idx]++;
    }
    if (wide_row[i+1] != 0) {
      C->row[idx][C->sz[idx]] = (re_t)wide_row[i+1];
      C->pos[idx][C->sz[idx]] = i+1;
      C->sz[idx]++;
    }
    if (wide_row[i+2] != 0) {
      C->row[idx][C->sz[idx]] = (re_t)wide_row[i+2];
      C->pos[idx][C->sz[idx]] = i+2;
      C->sz[idx]++;
    }
    if (wide_row[i+3] != 0) {
      C->row[idx][C->sz[idx]] = (re_t)wide_row[i+3];
      C->pos[idx][C->sz[idx]] = i+3;
      C->sz[idx]++;
    }
  }
  */
  i=0;
  for (; i<c_ncols; ++i) {
    if (wide_row[i] != 0) {
      C->row[idx][C->sz[idx]] = (re_t)wide_row[i];
      C->pos[idx][C->sz[idx]] = i;
      C->sz[idx]++;
    }
  }
  *C_in = C;
}

/**
 * \brief Copies entries from dense wide block wide_block to dense block
 * dense_block for elimination purposes.
 *
 * \param wide block wide_block
 *
 * \param dense block dense_block
 */
static inline void copy_wide_to_dense_block(re_l_t **wide_block, re_t **dense_block)
{
  re_t *db = *dense_block;
  bi_t i, j;
  //printf("IN wb %p\n",wide_block);
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
      //printf("wb[%d][%d] = %lu\n",i,j,wide_block[i][j]);
      if (wide_block[i][j] != 0) {
        goto not_zero;
      }
    }
  }
  // wide_block has only zeroes in it
  free(db);
  db  = NULL;
  *dense_block  = db;
  return;

not_zero:
  if (db == NULL)
    db = (re_t *)malloc(__GBLA_SIMD_BLOCK_SIZE_RECT * sizeof(re_t));

  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i)
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j)
      db[i*__GBLA_SIMD_BLOCK_SIZE+j]  = (re_t)wide_block[i][j];
  
  *dense_block  = db;
}

/**
 * \brief Modular reduction of wide block row resp. column
 *
 * \param wide block wide_block
 *
 * \param row in wide block row_idx
 *
 * \param modulus resp. field characteristic modulus
 */
static inline void modulo_wide_block_val(re_l_t **wide_block, const bi_t row_idx,
    const mod_t modulus) {
  bi_t i;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i)
		/* wide_block[row_idx][i]  = (re_l_t)MODP(wide_block[row_idx][i] , modulus); */
    wide_block[row_idx][i]  = (re_l_t)(wide_block[row_idx][i] % modulus);
}

/**
 * \brief Modular reduction of wide block before copying back
 *
 * \param wide block wide_block
 *
 * \param modulus resp. field characteristic modulus
 */
static inline void modulo_wide_block(re_l_t ***wide_block, const mod_t modulus) {
  re_l_t **wb = *wide_block;
  bi_t i, j;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i)
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j)
			/* wb[i][j]  = (re_l_t)MODP(wb[i][j] , modulus); */
      wb[i][j]  = (re_l_t)(wb[i][j] % modulus);
  *wide_block = wb;
}

/**
 * \brief Use hybrid blocks from A to update dense blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 *
 * \param hybrid block from A block_A
 *
 * \param dense block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_hybrid_dense_rectangular(dbl_t **block_A, const re_t *block_B,
  re_l_t **wide_block)
{
  bi_t i, j, k, l;
  register re_m_t a;
  //printf("block_A %p\n",block_A);
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    if (block_A[i] == NULL) {
      continue;
    } else {
  //printf("block_A[%d] %p\n",i,block_A[i]);
      for (j=0; j<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++j) {
        if (block_A[i][j].val == NULL) {
          continue;
        } else {
  //printf("block_A[%d][%d].val %p\n",i,j,block_A[i][j].val);
          for (k=0; k<__GBLA_SIMD_INNER_SIZE; ++k) {
            if (block_A[i][j].val[k] != 0) {
              a = block_A[i][j].val[k];
              if (a != 0) {
                //printf("%u - %d\n",block_A[i][j].val[k],j*__GBLA_SIMD_INNER_SIZE+k);
                for (l=0; l<__GBLA_SIMD_BLOCK_SIZE; ++l) {
                  wide_block[i][l] +=  a *
                    (re_l_t)block_B[(j*__GBLA_SIMD_INNER_SIZE+k)*__GBLA_SIMD_BLOCK_SIZE+l];
                }
              }
            }
          }
        }
      }
    }
  }
}

#if __GBLA_COLUMN_B
/**
 * \brief Use dense blocks from A to update sparse blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 *
 * \note B is stored by column.
 *
 * \param dense block from A block_A
 *
 * \param sparse block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_dense_sparse_rectangular(const re_t *block_A, const sbl_t *block_B,
  re_l_t **wide_block)
{
#if 1
  bi_t i, j, k;
  register re_m_t a1;
  register re_m_t b1, b2, b3, b4, b5, b6, b7, b8;
  for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
    register const bi_t *posB  = block_B->pos[j];
    k = 0;
#if 0
    if (block_B->sz[j]-k > 15) {
      for (; k<block_B->sz[j]-15; k=k+16) {
        /* register const re_m_t */ b1  = block_B->val[j][k];
        /* register const re_m_t */ b2  = block_B->val[j][k+1];
        /* register const re_m_t */ b3  = block_B->val[j][k+2];
        /* register const re_m_t */ b4  = block_B->val[j][k+3];
        /* register const re_m_t */ b5  = block_B->val[j][k+4];
        /* register const re_m_t */ b6  = block_B->val[j][k+5];
        /* register const re_m_t */ b7  = block_B->val[j][k+6];
        /* register const re_m_t */ b8  = block_B->val[j][k+7];
        /* register const re_m_t */ b9  = block_B->val[j][k+8];
        /* register const re_m_t */ b10 = block_B->val[j][k+9];
        /* register const re_m_t */ b11 = block_B->val[j][k+10];
        /* register const re_m_t */ b12 = block_B->val[j][k+11];
        /* register const re_m_t */ b13 = block_B->val[j][k+12];
        /* register const re_m_t */ b14 = block_B->val[j][k+13];
        /* register const re_m_t */ b15 = block_B->val[j][k+14];
        /* register const re_m_t */ b16 = block_B->val[j][k+15];
        for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; i=i+1) {
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k]];
          wide_block[i][j]  +=  (re_l_t)a1 * b1;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+1]];
          wide_block[i][j]  +=  (re_l_t)a1 * b2;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+2]];
          wide_block[i][j]  +=  (re_l_t)a1 * b3;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+3]];
          wide_block[i][j]  +=  (re_l_t)a1 * b4;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+4]];
          wide_block[i][j]  +=  (re_l_t)a1 * b5;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+5]];
          wide_block[i][j]  +=  (re_l_t)a1 * b6;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+6]];
          wide_block[i][j]  +=  (re_l_t)a1 * b7;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+7]];
          wide_block[i][j]  +=  (re_l_t)a1 * b8;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+8]];
          wide_block[i][j]  +=  (re_l_t)a1 * b9;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+9]];
          wide_block[i][j]  +=  (re_l_t)a1 * b10;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+10]];
          wide_block[i][j]  +=  (re_l_t)a1 * b11;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+11]];
          wide_block[i][j]  +=  (re_l_t)a1 * b12;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+12]];
          wide_block[i][j]  +=  (re_l_t)a1 * b13;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+13]];
          wide_block[i][j]  +=  (re_l_t)a1 * b14;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+14]];
          wide_block[i][j]  +=  (re_l_t)a1 * b15;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+15]];
          wide_block[i][j]  +=  (re_l_t)a1 * b16;

          /*
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k]];
          a2  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+1]];
          a3  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+2]];
          a4  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+3]];
          a5  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+4]];
          a6  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+5]];
          a7  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+6]];
          a8  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+7]];
          wide_block[i][j]  +=  (re_l_t)a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4
            + a5 * b5 + a6 * b6 + a7 * b7 + a8 * b8;
            */
        }
      }
    }
#endif
    if (block_B->sz[j]-k > 7) {
      for (; k<block_B->sz[j]-7; k=k+8) {
        /* register const re_m_t */ b1  = block_B->val[j][k];
        /* register const re_m_t */ b2  = block_B->val[j][k+1];
        /* register const re_m_t */ b3  = block_B->val[j][k+2];
        /* register const re_m_t */ b4  = block_B->val[j][k+3];
        /* register const re_m_t */ b5  = block_B->val[j][k+4];
        /* register const re_m_t */ b6  = block_B->val[j][k+5];
        /* register const re_m_t */ b7  = block_B->val[j][k+6];
        /* register const re_m_t */ b8  = block_B->val[j][k+7];
        for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; i=i+1) {
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k]];
          wide_block[i][j]  +=  (re_l_t)a1 * b1;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+1]];
          wide_block[i][j]  +=  (re_l_t)a1 * b2;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+2]];
          wide_block[i][j]  +=  (re_l_t)a1 * b3;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+3]];
          wide_block[i][j]  +=  (re_l_t)a1 * b4;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+4]];
          wide_block[i][j]  +=  (re_l_t)a1 * b5;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+5]];
          wide_block[i][j]  +=  (re_l_t)a1 * b6;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+6]];
          wide_block[i][j]  +=  (re_l_t)a1 * b7;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+7]];
          wide_block[i][j]  +=  (re_l_t)a1 * b8;
          /*
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k]];
          a2  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+1]];
          a3  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+2]];
          a4  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+3]];
          a5  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+4]];
          a6  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+5]];
          a7  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+6]];
          a8  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+7]];
          wide_block[i][j]  +=  (re_l_t)a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4
            + a5 * b5 + a6 * b6 + a7 * b7 + a8 * b8;
            */
        }
      }
    }
    if (block_B->sz[j]-k > 3) {
      for (; k<block_B->sz[j]-3; k=k+4) {
        /* register const re_m_t */ b1  = block_B->val[j][k];
        /* register const re_m_t */ b2  = block_B->val[j][k+1];
        /* register const re_m_t */ b3  = block_B->val[j][k+2];
        /* register const re_m_t */ b4  = block_B->val[j][k+3];
        for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k]];
          wide_block[i][j]  +=  (re_l_t)a1 * b1;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+1]];
          wide_block[i][j]  +=  (re_l_t)a1 * b2;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+2]];
          wide_block[i][j]  +=  (re_l_t)a1 * b3;
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k+3]];
          wide_block[i][j]  +=  (re_l_t)a1 * b4;
        }
      }
    }
    for (; k<block_B->sz[j]; ++k) {
      /* register const re_m_t */ b1  = block_B->val[j][k];
      for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
        a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + posB[k]];
        if (a1 != 0)
          wide_block[i][j]  += (re_l_t)a1 * b1;
      }
    }
  }
#else
  bi_t i, j, k, l, ctr;
  bi_t indices[__GBLA_SIMD_BLOCK_SIZE];
  register re_m_t b1,b2,b3,b4,b5,b6,b7,b8;
  for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
    for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; i=i+2) {
      ctr = 0;
      for (k=0; k<block_B->sz[j]; ++k) {
        if (block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k]] != 0 ||
            block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k]] != 0) {
          indices[ctr]  = k;
          ++ctr;
        }
      }
      l=0;
      /*
      if (ctr>7) {
        for (l=0; l<ctr-7; l=l+8) {
          b1 = block_B->val[j][indices[l]];
          b2 = block_B->val[j][indices[l+1]];
          b3 = block_B->val[j][indices[l+2]];
          b4 = block_B->val[j][indices[l+3]];
          b5 = block_B->val[j][indices[l+4]];
          b6 = block_B->val[j][indices[l+5]];
          b7 = block_B->val[j][indices[l+6]];
          b8 = block_B->val[j][indices[l+7]];
          wide_block[i][j]  +=
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l]]] * b1 +
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+1]]] * b2 +
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+2]]] * b3 +
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+3]]] * b4 +
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+4]]] * b5 +
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+5]]] * b6 +
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+6]]] * b7 +
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+7]]] * b8;
          wide_block[i+1][j]  +=
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l]]] * b1 +
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+1]]] * b2 +
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+2]]] * b3 +
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+3]]] * b4 +
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+4]]] * b5 +
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+5]]] * b6 +
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+6]]] * b7 +
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+7]]] * b8;
        }
      }
      */
      if (ctr>1) {
        for (l=0; l<ctr-1; l=l+2) {
          b1 = block_B->val[j][indices[l]];
          b2 = block_B->val[j][indices[l+1]];
          wide_block[i][j]  +=
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l]]] * b1 +
            (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+1]]] * b2;
          wide_block[i+1][j]  +=
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l]]] * b1 +
            (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l+1]]] * b2;
        }
      }
      for (; l<ctr; ++l) {
        b1 = block_B->val[j][indices[l]];
        wide_block[i][j]  +=
          (re_l_t)block_A[i*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l]]] * b1;
        wide_block[i+1][j]  +=
          (re_l_t)block_A[(i+1)*__GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][indices[l]]] * b1;
      }
    }
  }
#endif
}
#else

#define V1 1
#define V2 0
#define V3 0
/**
 * \brief Use dense blocks from A to update sparse blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 *
 * \note B is stored by row.
 *
 * \param dense block from A block_A
 *
 * \param sparse block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_dense_sparse_rectangular(const re_t *block_A, const sbl_t *block_B,
  re_l_t **wide_block)
{
#if V1
  bi_t i, j, k;
  register re_m_t a1, a2, a3, a4, a5, a6, a7, a8;
  register re_m_t b1, b2, b3, b4, b5, b6, b7, b8;
  for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
    k = 0;
    if (block_B->sz[j]-k > 7) {
      for (; k<block_B->sz[j]-7; k=k+8) {
        b1  = block_B->val[j][k];
        b2  = block_B->val[j][k+1];
        b3  = block_B->val[j][k+2];
        b4  = block_B->val[j][k+3];
        b5  = block_B->val[j][k+4];
        b6  = block_B->val[j][k+5];
        b7  = block_B->val[j][k+6];
        b8  = block_B->val[j][k+7];
        for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; i=i+1) {
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k]];
          a2  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+1]];
          a3  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+2]];
          a4  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+3]];
          a5  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+4]];
          a6  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+5]];
          a7  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+6]];
          a8  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+7]];
          wide_block[i][j]  +=  (re_l_t)a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4
            + a5 * b5 + a6 * b6 + a7 * b7 + a8 * b8;
          /*
          a1  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k]];
          a2  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+1]];
          a3  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+2]];
          a4  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+3]];
          a5  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+4]];
          a6  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+5]];
          a7  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+6]];
          a8  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+7]];
          wide_block[i+1][j]  +=  (re_l_t)a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4
            + a5 * b5 + a6 * b6 + a7 * b7 + a8 * b8;
            */
        }
      }
    }
    if (block_B->sz[j]-k > 3) {
      for (; k<block_B->sz[j]-3; k=k+4) {
        b1  = block_B->val[j][k];
        b2  = block_B->val[j][k+1];
        b3  = block_B->val[j][k+2];
        b4  = block_B->val[j][k+3];
        for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k]];
          a2  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+1]];
          a3  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+2]];
          a4  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k+3]];
          wide_block[i][j]  +=  (re_l_t)a1 * b1 + a2 * b2 + a3 * b3 + a4 * b4;
        }
      }
    }
    for (; k<block_B->sz[j]; ++k) {
      b1  = block_B->val[j][k];
      for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
        a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + block_B->pos[j][k]];
        if (a1 != 0)
          wide_block[i][j]  +=  a1 * b1;
      }
    }
  }
#endif
#if V2
  bi_t i, j, k;
  register re_m_t a1, a2, a3, a4, a5, a6, a7, a8;
  //register re_m_t a21, a22, a23, a24, a25, a26, a27, a28;
  register re_m_t b1, b2, b3, b4, b5, b6, b7, b8;
  for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
    k = 0;
    if (block_B->sz[j]-k > 7) {
      for (; k<block_B->sz[j]-7; k=k+8) {
        b1  = block_B->val[j][k];
        b2  = block_B->val[j][k+1];
        b3  = block_B->val[j][k+2];
        b4  = block_B->val[j][k+3];
        b5  = block_B->val[j][k+4];
        b6  = block_B->val[j][k+5];
        b7  = block_B->val[j][k+6];
        b8  = block_B->val[j][k+7];
        for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
          a1 = block_A[i * __GBLA_SIMD_BLOCK_SIZE + j];

          if (a1 != 0) {
            wide_block[i][block_B->pos[j][k]]  +=  (re_l_t)a1 * b1;
            wide_block[i][block_B->pos[j][k+1]]  +=  (re_l_t)a1 * b2;
            wide_block[i][block_B->pos[j][k+2]]  +=  (re_l_t)a1 * b3;
            wide_block[i][block_B->pos[j][k+3]]  +=  (re_l_t)a1 * b4;
            wide_block[i][block_B->pos[j][k+4]]  +=  (re_l_t)a1 * b5;
            wide_block[i][block_B->pos[j][k+5]]  +=  (re_l_t)a1 * b6;
            wide_block[i][block_B->pos[j][k+6]]  +=  (re_l_t)a1 * b7;
            wide_block[i][block_B->pos[j][k+7]]  +=  (re_l_t)a1 * b8;
          }
        }
      }
    }
    if (block_B->sz[j]-k > 3) {
      for (; k<block_B->sz[j]-3; k=k+4) {
        b1  = block_B->val[j][k];
        b2  = block_B->val[j][k+1];
        b3  = block_B->val[j][k+2];
        b4  = block_B->val[j][k+3];
        for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
          a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + j];
          if (a1 != 0) {
            wide_block[i][block_B->pos[j][k]]  +=  (re_l_t)a1 * b1;
            wide_block[i][block_B->pos[j][k+1]]  +=  (re_l_t)a1 * b2;
            wide_block[i][block_B->pos[j][k+2]]  +=  (re_l_t)a1 * b3;
            wide_block[i][block_B->pos[j][k+3]]  +=  (re_l_t)a1 * b4;
          }
        }
      }
    }
    for (; k<block_B->sz[j]; ++k) {
      b1  = block_B->val[j][k];
      for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
        a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + j];
        if (a1 != 0)
          wide_block[i][block_B->pos[j][k]]  +=  (re_l_t)a1 * b1;
      }
    }
  }
#endif
#if V3
  bi_t i, j, k, ctr;
  register re_m_t a1;
  register re_m_t a2;
  register re_m_t b1, b2, b3, b4, b5, b6, b7, b8;
  for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
      k = 0;
      if (block_B->sz[j]-k > 7) {
        for (; k<block_B->sz[j]-7; k=k+8) {
          b1  = block_B->val[j][k];
          b2  = block_B->val[j][k+1];
          b3  = block_B->val[j][k+2];
          b4  = block_B->val[j][k+3];
          b5  = block_B->val[j][k+4];
          b6  = block_B->val[j][k+5];
          b7  = block_B->val[j][k+6];
          b8  = block_B->val[j][k+7];
          for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; i=i+2) {
            ctr = 0;
            if (block_A[i * __GBLA_SIMD_BLOCK_SIZE + j] != 0) {
              a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + j];
              ctr = 1;
            }
            if (block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + j] != 0) {
              a2  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + j];
              if (ctr == 1)
                ctr = 3;
              else
                ctr = 2;
            }
            if (ctr == 0)
              continue;

            if (ctr == 3) {
              wide_block[i][block_B->pos[j][k]]    +=  (re_l_t)a1 * b1;
              wide_block[i][block_B->pos[j][k+1]]  +=  (re_l_t)a1 * b2;
              wide_block[i][block_B->pos[j][k+2]]  +=  (re_l_t)a1 * b3;
              wide_block[i][block_B->pos[j][k+3]]  +=  (re_l_t)a1 * b4;
              wide_block[i][block_B->pos[j][k+4]]  +=  (re_l_t)a1 * b5;
              wide_block[i][block_B->pos[j][k+5]]  +=  (re_l_t)a1 * b6;
              wide_block[i][block_B->pos[j][k+6]]  +=  (re_l_t)a1 * b7;
              wide_block[i][block_B->pos[j][k+7]]  +=  (re_l_t)a1 * b8;

              wide_block[i+1][block_B->pos[j][k]]    +=  (re_l_t)a2 * b1;
              wide_block[i+1][block_B->pos[j][k+1]]  +=  (re_l_t)a2 * b2;
              wide_block[i+1][block_B->pos[j][k+2]]  +=  (re_l_t)a2 * b3;
              wide_block[i+1][block_B->pos[j][k+3]]  +=  (re_l_t)a2 * b4;
              wide_block[i+1][block_B->pos[j][k+4]]  +=  (re_l_t)a2 * b5;
              wide_block[i+1][block_B->pos[j][k+5]]  +=  (re_l_t)a2 * b6;
              wide_block[i+1][block_B->pos[j][k+6]]  +=  (re_l_t)a2 * b7;
              wide_block[i+1][block_B->pos[j][k+7]]  +=  (re_l_t)a2 * b8;
              continue;
            }
            if (ctr == 1) {
              wide_block[i][block_B->pos[j][k]]    +=  (re_l_t)a1 * b1;
              wide_block[i][block_B->pos[j][k+1]]  +=  (re_l_t)a1 * b2;
              wide_block[i][block_B->pos[j][k+2]]  +=  (re_l_t)a1 * b3;
              wide_block[i][block_B->pos[j][k+3]]  +=  (re_l_t)a1 * b4;
              wide_block[i][block_B->pos[j][k+4]]  +=  (re_l_t)a1 * b5;
              wide_block[i][block_B->pos[j][k+5]]  +=  (re_l_t)a1 * b6;
              wide_block[i][block_B->pos[j][k+6]]  +=  (re_l_t)a1 * b7;
              wide_block[i][block_B->pos[j][k+7]]  +=  (re_l_t)a1 * b8;
              continue;
            }
            if (ctr == 2) {
              wide_block[i+1][block_B->pos[j][k]]    +=  (re_l_t)a2 * b1;
              wide_block[i+1][block_B->pos[j][k+1]]  +=  (re_l_t)a2 * b2;
              wide_block[i+1][block_B->pos[j][k+2]]  +=  (re_l_t)a2 * b3;
              wide_block[i+1][block_B->pos[j][k+3]]  +=  (re_l_t)a2 * b4;
              wide_block[i+1][block_B->pos[j][k+4]]  +=  (re_l_t)a2 * b5;
              wide_block[i+1][block_B->pos[j][k+5]]  +=  (re_l_t)a2 * b6;
              wide_block[i+1][block_B->pos[j][k+6]]  +=  (re_l_t)a2 * b7;
              wide_block[i+1][block_B->pos[j][k+7]]  +=  (re_l_t)a2 * b8;
              continue;
            }
          }
        }
      }
      if (block_B->sz[j]-k > 3) {
        for (; k<block_B->sz[j]-3; k=k+4) {
          b1  = block_B->val[j][k];
          b2  = block_B->val[j][k+1];
          b3  = block_B->val[j][k+2];
          b4  = block_B->val[j][k+3];
          for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; i=i+2) {
            ctr = 0;
            if (block_A[i * __GBLA_SIMD_BLOCK_SIZE + j] != 0) {
              a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + j];
              ctr = 1;
            }
            if (block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + j] != 0) {
              a2  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + j];
              if (ctr == 1)
                ctr = 3;
              else
                ctr = 2;
            }
            if (ctr == 0)
              continue;

            if (ctr == 3) {
              wide_block[i][block_B->pos[j][k]]    +=  (re_l_t)a1 * b1;
              wide_block[i][block_B->pos[j][k+1]]  +=  (re_l_t)a1 * b2;
              wide_block[i][block_B->pos[j][k+2]]  +=  (re_l_t)a1 * b3;
              wide_block[i][block_B->pos[j][k+3]]  +=  (re_l_t)a1 * b4;

              wide_block[i+1][block_B->pos[j][k]]    +=  (re_l_t)a2 * b1;
              wide_block[i+1][block_B->pos[j][k+1]]  +=  (re_l_t)a2 * b2;
              wide_block[i+1][block_B->pos[j][k+2]]  +=  (re_l_t)a2 * b3;
              wide_block[i+1][block_B->pos[j][k+3]]  +=  (re_l_t)a2 * b4;
              continue;
            }
            if (ctr == 1) {
              wide_block[i][block_B->pos[j][k]]    +=  (re_l_t)a1 * b1;
              wide_block[i][block_B->pos[j][k+1]]  +=  (re_l_t)a1 * b2;
              wide_block[i][block_B->pos[j][k+2]]  +=  (re_l_t)a1 * b3;
              wide_block[i][block_B->pos[j][k+3]]  +=  (re_l_t)a1 * b4;
              continue;
            }
            if (ctr == 2) {
              wide_block[i+1][block_B->pos[j][k]]    +=  (re_l_t)a2 * b1;
              wide_block[i+1][block_B->pos[j][k+1]]  +=  (re_l_t)a2 * b2;
              wide_block[i+1][block_B->pos[j][k+2]]  +=  (re_l_t)a2 * b3;
              wide_block[i+1][block_B->pos[j][k+3]]  +=  (re_l_t)a2 * b4;
              continue;
            }
          }
        }
      }
      for (; k<block_B->sz[j]; ++k) {
        b1  = block_B->val[j][k];
        for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; i=i+2) {
          ctr = 0;
          if (block_A[i * __GBLA_SIMD_BLOCK_SIZE + j] != 0) {
            a1  = block_A[i * __GBLA_SIMD_BLOCK_SIZE + j];
            ctr = 1;
          }
          if (block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + j] != 0) {
            a2  = block_A[(i+1) * __GBLA_SIMD_BLOCK_SIZE + j];
            if (ctr == 1)
              ctr = 3;
            else
              ctr = 2;
          }
          if (ctr == 0)
            continue;
        if (ctr == 3) {
          wide_block[i][block_B->pos[j][k]]   +=  (re_l_t)a1 * b1;
          wide_block[i+1][block_B->pos[j][k]] +=  (re_l_t)a2 * b1;
          continue;
        }
        if (ctr == 1) {
          wide_block[i][block_B->pos[j][k]]   +=  (re_l_t)a1 * b1;
          continue;
        }
        if (ctr == 2) {
          wide_block[i+1][block_B->pos[j][k]] +=  (re_l_t)a2 * b1;
          continue;
        }
      }
    }
}
#endif
}
#endif


#if __GBLA_COLUMN_B
/**
* \brief Use compressed intermediate sparse blocks from A to update sparse blocks in B.
* Delay modulus operations by storing results in wide blocks.
*
* \note B is stored by column.
*
* \param intermediate block from A block_A
*
 * \param sparse block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_intermediate_sparse_sparse_rectangular(const ibl_t *block_A, const sbl_t *block_B,
  re_l_t **wide_block)
{
#if 0
  register re_m_t a, b;
  register ci_t pos;
  bi_t posA[__GBLA_SIMD_BLOCK_SIZE];
  bi_t posB[__GBLA_SIMD_BLOCK_SIZE];
  bi_t i, j, k, l, m, ctr;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
STEP_0:
      //printf("%u || %u\n",i,j);
      ctr = 0;
      l   = 0;
      k   = 0;
      if (block_A->sz[i] > 0 && block_B->sz[j] > 0) {
STEP_1:
        //printf(">> %u || %u\n", k, l);
        if (block_A->pos[i][k] == block_B->pos[j][l]) {
          posA[ctr] = k;
          posB[ctr] = l;
          ++ctr;
          ++k;
          if (k==block_A->sz[i])
            goto STEP_0;
          ++l;
          if (l==block_B->sz[j]) {
            ++j;
            goto STEP_0;
          }
          goto STEP_1;
        } else {
          if (block_A->pos[i][k] < block_B->pos[j][l]) {
            ++k;
            while (k < block_A->sz[i] && block_A->pos[i][k] < block_B->pos[j][l])
              ++k;
            if (k == block_A->sz[i]) {
              for (m=0; m<ctr; ++m)
                wide_block[i][j]  +=  (re_l_t)block_A->row[i][posA[ctr]] * block_B->val[j][posB[ctr]];
              ++j;
              goto STEP_0;
            } else {
              goto STEP_1;
            }
          } else {
            ++l;
            while (l < block_B->sz[j] && block_A->pos[i][k] > block_B->pos[j][l])
              ++l;
            if (l == block_B->sz[j]) {
              for (m=0; m<ctr; ++m)
                wide_block[i][j]  +=  (re_l_t)block_A->row[i][posA[ctr]] * block_B->val[j][posB[ctr]];
              ++j;
              goto STEP_0;
            } else {
              goto STEP_1;
            }
          }
        }
      }
    }
  }
#else
  bi_t i, j, k, l, ctr = 0;
  register re_l_t wb;
  re_m_t a[4] = {0};
  re_m_t b[4] = {0};
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    j = 0;
STEP_0:
    for (; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
      wb  = wide_block[i][j];
      l   = 0;
      for (k=0; k<block_B->sz[j]; ++k) {
        //printf("%u | %u | %u | %u\n", i,j,k,l);
        while (l<block_A->sz[i] && block_A->pos[i][l] < block_B->pos[j][k])
          ++l;
        if (l == block_A->sz[i]) {
          wide_block[i][j]  = wb;
          ++j;
          goto STEP_0;
        }
        if (block_A->pos[i][l]  ==  block_B->pos[j][k]) {
          a[ctr]  = block_A->row[i][l];
          b[ctr]  = block_B->val[j][k];
          ++ctr;
          if (ctr == 4) {
            wb  +=  (re_l_t)a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
            ctr = 0;
          }
        }
      }
      // multiply some leftovers
      for (l=0; l<ctr; ++l)
        wb  +=  (re_l_t)a[l]*b[l];
      
      wide_block[i][j]  = wb;
    }
  }
#endif
}
#else
/**
 * \brief Use intermediate blocks from A to update sparse blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 *
 * \param intermediate block from A block_A
 *
 * \param sparse block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_intermediate_sparse_sparse_rectangular(const ibl_t *block_A, const sbl_t *block_B,
  re_l_t **wide_block)
{
  bi_t i, j, k;
  register re_m_t  a1;
  register re_m_t b1, b2, b3, b4, b5, b6, b7, b8;
  bi_t pos1;
  re_t *row_B1;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    j = 0;
    for (;j<block_A->sz[i]; ++j) {
      a1 = block_A->row[i][j];
      pos1 = block_A->pos[i][j];
      k = 0;
      row_B1 = block_B->val[pos1];
      if (block_B->sz[pos1] > 7) {
        for (; k<block_B->sz[pos1]-7; k=k+8) {
          b1  = (re_m_t)row_B1[k];
          b2  = (re_m_t)row_B1[k+1];
          b3  = (re_m_t)row_B1[k+2];
          b4  = (re_m_t)row_B1[k+3];
          b5  = (re_m_t)row_B1[k+4];
          b6  = (re_m_t)row_B1[k+5];
          b7  = (re_m_t)row_B1[k+6];
          b8  = (re_m_t)row_B1[k+7];
          wide_block[i][block_B->pos[pos1][k]]  +=
            (re_l_t)a1 * b1;
          wide_block[i][block_B->pos[pos1][k+1]]  +=
            (re_l_t)a1 * b2;
          wide_block[i][block_B->pos[pos1][k+2]]  +=
            (re_l_t)a1 * b3;
          wide_block[i][block_B->pos[pos1][k+3]]  +=
            (re_l_t)a1 * b4;
          wide_block[i][block_B->pos[pos1][k+4]]  +=
            (re_l_t)a1 * b5;
          wide_block[i][block_B->pos[pos1][k+5]]  +=
            (re_l_t)a1 * b6;
          wide_block[i][block_B->pos[pos1][k+6]]  +=
            (re_l_t)a1 * b7;
          wide_block[i][block_B->pos[pos1][k+7]]  +=
            (re_l_t)a1 * b8;
        }
      }
      if (block_B->sz[pos1] > 3) {
        for (; k<block_B->sz[pos1]-3; k=k+4) {
          b1  = (re_m_t)row_B1[k];
          b2  = (re_m_t)row_B1[k+1];
          b3  = (re_m_t)row_B1[k+2];
          b4  = (re_m_t)row_B1[k+3];
          wide_block[i][block_B->pos[pos1][k]]  +=
            (re_l_t)a1 * b1;
          wide_block[i][block_B->pos[pos1][k+1]]  +=
            (re_l_t)a1 * b2;
          wide_block[i][block_B->pos[pos1][k+2]]  +=
            (re_l_t)a1 * b3;
          wide_block[i][block_B->pos[pos1][k+3]]  +=
            (re_l_t)a1 * b4;
        }
      }
      for (; k<block_B->sz[pos1]; ++k) {
        wide_block[i][block_B->pos[pos1][k]]  +=
          (re_l_t)a1 * (re_m_t)row_B1[k];
      }
    }
  }
}
#endif
#if __GBLA_COLUMN_B
/**
 * \brief Use compressed sparse blocks from A to update sparse blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 *
 * \note B is stored by column.
 *
 * \param sparse block from A block_A
 *
 * \param sparse block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_sparse_sparse_rectangular(const sbl_t *block_A, const sbl_t *block_B,
  re_l_t **wide_block)
{
#if 0
  register re_m_t a, b;
  register ci_t pos;
  bi_t posA[__GBLA_SIMD_BLOCK_SIZE];
  bi_t posB[__GBLA_SIMD_BLOCK_SIZE];
  bi_t i, j, k, l, m, ctr;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
STEP_0:
      //printf("%u || %u\n",i,j);
      ctr = 0;
      l   = 0;
      k   = 0;
      if (block_A->sz[i] > 0 && block_B->sz[j] > 0) {
STEP_1:
        //printf(">> %u || %u\n", k, l);
        if (block_A->pos[i][k] == block_B->pos[j][l]) {
          posA[ctr] = k;
          posB[ctr] = l;
          ++ctr;
          ++k;
          if (k==block_A->sz[i])
            goto STEP_0;
          ++l;
          if (l==block_B->sz[j]) {
            ++j;
            goto STEP_0;
          }
          goto STEP_1;
        } else {
          if (block_A->pos[i][k] < block_B->pos[j][l]) {
            ++k;
            while (k < block_A->sz[i] && block_A->pos[i][k] < block_B->pos[j][l])
              ++k;
            if (k == block_A->sz[i]) {
              for (m=0; m<ctr; ++m)
                wide_block[i][j]  +=  (re_l_t)block_A->val[i][posA[ctr]] * block_B->val[j][posB[ctr]];
              ++j;
              goto STEP_0;
            } else {
              goto STEP_1;
            }
          } else {
            ++l;
            while (l < block_B->sz[j] && block_A->pos[i][k] > block_B->pos[j][l])
              ++l;
            if (l == block_B->sz[j]) {
              for (m=0; m<ctr; ++m)
                wide_block[i][j]  +=  (re_l_t)block_A->val[i][posA[ctr]] * block_B->val[j][posB[ctr]];
              ++j;
              goto STEP_0;
            } else {
              goto STEP_1;
            }
          }
        }
      }
    }
  }
#else
  bi_t i, j, k, l, ctr = 0;
  register re_l_t wb;
  re_m_t a[4] = {0};
  re_m_t b[4] = {0};
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    j = 0;
STEP_0:
    for (; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
      wb  = wide_block[i][j];
      l   = 0;
      for (k=0; k<block_B->sz[j]; ++k) {
        //printf("%u | %u | %u | %u\n", i,j,k,l);
        while (l<block_A->sz[i] && block_A->pos[i][l] < block_B->pos[j][k])
          ++l;
        if (l == block_A->sz[i]) {
          wide_block[i][j]  = wb;
          ++j;
          goto STEP_0;
        }
        if (block_A->pos[i][l]  ==  block_B->pos[j][k]) {
          a[ctr]  = block_A->val[i][l];
          b[ctr]  = block_B->val[j][k];
          ++ctr;
          if (ctr == 4) {
            wb  +=  (re_l_t)a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
            ctr = 0;
          }
        }
      }
      // multiply some leftovers
      for (l=0; l<ctr; ++l)
        wb  +=  (re_l_t)a[l]*b[l];
      
      wide_block[i][j]  = wb;
    }
  }
#endif
}
#else
/**
 * \brief Use compressed sparse blocks from A to update sparse blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 *
 * \param sparse block from A block_A
 *
 * \param sparse block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_sparse_sparse_rectangular(const sbl_t *block_A, const sbl_t *block_B,
  re_l_t **wide_block)
{
  bi_t i, j, k;
  register re_m_t  a1;
  register re_m_t b1, b2, b3, b4, b5, b6, b7, b8;
  bi_t pos1;
  re_t *val_B1;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    j = 0;
    for (;j<block_A->sz[i]; ++j) {
      a1 = block_A->val[i][j];
      pos1 = block_A->pos[i][j];
      k = 0;
      val_B1 = block_B->val[pos1];
      if (block_B->sz[pos1] > 7) {
        for (; k<block_B->sz[pos1]-7; k=k+8) {
          b1  = (re_m_t)val_B1[k];
          b2  = (re_m_t)val_B1[k+1];
          b3  = (re_m_t)val_B1[k+2];
          b4  = (re_m_t)val_B1[k+3];
          b5  = (re_m_t)val_B1[k+4];
          b6  = (re_m_t)val_B1[k+5];
          b7  = (re_m_t)val_B1[k+6];
          b8  = (re_m_t)val_B1[k+7];
          wide_block[i][block_B->pos[pos1][k]]  +=
            (re_l_t)a1 * b1;
          wide_block[i][block_B->pos[pos1][k+1]]  +=
            (re_l_t)a1 * b2;
          wide_block[i][block_B->pos[pos1][k+2]]  +=
            (re_l_t)a1 * b3;
          wide_block[i][block_B->pos[pos1][k+3]]  +=
            (re_l_t)a1 * b4;
          wide_block[i][block_B->pos[pos1][k+4]]  +=
            (re_l_t)a1 * b5;
          wide_block[i][block_B->pos[pos1][k+5]]  +=
            (re_l_t)a1 * b6;
          wide_block[i][block_B->pos[pos1][k+6]]  +=
            (re_l_t)a1 * b7;
          wide_block[i][block_B->pos[pos1][k+7]]  +=
            (re_l_t)a1 * b8;
        }
      }
      if (block_B->sz[pos1] > 3) {
        for (; k<block_B->sz[pos1]-3; k=k+4) {
          b1  = (re_m_t)val_B1[k];
          b2  = (re_m_t)val_B1[k+1];
          b3  = (re_m_t)val_B1[k+2];
          b4  = (re_m_t)val_B1[k+3];
          wide_block[i][block_B->pos[pos1][k]]  +=
            (re_l_t)a1 * b1;
          wide_block[i][block_B->pos[pos1][k+1]]  +=
            (re_l_t)a1 * b2;
          wide_block[i][block_B->pos[pos1][k+2]]  +=
            (re_l_t)a1 * b3;
          wide_block[i][block_B->pos[pos1][k+3]]  +=
            (re_l_t)a1 * b4;
        }
      }
      for (; k<block_B->sz[pos1]; ++k) {
        wide_block[i][block_B->pos[pos1][k]]  +=
          (re_l_t)a1 * (re_m_t)val_B1[k];
      }
    }
  }
}
#endif


/**
 * \brief Use compressed dense blocks from A to update dense blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 *
 * \param sparse block from A block_A
 *
 * \param dense block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_sparse_dense_rectangular(const sbl_t *block_A, const re_t *block_B,
  re_l_t **wide_block)
{
#ifdef NOSSE
  bi_t i, j, k;
  register re_t a, b, c, d, e, f, g, h;
  register bi_t pa, pb, pc, pd, pe, pf, pg, ph;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    j = 0;
    if (block_A->sz[i] > 7) {
      for (; j<block_A->sz[i]-7; j = j+8) {
        a = block_A->val[i][j];
        b = block_A->val[i][j+1];
        c = block_A->val[i][j+2];
        d = block_A->val[i][j+3];
        e = block_A->val[i][j+4];
        f = block_A->val[i][j+5];
        g = block_A->val[i][j+6];
        h = block_A->val[i][j+7];
        pa  = block_A->pos[i][j] * __GBLA_SIMD_BLOCK_SIZE;
        pb  = block_A->pos[i][j+1] * __GBLA_SIMD_BLOCK_SIZE;
        pc  = block_A->pos[i][j+2] * __GBLA_SIMD_BLOCK_SIZE;
        pd  = block_A->pos[i][j+3] * __GBLA_SIMD_BLOCK_SIZE;
        pe  = block_A->pos[i][j+4] * __GBLA_SIMD_BLOCK_SIZE;
        pf  = block_A->pos[i][j+5] * __GBLA_SIMD_BLOCK_SIZE;
        pg  = block_A->pos[i][j+6] * __GBLA_SIMD_BLOCK_SIZE;
        ph  = block_A->pos[i][j+7] * __GBLA_SIMD_BLOCK_SIZE;
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k = k+4) {
          wide_block[i][k]  +=
            (a * (re_l_t)block_B[pa + k] +
             b * (re_m_t)block_B[pb + k] +
             c * (re_m_t)block_B[pc + k] +
             d * (re_m_t)block_B[pd + k] +
             e * (re_m_t)block_B[pe + k] +
             f * (re_m_t)block_B[pf + k] +
             g * (re_m_t)block_B[pg + k] +
             h * (re_m_t)block_B[ph + k]);
          wide_block[i][k+1]  +=
            (a * (re_l_t)block_B[pa + k+1] +
             b * (re_m_t)block_B[pb + k+1] +
             c * (re_m_t)block_B[pc + k+1] +
             d * (re_m_t)block_B[pd + k+1] +
             e * (re_m_t)block_B[pe + k+1] +
             f * (re_m_t)block_B[pf + k+1] +
             g * (re_m_t)block_B[pg + k+1] +
             h * (re_m_t)block_B[ph + k+1]);
          wide_block[i][k+2]  +=
            (a * (re_l_t)block_B[pa + k+2] +
             b * (re_m_t)block_B[pb + k+2] +
             c * (re_m_t)block_B[pc + k+2] +
             d * (re_m_t)block_B[pd + k+2] +
             e * (re_m_t)block_B[pe + k+2] +
             f * (re_m_t)block_B[pf + k+2] +
             g * (re_m_t)block_B[pg + k+2] +
             h * (re_m_t)block_B[ph + k+2]);
          wide_block[i][k+3]  +=
            (a * (re_l_t)block_B[pa + k+3] +
             b * (re_m_t)block_B[pb + k+3] +
             c * (re_m_t)block_B[pc + k+3] +
             d * (re_m_t)block_B[pd + k+3] +
             e * (re_m_t)block_B[pe + k+3] +
             f * (re_m_t)block_B[pf + k+3] +
             g * (re_m_t)block_B[pg + k+3] +
             h * (re_m_t)block_B[ph + k+3]);
        }
      }
    }
    if (block_A->sz[i]-j > 3) {
      a = block_A->val[i][j];
      b = block_A->val[i][j+1];
      c = block_A->val[i][j+2];
      d = block_A->val[i][j+3];
      pa  = block_A->pos[i][j] * __GBLA_SIMD_BLOCK_SIZE;
      pb  = block_A->pos[i][j+1] * __GBLA_SIMD_BLOCK_SIZE;
      pc  = block_A->pos[i][j+2] * __GBLA_SIMD_BLOCK_SIZE;
      pd  = block_A->pos[i][j+3] * __GBLA_SIMD_BLOCK_SIZE;
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k = k+4) {
        wide_block[i][k]  +=
          (a * (re_l_t)block_B[pa + k] +
           b * (re_m_t)block_B[pb + k] +
           c * (re_m_t)block_B[pc + k] +
           d * (re_m_t)block_B[pd + k]);
          wide_block[i][k+1]  +=
            (a * (re_l_t)block_B[pa + k+1] +
             b * (re_m_t)block_B[pb + k+1] +
             c * (re_m_t)block_B[pc + k+1] +
             d * (re_m_t)block_B[pd + k+1]);
          wide_block[i][k+2]  +=
            (a * (re_l_t)block_B[pa + k+2] +
             b * (re_m_t)block_B[pb + k+2] +
             c * (re_m_t)block_B[pc + k+2] +
             d * (re_m_t)block_B[pd + k+2]);
          wide_block[i][k+3]  +=
            (a * (re_l_t)block_B[pa + k+3] +
             b * (re_m_t)block_B[pb + k+3] +
             c * (re_m_t)block_B[pc + k+3] +
             d * (re_m_t)block_B[pd + k+3]);
      }
      j = j+4;
    }
    for (;j<block_A->sz[i]; ++j) {
      a   = block_A->val[i][j];
      pa  = block_A->pos[i][j] * __GBLA_SIMD_BLOCK_SIZE;
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k = k+4) {
        wide_block[i][k]  +=
          a * (re_m_t)block_B[pa + k];
        wide_block[i][k+1]  +=
          a * (re_m_t)block_B[pa + k+1];
        wide_block[i][k+2]  +=
          a * (re_m_t)block_B[pa + k+2];
        wide_block[i][k+3]  +=
          a * (re_m_t)block_B[pa + k+3];
      }
    }
  }
#endif
#ifdef NOSSETEST
  bi_t i, j, k;
  register re_t a, b;
  register bi_t pa, pb;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {

    j = block_A->sz[i] % 2;
    j = j & 1 ? 1 : 0;
    while (j < block_A->sz[i]) {
      a = block_A->val[i][j];
      b = block_A->val[i][j+1];
      pa  = block_A->pos[i][j] * __GBLA_SIMD_BLOCK_SIZE;
      pb  = block_A->pos[i][j+1] * __GBLA_SIMD_BLOCK_SIZE;
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k = k+4) {
        wide_block[i][k]  +=
          (a * (re_l_t)block_B[pa + k] +
           b * (re_m_t)block_B[pb + k]);
        wide_block[i][k+1]  +=
          (a * (re_l_t)block_B[pa + k+1] +
           b * (re_m_t)block_B[pb + k+1]);
        wide_block[i][k+2]  +=
          (a * (re_l_t)block_B[pa + k+2] +
           b * (re_m_t)block_B[pb + k+2]);
        wide_block[i][k+3]  +=
          (a * (re_l_t)block_B[pa + k+3] +
           b * (re_m_t)block_B[pb + k+3]);
      }
      j +=  2;
    }

    printf("block_A->sz = %u\n",block_A->sz);
    // extra step we have to due due to above loop unrolling
    a = block_A->val[i][0];
    pa  = block_A->pos[i][0] * __GBLA_SIMD_BLOCK_SIZE;
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k = k+4) {
      wide_block[i][k]  +=
        (a * (re_l_t)block_B[pa + k]);
      wide_block[i][k+1]  +=
        (a * (re_l_t)block_B[pa + k+1]);
      wide_block[i][k+2]  +=
        (a * (re_l_t)block_B[pa + k+2]);
      wide_block[i][k+3]  +=
        (a * (re_l_t)block_B[pa + k+3]);
    }
  }
#endif
#ifdef NOSSE2
  bi_t i, j, k;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    j = 0;
    const bi_t bs = __GBLA_SIMD_BLOCK_SIZE;
    const re_t *bAv = block_A->val[i];
    const bi_t *bAp = block_A->pos[i];
    if (block_A->sz[i] > 7) {
      for (; j<block_A->sz[i]-7; j = j+8) {
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k = k+4) {
          wide_block[i][k]  +=
            (bAv[j] * (re_l_t)block_B[bAp[j]*bs + k] +
             bAv[j+1] * (re_m_t)block_B[bAp[j+1]*bs + k] +
             bAv[j+2] * (re_m_t)block_B[bAp[j+2]*bs + k] +
             bAv[j+3] * (re_m_t)block_B[bAp[j+3]*bs + k] +
             bAv[j+4] * (re_m_t)block_B[bAp[j+4]*bs + k] +
             bAv[j+5] * (re_m_t)block_B[bAp[j+5]*bs + k] +
             bAv[j+6] * (re_m_t)block_B[bAp[j+6]*bs + k] +
             bAv[j+7] * (re_m_t)block_B[bAp[j+7]*bs + k]);
          wide_block[i][k+1]  +=
            (bAv[j] * (re_l_t)block_B[bAp[j]*bs + k+1] +
             bAv[j+1] * (re_m_t)block_B[bAp[j+1]*bs + k+1] +
             bAv[j+2] * (re_m_t)block_B[bAp[j+2]*bs + k+1] +
             bAv[j+3] * (re_m_t)block_B[bAp[j+3]*bs + k+1] +
             bAv[j+4] * (re_m_t)block_B[bAp[j+4]*bs + k+1] +
             bAv[j+5] * (re_m_t)block_B[bAp[j+5]*bs + k+1] +
             bAv[j+6] * (re_m_t)block_B[bAp[j+6]*bs + k+1] +
             bAv[j+7] * (re_m_t)block_B[bAp[j+7]*bs + k+1]);
          wide_block[i][k+2]  +=
            (bAv[j] * (re_l_t)block_B[bAp[j]*bs + k+2] +
             bAv[j+1] * (re_m_t)block_B[bAp[j+1]*bs + k+2] +
             bAv[j+2] * (re_m_t)block_B[bAp[j+2]*bs + k+2] +
             bAv[j+3] * (re_m_t)block_B[bAp[j+3]*bs + k+2] +
             bAv[j+4] * (re_m_t)block_B[bAp[j+4]*bs + k+2] +
             bAv[j+5] * (re_m_t)block_B[bAp[j+5]*bs + k+2] +
             bAv[j+6] * (re_m_t)block_B[bAp[j+6]*bs + k+2] +
             bAv[j+7] * (re_m_t)block_B[bAp[j+7]*bs + k+2]);
          wide_block[i][k+3]  +=
            (bAv[j] * (re_l_t)block_B[bAp[j]*bs + k+3] +
             bAv[j+1] * (re_m_t)block_B[bAp[j+1]*bs + k+3] +
             bAv[j+2] * (re_m_t)block_B[bAp[j+2]*bs + k+3] +
             bAv[j+3] * (re_m_t)block_B[bAp[j+3]*bs + k+3] +
             bAv[j+4] * (re_m_t)block_B[bAp[j+4]*bs + k+3] +
             bAv[j+5] * (re_m_t)block_B[bAp[j+5]*bs + k+3] +
             bAv[j+6] * (re_m_t)block_B[bAp[j+6]*bs + k+3] +
             bAv[j+7] * (re_m_t)block_B[bAp[j+7]*bs + k+3]);
        }
      }
    }
    if (block_A->sz[i]-j > 3) {
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k = k+4) {
        wide_block[i][k]  +=
          (bAv[j] * (re_l_t)block_B[bAp[j]*bs + k] +
            bAv[j+1] * (re_m_t)block_B[bAp[j+1]*bs + k] +
            bAv[j+2] * (re_m_t)block_B[bAp[j+2]*bs + k] +
            bAv[j+3] * (re_m_t)block_B[bAp[j+3]*bs + k]);
        wide_block[i][k+1]  +=
          (bAv[j] * (re_l_t)block_B[bAp[j]*bs + k+1] +
            bAv[j+1] * (re_m_t)block_B[bAp[j+1]*bs + k+1] +
            bAv[j+2] * (re_m_t)block_B[bAp[j+2]*bs + k+1] +
            bAv[j+3] * (re_m_t)block_B[bAp[j+3]*bs + k+1]);
        wide_block[i][k+2]  +=
          (bAv[j] * (re_l_t)block_B[bAp[j]*bs + k+2] +
            bAv[j+1] * (re_m_t)block_B[bAp[j+1]*bs + k+2] +
            bAv[j+2] * (re_m_t)block_B[bAp[j+2]*bs + k+2] +
            bAv[j+3] * (re_m_t)block_B[bAp[j+3]*bs + k+2]);
        wide_block[i][k+3]  +=
          (bAv[j] * (re_l_t)block_B[bAp[j]*bs + k+3] +
            bAv[j+1] * (re_m_t)block_B[bAp[j+1]*bs + k+3] +
            bAv[j+2] * (re_m_t)block_B[bAp[j+2]*bs + k+3] +
            bAv[j+3] * (re_m_t)block_B[bAp[j+3]*bs + k+3] );
      }
      j = j+4;
    }
    for (;j<block_A->sz[i]; ++j) {
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k = k+4) {
        wide_block[i][k]  +=
          (bAv[j] * (re_m_t)block_B[bAp[j]*bs + k]);
        wide_block[i][k+1]  +=
          (bAv[j] * (re_m_t)block_B[bAp[j]*bs + k+1]);
        wide_block[i][k+2]  +=
          (bAv[j] * (re_m_t)block_B[bAp[j]*bs + k+2]);
        wide_block[i][k+3]  +=
          (bAv[j] * (re_m_t)block_B[bAp[j]*bs + k+3]);
      }
    }
  }
#endif
#ifdef AVX
  bi_t i, j, k;
  register __m256i a1, a2, a3;
  register __m256i b1, b2, b3, b4;
  register __m256i b5, b6, b7, b8;
  register __m256i b9, b10, b11, b12;
  register __m256i c1, c2, c3, c4;
  bi_t pa1, pa2, pa3;
  __m256i c[__GBLA_SIMD_BLOCK_SIZE][__GBLA_SIMD_BLOCK_SIZE/4];

  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k=k+4) {
      //printf("(%u,%u)\n", i, k/4);
      c[i][k/4] = _mm256_setr_epi64x(
          wide_block[i][k],
          wide_block[i][k+1],
          wide_block[i][k+2],
          wide_block[i][k+3]);
    }
    //free(wide_block[i]);
  }
  //free(wide_block);
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k=k+16) {
      // load B block
      b1  = _mm256_setr_epi64x(
          block_B[pa1 + k],
          block_B[pa1 + k+1],
          block_B[pa1 + k+2],
          block_B[pa1 + k+3]
          );
      b2  = _mm256_setr_epi64x(
          block_B[pa1 + k+4],
          block_B[pa1 + k+5],
          block_B[pa1 + k+6],
          block_B[pa1 + k+7]
          );
      b3  = _mm256_setr_epi64x(
          block_B[pa1 + k+8],
          block_B[pa1 + k+9],
          block_B[pa1 + k+10],
          block_B[pa1 + k+11]
          );
      b4  = _mm256_setr_epi64x(
          block_B[pa1 + k+12],
          block_B[pa1 + k+13],
          block_B[pa1 + k+14],
          block_B[pa1 + k+15]
          );
      /*
      if (block_A->sz[i] > 3) {
        for (j=0; j<block_A->sz[i]-2; j=j+3) {
          // load A block
          a1  = _mm256_set1_epi64x(block_A->val[i][j]);
          pa1 = block_A->pos[i][j] * __GBLA_SIMD_BLOCK_SIZE;
          a2  = _mm256_set1_epi64x(block_A->val[i][j+1]);
          pa2 = block_A->pos[i][j+1] * __GBLA_SIMD_BLOCK_SIZE;
          a3  = _mm256_set1_epi64x(block_A->val[i][j+2]);
          pa3 = block_A->pos[i][j+2] * __GBLA_SIMD_BLOCK_SIZE;


          // load B block
          b1  = _mm256_setr_epi64x(
              block_B[pa1 + k],
              block_B[pa1 + k+1],
              block_B[pa1 + k+2],
              block_B[pa1 + k+3]
              );
          b2  = _mm256_setr_epi64x(
              block_B[pa1 + k+4],
              block_B[pa1 + k+5],
              block_B[pa1 + k+6],
              block_B[pa1 + k+7]
              );
          b3  = _mm256_setr_epi64x(
              block_B[pa1 + k+8],
              block_B[pa1 + k+9],
              block_B[pa1 + k+10],
              block_B[pa1 + k+11]
              );
          b4  = _mm256_setr_epi64x(
              block_B[pa1 + k+12],
              block_B[pa1 + k+13],
              block_B[pa1 + k+14],
              block_B[pa1 + k+15]
              );
          b5  = _mm256_setr_epi64x(
              block_B[pa2 + k],
              block_B[pa2 + k+1],
              block_B[pa2 + k+2],
              block_B[pa2 + k+3]
              );
          b6  = _mm256_setr_epi64x(
              block_B[pa2 + k+4],
              block_B[pa2 + k+5],
              block_B[pa2 + k+6],
              block_B[pa2 + k+7]
              );
          b7  = _mm256_setr_epi64x(
              block_B[pa2 + k+8],
              block_B[pa2 + k+9],
              block_B[pa2 + k+10],
              block_B[pa2 + k+11]
              );
          b8  = _mm256_setr_epi64x(
              block_B[pa2 + k+12],
              block_B[pa2 + k+13],
              block_B[pa2 + k+14],
              block_B[pa2 + k+15]
              );
          b9  = _mm256_setr_epi64x(
              block_B[pa3 + k],
              block_B[pa3 + k+1],
              block_B[pa3 + k+2],
              block_B[pa3 + k+3]
              );
          b10 = _mm256_setr_epi64x(
              block_B[pa3 + k+4],
              block_B[pa3 + k+5],
              block_B[pa3 + k+6],
              block_B[pa3 + k+7]
              );
          b11 = _mm256_setr_epi64x(
              block_B[pa3 + k+8],
              block_B[pa3 + k+9],
              block_B[pa3 + k+10],
              block_B[pa3 + k+11]
              );
          b12 = _mm256_setr_epi64x(
              block_B[pa3 + k+12],
              block_B[pa3 + k+13],
              block_B[pa3 + k+14],
              block_B[pa3 + k+15]
              );
          c1 = _mm256_add_epi64(_mm256_add_epi64(c1, _mm256_mullo_epi32(a3,b9)), _mm256_add_epi64(_mm256_mullo_epi32(a1,b1), _mm256_mullo_epi32(a2,b5)));
          c2 = _mm256_add_epi64(_mm256_add_epi64(c2, _mm256_mullo_epi32(a3,b10)), _mm256_add_epi64(_mm256_mullo_epi32(a1,b2), _mm256_mullo_epi32(a2,b6)));
          c3 = _mm256_add_epi64(_mm256_add_epi64(c3, _mm256_mullo_epi32(a3,b11)), _mm256_add_epi64(_mm256_mullo_epi32(a1,b3), _mm256_mullo_epi32(a2,b7)));
          c4 = _mm256_add_epi64(_mm256_add_epi64(c4, _mm256_mullo_epi32(a3,b12)), _mm256_add_epi64(_mm256_mullo_epi32(a1,b4), _mm256_mullo_epi32(a2,b8)));
        }
      }
      */
    j = 0;

    for (j; j<block_A->sz[i]; ++j) {
      // load A block
      a1  = _mm256_set1_epi64x(block_A->val[i][j]);
      pa1 = block_A->pos[i][j] * __GBLA_SIMD_BLOCK_SIZE;
        
        

        //printf("--(%u,%u)--\n",i,k/4);
        c[i][k/4] = _mm256_add_epi64(c[i][k/4], _mm256_mul_epu32(a1,b1));
        c[i][k/4+1] = _mm256_add_epi64(c[i][k/4+1], _mm256_mul_epu32(a1,b2));
        c[i][k/4+2] = _mm256_add_epi64(c[i][k/4+2], _mm256_mul_epu32(a1,b3));
        c[i][k/4+3] = _mm256_add_epi64(c[i][k/4+3], _mm256_mul_epu32(a1,b4));
      }
      /*
      _mm256_storeu_si256((re_l_t *)wide_block[i] + k, c[i][k/4]);
      _mm256_storeu_si256((re_l_t *)wide_block[i] + k + 4, c[i][k/4+1]);
      _mm256_storeu_si256((re_l_t *)wide_block[i] + k + 8, c[i][k/4+2]);
      _mm256_storeu_si256((re_l_t *)wide_block[i] + k + 12, c[i][k/4+3]);
      */
    }
  }
  //init_wide_blocks(&wide_block);
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k=k+4) {
      _mm256_storeu_si256((re_l_t *)wide_block[i] + k, c[i][k/4]);
    }
  }
  /*
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
      printf("%lu ",wide_block[i][k]);
    }
    printf("\n");
  }
  printf("\n");
  */
#endif
}

/**
 * \brief Use compressed dense blocks from A to update dense blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 *
 * \param compressed dense block from A block_A
 *
 * \param dense block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_dense_rectangular(const re_t *block_A, const re_t *block_B,
  re_l_t **wide_block)
{
  bi_t i, j, k;
  register re_m_t a, b, c, d;
#if 1
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; j=j+4) {
      a = block_A[i*__GBLA_SIMD_BLOCK_SIZE+j];
      b = block_A[i*__GBLA_SIMD_BLOCK_SIZE+j+1];
      c = block_A[i*__GBLA_SIMD_BLOCK_SIZE+j+2];
      d = block_A[i*__GBLA_SIMD_BLOCK_SIZE+j+3];
      if (a != 0 || b != 0) {
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
          wide_block[i][k]  +=
            (a * (re_l_t)block_B[j*__GBLA_SIMD_BLOCK_SIZE + k] +
            b * (re_l_t)block_B[(j+1)*__GBLA_SIMD_BLOCK_SIZE + k] +
            c * (re_l_t)block_B[(j+2)*__GBLA_SIMD_BLOCK_SIZE + k] +
            d * (re_l_t)block_B[(j+3)*__GBLA_SIMD_BLOCK_SIZE + k]);
        }
      }
    }
  }
#else
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (j=0; j<__GBLA_SIMD_BLOCK_SIZE; ++j) {
      a = block_A[i*__GBLA_SIMD_BLOCK_SIZE+j];
      if (a != 0) {
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
          wide_block[i][k]  +=
            a * (re_m_t)block_B[j*__GBLA_SIMD_BLOCK_SIZE + k];
        }
      }
    }
  }
#endif
}

/**
 * \brief Uses hybrid blocks from A to update rows in wide_block.
 * Delay modulus operations by storing results in wide blocks.
 * 
 * \param compressed dense block from A block_A
 *
 * \param wide block storing the result wide_block
 *
 * \param modulus resp. field characteristic modulus
 */
static inline void red_hybrid_triangular(dbl_t **block_A,
    re_l_t **wide_block, mod_t modulus)
{
  bi_t i, j, k, l;
  register re_m_t a;
  //            printf("%u || %d | %d\n", block_A[3][0].val[0], 3,0);
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    // if we reach a zero row then the algorithm is in the very last triangular
    // block of A and the number of rows of A does not divide
    // __GBLA_SIMD_BLOCK_SIZE. Thus, this last triangular block is not filled
    // completely and once we reach the very first such row, we are done and can
    // return to the caller.
    if (block_A[i] == NULL) {
      return;
    }
    // we are doing j < i/__GBLA_SIMD_INNER_SIZE+1, but we have to distinguish:
    // 1. For j < i/__GBLA_SIMD_INNER_SIZE we can always go through the
    // full inner block
    // 2. For j =  i/__GBLA_SIMD_INNER_SIZE+1 we have
    // i%__GBLA_SIMD_INNER_SIZE as maximum for k
    for (j=0; j<i/__GBLA_SIMD_INNER_SIZE; ++j) {
      if (block_A[i][j].val == NULL) {
        continue;
      } else {
        for (k=0; k<__GBLA_SIMD_INNER_SIZE; ++k) {
          if (block_A[i][j].val[k] != 0) {
            a = block_A[i][j].val[k];
            if (a != 0) {
              for (l=0; l<__GBLA_SIMD_BLOCK_SIZE; ++l) {
                wide_block[i][l] +=  a *
                  wide_block[j*__GBLA_SIMD_INNER_SIZE+k][l];
              }
            }
          }
        }
      }
    }
    // now the last inner block in this line which we might only go through
    // until i%__GBLA_SIMD_INNER_SIZE
    if (block_A[i][j].val == NULL) {
      continue;
    } else {
      for (k=0; k<i%__GBLA_SIMD_INNER_SIZE; ++k) {
        if (block_A[i][j].val[k] != 0) {
          a = block_A[i][j].val[k];
          if (a != 0) {
            for (l=0; l<__GBLA_SIMD_BLOCK_SIZE; ++l) {
              wide_block[i][l] +=  a *
                wide_block[j*__GBLA_SIMD_INNER_SIZE+k][l];
            }
          }
        }
      }
    }
    modulo_wide_block_val(wide_block, i, modulus);
  }
}

/**
 * \brief Use hybrid blocks from A to update hybrid blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 *
 * \param hybrid block from A block_A
 *
 * \param hybrid block from B block_B
 *
 * \param wide block storing the result wide_block
 */
static inline void red_hybrid_rectangular(dbl_t **block_A, dbl_t **block_B,
  re_l_t **wide_block)
{
  
  bi_t i, j, k, l, m;
  register re_m_t a;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    if (block_A[i] == NULL) {
      continue;
    } else {
      for (j=0; j<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++j) {
        if (block_A[i][j].val == NULL) {
          continue;
        } else {
          for (k=0; k<__GBLA_SIMD_INNER_SIZE; ++k) {
            if (block_A[i][j].val[k] != 0) {
              a = block_A[i][j].val[k];
              if ((a != 0) && (block_B[j*__GBLA_SIMD_INNER_SIZE+k] != NULL)) {
                for (l=0; l<__GBLA_SIMD_INNER_BLOCKS_PER_ROW; ++l){
                  if (block_B[j*__GBLA_SIMD_INNER_SIZE+k][l].val != NULL) {
                    for (m=0; m<__GBLA_SIMD_INNER_SIZE; ++m) {
                      wide_block[i][l*__GBLA_SIMD_INNER_SIZE+m] +=  a *
                        (re_m_t)block_B[j*__GBLA_SIMD_INNER_SIZE+k][l].val[m];
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

/**
 * \brief Takes wide_row and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param column position in wide_row resp. corresponding row index in A idx
 */
static inline void update_wide_row(re_l_t *wide_row1, const sm_fl_t *A,
  const re_m_t multiplier1, const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2;
  register ci_t a_idx1, a_idx2;
  register const re_m_t m1  = multiplier1;
  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
  if (A->sz[row_idx] > 1) {
    for (; i<A->sz[row_idx]-1; i=i+2) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
  }
}

/**
 * \brief Takes two wide rows and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 *
 * \param wide row storing the result wide_row2
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param multiplier from C multiplier2
 *
 * \param column position in wide_row resp. corresponding row index in A pos
 */
static inline void update_wide_rows(re_l_t *wide_row1, re_l_t *wide_row2,
    const sm_fl_t *A, const re_m_t multiplier1, const re_m_t multiplier2,
    const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2;
  register ci_t a_idx1, a_idx2;
  register const re_m_t m1  = multiplier1;
  register const re_m_t m2  = multiplier2;

  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
  if (A->sz[row_idx] > 1) {
    for (; i<A->sz[row_idx]-1; i=i+2) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row2[a_idx2] += m2 * a_elt2;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
    wide_row2[a_idx1] += m2 * a_elt1;
  }
}

/**
 * \brief Takes three wide rows and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 *
 * \param wide row storing the result wide_row2
 *
 * \param wide row storing the result wide_row3
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param multiplier from C multiplier2
 *
 * \param multiplier from C multiplier3
 *
 * \param column position in wide_row resp. corresponding row index in A pos
 */
static inline void update_wide_rows_three(re_l_t *wide_row1, re_l_t *wide_row2,
    re_l_t *wide_row3, const sm_fl_t *A, const re_m_t multiplier1,
    const re_m_t multiplier2, const re_m_t multiplier3, const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2;
  register ci_t a_idx1, a_idx2;
  register const re_m_t m1  = multiplier1;
  register const re_m_t m2  = multiplier2;
  register const re_m_t m3  = multiplier3;
  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
  if (A->sz[row_idx] > 1) {
    for (; i<A->sz[row_idx]-1; i=i+2) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row2[a_idx2] += m2 * a_elt2;
      wide_row3[a_idx1] += m3 * a_elt1;
      wide_row3[a_idx2] += m3 * a_elt2;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
    wide_row2[a_idx1] += m2 * a_elt1;
    wide_row3[a_idx1] += m3 * a_elt1;
  }
}

/**
 * \brief Takes ten wide rows and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 *
 * \param wide row storing the result wide_row2
 *
 * \param wide row storing the result wide_row3
 *
 * \param wide row storing the result wide_row4
 *
 * \param wide row storing the result wide_row5
 *
 * \param wide row storing the result wide_row6
 *
 * \param wide row storing the result wide_row7
 *
 * \param wide row storing the result wide_row8
 *
 * \param wide row storing the result wide_row9
 *
 * \param wide row storing the result wide_row10
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param multiplier from C multiplier2
 *
 * \param multiplier from C multiplier3
 *
 * \param multiplier from C multiplier4
 *
 * \param multiplier from C multiplier5
 *
 * \param multiplier from C multiplier6
 *
 * \param multiplier from C multiplier7
 *
 * \param multiplier from C multiplier8
 *
 * \param multiplier from C multiplier9
 *
 * \param multiplier from C multiplier10
 *
 * \param column position in wide_row resp. corresponding row index in A pos
 */
static inline void update_wide_rows_ten(re_l_t *wide_row1, re_l_t *wide_row2,
    re_l_t *wide_row3, re_l_t *wide_row4, re_l_t *wide_row5, re_l_t *wide_row6,
    re_l_t *wide_row7, re_l_t *wide_row8, re_l_t *wide_row9, re_l_t *wide_row10,
    const sm_fl_t *A,
    const re_m_t multiplier1, const re_m_t multiplier2, const re_m_t multiplier3,
    const re_m_t multiplier4, const re_m_t multiplier5, const re_m_t multiplier6,
    const re_m_t multiplier7, const re_m_t multiplier8, const re_m_t multiplier9,
    const re_m_t multiplier10, const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2;
  register ci_t a_idx1, a_idx2;
  register const re_m_t m1  = multiplier1;
  register const re_m_t m2  = multiplier2;
  register const re_m_t m3  = multiplier3;
  register const re_m_t m4  = multiplier4;
  register const re_m_t m5  = multiplier5;
  register const re_m_t m6  = multiplier6;
  register const re_m_t m7  = multiplier7;
  register const re_m_t m8  = multiplier8;
  register const re_m_t m9  = multiplier9;
  register const re_m_t m10 = multiplier10;
  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
  if (A->sz[row_idx] > 1) {
    for (; i<A->sz[row_idx]-1; i=i+2) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row2[a_idx2] += m2 * a_elt2;
      wide_row3[a_idx1] += m3 * a_elt1;
      wide_row3[a_idx2] += m3 * a_elt2;
      wide_row4[a_idx1] += m4 * a_elt1;
      wide_row4[a_idx2] += m4 * a_elt2;
      wide_row5[a_idx1] += m5 * a_elt1;
      wide_row5[a_idx2] += m5 * a_elt2;
      wide_row6[a_idx1] += m6 * a_elt1;
      wide_row6[a_idx2] += m6 * a_elt2;
      wide_row7[a_idx1] += m7 * a_elt1;
      wide_row7[a_idx2] += m7 * a_elt2;
      wide_row8[a_idx1] += m8 * a_elt1;
      wide_row8[a_idx2] += m8 * a_elt2;
      wide_row9[a_idx1] += m9 * a_elt1;
      wide_row9[a_idx2] += m9 * a_elt2;
      wide_row10[a_idx1] += m10 * a_elt1;
      wide_row10[a_idx2] += m10 * a_elt2;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
    wide_row2[a_idx1] += m2 * a_elt1;
    wide_row3[a_idx1] += m3 * a_elt1;
    wide_row4[a_idx1] += m4 * a_elt1;
    wide_row5[a_idx1] += m5 * a_elt1;
    wide_row6[a_idx1] += m6 * a_elt1;
    wide_row7[a_idx1] += m7 * a_elt1;
    wide_row8[a_idx1] += m8 * a_elt1;
    wide_row9[a_idx1] += m9 * a_elt1;
    wide_row10[a_idx1] += m10 * a_elt1;
  }
}

/**
 * \brief Takes nine wide rows and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 *
 * \param wide row storing the result wide_row2
 *
 * \param wide row storing the result wide_row3
 *
 * \param wide row storing the result wide_row4
 *
 * \param wide row storing the result wide_row5
 *
 * \param wide row storing the result wide_row6
 *
 * \param wide row storing the result wide_row7
 *
 * \param wide row storing the result wide_row8
 *
 * \param wide row storing the result wide_row9
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param multiplier from C multiplier2
 *
 * \param multiplier from C multiplier3
 *
 * \param multiplier from C multiplier4
 *
 * \param multiplier from C multiplier5
 *
 * \param multiplier from C multiplier6
 *
 * \param multiplier from C multiplier7
 *
 * \param multiplier from C multiplier8
 *
 * \param multiplier from C multiplier9
 *
 * \param column position in wide_row resp. corresponding row index in A pos
 */
static inline void update_wide_rows_nine(re_l_t *wide_row1, re_l_t *wide_row2,
    re_l_t *wide_row3, re_l_t *wide_row4, re_l_t *wide_row5, re_l_t *wide_row6,
    re_l_t *wide_row7, re_l_t *wide_row8, re_l_t *wide_row9, const sm_fl_t *A,
    const re_m_t multiplier1, const re_m_t multiplier2, const re_m_t multiplier3,
    const re_m_t multiplier4, const re_m_t multiplier5, const re_m_t multiplier6,
    const re_m_t multiplier7, const re_m_t multiplier8, const re_m_t multiplier9,
    const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2;
  register ci_t a_idx1, a_idx2;
  register const re_m_t m1  = multiplier1;
  register const re_m_t m2  = multiplier2;
  register const re_m_t m3  = multiplier3;
  register const re_m_t m4  = multiplier4;
  register const re_m_t m5  = multiplier5;
  register const re_m_t m6  = multiplier6;
  register const re_m_t m7  = multiplier7;
  register const re_m_t m8  = multiplier8;
  register const re_m_t m9  = multiplier9;
  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
  if (A->sz[row_idx] > 1) {
    for (; i<A->sz[row_idx]-1; i=i+2) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row2[a_idx2] += m2 * a_elt2;
      wide_row3[a_idx1] += m3 * a_elt1;
      wide_row3[a_idx2] += m3 * a_elt2;
      wide_row4[a_idx1] += m4 * a_elt1;
      wide_row4[a_idx2] += m4 * a_elt2;
      wide_row5[a_idx1] += m5 * a_elt1;
      wide_row5[a_idx2] += m5 * a_elt2;
      wide_row6[a_idx1] += m6 * a_elt1;
      wide_row6[a_idx2] += m6 * a_elt2;
      wide_row7[a_idx1] += m7 * a_elt1;
      wide_row7[a_idx2] += m7 * a_elt2;
      wide_row8[a_idx1] += m8 * a_elt1;
      wide_row8[a_idx2] += m8 * a_elt2;
      wide_row9[a_idx1] += m9 * a_elt1;
      wide_row9[a_idx2] += m9 * a_elt2;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
    wide_row2[a_idx1] += m2 * a_elt1;
    wide_row3[a_idx1] += m3 * a_elt1;
    wide_row4[a_idx1] += m4 * a_elt1;
    wide_row5[a_idx1] += m5 * a_elt1;
    wide_row6[a_idx1] += m6 * a_elt1;
    wide_row7[a_idx1] += m7 * a_elt1;
    wide_row8[a_idx1] += m8 * a_elt1;
    wide_row9[a_idx1] += m9 * a_elt1;
  }
}

/**
 * \brief Takes eight wide rows and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 *
 * \param wide row storing the result wide_row2
 *
 * \param wide row storing the result wide_row3
 *
 * \param wide row storing the result wide_row4
 *
 * \param wide row storing the result wide_row5
 *
 * \param wide row storing the result wide_row6
 *
 * \param wide row storing the result wide_row7
 *
 * \param wide row storing the result wide_row8
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param multiplier from C multiplier2
 *
 * \param multiplier from C multiplier3
 *
 * \param multiplier from C multiplier4
 *
 * \param multiplier from C multiplier5
 *
 * \param multiplier from C multiplier6
 *
 * \param multiplier from C multiplier7
 *
 * \param multiplier from C multiplier8
 *
 * \param column position in wide_row resp. corresponding row index in A pos
 */
static inline void update_wide_rows_eight(re_l_t *wide_row1, re_l_t *wide_row2,
    re_l_t *wide_row3, re_l_t *wide_row4, re_l_t *wide_row5, re_l_t *wide_row6,
    re_l_t *wide_row7, re_l_t *wide_row8, const sm_fl_t *A, const re_m_t multiplier1,
    const re_m_t multiplier2, const re_m_t multiplier3, const re_m_t multiplier4,
    const re_m_t multiplier5, const re_m_t multiplier6, const re_m_t multiplier7,
    const re_m_t multiplier8, const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2;
  register ci_t a_idx1, a_idx2;
  register const re_m_t m1  = multiplier1;
  register const re_m_t m2  = multiplier2;
  register const re_m_t m3  = multiplier3;
  register const re_m_t m4  = multiplier4;
  register const re_m_t m5  = multiplier5;
  register const re_m_t m6  = multiplier6;
  register const re_m_t m7  = multiplier7;
  register const re_m_t m8  = multiplier8;
  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
  if (A->sz[row_idx] > 1) {
    for (; i<A->sz[row_idx]-1; i=i+2) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row2[a_idx2] += m2 * a_elt2;
      wide_row3[a_idx1] += m3 * a_elt1;
      wide_row3[a_idx2] += m3 * a_elt2;
      wide_row4[a_idx1] += m4 * a_elt1;
      wide_row4[a_idx2] += m4 * a_elt2;
      wide_row5[a_idx1] += m5 * a_elt1;
      wide_row5[a_idx2] += m5 * a_elt2;
      wide_row6[a_idx1] += m6 * a_elt1;
      wide_row6[a_idx2] += m6 * a_elt2;
      wide_row7[a_idx1] += m7 * a_elt1;
      wide_row7[a_idx2] += m7 * a_elt2;
      wide_row8[a_idx1] += m8 * a_elt1;
      wide_row8[a_idx2] += m8 * a_elt2;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
    wide_row2[a_idx1] += m2 * a_elt1;
    wide_row3[a_idx1] += m3 * a_elt1;
    wide_row4[a_idx1] += m4 * a_elt1;
    wide_row5[a_idx1] += m5 * a_elt1;
    wide_row6[a_idx1] += m6 * a_elt1;
    wide_row7[a_idx1] += m7 * a_elt1;
    wide_row8[a_idx1] += m8 * a_elt1;
  }
}

/**
 * \brief Takes seven wide rows and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 *
 * \param wide row storing the result wide_row2
 *
 * \param wide row storing the result wide_row3
 *
 * \param wide row storing the result wide_row4
 *
 * \param wide row storing the result wide_row5
 *
 * \param wide row storing the result wide_row6
 *
 * \param wide row storing the result wide_row7
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param multiplier from C multiplier2
 *
 * \param multiplier from C multiplier3
 *
 * \param multiplier from C multiplier4
 *
 * \param multiplier from C multiplier5
 *
 * \param multiplier from C multiplier6
 *
 * \param multiplier from C multiplier7
 *
 * \param column position in wide_row resp. corresponding row index in A pos
 */
static inline void update_wide_rows_seven(re_l_t *wide_row1, re_l_t *wide_row2,
    re_l_t *wide_row3, re_l_t *wide_row4, re_l_t *wide_row5, re_l_t *wide_row6,
    re_l_t *wide_row7, const sm_fl_t *A, const re_m_t multiplier1,
    const re_m_t multiplier2, const re_m_t multiplier3, const re_m_t multiplier4,
    const re_m_t multiplier5, const re_m_t multiplier6, const re_m_t multiplier7,
    const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2;
  register ci_t a_idx1, a_idx2;
  register const re_m_t m1  = multiplier1;
  register const re_m_t m2  = multiplier2;
  register const re_m_t m3  = multiplier3;
  register const re_m_t m4  = multiplier4;
  register const re_m_t m5  = multiplier5;
  register const re_m_t m6  = multiplier6;
  register const re_m_t m7  = multiplier7;
  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
  if (A->sz[row_idx] > 1) {
    for (; i<A->sz[row_idx]-1; i=i+2) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row2[a_idx2] += m2 * a_elt2;
      wide_row3[a_idx1] += m3 * a_elt1;
      wide_row3[a_idx2] += m3 * a_elt2;
      wide_row4[a_idx1] += m4 * a_elt1;
      wide_row4[a_idx2] += m4 * a_elt2;
      wide_row5[a_idx1] += m5 * a_elt1;
      wide_row5[a_idx2] += m5 * a_elt2;
      wide_row6[a_idx1] += m6 * a_elt1;
      wide_row6[a_idx2] += m6 * a_elt2;
      wide_row7[a_idx1] += m7 * a_elt1;
      wide_row7[a_idx2] += m7 * a_elt2;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
    wide_row2[a_idx1] += m2 * a_elt1;
    wide_row3[a_idx1] += m3 * a_elt1;
    wide_row4[a_idx1] += m4 * a_elt1;
    wide_row5[a_idx1] += m5 * a_elt1;
    wide_row6[a_idx1] += m6 * a_elt1;
    wide_row7[a_idx1] += m7 * a_elt1;
  }
}

/**
 * \brief Takes six wide rows and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 *
 * \param wide row storing the result wide_row2
 *
 * \param wide row storing the result wide_row3
 *
 * \param wide row storing the result wide_row4
 *
 * \param wide row storing the result wide_row5
 *
 * \param wide row storing the result wide_row6
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param multiplier from C multiplier2
 *
 * \param multiplier from C multiplier3
 *
 * \param multiplier from C multiplier4
 *
 * \param multiplier from C multiplier5
 *
 * \param multiplier from C multiplier6
 *
 * \param column position in wide_row resp. corresponding row index in A pos
 */
static inline void update_wide_rows_six(re_l_t *wide_row1, re_l_t *wide_row2,
    re_l_t *wide_row3, re_l_t *wide_row4, re_l_t *wide_row5, re_l_t *wide_row6,
    const sm_fl_t *A, const re_m_t multiplier1, const re_m_t multiplier2,
    const re_m_t multiplier3, const re_m_t multiplier4, const re_m_t multiplier5,
    const re_m_t multiplier6, const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2;
  register ci_t a_idx1, a_idx2;
  register const re_m_t m1  = multiplier1;
  register const re_m_t m2  = multiplier2;
  register const re_m_t m3  = multiplier3;
  register const re_m_t m4  = multiplier4;
  register const re_m_t m5  = multiplier5;
  register const re_m_t m6  = multiplier6;
  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
  if (A->sz[row_idx] > 1) {
    for (; i<A->sz[row_idx]-1; i=i+2) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row2[a_idx2] += m2 * a_elt2;
      wide_row3[a_idx1] += m3 * a_elt1;
      wide_row3[a_idx2] += m3 * a_elt2;
      wide_row4[a_idx1] += m4 * a_elt1;
      wide_row4[a_idx2] += m4 * a_elt2;
      wide_row5[a_idx1] += m5 * a_elt1;
      wide_row5[a_idx2] += m5 * a_elt2;
      wide_row6[a_idx1] += m6 * a_elt1;
      wide_row6[a_idx2] += m6 * a_elt2;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
    wide_row2[a_idx1] += m2 * a_elt1;
    wide_row3[a_idx1] += m3 * a_elt1;
    wide_row4[a_idx1] += m4 * a_elt1;
    wide_row5[a_idx1] += m5 * a_elt1;
    wide_row6[a_idx1] += m6 * a_elt1;
  }
}

/**
 * \brief Takes five wide rows and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 *
 * \param wide row storing the result wide_row2
 *
 * \param wide row storing the result wide_row3
 *
 * \param wide row storing the result wide_row4
 *
 * \param wide row storing the result wide_row5
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param multiplier from C multiplier2
 *
 * \param multiplier from C multiplier3
 *
 * \param multiplier from C multiplier4
 *
 * \param multiplier from C multiplier5
 *
 * \param column position in wide_row resp. corresponding row index in A pos
 */
static inline void update_wide_rows_five(re_l_t *wide_row1, re_l_t *wide_row2,
    re_l_t *wide_row3, re_l_t *wide_row4, re_l_t *wide_row5, const sm_fl_t *A,
    const re_m_t multiplier1, const re_m_t multiplier2, const re_m_t multiplier3,
    const re_m_t multiplier4, const re_m_t multiplier5, const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2;
  register ci_t a_idx1, a_idx2;
  register const re_m_t m1  = multiplier1;
  register const re_m_t m2  = multiplier2;
  register const re_m_t m3  = multiplier3;
  register const re_m_t m4  = multiplier4;
  register const re_m_t m5  = multiplier5;
  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
  if (A->sz[row_idx] > 1) {
    for (; i<A->sz[row_idx]-1; i=i+2) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row2[a_idx2] += m2 * a_elt2;
      wide_row3[a_idx1] += m3 * a_elt1;
      wide_row3[a_idx2] += m3 * a_elt2;
      wide_row4[a_idx1] += m4 * a_elt1;
      wide_row4[a_idx2] += m4 * a_elt2;
      wide_row5[a_idx1] += m5 * a_elt1;
      wide_row5[a_idx2] += m5 * a_elt2;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
    wide_row2[a_idx1] += m2 * a_elt1;
    wide_row3[a_idx1] += m3 * a_elt1;
    wide_row4[a_idx1] += m4 * a_elt1;
    wide_row5[a_idx1] += m5 * a_elt1;
  }
}

/**
 * \brief Takes four wide rows and updates it with multiples of a row from A
 * corresponding to the column position pos defined by the first non zero entry
 * in wide_row.
 *
 * \param wide row storing the result wide_row1
 *
 * \param wide row storing the result wide_row2
 *
 * \param wide row storing the result wide_row3
 *
 * \param wide row storing the result wide_row4
 * 
 * \param sparse block matrix to update wide row A
 *
 * \param multiplier from C multiplier1
 *
 * \param multiplier from C multiplier2
 *
 * \param multiplier from C multiplier3
 *
 * \param multiplier from C multiplier4
 *
 * \param column position in wide_row resp. corresponding row index in A pos
 */
static inline void update_wide_rows_four(re_l_t *wide_row1, re_l_t *wide_row2,
    re_l_t *wide_row3, re_l_t *wide_row4, const sm_fl_t *A, const re_m_t multiplier1,
    const re_m_t multiplier2, const re_m_t multiplier3, const re_m_t multiplier4,
    const ci_t idx)
{
  ci_t i;
  const ri_t row_idx  = A->nrows - idx - 1;
  register re_m_t a_elt1, a_elt2, a_elt3, a_elt4;
  register ci_t a_idx1, a_idx2, a_idx3, a_idx4;
  register const re_m_t m1  = multiplier1;
  register const re_m_t m2  = multiplier2;
  register const re_m_t m3  = multiplier3;
  register const re_m_t m4  = multiplier4;
  
  // we do not need to update with A->blocks[bir][bir].row[rib][0] since we
  // cancel out the element at this position in wide_row
  i = 0;
 /*
  if (A->sz[row_idx] > 7) {
    for (; i<A->sz[row_idx]-7; i=i+8) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      a_idx3  = A->pos[row_idx][i+2];
      a_elt3  = A->row[row_idx][i+2];
      a_idx4  = A->pos[row_idx][i+3];
      a_elt4  = A->row[row_idx][i+3];
      a_idx5  = A->pos[row_idx][i+4];
      a_elt5  = A->row[row_idx][i+4];
      a_idx6  = A->pos[row_idx][i+5];
      a_elt6  = A->row[row_idx][i+5];
      a_idx7  = A->pos[row_idx][i+6];
      a_elt7  = A->row[row_idx][i+6];
      a_idx8  = A->pos[row_idx][i+7];
      a_elt8  = A->row[row_idx][i+7];

      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row1[a_idx3] += m1 * a_elt3;
      wide_row1[a_idx4] += m1 * a_elt4;
      wide_row1[a_idx5] += m1 * a_elt5;
      wide_row1[a_idx6] += m1 * a_elt6;
      wide_row1[a_idx7] += m1 * a_elt7;
      wide_row1[a_idx8] += m1 * a_elt8;
      
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row2[a_idx2] += m2 * a_elt2;
      wide_row2[a_idx3] += m2 * a_elt3;
      wide_row2[a_idx4] += m2 * a_elt4;
      wide_row2[a_idx5] += m2 * a_elt5;
      wide_row2[a_idx6] += m2 * a_elt6;
      wide_row2[a_idx7] += m2 * a_elt7;
      wide_row2[a_idx8] += m2 * a_elt8;

      wide_row3[a_idx1] += m3 * a_elt1;
      wide_row3[a_idx2] += m3 * a_elt2;
      wide_row3[a_idx3] += m3 * a_elt3;
      wide_row3[a_idx4] += m3 * a_elt4;
      wide_row3[a_idx5] += m3 * a_elt5;
      wide_row3[a_idx6] += m3 * a_elt6;
      wide_row3[a_idx7] += m3 * a_elt7;
      wide_row3[a_idx8] += m3 * a_elt8;

      wide_row4[a_idx1] += m4 * a_elt1;
      wide_row4[a_idx2] += m4 * a_elt2;
      wide_row4[a_idx3] += m4 * a_elt3;
      wide_row4[a_idx4] += m4 * a_elt4;
      wide_row4[a_idx5] += m4 * a_elt5;
      wide_row4[a_idx6] += m4 * a_elt6;
      wide_row4[a_idx7] += m4 * a_elt7;
      wide_row4[a_idx8] += m4 * a_elt8;
    }
  }
  */
  if (A->sz[row_idx] > 3) {
    for (; i<A->sz[row_idx]-3; i=i+4) {
      a_idx1  = A->pos[row_idx][i];
      a_elt1  = A->row[row_idx][i];
      a_idx2  = A->pos[row_idx][i+1];
      a_elt2  = A->row[row_idx][i+1];
      a_idx3  = A->pos[row_idx][i+2];
      a_elt3  = A->row[row_idx][i+2];
      a_idx4  = A->pos[row_idx][i+3];
      a_elt4  = A->row[row_idx][i+3];

      wide_row1[a_idx1] += m1 * a_elt1;
      wide_row1[a_idx2] += m1 * a_elt2;
      wide_row1[a_idx3] += m1 * a_elt3;
      wide_row1[a_idx4] += m1 * a_elt4;
      
      wide_row2[a_idx1] += m2 * a_elt1;
      wide_row2[a_idx2] += m2 * a_elt2;
      wide_row2[a_idx3] += m2 * a_elt3;
      wide_row2[a_idx4] += m2 * a_elt4;

      wide_row3[a_idx1] += m3 * a_elt1;
      wide_row3[a_idx2] += m3 * a_elt2;
      wide_row3[a_idx3] += m3 * a_elt3;
      wide_row3[a_idx4] += m3 * a_elt4;

      wide_row4[a_idx1] += m4 * a_elt1;
      wide_row4[a_idx2] += m4 * a_elt2;
      wide_row4[a_idx3] += m4 * a_elt3;
      wide_row4[a_idx4] += m4 * a_elt4;
    }
  }
  for (; i<A->sz[row_idx]; ++i) {
    a_idx1  = A->pos[row_idx][i];
    a_elt1  = A->row[row_idx][i];
    wide_row1[a_idx1] += m1 * a_elt1;
    wide_row2[a_idx1] += m2 * a_elt1;
    wide_row3[a_idx1] += m3 * a_elt1;
    wide_row4[a_idx1] += m4 * a_elt1;
  }
}


/**
 * \brief Use sparse blocks from A to update dense blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 * 
 * \param sparse block from A block_A
 *
 * \param wide block storing the result wide_block
 *
 * \param modulus resp. field characteristic modulus
 */
static inline void red_sparse_triangular(const sbl_t *block_A,
  re_l_t **wide_block, const mod_t modulus)
{
  bi_t i, j, k;
#if 0
  register re_m_t a, b ,c, d, e, f, g, h;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (j=0; j<block_A->sz[i]-8; j = j+8) {
      a = block_A->val[i][j];
      b = block_A->val[i][j+1];
      c = block_A->val[i][j+2];
      d = block_A->val[i][j+3];
      e = block_A->val[i][j+4];
      f = block_A->val[i][j+5];
      g = block_A->val[i][j+6];
      h = block_A->val[i][j+7];
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
        wide_block[i][k]  +=
          (a * wide_block[block_A->pos[i][j]][k] +
          b * wide_block[block_A->pos[i][j+1]][k] +
          c * wide_block[block_A->pos[i][j+2]][k] +
          d * wide_block[block_A->pos[i][j+3]][k] +
          e * wide_block[block_A->pos[i][j+4]][k] +
          f * wide_block[block_A->pos[i][j+5]][k] +
          g * wide_block[block_A->pos[i][j+6]][k] +
          h * wide_block[block_A->pos[i][j+7]][k]);
      }
    }
    for (;j<block_A->sz[i]-1; ++j) {
      a = block_A->val[i][j];
      for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
        wide_block[i][k]  +=
          a * wide_block[block_A->pos[i][j]][k];
      }
    }
    modulo_wide_block_val(wide_block, i, modulus);
  }
#else
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    const re_t *bAv = block_A->val[i];
    const bi_t *bAp = block_A->pos[i];
    j = 0;
    if (block_A->sz[i] > 0) {
      if (block_A->sz[i] > 8) {
        for (; j<block_A->sz[i]-8; j = j+8) {
          for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k=k+4) {
            wide_block[i][k]  +=
              ((re_l_t)bAv[j] * wide_block[bAp[j]][k] +
               (re_m_t)bAv[j+1] * wide_block[bAp[j+1]][k] +
               (re_m_t)bAv[j+2] * wide_block[bAp[j+2]][k] +
               (re_m_t)bAv[j+3] * wide_block[bAp[j+3]][k] +
               (re_m_t)bAv[j+4] * wide_block[bAp[j+4]][k] +
               (re_m_t)bAv[j+5] * wide_block[bAp[j+5]][k] +
               (re_m_t)bAv[j+6] * wide_block[bAp[j+6]][k] +
               (re_m_t)bAv[j+7] * wide_block[bAp[j+7]][k]);
            wide_block[i][k+1]  +=
              ((re_l_t)bAv[j] * wide_block[bAp[j]][k+1] +
               (re_m_t)bAv[j+1] * wide_block[bAp[j+1]][k+1] +
               (re_m_t)bAv[j+2] * wide_block[bAp[j+2]][k+1] +
               (re_m_t)bAv[j+3] * wide_block[bAp[j+3]][k+1] +
               (re_m_t)bAv[j+4] * wide_block[bAp[j+4]][k+1] +
               (re_m_t)bAv[j+5] * wide_block[bAp[j+5]][k+1] +
               (re_m_t)bAv[j+6] * wide_block[bAp[j+6]][k+1] +
               (re_m_t)bAv[j+7] * wide_block[bAp[j+7]][k+1]);
            wide_block[i][k+2]  +=
              ((re_l_t)bAv[j] * wide_block[bAp[j]][k+2] +
               (re_m_t)bAv[j+1] * wide_block[bAp[j+1]][k+2] +
               (re_m_t)bAv[j+2] * wide_block[bAp[j+2]][k+2] +
               (re_m_t)bAv[j+3] * wide_block[bAp[j+3]][k+2] +
               (re_m_t)bAv[j+4] * wide_block[bAp[j+4]][k+2] +
               (re_m_t)bAv[j+5] * wide_block[bAp[j+5]][k+2] +
               (re_m_t)bAv[j+6] * wide_block[bAp[j+6]][k+2] +
               (re_m_t)bAv[j+7] * wide_block[bAp[j+7]][k+2]);
            wide_block[i][k+3]  +=
              ((re_l_t)bAv[j] * wide_block[bAp[j]][k+3] +
               (re_m_t)bAv[j+1] * wide_block[bAp[j+1]][k+3] +
               (re_m_t)bAv[j+2] * wide_block[bAp[j+2]][k+3] +
               (re_m_t)bAv[j+3] * wide_block[bAp[j+3]][k+3] +
               (re_m_t)bAv[j+4] * wide_block[bAp[j+4]][k+3] +
               (re_m_t)bAv[j+5] * wide_block[bAp[j+5]][k+3] +
               (re_m_t)bAv[j+6] * wide_block[bAp[j+6]][k+3] +
               (re_m_t)bAv[j+7] * wide_block[bAp[j+7]][k+3]);
          }
        }
      }
      if (block_A->sz[i]-j > 4) {
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k=k+4) {
          wide_block[i][k]  +=
            ((re_l_t)bAv[j] * wide_block[bAp[j]][k] +
             (re_m_t)bAv[j+1] * wide_block[bAp[j+1]][k] +
             (re_m_t)bAv[j+2] * wide_block[bAp[j+2]][k] +
             (re_m_t)bAv[j+3] * wide_block[bAp[j+3]][k]);
          wide_block[i][k+1]  +=
            ((re_l_t)bAv[j] * wide_block[bAp[j]][k+1] +
             (re_m_t)bAv[j+1] * wide_block[bAp[j+1]][k+1] +
             (re_m_t)bAv[j+2] * wide_block[bAp[j+2]][k+1] +
             (re_m_t)bAv[j+3] * wide_block[bAp[j+3]][k+1]);
          wide_block[i][k+2]  +=
            ((re_l_t)bAv[j] * wide_block[bAp[j]][k+2] +
             (re_m_t)bAv[j+1] * wide_block[bAp[j+1]][k+2] +
             (re_m_t)bAv[j+2] * wide_block[bAp[j+2]][k+2] +
             (re_m_t)bAv[j+3] * wide_block[bAp[j+3]][k+2]);
          wide_block[i][k+3]  +=
            ((re_l_t)bAv[j] * wide_block[bAp[j]][k+3] +
             (re_m_t)bAv[j+1] * wide_block[bAp[j+1]][k+3] +
             (re_m_t)bAv[j+2] * wide_block[bAp[j+2]][k+3] +
             (re_m_t)bAv[j+3] * wide_block[bAp[j+3]][k+3]);
        }
        j = j+4;
      }
      //printf("j %u | block_A->sz[%u] = %u\n", j, i, block_A->sz[i]);
      for (;j<block_A->sz[i]-1; ++j) {
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; k=k+4) {
          wide_block[i][k]  +=
            (re_m_t)bAv[j] * wide_block[bAp[j]][k];
          wide_block[i][k+1]  +=
            (re_m_t)bAv[j] * wide_block[bAp[j]][k+1];
          wide_block[i][k+2]  +=
            (re_m_t)bAv[j] * wide_block[bAp[j]][k+2];
          wide_block[i][k+3]  +=
            (re_m_t)bAv[j] * wide_block[bAp[j]][k+3];
        }
      }
      modulo_wide_block_val(wide_block, i, modulus);
    }
  }
#endif
}

/**
 * \brief Use compressed dense blocks from A to update dense blocks in B.
 * Delay modulus operations by storing results in wide blocks.
 * 
 * \param compressed dense block from A block_A
 *
 * \param wide block storing the result wide_block
 *
 * \param modulus resp. field characteristic modulus
 */
static inline void red_dense_triangular(const re_t *block_A,
  re_l_t **wide_block, const mod_t modulus)
{
  //printf("INDENSE\n");
  bi_t i, j, k;
  register re_m_t a;
  for (i=0; i<__GBLA_SIMD_BLOCK_SIZE; ++i) {
    for (j=0; j<i; ++j) {
      // Note: diagonal blocks in A are compressed, only lower triangular
      // entries are stored due to inversion of ordering of the elements when
      // storing them in src/mapping.h
      if (block_A[i + ((2*__GBLA_SIMD_BLOCK_SIZE-(j+1))*j)/2] != 0) {
        a = block_A[i + ((2*__GBLA_SIMD_BLOCK_SIZE-(j+1))*j)/2];
        //printf("%lu - %u | %u\n",a,i,j);
       // printf("%lu || %d | %d\n", a, i, j);
              //printf("%u || %d | %d\n", a, i,j);
        for (k=0; k<__GBLA_SIMD_BLOCK_SIZE; ++k) {
          /*
          if (i==83 && k==0) {
            printf("%lu !! %lu, %d\n", a, wide_block[83][0], j);
          }
          */
          wide_block[i][k]  +=  a * wide_block[j][k];
        }
        //printf("---> %lu\n",wide_block[83][0]);
      }
    }
    modulo_wide_block_val(wide_block, i, modulus);
  }
}

/**
 * \brief Copies entry from sparse block sparse_block to dense block dense_block
 * for elimination purposes.
 *
 * \param sparse block sparse_block
 *
 * \param dense block dense_block
 *
 * \param block height bheight
 *
 * \param block width bwidth
 */
void copy_sparse_to_dense_block(
		mbl_t *sparse_block, re_l_t **dense_block, ri_t bheight, ci_t bwidth) ;

/**
 * \brief Copies entry from dense block to sparse rows
 * for elimination purposes.
 *
 * \param dense block dense_block
 *
 * \param sparse block sparse_block
 *
 * \param block height bheight
 *
 * \param block width bwidth
 *
 * \param modulus resp. field characteristic modulus
 */
void copy_dense_block_to_sparse(
		re_l_t **dense_block, mbl_t **sparse_block, ri_t bheight, ci_t bwidth,
		mod_t modulus);

/**
 * \brief Computes a dense AXPY for one row
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param corresponding multiline block in B
 *
 * \param line index line_idx
 *
 * \param dense block width bwidth
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
void dense_scal_mul_sub_1_row_vect_array(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const mbl_t multiline,
		const bi_t line_idx,
		const bi_t  bwidth,
		re_l_t *dense_val1,
		re_l_t *dense_val2) ;

/**
 * \brief Computes a dense AXPY for one row using full multilines
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param corresponding multiline
 *
 * \param line index line_idx
 *
 * \param dense block width bwidth
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 *
 * \param offset for first nonzero coefficient in multiline
 */
void dense_scal_mul_sub_1_row_vect_array_multiline_var_size(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const ml_t multiline,
		const bi_t line_idx,
		re_l_t *dense_val1,
		re_l_t *dense_val2,
		const ci_t offset) ;

/**
 * \brief Computes a dense AXPY for two rows for multilines
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param corresponding multiline
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 *
 * \param offset for first nonzero coefficient in multiline
 */
void dense_scal_mul_sub_2_rows_vect_array_multiline_var_size(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const re_m_t Av1_col2,
		const re_m_t Av2_col2,
		const ml_t multiline,
		re_l_t *dense_val1,
		re_l_t *dense_val2,
		const ci_t offset1,
		const ci_t offset2) ;

/**
 * \brief Computes a sparse AXPY for one row for full multilines
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param corresponding multiline
 *
 * \param line index line_idx
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
void sparse_scal_mul_sub_1_row_vect_array_multiline(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const ml_t multiline,
		const bi_t line_idx,
		re_l_t *dense_val1,
		re_l_t *dense_val2) ;

/**
 * \brief Computes a sparse AXPY for one row
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param corresponding multiline block in B
 *
 * \param line index line_idx
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
void sparse_scal_mul_sub_1_row_vect_array(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const mbl_t multiline,
		const bi_t line_idx,
		re_l_t *dense_val1,
		re_l_t *dense_val2) ;

/**
 * \brief Computes a dense AXPY for two rows
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param corresponding multiline block in B
 *
 * \param dense block width bwidth
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
void dense_scal_mul_sub_2_rows_vect_array(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const re_m_t Av1_col2,
		const re_m_t Av2_col2,
		const mbl_t multiline,
		const bi_t  bwidth,
		re_l_t *dense_val1,
		re_l_t *dense_val2) ;

/**
 * \brief Computes a sparse AXPY for two rows
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param corresponding multiline block in B
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
void sparse_scal_mul_sub_2_rows_vect_array(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const re_m_t Av1_col2,
		const re_m_t Av2_col2,
		const mbl_t multiline,
		re_l_t *dense_val1,
		re_l_t *dense_val2) ;

/**
 * \brief Computes a sparse AXPY for two rows for multilines
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param corresponding multiline
 *
 * \param dense value 1 holder for delayed modulus dense_val1
 *
 * \param dense value 2 holder for delayed modulus dense_val2
 */
void sparse_scal_mul_sub_2_rows_vect_array_multiline(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const re_m_t Av1_col2,
		const re_m_t Av2_col2,
		const ml_t multiline,
		re_l_t *dense_val1,
		re_l_t *dense_val2) ;

/**
 * \brief Computes a dense AXPY for one dense row (in triangular A^-1B situation)
 *
 * \param value 1 from A Av1_col1
 *
 * \param value 2 from A Av2_col1
 *
 * \param block width bwidth
 *
 * \param dense source array dense_array_source
 *
 * \param dense array 1 holder for delayed modulus dense_array1
 *
 * \param dense array 2 holder for delayed modulus dense_array2
 */
void dense_scal_mul_sub_1_row_array_array(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const bi_t bwidth,
		const re_l_t *dense_array_source,
		re_l_t *dense_array1,
		re_l_t *dense_array2) ;

/**
 * \brief Computes a dense AXPY for two dense rows (in triangular A^-1B situation)
 *
 * \param value 1,1 from A Av1_col1
 *
 * \param value 2,1 from A Av2_col1
 *
 * \param value 1,2 from A Av1_col2
 *
 * \param value 2,2 from A Av2_col2
 *
 * \param block width bwidth
 *
 * \param dense source array dense_array_source1
 *
 * \param dense source array dense_array_source2
 *
 * \param dense array 1 holder for delayed modulus dense_array1
 *
 * \param dense array 2 holder for delayed modulus dense_array2
 */
void dense_scal_mul_sub_2_rows_array_array(
		const re_m_t Av1_col1,
		const re_m_t Av2_col1,
		const re_m_t Av1_col2,
		const re_m_t Av2_col2,
		const bi_t bwidth,
		const re_l_t *dense_array_source1,
		const re_l_t *dense_array_source2,
		re_l_t *dense_array1,
		re_l_t *dense_array2);

/**
 * \brief Modular reduction of dense row array
 *
 * \param dense row array to be reduced dense_array
 *
 * \param block width bwidth
 *
 * \param modulus resp. field characteristic modulus
 */
void red_dense_array_modular(re_l_t *dense_array, const bi_t bwidth, const mod_t modulus);

/**
 * \brief Gets first nonzero entry in multiline m at line index line_idx and
 * stores it in h resp. h_idx.
 *
 * \param multiline m
 *
 * \param line index line_idx
 *
 * \param storage for "head" of line h
 *
 * \param storage for "head index" of line h_idx
 *
 * \return index value corresponding to h from input matrix M
 */
mli_t get_head_multiline(const ml_t *m, const bi_t line_idx, re_t *h, mli_t *h_idx) ;

/**
 * \brief Gets first nonzero entry in multiline m at line index line_idx and
 * stores it in h resp. h_idx. Use hybrid representation of multiline.
 *
 * \param multiline m
 *
 * \param line index line_idx
 *
 * \param storage for "head" of line h
 *
 * \param storage for "head index" of line h_idx
 *
 * \return index value corresponding to h from input matrix M
 */
mli_t get_head_multiline_hybrid(const ml_t *m,
		const bi_t line_idx, re_t *h, mli_t *h_idx, const ci_t coldim) ;

/**
 * \brief Gets first nonzero entry in dense array. Reduces to zero on the go if
 * possible.
 *
 * \param dense array dense_array
 *
 * \param holder for nonzero head value
 *
 * \param column dimension coldim
 *
 * \param field characteristic modulus
 *
 * \return index value corresponding to head value index in dense array
 */
int get_head_dense_array(re_l_t *dense_array,
		re_t *val, const ci_t coldim, const mod_t modulus) ;

/**
 * \brief Computes inverse value of x modulo y:
 * We compute the inverse using the extended GCD. So we are only interested in x.
 * Note that internally we need signed types, but we return only unsigned type
 * re_t for x.
 *
 * \param x
 *
 * \param y
 */
void inverse_val(re_t *x, const mod_t modulus) ;

/**
 * \brief Normalizes dense array and returns index of head element in array
 *
 * \param dense array dense_array
 *
 * \param column dimension coldim
 *
 * \param field characteristic modulus
 *
 * \return index value corresponding to head value index in dense array
 */
int normalize_dense_array(re_l_t *dense_array,
		const ci_t coldim, const mod_t modulus);

/**
 * \brief Normalizes multiline vector.
 *
 * \param multiline m
 *
 * \param field characteristic modulus
 */
void normalize_multiline(ml_t *m, const ci_t coldim, const mod_t modulus) ;

/**
 * \brief Copies two dense arrays to a dense multiline for further processing.
 * \note Multiline m must have already memory allocated correspondingly.
 * Moreover all entries of m->val must be initialized to zero beforehand.
 *
 * \param dense array dense_1
 *
 * \param dense array dense_2
 *
 * \param minimal first position of non-zero entries in dense_1 and dense_2
 * start_pos
 *
 * \param multiline m
 *
 * \param array size resp. column dimension coldim
 */
void copy_dense_arrays_to_zero_dense_multiline(const re_l_t *dense_1,
		const re_l_t *dense_2, const int start_pos, ml_t *m, const ci_t coldim,
		const mod_t modulus) ;

/**
 * \brief Copies two dense arrays to a multiline for further processing.
 *
 * \param dense array dense_1
 *
 * \param dense array dense_2
 *
 * \param minimal first position of non-zero entries in dense_1 and dense_2
 * start_pos
 *
 * \param multiline m
 *
 * \param array size resp. column dimension coldim
 */
void copy_dense_arrays_to_multiline(const re_l_t *dense_1,
		const re_l_t *dense_2, const int start_pos, ml_t *m, const ci_t coldim,
		const mod_t modulus) ;

/**
 * \brief Copies one dense array to a dense multiline for further processing.
 * \note Multiline m must have already memory allocated correspondingly.
 * Moreover all entries of m->val must be initialized to zero beforehand.
 *
 * \param dense array dense_1
 *
 * \param position of first non-zero entry in array dense_1 start_pos
 *
 * \param multiline m
 *
 * \param array size resp. column dimension coldim
 */
void copy_dense_array_to_zero_dense_multiline(const re_l_t *dense_1,
		const int start_pos, ml_t *m, const ci_t coldim, const mod_t modulus) ;

#if 0
void copy_dense_arrays_to_dense_multiline(const re_l_t *dense_1,
		const re_l_t *dense_2, ml_t *m, const ci_t coldim, const mod_t modulus) {

	ci_t i;
	re_l_t tmp_1, tmp_2;
	for (i=0; i<coldim; ++i) {
		tmp_1 = MODP(dense_1[i] , modulus);
		tmp_2 = MODP(dense_2[i] , modulus);
		m->val[2*i]   = (re_t)tmp_1;
		m->val[2*i+1] = (re_t)tmp_2;
	}
}
#endif

/**
 * \brief Copies multiline to two dense arrays for further processing.
 *
 * \param multiline m
 *
 * \param dense array dense_1
 *
 * \param dense array dense_2
 *
 * \param array size resp. column dimension coldim
 */
void copy_multiline_to_dense_array(const ml_t m, re_l_t *dense_1,
		re_l_t *dense_2, const ci_t coldim) ;

/**
 * \brief Returns smallest row index in waiting list for echelonization.
 *
 * \param waiting list waiting
 *
 * \return 1 if waiting list is not empty, 0 else
 */
int get_smallest_waiting_row(wl_t *waiting_global,
		ri_t *wl_idx, ri_t *wl_lp) ;

#if GBLA_WITH_FFLAS
/**
 * \brief Copies meta data from block matrix representation B to dense matrix
 * representation A. Needs modulus, too.
 *
 * \param dense matrix representation A
 *
 * \param block multiline matrix representation B
 *
 * \param characteristic of underlying field modulus
 */
void copyMetaData(DNS *A, sbm_fl_t *B, mod_t modulus) {
	A->row  = (dimen_t)B->nrows;
	A->col  = (dimen_t)B->ncols;
	A->ld   = ALIGN(A->col);
	A->mod  = (elemt_t)modulus;
	A->nnz  = (index_t)B->nnz;
}
#endif

/**
 * \brief Reduces dense block block_B with rectangular sparse block block_A.
 *
 * \param sparse block block_A
 *
 * \param sparse block block_B
 *
 * \param dense block dense_B
 *
 * \param block height bheight
 *
 * \param invert scalars? inv_scalars
 *
 * \param characteristic of underlying field modulus
 */
void red_with_rectangular_block(mbl_t *block_A, mbl_t *block_B, re_l_t **dense_B,
		const ri_t bheight, const int inv_scalars, const mod_t modulus);

/**
 * \brief Reduces dense block block_B with triangular sparse block block_A.
 *
 * \param sparse block block_A
 *
 * \param dense block dense_B
 *
 * \param block height bheight
 *
 * \param characteristic of underlying field modulus
 *
 * \param invert scalars? inv_scalars
 */
void red_with_triangular_block(mbl_t *block_A, re_l_t **dense_B,
		const ri_t bheight, int inv_scalars, const mod_t modulus);

/**
 * \brief Elimination procedure which reduces the block submatrix A to the unit
 * matrix. Corresponding changes in block submatrix B are carried out, too.
 *
 * \param block submatrix A (left upper side)
 *
 * \param block submatrix B (right upper side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_block(sbm_fl_t **A, sbm_fl_t *B, const mod_t modulus, int nthrds);

/**
 * \brief Elimination procedure which reduces the sparse block submatrix A to
 * the unit matrix. Corresponding changes in dense block submatrix B are
 * carried out, too.
 *
 * \param sparse block submatrix A (left upper side)
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_sparse_dense_block(sb_fl_t **A, dbm_fl_t *B, const mod_t modulus,
    int nthrds);

/**
 * \brief Different block tasks when reducing sparse block submatrix A.
 *
 * \param sparse block submatrix A (left upper side)
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param column index of blocks in B block_col_idx_B
 *
 * \param number of block rows in A nbrows_A
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_sparse_dense_blocks_task(const sb_fl_t *A, dbm_fl_t *B,
    const ci_t block_col_idx_B, const ri_t nbrows_A, const mod_t modulus);

/**
 * \brief Elimination procedure which reduces the hybrid block submatrix A to
 * the unit matrix. Corresponding changes in dense block submatrix B are
 * carried out, too.
 *
 * \param hybrid block submatrix A (left upper side)
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_hybrid_dense_block(hbm_fl_t **A, dbm_fl_t *B, const mod_t modulus,
    int nthrds);

/**
 * \brief Different block tasks when reducing hybrid block submatrix A.
 *
 * \param hybrid block submatrix A (left upper side)
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param column index of blocks in B block_col_idx_B
 *
 * \param number of block rows in A nbrows_A
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_hybrid_dense_blocks_task(hbm_fl_t *A, dbm_fl_t *B,
    const ci_t block_col_idx_B, const ri_t nbrows_A, const mod_t modulus);

/**
 * \brief Elimination procedure which reduces the hybrid block submatrix A to
 * the unit matrix. Corresponding changes in dense block submatrix B are
 * carried out, too.
 *
 * \param hybrid block submatrix A (left upper side)
 *
 * \param hybrid block submatrix B (right upper side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_hybrid_block(hbm_fl_t **A, hbm_fl_t *B, const mod_t modulus, int nthrds);

/**
 * \brief Different block tasks when reducing hybrid block submatrix A.
 *
 * \param hybrid block submatrix A (left upper side)
 *
 * \param hybrid block submatrix B (right upper side)
 *
 * \param column index of blocks in B block_col_idx_B
 *
 * \param number of block rows in A nbrows_A
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_hybrid_blocks_task(hbm_fl_t *A, hbm_fl_t *B,
    const ci_t block_col_idx_B, const ri_t nbrows_A, const mod_t modulus);

/**
 * \brief Elimination procedure which reduces the dense block submatrix A to
 * the unit matrix. Corresponding changes in dense block submatrix B are
 * carried out, too.
 *
 * \param dense block submatrix A (left upper side)
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_dense_block(dbm_fl_t **A, dbm_fl_t *B, const mod_t modulus, int nthrds);

/**
 * \brief Different block tasks when reducing dense block submatrix A.
 *
 * \param dense block submatrix A (left upper side)
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param column index of blocks in B block_col_idx_B
 *
 * \param number of block rows in A nbrows_A
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_dense_blocks_task(const dbm_fl_t *A, dbm_fl_t *B,
    const ci_t block_col_idx_B, const ri_t nbrows_A, const mod_t modulus);

/**
 * \brief Different block tasks when reducing block submatrix A.
 *
 * \param block submatrix A (left upper side)
 *
 * \param block submatrix B (right upper side)
 *
 * \param number of block rows in A nbrows_A
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_A_blocks_task(sbm_fl_t *A, sbm_fl_t *B,
		const ri_t nbrows_A, const mod_t modulus);

/**
 * \brief Elimination procedure which reduces the intermediate block submatrix C to zero.
 * Corresponding changes in dense block submatrix D are carried out using sparse block
 * submatrix B in column style, too.
 *
 * \param sparse block submatrix B (right upper side)
 *
 * \param intermediate block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_intermediate_block(sb_fl_t *B, ibm_fl_t **C, dbm_fl_t *D,
    const mod_t modulus, const int nthrds);

/**
 * \brief Elimination procedure which reduces the dense block submatrix C to zero.
 * Corresponding changes in dense block submatrix D are carried out using sparse block
 * submatrix B in column style, too.
 *
 * \param sparse block submatrix B (right upper side)
 *
 * \param dense block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_dense_sparse_block(sb_fl_t *B, dbm_fl_t **C, dbm_fl_t *D,
    const mod_t modulus, const int nthrds);

/**
 * \brief Elimination procedure which reduces the sparse block submatrix C to zero.
 * Corresponding changes in dense block submatrix D are carried out using dense block
 * submatrixB, too.
 *
 * \param sparse block submatrix B (right upper side)
 *
 * \param sparse block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_sparse_block(sb_fl_t *B, sb_fl_t **C, dbm_fl_t *D,
    const mod_t modulus, const int nthrds);

/**
 * \brief Elimination procedure which reduces the dense block submatrix C to zero.
 * Corresponding changes in dense block submatrix D are carried out using dense block
 * submatrixB, too.
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param sparse block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_block(dbm_fl_t *B, sb_fl_t **C, dbm_fl_t *D,
    const mod_t modulus, const int nthrds);

/**
 * \brief Different tasks of elimination of C operating on ten five of C only.
 * Updating all row entries via corresponding multiples from A.
 *
 * \param sparse submatrix C (left lower side)
 *
 * \param sparse submatrix A (left upper side)
 *
 * \param first row index in C idx1
 *
 * \param second row index in C idx2
 *
 * \param second row index in C idx3
 *
 * \param second row index in C idx4
 *
 * \param second row index in C idx5
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_keep_A_tasks_five(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx1, const ri_t idx2, const ri_t idx3, const ri_t idx4,
    const ri_t idx5, const mod_t modulus);

/**
 * \brief Different tasks of elimination of C operating on ten six of C only.
 * Updating all row entries via corresponding multiples from A.
 *
 * \param sparse submatrix C (left lower side)
 *
 * \param sparse submatrix A (left upper side)
 *
 * \param first row index in C idx1
 *
 * \param second row index in C idx2
 *
 * \param second row index in C idx3
 *
 * \param second row index in C idx4
 *
 * \param second row index in C idx5
 *
 * \param second row index in C idx6
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_keep_A_tasks_six(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx1, const ri_t idx2, const ri_t idx3, const ri_t idx4,
    const ri_t idx5, const ri_t idx6, const mod_t modulus);

/**
 * \brief Different tasks of elimination of C operating on eight rows of C only.
 * Updating all row entries via corresponding multiples from A.
 *
 * \param sparse submatrix C (left lower side)
 *
 * \param sparse submatrix A (left upper side)
 *
 * \param first row index in C idx1
 *
 * \param second row index in C idx2
 *
 * \param second row index in C idx3
 *
 * \param second row index in C idx4
 *
 * \param second row index in C idx5
 *
 * \param second row index in C idx6
 *
 * \param second row index in C idx7
 *
 * \param second row index in C idx8
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_keep_A_tasks_eight(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx1, const ri_t idx2, const ri_t idx3, const ri_t idx4,
    const ri_t idx5, const ri_t idx6, const ri_t idx7, const ri_t idx8,
    const mod_t modulus);

/**
 * \brief Different tasks of elimination of C operating on ten rows of C only.
 * Updating all row entries via corresponding multiples from A.
 *
 * \param sparse submatrix C (left lower side)
 *
 * \param sparse submatrix A (left upper side)
 *
 * \param first row index in C idx1
 *
 * \param second row index in C idx2
 *
 * \param second row index in C idx3
 *
 * \param second row index in C idx4
 *
 * \param second row index in C idx5
 *
 * \param second row index in C idx6
 *
 * \param second row index in C idx7
 *
 * \param second row index in C idx8
 *
 * \param second row index in C idx9
 *
 * \param second row index in C idx10
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_keep_A_tasks_ten(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx1, const ri_t idx2, const ri_t idx3, const ri_t idx4,
    const ri_t idx5, const ri_t idx6, const ri_t idx7, const ri_t idx8,
    const ri_t idx9, const ri_t idx10, const mod_t modulus);

/**
 * \brief Different tasks of elimination of C operating on four rows of C only.
 * Updating all row entries via corresponding multiples from A.
 *
 * \param sparse submatrix C (left lower side)
 *
 * \param sparse submatrix A (left upper side)
 *
 * \param first row index in C idx1
 *
 * \param second row index in C idx2
 *
 * \param second row index in C idx3
 *
 * \param second row index in C idx4
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_keep_A_tasks_four(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx1, const ri_t idx2, const ri_t idx3, const ri_t idx4,
    const mod_t modulus);

/**
 * \brief Different tasks of elimination of C operating on three rows of C only.
 * Updating all row entries via corresponding multiples from A.
 *
 * \param sparse submatrix C (left lower side)
 *
 * \param sparse submatrix A (left upper side)
 *
 * \param first row index in C idx1
 *
 * \param second row index in C idx2
 *
 * \param second row index in C idx3
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_keep_A_tasks_three(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx1, const ri_t idx2, const ri_t idx3, const mod_t modulus);

/**
 * \brief Different tasks of elimination of C operating on two rows of C only.
 * Updating all row entries via corresponding multiples from A.
 *
 * \param sparse submatrix C (left lower side)
 *
 * \param sparse submatrix A (left upper side)
 *
 * \param first row index in C idx1
 *
 * \param second row index in C idx2
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_keep_A_tasks_double(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx1, const ri_t idx2, const mod_t modulus);

/**
 * \brief Different tasks of elimination of C operating on one row of C only.
 * Updating all row entries via corresponding multiples from A.
 *
 * \param sparse submatrix C (left lower side)
 *
 * \param sparse submatrix A (left upper side)
 *
 * \param row index in C idx
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_keep_A_tasks_single(sm_fl_t *C, const sm_fl_t *A,
    const ri_t idx, const mod_t modulus);

/**
 * \brief Elimination procedure which reduces the sparse submatrix C to zero.
 * For this an intermediate dense wide representation is used. Multiplies
 * corresponding elements from A and prepares C for updating D with B later on.
 *
 * \param sparse submatrix C (left lower side)
 *
 * \param sparse submatrix A (left upper side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_keep_A(sm_fl_t *C, sm_fl_t **A, const mod_t modulus,
    const int nthrds);

/**
 * \brief Different block tasks when reducing intermediate block submatrix C.
 *
 * \param sparse block submatrix in column style B (right upper side)
 *
 * \param intermediate block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param column index of blocks in D block_col_idx_D
 *
 * \param number of block rows in C nblock_rows_C
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_intermediate_blocks_task(sb_fl_t *B, ibm_fl_t *C, dbm_fl_t *D,
    const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
    const mod_t modulus);

/**
 * \brief Different block tasks when reducing dense block submatrix C.
 *
 * \param sparse block submatrix in column style B (right upper side)
 *
 * \param dense block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param column index of blocks in D block_col_idx_D
 *
 * \param number of block rows in C nblock_rows_C
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_dense_sparse_blocks_task(sb_fl_t *B, dbm_fl_t *C, dbm_fl_t *D,
    const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
    const mod_t modulus);

/**
 * \brief Different block tasks when reducing sparse block submatrix C.
 *
 * \param sparse block submatrix B (right upper side)
 *
 * \param sparse block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param column index of blocks in D block_col_idx_D
 *
 * \param number of block rows in C nblock_rows_C
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_sparse_blocks_task(sb_fl_t *B, sb_fl_t *C, dbm_fl_t *D,
    const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
    const mod_t modulus);

/**
 * \brief Different block tasks when reducing dense block submatrix C.
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param sparse block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param column index of blocks in D block_col_idx_D
 *
 * \param number of block rows in C nblock_rows_C
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_sparse_dense_blocks_task(dbm_fl_t *B, sb_fl_t *C, dbm_fl_t *D,
    const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
    const mod_t modulus);

/**
 * \brief Elimination procedure which reduces the dense block submatrix C to zero.
 * Corresponding changes in dense block submatrix D are carried out using dense block
 * submatrixB, too.
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param dense block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_dense_block(dbm_fl_t *B, dbm_fl_t **C, dbm_fl_t *D,
    const mod_t modulus, const int nthrds);

/**
 * \brief Different block tasks when reducing denes block submatrix C.
 *
 * \param dense block submatrix B (right upper side)
 *
 * \param dense block submatrix C (left lower side)
 *
 * \param dense block submatrix D (right lower side)
 *
 * \param column index of blocks in D block_col_idx_D
 *
 * \param number of block rows in C nblock_rows_C
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_dense_blocks_task(dbm_fl_t *B, dbm_fl_t *C, dbm_fl_t *D,
    const ci_t block_col_idx_D, const ri_t nbrows_C, const ci_t nbcols_C,
    const mod_t modulus);

/**
 * \brief Elimination procedure which reduces the block submatrix C to zero.
 * Corresponding changes in block submatrix D are carried out using B, too.
 *
 * \param block submatrix B (right upper side)
 *
 * \param block submatrix C (left lower side)
 *
 * \param block submatrix D (right lower side)
 *
 * \param inverse scalars? inv_scalars
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_block(sbm_fl_t *B, sbm_fl_t **C, sbm_fl_t *D, const int inv_scalars,
		const mod_t modulus, const int nthrds);

/**
 * \brief Different block tasks when reducing block submatrix C.
 *
 * \param block submatrix B (right upper side)
 *
 * \param block submatrix C (left lower side)
 *
 * \param block submatrix D (right lower side)
 *
 * \param number of block rows in C nblock_rows_C
 *
 * \param inverse scalars? inv_scalars
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_blocks_task(sbm_fl_t *B, sbm_fl_t *C, sbm_fl_t *D,
		const ri_t nbrows_C, const ci_t nbcols_C,
		const int inv_scalars, const mod_t modulus);

/**
 * \brief Elimination procedure which reduces the block submatrix D to an
 * upper triangular matrix. For this FFLAS-FFPACK procedures are used. Thus the
 * input matrix D has to be converted to DNS type first.
 *
 * \param block submatrix D (right lower side), input matrix
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return rank of D_red
 */
ri_t elim_fl_D_fflas_ffpack(sbm_fl_t *D, const mod_t modulus, int nthrds);

/**
 * \brief Elimination procedure which reduces the block submatrix D to an
 * upper triangular matrix. Note that the input matrix D will be removed
 * later on and the output matrix D_red is in multiline format for further
 * reduction steps.
 *
 * \param block submatrix D (right lower side), input matrix
 *
 * \param multiline submatrix D_red, output matrix
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return rank of D_red
 */
ri_t elim_fl_D_block(sbm_fl_t *D, sm_fl_ml_t *D_red, const mod_t modulus, int nthrds);

/**
 * \brief Reduces dense row modulo prime M->mod
 *
 * \param dense row submatrix A
 *
 * \param row index ridx
 */
static inline void red_mod_dense_row(dm_t *A, const ri_t ridx)
{
  ci_t i;
  for (i=A->row[ridx]->lead+1; i<A->ncols; ++i) {
    /* A->row[ridx]->val[i]  = MODP(A->row[ridx]->val[i], A->mod); */
    A->row[ridx]->val[i]  = (re_t)(A->row[ridx]->val[i] % A->mod);
  }
}

/**
 * \brief Normalizes dense row with index ridx in dense row matrix D
 *
 * \param dense row submatrix D
 *
 * \param row index ridx
 */
static inline void normalize_dense_row(dm_t *A, const ri_t ridx)
{
  ci_t i;
  const ci_t lead_idx  = A->row[ridx]->lead;
  /* re_t inv =  (re_t)MODP(A->row[ridx]->val[lead_idx], A->mod); */
  re_t inv =  (re_t)(A->row[ridx]->val[lead_idx] % A->mod);
  //printf("normalize %u\n", inv);
  A->row[ridx]->val[lead_idx]  = (re_l_t)1;
  if (inv == 1) {
    for (i=lead_idx+1; i<A->ncols; ++i) {
      /* A->row[ridx]->val[i] =   MODP(A->row[ridx]->val[i], A->mod); */
      A->row[ridx]->val[i] =   (re_t)(A->row[ridx]->val[i] % A->mod);
    }
  } else {
    inverse_val(&inv, A->mod);
    //printf("inv normalize %u\n", inv);
    register const re_t cinv  = inv;
    for (i=lead_idx+1; i<A->ncols; ++i) {
#if defined(GBLA_USE_UINT16) || defined(GBLA_USE_UINT32)
#else
      /* A->row[ridx]->val[i] =   MODP(A->row[ridx]->val[i], A->mod); */
      A->row[ridx]->val[i] =   (re_t)(A->row[ridx]->val[i] % A->mod);
#endif
      A->row[ridx]->val[i] *=  cinv;
      /* A->row[ridx]->val[i] =   MODP(A->row[ridx]->val[i], A->mod); */
      A->row[ridx]->val[i] =   (re_t)(A->row[ridx]->val[i] % A->mod);
    }
  }
}

/**
 * \brief Normalizes row of input matrix from Schreyer (Prym-Green
 * conjecture)
 *
 * \param sparse input matrix M
 */
static inline void normalize_schreyer_input_rows(sm_t *M)
{
  ri_t i;
  ci_t j;
  re_t inv;
  re_l_t tmp;
  for (i=0; i<M->nrows; ++i) {
    inv = (re_t)(M->rows[i][0] % M->mod);
    M->rows[i][0] = (re_t)1;
    inverse_val(&inv, M->mod);
    register const re_t cinv  = inv;
    for (j=1; j<M->rwidth[i]; ++j) {
      tmp = (re_l_t)M->rows[i][j];
      tmp *=  cinv;
      M->rows[i][j] =   (re_t)(tmp % M->mod);
    }
  }
}

/**
 * \brief Acronym for normalization of dense row with index curr_row_to_reduce
 * in dense row matrix D
 *
 * \param dense row matrix D
 *
 * \param index of row to be normalized curr_row_to_reduce
 *
 * \param index of pivot row
 */
static inline void save_pivot(dm_t *D, const ri_t curr_row_to_reduce, const ri_t new_piv_idx)
{
  ci_t i;
  //D->row[curr_row_to_reduce]->piv_val  = (re_t *)calloc(D->ncols, sizeof(re_t));
  normalize_dense_row(D, curr_row_to_reduce);
  if (D->row[curr_row_to_reduce]->lead < D->ncols) {
    D->row[new_piv_idx]->piv_val  = (re_t *)malloc(D->ncols * sizeof(re_t));
    // set all elements before lead to zero
    memset(D->row[new_piv_idx]->piv_val, 0, D->row[curr_row_to_reduce]->lead * sizeof(re_t));
    for (i=D->row[curr_row_to_reduce]->lead; i<D->ncols; ++i)
      D->row[new_piv_idx]->piv_val[i] = (re_t)D->row[curr_row_to_reduce]->val[i];
  }
  D->row[new_piv_idx]->piv_lead = D->row[curr_row_to_reduce]->lead;
  free(D->row[curr_row_to_reduce]->val);
  D->row[curr_row_to_reduce]->val  = NULL;
  //printf("NEW PIVOT %u -- lead %u\n",curr_row_to_reduce, D->row[curr_row_to_reduce]->lead);
  //for (int ii=0; ii<D->ncols; ++ii)
  //  printf("%lu ",D->row[curr_row_to_reduce]->val[ii]);
  //printf("\n");
}

/**
 * \brief Saves pivot values of dense row with index curr_row_to_reduce
 * in dense row matrix D
 *
 * \note Does not normalize the row.
 *
 * \param dense row matrix D
 *
 * \param index of row to be normalized curr_row_to_reduce
 *
 * \param index of pivot row
 */
static inline void save_pivot_without_normalization(dm_t *D, const ri_t curr_row_to_reduce, const ri_t new_piv_idx)
{
  ci_t i;
  //D->row[curr_row_to_reduce]->piv_val  = (re_t *)calloc(D->ncols, sizeof(re_t));
  // we only reduce modulo field charactistic, but do not normalize the row
  red_mod_dense_row(D, curr_row_to_reduce);
  if (D->row[curr_row_to_reduce]->lead < D->ncols) {
    D->row[new_piv_idx]->piv_val  = (re_t *)malloc(D->ncols * sizeof(re_t));
    // set all elements before lead to zero
    memset(D->row[new_piv_idx]->piv_val, 0, D->row[curr_row_to_reduce]->lead * sizeof(re_t));
    for (i=D->row[curr_row_to_reduce]->lead; i<D->ncols; ++i)
      D->row[new_piv_idx]->piv_val[i] = (re_t)D->row[curr_row_to_reduce]->val[i];
  }
  D->row[new_piv_idx]->piv_lead = D->row[curr_row_to_reduce]->lead;
  free(D->row[curr_row_to_reduce]->val);
  D->row[curr_row_to_reduce]->val  = NULL;
  //printf("NEW PIVOT %u -- lead %u\n",curr_row_to_reduce, D->row[curr_row_to_reduce]->lead);
  //for (int ii=0; ii<D->ncols; ++ii)
  //  printf("%lu ",D->row[curr_row_to_reduce]->val[ii]);
  //printf("\n");
}


/**
 * \brief Reduces dense row ri via row rj in D using precomputed inverted
 * multiplier mult. Gives special starting column position from
 *
 * \param dense row submatrix D
 *
 * \param row index ri
 *
 * \param row index rj
 *
 * \param precomputed inverted multiplier mult
 *
 * \param starting column position from
 *
 */
static inline void reduce_dense_row_from(dm_t *A, const ri_t ri, const ri_t rj,
    const re_t mult, const ci_t from)
{
  ci_t i;
  const re_t *reducers  = A->row[rj]->piv_val;
  
  i = from;
  
  //printf("i initially %u\n",i);
  // leading nonzero element has to become zero
  //assert(MODP(A->row[ri]->val[i-1] + mult * reducers[i-1], A->mod) == 0);
  //A->row[ri]->val[i-1]  = 0;

  if (A->ncols-i > 7) {
    for (; i<A->ncols-7; i=i+8) {
      A->row[ri]->val[i]    +=  (re_l_t)mult * reducers[i];
      A->row[ri]->val[i+1]    +=  (re_l_t)mult * reducers[i+1];
      A->row[ri]->val[i+2]    +=  (re_l_t)mult * reducers[i+2];
      A->row[ri]->val[i+3]    +=  (re_l_t)mult * reducers[i+3];
      A->row[ri]->val[i+4]    +=  (re_l_t)mult * reducers[i+4];
      A->row[ri]->val[i+5]    +=  (re_l_t)mult * reducers[i+5];
      A->row[ri]->val[i+6]    +=  (re_l_t)mult * reducers[i+6];
      A->row[ri]->val[i+7]    +=  (re_l_t)mult * reducers[i+7];
#if COUNT_REDS
      nreductions +=  8;
#endif
    }
  }
  if (A->ncols-i > 4) {
    A->row[ri]->val[i]    +=  (re_l_t)mult * reducers[i];
    A->row[ri]->val[i+1]    +=  (re_l_t)mult * reducers[i+1];
    A->row[ri]->val[i+2]    +=  (re_l_t)mult * reducers[i+2];
    A->row[ri]->val[i+3]    +=  (re_l_t)mult * reducers[i+3];
    i = i+4;
#if COUNT_REDS
    nreductions +=  4;
#endif
  }
  for (; i<A->ncols; ++i) {
    A->row[ri]->val[i]    +=  (re_l_t)mult * reducers[i];
#if COUNT_REDS
    nreductions +=  1;
#endif
  }
  // search new lead
  update_lead_of_row(A, ri);
  // if we get here then the row is completely zero
  if (A->row[ri]->lead == A->ncols) {
    free(A->row[ri]->val);
    A->row[ri]->val = NULL;
  }
}

/**
 * \brief Reduces dense row ri via rows rj1 and rj2 in D using precomputed inverted
 * multiplier mult
 *
 * \param dense row submatrix D
 *
 * \param row index ri
 *
 * \param row index rj1
 *
 * \param row index rj2
 *
 * \param precomputed inverted multiplier mult1 for pivot row rj
 *
 * \param precomputed inverted multiplier mult2 for pivot row rj2
 *
 */
static inline void reduce_dense_row_twice_from(dm_t *A, const ri_t ri, const ri_t rj1,
    const ri_t rj2, const re_t mult1, const re_t mult2, const ci_t from)
{
  ci_t i;
  const re_t *reducers1 = A->row[rj1]->piv_val;
  const re_t *reducers2 = A->row[rj2]->piv_val;
  
  i = from;
  
  if (A->ncols-i > 7) {
    for (; i<A->ncols-7; i=i+8) {
      A->row[ri]->val[i]    +=
        (re_l_t)mult1 * reducers1[i] + (re_l_t)mult2 * reducers2[i];
      A->row[ri]->val[i+1]    +=
        (re_l_t)mult1 * reducers1[i+1] + (re_l_t)mult2 * reducers2[i+1];
      A->row[ri]->val[i+2]    +=
        (re_l_t)mult1 * reducers1[i+2] + (re_l_t)mult2 * reducers2[i+2];
      A->row[ri]->val[i+3]    +=
        (re_l_t)mult1 * reducers1[i+3] + (re_l_t)mult2 * reducers2[i+3];
      A->row[ri]->val[i+4]    +=
        (re_l_t)mult1 * reducers1[i+4] + (re_l_t)mult2 * reducers2[i+4];
      A->row[ri]->val[i+5]    +=
        (re_l_t)mult1 * reducers1[i+5] + (re_l_t)mult2 * reducers2[i+5];
      A->row[ri]->val[i+6]    +=
        (re_l_t)mult1 * reducers1[i+6] + (re_l_t)mult2 * reducers2[i+6];
      A->row[ri]->val[i+7]    +=
        (re_l_t)mult1 * reducers1[i+7] + (re_l_t)mult2 * reducers2[i+7];
#if COUNT_REDS
      nreductions +=  8;
#endif
    }
  }
  if (A->ncols-i > 4) {
      A->row[ri]->val[i]    +=
        (re_l_t)mult1 * reducers1[i] + (re_l_t)mult2 * reducers2[i];
      A->row[ri]->val[i+1]    +=
        (re_l_t)mult1 * reducers1[i+1] + (re_l_t)mult2 * reducers2[i+1];
      A->row[ri]->val[i+2]    +=
        (re_l_t)mult1 * reducers1[i+2] + (re_l_t)mult2 * reducers2[i+2];
      A->row[ri]->val[i+3]    +=
        (re_l_t)mult1 * reducers1[i+3] + (re_l_t)mult2 * reducers2[i+3];
    i = i+4;
#if COUNT_REDS
    nreductions +=  4;
#endif
  }
  for (; i<A->ncols; ++i) {
      A->row[ri]->val[i]    +=
        (re_l_t)mult1 * reducers1[i] + (re_l_t)mult2 * reducers2[i];
#if COUNT_REDS
    nreductions +=  1;
#endif
  }
  // search new lead
  update_lead_of_row(A, ri);
  // if we get here then the row is completely zero
  if (A->row[ri]->lead == A->ncols) {
    free(A->row[ri]->val);
    A->row[ri]->val = NULL;
  }
}

/**
 * \brief Reduces dense row ri via rows rj1, rj2 and rj3 in D using precomputed inverted
 * multiplier mult
 *
 * \param dense row submatrix D
 *
 * \param row index ri
 *
 * \param row index rj1
 *
 * \param row index rj2
 *
 * \param row index rj3
 *
 * \param precomputed inverted multiplier mult1 for pivot row rj
 *
 * \param precomputed inverted multiplier mult2 for pivot row rj2
 *
 * \param precomputed inverted multiplier mult3 for pivot row rj3
 *
 */
static inline void reduce_dense_row_three(dm_t *A, const ri_t ri, const ri_t rj1,
    const ri_t rj2, const ri_t rj3, const re_t mult1, const re_t mult2,
    const re_t mult3)
{
  ci_t i;
  const re_t *reducers1 = A->row[rj1]->piv_val;
  const re_t *reducers2 = A->row[rj2]->piv_val;
  const re_t *reducers3 = A->row[rj3]->piv_val;
  
  i = A->row[rj3]->piv_lead + 1;
  //A->row[ri]->val[i-1]  = 0;
  
  //printf("i initially %u\n",i);
  // leading nonzero element has to become zero
  A->row[ri]->val[i-1]  = 0;
  if (A->ncols-i > 7) {
    for (; i<A->ncols-7; i=i+8) {
      A->row[ri]->val[i]    +=
        (re_l_t)mult1 * reducers1[i] + (re_l_t)mult2 * reducers2[i]
        + (re_l_t)mult3 * reducers3[i];
      A->row[ri]->val[i+1]    +=
        (re_l_t)mult1 * reducers1[i+1] + (re_l_t)mult2 * reducers2[i+1]
        + (re_l_t)mult3 * reducers3[i+1];
      A->row[ri]->val[i+2]    +=
        (re_l_t)mult1 * reducers1[i+2] + (re_l_t)mult2 * reducers2[i+2]
        + (re_l_t)mult3 * reducers3[i+2];
      A->row[ri]->val[i+3]    +=
        (re_l_t)mult1 * reducers1[i+3] + (re_l_t)mult2 * reducers2[i+3]
        + (re_l_t)mult3 * reducers3[i+3];
      A->row[ri]->val[i+4]    +=
        (re_l_t)mult1 * reducers1[i+4] + (re_l_t)mult2 * reducers2[i+4]
        + (re_l_t)mult3 * reducers3[i+4];
      A->row[ri]->val[i+5]    +=
        (re_l_t)mult1 * reducers1[i+5] + (re_l_t)mult2 * reducers2[i+5]
        + (re_l_t)mult3 * reducers3[i+5];
      A->row[ri]->val[i+6]    +=
        (re_l_t)mult1 * reducers1[i+6] + (re_l_t)mult2 * reducers2[i+6]
        + (re_l_t)mult3 * reducers3[i+6];
      A->row[ri]->val[i+7]    +=
        (re_l_t)mult1 * reducers1[i+7] + (re_l_t)mult2 * reducers2[i+7]
        + (re_l_t)mult3 * reducers3[i+7];
#if COUNT_REDS
      nreductions +=  8;
#endif
    }
  }
  if (A->ncols-i > 4) {
      A->row[ri]->val[i]    +=
        (re_l_t)mult1 * reducers1[i] + (re_l_t)mult2 * reducers2[i]
        + (re_l_t)mult3 * reducers3[i];
      A->row[ri]->val[i+1]    +=
        (re_l_t)mult1 * reducers1[i+1] + (re_l_t)mult2 * reducers2[i+1]
        + (re_l_t)mult3 * reducers3[i+1];
      A->row[ri]->val[i+2]    +=
        (re_l_t)mult1 * reducers1[i+2] + (re_l_t)mult2 * reducers2[i+2]
        + (re_l_t)mult3 * reducers3[i+2];
      A->row[ri]->val[i+3]    +=
        (re_l_t)mult1 * reducers1[i+3] + (re_l_t)mult2 * reducers2[i+3]
        + (re_l_t)mult3 * reducers3[i+3];
    i = i+4;
#if COUNT_REDS
    nreductions +=  4;
#endif
  }
  for (; i<A->ncols; ++i) {
      A->row[ri]->val[i]    +=
        (re_l_t)mult1 * reducers1[i] + (re_l_t)mult2 * reducers2[i]
        + (re_l_t)mult3 * reducers3[i];
#if COUNT_REDS
    nreductions +=  1;
#endif
  }
  // search new lead
  update_lead_of_row(A, ri);
  // if we get here then the row is completely zero
  if (A->row[ri]->lead == A->ncols) {
    free(A->row[ri]->val);
    A->row[ri]->val = NULL;
  }
}

/**
 * \brief Reduces dense row ri via rows rj1 and rj2 in D using precomputed inverted
 * multiplier mult
 *
 * \param dense row submatrix D
 *
 * \param row index ri
 *
 * \param row index rj1
 *
 * \param row index rj2
 *
 * \param precomputed inverted multiplier mult1 for pivot row rj
 *
 * \param precomputed inverted multiplier mult2 for pivot row rj2
 *
 */
static inline void reduce_dense_row_twice(dm_t *A, const ri_t ri, const ri_t rj1,
    const ri_t rj2, const re_t mult1, const re_t mult2)
{
  ci_t i;
  const re_t *reducers1 = A->row[rj1]->piv_val;
  const re_t *reducers2 = A->row[rj2]->piv_val;
  
  i = A->row[rj2]->piv_lead + 1;
  
  //printf("i initially %u\n",i);
  // leading nonzero element has to become zero
  //assert(MODP(A->row[ri]->val[i-1] + mult * reducers[i-1], A->mod) == 0);
  //A->row[ri]->val[i-1]  = 0;

  if (A->ncols-i > 7) {
    for (; i<A->ncols-7; i=i+8) {
      A->row[ri]->val[i]    +=
        (re_l_t)mult1 * reducers1[i] + (re_l_t)mult2 * reducers2[i];
      A->row[ri]->val[i+1]    +=
        (re_l_t)mult1 * reducers1[i+1] + (re_l_t)mult2 * reducers2[i+1];
      A->row[ri]->val[i+2]    +=
        (re_l_t)mult1 * reducers1[i+2] + (re_l_t)mult2 * reducers2[i+2];
      A->row[ri]->val[i+3]    +=
        (re_l_t)mult1 * reducers1[i+3] + (re_l_t)mult2 * reducers2[i+3];
      A->row[ri]->val[i+4]    +=
        (re_l_t)mult1 * reducers1[i+4] + (re_l_t)mult2 * reducers2[i+4];
      A->row[ri]->val[i+5]    +=
        (re_l_t)mult1 * reducers1[i+5] + (re_l_t)mult2 * reducers2[i+5];
      A->row[ri]->val[i+6]    +=
        (re_l_t)mult1 * reducers1[i+6] + (re_l_t)mult2 * reducers2[i+6];
      A->row[ri]->val[i+7]    +=
        (re_l_t)mult1 * reducers1[i+7] + (re_l_t)mult2 * reducers2[i+7];
#if COUNT_REDS
      nreductions +=  8;
#endif
    }
  }
  if (A->ncols-i > 4) {
      A->row[ri]->val[i]    +=
        (re_l_t)mult1 * reducers1[i] + (re_l_t)mult2 * reducers2[i];
      A->row[ri]->val[i+1]    +=
        (re_l_t)mult1 * reducers1[i+1] + (re_l_t)mult2 * reducers2[i+1];
      A->row[ri]->val[i+2]    +=
        (re_l_t)mult1 * reducers1[i+2] + (re_l_t)mult2 * reducers2[i+2];
      A->row[ri]->val[i+3]    +=
        (re_l_t)mult1 * reducers1[i+3] + (re_l_t)mult2 * reducers2[i+3];
    i = i+4;
#if COUNT_REDS
    nreductions +=  4;
#endif
  }
  for (; i<A->ncols; ++i) {
    A->row[ri]->val[i]    +=
      (re_l_t)mult1 * reducers1[i] + (re_l_t)mult2 * reducers2[i];
#if COUNT_REDS
    nreductions +=  1;
#endif
  }
  // search new lead
  update_lead_of_row(A, ri);
  // if we get here then the row is completely zero
  if (A->row[ri]->lead == A->ncols) {
    free(A->row[ri]->val);
    A->row[ri]->val = NULL;
  }
}

/**
 * \brief Reduces dense row ri via row rj in D using precomputed inverted
 * multiplier mult
 *
 * \param dense row submatrix D
 *
 * \param row index ri
 *
 * \param row index rj
 *
 * \param precomputed inverted multiplier mult
 *
 */
static inline void reduce_dense_row(dm_t *A, const ri_t ri, const ri_t rj, const re_t mult)
{
  ci_t i;
  const re_t *reducers  = A->row[rj]->piv_val;
  
  i = A->row[rj]->piv_lead + 1;
  
  //printf("i initially %u\n",i);
  // leading nonzero element has to become zero
  //A->row[ri]->val[i-1]  = 0;


  if (A->ncols-i > 7) {
    for (; i<A->ncols-7; i=i+8) {
      A->row[ri]->val[i]    +=  (re_l_t)mult * reducers[i];
      A->row[ri]->val[i+1]    +=  (re_l_t)mult * reducers[i+1];
      A->row[ri]->val[i+2]    +=  (re_l_t)mult * reducers[i+2];
      A->row[ri]->val[i+3]    +=  (re_l_t)mult * reducers[i+3];
      A->row[ri]->val[i+4]    +=  (re_l_t)mult * reducers[i+4];
      A->row[ri]->val[i+5]    +=  (re_l_t)mult * reducers[i+5];
      A->row[ri]->val[i+6]    +=  (re_l_t)mult * reducers[i+6];
      A->row[ri]->val[i+7]    +=  (re_l_t)mult * reducers[i+7];
#if COUNT_REDS
      nreductions +=  8;
#endif
    }
  }
  if (A->ncols-i > 4) {
    A->row[ri]->val[i]    +=  (re_l_t)mult * reducers[i];
    A->row[ri]->val[i+1]    +=  (re_l_t)mult * reducers[i+1];
    A->row[ri]->val[i+2]    +=  (re_l_t)mult * reducers[i+2];
    A->row[ri]->val[i+3]    +=  (re_l_t)mult * reducers[i+3];
    i = i+4;
#if COUNT_REDS
    nreductions +=  4;
#endif
  }
  for (; i<A->ncols; ++i) {
    A->row[ri]->val[i]    +=  (re_l_t)mult * reducers[i];
#if COUNT_REDS
    nreductions +=  1;
#endif
  }
  // search new lead
  update_lead_of_row(A, ri);
  // if we get here then the row is completely zero
  if (A->row[ri]->lead == A->ncols) {
    free(A->row[ri]->val);
    A->row[ri]->val = NULL;
  }
}

/**
 * \brief Reduces dense row ri from matrix B via row rj from matrix D using precomputed inverted
 * multiplier mult
 *
 * \param dense row submatrix B
 *
 * \param dense row submatrix D
 *
 * \param row index ri
 *
 * \param row index rj
 *
 * \param precomputed inverted multiplier mult
 *
 */
static inline void reduce_dense_row_of_B_by_D(dm_t *B, const dm_t *D, const ri_t ri,
    const ri_t rj, const re_t mult)
{
  ci_t i;
  const re_t *reducers  = D->row[rj]->piv_val;
  
  i = B->row[ri]->lead + 1;
  
  //printf("i initially %u / %u for %u and %u\n",i, B->ncols, ri, rj);
  // leading nonzero element has to become zero
  //printf("multiplier %u \n",mult);

  if (B->ncols-i > 7) {
    for (; i<B->ncols-7; i=i+8) {
      B->row[ri]->val[i]    +=  (re_l_t)mult * reducers[i];
      B->row[ri]->val[i+1]    +=  (re_l_t)mult * reducers[i+1];
      B->row[ri]->val[i+2]    +=  (re_l_t)mult * reducers[i+2];
      B->row[ri]->val[i+3]    +=  (re_l_t)mult * reducers[i+3];
      B->row[ri]->val[i+4]    +=  (re_l_t)mult * reducers[i+4];
      B->row[ri]->val[i+5]    +=  (re_l_t)mult * reducers[i+5];
      B->row[ri]->val[i+6]    +=  (re_l_t)mult * reducers[i+6];
      B->row[ri]->val[i+7]    +=  (re_l_t)mult * reducers[i+7];
#if COUNT_REDS
      nreductions +=  8;
#endif
    }
  }
  if (B->ncols-i > 4) {
    B->row[ri]->val[i]    +=  (re_l_t)mult * reducers[i];
    B->row[ri]->val[i+1]    +=  (re_l_t)mult * reducers[i+1];
    B->row[ri]->val[i+2]    +=  (re_l_t)mult * reducers[i+2];
    B->row[ri]->val[i+3]    +=  (re_l_t)mult * reducers[i+3];
    i = i+4;
#if COUNT_REDS
    nreductions +=  4;
#endif
  }
  for (; i<B->ncols; ++i) {
    B->row[ri]->val[i]    +=  (re_l_t)mult * reducers[i];
    //printf("%lu | ",B->row[ri]->val[i]);
#if COUNT_REDS
    nreductions +=  1;
#endif
  }
  //printf("\n");
  // search new lead
  update_lead_of_row(B, ri);
  // if we get here then the row is completely zero
  if (B->row[ri]->lead == B->ncols) {
    free(B->row[ri]->val);
    B->row[ri]->val = NULL;
  }
}

/**
 * \brief Does a sequential pre elimination of D in order to start a parallel
 * version with known pivots later on. It completely reduces D, i.e. if nthrds =
 * 1 the resulting D is completely reduced.
 *
 * \param dense row submatrix D
 *
 * \param global last pivot up to which to reduce sequentially global_last_piv
 */
void pre_elim_sequential_completely(dm_t *D, const ri_t global_last_piv);

/**
 * \brief Does a sequential pre elimination of D in order to start a parallel
 * version with known pivots later on.
 *
 * \param dense row submatrix D
 *
 * \param global last pivot up to which to reduce sequentially global_last_piv
 *
 * \param number of threads nthrds
 */
void pre_elim_sequential(dm_t *D, const ri_t global_last_piv, const int nthrds);

/**
 * \brief Does a sequential pre elimination of D in order to start a parallel
 * version with known pivots later on.
 *
 * \param dense row submatrix D
 *
 * \param global last pivot up to which to reduce sequentially global_last_piv
 *
 * \param number of threads nthrds
 */
void pre_elim_sequential_test(dm_t *D, const ri_t global_last_piv, const int nthrds);

/**
 * \brief After we have reduced a subset of rows of D we reduce them upwards
 * from index last_row to 0
 *
 * \param dense row submatrix D
 *
 * \param row up to which we have reduced D already
 */
void reduce_upwards(dm_t *D, const ri_t last_row);

/**
 * \brief Elimination procedure which reduces the dense row submatrix D to an
 * upper triangular matrix. Uses a structured Gaussian Elimination. Returns
 * the rank of D.
 *
 * \note Assumes D->nrows > 0.
 *
 * \param dense row submatrix D
 *
 * \param number of threads nthrds
 *
 * \return rank of D
 */
ri_t elim_fl_dense_D_test(dm_t *D, const int nthreads);

/**
 * \brief Elimination procedure which reduces the dense row submatrix D
 * to a completely reduced upper triangular matrix. Uses a structured
 * Gaussian Elimination. Returns the rank of D.
 *
 * \note Assumes D->nrows > 0.
 *
 * \param dense row submatrix D
 *
 * \param number of threads nthrds
 *
 * \return rank of D
 */
ri_t elim_fl_dense_D_completely(dm_t *D, const int nthreads);

/**
 * \brief Elimination procedure which reduces the dense row submatrix D to an
 * upper triangular matrix. Uses a structured Gaussian Elimination. Returns
 * the rank of D.
 *
 * \note Assumes D->nrows > 0.
 *
 * \param dense row submatrix D
 *
 * \param number of threads nthrds
 *
 * \return rank of D
 */
ri_t elim_fl_dense_D(dm_t *D, const int nthreads);

/**
 * \brief Parallel tasks doing the complete/reduced Gaussian Elimination of D.
 *
 * \param dense row submatrix D
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_dense_D_completely_tasks(dm_t *D);

/**
 * \brief Parallel tasks doing the Gaussian Elimination of D.
 *
 * \param dense row submatrix D
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_dense_D_tasks(dm_t *D);

/**
 * \brief Elimination procedure which reduces the multiline matrix C to zero
 * carrying out corresponding computations by A and B
 *
 * \param multiline submatrix C (left lower side)
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param characteristic of underlying field modulus
 *
 * \param number of threads nthrds
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_ml(sm_fl_ml_t *C, sm_fl_ml_t *A, mod_t modulus, int nthrds);

/**
 * \brief Executing corresponding multiples from A in multilines of C in order
 * to prepare C for the reduction of D by multiples of B
 *
 * \param multiline submatrix C (left lower side)
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param characteristic of underlying field modulus
 *
 * \return 0 if success, 1 if failure
 */
int elim_fl_C_ml_task(sm_fl_ml_t *C, sm_fl_ml_t *A, ri_t row_idx, mod_t modulus);

/**
 * \brief Echelonizes the multiline rows of A from row index from up to row
 * index to sequential. This is done in order to prepare a basis for the
 * upcoming parallel structure echelonization
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param start row index from
 *
 * \param end row index to
 *
 * \param field characteristic modulus
 *
 * \return number of real pivots found
 */
ri_t echelonize_rows_sequential(sm_fl_ml_t *A, const ri_t from, const ri_t to,
		const mod_t modulus);

/**
 * \brief Echelonizes the multiline rows of A from row index from up to row
 * index to sequential. This is done in order to prepare a basis for the
 * upcoming parallel structure echelonization
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param field characteristic modulus
 *
 * \return 0 if success, 1 if failure
 */
int echelonize_rows_task(sm_fl_ml_t *A, const ri_t ml_ncols,
		/* ri_t global_next_row_to_reduce, ri_t global_last_piv, */
		/* wl_t *waiting_global, */
		const mod_t modulus
		/* , omp_lock_t echelonize_lock */
		);

/**
 * \brief Echelonizes one multiline row (represented by dense_array_1 and
 * dense_array_2) by the already known pivot rows A->ml[first_piv] up to
 * A->ml[last_piv].
 *
 * \param multiline submatrix A (left upper side)
 *
 * \param dense array for first part of multiline row dense_array_1
 *
 * \param dense array for second part of multiline row dense_array_2
 *
 * \param first known pivot first_piv
 *
 * \param last known pivot last_piv
 *
 * \param field characteristic modulus
 */
void echelonize_one_row(sm_fl_ml_t *A, re_l_t *dense_array_1,
		re_l_t *dense_array_2,
		const ri_t first_piv, const ri_t last_piv,
		const mod_t modulus);

static inline void copy_piv_to_val(dm_t *D, const ri_t idx)
{
  // if D->row[idx]->init_val == NULL then D->row[idx] was already partly reduced
  // by some lower pivots. So we have the row already stored in large re_l_t
  // representation in D->row[idx]->val and can just use it
  if (D->row[idx]->piv_val != NULL) {
    ci_t i;
    D->row[idx]->val  = (re_l_t *)malloc(D->ncols * sizeof(re_l_t));
    for (i=0; i<D->ncols; ++i)
      D->row[idx]->val[i] = (re_l_t)D->row[idx]->piv_val[i];
    free(D->row[idx]->piv_val);
    D->row[idx]->piv_val  = NULL;
    D->row[idx]->lead     = D->row[idx]->piv_lead;
  }
}

static inline void copy_to_val(dm_t *D, const ri_t idx)
{
  // if D->row[idx]->init_val == NULL then D->row[idx] was already partly reduced
  // by some lower pivots. So we have the row already stored in large re_l_t
  // representation in D->row[idx]->val and can just use it
  if (D->row[idx]->init_val != NULL) {
    ci_t i;
    D->row[idx]->val  = (re_l_t *)malloc(D->ncols * sizeof(re_l_t));
    for (i=0; i<D->ncols; ++i)
      D->row[idx]->val[i] = (re_l_t)D->row[idx]->init_val[i];
    free(D->row[idx]->init_val);
    D->row[idx]->init_val  = NULL;
  }
}

/**
 * \brief Reduces dense row of index curr_row_to_reduce in matrix B with pivots that are
 * already interreduced from matrix D. we know that pivots of row index i have zero
 * entries for lead positions of pivots of index > i.
 *
 * \param dense row submatrix B
 *
 * \param dense row submatrix D
 *
 * \param row index curr_row_to_reduce
 */
static inline void reduce_B_by_D(dm_t *B, const dm_t *D, const ri_t curr_row_to_reduce)
{
  ri_t i;
  re_t mult1;
  i = 0;
  while (i<D->rank) {
    mult1 = 0;
    if (D->row[i]->piv_lead == B->row[curr_row_to_reduce]->lead)
      /* mult1  = MODP(B->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
      mult1  = (re_t)(B->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
    //  printf("pos: %u <= %u | lead of row %u = %lu\n",D->row[icurr_row_to_reduce,D->row[curr_row_to_reduce]->val[D->row[i]->lead]);
    //printf("mult for %u = %u\n",i,mult);
    if (mult1 != 0) {
      mult1  = (re_t)(D->mod - mult1);
      B->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
      //printf("inverted mult for %u = %u\n",i,mult);
      // also updates lead for row i
      //printf("lead in %u (from_row %u - reduced by %u) ",
      //    D->row[curr_row_to_reduce]->lead, from_row, i);
      reduce_dense_row_of_B_by_D(B, D, curr_row_to_reduce, i, mult1);
      //printf("--> out %u\n", D->row[curr_row_to_reduce]->lead);
      // if reduced row i is zero row then swap row down and get a new
      // row from the bottom
      if (B->row[curr_row_to_reduce]->val == NULL) {
        return;
      }
    }
    i++;
  }
  /*
  printf("reduction done: ");
  for (int ii=0; ii<B->ncols; ++ii) {
    printf("%lu ", B->row[curr_row_to_reduce]->val[ii]);
  }
  printf("\n");
  */
  save_pivot_without_normalization(B, curr_row_to_reduce, curr_row_to_reduce);
  /*
  for (int ii=0; ii<B->ncols; ++ii) {
    printf("%u ", B->row[curr_row_to_reduce]->piv_val[ii]);
  }
  printf("\n");
  */
}

/**
 * \brief Reduces dense row of index curr_row_to_reduce in matrix D with pivots that are
 * already interreduced from matrix D. we know that pivots of row index i have zero
 * entries for lead positions of pivots of index > i.
 *
 * \param dense row submatrix D
 *
 * \param row index curr_row_to_reduce
 */
static inline void completely_reduce_D(dm_t *D, const ri_t curr_row_to_reduce)
{
  ri_t i;
  re_t mult1;
  i = curr_row_to_reduce+1;
  while (i<D->rank) {
    mult1 = 0;
    if (D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] != 0)
      /* mult1  = MODP(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
      mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
    //  printf("pos: %u <= %u | lead of row %u = %lu\n",D->row[icurr_row_to_reduce,D->row[curr_row_to_reduce]->val[D->row[i]->lead]);
    //printf("mult for %u = %u\n",i,mult);
    if (mult1 != 0) {
      mult1  = (re_t)(D->mod - mult1);
      D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
      //printf("inverted mult for %u = %u\n",i,mult);
      // also updates lead for row i
      //printf("lead in %u (from_row %u - reduced by %u) ",
      //    D->row[curr_row_to_reduce]->lead, from_row, i);
      reduce_dense_row(D, curr_row_to_reduce, i, mult1);
      //reduce_dense_row_of_B_by_D(B, D, curr_row_to_reduce, i, mult1);
      //printf("--> out %u\n", D->row[curr_row_to_reduce]->lead);
      // if reduced row i is zero row then swap row down and get a new
      // row from the bottom
    }
    i++;
  }
  /*
  printf("reduction done: ");
  for (int ii=0; ii<B->ncols; ++ii) {
    printf("%lu ", B->row[curr_row_to_reduce]->val[ii]);
  }
  printf("\n");
  */
  save_pivot_without_normalization(D, curr_row_to_reduce, curr_row_to_reduce);
  /*
  for (int ii=0; ii<B->ncols; ++ii) {
    printf("%u ", B->row[curr_row_to_reduce]->piv_val[ii]);
  }
  printf("\n");
  */
}

/**
 * \brief Reduces dense row of index curr_row_to_reduce with pivots that are
 * already interreduced, i.e. we know that pivots of row index i have zero
 * entries for lead positions of pivots of index > i.
 *
 * \param dense row submatrix D
 *
 * \param row index curr_row_to_reduce
 *
 * \param row index from_row
 *
 * \param row index local_last_piv
 *
 */
static inline void reduce_dense_row_pre_elim(dm_t *D, const ri_t curr_row_to_reduce,
    const ri_t from_row, const ri_t local_last_piv)
{
  ri_t i, j;
  re_t mult1, mult2, mult3;
  i = from_row;
  if (local_last_piv > 0) {
    while (i<local_last_piv-1) {
      //printf("A crr[%u]->lead %u || %u i[%u]->plead\n", curr_row_to_reduce, D->row[curr_row_to_reduce]->lead, D->row[i]->piv_lead, i);
      if (D->row[i]->piv_lead >= D->row[curr_row_to_reduce]->lead) {
        if (D->row[i]->piv_lead == D->row[curr_row_to_reduce]->lead)
          mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead];
        else
          /* mult1  = MODP(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
          mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
        if (mult1 != 0) {
          mult1  = (re_t)(D->mod - mult1);
          D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
          for (j=D->row[i]->piv_lead+1; j<D->row[i+1]->piv_lead; ++j)
            D->row[curr_row_to_reduce]->val[j]  += (re_l_t)mult1 * D->row[i]->piv_val[j];
          /* mult2  = MODP(D->row[curr_row_to_reduce]->val[D->row[i+1]->piv_lead], D->mod); */
          mult2  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i+1]->piv_lead] % D->mod);
          if (mult2 != 0) {
            mult2 = (re_t)(D->mod - mult2);
            D->row[curr_row_to_reduce]->val[D->row[i+1]->piv_lead]  = 0;
            for (j=D->row[i+1]->piv_lead+1; j<D->row[i+2]->piv_lead; ++j)
              D->row[curr_row_to_reduce]->val[j]  +=
                ((re_l_t)mult1 * D->row[i]->piv_val[j] + (re_l_t)mult2 * D->row[i+1]->piv_val[j]);
            /* mult3  = MODP(D->row[curr_row_to_reduce]->val[D->row[i+2]->piv_lead], D->mod); */
            mult3  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i+2]->piv_lead] % D->mod);
            if (mult3 != 0) {
              mult3 = (re_t)(D->mod - mult3);
              D->row[curr_row_to_reduce]->val[D->row[i+2]->piv_lead]  = 0;
              reduce_dense_row_three(D, curr_row_to_reduce, i, i+1, i+2, mult1, mult2, mult3);
              // if reduced row i is zero row then swap row down and get a new
              // row from the bottom
              if (D->row[curr_row_to_reduce]->val == NULL) {
                return;
              }
            } else {
              reduce_dense_row_twice_from(D, curr_row_to_reduce, i, i+1, mult1, mult2, D->row[i+2]->piv_lead + 1);
            }
          } else {
            reduce_dense_row_from(D, curr_row_to_reduce, i, mult1, D->row[i+1]->piv_lead + 1);
            // if reduced row i is zero row then swap row down and get a new
            // row from the bottom
            if (D->row[curr_row_to_reduce]->val == NULL) {
              return;
            }
            i = i+2;
            continue;
          }
        } else {
          i++;
          continue;
        }
        i = i+3;
      } else {
        i++;
      }
    }
    if (i == local_last_piv-1) {
      if (D->row[i]->piv_lead >= D->row[curr_row_to_reduce]->lead) {
        if (D->row[i]->piv_lead == D->row[curr_row_to_reduce]->lead)
          mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead];
        else
          /* mult1  = MODP(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
          mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
        //  printf("pos: %u <= %u | lead of row %u = %lu\n",D->row[icurr_row_to_reduce,D->row[curr_row_to_reduce]->val[D->row[i]->lead]);
        //printf("mult for %u = %u\n",i,mult);
        if (mult1 != 0) {
          mult1  = (re_t)(D->mod - mult1);
          D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
          //printf("inverted mult for %u = %u\n",i,mult);
          // also updates lead for row i
          //printf("lead in %u (from_row %u - reduced by %u) ",
          //    D->row[curr_row_to_reduce]->lead, from_row, i);
          reduce_dense_row(D, curr_row_to_reduce, i, mult1);
          //printf("--> out %u\n", D->row[curr_row_to_reduce]->lead);
          // if reduced row i is zero row then swap row down and get a new
          // row from the bottom
          if (D->row[curr_row_to_reduce]->val == NULL) {
            return;
          }
        }
      }
      i++;
    }
    if (i == local_last_piv) {
      if (D->row[i]->piv_lead >= D->row[curr_row_to_reduce]->lead) {
        if (D->row[i]->piv_lead == D->row[curr_row_to_reduce]->lead)
          mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead];
        else
          /* mult1  = MODP(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
          mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
        //  printf("pos: %u <= %u | lead of row %u = %lu\n",D->row[icurr_row_to_reduce,D->row[curr_row_to_reduce]->val[D->row[i]->lead]);
        //printf("mult for %u = %u\n",i,mult);
        if (mult1 != 0) {
          mult1  = (re_t)(D->mod - mult1);
          D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
          //printf("inverted mult for %u = %u\n",i,mult);
          // also updates lead for row i
          //printf("lead in %u (from_row %u - reduced by %u) ",
          //    D->row[curr_row_to_reduce]->lead, from_row, i);
          reduce_dense_row(D, curr_row_to_reduce, i, mult1);
          //printf("--> out %u\n", D->row[curr_row_to_reduce]->lead);
          // if reduced row i is zero row then swap row down and get a new
          // row from the bottom
          if (D->row[curr_row_to_reduce]->val == NULL) {
            return;
          }
        }
      }
    }
  } else {
    if (i == local_last_piv) {
      if (D->row[i]->piv_lead >= D->row[curr_row_to_reduce]->lead) {
        if (D->row[i]->piv_lead == D->row[curr_row_to_reduce]->lead)
          mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead];
        else
          /* mult1  = MODP(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
          mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
        //  printf("pos: %u <= %u | lead of row %u = %lu\n",D->row[icurr_row_to_reduce,D->row[curr_row_to_reduce]->val[D->row[i]->lead]);
        //printf("mult for %u = %u\n",i,mult);
        if (mult1 != 0) {
          mult1  = (re_t)(D->mod - mult1);
          D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
          //printf("inverted mult for %u = %u\n",i,mult);
          // also updates lead for row i
          //printf("lead in %u (from_row %u - reduced by %u) ",
          //    D->row[curr_row_to_reduce]->lead, from_row, i);
          reduce_dense_row(D, curr_row_to_reduce, i, mult1);
          //printf("--> out %u\n", D->row[curr_row_to_reduce]->lead);
          // if reduced row i is zero row then swap row down and get a new
          // row from the bottom
          if (D->row[curr_row_to_reduce]->val == NULL) {
            return;
          }
        }
      }
    }
  }
}

/**
 * \brief Reduces dense row of index curr_row_to_reduce with pivots that are
 * no necessarily interreduced, i.e. it is possible hat pivots of row index i
 * have non-zero entries for lead positions of pivots of index > i.
 *
 * \param dense row submatrix D
 *
 * \param row index curr_row_to_reduce
 *
 * \param row index from_row
 *
 * \param row index local_last_piv
 *
 */
static inline void reduce_dense_row_general(dm_t *D, const ri_t curr_row_to_reduce,
    const ri_t from_row, const ri_t local_last_piv)
{
  ri_t i, j;
  re_t mult1, mult2, mult3;
  i = from_row;
  if (local_last_piv > 0) {
    while (i<local_last_piv-1) {
      //printf("crr[%u]->lead %u || %u i[%u]->plead | %u i+1[%u]->plead | %u i+2[%u]->plead\n", curr_row_to_reduce, D->row[curr_row_to_reduce]->lead, D->row[i]->piv_lead, i, D->row[i+1]->piv_lead, i+1, D->row[i+2]->piv_lead, i+2);
      if (D->row[i]->piv_lead >= D->row[curr_row_to_reduce]->lead) {
        if (D->row[i]->piv_lead == D->row[curr_row_to_reduce]->lead)
          mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead];
        else
          /* mult1  = MODP(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
          mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
        if (mult1 != 0) {
          mult1  = (re_t)(D->mod - mult1);
          D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
          for (j=D->row[i]->piv_lead+1; j<D->row[i+1]->piv_lead+1; ++j)
            D->row[curr_row_to_reduce]->val[j]  += (re_l_t)mult1 * D->row[i]->piv_val[j];
          /* mult2  = MODP(D->row[curr_row_to_reduce]->val[D->row[i+1]->piv_lead], D->mod); */
          mult2  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i+1]->piv_lead] % D->mod);
          if (mult2 != 0) {
            mult2 = (re_t)(D->mod - mult2);
            D->row[curr_row_to_reduce]->val[D->row[i+1]->piv_lead]  = 0;
            for (j=D->row[i+1]->piv_lead+1; j<D->row[i+2]->piv_lead+1; ++j)
              D->row[curr_row_to_reduce]->val[j]  +=
                ((re_l_t)mult1 * D->row[i]->piv_val[j] + (re_l_t)mult2 * D->row[i+1]->piv_val[j]);
            /* mult3  = MODP(D->row[curr_row_to_reduce]->val[D->row[i+2]->piv_lead], D->mod); */
            mult3  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i+2]->piv_lead] % D->mod);
            if (mult3 != 0) {
              mult3 = (re_t)(D->mod - mult3);
              D->row[curr_row_to_reduce]->val[D->row[i+2]->piv_lead]  = 0;
              reduce_dense_row_three(D, curr_row_to_reduce, i, i+1, i+2, mult1, mult2, mult3);
              // if reduced row i is zero row then swap row down and get a new
              // row from the bottom
              if (D->row[curr_row_to_reduce]->val == NULL) {
                return;
              }
            } else {
              reduce_dense_row_twice_from(D, curr_row_to_reduce, i, i+1, mult1, mult2, D->row[i+2]->piv_lead + 1);
              if (D->row[curr_row_to_reduce]->val == NULL) {
                return;
              }
            }
          } else {
            reduce_dense_row_from(D, curr_row_to_reduce, i, mult1, D->row[i+1]->piv_lead + 1);
            // if reduced row i is zero row then swap row down and get a new
            // row from the bottom
            if (D->row[curr_row_to_reduce]->val == NULL) {
              return;
            }
            i = i+2;
            continue;
          }
        } else {
          i++;
          continue;
        }
        i = i+3;
      } else {
        i++;
      }
    }
    if (i == local_last_piv-1) {
      //printf("crr[%u]->lead %u || %u i[%u]->plead\n", curr_row_to_reduce, D->row[curr_row_to_reduce]->lead, D->row[i]->piv_lead, i);
      if (D->row[i]->piv_lead >= D->row[curr_row_to_reduce]->lead) {
        if (D->row[i]->piv_lead == D->row[curr_row_to_reduce]->lead)
          mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead];
        else
          /* mult1  = MODP(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
          mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
        //  printf("pos: %u <= %u | lead of row %u = %lu\n",D->row[icurr_row_to_reduce,D->row[curr_row_to_reduce]->val[D->row[i]->lead]);
        //printf("mult for %u = %u\n",i,mult);
        if (mult1 != 0) {
          mult1  = (re_t)(D->mod - mult1);
          D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
          //printf("inverted mult for %u = %u\n",i,mult);
          // also updates lead for row i
          //printf("lead in %u (from_row %u - reduced by %u) ",
          //    D->row[curr_row_to_reduce]->lead, from_row, i);
          reduce_dense_row(D, curr_row_to_reduce, i, mult1);
          //printf("--> out %u\n", D->row[curr_row_to_reduce]->lead);
          // if reduced row i is zero row then swap row down and get a new
          // row from the bottom
          if (D->row[curr_row_to_reduce]->val == NULL) {
            return;
          }
        }
      }
      i++;
    }
    if (i == local_last_piv) {
      //printf("crr[%u]->lead %u || %u i[%u]->plead\n", curr_row_to_reduce, D->row[curr_row_to_reduce]->lead, D->row[i]->piv_lead, i);
      if (D->row[i]->piv_lead >= D->row[curr_row_to_reduce]->lead) {
        if (D->row[i]->piv_lead == D->row[curr_row_to_reduce]->lead)
          mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead];
        else
          /* mult1  = MODP(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
          mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
        //  printf("pos: %u <= %u | lead of row %u = %lu\n",D->row[icurr_row_to_reduce,D->row[curr_row_to_reduce]->val[D->row[i]->lead]);
        if (mult1 != 0) {
          mult1  = (re_t)(D->mod - mult1);
          D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
          //printf("inverted mult for %u = %u\n",i,mult);
          // also updates lead for row i
          //printf("lead in %u (from_row %u - reduced by %u) ",
          //    D->row[curr_row_to_reduce]->lead, from_row, i);
          reduce_dense_row(D, curr_row_to_reduce, i, mult1);
          //printf("--> out %u\n", D->row[curr_row_to_reduce]->lead);
          // if reduced row i is zero row then swap row down and get a new
          // row from the bottom
          if (D->row[curr_row_to_reduce]->val == NULL) {
            return;
          }
        }
      }
    }
  } else {
    if (i == local_last_piv) {
      if (D->row[i]->piv_lead >= D->row[curr_row_to_reduce]->lead) {
        if (D->row[i]->piv_lead == D->row[curr_row_to_reduce]->lead)
          mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead];
        else
          /* mult1  = MODP(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead], D->mod); */
          mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead] % D->mod);
        //  printf("pos: %u <= %u | lead of row %u = %lu\n",D->row[icurr_row_to_reduce,D->row[curr_row_to_reduce]->val[D->row[i]->lead]);
        //printf("mult for %u = %u\n",i,mult);
        if (mult1 != 0) {
          mult1  = (re_t)(D->mod - mult1);
          D->row[curr_row_to_reduce]->val[D->row[i]->piv_lead]  = 0;
          //printf("inverted mult for %u = %u\n",i,mult);
          // also updates lead for row i
          //printf("lead in %u (from_row %u - reduced by %u) ",
          //    D->row[curr_row_to_reduce]->lead, from_row, i);
          reduce_dense_row(D, curr_row_to_reduce, i, mult1);
          //printf("--> out %u\n", D->row[curr_row_to_reduce]->lead);
          // if reduced row i is zero row then swap row down and get a new
          // row from the bottom
          if (D->row[curr_row_to_reduce]->val == NULL) {
            return;
          }
        }
      }
    }
  }
}

/**
 * \brief Reduces dense row of index curr_row_to_reduce with pivots of index
 * from_row to local_last_piv in the task based reduction of D. Tries to reduce
 * a row with two pivots at once. This is the complete interreduction of already
 * computed pivots.
 *
 * \param dense row submatrix D
 *
 * \param row index curr_row_to_reduce
 *
 * \param row index from_row
 *
 * \param row index local_last_piv
 *
 */
static inline void completely_reduce_dense_row_task_new(dm_t *D, const ri_t curr_row_to_reduce,
    const ri_t from_row, const ri_t local_last_piv)
{
  // can we assume that all pivots up to local_last_piv have been fully reduced
  // with respect to all other pivots?

  // in very few cases it is possible that we have already a new pivot stored in
  // curr_row_to_reduce's piv_val, but we still have some not fully reduced data
  // in curr_row_to_reduce's val. In this cases copy_piv_to_val would overwrite
  // the data in curr_row_to_reduce's val. So we have to store it temporarily
  // and copy it back afterwards.
  // One example where this is happening is f4/eco/14/mat10.gbm.gz when the
  // number of threads used is between 3 - 5.
  re_l_t *tmp_val = NULL;
  ri_t tmp_lead = 0;
  if (D->row[curr_row_to_reduce]->val != NULL) {
    tmp_val   = D->row[curr_row_to_reduce]->val;
    tmp_lead  = D->row[curr_row_to_reduce]->lead; 
  }
  copy_piv_to_val(D, curr_row_to_reduce);
  reduce_dense_row_pre_elim(D, curr_row_to_reduce, from_row, local_last_piv);
  save_pivot(D, curr_row_to_reduce, curr_row_to_reduce);
  if (tmp_val != NULL) {
    D->row[curr_row_to_reduce]->val   = tmp_val;
    D->row[curr_row_to_reduce]->lead  = tmp_lead;
  }
}

/**
 * \brief Reduces dense row of index curr_row_to_reduce with pivots of index
 * from_row to local_last_piv in the task based reduction of D. Tries to reduce
 * a row with three pivots at once. This is the special call for the sequential
 * pre-elimination since there global_pre_elim == 0.
 *
 * \param dense row submatrix D
 *
 * \param row index curr_row_to_reduce
 *
 * \param row index from_row
 *
 * \param row index local_last_piv
 *
 */
static inline void reduce_dense_row_task_sequential(dm_t *D, const ri_t curr_row_to_reduce,
    const ri_t from_row, const ri_t local_last_piv)
{
  // can we assume that all pivots up to local_last_piv have been fully reduced
  // with respect to all other pivots?

  //printf("CURR ROW TO REDUCE %u ( %u ) - %u -- %u || %u\n", curr_row_to_reduce,D->row[curr_row_to_reduce]->lead,from_row,local_last_piv, global_pre_elim);
  copy_to_val(D, curr_row_to_reduce);
  
  if (D->row[curr_row_to_reduce]->val == NULL)
    return;

  reduce_dense_row_general(D, curr_row_to_reduce, from_row, local_last_piv);
}

/**
 * \brief Reduces dense row of index curr_row_to_reduce with pivots of index
 * from_row to local_last_piv in the task based reduction of D. Tries to reduce
 * a row with three pivots at once.
 *
 * \param dense row submatrix D
 *
 * \param row index curr_row_to_reduce
 *
 * \param row index from_row
 *
 * \param row index local_last_piv
 *
 */
static inline void reduce_dense_row_task_new(dm_t *D, const ri_t curr_row_to_reduce,
    const ri_t from_row, const ri_t local_last_piv)
{
  // can we assume that all pivots up to local_last_piv have been fully reduced
  // with respect to all other pivots?

  //printf("CURR ROW TO REDUCE %u ( %u ) - %u -- %u || %u\n", curr_row_to_reduce,D->row[curr_row_to_reduce]->lead,from_row,local_last_piv, global_pre_elim);
  copy_to_val(D, curr_row_to_reduce);
  
  if (D->row[curr_row_to_reduce]->val == NULL)
    return;

  //reduce_dense_row_general(D, curr_row_to_reduce, from_row, local_last_piv);
  reduce_dense_row_pre_elim(D, curr_row_to_reduce, from_row, global_pre_elim);
  reduce_dense_row_general(D, curr_row_to_reduce, global_pre_elim+1, local_last_piv);
}

/**
 * \brief Reduces dense row of index curr_row_to_reduce with pivots of index
 * from_row to local_last_piv in the task based reduction of D. Tries to reduce
 * a row with two pivots at once.
 *
 * \param dense row submatrix D
 *
 * \param row index curr_row_to_reduce
 *
 * \param row index from_row
 *
 * \param row index local_last_piv
 *
 */
static inline void reduce_dense_row_task(dm_t *D, const ri_t curr_row_to_reduce,
    const ri_t from_row, const ri_t local_last_piv)
{
  // can we assume that all pivots up to local_last_piv have been fully reduced
  // with respect to all other pivots?

  ri_t i, j;
  re_t mult1, mult2;
  copy_to_val(D, curr_row_to_reduce);
  i = from_row;
  while (i<local_last_piv) {
    if (D->row[i]->lead >= D->row[curr_row_to_reduce]->lead) {
      if (D->row[i]->lead == D->row[curr_row_to_reduce]->lead)
        mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->lead];
      else
        /* mult1  = (re_t)MODP(D->row[curr_row_to_reduce]->val[D->row[i]->lead], D->mod); */
        mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->lead] % D->mod);
      if (mult1 != 0) {
        mult1  = (re_t)(D->mod - mult1);
        D->row[curr_row_to_reduce]->val[D->row[i]->lead]  = 0;
        for (j=D->row[i]->lead+1; j<D->row[i+1]->lead+1; ++j)
          D->row[curr_row_to_reduce]->val[j]  += (re_l_t)mult1 * D->row[i]->piv_val[j];
        /* mult2  = (re_t)MODP(D->row[curr_row_to_reduce]->val[D->row[i+1]->lead], D->mod); */
        mult2  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i+1]->lead] % D->mod);
        if (mult2 != 0) {
          mult2 = (re_t)(D->mod - mult2);
          reduce_dense_row_twice(D, curr_row_to_reduce, i, i+1, mult1, mult2);
          // if reduced row i is zero row then swap row down and get a new
          // row from the bottom
          if (D->row[curr_row_to_reduce]->val == NULL) {
            return;
          }
        } else {
          reduce_dense_row_from(D, curr_row_to_reduce, i, mult1, D->row[i+1]->lead + 1);
          // if reduced row i is zero row then swap row down and get a new
          // row from the bottom
          if (D->row[curr_row_to_reduce]->val == NULL) {
            return;
          }
        }
      } else {
        i++;
        continue;
      }
      i = i+2;
    } else {
      i++;
    }
  }
  // if local_last_piv % 2 = 1 there is one pivot left
  if (i == local_last_piv) {
    if (D->row[i]->lead >= D->row[curr_row_to_reduce]->lead) {
      if (D->row[i]->lead == D->row[curr_row_to_reduce]->lead)
        mult1  = (re_t)D->row[curr_row_to_reduce]->val[D->row[i]->lead];
      else
        /* mult1  = (re_t)MODP(D->row[curr_row_to_reduce]->val[D->row[i]->lead], D->mod); */
        mult1  = (re_t)(D->row[curr_row_to_reduce]->val[D->row[i]->lead] % D->mod);
      //  printf("pos: %u <= %u | lead of row %u = %lu\n",D->row[icurr_row_to_reduce,D->row[curr_row_to_reduce]->val[D->row[i]->lead]);
      //printf("mult for %u = %u\n",i,mult);
      if (mult1 != 0) {
        mult1  = (re_t)(D->mod - mult1);
        //printf("inverted mult for %u = %u\n",i,mult);
        // also updates lead for row i
        //printf("lead in %u (from_row %u - reduced by %u) ",
        //    D->row[curr_row_to_reduce]->lead, from_row, i);
        reduce_dense_row(D, curr_row_to_reduce, i, mult1);
        //printf("--> out %u\n", D->row[curr_row_to_reduce]->lead);
        // if reduced row i is zero row then swap row down and get a new
        // row from the bottom
        if (D->row[curr_row_to_reduce]->val == NULL) {
          return;
        }
      }
    }
  }
}

/**
 * \brief Updates last pivot entry in global waiting list if a new pivot had to
 * be sorted in the list of pivots. If lp in waiting list is bigger than
 * from_row then we exchange them.
 *
 * \param waiting list wl
 *
 * \param row index of the inserted last pivot
 */
static inline void update_from_row_in_waiting_list(wl_t *waiting_global, const ri_t from_row)
{
  ri_t i;
  for (i=0; i<waiting_global->sz; ++i) {
    if (waiting_global->list[i].lp > from_row) {
      waiting_global->list[i].lp  = from_row;
    }
  }
}

/**
 * \brief Adds a new row to end of waiting list
 *
 * \param waiting list wl
 *
 * \param row index to be added row_idx
 *
 * \param row index of last pivot row_idx was reduced by last_piv_reduced_by
 */
static inline void push_row_to_waiting_list(wl_t *waiting_global, const ri_t row_idx,
    const ri_t last_piv_reduced_by)
{
#if DEBUG_ECHELONIZE
  int tid = omp_get_thread_num();
  printf("BEFORE PUSH\n");
  for (int ii=0; ii<waiting_global->sz; ++ii) {
    printf("T(%d) %d . %d\n",tid,waiting_global->list[ii].idx,waiting_global->list[ii].lp);
  }
#endif
  waiting_global->list[waiting_global->sz].idx  = row_idx;
  waiting_global->list[waiting_global->sz].lp   = last_piv_reduced_by;
#if DEBUG_ECHELONIZE
  printf("rowidx %u -- lprb %u\n", row_idx, last_piv_reduced_by);
#endif
  waiting_global->sz++;
#if DEBUG_ECHELONIZE
  printf("AFTER PUSH\n");
  for (int ii=0; ii<waiting_global->sz; ++ii) {
    printf("T(%d) %d . %d\n",tid,waiting_global->list[ii].idx,waiting_global->list[ii].lp);
  }
#endif
}

static inline void realloc_rows_ml(sm_fl_ml_t *A, const mli_t mli,
    const bi_t init_buffer_A, mli_t *buffer_A) {
  *buffer_A +=  init_buffer_A;
  A->ml[mli].idx = (mli_t*) realloc(A->ml[mli].idx, (*buffer_A) * sizeof(mli_t));
  A->ml[mli].val = (re_t*)  realloc(A->ml[mli].val, 2 * (*buffer_A) * sizeof(re_t));
}

static inline void realloc_block_rows(sbm_fl_t *A, const ri_t rbi, const ci_t bir,
    const bi_t lib, const bi_t init_buffer_A, bi_t *buffer_A) {
  *buffer_A +=  init_buffer_A;
  A->blocks[rbi][bir][lib].idx = realloc(
      A->blocks[rbi][bir][lib].idx,
      (*buffer_A) * sizeof(bi_t));
  A->blocks[rbi][bir][lib].val = realloc(
      A->blocks[rbi][bir][lib].val,
      2 * (*buffer_A) * sizeof(re_t));
}

static inline void insert_row_data_ml_1_1(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi1, const ci_t i1) {
  A->ml[mli].idx[A->ml[mli].sz]       = eil;
  A->ml[mli].val[2*A->ml[mli].sz]     = M->rows[bi1][i1];
  A->ml[mli].val[(2*A->ml[mli].sz)+1] = 0;
  A->ml[mli].sz++;
}

static inline void insert_row_data_ml_1_2(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi2, const ci_t i2) {
  /* printf("sz %d -- %d\n", A->ml[mli].sz, eil); */
  A->ml[mli].idx[A->ml[mli].sz]       = eil;
  A->ml[mli].val[2*A->ml[mli].sz]     = 0;
  A->ml[mli].val[(2*A->ml[mli].sz)+1] = M->rows[bi2][i2];
  A->ml[mli].sz++;
}

static inline void insert_row_data_ml_2(sm_fl_ml_t *A, const sm_t *M,
    const mli_t mli, const ci_t eil, const ci_t bi1, const ci_t i1,
    const ci_t bi2, const ci_t i2) {
  A->ml[mli].idx[A->ml[mli].sz]       = eil;
  A->ml[mli].val[2*A->ml[mli].sz]     = M->rows[bi1][i1];
  A->ml[mli].val[(2*A->ml[mli].sz)+1] = M->rows[bi2][i2];
  A->ml[mli].sz++;
}


static inline void insert_block_row_data_ml_1_1(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi1, const ci_t i1) {
  A->blocks[rbi][bir][lib].idx[A->blocks[rbi][bir][lib].sz]       = eil;
  A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz]     = M->rows[bi1][i1];
  A->blocks[rbi][bir][lib].val[(2*A->blocks[rbi][bir][lib].sz)+1] = 0;
  A->blocks[rbi][bir][lib].sz++;
}

static inline void insert_block_row_data_ml_1_2(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi2, const ci_t i2) {
  A->blocks[rbi][bir][lib].idx[A->blocks[rbi][bir][lib].sz]       = eil;
  A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz]     = 0;
  A->blocks[rbi][bir][lib].val[(2*A->blocks[rbi][bir][lib].sz)+1] = M->rows[bi2][i2];
  A->blocks[rbi][bir][lib].sz++;
}

static inline void insert_block_row_data_ml_2(sbm_fl_t *A, const sm_t *M,
    const bi_t rbi, const bi_t bir, const bi_t lib, const bi_t eil,
    const ci_t bi1, const ci_t i1, const ci_t bi2, const ci_t i2) {
  A->blocks[rbi][bir][lib].idx[A->blocks[rbi][bir][lib].sz]       = eil;
  A->blocks[rbi][bir][lib].val[2*A->blocks[rbi][bir][lib].sz]     = M->rows[bi1][i1];
  A->blocks[rbi][bir][lib].val[(2*A->blocks[rbi][bir][lib].sz)+1] = M->rows[bi2][i2];
  A->blocks[rbi][bir][lib].sz++;
}


void swap_block_data(sbm_fl_t *A, const ci_t clA, const bi_t rbi,
    const bi_t cvb);
/**
 * \brief Restores the two dense arrays in a multiline row in D after
 * echelonizing the corresponding multiline row
 *
 * \param multiline row in D ml
 *
 * \param dense array for first part of multiline row dense_array_1
 *
 * \param dense array for second part of multiline row dense_array_2
 *
 * \param columns dimension coldim
 *
 * \param field characteristic modulus
 *
 * \param reduce or just store in multiline row? reduce
 * if 1 then reduction is done, else we only store the multiline row
 */
void save_back_and_reduce(ml_t *ml, re_l_t *dense_array_1,
		re_l_t *dense_array_2, const ci_t coldim, const mod_t modulus,
		const int reduce);

/*  global variables for echelonization of D */
extern  omp_lock_t echelonize_lock;
extern  ri_t global_next_row_to_reduce;
extern  ri_t global_last_piv;
extern  ri_t global_last_row_fully_reduced;
extern  wl_t waiting_global;

extern ri_t global_initial_D_rank;
extern ri_t global_piv_lead_drop;

/*  global variables for reduction of C */
extern  omp_lock_t reduce_C_lock;
extern  ri_t reduce_C_next_col_to_reduce;

/*  global variables for reduction of C */
extern  omp_lock_t reduce_A_lock;
extern  ri_t reduce_A_next_col_to_reduce;

extern ri_t pivs_in_use;
extern ri_t sorting_pivs;
extern ri_t last_pivot_idx;

#endif /* GBLA_ELIMINATION_H */

/* vim:sts=2:sw=2:ts=2:
*/
