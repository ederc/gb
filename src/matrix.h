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
mat_t *initialize_gbla_matrix(const spd_t *spd, const gb_t *basis);

/**
 * \brief Frees gbla matrix after reduction.
 *
 * \param gbla matrix mat
 */
void free_gbla_matrix(mat_t *mat);

/**
 * \brief Initializes a dense block row for buffering values when converting a
 * polynomial to a row in the gbla matrix.
 *
 * \param number of blocks nb
 *
 * \param block size bs
 *
 * \return dense block row
 */
dbr_t *initialize_dense_block_row(const nelts_t nb, const bi_t bs);

/**
 * \brief Frees dense block row.
 *
 * \param dense block row dbr
 *
 * \param number of blocks nb
 */
void free_dense_block_row(dbr_t *dbr, const nelts_t nb);

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
void allocate_sparse_block(sb_fl_t *A, const nelts_t rbi, const nelts_t bir,
    const bi_t bs);

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
void allocate_densee_block(sb_fl_t *A, const nelts_t rbi, const nelts_t bir,
    const bi_t bs);

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
 * \param polynomial index in basis pi
 *
 * \param hash position of multiplier mul
 *
 * \param first column on the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param intermediate groebner basis basis
 *
 * \param hash table ht
 */
inline void store_in_buffer(dbr_t *dbr, const nelts_t pi, const hash_t mul,
    const nelts_t fr, const bi_t bs, const gb_t *basis, const mp_cf4_ht_t *ht)
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
  if (basis->nt[pi]>3) {
    for (j=0; j<basis->nt[pi]-3; j=j+4) {
      hp  = find_in_hash_table_product(mul,basis->eh[pi][j], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j, basis->eh[pi][j], cp, basis->cf[pi][j]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = basis->cf[pi][j];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = basis->cf[pi][j];
        dbr->ctr[fbr+cp/bs]++;
      }
      hp  = find_in_hash_table_product(mul, basis->eh[pi][j+1], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+1, basis->eh[pi][j+1], cp, basis->cf[pi][j+1]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = basis->cf[pi][j+1];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = basis->cf[pi][j+1];
        dbr->ctr[fbr+cp/bs]++;
      }
      hp  = find_in_hash_table_product(mul, basis->eh[pi][j+2], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+2, basis->eh[pi][j+2], cp, basis->cf[pi][j+2]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = basis->cf[pi][j+2];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = basis->cf[pi][j+2];
        dbr->ctr[fbr+cp/bs]++;
      }
      hp  = find_in_hash_table_product(mul, basis->eh[pi][j+3], ht);
      cp  = ht->idx[hp];
#if MATRIX_DEBUG
      printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j+3, basis->eh[pi][j+3], cp, basis->cf[pi][j+3]);
#endif
      if (cp<fr) {
        dbr->cf[cp/bs][cp%bs]  = basis->cf[pi][j+3];
        dbr->ctr[cp/bs]++;
      } else {
        cp = cp - fr;
        dbr->cf[fbr+cp/bs][cp%bs]  = basis->cf[pi][j+3];
        dbr->ctr[fbr+cp/bs]++;
      }
    }
  }
  tmp = j;
  for (j=tmp; j<basis->nt[pi]; ++j) {
    hp  = find_in_hash_table_product(mul, basis->eh[pi][j], ht);
    cp  = ht->idx[hp];
#if MATRIX_DEBUG
    printf("fr %u | hp %u | eh[%u][%u] %u | cp %u | cf %u\n", fr, hp, pi, j, basis->eh[pi][j], cp, basis->cf[pi][j]);
#endif
    if (cp<fr) {
      dbr->cf[cp/bs][cp%bs]  = basis->cf[pi][j];
      dbr->ctr[cp/bs]++;
    } else {
      cp = cp - fr;
      dbr->cf[fbr+cp/bs][cp%bs]  = basis->cf[pi][j];
      dbr->ctr[fbr+cp/bs]++;
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
 *
 * \param symbolic preprocessing monomials col
 */
inline void generate_row_blocks(sb_fl_t * A, dbm_fl_t *B, const nelts_t rbi,
    const nelts_t nr, const nelts_t fr, const bi_t bs, const nelts_t ncb,
    const gb_t *basis, const sel_t *sel, const pre_t *col)
{
  nelts_t i;
  // get new row index in block rib
  bi_t rib;
  // polynomial index in basis
  nelts_t pi;
  // multiplier
  hash_t mul;

  // preallocate buffer to store row in dense format
  dbr_t *dbr  = initialize_dense_block_row(ncb, bs);

  nelts_t min = (rbi+1)*bs > nr ? nr : (rbi+1)*bs;
  // for each row we allocate memory in the sparse, left side and go through the
  // polynomials and add corresponding entries in the matrix
  for (i=rbi*bs; i<min; ++i) {
    // zero out buffer data
    reset_buffer(dbr, ncb, bs);

    rib = i % bs;
    pi  = sel->mpp[i].idx;
    mul = sel->mpp[i].mul;

    store_in_buffer(dbr, pi, mul, fr, bs, basis, ht);
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
void store_in_matrix(sb_fl_t *A, dbm_fl_t *B, const dbr_t *dbr, const nelts_t  rbi,
    const nelts_t rib, const nelts_t ncb, const nelts_t fr, const bi_t bs,
    const coeff_t mod);

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
void write_to_sparse_row(sb_fl_t *A, const coeff_t *cf, const nelts_t rbi,
    const nelts_t bir, const nelts_t rib, const nelts_t sz,
    const bi_t bs, const coeff_t mod);

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
void write_to_dense_row(dbm_fl_t *A, const coeff_t *cf, const nelts_t rbi,
    const nelts_t bir, const nelts_t rib, const bi_t bs);

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
ci_t get_number_of_left_column_blocks(const pre_t *col, const nelts_t bs);

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
ci_t get_number_of_right_column_blocks(const pre_t *col, const nelts_t bs);

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
ri_t get_number_of_row_blocks(const sel_t *sel, const nelts_t bs);

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
 * \param symbplic preprocessing data spd
 *
 * \param number of threads nthreads
 *
 * \return gbla matrix mat
 */
mat_t *generate_gbla_matrix(const gb_t *basis, const spd_t *spd,
    const int nthreads);

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
 * \brief Resets buffer to all entries zero once the a row is done.
 *
 * \param dense block row buffer dbr
 *
 * \param number of column blocks ncb
 *
 * \param block size bs
 */
void reset_buffer(dbr_t *dbr, const nelts_t ncb, const bi_t bs);
#endif
