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

#define MATRIX_DEBUG 0

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
 * \param first row for the righthand side in matrix fr
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
void generate_row_blocks(sb_fl_t * A, dbm_fl_t *B, const nelts_t rbi,
    const nelts_t nr, const nelts_t fr, const bi_t bs, const nelts_t ncb,
    const gb_t *basis, const sel_t *sel, const pre_t *col);

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
 * \param first row for the righthand side in matrix fr
 *
 * \param block size bs
 *
 * \param intermediate groebner basis basis
 *
 * \param hash table ht
 */
void store_in_buffer(dbr_t *dbr, const nelts_t pi, const hash_t mul,
    const nelts_t fr, const bi_t bs, const gb_t *basis, const mp_cf4_ht_t *ht);

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
 * \param first row for the righthand side in matrix fr
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
#endif
