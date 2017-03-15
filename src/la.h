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
 * \file la.h
 * \brief Implementation of the linear algebra parts.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_LA_H
#define GB_LA_H

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <config.h>
#include "types.h"

static inline void copy_first_row_from_dense(nelts_t *len, bs_t *pos,
  cf_t *val, const mat_gb_block_t *bl, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  len[0]  = 0;
  len[1]  = len[0];
  
  for (i=0; i<meta->bs; ++i) {
    if (bl->val[i] != 0) {
      pos[len[1]] = (bs_t)i;
      val[len[1]] = bl->val[i];
      len[1]++;
    }
  }
}

static inline void copy_first_row_from_sparse(nelts_t *len, bs_t *pos,
  cf_t *val, const mat_gb_block_t *bl)
{
  len[0]  = bl->len[0];
  len[1]  = bl->len[1];
  memcpy(pos, bl->pos, (len[1]-len[0]) * sizeof(bs_t));
  memcpy(val, bl->val, (len[1]-len[0]) * sizeof(cf_t));
}

static inline void load_dense_row_for_update_from_dense(bf_t *dr,
    const nelts_t idx, const mat_gb_block_t *bl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  memset(dr, 0, meta->bs * sizeof(bf_t));

  for (i=0; i<meta->bs; ++i)
    dr[i]  = (bf_t)bl->val[idx*meta->bs+i];
}

static inline void load_dense_row_for_update_from_sparse(bf_t *dr,
    const nelts_t idx, const mat_gb_block_t *bl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  memset(dr, 0, meta->bs * sizeof(bf_t));

  for (i=bl->len[idx]; i<bl->len[idx+1]; ++i)
    dr[bl->pos[i]]  = (bf_t)bl->val[i];
}

static inline void  update_dense_row(bf_t *dr, const nelts_t idx,
    const nelts_t *len, const bs_t *pos, const cf_t *val,
    const mat_gb_block_t *fbl, const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  for (i=fbl->len[idx]+1; i<fbl->len[idx+1]; ++i) {
    const bf_t mul    = (bf_t)fbl->val[i];
    const nelts_t ri  = meta->bs - fbl->pos[i] - 1;

    for (j=len[ri]; j<len[ri+1]; ++j)
      dr[pos[j]] +=  mul * val[j];
  }
}

static inline void write_updated_row_to_sparse_format(nelts_t *len,
    bs_t *pos, cf_t *val, const bf_t *dr, const nelts_t idx,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i;
  bf_t tmp;

  len[idx+1]  = len[idx];
  for (i=0; i<meta->bs; ++i) {
    if (dr[i] != 0) {
      tmp = dr[i] % meta->mod;
      if (tmp != 0) {
        pos[len[idx+1]] = (bs_t)i;
        val[len[idx+1]] = (cf_t)tmp;
        len[idx+1]++;
      }
    }
  }
}

static inline void set_updated_block(mat_gb_block_t *bl, nelts_t *len,
    bs_t *pos, cf_t *val)
{
  free(bl->len);
  free(bl->pos);
  free(bl->val);

  bl->len  = len;
  bl->pos  = pos;
  bl->val  = val;
}

static inline void update_single_block_dense(mat_gb_block_t *mat,
    const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  const nelts_t bs_square = (nelts_t)meta->bs * meta->bs;
  /* first block used as lookup table for updates */
  const mat_gb_block_t *fbl  = mat;
  /* block to be updated */
  mat_gb_block_t *ubl  = mat+idx;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  nelts_t *new_len  = (nelts_t *)malloc((meta->bs+1) * sizeof(nelts_t));
  bs_t *new_pos     = (bs_t *)malloc(bs_square * sizeof(bs_t));
  cf_t *new_val     = (cf_t *)malloc(bs_square * sizeof(cf_t));

  /* copy first row, it is not changed */
  copy_first_row_from_dense(new_len, new_pos, new_val, ubl, meta);

  for (i=1; i<meta->bs; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update_from_dense(dr, i, ubl, meta);

    /* find corresponding row and multiplier */
    update_dense_row(dr, i, new_len, new_pos, new_val, fbl, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(new_len, new_pos, new_val, dr, i, meta);
  }
  set_updated_block(ubl, new_len, new_pos, new_val);
}

static inline void update_single_block_sparse(mat_gb_block_t *mat,
    const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  const nelts_t bs_square = (nelts_t)meta->bs * meta->bs;
  /* first block used as lookup table for updates */
  const mat_gb_block_t *fbl  = mat;
  /* block to be updated */
  mat_gb_block_t *ubl  = mat+idx;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  nelts_t *new_len  = (nelts_t *)malloc((meta->bs+1) * sizeof(nelts_t));
  bs_t *new_pos     = (bs_t *)malloc(bs_square * sizeof(bs_t));
  cf_t *new_val     = (cf_t *)malloc(bs_square * sizeof(cf_t));

  /* copy first row, it is not changed */
  copy_first_row_from_sparse(new_len, new_pos, new_val, ubl);

  for (i=1; i<meta->bs; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update_from_sparse(dr, i, ubl, meta);

    /* find corresponding row and multiplier */
    update_dense_row(dr, i, new_len, new_pos, new_val, fbl, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(new_len, new_pos, new_val, dr, i, meta);
  }
  set_updated_block(ubl, new_len, new_pos, new_val);
}

static inline void update_single_block(mat_gb_block_t *mat,
    const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  if (mat[idx].len != NULL)
    update_single_block_sparse(mat, idx, meta);
  else
    update_single_block_dense(mat, idx, meta);
}

static inline void update_upper_row_block(mat_gb_block_t *mat,
    const mat_gb_meta_data_t *meta, const int t)
{
  nelts_t i;

  #pragma omp parallel num_threads(t)
  {
    #pragma omp single nowait
    {
      /* the first block is used for updated the remaining ones,
       * i.e. we start at block index i=1 */
      for (i=1; i<meta->ncb_AC; ++i) {
        #pragma omp task
        update_single_block(mat, i, meta);
      }
    }
  }
  /* check density of blocks */
  adjust_block_row_types(mat, meta);
}


#endif
