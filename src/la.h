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

static inline void copy_first_row_from_dense(mat_gb_block_t *bl,
    const mat_gb_block_t *obl, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  bl->len[0]  = 0;
  bl->len[1]  = bl->len[0];
  
  for (i=0; i<meta->bs; ++i) {
    if (obl->val[i] != 0) {
      bl->pos[bl->len[1]] = (bs_t)i;
      bl->val[bl->len[1]] = obl->val[i];
      bl->len[1]++;
    }
  }
}

static inline void copy_first_row_from_sparse(mat_gb_block_t *bl,
    const mat_gb_block_t *obl)
{
  bl->len[0]  = obl->len[0];
  bl->len[1]  = obl->len[1];
  memcpy(bl->pos, obl->pos, (bl->len[1]-bl->len[0]) * sizeof(bs_t));
  memcpy(bl->val, obl->val, (bl->len[1]-bl->len[0]) * sizeof(cf_t));
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

static inline void write_updated_row_to_sparse_format(mat_gb_block_t *bl,
    const bf_t *dr, const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;
  bf_t tmp;

  bl->len[idx+1]  = bl->len[idx];
  for (i=0; i<meta->bs; ++i) {
    if (dr[i] != 0) {
      tmp = dr[i] % meta->mod;
      if (tmp != 0) {
        bl->pos[bl->len[idx+1]] = (bs_t)i;
        bl->val[bl->len[idx+1]] = (cf_t)tmp;
        bl->len[idx+1]++;
      }
    }
  }
}

static inline void set_updated_block(mat_gb_block_t **bl, mat_gb_block_t* new_bl)
{

  free((*bl)->len);
  free((*bl)->pos);
  free((*bl)->val);
  free((*bl));

  *bl = new_bl;
}

static inline void update_dense_row(bf_t *dr, const nelts_t idx,
    const mat_gb_block_t *bl, const mat_gb_block_t *mbl,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i, j;

  for (i=mbl->len[idx]+1; i<mbl->len[idx+1]; ++i) {
    const bf_t mul    = (bf_t)mbl->val[i];
    const nelts_t ri  = meta->bs - mbl->pos[i] - 1;

    for (j=bl->len[ri]; j<bl->len[ri+1]; ++j)
      dr[bl->pos[j]] +=  mul * bl->val[j];
  }
}

static inline void update_single_block_dense(mat_gb_block_t *mat,
    const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* first block used as lookup table for updates */
  const mat_gb_block_t *fbl  = mat;
  /* block to be updated */
  mat_gb_block_t *ubl  = mat+idx;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = initialize_mat_gb_block(meta);

  /* copy first row, it is not changed */
  copy_first_row_from_dense(nbl, ubl, meta);

  for (i=1; i<meta->bs; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update_from_dense(dr, i, ubl, meta);

    /* find corresponding row and multiplier */
    update_dense_row(dr, i, nbl, fbl, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(nbl, dr, i, meta);
  }

  set_updated_block(&ubl, nbl);
}

static inline void update_single_block_sparse(mat_gb_block_t *mat,
    const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* first block used as lookup table for updates */
  const mat_gb_block_t *fbl  = mat;
  /* block to be updated */
  mat_gb_block_t *ubl  = mat+idx;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = initialize_mat_gb_block(meta);

  /* copy first row, it is not changed */
  copy_first_row_from_sparse(nbl, ubl);

  for (i=1; i<meta->bs; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update_from_sparse(dr, i, ubl, meta);

    /* find corresponding row and multiplier */
    update_dense_row(dr, i, nbl, fbl, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(nbl, dr, i, meta);
  }

  set_updated_block(&ubl, nbl);
}

static inline void update_single_block(mat_gb_block_t *mat,
    const nelts_t idx, const mat_gb_meta_data_t *meta)
{
  if (mat[idx].len != NULL) {
    update_single_block_sparse(mat, idx, meta);
  } else {
    if (mat[idx].val != NULL)
      update_single_block_dense(mat, idx, meta);
    else
      return;
  }
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
  /* adjust_block_row_types(mat, meta); */
}

static inline void sparse_update_lower_block_by_upper_block(mat_gb_block_t *l,
    const mat_gb_block_t *u, const nelts_t rbi, const nelts_t cbi,
    const mat_gb_meta_data_t *meta)
{
  nelts_t i;

  /* first multiplier block used as lookup table for updates */
  const mat_gb_block_t *mbl  = l + rbi*meta->ncb_AC;
  /* block to be updated */
  mat_gb_block_t *ubl  = l + rbi*meta->ncb_AC+cbi;
  /* row to be updated in dense format */
  bf_t *dr  = (bf_t *)malloc(meta->bs * sizeof(bf_t));

  mat_gb_block_t *nbl = initialize_mat_gb_block(meta);

  for (i=0; i<meta->bs; ++i) {
    /* load dense row to be updated */
    load_dense_row_for_update_from_sparse(dr, i, ubl, meta);

    /* find corresponding row and multiplier */
    update_dense_row(dr, i, u, mbl, meta);

    /* write updated row to new storage holders */
    write_updated_row_to_sparse_format(nbl, dr, i, meta);
  }

  set_updated_block(&ubl, nbl);
}

static inline void update_lower_by_upper_row_block(mat_gb_block_t *l,
    const mat_gb_block_t *u, const mat_gb_meta_data_t *meta, const int t)
{
  nelts_t i, j;

#pragma omp parallel num_threads(t)
  {
#pragma omp single nowait
    {
      /* the first block is used for updated the remaining ones,
       * i.e. we start at block index i=1 */
      for (i=0; i<meta->nrb_CD; ++i) {
        /* need to look at the first block in each row in
         * order to decide which algorithm to be chosen */
        if (l[i+meta->ncb_AC].len != NULL) {
          for (j=1; j<meta->ncb_AC; ++j) {
#pragma omp task
            {
              sparse_update_lower_block_by_upper_block(l, u, i, j, meta);
            }
          }
        } else {
          continue;
        }
      }
    }
/* at the moment we only work with sparse blocks */
#if 0 
#pragma omp single nowait
    {
      for (i=0; i<meta->nrb_CD; ++i) {
#pragma omp task
        {
          /* check density of blocks */
          adjust_block_row_types(l+i*meta->ncb_AC, meta);
        }
      }
    }
#endif
  }
}
#endif
