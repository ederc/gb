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

static inline void update_single_block(mat_gb_block_t *mat,
    const nelts_t idx, const mat_gb_meta_data_t *meta)
{
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

}


#endif
