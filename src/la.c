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
 * \file la.c
 * \brief Implementation of linear algebra.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "data.h"

static inline void normalize_matrix_row(
    val_t *row
    )
{
  int32_t i;

  const len_t len   = row[0];
  const int64_t inv = (int64_t)mod_p_inverse_32(row[3], fc);
  
  i = (len-2)/2;
  i = i & 1 ? 5 : 3;
  for (; i < len; i = i+4) {
    row[i]    = (val_t)(row[i] * inv) % fc;
    row[i+2]  = (val_t)(row[i+2] * inv) % fc;
  }
  /* probably we left out the first coefficient */
  row[3]  = 1;
}
