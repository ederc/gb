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
 * \file tools.c
 * \brief Implementation of smaller tools
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

static inline int32_t mod_p_inverse_32(
    const int32_t val,
    const int32_t p
    )
{
  int32_t a, b, c, d, e, f;
  a =   p;
  b =   val % p;
  /* if b < 0 we shift correspondingly */
  b +=  (b >> 31) & p;
  c =   1;
  d =   0;

  while (b != 0) {
    f = b;
    e = a/f;
    b = a - e*f;
    a = f;
    f = c;
    c = d - e*f;
    d = f;
  }

  /* if d < 0 we shift correspondingly */
  d +=  (d >> 31) & p;

  return d;
}

static inline void free_matrix(
    val_t **mat
    )
{
  int32_t i;

  for (i = 0; i < nrall; ++i) {
    free(mat[i]);
  }
  free(mat);
}
