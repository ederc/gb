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
 * \file io.c
 * \brief Input and output handling
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

static val_t **import_julia_data(
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const int32_t nr_gens
    )
{
  int32_t i, j;
  int32_t off = 0; /* offset in arrays */
  val_t **mat = malloc((unsigned long)nr_gens * sizeof(val_t *));
  
  for (i = 0; i < nr_gens; ++i) {
    /* each matrix row has the following structure:
     * [length | eh1 | cf1 | eh2 | cf2 | .. | ehl | cfl | piv?]
     * where piv? is a label for being a known or an unknown pivot */
    mat[i]    = malloc((2*(unsigned long)lens[i]+2) * sizeof(val_t));
    mat[i][0] = 2*lens[i]+2;
    for (j = off; j < off+lens[i]; ++j) {
      mat[i][2*(j-off)+1] = insert_in_local_hash_table(exps+(nvars*j));
      mat[i][2*(j-off)+2] = cfs[j];
    }
    /* mark initial generators, they have to be added to the basis first */
    mat[i][2*lens[i]+1] = 1;
    off +=  lens[i];
  }

  return mat;
}
