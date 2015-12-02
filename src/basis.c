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
 * \file basis.h
 * \brief Implementation of handling of groebner basis.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "basis.h"

inline void free_basis_dynamic_data(gb_t *basis)
{
  if (basis) {
    nelts_t i, j;
    for (i=0; i<basis->nvars; ++i) {
      free(basis->vnames[i]);
    }
    free(basis->vnames);
    free(basis->deg);
    free(basis->nt);
    for (i=0; i<basis->load; ++i) {
      free(basis->cf[i]);
      free(basis->eh[i]);
    }
    free(basis->cf);
    free(basis->eh);
  }
}
