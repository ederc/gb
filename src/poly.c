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
 * \file poly.h
 * \brief Implementation of handling of polynomials and groebner bases.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "poly.h"
/*
gb_t *initialize_basis(const nelts_t size)
{
  gb_t *basis = (gb_t *)malloc(sizeof(gb_t));

  basis->size   = size;
  basis->load   = 0;
  basis->nv  = nvars;
  basis->mod    = modulus;
  basis->nt     = (nelts_t *)malloc(basis->size * sizeof(nelts_t));
  basis->deg    = (deg_t *)malloc(basis->size * sizeof(deg_t));
  basis->cf     = (coeff_t **)malloc(basis->size * sizeof(coeff_t *));
  basis->eh     = (hash_t **)malloc(basis->size * sizeof(hash_t *));

  return basis;
}
*/
inline void free_basis(gb_t *basis)
{
  if (basis) {
    nelts_t i;
    for (i=0; i<basis->nv; ++i) {
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

  free(basis);
  basis = NULL;
}

void add_new_element_to_basis_grevlex(gb_t *basis, const mat_t *mat,
    const nelts_t ri, const spd_t *spd, const mp_cf4_ht_t *ht)
{
  nelts_t i;

  if (basis->load == basis->size)
    enlarge_basis(basis, 2*basis->size);

  // maximal size is DR->ncols - row's piv lead
  nelts_t ms  = mat->DR->ncols - mat->DR->row[ri]->piv_lead;
  // use shorter names in here
  basis->cf[basis->load]  = (coeff_t *)malloc(ms * sizeof(coeff_t)); 
  basis->eh[basis->load]  = (hash_t *)malloc(ms * sizeof(hash_t)); 
  coeff_t *cf = basis->cf[basis->load];
  hash_t *eh  = basis->eh[basis->load];
  
  nelts_t ctr = 0;
  deg_t deg   = 0;
  for (i=mat->DR->row[ri]->piv_lead; i<mat->DR->ncols; ++i) {
    if (mat->DR->row[ri]->piv_val[i] != 0) {
      cf[ctr] = mat->DR->row[ri]->piv_val[i];
      eh[ctr] = spd->col->hpos[i];
      deg = ht->deg[eh[ctr]] > deg ? ht->deg[eh[ctr]] : deg; 
      ctr++;
    }
  }
  basis->nt[basis->load]  = ctr;
  basis->deg[basis->load] = deg;

  // realloc memory to the correct number of terms
  cf  = realloc(cf, basis->nt[basis->load] * sizeof(coeff_t)); 
  eh  = realloc(eh, basis->nt[basis->load] * sizeof(hash_t)); 
  basis->load++;
}

void enlarge_basis(gb_t *basis, const nelts_t size)
{
  basis->size = size;
  basis->nt   = realloc(basis->nt, basis->size * sizeof(nelts_t));
  basis->deg  = realloc(basis->deg, basis->size * sizeof(deg_t));
  basis->cf   = realloc(basis->cf, basis->size * sizeof(coeff_t *));
  basis->eh   = realloc(basis->eh, basis->size * sizeof(hash_t *));
}
