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
inline gb_t *initialize_basis(const gb_t *input)
{
  gb_t *basis = (gb_t *)malloc(sizeof(gb_t));

  // get meta data from input
  basis->size   = 3 * input->size;
  basis->load   = 0;
  basis->nv     = input->nv;
  basis->mod    = input->mod;
  basis->vnames = input->vnames;
  
  // allocate memory for elements
  basis->nt     = (nelts_t *)malloc(basis->size * sizeof(nelts_t));
  basis->deg    = (deg_t *)malloc(basis->size * sizeof(deg_t));
  basis->red    = (red_t *)malloc(basis->size * sizeof(red_t));
  basis->cf     = (coeff_t **)malloc(basis->size * sizeof(coeff_t *));
  basis->eh     = (hash_t **)malloc(basis->size * sizeof(hash_t *));
  printf("bred %p\n", basis->red);

  // initialize element at position 0 to NULL for faster divisibility checks
  // later on
  basis->cf[0]  = NULL;
  basis->eh[0]  = NULL;
  basis->load++;

  return basis;
}

inline void free_basis(gb_t *basis)
{
  if (basis) {
    nelts_t i;
    for (i=0; i<basis->nv; ++i) {
      free(basis->vnames[i]);
    }
    free(basis->vnames);
    free(basis->red);
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

hash_t add_new_element_to_basis_grevlex(gb_t *basis, const mat_t *mat,
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
  
  nelts_t ctr = 0;
  deg_t deg   = 0;
  nelts_t fc  = spd->col->nlm + mat->DR->row[ri]->piv_lead;
#if POLY_DEBUG
  printf("new lm: ");
#if !__GB_HAVE_SSE2
  for (int ii=0; ii<ht->nv; ++ii)
    printf("%u ",ht->exp[spd->col->hpos[fc]][ii]);
#endif
  printf(" %u  (%u)\n",ht->val[spd->col->hpos[fc]], spd->col->hpos[fc]);
#endif
  hash_t hv  = ht->val[spd->col->hpos[fc]];

  for (i=mat->DR->row[ri]->piv_lead; i<mat->DR->ncols; ++i) {
    if (mat->DR->row[ri]->piv_val[i] != 0) {
      basis->cf[basis->load][ctr] = mat->DR->row[ri]->piv_val[i];
      // note that we have to adjust the position via shifting it by
      // spd->col->nlm since DR is on the righthand side of the matrix
      basis->eh[basis->load][ctr] = spd->col->hpos[spd->col->nlm+i];
#if POLY_DEBUG
    printf("%u|",basis->cf[basis->load][ctr]);
#if !__GB_HAVE_SSE2
  for (int ii=0; ii<ht->nv; ++ii)
    printf("%u",ht->exp[basis->eh[basis->load][ctr]][ii]);
  printf("  ");
#endif
#endif
      deg = ht->deg[basis->eh[basis->load][ctr]] > deg ?
        ht->deg[basis->eh[basis->load][ctr]] : deg;
      ctr++;
    }
  }
#if POLY_DEBUG
  printf("\n");
  printf("deg: %u\n", deg);
#endif
  basis->nt[basis->load]  = ctr;
  basis->deg[basis->load] = deg;
  basis->red[basis->load] = NOT_REDUNDANT;


  // realloc memory to the correct number of terms
  basis->cf[basis->load]  = realloc(basis->cf[basis->load],
    basis->nt[basis->load] * sizeof(coeff_t));
  basis->eh[basis->load]  = realloc(basis->eh[basis->load],
      basis->nt[basis->load] * sizeof(hash_t));

  basis->load++;

  return hv;
}

inline void track_redundant_elements_in_basis(gb_t *basis)
{
  nelts_t i;
  // check for redundancy of other elements in basis
  for (i=basis->st; i<basis->load-2; ++i) {
    if (basis->red[i] == NOT_REDUNDANT) {
      if (check_monomial_division(basis->eh[i][0], basis->eh[basis->load-1][0], ht)) {
        basis->red[i] = REDUNDANT;
        basis->nred++;
      }
    }
  }

}

inline void enlarge_basis(gb_t *basis, const nelts_t size)
{
  basis->size = size;
  basis->nt   = realloc(basis->nt, basis->size * sizeof(nelts_t));
  basis->deg  = realloc(basis->deg, basis->size * sizeof(deg_t));
  basis->red  = realloc(basis->red, basis->size * sizeof(red_t));
  basis->cf   = realloc(basis->cf, basis->size * sizeof(coeff_t *));
  basis->eh   = realloc(basis->eh, basis->size * sizeof(hash_t *));
}
