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
inline gb_t *initialize_basis(const int ordering, const int nlines,
    const nvars_t nvars, char **vnames, const mod_t mod,
    const int simplify, const uint64_t fl)
{
  nvars_t i;
  gb_t *basis = (gb_t *)malloc(sizeof(gb_t));

  basis->ord  = (ord_t)ordering;
  basis->nred = 0;
  basis->hom  = 0;
  basis->mtl  = 0;
  // #generators of the input system = nlines - 2 since the first line has the
  // variable names and second line is the field modulus. Then we add 1 since we
  // keep the element at position 0 NULL for faster divisibility checks
  basis->load = nlines -2 +1;
  basis->st   = basis->load;

  basis->size   = 3 * basis->load;
  basis->nv     = nvars;
  basis->vnames = vnames;
  basis->mod    = mod;
  
  // allocate memory for elements
  basis->nt     = (nelts_t *)malloc(basis->size * sizeof(nelts_t));
  basis->deg    = (deg_t *)malloc(basis->size * sizeof(deg_t));
  basis->red    = (red_t *)malloc(basis->size * sizeof(red_t));
  basis->cf     = (coeff_t **)malloc(basis->size * sizeof(coeff_t *));
  basis->eh     = (hash_t **)malloc(basis->size * sizeof(hash_t *));

  if (simplify == 1) {
    basis->sf = (sf_t *)malloc(basis->size * sizeof(sf_t));
    // initialize dummy values for input polynomials
    for (i=0; i<basis->st; ++i) {
      basis->sf[i].load = basis->sf[i].size = 0;
      basis->sf[i].mul  = NULL;
      basis->sf[i].idx  = NULL;
    }
  } else {
    basis->sf = NULL;
  }

  // enter meta information from input file
  basis->fs   = (double) fl / 1024;
  basis->fsu  = "KB";
  if (basis->fs > 1000) {
    basis->fs   = basis->fs / 1024;
    basis->fsu  = "MB";
  }
  if (basis->fs > 1000) {
    basis->fs   = basis->fs / 1024;
    basis->fsu  = "GB";
  }
  for (i=0; i<basis->nv; ++i) {
    // calculates the maximal term length: here we take the longes variable
    // name, later we add 10 for coefficient and "+")
    if (basis->mtl < strlen(basis->vnames[i]))
        basis->mtl  = strlen(basis->vnames[i]);
  }
  // add to max. term length 5 for "^{exp}"
  basis->mtl  +=  5;
  // multiply max. term length by number of variables
  basis->mtl  *=  basis->nv;
  // now add maximal coefficient length to mtl (max. term length)
  basis->mtl  +=  COEFFICIENT_CHAR_LENGTH;
  return basis;
}

inline gb_t *initialize_simplifier_list(const gb_t *basis)
{
  gb_t *sf = (gb_t *)malloc(sizeof(gb_t));

  // get meta data from input
  sf->size    = basis->size; // TODO: This is too small, what is a good factor?
  sf->load    = 0;
  sf->nv      = basis->nv;
  sf->mod     = basis->mod;
  sf->vnames  = NULL;

  // allocate memory for elements
  sf->nt  = (nelts_t *)malloc(sf->size * sizeof(nelts_t));
  sf->deg = (deg_t *)malloc(sf->size * sizeof(deg_t));
  sf->red = (red_t *)malloc(sf->size * sizeof(red_t));
  sf->cf  = (coeff_t **)malloc(sf->size * sizeof(coeff_t *));
  sf->eh  = (hash_t **)malloc(sf->size * sizeof(hash_t *));
  sf->sf  = NULL;

  // initialize element at position 0 to NULL for faster divisibility checks
  // later on
  sf->cf[0]  = NULL;
  sf->eh[0]  = NULL;
  sf->load++;

  return sf;
}

void add_new_element_to_simplifier_list_grevlex(gb_t *basis, gb_t *sf,
    const dm_t *B, const nelts_t ri, const spd_t *spd, const mp_cf4_ht_t *ht)
{
  nelts_t i;

  if (sf->load == sf->size)
    enlarge_basis(sf, 2*sf->size);

  // maximal size is B->ncols + 1 (for A)
  nelts_t ms  = B->ncols + 1;
  // use shorter names in here
  sf->cf[sf->load] = (coeff_t *)malloc(ms * sizeof(coeff_t));
  sf->eh[sf->load] = (hash_t *)malloc(ms * sizeof(hash_t));

  nelts_t ctr = 0;
  deg_t deg   = 0;

  // first we add the term from A
  sf->cf[sf->load][ctr] = 1;
  sf->eh[sf->load][ctr] = spd->col->hpos[ri];
  deg = ht->deg[sf->eh[sf->load][ctr]] > deg ?
    ht->deg[sf->eh[sf->load][ctr]] : deg;
  ctr++;

#if POLY_DEBUG
  printf("new simplifier lm: ");
#if !__GB_HAVE_SSE2
  for (int ii=0; ii<ht->nv; ++ii)
    printf("%u ",ht->exp[spd->col->hpos[ri]][ii]);
#endif
  printf(" %u  (%u)\n",ht->val[spd->col->hpos[ri]], spd->col->hpos[ri]);
#endif
  // now we do B
  for (i=0; i<B->ncols; ++i) {
    if (B->row[ri]->init_val[i] != 0) {
      sf->cf[sf->load][ctr] = B->row[ri]->init_val[i];
      // note that we have to adjust the position via shifting it by
      // spd->col->nlm since DR is on the righthand side of the matrix
      sf->eh[sf->load][ctr] = spd->col->hpos[spd->col->nlm+i];
#if POLY_DEBUG
    printf("%u|",sf->cf[sf->load][ctr]);
#if !__GB_HAVE_SSE2
  for (int ii=0; ii<ht->nv; ++ii)
    printf("%u",ht->exp[sf->eh[sf->load][ctr]][ii]);
  printf("  ");
#endif
#endif
      deg = ht->deg[sf->eh[sf->load][ctr]] > deg ?
        ht->deg[sf->eh[sf->load][ctr]] : deg;
      ctr++;
    }
  }
#if POLY_DEBUG
  printf("\n");
  printf("deg: %u\n", deg);
#endif
  sf->nt[sf->load]  = ctr;
  sf->deg[sf->load] = deg;


  // realloc memory to the correct number of terms
  sf->cf[sf->load]  = realloc(sf->cf[sf->load],
    sf->nt[sf->load] * sizeof(coeff_t));
  sf->eh[sf->load]  = realloc(sf->eh[sf->load],
      sf->nt[sf->load] * sizeof(hash_t));

#if POLY_DEBUG
  printf("new simplifier: %u\n",ri);
  for (int ii=0; ii<sf->nt[sf->load]; ++ii)
    printf("%lu ", sf->eh[sf->load][ii]);
  printf("\n");
#endif
  sf->load++;
  // now we have to link the simplifier with the corresponding element in basis
  
  // get index of element in basis
  link_simplifier_to_basis(basis, sf, spd, ri);
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

  if (basis->sf != NULL) {
    basis->sf[basis->load].size = 3;
    basis->sf[basis->load].load = 0;
    basis->sf[basis->load].mul  = (hash_t *)malloc(basis->sf[basis->load].size * sizeof(hash_t));
    basis->sf[basis->load].idx  = (nelts_t *)malloc(basis->sf[basis->load].size * sizeof(nelts_t));
  }
  basis->load++;

  return hv;
}
