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
inline gb_t *initialize_basis(const int order, const int nlines,
    const nvars_t nvars, char **vnames, const mod_t mod,
    const int simplify, const uint64_t fl)
{
  nvars_t i;
  gb_t *basis = (gb_t *)malloc(sizeof(gb_t));

  basis->ord  = (ord_t)order;
  basis->nred = 0;
  // per default homogeneous, easier to check inhomogeneity when reading input data
  basis->init_hom = 1;
  basis->hom      = 1;
  basis->mtl      = 0;
  // #generators of the input system = nlines - 2 since the first line has the
  // variable names and second line is the field modulus. Then we add 1 since we
  // keep the element at position 0 NULL for faster divisibility checks
  basis->load     = nlines -2 +1;
  basis->load_ls  = basis->load;
  basis->st       = basis->load;

  basis->has_unit = 0;

  basis->size   = 3 * basis->load;
  // we add one extra variable space for possible homogenizations during the
  // computation
  basis->rnv    = nvars;
  basis->nv     = nvars+1;
  basis->vnames = vnames;
  basis->mod    = mod;
  
  // allocate memory for elements
  basis->p  = (poly_t **)malloc(basis->size * sizeof(poly_t *));
  basis->sl = simplify;

  if (basis->sl > 0) {
    basis->sf = (sf_t *)malloc(basis->size * sizeof(sf_t));
    // initialize dummy values for input polynomials
    for (i=0; i<basis->st; ++i) {
      basis->sf[i].load = basis->sf[i].size = 0;
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
  for (i=0; i<basis->rnv; ++i) {
    // calculates the maximal term length: here we take the longes variable
    // name, later we add 10 for coefficient and "+")
    if (basis->mtl < strlen(basis->vnames[i]))
        basis->mtl  = strlen(basis->vnames[i]);
  }
  // add to max. term length 5 for "^{exp}"
  basis->mtl  +=  5;
  // multiply max. term length by number of variables
  basis->mtl  *=  basis->rnv;
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
  sf->p   = (poly_t **)malloc(sf->size * sizeof(poly_t *));
  sf->sf  = NULL;

  // initialize element at position 0 to NULL for faster divisibility checks
  // later on
  sf->p[0]      = (poly_t *)malloc(sizeof(poly_t));
  sf->p[0]->cf  = NULL;
  sf->p[0]->eh  = NULL;
  sf->load++;

  return sf;
}

void add_new_element_to_simplifier_list(gb_t *basis, gb_t *sf,
    const dm_t *B, const nelts_t ri, const spd_t *spd, const mp_cf4_ht_t *ht)
{
  nelts_t i;

  if (sf->load == sf->size)
    enlarge_basis(sf, 2*sf->size);

  // maximal size is B->ncols + 1 (for A)
  nelts_t ms  = B->ncols + 1;

  // allocate memory for new poly
  sf->p[sf->load]   = (poly_t *)malloc(sizeof(poly_t));
  poly_t * const p  = sf->p[sf->load];
  // use shorter names in here
  p->cf = (coeff_t *)malloc(ms * sizeof(coeff_t));
  p->eh = (hash_t *)malloc(ms * sizeof(hash_t));

  nelts_t ctr = 0;
  deg_t deg   = 0;

  // first we add the term from A
  p->cf[ctr] = 1;
  p->eh[ctr] = spd->col->hpos[ri];
  deg = ht->deg[p->eh[ctr]] > deg ?
    ht->deg[p->eh[ctr]] : deg;
  ctr++;

#if POLY_DEBUG
  printf("row %u -- new simplifier lm: ", ri);
#if !__GB_HAVE_SSE2
  for (int ii=0; ii<ht->nv; ++ii)
   // printf("%u ",ht->exp[spd->col->hpos[ri]][ii]);
    printf("%u ",ht->exp[sf->eh[sf->load][0]][ii]);
#else
    exp_t expa[ht->nev * ht->vl] __attribute__ ((aligned (16)));
    exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
    for (int ii=0; ii<ht->nev; ++ii) {
      _mm_store_si128((exp_v *)tmp, ht->ev[sf->eh[sf->load][0]][ii]);
      memcpy(expa+(ii*ht->vl), tmp, ht->vl*sizeof(exp_t));
    }
  for (int ii=0; ii<ht->nv; ++ii)
   // printf("%u ",ht->exp[spd->col->hpos[ri]][ii]);
    printf("%u ",expa[ii]);
#endif
  printf(" %u  (%u)\n",ht->val[spd->col->hpos[ri]], spd->col->hpos[ri]);
  printf("for the basis lm:  ");
#if !__GB_HAVE_SSE2
  printf("ri %u -- bi %u\n", ri, spd->selu->mpp[ri].bi);
  for (int ii=0; ii<ht->nv; ++ii)
    printf("%u ",ht->exp[basis->eh[spd->selu->mpp[ri].bi][0]][ii]);
#else
    exp_t expa[ht->nev * ht->vl] __attribute__ ((aligned (16)));
    exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
    for (int ii=0; ii<ht->nev; ++ii) {
      _mm_store_si128((exp_v *)tmp, ht->ev[basis->eh[spd->selu->mpp[ri].bi][0]][ii]);
      memcpy(expa+(ii*ht->vl), tmp, ht->vl*sizeof(exp_t));
    }
  for (int ii=0; ii<ht->nv; ++ii)
   // printf("%u ",ht->exp[spd->col->hpos[ri]][ii]);
    printf("%u ",expa[ii]);
#endif
  printf("\n");
#endif
  // now we do B
  //
  re_t *values;
  nelts_t si  = 0;
  if (basis->sl == 1) {
    values  = B->row[ri]->init_val;
  } else {
    values  = B->row[ri]->piv_val;
    si      = B->row[ri]->piv_lead;
  }
  if (values != NULL) {
    for (i=si; i<B->ncols; ++i) {
      if (values[i] != 0) {
        p->cf[ctr] = values[i];
        // note that we have to adjust the position via shifting it by
        // spd->col->nlm since DR is on the righthand side of the matrix
        p->eh[ctr] = spd->col->hpos[spd->col->nlm+i];
#if POLY_DEBUG
        printf("%u|",p->cf[ctr]);
#if !__GB_HAVE_SSE2
        for (int ii=0; ii<ht->nv; ++ii)
          printf("%u",ht->exp[p->eh[ctr]][ii]);
        printf("  ");
#endif
#endif
        deg = ht->deg[p->eh[ctr]] > deg ?
          ht->deg[p->eh[ctr]] : deg;
        ctr++;
      }
    }
  }
#if POLY_DEBUG
  printf("\n");
  printf("deg: %u\n", deg);
#endif
  p->nt   = ctr;
  p->deg  = deg;


  // realloc memory to the correct number of terms
  p->cf  = realloc(p->cf, p->nt * sizeof(coeff_t));
  p->eh  = realloc(p->eh, p->nt * sizeof(hash_t));

#if POLY_DEBUG
  printf("new simplifier: %u\n",ri);
  for (int ii=0; ii<p->nt; ++ii)
    printf("%lu ", p->eh[ii]);
  printf("\n");
#endif
  sf->load++;
  // now we have to link the simplifier with the corresponding element in basis
  
  // get index of element in basis
  link_simplifier_to_basis(basis, sf, spd, ri);
}

int add_new_element_to_basis(gb_t *basis, const mat_t *mat,
    const nelts_t ri, const spd_t *spd, const mp_cf4_ht_t *ht)
{
  // get position of lead term in this row
  const nelts_t fc  = spd->col->nlm + mat->DR->row[ri]->piv_lead;

  // check if we have found a unit in the basis
  hash_t hv  = ht->val[spd->col->hpos[fc]];
  if (hv == 0)
    return 0;

  // check next if this element might be redundant: this is only possible if
  // the input elements are not homogeneous. in this situation we might have
  // several new elements from D which have lead terms that divide each other.
  // if all polynomials are homogeneous this cannot happen since then such a
  // lead term divisibility must have been found already in the linear algebra
  // reduction process.
  if (basis->hom == 0 &&
      check_new_element_for_redundancy(spd->col->hpos[fc], basis) != 0) {
    return -1;
  }

  // if not redundandant
  nelts_t i;

  if (basis->load == basis->size)
    enlarge_basis(basis, 2*basis->size);

  // maximal size is DR->ncols - row's piv lead
  nelts_t ms  = mat->DR->ncols - mat->DR->row[ri]->piv_lead;

  // allocate memory for new poly
  basis->p[basis->load] = (poly_t *)malloc(sizeof(poly_t));
  poly_t * const p      = basis->p[basis->load];
  // use shorter names in here
  p->cf = (coeff_t *)malloc(ms * sizeof(coeff_t));
  p->eh = (hash_t *)malloc(ms * sizeof(hash_t));

  nelts_t ctr = 0;
  deg_t deg   = 0;
#if POLY_DEBUG
  printf("new lm from row %u (basis element %u): ", ri, basis->load);
#if !__GB_HAVE_SSE2
  for (int ii=0; ii<ht->nv; ++ii)
    printf("%u ",ht->exp[spd->col->hpos[fc]][ii]);
#else
    exp_t expa[ht->nev * ht->vl] __attribute__ ((aligned (16)));
    exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
    for (int ii=0; ii<ht->nev; ++ii) {
      _mm_store_si128((exp_v *)tmp, ht->ev[spd->col->hpos[fc]][ii]);
      memcpy(expa+(ii*ht->vl), tmp, ht->vl*sizeof(exp_t));
    }
  for (int ii=0; ii<ht->nv; ++ii)
   // printf("%u ",ht->exp[spd->col->hpos[ri]][ii]);
    printf("%u ",expa[ii]);
#endif
  printf(" %u  (%u)\n",ht->val[spd->col->hpos[fc]], spd->col->hpos[fc]);
#endif

  for (i=mat->DR->row[ri]->piv_lead; i<mat->DR->ncols; ++i) {
    if (mat->DR->row[ri]->piv_val[i] != 0) {
      p->cf[ctr]  = mat->DR->row[ri]->piv_val[i];
      // note that we have to adjust the position via shifting it by
      // spd->col->nlm since DR is on the righthand side of the matrix
      p->eh[ctr] = spd->col->hpos[spd->col->nlm+i];
#if POLY_DEBUG
    printf("%u|",basis->p[basis->load].cf[ctr]);
#if !__GB_HAVE_SSE2
  for (int ii=0; ii<ht->nv; ++ii)
    printf("%u",ht->exp[basis->eh[basis->load][ctr]][ii]);
  printf("  ");
#else
    exp_t expa[ht->nev * ht->vl] __attribute__ ((aligned (16)));
    exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
    for (int ii=0; ii<ht->nev; ++ii) {
      _mm_store_si128((exp_v *)tmp, ht->ev[basis->eh[basis->load][ctr]][ii]);
      memcpy(expa+(ii*ht->vl), tmp, ht->vl*sizeof(exp_t));
    }
  for (int ii=0; ii<ht->nv; ++ii)
   // printf("%u ",ht->exp[spd->col->hpos[ri]][ii]);
    printf("%u",expa[ii]);
  printf("  ");
#endif
#endif
      deg = ht->deg[p->eh[ctr]] > deg ? ht->deg[p->eh[ctr]] : deg;
      ctr++;
    }
  }
#if POLY_DEBUG
  printf("\n");
  printf("deg: %u\n", deg);
  printf("# terms = 1 + %u\n",ctr-1);
#endif
  p->nt  = ctr;
  p->deg = deg;
  p->red = 0;


  // realloc memory to the correct number of terms
  p->cf = realloc(p->cf, p->nt * sizeof(coeff_t));
  p->eh = realloc(p->eh, p->nt * sizeof(hash_t));

  if (basis->sf != NULL) {
    basis->sf[basis->load].size = 3;
    basis->sf[basis->load].load = 0;
    basis->sf[basis->load].idx  = (nelts_t *)malloc(basis->sf[basis->load].size * sizeof(nelts_t));
  }
  basis->load++;

  return 1;
}
