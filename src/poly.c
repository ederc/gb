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
    const long max_spairs, const int64_t fl)
{
  nvars_t i;
  gb_t *basis = (gb_t *)malloc(sizeof(gb_t));

  basis->ord  = (ord_t)order;
  basis->nred = 0;
  /* per default homogeneous, easier to check inhomogeneity when reading input data */
  basis->init_hom = 1;
  basis->hom      = 1;
  basis->mtl      = 0;
  basis->max_sel  = (nelts_t)max_spairs;
  /* #generators of the input system = nlines - 2 since the first line has the
   * variable names and second line is the field modulus. Then we add 1 since we
   * keep the element at position 0 NULL for faster divisibility checks */
  basis->load     = (nelts_t)nlines -2 +1;
  basis->load_ls  = basis->load;
  basis->st       = basis->load;

  basis->has_unit = 0;

  basis->size   = 3 * basis->load;
  /* we add one extra variable space for possible homogenizations during the
   * computation */
  basis->rnv    = nvars;
  basis->nv     = nvars+1;
  basis->vnames = vnames;
  basis->mod    = mod;
  
  /* allocate memory for elements */
  basis->deg    = (deg_t *)malloc(basis->size * sizeof(deg_t));
  basis->red    = (red_t *)malloc(basis->size * sizeof(red_t));
  /* for the zero and the initial elements we set the redundancy value to zero */
  memset(basis->red, 0, basis->st * sizeof(red_t));
  basis->p      = (poly_t **)malloc(basis->size * sizeof(poly_t *));

  /* enter meta information from input file */
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
    /* calculates the maximal term length: here we take the longes variable
     * name, later we add 10 for coefficient and "+") */
    if (basis->mtl < strlen(basis->vnames[i]))
        basis->mtl  = strlen(basis->vnames[i]);
  }
  /* add to max. term length 5 for "^{exp}" */
  basis->mtl  +=  5;
  /* multiply max. term length by number of variables */
  basis->mtl  *=  basis->rnv;
  /* now add maximal coefficient length to mtl (max. term length) */
  basis->mtl  +=  COEFFICIENT_CHAR_LENGTH;
  return basis;
}


int add_new_element_to_basis(gb_t *basis, src_t *row,
    const pre_t *mon, const ht_t *ht)
{
  /* check if we have found a unit in the basis */
  if (ht->val[mon->hash[row[2]]] == 0)
    return 0;

  /* check next if this element might be redundant: this is only possible if
   * the input elements are not homogeneous. in this situation we might have
   * several new elements from D which have lead terms that divide each other.
   * if all polynomials are homogeneous this cannot happen since then such a
   * lead term divisibility must have been found already in the linear algebra
   * reduction process. */
  if (basis->hom == 0 &&
      check_new_element_for_redundancy(mon->hash[row[2]], basis, ht) != 0) {
    return -1;
  }

  /* if not redundandant */
  deg_t deg = 0;
  for (size_t i = 2; i < row[1]; i += 2) {
    /* change column position to corresponding exponent hash */
    row[i]  = mon->hash[row[i]];

    /* compute degree of polynomial */
    deg = ht->deg[row[i]] > deg ?
      ht->deg[row[i]] : deg;
  }
  basis->deg[basis->load] = deg;
  basis->p[basis->load]   = row;
  basis->red[basis->load] = 0;

#if POLY_DEBUG
  printf("new element (%u) ", basis->load);
  for (size_t i = 2; i < basis->p[basis->load][1]; i = i+2)
    printf("%u %u | ", basis->p[basis->load][i+1], basis->p[basis->load][i]);
  printf("\n")
#endif

  basis->load++;
  return 1;
}

int add_new_element_to_basis_all_pivs(gb_t *basis, src_t *row,
    const pre_t *mon, const ht_t *ht)
{
  /* check if we have found a unit in the basis */
  if (ht->val[mon->hash[row[2]]] == 0)
    return 0;

  /* check next if this element might be redundant: this is only possible if
   * the input elements are not homogeneous. in this situation we might have
   * several new elements from D which have lead terms that divide each other.
   * if all polynomials are homogeneous this cannot happen since then such a
   * lead term divisibility must have been found already in the linear algebra
   * reduction process. */
  if (basis->hom == 0 &&
      check_new_element_for_redundancy(mon->hash[row[2]], basis, ht) != 0) {
    return -1;
  }

  /* if not redundandant */
  deg_t deg = 0;
  for (size_t i = 2; i < row[1]; i += 2) {
    /* change column position to corresponding exponent hash */
    row[i]  = mon->hash[row[i]];

    /* compute degree of polynomial */
    deg = ht->deg[row[i]] > deg ?
      ht->deg[row[i]] : deg;
  }
  basis->deg[basis->load] = deg;
  basis->p[basis->load]   = row;
  basis->red[basis->load] = 0;

  basis->load++;

  return 1;
}
