/* gb: Gröbner Basis
 * Copyright (C) 2018 Christian Eder <ederc@mathematik.uni-kl.de>
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
 * \brief Implementation of handling of polynomials and groebner bases.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_POLY_H
#define GB_POLY_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <config.h>
#include "types.h"
#include "hash.h"

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#ifndef POLY_DEBUG
#define POLY_DEBUG 0
#endif

#define COEFFICIENT_CHAR_LENGTH 10


/**
 * \brief Initializes basis by taking meta data from input elements.
 *
 * \note Input elements are not directly added to the basis, they are added as
 * spairs to the pair set and enter basis afterwards.
 *
 * \param given monomial order order
 *
 * \param number of lines of input files nlines
 *
 * \param number of variables nvars
 *
 * \param variable names vnames
 *
 * \param modulo resp. field characteristic mod
 *
 * \param maximal number of spairs handled at once
 *
 * \param file length of input file fl
 *
 * \return intermediate groebner basis basis
 */
gb_t *initialize_basis(const int order, const int nlines,
    const nvars_t nvars, char **vnames, const mod_t mod,
    const long max_spairs, const int64_t fl);


/**
 * \brief Adds new element from reduced D part of gbla matrix to basis.
 *
 * \note Also enlarges basis if needed.
 *
 * \note Symbolic preprocessing data is needed to convert a matrix row back to a
 * polynomial with hashed exponents.
 *
 * \param intermediate groebner basis basis
 *
 * \param reduced gbla matrix mat
 *
 * \param row index in reduced D part ri
 *
 * \param symbolic preprocessing data spd
 *
 * \param hash table ht
 *
 * \return integer that indicates the following situation:
 * 0:   the new element would be a constant
 * 1:   the new element is a added to the basis
 * -1:  the new element is not added to the basis since it is redundant
 */
int add_new_element_to_basis(gb_t *basis, src_t *row,
    const pre_t *mon, const ht_t *ht);

/**
 * \brief Frees dynamically allocated memory from groebner basis. Sets the basis
 * to NULL.
 *
 * \param groebner basis basis
 */
static inline void free_basis(gb_t **basis_in)
{
  gb_t *basis = *basis_in;
  if (basis) {
    nelts_t i;

    for (i=0; i<basis->rnv; ++i) {
      free(basis->vnames[i]);
    }
    free(basis->vnames);
    free(basis->red);
    free(basis->deg);
    for (i=0; i<basis->load; ++i) {
      free(basis->p[i]);
    }
    free(basis->p);
  }

  free(basis);
  basis     = NULL;
  *basis_in = basis;
}

/**
 * \brief Enlarges groebner basis to given new size.
 *
 * \param intermediate groebner basis basis
 *
 * \param new size size
 */
static inline void enlarge_basis(gb_t *basis, const nelts_t size)
{
  basis->size = size;
  basis->deg  = realloc(basis->deg, basis->size * sizeof(deg_t));
  basis->red  = realloc(basis->red, basis->size * sizeof(red_t));
  basis->p    = realloc(basis->p, basis->size * sizeof(poly_t *));
}

/**
 * \brief Tracks and labels redundant elements in basis. In particular, we check
 * if the lead monomial of the last basis element (currently added to the basis)
 * divides the lead monomial of other elements in the basis. Those elements are
 * then labeled to be redundant.
 *
 * \param intermediate groebner basis basis
 */
static inline void track_redundant_elements_in_basis(gb_t *basis, const ht_t *ht)
{
  nelts_t i;
  /* check for redundancy of other elements in basis */
  for (i=basis->st; i<basis->load-1; ++i) {
    if (basis->red[i] == 0) {
      if (check_monomial_division(
            basis->p[i][2], basis->p[basis->load-1][2], ht)) {
        basis->red[i] = basis->load-1;
        basis->nred++;
      }
    }
  }
}

/**
 * \brief Checks if the new element is probably already redundant. This is
 * possible if elements of lower degree were added to basis in the same
 * reduction step and one of these elements has a lead term dividing the lead
 * term of this new, higher-degree element.
 *
 * \note We only have to check with the elements already added in this step of
 * the algorithm. Elements that have been in the basis earlier must have reduced
 * this lead term.
 *
 * \param hash of new lead term hash
 *
 * \param intermediate groebner basis basis
 *
 * \return 0 if not redundant, =/= 0 else
 */
static inline int check_new_element_for_redundancy(hash_t hash, const gb_t *basis,
    const ht_t *ht)
{
  /* check for redundancy of other elements in basis */
  for (nelts_t i=basis->load_ls; i<basis->load; ++i) {
    if (basis->red[i] == 0 && check_monomial_division(hash, basis->p[i][2], ht) == 1)
      return 1;
  }
  return 0;
}
#endif
