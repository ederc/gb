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
 * \brief Implementation of handling of polynomials and groebner bases.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_POLY_H
#define GB_POLY_H

#include "gb_config.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
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
 * \param switch for possible simplification simplify
 *
 * \param file length of input file fl
 *
 * \return intermediate groebner basis basis
 */
gb_t *initialize_basis(const int order, const int nlines,
    const nvars_t nvars, char **vnames, const mod_t mod,
    const int simplify, const uint64_t fl);

/**
 * \brief Initializes simplifier list, takes meta data from intermediate
 * groebner basis and allocates memory correspondingly for list entries.
 *
 * \param input elements input
 *
 * \return intermediate groebner basis basis
 */
gb_t *initialize_simplifier_list(const gb_t *basis);

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
 * \return hash value of new lead monomial
 */
hash_t add_new_element_to_basis_grevlex(gb_t *basis, const mat_t *mat,
    const nelts_t ri, const spd_t *spd, const mp_cf4_ht_t *ht);

/**
 * \brief Adds new element from reduced B part of gbla matrix to simplifier list.
 *
 * \note Also enlarges list if needed.
 *
 * \note Symbolic preprocessing data is needed to convert a matrix row back to a
 * polynomial with hashed exponents.
 *
 * \param intermediate groebner basis basis
 *
 * \param simplifier list sf
 *
 * \param part B of gbla matrix after reduction B
 *
 * \param row index in reduced D part ri
 *
 * \param symbolic preprocessing data spd
 *
 * \param hash table ht
 */
void add_new_element_to_simplifier_list_grevlex(gb_t *basis, gb_t *sf,
    const dm_t *B, const nelts_t ri, const spd_t *spd, const mp_cf4_ht_t *ht);

/**
 * \brief Frees dynamically allocated memory from simplifier list. Sets the list
 * to NULL.
 *
 * \param simplifier list sf
 */
static inline void free_simplifier_list(gb_t **sf_in)
{
  gb_t *sf  = *sf_in;
  if (sf) {
    nelts_t i;
    free(sf->vnames);
    free(sf->red);
    free(sf->deg);
    free(sf->nt);
    for (i=0; i<sf->load; ++i) {
      free(sf->cf[i]);
      free(sf->eh[i]);
    }
    free(sf->cf);
    free(sf->eh);
  }

  free(sf);
  sf      = NULL;
  *sf_in  = sf;
}


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

    if (basis->sf != NULL) {
      for (i=basis->st; i<basis->load; ++i) {
        free(basis->sf[i].idx);
      }
      free(basis->sf);
      basis->sf = NULL;
    }

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
  basis->nt   = realloc(basis->nt, basis->size * sizeof(nelts_t));
  basis->deg  = realloc(basis->deg, basis->size * sizeof(deg_t));
  basis->red  = realloc(basis->red, basis->size * sizeof(red_t));
  basis->cf   = realloc(basis->cf, basis->size * sizeof(coeff_t *));
  basis->eh   = realloc(basis->eh, basis->size * sizeof(hash_t *));
  if (basis->sf != NULL)
    basis->sf   = realloc(basis->sf, basis->size * sizeof(sf_t));
}

/**
 * \brief Initializes simplifier link between elements of intermediate groebner
 * basis and simplifier list.
 *
 * \param intermediate groebner basis basis
 */
static inline void initialize_simplifier_link(gb_t *basis)
{
  nelts_t i;
  //basis->sf = (sf_t *)malloc(basis->size * sizeof(sf_t));
  for (i=0; i<basis->size; ++i) {
    basis->sf[i].size = 3;
    basis->sf[i].load = 0;
    basis->sf[i].idx  = (nelts_t *)malloc(basis->sf[i].size * sizeof(nelts_t));
  }
}

/**
 * \brief Tracks and labels redundant elements in basis. In particular, we check
 * if the lead monomial of the last basis element (currently added to the basis)
 * divides the lead monomial of other elements in the basis. Those elements are
 * then labeled to be redundant.
 *
 * \param intermediate groebner basis basis
 */
static inline void track_redundant_elements_in_basis(gb_t *basis)
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

/**
 * \brief Connects simplifier entries for an element in the intermediate
 * groebner basis with entry in the global simplifier list.
 *
 * \param intermediate groebner basis basis
 *
 * \param simplifier list sf
 *
 * \param symbolic preprocessing data spd
 *
 * \param row index from gbla matrix ri
 */
static inline void link_simplifier_to_basis(gb_t *basis, const gb_t *sf,
    const spd_t *spd, const ri_t ri)
{
  // get index of basis element
  const nelts_t bi  = spd->selu->mpp[ri].bi;

  // enlarge array if needed
  if (basis->sf[bi].load  ==  basis->sf[bi].size) {
    basis->sf[bi].idx   =   realloc(basis->sf[bi].idx,
                              2 * basis->sf[bi].size * sizeof(nelts_t));
    basis->sf[bi].size  *=  2;
  }
  // insert link to simplifier list
  basis->sf[bi].idx[basis->sf[bi].load] = sf->load-1;
  basis->sf[bi].load++;
}
#endif
