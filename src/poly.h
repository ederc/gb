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


/**
 * \brief Initializes basis by taking meta data from input elements.
 *
 * \note Input elements are not directly added to the basis, they are added as
 * spairs to the pair set and enter basis afterwards.
 *
 * \param input elements input
 *
 * \return intermediate groebner basis basis
 */
gb_t *initialize_basis(const gb_t *input);

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
 * \brief Frees dynamically allocated memory from simplifier list. Sets the list
 * to NULL.
 *
 * \param simplifier list sf
 */
void free_simplifier_list(gb_t *sf);

/**
 * \brief Frees dynamically allocated memory from groebner basis. Sets the basis
 * to NULL.
 *
 * \param groebner basis basis
 */
void free_basis(gb_t *basis);

/**
 * \brief Enlarges groebner basis to given new size.
 *
 * \param intermediate groebner basis basis
 *
 * \param new size size
 */
void enlarge_basis(gb_t *basis, nelts_t size);

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
void add_new_element_to_simplifier_list_grevlex(gb_t *sf, const dm_t *B,
    const nelts_t ri, const spd_t *spd, const mp_cf4_ht_t *ht);

/**
 * \brief Tracks and labels redundant elements in basis. In particular, we check
 * if the lead monomial of the last basis element (currently added to the basis)
 * divides the lead monomial of other elements in the basis. Those elements are
 * then labeled to be redundant.
 *
 * \param intermediate groebner basis basis
 */
void track_redundant_elements_in_basis(gb_t *basis);
#endif
