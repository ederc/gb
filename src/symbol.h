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
 * \file symbol.h
 * \brief Implementation of the symbolic pre- and postprocessing in gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_SYMBOL_H
#define GB_SYMBOL_H

#include "gb_config.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include "types.h"
#include "hash.h"
#include "spairs.h"

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#define SYMBOL_DEBUG 0
#define __GB_SYM_LIST_LEN   1000

/**
 * \brief Symbolic preprocessing searching for all possible reducers of all
 * ocurring monomials. These data are then used to construct the matrices for
 * gbla.
 *
 * \param current pair set ps
 *
 * \param intermediate groebner basis basis
 *
 * \return full information from symbolic preprocessing for generationg corresponding
 * matrices out of polynomial data
 */
spd_t *symbolic_preprocessing(ps_t *ps, gb_t *basis);

/**
 * \brief Enters the lower order monomials of the selected spair generators to
 * the preprocessing hash list.
 *
 * \param selected list of elements to be used in the next round sel
 *
 * \param intermediate groebner basis basis
 *
 * \param preprocessing hash list mon
 */
void enter_spairs_to_preprocessing_hash_list(sel_t *sel, const gb_t *basis, pre_t *mon);

/**
 * \brief Enters one monomial (h1*h2) to preprocessing hash list.
 *
 * \param monomial hash position h1
 *
 * \param monomial hash position h2
 *
 * \param preprocessing hash list mon
 */
void enter_monomial_to_preprocessing_hash_list(const hash_t h1, const hash_t h2, pre_t *mon);

/**
 * \brief Initializes a hash list for symbolic preprocessing.
 *
 * \param size of hash list size
 *
 * \return hash list
 */
pre_t *init_preprocessing_hash_list(const nelts_t size);

/**
 * \brief Enlarges hash list for symbolic preprocessing to new_size.
 *
 * \param preprocessing hash list hl
 *
 * \param new size of hash list size
 */
void enlarge_preprocessing_hash_list(pre_t *hl, const nelts_t size);

/**
 * \brief Frees hash list for symbolic preprocessing.
 *
 * \param preprocessing hash list hl
 */
void free_preprocessing_hash_list(pre_t *hl);

/**
 * \brief Frees data structure for symbolic preprocessing.
 *
 * \param symbolic preprocessing data structure spd
 */
void free_symbolic_preprocessing_data(spd_t *spd);
#endif
