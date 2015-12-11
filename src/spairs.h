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
 * \file spairs.h
 * \brief Implementation of handling of pair sets.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_SPAIRS_H
#define GB_SPAIRS_H

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

#define SPAIRS_DEBUG  0

/**
 * \brief Initialize pair set
 *
 * \param groebner basis gb
 *
 * \param hash table ht
 */
ps_t *init_pair_set(gb_t *basis, mp_cf4_ht_t *ht);

/**
 * \brief Enlarge pair set ps to size new_size
 *
 * \param pair set ps
 *
 * \param new size new_size
 */
void enlarge_pair_set(ps_t *ps, nelts_t new_size);

/**
 * \brief Frees dynamically allocated memory from pair set
 *
 * \param pair set ps
 */
void free_pair_set_dynamic_data(ps_t *ps);

/**
 * \brief Generates spair given by gen1 and gen2
 *
 * \param first generator gen1
 *
 * \param second generator gen2
 *
 * \param current groebner basis basis
 *
 * \param hash table ht
 *
 * \return generated spair
 */
spair_t *generate_spair(nelts_t gen1, nelts_t gen2, gb_t *basis, mp_cf4_ht_t *ht);

/**
 * \brief Comparison implementation for qsort. Sorts pair set w.r.t. graded
 * reverse lexicographical order grevlex.
 *
 * \param element to be compared a
 *
 * \param element to be compared b
 *
 * \returns corresponding integer for qsort
 */
int cmp_spairs_grevlex(const void *a, const void *b);

/**
 * \brief Sorts pair set w.r.t. graded reverse lexicographical order
 * grevlex using qsort.
 *
 * \param pair set to be sorted ps
 *
 */
void sort_pair_set_by_lcm_grevlex(ps_t *ps);

/**
 * \brief Updates pair set including Gebauer-Moeller criteria checks
 *
 * \param pair set ps
 *
 * \param intermediate groebner basis gb
 *
 * \param index of new element in gb idx
 */
void update_pair_set(ps_t *ps, gb_t *basis, nelts_t idx);

/**
 * \brief Gebauer-Moeller checks for product and chain criterion
 *
 * \param pair set ps
 *
 * \param hash position of newly added basis element hash
 *
 * \param index in basis of newly added basis element idx
 */
void gebauer_moeller(ps_t *ps, hash_t hash, nelts_t idx);

/**
 * \brief Remove spairs detected by either product or chain criterion
 *
 * \param pair set ps
 *
 * \param index of basis element the new pairs were generated with idx
 *
 * \return number of removed pairs
 */
nelts_t remove_detected_pairs(ps_t *ps, nelts_t idx);

/**
 * \brief Selects pairs by lowest degree (normal selection strategy) and returns
 * a set of this selection that is later on filled with corresponding lower term
 * reducers.
 *
 * \note The function also already marks the hash values for the lcms of the
 * spairs with "2". So we know that there are 2 polynomials in the upcoming
 * matrix that hit that monomial. We use this in the following symbolic
 * preprocessing not to handle this monomial again and also, in the matrix
 * construction later on, for know the number of nonzero entries in the column
 * corresponding to this monomial. This will help splicing the matrix for gbla.
 *
 * \param pair set ps
 *
 * \param intermediate grobner basis basis
 *
 * \return selection set sel
 */
sel_t *select_pairs_by_minimal_degree(ps_t *ps, gb_t *basis);

/**
 * \brief Adds generator gen of the corresponding spair with least common
 * multiple lcm to selection list sel.
 *
 * \param intermediate groebner basis basis
 *
 * \param selection list sel
 *
 * \param least common multiple of spair lcm
 *
 * \param a generator of spair gen
 */
void add_spair_generator_to_selection(gb_t *basis, sel_t *sel,
    const hash_t lcm, const nelts_t gen);

/**
 * \brief Checks storage for multipliers and enlarges storage for given entry if
 * necessary
 *
 * \param selection set sel
 *
 *  \param new size new_size
 *
 *  \param index of entry idx
 */
void check_enlargement_mul_in_selection(sel_t *sel, nelts_t new_size, nelts_t idx);

/**
 * \brief Enlarges selection set to new_size
 *
 * \param selection set sel
 *
 *  \param new size new_size
 */
void enlarge_selection(sel_t *sel, nelts_t new_size);

/**
 * \brief Initializes selection for next reduction step of size size
 *
 * \param size of spairs to generate selection out of size
 *
 * \return selection set sel
 */
sel_t *init_selection(nelts_t size);

/**
 * \brief Frees selection set
 *
 * \param selection set
 */
void free_selection(sel_t *sel);
#endif
