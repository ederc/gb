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
 * \param selection set sel
 *
 * \param hash list of monomials mon
 */
void select_pairs_by_minimal_degree(ps_t *ps, gb_t *basis, sel_t *sel, pre_t *mon);

/**
 * \brief Enters the lower order monomials of the selected spair generators to
 * the preprocessing hash list.
 *
 * \param intermediate groebner basis basis
 *
 * \param selected list of elements to be used in the next round sel
 *
 * \param preprocessing hash list mon
 */
void enter_spairs_to_preprocessing_hash_list(const gb_t *basis, sel_t *sel, pre_t *mon);

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

/**
 * \brief Comparison function of monomials for quicksort. Compares the property
 * of being lead monomial or not.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
int cmp_symbolic_preprocessing_monomials_by_lead(const void *a, const void *b);

/**
 * \brief Comparison function of monomials for quicksort. Compares w.r.t. the
 * given monomial order grevlex.
 *
 * \note This is used in symbolic preprocessing, thus we know already that a =/=
 * b and we do not have to check this again.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
int cmp_symbolic_preprocessing_monomials_by_grevlex(const void *a, const void *b);

/**
 * \brief Sorts columns resp. monomials found by symbolic preprocessing to get
 * two different parts: Monomials which are lead monomials and monomials which
 * are non lead monomials.
 *
 * \param symbolic preprocessing data spd
 */
void sort_columns_by_lead(spd_t *spd);

/**
 * \brief Sorts lead monomials found by symbolic preprocessing w.r.t. the
 * given monomial order grevlex
 *
 * \note The list of columns resp. monomials was already presorted to
 * distinguish lead ones (at the beginning of the list) from non lead ones (at
 * the end of the list).
 *
 * \param symbolic preprocessing data spd
 */
void sort_lead_columns_by_grevlex(spd_t *spd);

/**
 * \brief Sorts non lead monomials found by symbolic preprocessing w.r.t. the
 * given monomial order grevlex
 *
 * \note The list of columns resp. monomials was already presorted to
 * distinguish lead ones (at the beginning of the list) from non lead ones (at
 * the end of the list).
 *
 * \param symbolic preprocessing data spd
 */
void sort_non_lead_columns_by_grevlex(spd_t *spd);

/**
 * \brief Sorts already presorted list of monomials resp. columns w.r.t. the
 * given monomial order grevlex. The list is already presorted by lead and non
 * lead monomials. Those monomials correspond then to the columns of the gbla
 * matrix generated later on
 *
 * \param symbolic preprocessing data spd
 *
 * \param number of threads to use in parallel nthreads
 */
void sort_presorted_columns_by_grevlex(spd_t *spd, int nthreads);
#endif
