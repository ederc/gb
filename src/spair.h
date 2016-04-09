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
 * \file spair.h
 * \brief Implementation of handling of pair sets.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_SPAIR_H
#define GB_SPAIR_H

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

#ifndef SPAIR_DEBUG
#define SPAIR_DEBUG  0
#endif

/**
 * \brief Initialize pair set
 *
 * \param input elements input
 *
 * \param hash table ht
 */
ps_t *initialize_pair_set(const gb_t *input, mp_cf4_ht_t *ht);

/**
 * \brief Enters the input elements as sepcial spairs to the pair set:
 * Their second generator is zero. Thus, when handling these pairs we
 * automatically interreduce the input.
 *
 * \note We do not need to check for pair set enlargement since the number of
 * added pairs here is fixed and smaller than the initially allocated size.
 *
 * \param pair set ps
 *
 * \param input elements input
 */
void enter_input_elements_to_pair_set(ps_t *ps, const gb_t *input);

/**
 * \brief Enlarge pair set ps to size new_size
 *
 * \param pair set ps
 *
 * \param new size new_size
 */
void enlarge_pair_set(ps_t *ps, const nelts_t new_size);

/**
 * \brief Frees dynamically allocated memory from pair set
 *
 * \param pair set ps
 */
static inline void free_pair_set(ps_t **ps_in)
{
  ps_t *ps  = *ps_in;
  free(ps->pairs);
  free(ps);
  ps      = NULL;
  *ps_in  = ps;
}


/**
 * \brief Generates spair given by one input element. The first generator is 0.
 *
 * \note We use the second generator since then it is easier in symbolic
 * preprocessing: For all spairs, intial ones or normal ones, gen2 has to go to
 * selection sell (lower one). So not so many conditionals have to be done.
 *
 * \param second generator gen2
 *
 * \param input elements input
 *
 * \param hash table ht
 *
 * \return generated spair
 */
spair_t *generate_input_element_spair(const nelts_t gen1, const gb_t *input,
    mp_cf4_ht_t *ht);

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
spair_t *generate_spair(const nelts_t gen1, const nelts_t gen2, const gb_t *basis,
    mp_cf4_ht_t *ht);

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
void update_pair_set(ps_t *ps, const gb_t *basis, const nelts_t idx);

/**
 * \brief Gebauer-Moeller checks for product and chain criterion
 *
 * \param pair set ps
 *
 * \param intermediate groebner basis basis
 *
 * \param index in basis of newly added basis element idx
 */
void gebauer_moeller(ps_t *ps, const gb_t *basis,  const nelts_t idx);

/**
 * \brief Remove spairs detected by either product or chain criterion
 *
 * \param pair set ps
 *
 * \param intermediate groebner basis basis
 *
 * \param index of basis element the new pairs were generated with idx
 *
 * \return number of removed pairs
 */
nelts_t remove_detected_pairs(ps_t *ps, const gb_t *basis, const nelts_t idx);

/**
 * \brief Selects pairs by lowest degree (normal selection strategy) and returns
 * the index of the last pair in the pair list
 *
 * \note This function also sorts the pair set correspondingly.
 *
 * \param pair set ps
 *
 * \return last index of pair selection in pair list
 */
nelts_t get_pairs_by_minimal_degree(ps_t *ps);

/**
 * \brief Adds generator gen of the corresponding spair with least common
 * multiple lcm to selection list sel.
 *
 * \param selection list sel
 *
 * \param intermediate groebner basis basis
 *
 * \param least common multiple of spair lcm
 *
 * \param a generator of spair gen
 */
void add_spair_generator_to_selection(sel_t *sel, const gb_t *basis,
    const hash_t lcm, const nelts_t gen);

/**
 * \brief Adjusts selection set size to new_size
 *
 * \note This function is not only used for enlargements, but also
 * to cut down memory at the end of symbolic preprocessing
 *
 * \param selection set sel
 *
 *  \param new size new_size
 */
static inline void adjust_size_of_selection(sel_t *sel, const nelts_t new_size)
{
  sel->size = new_size;
  sel->mpp  = realloc(sel->mpp, sel->size * sizeof(mpp_t));
}

/**
 * \brief Initializes selection for next reduction step of size size
 *
 * \param size of spairs to generate selection out of size
 *
 * \return selection set sel
 */
static inline sel_t *init_selection(const nelts_t size)
{
  sel_t *sel  = (sel_t *)malloc(sizeof(sel_t));
  sel->size   = size;
  sel->load   = 0;

  sel->mpp    = (mpp_t *)malloc(size * sizeof(mpp_t));

  return sel;
}

/**
 * \brief Frees selection set
 *
 * \param selection set
 */
static inline void free_selection(sel_t **sel_in)
{
  sel_t *sel  = *sel_in;
  free(sel->mpp);
  free(sel);
  sel     = NULL;
  *sel_in = sel;
}
#endif
