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
#include "spair.h"

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#ifndef SYMBOL_DEBUG
#define SYMBOL_DEBUG 0
#endif

/**
 * \brief Symbolic preprocessing searching for all possible reducers of all
 * ocurring monomials. These data are then used to construct the matrices for
 * gbla.
 *
 * \param current pair set ps
 *
 * \param intermediate groebner basis basis
 *
 * \param simplifier list sf
 *
 * \return full information from symbolic preprocessing for generationg corresponding
 * matrices out of polynomial data
 */
spd_t *symbolic_preprocessing(ps_t *ps, const gb_t *basis, const gb_t *sf);

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
 * \brief Enters one monomial h1 to preprocessing hash list.
 *
 * \param monomial hash position h1
 *
 * \param preprocessing hash list mon
 */
void enter_not_multiplied_monomial_to_preprocessing_hash_list(const hash_t h1,
    pre_t *mon);

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
 * \brief Adjusts size of hash list for symbolic preprocessing to new_size.
 *
 * \note It is not only used to enlarge the size, but also for cutting down
 * memory once symbolic preprocessing is done.
 *
 * \param preprocessing hash list hl
 *
 * \param new size of hash list size
 */
void adjust_size_of_preprocessing_hash_list(pre_t *hl, const nelts_t size);

/**
 * \brief Frees hash list for symbolic preprocessing.
 *
 * \param preprocessing hash list hl
 */
static inline void free_preprocessing_hash_list(pre_t **hl_in)
{
  pre_t *hl = *hl_in;
  free(hl->hpos);
  free(hl);
  hl      = NULL;
  *hl_in  = hl;
}


/**
 * \brief Frees data structure for symbolic preprocessing.
 *
 * \param symbolic preprocessing data structure spd
 */
static inline void free_symbolic_preprocessing_data(spd_t **spd_in)
{
  spd_t *spd = *spd_in;
  free_preprocessing_hash_list(&(spd->col));
  free_selection(&(spd->selu));
  free_selection(&(spd->sell));
  free(spd);
  spd     = NULL;
  *spd_in = spd;
}


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
 * \brief Comparison function of monomials for quicksort. Compares w.r.t. the
 * given monomial order grevlex and inverts the order.
 *
 * \note This is used in symbolic preprocessing, thus we know already that a =/=
 * b and we do not have to check this again.
 *
 * \note We invert grevlex since we have to align the A and C parts of the gbla
 * matrix in order to efficiently compute.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
int cmp_symbolic_preprocessing_monomials_by_inverse_grevlex(const void *a, const void *b);

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
 * \brief Sorts lead monomials found by symbolic preprocessing w.r.t. the
 * given monomial order grevlex and inverts the order
 *
 * \note The list of columns resp. monomials was already presorted to
 * distinguish lead ones (at the beginning of the list) from non lead ones (at
 * the end of the list).
 *
 * \note We invert grevlex since we have to align the A and C parts of the gbla
 * matrix in order to efficiently compute.
 *
 * \param symbolic preprocessing data spd
 */
void sort_lead_columns_by_inverse_grevlex(spd_t *spd);

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
void sort_presorted_columns_by_grevlex(spd_t *spd, const int nthreads);

/**
 * \brief Sets the index entry of the hash values to the corresponding columns
 * for the upcoming gbla matrix. The relation between hash values of monomials
 * and columns is given by the symbolic preprocessing data.
 *
 * \param hash table ht
 *
 * \param symbolic preprocessing data spd
 */
void set_column_index_in_hash_table(mp_cf4_ht_t *ht, const spd_t *spd);

/**
 * \brief Implements the comparison function for quicksort used in the function
 * sort_selection_by_column_index(). Takes multiplied lead monomials and sorts
 * corresponding to the predefined column index that is stored in the idx entry
 * of the hash table.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns ht->idx[a.mlm] - ht->idx[b.mlm]
 */
int cmp_monomial_polynomial_pair(const void *a, const void *b);

/**
 * \brief Sorts upper and lower selection of polynomials from preprocessing to
 * by the position of the corresponding multiplied lead monomial w.r.t. to the
 * predefined column order of the gbla matrix to be generated next.
 *
 * \note Both lists, upper and lower, can be sorted in parallel, thus the
 * implementation is done using open mp tasks.
 *
 * \param symbolic data structure spd
 *
 * \param hash table ht
 *
 * \param number of threads nthreads
 */
void sort_selection_by_column_index(spd_t *spd, const mp_cf4_ht_t *ht,
  const int nthreads);

/**
 * \brief Selects pairs due to the sorting and selection done by the previous
 * get_pairs_by_*() procedure and returns
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
 * \note There are two selection lists, one for the upper and one for the lower
 * part of the gbla matrix to be constructed. By default, the first generator of
 * an spair is added to the upper part, the second generator to the lower part.
 * All upcoming reducers are also added to the upper selection list.
 *
 * \param pair set ps
 *
 * \param selection set for upper part of gbla matrix sel_upp
 *
 * \param selection set for lower part of gbla matrix sel_low
 *
 * \param hash list of monomials mon
 *
 * \param intermediate grobner basis basis
 *
 * \param simplifier list sf
 *
 * \param last index of pair selection in pair list idx
 */
static inline void select_pairs(ps_t *ps, sel_t *selu, sel_t *sell, pre_t *mon,
    const gb_t *basis, const gb_t *sf, const nelts_t nsel)
{
  nelts_t i, j, k;
  spair_t *sp;

  hash_t mul  = 0;
  nelts_t nt  = 0;
  hash_t *pol = NULL;

  nelts_t have_sf = 0;

  // wVe do not need to check for size problems in sel du to above comment: we
  // have allocated basis->load slots, so enough for each possible element from
  // the basis
#if SYMBOL_DEBUG
  printf("selected pairs in this step of the algorithm:\n");
#endif
  for (i=0; i<nsel; ++i) {
    sp = ps->pairs[i];
#if SYMBOL_DEBUG
    printf("gen1 %u -- gen2 %u -- lcm %u\n", sp->gen1, sp->gen2, sp->lcm);
#endif
    // We have to distinguish between usual spairs and spairs consisting of one
    // initial input element: The latter ones have only 1 generator (gen2 has
    // the generator, gen1 is zero) and the
    // corresponding lead monomial is not part of the basis. Thus we cannot set
    // ht->idx[*] = 2 for these. The whole data for these special spairs is in
    // basis before basis->st, whereas the data of the usual spairs is
    // completely in basis starting at position basis->st.

    // so we can always put gen2 in the lower selection list sell
    add_spair_generator_to_selection(sell, basis, sp->lcm, sp->gen2);
    j = sell->load-1;
    /*
    if (basis->sf != NULL) {
      nelts_t l     = 0;
      hash_t sf_mul = 0;
      const nelts_t load  = basis->sf[sell->mpp[j].idx].load;
      for (l=0; l<load; ++l) {
        // we start searching from the end of the list since those elements
        // might be best reduced
        //sf_mul = monomial_division(sell->mpp[j].mul, basis->sf[sell->mpp[j].idx].mul[load-1-l], ht);
        sf_mul = monomial_division(sell->mpp[j].mlm, sf->eh[basis->sf[sell->mpp[j].idx].idx[load-1-l]][0], ht);
        if (sf_mul != 0) {
          printf("-- SIMPLIFY --\n");
          for (int ii=0; ii<basis->nv; ++ii)
            printf("%u ",ht->exp[sell->mpp[j].mul][ii]);
          printf("\n");
          for (int ii=0; ii<basis->nv; ++ii)
            printf("%u ",ht->exp[basis->eh[sell->mpp[j].idx][0]][ii]);
          printf("\n - - -\n");
          for (int ii=0; ii<basis->nv; ++ii)
            printf("%u ",ht->exp[sf_mul][ii]);
          printf("\n");
          for (int ii=0; ii<basis->nv; ++ii)
            printf("%u ",ht->exp[sf->eh[basis->sf[sell->mpp[j].idx].idx[load-1-l]][0]][ii]);
          printf("\n");
          have_sf = 1;
          mul = sf_mul;
          nt  = sf->nt[basis->sf[sell->mpp[j].idx].idx[load-1-l]];
          pol = sf->eh[basis->sf[sell->mpp[j].idx].idx[load-1-l]];
          printf("nt: %u\n", nt);
          break; 
        }
      }
    }
    */
    if (have_sf == 0) {
      mul = sell->mpp[j].mul;
      nt  = basis->nt[sell->mpp[j].idx];
      pol = basis->eh[sell->mpp[j].idx];
    }
    for (k=1; k<nt; ++k)
      enter_monomial_to_preprocessing_hash_list(mul, pol[k], mon);
    /*
    for (k=1; k<basis->nt[sell->mpp[j].idx]; ++k)
      enter_monomial_to_preprocessing_hash_list(sell->mpp[j].mul,
          basis->eh[sell->mpp[j].idx][k], mon);
    */
    // now we distinguish cases for gen1
    if (sp->gen1 == 0) {
      if (ht->idx[sp->lcm] == 0) {
        mon->hpos[mon->load]  = sp->lcm;
        ht->idx[sp->lcm]      = 1;
        mon->load++;
        if (mon->load == mon->size)
          adjust_size_of_preprocessing_hash_list(mon, 2*mon->size);
      }
    } else {
      // first generator for upper part of gbla matrix if there is no other pair
      // with the same lcm
      if (ht->idx[sp->lcm] == 0) {
        mon->hpos[mon->load]  = sp->lcm;
        ht->idx[sp->lcm]      = 2;
#if SYMBOL_DEBUG
      printf("hpos[%u] = %u\n", mon->load, mon->hpos[mon->load]);
      for (int ii=0; ii<ht->nv; ++ii)
        printf("%u ", ht->exp[sp->lcm][ii]);
      printf("\n");
#endif
        mon->nlm++;
        mon->load++;
        if (mon->load == mon->size)
          adjust_size_of_preprocessing_hash_list(mon, 2*mon->size);
        add_spair_generator_to_selection(selu, basis, sp->lcm, sp->gen1);
        j = selu->load-1;
        for (k=1; k<basis->nt[selu->mpp[j].idx]; ++k)
          enter_monomial_to_preprocessing_hash_list(selu->mpp[j].mul,
              basis->eh[selu->mpp[j].idx][k], mon);
      } else {
        add_spair_generator_to_selection(sell, basis, sp->lcm, sp->gen1);
        j = sell->load-1;
        for (k=1; k<basis->nt[sell->mpp[j].idx]; ++k)
          enter_monomial_to_preprocessing_hash_list(sell->mpp[j].mul,
              basis->eh[sell->mpp[j].idx][k], mon);
      }
    }
    // remove the selected pair from the pair set
    free(sp);
  }

  // adjust pair set after removing the bunch of selected pairs
  k = 0;
  for (i=nsel; i<ps->load; ++i) {
    ps->pairs[k] = ps->pairs[i];
    k++;
  }
  ps->load  = k;
}

#endif
