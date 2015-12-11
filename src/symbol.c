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
 * \file symbol.c
 * \brief Implementation of the symbolic pre- and postprocessing in gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "symbol.h"

spd_t *symbolic_preprocessing(ps_t *ps, gb_t *basis)
{
  nelts_t i, j, k, idx, last_div;
  hash_t hash_pos, hash_div = 0;

  // list of polynomials and their multipliers
  sel_t *sel  = select_pairs_by_minimal_degree(ps, basis);
  // list of monomials that appear in the matrix
  pre_t *mon  = init_preprocessing_hash_list(__GB_SYM_LIST_LEN);

  // enter selected spairs without their lead terms first
  enter_spairs_to_preprocessing_hash_list(sel, basis, mon);

  // we use mon as LIFO: last in, first out. thus we can easily remove and add
  // new elements to mon
  idx = 0;
  while (idx < mon->load) {
    hash_pos  = mon->hpos[idx];
    last_div  = ht->div[hash_pos];

    // takes last element in basis and test for division, goes down until we
    // reach the last known divisor last_div
    i = basis->load-1;
    while (i != last_div) {
      hash_div  = monomial_division(hash_pos, basis->eh[i][0], ht);
      if (hash_div != 0)
        break;
      i--;
    }
    // only if i > 0 we have found a reducer
    if (i != 0) {
      ht->div[hash_pos]  = i;
      for (j=0; j<sel->mload[i]; ++j)
        if (hash_div  ==  sel->mul[i][j])
          break;
      // if multiple is not already in the selected list
      // we have found another element with such a monomial, since we do not
      // take care of the lead monomial below when entering the other lower
      // order monomials, we have to adjust the idx for this given monomial
      // here.
      ht->idx[hash_pos]++;
      if (j == sel->mload[i]) {
        // check for enlarging
        check_enlargement_mul_in_selection(sel, 2*sel->msize[i], i);
        sel->mul[i][sel->mload[i]]  = hash_div;
        sel->mload[i]++;
        sel->load++;

        // now add new monomials to preprocessing hash list
        for (k=1; k<basis->nt[i]; ++k)
          enter_monomial_to_preprocessing_hash_list(sel->mul[i][sel->mload[i]-1],
              basis->eh[i][k], mon);

      }
    }
    idx++;
  }

  // next we store the information needed to construct the GBLA matrix in the
  // following
  spd_t *mat  = (spd_t *)malloc(sizeof(spd_t));
  mat->sel  = sel;
  mat->col  = mon;

  return mat;
}

inline void enter_monomial_to_preprocessing_hash_list(const hash_t h1,
    const hash_t h2, pre_t *mon)
{
  hash_t pos = check_in_hash_table_product(h1, h2, ht);
  ht->idx[pos]++;
  // only in this case we have this monomial hash for the first time,
  // otherwise it has already been taken care of
  if (ht->idx[pos] == 1) {
    mon->hpos[mon->load]  = pos;
    mon->load++;
    if (mon->load == mon->size)
      enlarge_preprocessing_hash_list(mon, 2*mon->size);
  }
}

inline void enter_spairs_to_preprocessing_hash_list(sel_t *sel, const gb_t *basis,
    pre_t *mon)
{
  nelts_t i, j, k;

  // add lower order terms of spairs that are in sel already
  for (i=1; i<basis->load; ++i) {
    for (j=0; j<sel->mload[i]; ++j) {
      // we have already taken care of the lead terms, thus we start with k=1
      for (k=1; k<basis->nt[i]; ++k) {
        enter_monomial_to_preprocessing_hash_list(sel->mul[i][j], basis->eh[i][k], mon);
      }
    }
  }
}

inline pre_t *init_preprocessing_hash_list(const nelts_t size)
{
  // allocate a list for hashes of monomials to be checked in the symbolic
  // preprocessing
  pre_t *mon  = (pre_t *)malloc(sizeof(pre_t));
  mon->hpos   = (hash_t *)malloc(__GB_SYM_LIST_LEN * sizeof(hash_t));
  mon->size   = __GB_SYM_LIST_LEN;
  mon->load   = 0;

  return mon;
}

inline void enlarge_preprocessing_hash_list(pre_t *hl, const nelts_t size)
{
  hl->hpos  = realloc(hl->hpos, size * sizeof(hash_t));
  hl->size  = size;
}

inline void free_preprocessing_hash_list(pre_t *hl)
{
  free(hl->hpos);
  free(hl);
  hl  = NULL;
}

inline void free_symbolic_preprocessing_data(spd_t *spd)
{
  free_preprocessing_hash_list(spd->col);
  free_selection(spd->sel);
  free(spd);
  spd = NULL;
}

