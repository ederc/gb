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

spd_t *symbolic_preprocessing(ps_t *ps, const gb_t *basis, const gb_t *sf)
{
  nelts_t i, idx, last_div, nsel;
  hash_t hash_pos;

  // clears hash table index: there we store during symbolic preprocessing a 2
  // if it is a lead monomial and 1 if it is not a lead monomial. all other
  // entries keep 0, thus they are not part of this reduction step
  clear_hash_table_idx(ht);

  nsel  = ht->sort.get_pairs_by_minimal_degree(ps);
  
  meta_data->sel_pairs  = nsel;
  meta_data->curr_deg   = ps->pairs[0]->deg;

  // list of monomials that appear in the matrix
  pre_t *mon      = init_preprocessing_hash_list(2*nsel);
  // the lower part of the gbla matrix resp. the selection is fixed:
  // those are just the second generators of the spairs, thus we need nsel
  // places.
  sel_t *sel_low  = init_selection(nsel);
  sel_t *sel_upp  = init_selection(5*nsel);
  sel_upp->deg    = ps->pairs[0]->deg;
  sel_low->deg    = ps->pairs[0]->deg;
  // list of polynomials and their multipliers
  select_pairs(ps, sel_upp, sel_low, mon, basis, sf, nsel);

  // we use mon as LIFO: last in, first out. thus we can easily remove and add
  // new elements to mon
  idx = 0;
  nelts_t hio;
  hash_t h, ho;
  while (idx < mon->load) {
    hash_pos  = mon->hpos[idx];
    // only if not already a lead monomial, e.g. if coming from spair
    if (ht->idx[hash_pos] != 2) {
      last_div  = ht->div[hash_pos];

      // takes last element in basis and test for division, goes down until we
      // reach the last known divisor last_div
      hio = 0;
      ho  = 0;
      i   = last_div == 0 ? basis->st : last_div;
#if 1
      // max value for an unsigned data type in order to ensure that the first
      // polynomial is taken
      nelts_t nto = -1;
      while (i<basis->load) {
        if (basis->red[i] == 0) {
          if (basis->nt[i] < nto) {
            h = monomial_division(hash_pos, basis->eh[i][0], ht);
            if ((h != 0)) {
              hio = i;
              nto = basis->nt[i];
              ho  = h;
            }
          }
        }
        i++;
      }
#else
      nelts_t b = i;
      // max value for an unsigned data type in order to ensure that the first
      // polynomial is taken
      for (i=basis->load; i>b; --i) {
        h = monomial_division(hash_pos, basis->eh[i-1][0], ht);
        if ((h != 0)) {
          hio = i-1;
          ho  = h;
        }
      }
#endif
      if (hio > 0) {
        mon->nlm++;
        //printf("this is the reducer finally taken %3u\n",hio);
        ht->div[hash_pos]  = hio;
        // if multiple is not already in the selected list
        // we have found another element with such a monomial, since we do not
        // take care of the lead monomial below when entering the other lower
        // order monomials, we have to adjust the idx for this given monomial
        // here.
        // we have reducer, i.e. the monomial is a leading monomial (important for
        // splicing matrix later on
        ht->idx[hash_pos] = 2;
        if (sel_upp->load == sel_upp->size)
          adjust_size_of_selection(sel_upp, 2*sel_upp->size);
        sel_upp->mpp[sel_upp->load].bi  = hio;
        sel_upp->mpp[sel_upp->load].mlm = hash_pos;
        sel_upp->mpp[sel_upp->load].mul = ho;
        sel_upp->mpp[sel_upp->load].nt  = basis->nt[hio];
        sel_upp->mpp[sel_upp->load].eh  = basis->eh[hio];
        sel_upp->mpp[sel_upp->load].cf  = basis->cf[hio];
        sel_upp->load++;

        if (basis->sf != NULL)
          try_to_simplify(&sel_upp->mpp[sel_upp->load-1], basis, sf);
        // now add new monomials to preprocessing hash list
        enter_monomial_to_preprocessing_hash_list(sel_upp->mpp[sel_upp->load-1],
            mon, ht);
      }
    }
    idx++;
  }

  // next we store the information needed to construct the GBLA matrix in the
  // following
  spd_t *mat  = (spd_t *)malloc(sizeof(spd_t));

  // adjust memory
  adjust_size_of_selection(sel_upp, sel_upp->load);
  adjust_size_of_preprocessing_hash_list(mon, mon->load);
#if SYMBOL_DEBUG
  for (int ii=0; ii<sel_low->load; ++ii) {
    for (int jj=0; jj<ht->nv; ++jj) {
      printf("%u ",ht->exp[sel_low->mpp[ii].mul][jj]);
    }
    printf(" || ");
    for (int jj=0; jj<ht->nv; ++jj) {
      printf("%u ",ht->exp[basis->eh[sel_low->mpp[ii].bi][0]][jj]);
    }
    printf(" ||| ");
    for (int jj=0; jj<ht->nv; ++jj) {
      printf("%u ",ht->exp[sel_low->mpp[ii].mul][jj] + ht->exp[basis->eh[sel_low->mpp[ii].bi][0]][jj]);
    }
    printf("\n");
  }
#endif
  mat->selu = sel_upp;
  mat->sell = sel_low;
  mat->col  = mon;

  return mat;
}
