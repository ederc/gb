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

#define SPARSE_REDUCERS 1

spd_t *symbolic_preprocessing(ps_t *ps, const gb_t *basis, const gb_t *sf)
{
  nelts_t i, idx, last_div, nsel;
  hash_t hash_pos, hash_div = 0;

  // clears hash table index: there we store during symbolic preprocessing a 2
  // if it is a lead monomial and 1 if it is not a lead monomial. all other
  // entries keep 0, thus they are not part of this reduction step
  clear_hash_table_idx(ht);

  nsel  = get_pairs_by_minimal_degree(ps);

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
  nelts_t hi, hio;
  hash_t h, ho, ntmin;
  while (idx < mon->load) {
    hash_pos  = mon->hpos[idx];
    // only if not already a lead monomial, e.g. if coming from spair
    if (ht->idx[hash_pos] != 2) {
      last_div  = ht->div[hash_pos];

      // takes last element in basis and test for division, goes down until we
      // reach the last known divisor last_div
      hio   = 0;
      ho    = 0;
      ntmin = 0;
      i   = last_div == 0 ? basis->st : last_div;
#if SPARSE_REDUCERS
      for (; i<basis->load; ++i) {
        h = monomial_division(hash_pos, basis->eh[i][0], ht);
        if ((h != 0)) {
          hi  = i;
          if (ho != 0) {
            if (basis->nt[hi] < ntmin + ntmin/4) {
              hio   = hi;
              ho    = h;
              ntmin = basis->nt[hi] < ntmin ? basis->nt[hi] : ntmin;
            }
          } else {
            hio   = hi;
            ho    = h;
            ntmin = basis->nt[hio];
          }
        }
      }
#else
      for (; i<basis->load; ++i) {
        h = monomial_division(hash_pos, basis->eh[i][0], ht);
        if ((h != 0)) {
          hi  = i;
          if (ho != 0) {
            //if (basis->nt[hi] < basis->nt[hio] + basis->nt[hio]/2) {
            if (basis->nt[hi] < basis->nt[hio]) {
              hio = hi;
              ho  = h;
            }
          } else {
            hio = hi;
            ho  = h;
          }
        }
      }
#endif
      /*
      i = basis->load-1;
      while (i != last_div) {
        hash_div  = monomial_division(hash_pos, basis->eh[i][0], ht);
        if (hash_div != 0)
          break;
        i--;
      }
      */
      // only if i >= basis->st we have found a reducer.
      // note: all reducers are added to the upper selection list!
      /*
      if (i >= basis->st) {
        if (i == last_div)
          hash_div  = monomial_division(hash_pos, basis->eh[i][0], ht);
        mon->nlm++;
        printf("this is the reducer finally taken %3u\n",i);
        ht->div[hash_pos]  = i;
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
        sel_upp->mpp[sel_upp->load].bi  = i;
        sel_upp->mpp[sel_upp->load].mlm = hash_pos;
        sel_upp->mpp[sel_upp->load].mul = hash_div;
        sel_upp->mpp[sel_upp->load].nt  = basis->nt[i];
        sel_upp->mpp[sel_upp->load].eh  = basis->eh[i];
        sel_upp->mpp[sel_upp->load].cf  = basis->cf[i];
        sel_upp->load++;

        if (basis->sf != NULL)
          try_to_simplify(sel_upp->mpp[sel_upp->load-1], basis, sf);
        printf("this is the simplified reducer finally taken with %3u terms\n", sel_upp->mpp[sel_upp->load-1].nt);

        // now add new monomials to preprocessing hash list
        enter_monomial_to_preprocessing_hash_list(sel_upp->mpp[sel_upp->load-1],
          mon, ht);
      }
      */
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

        if (basis->sf != NULL) {
          try_to_simplify(&sel_upp->mpp[sel_upp->load-1], basis, sf);
          /*
          if (ho != sel_upp->mpp[sel_upp->load-1].mul) {
            printf("simplify changes %3u with %3u terms to  %3u terms || mul %7lu to %7lu\n", hio, basis->nt[hio], sel_upp->mpp[sel_upp->load-1].nt, ho, sel_upp->mpp[sel_upp->load-1].mul);
          }
          */
        }
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

