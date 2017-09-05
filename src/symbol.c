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

void symbolic_preprocessing(ps_t *ps, smc_t *AB, smc_t *CD,
    pre_t *mon, const gb_t *basis)
{
  nelts_t i, idx, last_div, nsel;
  hash_t hash_pos;

  nsel  = ht->sort.get_pairs_by_minimal_degree(ps);
  
  /* allocate memory for matrices */
  AB->row = realloc(AB->row, 2 * nsel * sizeof(src_t *));
  AB->ncl = AB->ncr = AB->rk  = 0;
  CD->row = realloc(CD->row, 2 * nsel * sizeof(src_t *));
  CD->ncl = CD->ncr = CD->rk  = 0;
  AB->nr  = CD->nr  = 2 * nsel;

  /* check if limit for spair handling is set */
  if (basis->max_sel != 0) {
    while (nsel > basis->max_sel)
      nsel  = nsel / 2;
    /* nsel  = nsel < basis->max_sel ? nsel : nsebasis->max_sel; */
    /* try to keep all spairs of same lcm together */
    while (nsel != ps->load && ps->pairs[nsel-1].lcm == ps->pairs[nsel].lcm)
      nsel++;
  }
  
  meta_data->sel_pairs  = nsel;
  meta_data->curr_deg   = ps->pairs[0].deg;

  /* list of polynomials and their multipliers */
  select_pairs(ps, AB, CD, mon, basis, nsel);

  /* fix lower matrix size, won't change anymore */
  CD->nr  = CD->rk;
  CD->row = realloc(CD->row, CD->nr * sizeof(src_t *));

  /* we use mon as LIFO: last in, first out. thus we can easily remove and add
   * new elements to mon */
  idx = 0;
  i = 0;
  nelts_t hio;
  /* hash_t h, ho; */
  while (idx < mon->load) {
    hash_pos  = mon->hash[idx];
    /* only if not already a lead monomial, e.g. if coming from spair */
    if (ht->idx[hash_pos] != 2) {
      last_div  = ht->div[hash_pos];

      /* takes last element in basis and test for division, goes down until we
       * reach the last known divisor last_div */
      hio = 0;
      /* ho  = 0; */
      if (last_div > 0 && basis->red[last_div] == 0) {
        hio = last_div;
        goto done;
      }
      i   = last_div == 0 ? basis->st : last_div+1;

      while (i<basis->load) {
        if (basis->red[i] == 0 && check_monomial_division(hash_pos, basis->p[i][2], ht)) {
          hio = i;
          break;
        }
        i++;
      }
      if (hio > 0) {
        done:
        ht->div[hash_pos]  = hio;
        /* if multiple is not already in the selected list
         * we have found another element with such a monomial, since we do not
         * take care of the lead monomial below when entering the other lower
         * order monomials, we have to adjust the idx for this given monomial
         * here.
         * we have reducer, i.e. the monomial is a leading monomial (important for
         * splicing matrix later on */
        if (AB->rk == AB->nr) {
          AB->row =   realloc(AB->row, 2 * AB->rk * sizeof(src_t *));
          AB->nr  *=  2;
        }
        AB->row[AB->rk] = (src_t *)malloc(basis->p[hio][1] * sizeof(src_t));
        memcpy(AB->row[AB->rk], basis->p[hio],
            basis->p[hio][1] * sizeof(src_t));
        hash_t mul = get_multiplier(hash_pos, basis->p[hio][2], ht);
        /* printf("mul %u\n", mul); */
        for (size_t i = 2; i < AB->row[AB->rk][1]; i = i+2) {
          AB->row[AB->rk][i]  = check_in_hash_table_product(
              mul, AB->row[AB->rk][i], ht);
          if (ht->idx[AB->row[AB->rk][i]] == 0) {
            ht->idx[AB->row[AB->rk][i]] = 1;
            add_to_monomial_list(mon, AB->row[AB->rk][i]);
          }
        }
        AB->rk++;
        ht->idx[hash_pos] = 2;
        mon->nlm++;
      }
    }
    idx++;
  }
}
