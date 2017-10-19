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
  nelts_t i, idx, nsel;
  hash_t hash_pos;

  nsel  = ht->sort.get_pairs_by_minimal_degree(ps);

  /* allocate memory for matrices */
  AB->row = realloc(AB->row, 20 * nsel * sizeof(src_t *));
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
  
#if NOT_ADDING_MUL_TO_HT
  hash_t mul_hash;
  deg_t mul_deg;
  exp_t *mul_exp  = (exp_t *)malloc(ht->nv * sizeof(exp_t));
#endif
  /* hash_t h, ho; */
  while (idx < mon->load) {
    hash_pos  = mon->hash[idx];

    if (ht->idx[hash_pos] == 2) {
      idx++;
      continue;
    }
    i = ht->div[hash_pos] > 0 ? ht->div[hash_pos] : basis->st;

    /* printf("i %u last divisor checked with to start sc\n", i); */
    while (i<basis->load) {
      /* printf("i %u\n", i); */
      if (check_monomial_division(hash_pos, basis->p[i][2], ht)) {
        while (basis->red[i] != 0)
          i = basis->red[i];
        ht->div[hash_pos] = i;
        goto done;
      }
      i++;
    }
    ht->div[hash_pos] = i;
    idx++;
    continue;

done:
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
    AB->row[AB->rk] = (src_t *)malloc(basis->p[i][1] * sizeof(src_t));
    memcpy(AB->row[AB->rk], basis->p[i],
        basis->p[i][1] * sizeof(src_t));
#if NOT_ADDING_MUL_TO_HT
    /* printf("hash values: %d - %d = %d\n", ht->val[hash_pos], ht->val[basis->p[i][2]], mul_hash); */
    mul_hash = ht->val[hash_pos] - ht->val[basis->p[i][2]];
    for (size_t j = 0; j < ht->nv; ++j) {
      mul_exp[j]  = ht->exp[hash_pos][j] - ht->exp[basis->p[i][2]][j];
    }
    /* printf("--mul--\n");
     * for (size_t j = 0; j < ht->nv; ++j) {
     *   printf("%u", mul_exp[j]);
     * }
     * printf("\n");
     * printf("--div--\n");
     * for (size_t j = 0; j < ht->nv; ++j) {
     *   printf("%u", ht->exp[basis->p[i][2]][j]);
     * }
     * printf("\n");
     * printf("--elt--\n");
     * for (size_t j = 0; j < ht->nv; ++j) {
     *   printf("%u", ht->exp[hash_pos][j]);
     * }
     * printf("\n"); */
    mul_deg = 0;
    for (size_t j = 0; j < ht->nv; ++j) {
      mul_deg +=  mul_exp[j];
    }
#else
    hash_t mul = get_multiplier(hash_pos, basis->p[i][2], ht);
#endif

    /* hash_t mul = get_multiplier_after_poly(
     *     hash_pos, basis->p[i][2], basis->p[i][1], ht); */
    /* printf("mul %u\n", mul); */
    /* printf("NEW POLY STARTS %u\n", i); */
    for (size_t j = 2; j < AB->row[AB->rk][1]; j = j+2) {
#if NOT_ADDING_MUL_TO_HT
      AB->row[AB->rk][j]  = check_in_hash_table_product_special(
          AB->row[AB->rk][j], mul_hash, mul_deg, mul_exp, ht);
      /* printf("b %d\n", AB->row[AB->rk][j]); */
#else
      AB->row[AB->rk][j]  = check_in_hash_table_product(
          mul, AB->row[AB->rk][j], ht);
#endif
      if (ht->idx[AB->row[AB->rk][j]] == 0) {
        ht->idx[AB->row[AB->rk][j]] = 1;
        add_to_monomial_list(mon, AB->row[AB->rk][j]);
        /* printf("b mon->hash[%d] %d\n", mon->load, mon->hash[mon->load]); */
      }
    }
    AB->rk++;
    ht->idx[hash_pos] = 2;
    mon->nlm++;
    idx++;
  }
#if NOT_ADDING_MUL_TO_HT
  free(mul_exp);
#endif
}
