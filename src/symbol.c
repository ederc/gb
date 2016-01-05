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

  // clears hash table index: there we store during symbolic preprocessing a 2
  // if it is a lead monomial and 1 if it is not a lead monomial. all other
  // entries keep 0, thus they are not part of this reduction step
  clear_hash_table_idx(ht);

  // list of monomials that appear in the matrix
  pre_t *mon  = init_preprocessing_hash_list(__GB_SYM_LIST_LEN);
  sel_t *sel  = init_selection(basis->load);
  // list of polynomials and their multipliers
  select_pairs_by_minimal_degree(ps, basis, sel, mon);

  // enter selected spairs without their lead terms first
  enter_spairs_to_preprocessing_hash_list(basis, sel, mon);

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
      // we have reducer, i.e. the monomial is a leading monomial (important for
      // splicing matrix later on
      ht->idx[mon->hpos[idx]] = 2;
      mon->nlm++;
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

void select_pairs_by_minimal_degree(ps_t *ps, gb_t *basis, sel_t *sel, pre_t *mon)
{
  deg_t dmin  = ps->pairs[0]->deg;
  nelts_t i = 0;
  nelts_t nsel;
  spair_t *sp;

  // we assume here that the pair set is already sorted by degree of the lcms
  // (in particular, we assume grevlex ordering)
  while (i < ps->load && ps->pairs[i]->deg == dmin)
    i++;
  nsel  = i;

  sel->deg  = dmin;
  // wVe do not need to check for size problems in sel du to above comment: we
  // have allocated basis->load slots, so enough for each possible element from
  // the basis
#if SYMBOL_DEBUG
  printf("selected pairs in this step of the algorithm:\n");
#endif
  for (i=0; i<nsel; ++i) {
    sp = ps->pairs[i];
#if SYMBOL_DEBUG
    printf("gen1 %u -- gen2 %u\n", sp->gen1, sp->gen2);
#endif
    // first generator
    add_spair_generator_to_selection(basis, sel, sp->lcm, sp->gen1);
    // second generator
    add_spair_generator_to_selection(basis, sel, sp->lcm, sp->gen2);
    sel->nsp++;
    // corresponds to lcm of spair, tracking this information by setting ht->idx
    // to 1 keeping track that this monomial is a lead monomial
    ht->idx[ht->lut[sp->lcm]] = 2;
    mon->hpos[mon->load]      = ht->lut[sp->lcm];;
    mon->load++;
    mon->nlm++;

    // remove the selected pair from the pair set
    free(sp);
  }

  // adjust pair set after removing the bunch of selected pairs
  nelts_t k = 0;
  for (i=nsel; i<ps->load; ++i) {
    ps->pairs[k] = ps->pairs[i];
    k++;
  }
  ps->load  = k;
}

inline void enter_monomial_to_preprocessing_hash_list(const hash_t h1,
    const hash_t h2, pre_t *mon)
{
  hash_t pos = check_in_hash_table_product(h1, h2, ht);
  // only in this case we have this monomial hash for the first time,
  // otherwise it has already been taken care of
  if (ht->idx[pos] == 0) {
    ht->idx[pos]++;
    mon->hpos[mon->load]  = pos;
#if SYMBOL_DEBUG
    for (int i=0; i<ht->nvars; ++i)
      printf("%u ",ht->exp[h1][i]);
    printf("\n");
    for (int i=0; i<ht->nvars; ++i)
      printf("%u ",ht->exp[h2][i]);
    printf("\n");
    for (int i=0; i<ht->nvars; ++i)
      printf("%u ",ht->exp[pos][i]);
    printf("\n");
    printf("new mon %u + %u == %u\n", h1,h2,mon->hpos[mon->load]);
#endif
    mon->load++;
    if (mon->load == mon->size)
      enlarge_preprocessing_hash_list(mon, 2*mon->size);
  }
}

inline void enter_spairs_to_preprocessing_hash_list(const gb_t *basis, sel_t *sel,
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
  mon->nlm    = 0;
  
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

inline int cmp_monomials_by_lead(const void *a, const void *b)
{
  hash_t h1 = *((hash_t *)a);
  hash_t h2 = *((hash_t *)b);

  return (ht->idx[h2] - ht->idx[h1]);
}

inline void sort_columns_by_lead(spd_t *spd)
{
  qsort(spd->col->hpos, spd->col->load, sizeof(hash_t), cmp_monomials_by_lead);
}

inline void sort_lead_columns(spd_t *spd)
{
}

inline void sort_non_lead_columns(spd_t *spd)
{
}
