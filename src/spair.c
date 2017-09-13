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
 * \file spair.c
 * \brief Implementation of handling of pair sets.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "spair.h"
inline ps_t *initialize_pair_set(const gb_t *basis)
{
  ps_t *ps  = (ps_t *)malloc(sizeof(ps_t));
  ps->size  = 2 * basis->size;
  ps->pairs = (spair_t *)malloc(ps->size * sizeof(spair_t));
  ps->load  = 0;

  /* enter input elements as special spairs */
  enter_input_elements_to_pair_set(ps, basis);

  return ps;
}

void enter_input_elements_to_pair_set(ps_t *ps, const gb_t *basis)
{
  nelts_t i;

  for (i=1; i<basis->st; ++i) {
    generate_input_element_spair(ps, i, basis, ht);
  }
}


void gebauer_moeller(ps_t *ps, const gb_t *basis, const nelts_t idx)
{
  nelts_t pos1, pos2;
  /* current length can be computed already, need to adjust by the starting
   * position in basis */
  int cur_len = (int)(ps->load + idx - basis->st);
  /* printf("curlen %u for idx %u | %u\n",cur_len,idx, idx-basis->st); */
  const hash_t hash     = basis->eh[idx][0];
  int i, j; /* we need ints to cover cases where i=0 and j=i-1 */
  const int load  = (int)ps->load;

  spair_t *sp = ps->pairs;

  for (i=0; i<load; ++i) {
    /* do not check on initial spairs */
    /* See note on gb_t in src/types.h why we adjust position by -basis->st. */
    pos1  = sp[i].gen1 - basis->st;
    pos2  = sp[i].gen2 - basis->st;
    if (sp[i].lcm != sp[ps->load+pos1].lcm &&
        sp[i].lcm != sp[ps->load+pos2].lcm &&
        check_monomial_division(sp[i].lcm, hash, ht) != 0) {
      sp[i].crit  = CHAIN_CRIT;
    }
  }

  /* next: sort new pairs */
  qsort(sp+ps->load, idx-basis->st, sizeof(spair_t), ht->sort.compare_spairs);


  /* second step: remove new pairs by themselves w.r.t the chain criterion */

  int next;
  for (j = load; j < cur_len; ++j) {
    next  = j;
    if (sp[j].crit != CHAIN_CRIT) {
      i = j+1;
      while (i < cur_len && sp[i].lcm == sp[j].lcm)
        ++i;
      next  = i-1;
      while (i < cur_len) {
        /* if (check_monomial_division(sp[i].lcm, sp[j].lcm, ht) != 0) { */
        if (sp[i].crit == NO_CRIT &&
            check_monomial_division(sp[i].lcm, sp[j].lcm, ht) != 0) {
          sp[i].crit  = CHAIN_CRIT;
        }
        ++i;
      }
    }
    j = next;
  }

/*   nelts_t new_len = remove_chain_crit_pairs(ps, cur_len);
 *   [> fix pair set and remove detected pairs <]
 *   meta_data->ncrit_last   =   cur_len - (int)new_len;
 *   meta_data->ncrit_total  +=  meta_data->ncrit_last;
 *   [> printf("old cur_len %d | new cur_len %u\n", cur_len, new_len); <]
 *   cur_len = (int)new_len;
 *  */

  /* third step */
  for (i = load; i < cur_len; ++i) {
    if (sp[i].crit == CHAIN_CRIT) {
      continue;
    }
    const hash_t lcm  = sp[i].lcm;
    if (sp[i].crit == NO_CRIT) {
      j = i+1;
      for (j = i+1; j < cur_len; ++j) {
        if (sp[j].lcm == lcm) {
          sp[j].crit = CHAIN_CRIT;
        } else {
          break;
        }
      }
    }
  }

  /* fix pair set and remove detected pairs */
  meta_data->ncrit_last   =   remove_detected_pairs(ps, cur_len-ps->load);
  meta_data->ncrit_total  +=  meta_data->ncrit_last;
}

inline nelts_t remove_chain_crit_pairs(ps_t *ps, const nelts_t cur_len)
{
  /* current length can be computed already, need to adjust by the starting
   * position in basis */
  nelts_t i, j;

  j = ps->load;
  for (i = ps->load; i < cur_len; ++i) {
    if (ps->pairs[i].crit == CHAIN_CRIT) {
      continue;
    }
    ps->pairs[j++]  = ps->pairs[i];
  }
  return j;
}
inline nelts_t remove_detected_pairs(ps_t *ps, const nelts_t ctr)
{
  /* current length can be computed already, need to adjust by the starting
   * position in basis */
  const nelts_t cur_len = ps->load + ctr;
  nelts_t i, j, nremoved;

  j         = 0;
  nremoved  = 0;
  for (i=0; i<cur_len; ++i) {
    if (ps->pairs[i].crit != NO_CRIT) {
#if SPAIR_DEBUG
      printf("REMOVED %u (%u,%u) -> %u\n",i, ps->pairs[i].gen1, ps->pairs[i].gen2, ps->pairs[i].crit);
#endif
      nremoved++;
      /* printf("%p %u\n", ps->pairs[i], i); */
      /* free(ps->pairs[i]);
       * ps->pairs[i]  = NULL; */
      continue;
    }
    ps->pairs[j++]  = ps->pairs[i];
  }
  ps->load  = j;

  return nremoved;
}


inline void enlarge_pair_set(ps_t *ps, const nelts_t new_size)
{
  ps->pairs = (spair_t *)realloc(ps->pairs, new_size * sizeof(spair_t));
  ps->size  = new_size;
}

inline void generate_input_element_spair(ps_t *ps, const nelts_t gen2, const gb_t *basis, ht_t *ht)
{
  /* spair_t *sp = (spair_t *)malloc(sizeof(spair_t)); */
  /* sp->gen1  = 0; */
  spair_t *sp = ps->pairs + ps->load;
  
  sp->gen1    = gen2;
  sp->gen2    = gen2;
  sp->lcm     = basis->eh[gen2][0];
  sp->deg     = ht->deg[sp->lcm];
  sp->crit    = NO_CRIT;

  ps->load++;
}


inline void generate_spair(ps_t *ps, const nelts_t gen1,
    const nelts_t gen2, const gb_t *basis, ht_t *ht)
{
  spair_t *sp = ps->pairs + ps->load + gen2 - basis->st;
  /* we have to fix the positions where the new basis element is put (gen2),
   * since we are trying to remove as much as possible useless elements in
   * select_pairs(). if we would dynamically adjust the positioning (as done in
   * the below commented out code) we could no longer track this correctly. */
  sp->gen1  = gen2;
  sp->gen2  = gen1;

  /* if (basis->nt[gen1] < basis->nt[gen2]) {
   *   sp->gen1  = gen1;
   *   sp->gen2  = gen2;
   * } else {
   *   sp->gen1  = gen2;
   *   sp->gen2  = gen1;
   * } */

  sp->lcm   = get_lcm(basis->eh[gen1][0], basis->eh[gen2][0], ht);

  sp->deg   = ht->deg[sp->lcm];

  /* if one of the generators is redundant we can stop already here and mark it
   * with the CHAIN_CRIT in order to remove it later on */
  /* else */
  if (basis->red[gen2] > 0) {
    /* printf("%u - %u || %u\n", gen2, gen1, ps->load + gen2 - basis->st); */
    sp->crit  = CHAIN_CRIT;
    return;
  }
  /* check for product criterion and mark correspondingly, i.e. we set sp->deg=0 */
  if (sp->deg == ht->deg[basis->eh[gen1][0]] + ht->deg[basis->eh[gen2][0]]) {
    sp->crit  = PROD_CRIT;
    return;
  }
  sp->crit  = NO_CRIT;
  return;
}

inline void add_spair_generator_to_selection(sel_t *sel, const gb_t *basis,
    const hash_t lcm, const nelts_t gen)
{
  hash_t mul;
  mul = get_multiplier(lcm, basis->eh[gen][0], ht);
  /* if (sel->load == sel->size)
   *   adjust_size_of_selection(sel, 2*sel->size); */
  sel->mpp[sel->load].mlm = lcm;
  sel->mpp[sel->load].mul = mul;
  sel->mpp[sel->load].bi  = gen;
  sel->mpp[sel->load].sf  = 0;
#if 0
  sel->mpp[sel->load].eh  = basis->eh[gen];
  sel->mpp[sel->load].cf  = basis->cf[gen];
  sel->mpp[sel->load].nt  = basis->nt[gen];
#endif
  sel->load++;
}
