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
 * \file spairs.c
 * \brief Implementation of handling of pair sets.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "spairs.h"

#define SPAIRS_DEBUG  0

inline ps_t *init_pair_set(gb_t *basis, mp_cf4_ht_t *ht)
{
  nelts_t i;

  ps_t *ps  = (ps_t *)malloc(sizeof(ps_t));
  ps->size  = 2 * basis->size;
  ps->pairs = (spair_t **)malloc(ps->size * sizeof(spair_t *));
  ps->load  = 0;

  // generate spairs with the initial elements in basis
  for (i=1; i<basis->load; ++i) {
    update_pair_set(ps, basis, i);
#if META_DATA_DEBUG
    printf("criteria applied %u / %u\n", meta_data->ncrit_last, meta_data->ncrit_total);
#endif
  }
  // sort pair set by lcms
  sort_pair_set_by_lcm_grevlex(ps);
  return ps;
}

inline void update_pair_set(ps_t *ps, gb_t *basis, nelts_t idx)
{
  nelts_t i;
  // we get maximal idx-1 new pairs
  if (ps->size <= ps->load + (idx-1))
    enlarge_pair_set(ps, 2*ps->size);
  for (i=0; i<idx; ++i) {
    ps->pairs[ps->load+i] = generate_spair(idx, i, basis, ht);
#if SPAIRS_DEBUG
    printf("pair %u %u | %u\n",idx,i,tmp_pairs[i]->deg);
#endif
  }
  // we do not update ps->load at the moment in order to be able to distinguish
  // old and new pairs for the gebauer-moeller update following

  // check product and chain criterion in gebauer moeller style
  // note that we have already marked the pairs for which the product criterion
  // applies in generate_spair()
  gebauer_moeller(ps, basis->eh[idx][0], idx);

  // fix pair set and remove detected pairs
  meta_data->ncrit_last   =   remove_detected_pairs(ps, idx);
  meta_data->ncrit_total  +=  meta_data->ncrit_last;
}

void gebauer_moeller(ps_t *ps, hash_t hash, nelts_t idx)
{
  nelts_t gen1, gen2;
  nelts_t cur_len = ps->load + (idx - 1);
  int i, j; // we need ints to cover cases where i=0 and j=i-1

  // first step: remove elements already in ps due to chain criterion with new
  // pairs in new_pairs
  for (i=0; i<ps->load; ++i) {
    gen1  = ps->pairs[i]->gen1;
    gen2  = ps->pairs[i]->gen2;
    if (ps->pairs[i]->lcm != ps->pairs[ps->load+gen1]->lcm &&
        ps->pairs[i]->lcm != ps->pairs[ps->load+gen2]->lcm &&
        monomial_division(ps->pairs[i]->lcm, hash, ht)) {
      ps->pairs[i]->crit  = CHAIN_CRIT;
    }
  }

  // next: sort new pairs
  qsort(ps->pairs+ps->load, idx-1, sizeof(spair_t **), cmp_spairs_grevlex);
  
  // second step: remove new pairs by themselves w.r.t the chain criterion
  for (i=ps->load; i<cur_len; ++i) {
    if (ps->pairs[i]->crit != NO_CRIT)
      continue;
    for (j=ps->load; j<i; ++j) {
      if (i==j) // smaller lcm eliminated j
        continue;
      if (ps->pairs[j]->lcm == ps->pairs[i]->lcm)
        ps->pairs[j]->crit  = CHAIN_CRIT;
    }
  }
  // third step: remove new pairs via product criterion
  for (i=ps->load; i<cur_len; ++i) {
    if (ps->pairs[i]->crit == PROD_CRIT) {
      // eliminate all new pairs with this lcm
      for (j=ps->load; j<cur_len; ++j) {
        if (ps->pairs[j]->lcm == ps->pairs[i]->lcm) {
          ps->pairs[j]->crit  = CHAIN_CRIT;
        }
      }
    } else { // earlier pairs my eliminate this pair
      if (i > ps->load) {
        for (j=i-1; j>=(int)ps->load; --j) {
          if (ps->pairs[j]->lcm == ps->pairs[i]->lcm) {
            ps->pairs[i]->crit  = CHAIN_CRIT;
            break;
          }
        }
      }
    }
  }
}

inline nelts_t remove_detected_pairs(ps_t *ps, nelts_t idx)
{
  nelts_t cur_len = ps->load + (idx-1);
  nelts_t i, j, nremoved;

  j         = 0;
  nremoved  = 0;
  for (i=0; i<cur_len; ++i) {
    if (ps->pairs[i]->crit != NO_CRIT) {
      nremoved++;
      free(ps->pairs[i]);
      ps->pairs[i]  = NULL;
      continue;
    }
    ps->pairs[j++]  = ps->pairs[i];
  }
  ps->load  = j;

  return nremoved;
}


inline void enlarge_pair_set(ps_t *ps, nelts_t new_size)
{
  ps->pairs = (spair_t **)realloc(ps->pairs, new_size * sizeof(spair_t *));
  ps->size  = new_size;
}

inline void free_pair_set_dynamic_data(ps_t *ps)
{
  free(ps->pairs);
}

inline spair_t *generate_spair(nelts_t gen1, nelts_t gen2, gb_t *basis, mp_cf4_ht_t *ht)
{
  spair_t *sp = (spair_t *)malloc(sizeof(spair_t));
  sp->gen1  = gen1;
  sp->gen2  = gen2;
  sp->lcm   = get_lcm(basis->eh[gen1][0], basis->eh[gen2][0], ht);
  sp->deg   = ht->deg[sp->lcm];
  sp->crit  = NO_CRIT;

  // check for product criterion and mark correspondingly, i.e. we set sp->deg=0
  if (sp->deg == ht->deg[gen1] + ht->deg[gen2])
    sp->crit  = PROD_CRIT;

  return sp;
}

inline int cmp_spairs_grevlex(const void *a, const void *b) {
  spair_t *spa  = *((spair_t **)a);
  spair_t *spb  = *((spair_t **)b);
#if SPAIRS_DEBUG
  printf("%p | %p\n", spa, spb);
  printf("nvars %u\n",ht->nvars);
  printf("%u | %u\n",spa->lcm, spb->lcm);
  printf("%u | %u\n",spa->deg, spb->deg);
#endif
  if (spa->lcm != spb->lcm) {
    // compare degree
    if (spa->deg > spb->deg) {
      return 1;
    } else {
      if (spa->deg != spb->deg)
        return -1;
    }
    // compare reverse lexicographical
    nvars_t i;
    exp_t *expa = ht->exp[spa->lcm];
    exp_t *expb = ht->exp[spb->lcm];
    for (i=ht->nvars-1; i>=0; --i) {
      if (expa[i] < expb[i]) {
        return 1;
      } else {
        if (expa[i] != expb[i])
          return -1;
      }
    }
    return 0;
  } else {
    // we check for spairs labeled by the product criterion, those are moved to
    // the end in order to have an efficient Gebauer-Moeller implementation: If
    // an spair was labeled for product criterion we set the corresponding crit
    // enum to one, otherwise it is zero
    if (spa->deg != spb->deg) {
      return (spa->crit > spb->crit) ? -1 : 1;
    } else {
      // both have the same lcms and are not detected by the product criterion,
      // then we break ties by the first generator
      if (spa->gen1 != spb->gen1) {
        return (spa->gen1 < spb->gen1) ? -1 : 1;
      } else {
        return 0;
      }
    }
  }
}

inline void sort_pair_set_by_lcm_grevlex(ps_t *ps)
{
  qsort(ps->pairs, ps->load, sizeof(spair_t **), cmp_spairs_grevlex);
}
