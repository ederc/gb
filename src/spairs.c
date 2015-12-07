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
  nelts_t i,j;

  ps_t *ps  = (ps_t *)malloc(sizeof(ps_t));
  ps->size  = 2 * basis->size;
  ps->pairs = (spair_t **)malloc(ps->size * sizeof(spair_t *));
  ps->load  = 0;

  // generate spairs with the initial elements in basis
  for (i=1; i<basis->load; ++i) {
    for (j=0; j<i; ++j) {
      if (ps->size == ps->load)
        enlarge_pair_set(ps, 2*ps->size);
      ps->pairs[ps->load] = generate_spair(i, j, basis, ht);
#if SPAIRS_DEBUG
      printf("pair %u %u | %u\n",i,j,ps->pairs[ps->load]->deg);
#endif
      ps->load++;
    }
  }
  // sort pair set by lcms
  sort_pair_set_by_lcm_grevlex(ps);
  return ps;
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
  spair_t *sp = (spair_t *)malloc(sizeof(spair_t));;
  sp->gen1  = gen1;
  sp->gen2  = gen2;
  sp->lcm   = get_lcm(basis->eh[gen1][0], basis->eh[gen2][0], ht);
  sp->deg   = ht->deg[sp->lcm];
  // check for product criterion and mark correspondingly
  //if (ht->deg[sp.lcm] == ht->deg[gen1] + ht->deg[gen2])
    // TODOa

  return sp;
}

inline int cmp_spairs_grevlex(const void *a, const void *b) {
  spair_t *spa  = *((spair_t **)a);
  spair_t *spb  = *((spair_t **)b);
  //spair_t *spa    = *spap;
  //spair_t *spb    = *spbp;
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
    // an spair was labeled for product criterion we set the corresponding deg
    // of the lcm to zero
    if (spa->deg != spb->deg) {
      return (spa->deg < spb->deg) ? -1 : 1;
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
