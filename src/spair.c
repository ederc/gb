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
  ps->pairs = (spair_t **)malloc(ps->size * sizeof(spair_t *));
  ps->load  = 0;

  /* enter input elements as special spairs */
  enter_input_elements_to_pair_set(ps, basis);

  return ps;
}

void enter_input_elements_to_pair_set(ps_t *ps, const gb_t *basis)
{
  nelts_t i;

  for (i=1; i<basis->st; ++i) {
    ps->pairs[ps->load] = generate_input_element_spair(i, basis, ht);
    ps->load++;
  }
}

inline void update_pair_set(ps_t *ps, const gb_t *basis, const nelts_t idx)
{
  nelts_t i;

  /* we get maximal idx-1 new pairs */
  if (ps->size <= ps->load + (idx-1))
    enlarge_pair_set(ps, 2*ps->size);
  /* generate spairs with the initial elements in basis
   * See note on gb_t in src/types.h why we start at position 1 here. */
  for (i=basis->st; i<idx; ++i) {
    ps->pairs[ps->load+i-basis->st] = generate_spair(idx, i, basis, ht);
#if SPAIR_DEBUG
    printf("pair %u, %u + %u | %u\n",idx,i,ps->load,ps->pairs[ps->load+i-basis->st]->deg);
#endif
  }
  /* we do not update ps->load at the moment in order to be able to distinguish
   * old and new pairs for the gebauer-moeller update following */

  /* check product and chain criterion in gebauer moeller style
   * note that we have already marked the pairs for which the product criterion
   * applies in generate_spair() */
  if (idx > basis->st)
    gebauer_moeller(ps, basis, idx);

  /* fix pair set and remove detected pairs */
  meta_data->ncrit_last   =   remove_detected_pairs(ps, idx-basis->st);
  meta_data->ncrit_total  +=  meta_data->ncrit_last;
}

void gebauer_moeller(ps_t *ps, const gb_t *basis, const nelts_t idx)
{
  nelts_t pos1, pos2;
  /* current length can be computed already, need to adjust by the starting
   * position in basis */
  const int cur_len = (int)(ps->load + idx - basis->st);
  /* printf("curlen %u for idx %u | %u\n",cur_len,idx, idx-basis->st); */
  const hash_t hash     = basis->eh[idx][0];
  int i, j; /* we need ints to cover cases where i=0 and j=i-1 */
  const int load  = (int)ps->load;

  /* printf("idx %u | psl %u | cl %u\n", idx, ps->load, cur_len); */
  /* first step: remove elements already in ps due to chain criterion with new
   * pairs in new_pairs */
  for (i=0; i<load; ++i) {
    if ((basis->red[ps->pairs[i]->gen2] != 0) ||
        (basis->red[ps->pairs[i]->gen1] != 0 && basis->red[ps->pairs[i]->gen1] != ps->pairs[i]->gen2)
       )
      ps->pairs[i]->crit  = CHAIN_CRIT;
  }
  for (i=0; i<load; ++i) {
    /* do not check on initial spairs */
    if (ps->pairs[i]->crit == NO_CRIT && ps->pairs[i]->gen1 != ps->pairs[i]->gen2) {
      /* See note on gb_t in src/types.h why we adjust position by -basis->st. */
      pos1  = ps->pairs[i]->gen1 - basis->st;
      pos2  = ps->pairs[i]->gen2 - basis->st;
      if (ps->pairs[i]->lcm != ps->pairs[ps->load+pos1]->lcm &&
          ps->pairs[i]->lcm != ps->pairs[ps->load+pos2]->lcm &&
          check_monomial_division(ps->pairs[i]->lcm, hash, ht) != 0) {
        ps->pairs[i]->crit  = CHAIN_CRIT;
#if SPAIR_DEBUG
        printf("CC for (%u,%u)\n",pos1+1, pos2+1);
#endif
      }
    }
  }

  /* next: sort new pairs */
  qsort(ps->pairs+ps->load, idx-basis->st, sizeof(spair_t **), ht->sort.compare_spairs);

  /* second step: remove new pairs by themselves w.r.t the chain criterion */
  for (i=load; i<cur_len; ++i) {
    if (ps->pairs[i]->crit != NO_CRIT)
      continue;
    for (j=load; j<i; ++j) {
      if (ps->pairs[j]->crit != NO_CRIT) /* smaller lcm eliminated j */
        continue;
      /* if (ps->pairs[i]->lcm == ps->pairs[j]->lcm) { */
      if (check_monomial_division(ps->pairs[i]->lcm, ps->pairs[j]->lcm, ht) != 0 && ps->pairs[i]->lcm != ps->pairs[j]->lcm
) {
      /* if (ps->pairs[i]->lcm != ps->pairs[j]->lcm &&
       *     check_monomial_division(ps->pairs[i]->lcm, ps->pairs[j]->lcm, ht) != 0) { */
        ps->pairs[i]->crit  = CHAIN_CRIT;
#if SPAIR_DEBUG
        printf("2CC for (%u,%u) by (%u,%u)\n",ps->pairs[i]->gen1, ps->pairs[i]->gen2,ps->pairs[j]->gen1, ps->pairs[j]->gen2);
#endif
        break;
      }
    }
  }
#if 0
  for (i=ps->load; i<cur_len; ++i) {
    switch (ps->pairs[i]->crit) {
      case CHAIN_CRIT:
        continue;
        break;
      case PROD_CRIT:
        for (j=i+1; j<cur_len; ++j) {
          if (ps->pairs[j]->lcm == ps->pairs[i]->lcm) {
            if (ps->pairs[j]->crit == NO_CRIT) {
              ps->pairs[j]->crit  = CHAIN_CRIT;
            }
          } else {
            break;
          }
        }
        for (j=i-1; j>ps->load-1; --j) {
          if (ps->pairs[j]->lcm == ps->pairs[i]->lcm) {
            if (ps->pairs[j]->crit == NO_CRIT) {
              ps->pairs[j]->crit  = CHAIN_CRIT;
            }
          } else {
            break;
          }
        }
        break;
      case NO_CRIT:
        j = i-1 > 0 ? i-1 : 0;
        while (j>ps->load-1) {
          if (ps->pairs[j]->lcm == ps->pairs[i]->lcm) {
            ps->pairs[i]->crit  = CHAIN_CRIT;
            j--;
          } else {
            break;
          }
        }
        break;
      default:
        break;
    }
  }
#else
  /* third step: remove new pairs that still have the same lcm */
  hash_t lcm    = ps->pairs[ps->load]->lcm;
  int start = (int)ps->load;
  int end   = 0;
  i         = start+1;
  while (i<cur_len) {
    while (i<cur_len && lcm == ps->pairs[i]->lcm)
      ++i;
    end = i;
    if (end-start>1) {
      /* if we have a PROD_CRIT element we can remove all others, otherwise we
       * have to keep at least the first one: we set start to start+1, if we
       * find a PROD_CRIT element we decrement start again and then we remove
       * NO_CRIT elements with the same lcm */
      start++;
      for (i=start; i<end; ++i) {
        if (ps->pairs[i]->crit == PROD_CRIT) {
          start--;
          break;
        }
      }
      for (i=start; i<end; ++i)
        if (ps->pairs[i]->crit == NO_CRIT)
          ps->pairs[i]->crit  = CHAIN_CRIT;
    }
    if (end == cur_len)
      break;
    lcm   = ps->pairs[end]->lcm;
    start = end;
    i     = start+1;
    continue;
  }
#endif
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
    if (ps->pairs[i]->crit != NO_CRIT) {
#if SPAIR_DEBUG
      printf("REMOVED (%u,%u)\n",ps->pairs[i]->gen1, ps->pairs[i]->gen2);
#endif
      nremoved++;
      /* printf("%p %u\n", ps->pairs[i], i); */
      free(ps->pairs[i]);
      ps->pairs[i]  = NULL;
      continue;
    }
    ps->pairs[j++]  = ps->pairs[i];
  }
  ps->load  = j;

  return nremoved;
}


inline void enlarge_pair_set(ps_t *ps, const nelts_t new_size)
{
  ps->pairs = (spair_t **)realloc(ps->pairs, new_size * sizeof(spair_t *));
  ps->size  = new_size;
}

inline spair_t *generate_input_element_spair(const nelts_t gen2, const gb_t *basis, mp_cf4_ht_t *ht)
{
  spair_t *sp = (spair_t *)malloc(sizeof(spair_t));
  /* sp->gen1  = 0; */
  sp->gen1  = gen2;
  sp->gen2  = gen2;
  sp->lcm   = basis->eh[gen2][0];
  sp->nt    = basis->nt[gen2];
  sp->deg   = ht->deg[sp->lcm];
  sp->crit  = NO_CRIT;

  return sp;
}

inline spair_t *generate_spair(const nelts_t gen1, const nelts_t gen2, const gb_t *basis, mp_cf4_ht_t *ht)
{
  spair_t *sp = (spair_t *)malloc(sizeof(spair_t));
  /* we have to fix the positions where the new basis element is put (gen2),
   * since we are trying to remove as much as possible useless elements in
   * select_pairs(). if we would dynamically adjust the positioning (as done in
   * the below commented out code) we could no longer track this correctly. */
  sp->gen1  = gen2;
  sp->gen2  = gen1;
  /*
  if (basis->nt[gen1] < basis->nt[gen2]) {
    sp->gen1  = gen1;
    sp->gen2  = gen2;
  } else {
    sp->gen1  = gen2;
    sp->gen2  = gen1;
  }
  */
  sp->lcm   = get_lcm(basis->eh[gen1][0], basis->eh[gen2][0], ht);
  sp->nt    = basis->nt[gen1] + basis->nt[gen2];
  sp->deg   = ht->deg[sp->lcm];
  
  /* if one of the generators is redundant we can stop already here and mark it
   * with the CHAIN_CRIT in order to remove it later on */
  if (basis->red[gen2] > 0) {
    /* sp->crit  = CHAIN_CRIT; */
    sp->crit  = PROD_CRIT;
    return sp;
  }
  /* check for product criterion and mark correspondingly, i.e. we set sp->deg=0 */
  if (sp->deg == ht->deg[basis->eh[gen1][0]] + ht->deg[basis->eh[gen2][0]]) {
    sp->crit  = PROD_CRIT;
    return sp;
  }
  /* else */
  sp->crit  = NO_CRIT;
  return sp;
}

inline void add_spair_generator_to_selection(sel_t *sel, const gb_t *basis,
    const hash_t lcm, const nelts_t gen)
{
  hash_t mul;
  mul = get_multiplier(lcm, basis->eh[gen][0], ht);
  if (sel->load == sel->size)
    adjust_size_of_selection(sel, 2*sel->size);
  sel->mpp[sel->load].mlm = lcm;
  sel->mpp[sel->load].mul = mul;
  sel->mpp[sel->load].bi  = gen;
  sel->mpp[sel->load].eh  = basis->eh[gen];
  sel->mpp[sel->load].cf  = basis->cf[gen];
  sel->mpp[sel->load].nt  = basis->nt[gen];
  sel->load++;
}
