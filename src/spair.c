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
inline ps_t *initialize_pair_set(const gb_t *basis, mp_cf4_ht_t *ht)
{
  ps_t *ps  = (ps_t *)malloc(sizeof(ps_t));
  ps->size  = 2 * basis->size;
  ps->pairs = (spair_t **)malloc(ps->size * sizeof(spair_t *));
  ps->load  = 0;

  // enter input elements as special spairs
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

  nelts_t ctr = 0;
  // we get maximal idx-1 new pairs
  if (ps->size <= ps->load + (idx-1))
    enlarge_pair_set(ps, 2*ps->size);
  // generate spairs with the initial elements in basis
  // See note on gb_t in src/types.h why we start at position 1 here.
  for (i=basis->st; i<idx; ++i) {
    ps->pairs[ps->load+i-basis->st] = generate_spair(idx, i, basis, ht);
    ctr++;
#if SPAIR_DEBUG
    printf("pair %u, %u + %u | %u\n",idx,i,ps->load,ps->pairs[ps->load+i-basis->st]->deg);
#endif
  }
  // we do not update ps->load at the moment in order to be able to distinguish
  // old and new pairs for the gebauer-moeller update following

  // check product and chain criterion in gebauer moeller style
  // note that we have already marked the pairs for which the product criterion
  // applies in generate_spair()
  gebauer_moeller(ps, basis, idx, ctr);

  // fix pair set and remove detected pairs
  meta_data->ncrit_last   =   remove_detected_pairs(ps, basis, ctr);
  meta_data->ncrit_total  +=  meta_data->ncrit_last;
}

void gebauer_moeller(ps_t *ps, const gb_t *basis, const nelts_t idx, const nelts_t ctr)
{
  nelts_t pos1, pos2;
  // current length can be computed already, need to adjust by the starting
  // position in basis
  const nelts_t cur_len = ps->load + ctr;
  const hash_t hash     = basis->eh[idx][0];
  int i, j; // we need ints to cover cases where i=0 and j=i-1

  //printf("idx %u | psl %u | cl %u\n", idx, ps->load, cur_len);
  // first step: remove elements already in ps due to chain criterion with new
  // pairs in new_pairs
  for (i=0; i<ps->load; ++i) {
    // do not check on initial spairs
    if (ps->pairs[i]->gen1 != 0) {
      // See note on gb_t in src/types.h why we adjust position by -basis->st.
      pos1  = ps->pairs[i]->gen1 - basis->st;
      pos2  = ps->pairs[i]->gen2 - basis->st;
      //printf("gen1 %u | gen2 %u | pos1 %u | pos2 %u\n", ps->pairs[i]->gen1, ps->pairs[i]->gen2, pos1, pos2);
      if (ps->pairs[i]->lcm != ps->pairs[ps->load+pos1]->lcm &&
          ps->pairs[i]->lcm != ps->pairs[ps->load+pos2]->lcm &&
          monomial_division(ps->pairs[i]->lcm, hash, ht)) {
        ps->pairs[i]->crit  = CHAIN_CRIT;
#if SPAIR_DEBUG
        printf("CC for (%u,%u)\n",pos1+1, pos2+1);
#endif
      }
    }
  }

  // next: sort new pairs
  qsort(ps->pairs+ps->load, ctr, sizeof(spair_t **), cmp_spairs_grevlex);
  
  // second step: remove new pairs by themselves w.r.t the chain criterion
  for (i=ps->load; i<cur_len; ++i) {
    if (ps->pairs[i]->crit != NO_CRIT)
      continue;
    for (j=ps->load; j<i; ++j) {
      if (i==j || ps->pairs[j] != NO_CRIT) // smaller lcm eliminated j
        continue;
      if (ps->pairs[j]->lcm == ps->pairs[i]->lcm) {
        ps->pairs[j]->crit  = CHAIN_CRIT;
#if SPAIR_DEBUG
        printf("2CC for (%u,%u)\n",ps->pairs[j]->gen1, ps->pairs[j]->gen2);
#endif
      }
    }
  }
  // third step: remove new pairs via product criterion
  for (i=ps->load; i<cur_len; ++i) {
    if (ps->pairs[i]->crit == PROD_CRIT) {
      // eliminate all new pairs with this lcm
      for (j=ps->load; j<cur_len; ++j) {
        if (ps->pairs[j]->lcm == ps->pairs[i]->lcm) {
          ps->pairs[j]->crit  = CHAIN_CRIT;
#if SPAIR_DEBUG
          printf("3CC for (%u,%u)\n",ps->pairs[j]->gen1, ps->pairs[j]->gen2);
#endif
        }
      }
    } else { // earlier pairs may eliminate this pair
      if (i > ps->load) {
        for (j=i-1; j>=(int)ps->load; --j) {
          if (ps->pairs[j]->lcm == ps->pairs[i]->lcm) {
            ps->pairs[i]->crit  = CHAIN_CRIT;
#if SPAIR_DEBUG
            printf("4CC for (%u,%u)\n",ps->pairs[i]->gen1, ps->pairs[i]->gen2);
#endif
            break;
          }
        }
      }
    }
  }
}

inline nelts_t remove_detected_pairs(ps_t *ps, const gb_t *basis, const nelts_t ctr)
{
  // current length can be computed already, need to adjust by the starting
  // position in basis
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
      //printf("%p %u\n", ps->pairs[i], i);
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

inline void free_pair_set(ps_t *ps)
{
  free(ps->pairs);
  free(ps);
  ps  = NULL;
}

inline spair_t *generate_input_element_spair(const nelts_t gen2, const gb_t *basis, mp_cf4_ht_t *ht)
{
  spair_t *sp = (spair_t *)malloc(sizeof(spair_t));
  sp->gen1  = 0;
  sp->gen2  = gen2;
  sp->lcm   = basis->eh[gen2][0];
  sp->deg   = ht->deg[sp->lcm];
  sp->crit  = NO_CRIT;

  return sp;
}

inline spair_t *generate_spair(const nelts_t gen1, const nelts_t gen2, const gb_t *basis, mp_cf4_ht_t *ht)
{
  spair_t *sp = (spair_t *)malloc(sizeof(spair_t));
  // number of terms in polynomials decides which one is going to the upper part
  // of the gbla matrix (sp->gen1) and which one to the lower part (sp->gen2)
  if (basis->nt[gen1] < basis->nt[gen2]) {
    sp->gen1  = gen1;
    sp->gen2  = gen2;
  } else {
    sp->gen1  = gen2;
    sp->gen2  = gen1;
  }
  sp->lcm   = get_lcm(basis->eh[gen1][0], basis->eh[gen2][0], ht);
  sp->deg   = ht->deg[sp->lcm];
  
  // if one of the generators is redundant we can stop already here and mark it
  // with the CHAIN_CRIT in order to remove it later on
  /*
  if (basis->red[gen1] | basis->red[gen2]) {
    sp->crit  = CHAIN_CRIT;
    return sp;
  }
  */
  sp->crit  = NO_CRIT;
  // check for product criterion and mark correspondingly, i.e. we set sp->deg=0
  if (sp->deg == ht->deg[basis->eh[gen1][0]] + ht->deg[basis->eh[gen2][0]])
    sp->crit  = PROD_CRIT;

  return sp;
}

inline int cmp_spairs_grevlex(const void *a, const void *b)
{
  spair_t *spa  = *((spair_t **)a);
  spair_t *spb  = *((spair_t **)b);
#if SPAIR_DEBUG
  printf("%p | %p\n", spa, spb);
  printf("nvars %u\n",ht->nv);
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
    // Note that this loop only works since we assume that spa->lcm =/=
    // spb->lcm. Otherwise i might get to zero and decremented again, which
    // means that we would get into an infinite loop as nvars_t is unsigned.
    for (i=ht->nv-1; i>=0; --i) {
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

nelts_t get_pairs_by_minimal_degree(ps_t *ps)
{
  // sort pair set by lcms
  sort_pair_set_by_lcm_grevlex(ps);

  deg_t dmin  = ps->pairs[0]->deg;
  nelts_t i = 0;

  // we assume here that the pair set is already sorted by degree of the lcms
  // (in particular, we assume grevlex ordering)
  while (i < ps->load && ps->pairs[i]->deg == dmin)
    i++;

  return i;
}


inline void sort_pair_set_by_lcm_grevlex(ps_t *ps)
{
  qsort(ps->pairs, ps->load, sizeof(spair_t **), cmp_spairs_grevlex);
}

inline void add_spair_generator_to_selection(sel_t *sel, const gb_t *basis,
    const hash_t lcm, const nelts_t gen)
{
  hash_t mul;
  mul = get_multiplier(lcm, basis->eh[gen][0], ht);
  if (sel->load == sel->size)
    adjust_size_of_selection(sel, 2*sel->size);
  sel->mpp[sel->load].mlm  = lcm;
  sel->mpp[sel->load].mul  = mul;
  sel->mpp[sel->load].idx  = gen;
  sel->load++;
}

inline sel_t *init_selection(const nelts_t size)
{
  sel_t *sel  = (sel_t *)malloc(sizeof(sel_t));
  sel->size   = size;
  sel->load   = 0;

  sel->mpp    = (mpp_t *)malloc(size * sizeof(mpp_t));

  return sel;
}

inline void adjust_size_of_selection(sel_t *sel, const nelts_t new_size)
{
  sel->size = new_size;
  sel->mpp  = realloc(sel->mpp, sel->size * sizeof(mpp_t));
}

inline void free_selection(sel_t *sel)
{
  free(sel->mpp);
  free(sel);
  sel = NULL;
}
