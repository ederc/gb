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
 * \file spair.h
 * \brief Implementation of handling of pair sets.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_SPAIR_H
#define GB_SPAIR_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include <config.h>
#include "types.h"
#include "hash.h"
#include <cli/io.h>

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#ifndef SPAIR_DEBUG
#define SPAIR_DEBUG  0
#endif

/**
 * \brief Initialize pair set
 *
 * \param input elements input
 *
 * \param hash table ht
 */
ps_t *initialize_pair_set(const gb_t *input);

/**
 * \brief Enters the input elements as sepcial spairs to the pair set:
 * Their second generator is zero. Thus, when handling these pairs we
 * automatically interreduce the input.
 *
 * \note We do not need to check for pair set enlargement since the number of
 * added pairs here is fixed and smaller than the initially allocated size.
 *
 * \param pair set ps
 *
 * \param input elements input
 */
void enter_input_elements_to_pair_set(ps_t *ps, const gb_t *input);

/**
 * \brief Enlarge pair set ps to size new_size
 *
 * \param pair set ps
 *
 * \param new size new_size
 */
void enlarge_pair_set(ps_t *ps, const nelts_t new_size);

/**
 * \brief Frees dynamically allocated memory from pair set
 *
 * \param pair set ps
 */
static inline void free_pair_set(ps_t **ps_in)
{
  ps_t *ps  = *ps_in;
  free(ps->pairs);
  free(ps);
  ps      = NULL;
  *ps_in  = ps;
}


/* updates many pairs at once starting from the first index fidx
 * up to the current basis->load */
void update_pair_set_many(ps_t *ps, const gb_t *basis, const nelts_t fidx);

/**
 * \brief Gebauer-Moeller checks for product and chain criterion
 *
 * \param pair set ps
 *
 * \param intermediate groebner basis basis
 *
 * \param index in basis of newly added basis element idx
 */
void gebauer_moeller(ps_t *ps, const gb_t *basis,  const nelts_t idx);

/**
 * \brief Remove spairs detected by either product or chain criterion
 *
 * \param pair set ps
 *
 * \param index of basis element the new pairs were generated with idx
 *
 * \return number of removed pairs
 */
nelts_t remove_detected_pairs(ps_t *ps, const nelts_t idx);

/**
 * \brief Generates spair given by one input element. The first generator is 0.
 *
 * \note We use the second generator since then it is easier in symbolic
 * preprocessing: For all spairs, intial ones or normal ones, gen2 has to go to
 * selection sell (lower one). So not so many conditionals have to be done.
 *
 * \param second generator gen2
 *
 * \param input elements input
 *
 * \param hash table ht
 *
 * \return generated spair
 */
void generate_input_element_spair(ps_t *ps, const nelts_t gen1, const gb_t *input,
    ht_t *ht);

/**
 * \brief Generates spair given by gen1 and gen2
 *
 * \param first generator gen1
 *
 * \param second generator gen2
 *
 * \param current groebner basis basis
 *
 * \param hash table ht
 *
 * \return generated spair
 */
static inline void generate_all_spairs(ps_t *ps, const nelts_t new,
    const gb_t *basis, ht_t *ht)
{
  const hash_t hash_new  = basis->p[new][2];
  for (size_t i = basis->st; i < new; ++i) {
    spair_t *sp = ps->pairs + ps->load + i - basis->st;
    /* we have to fix the positions where the new basis element is put (gen2),
     * since we are trying to remove as much as possible useless elements in
     * select_pairs(). if we would dynamically adjust the positioning (as done in
     * the below commented out code) we could no longer track this correctly. */
    sp->gen1  = (nelts_t)i;
    sp->gen2  = new;

    /* if (basis->nt[gen1] < basis->nt[gen2]) {
     *   sp->gen1  = gen1;
     *   sp->gen2  = gen2;
     * } else {
     *   sp->gen1  = gen2;
     *   sp->gen2  = gen1;
     * } */

    sp->lcm   = get_lcm(basis->p[i][2], hash_new, ht);
    /* if (ht->ld[sp->lcm] == 0)
     *   ht->ld[sp->lcm] = new; */

    sp->deg   = ht->deg[sp->lcm];

    /* if one of the generators is redundant we can stop already here and mark it
     * with the CHAIN_CRIT in order to remove it later on */
    /* else */
    if (basis->red[i] > 0) {
      /* printf("%u - %u || %u\n", gen2, gen1, ps->load + gen2 - basis->st); */
      sp->crit  = CHAIN_CRIT;
      /* eliminate(ps->spt[new], i); */
      continue;
    }
    /* check for product criterion and mark correspondingly, i.e. we set sp->deg=0 */
    if (sp->deg == ht->deg[basis->p[i][2]] + ht->deg[hash_new]) {
      sp->crit  = PROD_CRIT;
      /* eliminate(ps->spt[new], i); */
      continue;
    }
    sp->crit  = NO_CRIT;
  }
}
static inline void generate_spair(ps_t *ps, const nelts_t gen1,
    const nelts_t gen2, const gb_t *basis, ht_t *ht)
{
  spair_t *sp = ps->pairs + ps->load + gen2 - basis->st;
  /* we have to fix the positions where the new basis element is put (gen2),
   * since we are trying to remove as much as possible useless elements in
   * select_pairs(). if we would dynamically adjust the positioning (as done in
   * the below commented out code) we could no longer track this correctly. */
  sp->gen1  = gen2;
  sp->gen2  = gen1;

  sp->lcm   = get_lcm(basis->p[gen1][2], basis->p[gen2][2], ht);
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
  if (sp->deg == ht->deg[basis->p[gen1][2]] + ht->deg[basis->p[gen2][2]]) {
    sp->crit  = PROD_CRIT;
    return;
  }
  sp->crit  = NO_CRIT;
  return;
}

/**
 * \brief Updates pair set including Gebauer-Moeller criteria checks
 *
 * \param pair set ps
 *
 * \param intermediate groebner basis gb
 *
 * \param index of new element in gb idx
 */
static inline void update_pair_set(ps_t *ps, const gb_t *basis, const nelts_t idx)
{
  /* we get maximal idx-1 new pairs */
  if (ps->size <= ps->load + (idx-1))
    enlarge_pair_set(ps, 2*ps->size);
  /* generate spairs with the initial elements in basis
   * See note on gb_t in src/types.h why we start at position 1 here. */
  /* gettimeofday(&t_start, NULL); */
  generate_all_spairs(ps, idx, basis, ht);
/*   for (i=basis->st; i<idx; ++i) {
 *     generate_spair(ps, idx, i, basis, ht);
 * #if SPAIR_DEBUG
 *     printf("pair %u, %u + %u | %u\n",idx,i,ps->load,ps->pairs[ps->load+i-basis->st]->deg);
 * #endif
 *   } */
  /* t_gen_spair +=  walltime(t_start); */
  /* printf("gen spairs %9.3f sec\n",
   *     t_gen_spair / (1000000)); */

  /* we do not update ps->load at the moment in order to be able to distinguish
   * old and new pairs for the gebauer-moeller update following */

  /* check product and chain criterion in gebauer moeller style
   * note that we have already marked the pairs for which the product criterion
   * applies in generate_spair() */
  if (idx > basis->st)
    gebauer_moeller(ps, basis, idx);

}

/**
 * \brief Adjusts selection set size to new_size
 *
 * \note This function is not only used for enlargements, but also
 * to cut down memory at the end of symbolic preprocessing
 *
 * \param selection set sel
 *
 *  \param new size new_size
 */
static inline void adjust_size_of_selection(sel_t *sel, const nelts_t new_size)
{
  sel->size = new_size;
}

/**
 * \brief Initializes selection for next reduction step of size size
 *
 * \param size of spairs to generate selection out of size
 *
 * \return selection set sel
 */
static inline sel_t *init_selection(const nelts_t size)
{
  sel_t *sel  = (sel_t *)malloc(sizeof(sel_t));
  sel->size   = size;
  sel->load   = 0;


  return sel;
}

/**
 * \brief Frees selection set
 *
 * \param selection set
 */
static inline void free_selection(sel_t **sel_in)
{
  sel_t *sel  = *sel_in;
  free(sel);
  sel     = NULL;
  *sel_in = sel;
}

/**
 * \brief Comparison implementation for qsort. Sorts pair set w.r.t. graded
 * reverse lexicographical order grevlex.
 *
 * \param element to be compared a
 *
 * \param element to be compared b
 *
 * \returns corresponding integer for qsort
 */
static inline int cmp_spairs_by_grevlex(const void *a, const void *b)
{
  const spair_t *spa  = ((spair_t *)a);
  const spair_t *spb  = ((spair_t *)b);
#if SPAIR_DEBUG
  printf("GEN %u | %u || %u | %u\n",spa->gen1, spa->gen2, spb->gen1, spb->gen2);
  printf("LCM %lu | %lu\n",spa->lcm, spb->lcm);
  printf("DEG %u | %u\n",spa->deg, spb->deg);
#endif
  if (spa->lcm != spb->lcm) {
    /* compare degree */
    if (spa->deg > spb->deg) {
      return 1;
    } else {
      if (spa->deg != spb->deg)
        return -1;
    }
    /* compare reverse lexicographical
     * NOTE: for graded reverse lexicographical ordering we store the exponents
     * ht->exp and ht->ev in reverse order => we can use memcmp() for reverse
     * lex comparison */
    exp_t *expa = ht->exp + (ht->nv * spa->lcm);
    exp_t *expb = ht->exp + (ht->nv * spb->lcm);
    /*
    for (int ii=0; ii<ht->nv; ++ii)
      printf("%u | ",expa[ii]);
    printf("\n");
    for (int ii=0; ii<ht->nv; ++ii)
      printf("%u | ",expb[ii]);
    printf("\n");
    printf("%d\n", memcmp(expb,expa, sizeof(exp_t) * ht->nv));
    */
    return memcmp(expb,expa, sizeof(exp_t) * ht->nv);
  } else {
    return 0;
  }
}

/**
 * \brief Comparison implementation for qsort. Sorts pair set w.r.t. the
 * lexicographical order lex.
 *
 * \note We first sort by degree, then by lex.
 *
 * \param element to be compared a
 *
 * \param element to be compared b
 *
 * \returns corresponding integer for qsort
 */
static inline int cmp_spairs_by_deg_lex(const void *a, const void *b)
{
  const spair_t *spa  = *((spair_t **)a);
  const spair_t *spb  = *((spair_t **)b);
#if SPAIR_DEBUG
  printf("GEN %u | %u || %u | %u\n",spa->gen1, spa->gen2, spb->gen1, spb->gen2);
  printf("LCM %lu | %lu\n",spa->lcm, spb->lcm);
  printf("DEG %u | %u\n",spa->deg, spb->deg);
#endif
  if (spa->lcm != spb->lcm) {
    /* compare degree */
    if (spa->deg > spb->deg) {
      return 1;
    } else {
      if (spa->deg != spb->deg)
        return -1;
    }
    /* compare lexicographical */
    exp_t *expa = ht->exp + (ht->nv * spa->lcm);
    exp_t *expb = ht->exp + (ht->nv * spb->lcm);
    /*
    for (int ii=0; ii<ht->nv; ++ii)
      printf("%u | ",expa[ii]);
    printf("\n");
    for (int ii=0; ii<ht->nv; ++ii)
      printf("%u | ",expb[ii]);
    printf("\n");
    printf("%d\n", memcmp(expb,expa, sizeof(exp_t) * ht->nv));
    */
    return memcmp(expa,expb, sizeof(exp_t) * ht->nv);
  } else {
    return 0;
  }
}

/**
 * \brief Sorts pair set w.r.t. the lexicographical order
 * lex using qsort.
 *
 * \note We first sort by degree, then by lex.
 *
 * \param pair set to be sorted ps
 *
 */
static inline void sort_pair_set_by_lcm_deg_lex(ps_t *ps)
{
  qsort(ps->pairs, ps->load, sizeof(ps), cmp_spairs_by_deg_lex);
}

/**
 * \brief Sorts pair set w.r.t. graded reverse lexicographical order
 * grevlex using qsort.
 *
 * \param pair set to be sorted ps
 *
 */
static inline void sort_pair_set_by_lcm_grevlex(ps_t *ps)
{
  qsort(ps->pairs, ps->load, sizeof(spair_t), cmp_spairs_by_grevlex);
}

/**
 * \brief Selects pairs by lowest degree (normal selection strategy) and returns
 * the index of the last pair in the pair list. Ties are broken using the
 * lexicographical order lex.
 *
 * \note This function also sorts the pair set correspondingly.
 *
 * \param pair set ps
 *
 * \return last index of pair selection in pair list
 */
static inline nelts_t get_pairs_by_minimal_degree_lex(ps_t *ps)
{
  /* sort pair set by lcms */
  sort_pair_set_by_lcm_deg_lex(ps);

  nelts_t i   = 0;
  deg_t dmin  = ps->pairs[0].deg;

  /* we assume here that the pair set is already sorted by degree of the lcms
   * (in particular, we assume grevlex ordering) */
  while (i < ps->load && ps->pairs[i].deg == dmin)
    i++;

  return i;
}


/**
 * \brief Selects pairs by lowest degree (normal selection strategy) and returns
 * the index of the last pair in the pair list. Ties are broken using the graded
 * reverse lexicographical order grevlex.
 *
 * \note This function also sorts the pair set correspondingly.
 *
 * \param pair set ps
 *
 * \return last index of pair selection in pair list
 */
static inline nelts_t get_pairs_by_minimal_degree_grevlex(ps_t *ps)
{
  /* sort pair set by lcms */
  sort_pair_set_by_lcm_grevlex(ps);

  deg_t dmin  = ps->pairs[0].deg;
  nelts_t i   = 0;

  /* we assume here that the pair set is already sorted by degree of the lcms
   * (in particular, we assume grevlex ordering) */
  while (i < ps->load && ps->pairs[i].deg == dmin)
    i++;

  return i;
}
#endif
