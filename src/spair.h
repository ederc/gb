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

#include "gb_config.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include "types.h"
#include "hash.h"

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
ps_t *initialize_pair_set(const gb_t *input, mp_cf4_ht_t *ht);

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
spair_t *generate_input_element_spair(const nelts_t gen1, const gb_t *input,
    mp_cf4_ht_t *ht);

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
spair_t *generate_spair(const nelts_t gen1, const nelts_t gen2, const gb_t *basis,
    mp_cf4_ht_t *ht);

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
void gebauer_moeller_new(ps_t *ps, const gb_t *basis);

/**
 * \brief Remove spairs detected by either product or chain criterion
 *
 * \param pair set ps
 *
 * \param intermediate groebner basis basis
 *
 * \param index of basis element the new pairs were generated with idx
 *
 * \return number of removed pairs
 */
static inline nelts_t remove_detected_pairs(ps_t *ps, const gb_t *basis, const nelts_t ctr)
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

/**
 * \brief Updates pair set including Gebauer-Moeller criteria checks
 *
 * \param pair set ps
 *
 * \param intermediate groebner basis gb
 *
 * \param index of new element in gb idx
 */
void update_pair_set(ps_t *ps, const gb_t *basis, const nelts_t idx);

static inline void update_pair_set_new(ps_t *ps, const gb_t *basis)
{
  const nelts_t fidx  = basis->load_ls; // index of first new element
  const nelts_t lidx  = basis->load-1;  // index of last new element
  nelts_t i, j;
  int sum = 0;
  for (j=fidx; j<lidx+1; ++j)
    sum +=  (int)(j-2-fidx) > 0 ? j-2-fidx : 0;

  nelts_t end = 0;
  //printf("idx %u\n", i);
  // we get maximal (lidx-1)*(basis->load_ls-basis->st)+sum new pairs
  while (ps->size <= ps->load + (lidx-1)*(basis->load_ls-basis->st)+sum)
    enlarge_pair_set(ps, 2*ps->size);
  for (i=fidx; i<lidx+1; ++i) {
    // generate spairs with the initial elements in basis
    // See note on gb_t in src/types.h why we start at position 1 here.
    for (j=basis->st; j<i; ++j) {
      ps->pairs[ps->load+end] = generate_spair(i, j, basis, ht);
#if SPAIR_DEBUG
      printf("pair %p %u, %u + %u | %u\n",ps->pairs[ps->load+j-basis->st],i,j,ps->load,ps->pairs[ps->load+j-basis->st]->deg);
#endif
      end++;
    }
    // we do not update ps->load at the moment in order to be able to distinguish
    // old and new pairs for the gebauer-moeller update following

    // check product and chain criterion in gebauer moeller style
    // note that we have already marked the pairs for which the product criterion
    // applies in generate_spair()
  }
  //printf("%u new elements, %u new pairs, %u length of pair set\n", lidx+1-fidx, end, ps->load+end);
  gebauer_moeller_new(ps, basis);
  // fix pair set and remove detected pairs
  meta_data->ncrit_last   =   remove_detected_pairs(ps, basis, end);
  meta_data->ncrit_total  +=  meta_data->ncrit_last;
}


/**
 * \brief Adds generator gen of the corresponding spair with least common
 * multiple lcm to selection list sel.
 *
 * \param selection list sel
 *
 * \param intermediate groebner basis basis
 *
 * \param least common multiple of spair lcm
 *
 * \param a generator of spair gen
 */
void add_spair_generator_to_selection(sel_t *sel, const gb_t *basis,
    const hash_t lcm, const nelts_t gen);

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
  sel->mpp  = realloc(sel->mpp, sel->size * sizeof(mpp_t));
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

  sel->mpp    = (mpp_t *)malloc(size * sizeof(mpp_t));

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
  free(sel->mpp);
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
  spair_t *spa  = *((spair_t **)a);
  spair_t *spb  = *((spair_t **)b);
//#if SPAIR_DEBUG
  //printf("%p | %p\n", spa, spb);
  //printf("nvars %u\n",ht->nv);
  /*
  printf("GEN %u | %u || %u | %u\n",spa->gen1, spa->gen2, spb->gen1, spb->gen2);
  printf("LCM %lu | %lu\n",spa->lcm, spb->lcm);
  printf("DEG %u | %u\n",spa->deg, spb->deg);
  */
//#endif
  if (spa->lcm != spb->lcm) {
    // compare degree
    if (spa->deg > spb->deg) {
      return 1;
    } else {
      if (spa->deg != spb->deg)
        return -1;
    }
    // compare reverse lexicographical
    // NOTE: for graded reverse lexicographical ordering we store the exponents
    // ht->exp and ht->ev in reverse order => we can use memcmp() for reverse
    // lex comparison
#if __GB_HAVE_SSE2
    nvars_t i;
    exp_t expa[ht->nev * ht->vl] __attribute__ ((aligned (16)));
    exp_t expb[ht->nev * ht->vl] __attribute__ ((aligned (16)));
    exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
    for (i=0; i<ht->nev; ++i) {
      _mm_store_si128((exp_v *)tmp, ht->ev[spa->lcm][i]);
      memcpy(expa+(i*ht->vl), tmp, ht->vl*sizeof(exp_t));
      _mm_store_si128((exp_v *)tmp, ht->ev[spb->lcm][i]);
      memcpy(expb+(i*ht->vl), tmp, ht->vl*sizeof(exp_t));
    }
#else
    exp_t *expa = ht->exp[spa->lcm];
    exp_t *expb = ht->exp[spb->lcm];
#endif
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
    // both have the same lcms and are not detected by the product criterion,
    // then we break ties by the overall number of terms
    if (spa->nt != spb->nt) {
      if (spa->nt < spb->nt) {
        return -1;
      } else {
        if (spa->nt != spb->nt) {
          return 1;
        } else {
          if (spa->gen1 < spb->gen1 || spa->gen2 < spb->gen2) {
            return -1;
          }
        }
      }
    }
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
  spair_t *spa  = *((spair_t **)a);
  spair_t *spb  = *((spair_t **)b);
//#if SPAIR_DEBUG
  //printf("%p | %p\n", spa, spb);
  //printf("nvars %u\n",ht->nv);
  /*
  printf("GEN %u | %u || %u | %u\n",spa->gen1, spa->gen2, spb->gen1, spb->gen2);
  printf("LCM %lu | %lu\n",spa->lcm, spb->lcm);
  printf("DEG %u | %u\n",spa->deg, spb->deg);
  */
//#endif
  if (spa->lcm != spb->lcm) {
    // compare degree
    if (spa->deg > spb->deg) {
      return 1;
    } else {
      if (spa->deg != spb->deg)
        return -1;
    }
    // compare lexicographical
#if __GB_HAVE_SSE2
    nvars_t i;
    exp_t expa[ht->nev * ht->vl] __attribute__ ((aligned (16)));
    exp_t expb[ht->nev * ht->vl] __attribute__ ((aligned (16)));
    exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
    for (i=0; i<ht->nev; ++i) {
      _mm_store_si128((exp_v *)tmp, ht->ev[spa->lcm][i]);
      memcpy(expa+(i*ht->vl), tmp, ht->vl*sizeof(exp_t));
      _mm_store_si128((exp_v *)tmp, ht->ev[spb->lcm][i]);
      memcpy(expb+(i*ht->vl), tmp, ht->vl*sizeof(exp_t));
    }
#else
    exp_t *expa = ht->exp[spa->lcm];
    exp_t *expb = ht->exp[spb->lcm];
#endif
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
    // both have the same lcms and are not detected by the product criterion,
    // then we break ties by the overall number of terms
    if (spa->nt != spb->nt) {
      if (spa->nt < spb->nt) {
        return -1;
      } else {
        if (spa->nt != spb->nt) {
          return 1;
        } else {
          if (spa->gen1 < spb->gen1 || spa->gen2 < spb->gen2) {
            return -1;
          }
        }
      }
    }
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
  qsort(ps->pairs, ps->load, sizeof(spair_t **), cmp_spairs_by_deg_lex);
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
  qsort(ps->pairs, ps->load, sizeof(spair_t **), cmp_spairs_by_grevlex);
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
  // sort pair set by lcms
  sort_pair_set_by_lcm_deg_lex(ps);

  nelts_t i   = 0;
  deg_t dmin  = ps->pairs[0]->deg;

  // we assume here that the pair set is already sorted by degree of the lcms
  // (in particular, we assume grevlex ordering)
  while (i < ps->load && ps->pairs[i]->deg == dmin)
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
  // sort pair set by lcms
  sort_pair_set_by_lcm_grevlex(ps);

  deg_t dmin  = ps->pairs[0]->deg;
  nelts_t i   = 0;

  // we assume here that the pair set is already sorted by degree of the lcms
  // (in particular, we assume grevlex ordering)
  while (i < ps->load && ps->pairs[i]->deg == dmin)
    i++;

  return i;
}
#endif
