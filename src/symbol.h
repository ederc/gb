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
 * \file symbol.h
 * \brief Implementation of the symbolic pre- and postprocessing in gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_SYMBOL_H
#define GB_SYMBOL_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#if defined(_OPENMP)
#include <omp.h>
#endif
#include <config.h>
#include "types.h"
#include "hash.h"
#include "spair.h"

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#ifndef SYMBOL_DEBUG
#define SYMBOL_DEBUG 0
#endif

#ifndef HASH_CHECK
#define HASH_CHECK 0
#endif

unsigned long load_global = 0;

/**
 * \brief Symbolic preprocessing searching for all possible reducers of all
 * ocurring monomials. These data are then used to construct the matrices for
 * gbla.
 */
void symbolic_preprocessing(ps_t *ps, smc_t *AB, smc_t *CD,
    pre_t *mon, const gb_t *basis);

/**
 * \brief Adjusts size of hash list for symbolic preprocessing to new_size.
 *
 * \note It is not only used to enlarge the size, but also for cutting down
 * memory once symbolic preprocessing is done.
 *
 * \param preprocessing hash list hl
 *
 * \param new size of hash list size
 */
static inline void adjust_size_of_preprocessing_hash_list(pre_t *hl, const nelts_t size)
{
  /* printf("hl->size %u --> %u\n", hl->size, size); */
  hl->hash  = realloc(hl->hash, size * sizeof(hash_t));
  hl->size  = size;
}

/**
 * \brief Initializes a hash list for symbolic preprocessing.
 *
 * \param size of hash list size
 *
 * \return hash list
 */
static inline pre_t *init_preprocessing_hash_list(const nelts_t size)
{
  /* allocate a list for hashes of monomials to be checked in the symbolic
   * preprocessing */
  pre_t *mon  = (pre_t *)malloc(sizeof(pre_t));
  mon->hash   = (hash_t *)malloc(size * sizeof(hash_t));
  mon->size   = size;
  mon->load   = 0;
  mon->nlm    = 0;
  
  return mon;
}


/**
 * \brief Frees hash list for symbolic preprocessing.
 *
 * \param preprocessing hash list hl
 */
static inline void free_preprocessing_hash_list(pre_t **hl_in)
{
  pre_t *hl = *hl_in;
  free(hl->hash);
  free(hl);
  hl      = NULL;
  *hl_in  = hl;
}


/**
 * \brief Frees data structure for symbolic preprocessing.
 *
 * \param symbolic preprocessing data structure spd
 */
static inline void free_symbolic_preprocessing_data(spd_t **spd_in)
{
  spd_t *spd = *spd_in;
  free_preprocessing_hash_list(&(spd->col));
  free_selection(&(spd->selu));
  free_selection(&(spd->sell));
  free(spd);
  spd     = NULL;
  *spd_in = spd;
}


/**
 * \brief Comparison function of monomials for quicksort. Compares the property
 * of being lead monomial or not.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
static inline int cmp_symbolic_preprocessing_monomials_by_lead(const void *a,
    const void *b)
{
  hash_t h1 = *((hash_t *)a);
  hash_t h2 = *((hash_t *)b);

  return (int)(ht->idx[h2] - ht->idx[h1]);
}

/**
 * \brief Comparison function of monomials for quicksort. Compares w.r.t. the
 * given monomial order lex.
 *
 * \note This is used in symbolic preprocessing, thus we know already that a =/=
 * b and we do not have to check this again.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
static inline int cmp_symbolic_preprocessing_monomials_by_lex(const void *a,
    const void *b)
{
  const hash_t ha = *((hash_t *)a);
  const hash_t hb = *((hash_t *)b);

  /* else we have to check lexicographical */
  return memcmp(ht->exp+(ht->nv * hb), ht->exp+(ht->nv * ha), sizeof(exp_t) * ht->nv);
}

/**
 * \brief Comparison function of monomials for quicksort. Compares w.r.t. the
 * given monomial order lex and inverts the order.
 *
 * \note This is used in symbolic preprocessing, thus we know already that a =/=
 * b and we do not have to check this again.
 *
 * \note We invert lex since we have to align the A and C parts of the gbla
 * matrix in order to efficiently compute.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
static inline int cmp_symbolic_preprocessing_monomials_by_inverse_lex(const void *a,
    const void *b)
{
  const hash_t ha = *((hash_t *)a);
  const hash_t hb = *((hash_t *)b);

  /* else we have to check lexicographical */
  return memcmp(ht->exp+(ht->nv * ha), ht->exp+(ht->nv * hb), sizeof(exp_t) * ht->nv);
}

/**
 * \brief Comparison function of monomials for quicksort. Compares w.r.t. the
 * given monomial order grevlex.
 *
 * \note This is used in symbolic preprocessing, thus we know already that a =/=
 * b and we do not have to check this again.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
static inline int cmp_symbolic_preprocessing_monomials_by_grevlex(const void *a,
    const void *b)
{
  const hash_t ha = *((hash_t *)a);
  const hash_t hb = *((hash_t *)b);

  /* compare degree first */
  if (ht->deg[hb] != ht->deg[ha])
    return (int)(ht->deg[hb]-ht->deg[ha]);

  /* else we have to check reverse lexicographical
   * NOTE: We store the exponents in reverse order in ht->exp and ht->ev
   * => we can use memcmp() here and still get reverse lexicographical ordering */
  return memcmp(ht->exp+(ht->nv * ha), ht->exp+(ht->nv * hb), sizeof(exp_t) * ht->nv);
}

/**
 * \brief Comparison function of monomials for quicksort. Compares w.r.t. the
 * given monomial order grevlex and inverts the order.
 *
 * \note This is used in symbolic preprocessing, thus we know already that a =/=
 * b and we do not have to check this again.
 *
 * \note We invert grevlex since we have to align the A and C parts of the gbla
 * matrix in order to efficiently compute.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
static inline int cmp_symbolic_preprocessing_monomials_by_inverse_grevlex(const void *a,
    const void *b)
{
  const hash_t ha = *((hash_t *)a);
  const hash_t hb = *((hash_t *)b);

  /* compare degree first */
  if (ht->deg[hb] != ht->deg[ha])
    return (int)(ht->deg[ha]-ht->deg[hb]);

  /* else we have to check reverse lexicographical */
  return memcmp(ht->exp+(ht->nv * hb), ht->exp+(ht->nv * ha), sizeof(exp_t) * ht->nv);
}

/**
 * \brief Comparison function of polynomials for quicksort. Compares w.r.t. the
 * given monomial order lex.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
static inline int cmp_polynomials_by_lex(const void *a,
    const void *b)
{
  const poly_t *pa = *((poly_t **)a);
  const poly_t *pb = *((poly_t **)b);

  const hash_t ha = pa[2];
  const hash_t hb = pb[2];

  return memcmp(ht->exp+(ht->nv * hb), ht->exp+(ht->nv * ha), sizeof(exp_t) * ht->nv);
}

/**
 * \brief Comparison function of polynomials for quicksort. Compares w.r.t. the
 * given monomial order grevlex.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
static inline int cmp_polynomials_by_grevlex(const void *a,
    const void *b)
{
  const poly_t *pa = *((poly_t **)a);
  const poly_t *pb = *((poly_t **)b);

  const hash_t ha = pa[2];
  const hash_t hb = pb[2];

  /* compare degree first */
  if (ht->deg[hb] != ht->deg[ha])
    return (int)(ht->deg[hb]-ht->deg[ha]);

  /* else we have to check reverse lexicographical
   * NOTE: We store the exponents in reverse order in ht->exp and ht->ev
   * => we can use memcmp() here and still get reverse lexicographical ordering */
  return memcmp(ht->exp+(ht->nv * ha), ht->exp+(ht->nv * hb), sizeof(exp_t) * ht->nv);
}

/**
 * \brief Comparison function of polynomials for quicksort. Compares inverted 
 * w.r.t. the given monomial order lex.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
static inline int cmp_polynomials_by_inverse_lex(const void *a,
    const void *b)
{
  const poly_t *pa = *((poly_t **)a);
  const poly_t *pb = *((poly_t **)b);

  const hash_t ha = pa[2];
  const hash_t hb = pb[2];

  return memcmp(ht->exp+(ht->nv * ha), ht->exp+(ht->nv * hb), sizeof(exp_t) * ht->nv);
}

/**
 * \brief Comparison function of polynomials for quicksort. Compares inverted
 * w.r.t. the given monomial order grevlex.
 *
 * \param value a
 *
 * \param value b
 *
 * \returns negative value if a is non lead and b is lead; 0 if both are lead or
 * both are non lead; a positive value if a is lead and b is non lead
 */
static inline int cmp_polynomials_by_inverse_grevlex(const void *a,
    const void *b)
{
  const poly_t *pa = *((poly_t **)a);
  const poly_t *pb = *((poly_t **)b);

  const hash_t ha = pa[2];
  const hash_t hb = pb[2];

  /* compare degree first */
  if (ht->deg[hb] != ht->deg[ha])
    return (int)(ht->deg[ha]-ht->deg[hb]);

  /* else we have to check reverse lexicographical
   * NOTE: We store the exponents in reverse order in ht->exp and ht->ev
   * => we can use memcmp() here and still get reverse lexicographical ordering */
  return memcmp(ht->exp+(ht->nv * hb), ht->exp+(ht->nv * ha), sizeof(exp_t) * ht->nv);
}

/**
 * \brief Sorts columns resp. monomials found by symbolic preprocessing to get
 * two different parts: Monomials which are lead monomials and monomials which
 * are non lead monomials.
 */
static inline void sort_columns_by_lead(pre_t *mon)
{
  qsort(mon->hash, mon->load, sizeof(hash_t),
      cmp_symbolic_preprocessing_monomials_by_lead);
}

/**
 * \brief Sorts lead monomials found by symbolic preprocessing w.r.t. the
 * given monomial order lex
 *
 * \note The list of columns resp. monomials was already presorted to
 * distinguish lead ones (at the beginning of the list) from non lead ones (at
 * the end of the list).
 *
 * \param symbolic preprocessing data spd
 */
static inline void sort_lead_columns_by_lex(pre_t *mon)
{
  if (mon->nlm != 0) {
    /* sort the start of spd->col, i.e. the lead monomial list */
    qsort(mon->hash, mon->nlm, sizeof(hash_t),
        cmp_symbolic_preprocessing_monomials_by_lex);
  }
}

/**
 * \brief Sorts lead monomials found by symbolic preprocessing w.r.t. the
 * given monomial order grevlex
 *
 * \note The list of columns resp. monomials was already presorted to
 * distinguish lead ones (at the beginning of the list) from non lead ones (at
 * the end of the list).
 *
 * \param symbolic preprocessing data spd
 */
static inline void sort_lead_columns_by_grevlex(pre_t *mon)
{
  if (mon->nlm != 0) {
    /* sort the start of spd->col, i.e. the lead monomial list */
    qsort(mon->hash, mon->nlm, sizeof(hash_t),
        cmp_symbolic_preprocessing_monomials_by_grevlex);
  }
}

static inline void sort_columns_by_grevlex(pre_t *mon)
{
  qsort(mon->hash, mon->load, sizeof(hash_t),
      cmp_symbolic_preprocessing_monomials_by_grevlex);
}

static inline void sort_columns_by_lex(pre_t *mon)
{
  qsort(mon->hash, mon->load, sizeof(hash_t),
      cmp_symbolic_preprocessing_monomials_by_lex);
}

/**
 * \brief Sorts lead monomials found by symbolic preprocessing w.r.t. the
 * given monomial order lex and inverts the order
 *
 * \note The list of columns resp. monomials was already presorted to
 * distinguish lead ones (at the beginning of the list) from non lead ones (at
 * the end of the list).
 *
 * \note We invert lex since we have to align the A and C parts of the gbla
 * matrix in order to efficiently compute.
 *
 * \param symbolic preprocessing data spd
 */
static inline void sort_lead_columns_by_inverse_lex(pre_t *mon)
{
  if (mon->nlm != 0) {
    /* sort the start of spd->col, i.e. the lead monomial list */
    qsort(mon->hash, mon->nlm, sizeof(hash_t),
        cmp_symbolic_preprocessing_monomials_by_inverse_lex);
  }
}

/**
 * \brief Sorts lead monomials found by symbolic preprocessing w.r.t. the
 * given monomial order grevlex and inverts the order
 *
 * \note The list of columns resp. monomials was already presorted to
 * distinguish lead ones (at the beginning of the list) from non lead ones (at
 * the end of the list).
 *
 * \note We invert grevlex since we have to align the A and C parts of the gbla
 * matrix in order to efficiently compute.
 *
 * \param symbolic preprocessing data spd
 */
static inline void sort_lead_columns_by_inverse_grevlex(pre_t *mon)
{
  if (mon->nlm != 0) {
    /* sort the start of spd->col, i.e. the lead monomial list */
    qsort(mon->hash, mon->nlm, sizeof(hash_t),
        cmp_symbolic_preprocessing_monomials_by_inverse_grevlex);
  }
}

/**
 * \brief Sorts non lead monomials found by symbolic preprocessing w.r.t. the
 * given monomial order lex
 *
 * \note The list of columns resp. monomials was already presorted to
 * distinguish lead ones (at the beginning of the list) from non lead ones (at
 * the end of the list).
 *
 * \param symbolic preprocessing data spd
 */
static inline void sort_non_lead_columns_by_lex(pre_t *mon)
{
  /* sort the end of spd->col, i.e. the non lead monomial list */
  qsort(mon->hash+mon->nlm, (mon->load - mon->nlm),
      sizeof(hash_t), cmp_symbolic_preprocessing_monomials_by_lex);
}

/**
 * \brief Sorts non lead monomials found by symbolic preprocessing w.r.t. the
 * given monomial order grevlex
 *
 * \note The list of columns resp. monomials was already presorted to
 * distinguish lead ones (at the beginning of the list) from non lead ones (at
 * the end of the list).
 *
 * \param symbolic preprocessing data spd
 */
static inline void sort_non_lead_columns_by_grevlex(pre_t *mon)
{
  /* sort the end of spd->col, i.e. the non lead monomial list */
  qsort(mon->hash+mon->nlm, (mon->load - mon->nlm),
      sizeof(hash_t), cmp_symbolic_preprocessing_monomials_by_grevlex);
}

/**
 * \brief Sorts already presorted list of monomials resp. columns w.r.t. the
 * given monomial order lex. The list is already presorted by lead and non
 * lead monomials. Those monomials correspond then to the columns of the gbla
 * matrix generated later on.
 */
static inline void sort_presorted_columns_by_lex(pre_t *mon,
    const int nthreads)
{
  const int t = 2<nthreads ? 2 : nthreads;
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single
    {
      #pragma omp task
      sort_lead_columns_by_lex(mon);
      #pragma omp task
      sort_non_lead_columns_by_lex(mon);
      #pragma omp taskwait
    }
  }
}

/**
 * \brief Sorts already presorted list of monomials resp. columns w.r.t. the
 * given monomial order grevlex. The list is already presorted by lead and non
 * lead monomials. Those monomials correspond then to the columns of the gbla
 * matrix generated later on.
 */
static inline void sort_presorted_columns_by_grevlex(pre_t *mon,
    const int nthreads)
{
  const int t = 2<nthreads ? 2 : nthreads;
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single
    {
      #pragma omp task
      sort_lead_columns_by_grevlex(mon);
      #pragma omp task
      sort_non_lead_columns_by_grevlex(mon);
      #pragma omp taskwait
    }
  }
}

/**
 * \brief Sorts already presorted list of monomials resp. columns w.r.t. the
 * given monomial order lex. The list is already presorted by lead and non
 * lead monomials. Those monomials correspond then to the columns of the gbla
 * matrix generated later on.
 *
 * \note On the lefthand side of the gbla matrix we invert the sorting of the
 * columns due to gbla's internal ordering for reducing A later on.
 */
static inline void sort_presorted_columns_by_lex_invert_left_side(pre_t *mon,
    const int nthreads)
{
  const int t = 2<nthreads ? 2 : nthreads;
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single
    {
      #pragma omp task
      sort_lead_columns_by_inverse_lex(mon);
      #pragma omp task
      sort_non_lead_columns_by_lex(mon);
      #pragma omp taskwait
    }
  }
}

/**
 * \brief Sorts already presorted list of monomials resp. columns w.r.t. the
 * given monomial order grevlex. The list is already presorted by lead and non
 * lead monomials. Those monomials correspond then to the columns of the gbla
 * matrix generated later on.
 *
 * \note On the lefthand side of the gbla matrix we invert the sorting of the
 * columns due to gbla's internal ordering for reducing A later on.
 */
static inline void sort_presorted_columns_by_grevlex_invert_left_side(pre_t *mon,
    const int nthreads)
{
  const int t = 2<nthreads ? 2 : nthreads;
  #pragma omp parallel num_threads(t)
  {
    #pragma omp single
    {
      #pragma omp task
      sort_lead_columns_by_inverse_grevlex(mon);
      #pragma omp task
      sort_non_lead_columns_by_grevlex(mon);
      #pragma omp taskwait
    }
  }
}

/**
 * \brief Sets the index entry of the hash values to the corresponding columns
 * for the upcoming gbla matrix. The relation between hash values of monomials
 * and columns is given by the symbolic preprocessing data.
 */
static inline void set_column_index_in_hash_table(ht_t *ht, const pre_t *mon)
{
  for (size_t i = 0; i < mon->load; ++i)
    ht->idx[mon->hash[i]]  = (ht_size_t)i;
}

/**
 * \brief Selects pairs due to the sorting and selection done by the previous
 * get_pairs_by_*() procedure and returns
 * a set of this selection that is later on filled with corresponding lower term
 * reducers.
 *
 * \note The function also already marks the hash values for the lcms of the
 * spairs with "2". So we know that there are 2 polynomials in the upcoming
 * matrix that hit that monomial. We use this in the following symbolic
 * preprocessing not to handle this monomial again and also, in the matrix
 * construction later on, to know the number of nonzero entries in the column
 * corresponding to this monomial. This will help splicing the matrix for gbla.
 *
 * \note There are two selection lists, one for the upper and one for the lower
 * part of the gbla matrix to be constructed. By default, the first generator of
 * an spair is added to the upper part, the second generator to the lower part.
 * All upcoming reducers are also added to the upper selection list.
 */

static inline void add_to_monomial_list(pre_t *mon, const hash_t h)
{
  if (mon->load >= mon->size) {
    mon->hash =   realloc(mon->hash, 2 * mon->size * sizeof(hash_t));
    mon->size *=  2;
  }
  mon->hash[mon->load] = h;
  mon->load++;
}

static inline void select_pairs(ps_t *ps, smc_t *AB, smc_t *CD,
    pre_t *mon, const gb_t *basis, const nelts_t nsel)
{
  nelts_t i, j ,k, l;

  /* for each lcm we first detect which elements have to be added to the
   * symbolic preprocessing step, i.e. we remove duplicates */
  i = 0;
  hash_t lcm    = 0;
  nelts_t *gens = (nelts_t *)malloc(2 * nsel * sizeof(nelts_t));
  nelts_t load  = 0;
  nelts_t tmp, min_nt_pos;
#define NOT_ADDING_MUL_TO_HT 1
#if NOT_ADDING_MUL_TO_HT
  hash_t mul_hash;
  deg_t mul_deg;
  exp_t *mul_exp  = (exp_t *)malloc(ht->nv * sizeof(exp_t));
#endif
  while (i<nsel) {
    memset(gens, 0, 2 * nsel * sizeof(nelts_t));
    lcm     = ps->pairs[i].lcm;
    gens[0] = ps->pairs[i].gen1;
    j = i;
    while (j<nsel && ps->pairs[j].lcm == lcm) {
      /* check first generator */
      for (k=0; k<load; ++k) {
        if (ps->pairs[j].gen1 == gens[k]) {
          break;
        }
      }
      if (k == load) {
        gens[load]  = ps->pairs[j].gen1;
        load++;
      }
      /* check second generator */
      for (k=0; k<load; ++k) {
        if (ps->pairs[j].gen2 == gens[k]) {
          break;
        }
      }
      if (k == load) {
        gens[load]  = ps->pairs[j].gen2;
        load++;
      }
      j++;
    }

    /* TODO:
     * does it really makes sense to sort the elements such that the known pivot
     * is the corresponding sparsest row? or is the benefit from this not enough
     * to compensate the drawback of searching for the element? */
    /* put generator with minimal number of terms on first position, i.e.
     * sparsest row of given lcm is used as known pivot */
    min_nt_pos  = 0;
    for (nelts_t l = 1; l < load; ++l) {
      if (basis->p[gens[min_nt_pos]][1] > basis->p[gens[l]][1])
        min_nt_pos  = l;
    }
    /* put generator with minimal number of terms on first position */
    tmp               = gens[0];
    gens[0]           = gens[min_nt_pos];
    gens[min_nt_pos]  = tmp;
    i     = j;
    /* now we have handled all poairs of the given lcm, we add them to the
     * symbolic preprocessing step in the following */
    k = 0;

    /* add to upper matrix AB */
    exp_t *ev = ht->exp + (ht->nv * lcm);
    exp_t *ed = NULL;
    if (load > 1) {
      AB->row[AB->rk] = (src_t *)malloc(basis->p[gens[k]][1] * sizeof(src_t));
      memcpy(AB->row[AB->rk], basis->p[gens[k]],
          basis->p[gens[k]][1] * sizeof(src_t));
#if NOT_ADDING_MUL_TO_HT
      mul_hash = ht->val[lcm] - ht->val[basis->p[gens[k]][2]];
      ed  = ht->exp + (ht->nv * basis->p[gens[k]][2]);
      for (size_t j = 0; j < ht->nv; ++j) {
        mul_exp[j]  = ev[j] - ed[j];
      }
      mul_deg = 0;
      for (size_t j = 0; j < ht->nv; ++j) {
        mul_deg +=  mul_exp[j];
      }
#else
      hash_t mul = get_multiplier(lcm, basis->p[gens[k]][2], ht);
#endif
      /* printf("adds %u with mul %u to AB at row %u || lcm %u\n", gens[k], mul, AB->rk, lcm); */
      for (size_t i = 2; i < AB->row[AB->rk][1]; i = i+2) {
#if NOT_ADDING_MUL_TO_HT
        AB->row[AB->rk][i]  = check_in_hash_table_product_special(
            AB->row[AB->rk][i], mul_hash, mul_deg, mul_exp, ht);
      /* printf("a1 %d\n", AB->row[AB->rk][i]); */
#else
        AB->row[AB->rk][i]  = check_in_hash_table_product(
            mul, AB->row[AB->rk][i], ht);
#endif
        if (ht->idx[AB->row[AB->rk][i]] == 0) {
          ht->idx[AB->row[AB->rk][i]] = 1;
          add_to_monomial_list(mon, AB->row[AB->rk][i]);
          /* printf("a mon->hash[%d] %d\n", mon->load, mon->hash[mon->load]); */
        }
      }
      AB->rk++;
      mon->nlm++;
      /* mark lcm term, only here, otherwise there might not be a second 
       * polynomial as reducer */
      ht->idx[lcm]  = 2;
      k++;
    }
    /* add to lower matrix CD */
    for (l=k; l<load; l++) {
      ed  = ht->exp + (ht->nv * basis->p[gens[l]][2]);
      CD->row[CD->rk] = (src_t *)malloc(basis->p[gens[l]][1] * sizeof(src_t));
      memcpy(CD->row[CD->rk], basis->p[gens[l]],
          basis->p[gens[l]][1] * sizeof(src_t));
#if NOT_ADDING_MUL_TO_HT
      mul_hash = ht->val[lcm] - ht->val[basis->p[gens[l]][2]];
      /* printf("LCM ");
       * for (size_t k = 0; k < ht->nv; ++k) {
       *   printf("%d ", ht->exp[lcm][k]);
       * }
       * printf("\n");
       * printf("POL %u ", gens[l]);
       * for (size_t k = 0; k < ht->nv; ++k) {
       *   printf("%d ", ht->exp[basis->p[gens[l]][2]][k]);
       * }
       * printf("\n"); */
      for (size_t j = 0; j < ht->nv; ++j) {
        mul_exp[j]  = ev[j] - ed[j];
      }
      mul_deg = 0;
      for (size_t j = 0; j < ht->nv; ++j) {
        mul_deg +=  mul_exp[j];
      }
#else
      hash_t mul = get_multiplier(lcm, basis->p[gens[l]][2], ht);
#endif
      /* printf("adds %u with mul %u to CD at row %u || lcm %u\n", gens[l], mul, CD->rk, lcm); */
      for (size_t i = 2; i < CD->row[CD->rk][1]; i = i+2) {
#if NOT_ADDING_MUL_TO_HT
        CD->row[CD->rk][i]  = check_in_hash_table_product_special(
            CD->row[CD->rk][i], mul_hash, mul_deg, mul_exp, ht);
      /* printf("c %d\n", CD->row[CD->rk][i]); */
#else
        CD->row[CD->rk][i]  = check_in_hash_table_product(
            mul, CD->row[CD->rk][i], ht);
#endif
        if (ht->idx[CD->row[CD->rk][i]] == 0) {
          ht->idx[CD->row[CD->rk][i]] = 1;
          add_to_monomial_list(mon, CD->row[CD->rk][i]);
        }
      }
      CD->rk++;
    }
    /* set data for next lcm round */
    load  = 0;
  }
  free(gens);

#if NOT_ADDING_MUL_TO_HT
  free(mul_exp);
#endif

  /* adjust pair set after removing the bunch of selected pairs */
  k = 0;
  /* for (i=0; i<nsel; ++i)
   *   free(ps->pairs[i]); */
  for (i=nsel; i<ps->load; ++i) {
    ps->pairs[k] = ps->pairs[i];
    k++;
  }
  ps->load  = k;
}
#endif
