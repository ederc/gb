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

spd_t *symbolic_preprocessing(ps_t *ps, const gb_t *basis, const gb_t *sf)
{
  nelts_t i, k, idx, last_div, nsel;
  hash_t hash_pos, hash_div = 0;

  // clears hash table index: there we store during symbolic preprocessing a 2
  // if it is a lead monomial and 1 if it is not a lead monomial. all other
  // entries keep 0, thus they are not part of this reduction step
  clear_hash_table_idx(ht);

  nsel  = get_pairs_by_minimal_degree(ps);

  // list of monomials that appear in the matrix
  pre_t *mon      = init_preprocessing_hash_list(2*nsel);
  // the lower part of the gbla matrix resp. the selection is fixed:
  // those are just the second generators of the spairs, thus we need nsel
  // places.
  sel_t *sel_low  = init_selection(nsel);
  sel_t *sel_upp  = init_selection(5*nsel);
  sel_upp->deg    = ps->pairs[0]->deg;
  sel_low->deg    = ps->pairs[0]->deg;
  // list of polynomials and their multipliers
  select_pairs(ps, sel_upp, sel_low, mon, basis, sf, nsel);

  // we use mon as LIFO: last in, first out. thus we can easily remove and add
  // new elements to mon
  idx = 0;
  while (idx < mon->load) {
    hash_pos  = mon->hpos[idx];
    // only if not already a lead monomial, e.g. if coming from spair
    if (ht->idx[hash_pos] != 2) {
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
      // only if i > 0 we have found a reducer.
      // note: all reducers are added to the upper selection list!
      if (i >= basis->st) {
        if (i == last_div)
          hash_div  = monomial_division(hash_pos, basis->eh[i][0], ht);
        mon->nlm++;
        ht->div[hash_pos]  = i;
        // if multiple is not already in the selected list
        // we have found another element with such a monomial, since we do not
        // take care of the lead monomial below when entering the other lower
        // order monomials, we have to adjust the idx for this given monomial
        // here.
        // we have reducer, i.e. the monomial is a leading monomial (important for
        // splicing matrix later on
        ht->idx[hash_pos] = 2;
        if (sel_upp->load == sel_upp->size)
          adjust_size_of_selection(sel_upp, 2*sel_upp->size);
        sel_upp->mpp[sel_upp->load].mlm = hash_pos;
        sel_upp->mpp[sel_upp->load].mul = hash_div;
        sel_upp->mpp[sel_upp->load].idx = i;
        sel_upp->load++;

        // now add new monomials to preprocessing hash list
        for (k=1; k<basis->nt[i]; ++k) {
          enter_monomial_to_preprocessing_hash_list(sel_upp->mpp[sel_upp->load-1].mul,
              basis->eh[sel_upp->mpp[sel_upp->load-1].idx][k], mon);
        }
      }
    }
    idx++;
  }

  // next we store the information needed to construct the GBLA matrix in the
  // following
  spd_t *mat  = (spd_t *)malloc(sizeof(spd_t));

  // adjust memory
  adjust_size_of_selection(sel_upp, sel_upp->load);
  adjust_size_of_preprocessing_hash_list(mon, mon->load);

  mat->selu = sel_upp;
  mat->sell = sel_low;
  mat->col  = mon;

  return mat;
}

inline void enter_not_multiplied_monomial_to_preprocessing_hash_list(const hash_t h1,
    pre_t *mon)
{
  hash_t pos  = h1;;
  // only in this case we have this monomial hash for the first time,
  // otherwise it has already been taken care of
  if (ht->idx[pos] == 0) {
    ht->idx[pos]++;
    mon->hpos[mon->load]  = pos;
#if SYMBOL_DEBUG
    for (int i=0; i<ht->nv; ++i)
      printf("%u ",ht->exp[h1][i]);
    printf("\n");
    for (int i=0; i<ht->nv; ++i)
      printf("%u ",ht->exp[pos][i]);
    printf("\n");
    printf("new mon %u == %u\n", h1,mon->hpos[mon->load]);
#endif
    mon->load++;
    if (mon->load == mon->size)
      adjust_size_of_preprocessing_hash_list(mon, 2*mon->size);
  }
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
    for (int i=0; i<ht->nv; ++i)
      printf("%u ",ht->exp[h1][i]);
    printf("\n");
    for (int i=0; i<ht->nv; ++i)
      printf("%u ",ht->exp[h2][i]);
    printf("\n");
    for (int i=0; i<ht->nv; ++i)
      printf("%u ",ht->exp[pos][i]);
    printf("\n");
    printf("2 new mon %u + %u == %u\n", h1,h2,mon->hpos[mon->load]);
#endif
    mon->load++;
    if (mon->load == mon->size)
      adjust_size_of_preprocessing_hash_list(mon, 2*mon->size);
  }
}

inline pre_t *init_preprocessing_hash_list(const nelts_t size)
{
  // allocate a list for hashes of monomials to be checked in the symbolic
  // preprocessing
  pre_t *mon  = (pre_t *)malloc(sizeof(pre_t));
  mon->hpos   = (hash_t *)malloc(size * sizeof(hash_t));
  mon->size   = size;
  mon->load   = 0;
  mon->nlm    = 0;
  
  return mon;
}

inline void adjust_size_of_preprocessing_hash_list(pre_t *hl, const nelts_t size)
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
  free_selection(spd->selu);
  free_selection(spd->sell);
  free(spd);
  spd = NULL;
}

inline int cmp_symbolic_preprocessing_monomials_by_lead(const void *a,
    const void *b)
{
  hash_t h1 = *((hash_t *)a);
  hash_t h2 = *((hash_t *)b);

  return (ht->idx[h2] - ht->idx[h1]);
}

inline int cmp_symbolic_preprocessing_monomials_by_grevlex(const void *a,
    const void *b)
{
  hash_t ha = *((hash_t *)a);
  hash_t hb = *((hash_t *)b);

  // compare degree first
  if (ht->deg[hb] > ht->deg[ha]) {
    return 1;
  } else {
    if (ht->deg[ha] > ht->deg[hb])
      return -1;
  }

  // else we have to check reverse lexicographical
  // NOTE: We store the exponents in reverse order in ht->exp and ht->ev
  // => we can use memcmp() here and still get reverse lexicographical ordering
#if __GB_HAVE_SSE2
  nvars_t i;
  exp_t expa[ht->nev * ht->vl];
  exp_t expb[ht->nev * ht->vl];
  for (i=0; i<ht->nev; ++i) {
    _mm_storeu_si128((exp_v *)expa + i*ht->vl, ht->ev[ha][i]);
    _mm_storeu_si128((exp_v *)expb + i*ht->vl, ht->ev[hb][i]);
  }
  return memcmp(expa, expb, ht->nv);
#else
  return memcmp(ht->exp[ha], ht->exp[hb], ht->nv);
#endif
}

inline int cmp_symbolic_preprocessing_monomials_by_inverse_grevlex(const void *a,
    const void *b)
{
  hash_t ha = *((hash_t *)a);
  hash_t hb = *((hash_t *)b);

  // compare degree first
  if (ht->deg[hb] > ht->deg[ha]) {
    return -1;
  } else {
    if (ht->deg[ha] > ht->deg[hb])
      return 1;
  }

  // else we have to check reverse lexicographical
#if __GB_HAVE_SSE2
  nvars_t i;
  exp_t expa[ht->nev * ht->vl];
  exp_t expb[ht->nev * ht->vl];
  for (i=0; i<ht->nev; ++i) {
    _mm_storeu_si128((exp_v *)expa + i*ht->vl, ht->ev[ha][i]);
    _mm_storeu_si128((exp_v *)expb + i*ht->vl, ht->ev[hb][i]);
  }
  return memcmp(expb, expa, ht->nv);
#else
  return memcmp(ht->exp[hb], ht->exp[ha], ht->nv);
#endif
}

inline void sort_columns_by_lead(spd_t *spd)
{
  qsort(spd->col->hpos, spd->col->load, sizeof(hash_t),
      cmp_symbolic_preprocessing_monomials_by_lead);
}

inline void sort_lead_columns_by_inverse_grevlex(spd_t *spd)
{
  if (spd->col->nlm != 0) {
    // sort the start of spd->col, i.e. the lead monomial list
    qsort(spd->col->hpos, spd->col->nlm, sizeof(hash_t),
        cmp_symbolic_preprocessing_monomials_by_inverse_grevlex);
  }
}

inline void sort_non_lead_columns_by_grevlex(spd_t *spd)
{
  // sort the end of spd->col, i.e. the non lead monomial list
  qsort(spd->col->hpos+spd->col->nlm, (spd->col->load - spd->col->nlm),
      sizeof(hash_t), cmp_symbolic_preprocessing_monomials_by_grevlex);
}

inline void sort_presorted_columns_by_grevlex(spd_t *spd, const int nthreads)
{
  #pragma omp parallel num_threads(nthreads)
  {
    #pragma omp single
    {
      #pragma omp task
      sort_lead_columns_by_inverse_grevlex(spd);
      #pragma omp task
      sort_non_lead_columns_by_grevlex(spd);
      #pragma omp taskwait
    }
  }
}

inline void set_column_index_in_hash_table(mp_cf4_ht_t *ht, const spd_t *spd)
{
  nelts_t i;

  for (i=0; i<spd->col->load; ++i)
    ht->idx[spd->col->hpos[i]]  = i;
}

inline int cmp_monomial_polynomial_pair(const void *a, const void *b)
{
  hash_t h1 = ((mpp_t *)a)->mlm;
  hash_t h2 = ((mpp_t *)b)->mlm;

  return (ht->idx[h1] - ht->idx[h2]);

}

inline void sort_selection_by_column_index(spd_t *spd, const mp_cf4_ht_t *ht,
  const int nthreads)
{
  #pragma omp parallel num_threads(nthreads)
  {
    #pragma omp single
    {
      // upper selection
      #pragma omp task
      qsort(spd->selu->mpp, spd->selu->load, sizeof(mpp_t),
          cmp_monomial_polynomial_pair);
      // lower selection
      #pragma omp task
      qsort(spd->sell->mpp, spd->sell->load, sizeof(mpp_t),
          cmp_monomial_polynomial_pair);
      #pragma omp taskwait
    }
  }
}
