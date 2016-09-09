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
 * \file gb.h
 * \brief Input/output routines for matrices
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_GB_H
#define GB_GB_H

#include <math.h>
#include <unistd.h>
#include <limits.h>
#include <getopt.h>
#include "io.h"
#include <src/poly.h>
#include <src/spair.h>
#include <src/symbol.h>
#include <src/matrix.h>

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#ifndef GB_DEBUG
#define GB_DEBUG 0
#endif

/**
 * \brief Prints help for gb call.
 */
void print_help();

/**
 * \brief Simplify functions are stored a structure of function pointers in the basis
 *
 * \param intermediate groebner basis basis
 */
static inline void set_simplify_functions(mp_cf4_ht_t *ht, const gb_t *basis)
{
  switch (basis->sl) {
    // graded reverse lexicographical order
    case 0:
      ht->sf.simplify = no_simplify;
      break;
    // lexicographical order
    case 1:
      ht->sf.simplify = try_to_simplify;
      break;
    case 2:
      ht->sf.simplify = try_to_simplify;
      break;
    default:
      abort ();
  }
}

/**
 * \brief Sort functions are stored a structure of function pointers in the hash
 * table. Depening on the chosen monomial order the functions pointers are set.
 *
 * \param hash table ht
 */
static inline void set_sort_functions_depending_on_monomial_order(mp_cf4_ht_t *ht, const ord_t ord)
{
  switch (ord) {
    // graded reverse lexicographical order
    case 0:
      ht->sort.get_pairs_by_minimal_degree              = get_pairs_by_minimal_degree_grevlex;
      ht->sort.sort_presorted_columns                   = sort_presorted_columns_by_grevlex;
      ht->sort.sort_presorted_columns_invert_left_side  = sort_presorted_columns_by_grevlex_invert_left_side;
      ht->sort.compare_spairs                           = cmp_spairs_by_grevlex;
      ht->sort.compare_monomials                        = cmp_symbolic_preprocessing_monomials_by_grevlex;
      ht->sort.compare_monomials_inverse                = cmp_symbolic_preprocessing_monomials_by_inverse_grevlex;
      ht->sort.compare_polynomials                      = cmp_polynomials_by_grevlex;
      ht->sort.compare_polynomials_inverse              = cmp_polynomials_by_inverse_grevlex;
      break;
    // lexicographical order
    case 1:
      ht->sort.get_pairs_by_minimal_degree              = get_pairs_by_minimal_degree_lex;
      ht->sort.sort_presorted_columns                   = sort_presorted_columns_by_lex;
      ht->sort.sort_presorted_columns_invert_left_side  = sort_presorted_columns_by_lex_invert_left_side;
      ht->sort.compare_spairs                           = cmp_spairs_by_deg_lex;
      ht->sort.compare_monomials                        = cmp_symbolic_preprocessing_monomials_by_lex;
      ht->sort.compare_monomials_inverse                = cmp_symbolic_preprocessing_monomials_by_inverse_lex;
      ht->sort.compare_polynomials                      = cmp_polynomials_by_lex;
      ht->sort.compare_polynomials_inverse              = cmp_polynomials_by_inverse_lex;
      break;
    default:
      abort ();
  }
}

/**
 * \brief Updates basis and pair set after reducing current gbla matrix.
 *
 * \param intermediate groebner basis basis
 *
 * \param pair set ps
 *
 * \param symbolic preprocessing data structure spd
 *
 * \param already reduced gbla matrix mat
 *
 * \param hash table ht
 *
 * \param rank of reduced D part of gbla matrix rankDR
 *
 * \return returns 1 if we have added the constant 1 to the groebner basis, i.e.
 * then the computation is done; else it returns 0.
 */
static inline int update_basis(gb_t *basis, ps_t *ps, spd_t *spd, const mat_t *mat,
    const mp_cf4_ht_t *ht,  const ri_t rankDR)
{
  ri_t i;
  int res;
  for (i=0; i<rankDR; ++i) {
    // add lowest row first, it has the smallest new lead monomial
    res = add_new_element_to_basis(basis, mat, rankDR-1-i, spd, ht);
    // if hash value 0 is new lead monomial we are done, since we have found a
    // unit in the basis, i.e. basis = { 1 }
    if (res == -1)
      continue;
    if (res == 0)
      return 1;
    //printf("psl before generating with row %u: %u\n", rankDR-1-i, ps->load);
    update_pair_set(ps, basis, basis->load-1);
    //printf("psl after: %u\n", ps->load);
    // if elements are homogeneous we compute by degree, thus no redundancy can
    // appear
    if (basis->hom == 0)
      track_redundant_elements_in_basis(basis);
  }
  // track load of basis at the end of this step
  basis->load_ls  = basis->load;
  return 0;
}


/**
 * \brief Updates basis and pair set after reducing current gbla matrix.
 * Moreover, it adds simplifier to simplification list for further optimizations
 * in upcoming symbolic preprocessing steps.
 *
 * \note If nthreads > 1 it tries to update the basis and to add simplifier in
 * parallel.
 *
 * \param intermediate groebner basis basis
 *
 * \param pair set ps
 *
 * \param symbolic preprocessing data structure spd
 *
 * \param already reduced gbla matrix mat
 *
 * \param hash table ht
 *
 * \param rank of reduced D part of gbla matrix rankDR
 *
 * \return returns 1 if we have added the constant 1 to the groebner basis, i.e.
 * then the computation is done; else it returns 0.
 */
int update_basis_and_add_simplifier(gb_t *basis, gb_t *sf, ps_t *ps,
    spd_t *spd, mat_t *mat, const mp_cf4_ht_t *ht,  const ri_t rankDR,
    const int nthreads);

/**
 * \brief Adds simplifier elements to sf list for further exchanges as better
 * spair generators or reducers during upcoming symbolic preprocessing. For this
 * we use the rows from AB in mat.
 *
 * \param intermediate groebner basis basis
 * \param simplifier list sf
 *
 * \param reduced gbla matrix mat
 *
 * \param symbolic preprocessing data structure spd
 *
 * \param hash table ht
 */
void add_simplifier(gb_t *basis, gb_t *sf, mat_t *mat, const spd_t *spd,
    const mp_cf4_ht_t *ht);
#endif
