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
int update_basis(gb_t *basis, ps_t *ps, spd_t *spd, const mat_t *mat,
    const mp_cf4_ht_t *ht,  const ri_t rankDR);

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
void add_simplifier_grevlex(gb_t *basis, gb_t *sf, mat_t *mat, const spd_t *spd,
    const mp_cf4_ht_t *ht);
#endif
