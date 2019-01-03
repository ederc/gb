/* gb: Gröbner Basis
 * Copyright (C) 2018 Christian Eder <ederc@mathematik.uni-kl.de>
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
 * \file convert.h
 * \brief Implementation of hash value / column index conversions
 * for the exponent hashes of polynomials resp. rows.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_CONVERT_H
#define GB_CONVERT_H

#include "data.h"
#include "basis.h"

hd_t **convert_hashes_to_columns(
        mat_t *mat,
        ht_t *ht,
        stat_t *st
        );

void convert_sparse_matrix_rows_to_basis_elements(
        mat_t *mat,
        bs_t *bs,
        ht_t *ht,
        hd_t **hcm,
        stat_t *st
        );
#endif
