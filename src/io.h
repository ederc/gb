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
 * \file io.h
 * \brief Input and output handling
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_IO_H
#define GB_IO_H

#include "data.h"
#include "order.h"

void set_function_pointers(
        const stat_t *st
        );
int32_t check_and_set_meta_data(
        stat_t *st;
        const int32_t *lens,
        const int32_t *cfs,
        const int32_t *exps,
        const int32_t field_char,
        const int32_t mon_order,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t regenerate_ht,
        const int32_t la_option,
        const int32_t info_level
        );

void import_julia_data_16(
        bs_t *bs,
        ht_t *ht,
        mat_t *mat,
        const int32_t *const lens,
        void *cfs_julia,
        const int32_t *const exps,
        const int32_t nr_gens
        );

void import_julia_data_32(
        bs_t *bs,
        ht_t *ht,
        mat_t *mat,
        const int32_t *const lens,
        void *cfs_julia,
        const int32_t *const exps,
        const int32_t nr_gens
        );
#endif
