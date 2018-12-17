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
 * \file symbol.h
 * \brief Symbolic preprocessing routines
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_SYMBOL_H
#define GB_SYMBOL_H

#include "data.h"
#include "hash.h"
#include "time.h"

mat_t *symbolic_preprocessing(
        mat_t *mat,
        ht_t *ht,
        const bs_t *const bs,
        stat_t *st
        );

mat_t *select_spairs_by_minimal_degree(
        mat_t *mat,
        ps_t *psl,
        ht_t *ht,
        const bs_t *const bs,
        stat_t *st
        );

int find_multiplied_reducer(
        mat_t *mat,
        ht_t *ht,
        hd_t *m,
        const bs_t *const bs
        );
#endif
