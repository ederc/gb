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
 * \file la_ff.h
 * \brief Implementation of the linear algebra parts over finite fields.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_LA_FF_H
#define GB_LA_FF_H

#include "data.h"
#include "hash.h"
#include "time.h"

inline mat_t *initialize_matrix(
        stat_t *st
        )
{
    mat_t *mat  = (mat_t *)calloc(1, sizeof(mat_t));

    return mat;
}

inline void free_matrix(
        mat_t **matp
        )
{
    mat_t *mat  = *matp;
    free(mat);
    *matp = mat;
}

#endif
