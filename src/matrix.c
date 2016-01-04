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
 * \file matrix.c
 * \brief Implementation of the construction and conversion from and to groebner
 * basis matrices.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "matrix.h"
/*
mat_t *initialize_gbla_matrix(spd_t *spd)
{
  mat_t mat = (mat_t *)malloc(sizeof(mat_t));

  sb_fl_t *A  = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  dbm_fl_t *B = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
  sb_fl_t *C  = (sb_fl_t *)malloc(sizeof(sb_fl_t));
  dbm_fl_t *D = (dbm_fl_t *)malloc(sizeof(dbm_fl_t));
}
*/
