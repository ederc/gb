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
 * \file statistics.c
 * \brief Global data for covering computational statistics
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include <stdint.h>

static double select_time = 0;
static double symbol_time = 0;
static double matrix_time = 0;
static double update_time = 0;
static double polmat_time = 0;
static double reduce_time = 0;

static double select_cput = 0;
static double symbol_cput = 0;
static double matrix_cput = 0;
static double update_cput = 0;
static double polmat_cput = 0;
static double reduce_cput = 0;

static int64_t num_reduced = 0;
static int64_t num_useless = 0;
static int64_t num_rowsred = 0;
static int64_t num_zerored = 0;

static void initialize_statistics(
    void
    )
{
	num_reduced = 0;
  num_useless = 0;
  num_rowsred = 0;
  num_zerored = 0;

	select_time = 0;
  symbol_time = 0;
  matrix_time = 0;
  update_time = 0;
  polmat_time = 0;
  reduce_time = 0;

	select_cput = 0;
  symbol_cput = 0;
  matrix_cput = 0;
  update_cput = 0;
  polmat_cput = 0;
  reduce_cput = 0;
}
