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
 * \file data.c
 * \brief General and global data
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include "data.h";

static void initialize()
{
  bload = 0;
  bsize = 64;

  pload = 0;
  psize = 192;

  gb  = malloc(bsize * sizeof(poly_t));
  ps  = malloc(psize * sizeof(spair_t));

  /* statistics */
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
