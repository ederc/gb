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
#include "io.h"
#include <src/poly.h>
#include <src/spair.h>
#include <src/symbol.h>
#include <src/matrix.h>

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

/**
 * \brief Prints help for gb call.
 */
void print_help();


#endif
