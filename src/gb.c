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
 * \file gb.c
 * \brief Overall GB library file
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "stat.h"     /* computational statistics */
#include "time.h"     /* timing business */
#include "tools.h"    /* tools like inversion mod p, etc. */
#include "hash.h"     /* hash table stuff */
#include "order.h"    /* order and comparison procedures */
#include "basis.h"    /* basis and polynomial handling */
#include "la_ff.h"    /* finite field linear algebra */
#include "update.h"   /* update process and pairset handling */
#include "convert.h"  /* conversion between hashes and column indices*/
#include "symbol.h"   /* symbolic preprocessing */
#include "io.h"       /* input and output data handling */
#include "f4.h"       /* implemenation of f4 algorithm */
