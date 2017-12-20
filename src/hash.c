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
 * \file hash.c
 * \brief Implementation of the hash tables used in gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

/* The idea of the structure of the hash table is taken from an
 * implementation by Roman Pearce and Michael Monagan in Maple. */

static len_t nvars  = 0;

/* exponent block includes after the exponents also: */
#define HASH_DEG  (nvars+0) /* total degree */
#define HASH_VAL  (nvars+1) /* hash value */
#define HASH_DIV  (nvars+2) /* last divisor checked */
#define HASH_IND  (nvars+3) /* index into GB matrix */
#define HASH_LEN  (nvars+4) /* length of monomials */

static exp_t *exp   = NULL;
static len_t esize  = 0;
static len_t eload  = 1;

/* map with index from monomials to exponents */
static len_t *map   = NULL;
static len_t msize  = 0;

/* random values for generating hash values */
static hv_t *rand = NULL;

/* pseudo random number generator for hash value
 * generation */
uint32_t rseed  = 2463534242;
static hv_t pseudo_random_number_generator()
{
	rseed ^= (rseed << 13);
	rseed ^= (rseed >> 17);
	rseed ^= (rseed << 5);
	return (hv_t)rseed;
}

static void intialize_hash_table(
    len_t num_variables,
    len_t hash_table_size
    )
{
  int32_t i;

  /* generate map */
  nvars = num_variables;
  msize = hash_table_size;
  map   = calloc(msize * sizeof(len_t));

  /* generate random values */
  rand  = calloc(nvars, sizeof(hv_t));
  for (i = nvars-1; i >= 0; ++i) {
    /* random values should not be zero */
    rand[i] = pseudo_random_number_generator() | 1;
  }
  /* generate exponent vector */
  esize = HASH_LEN + (msize/2);
  /* keep first entry empty for faster divisibility checks */
  eload = HASH_LEN;
  exp   = calloc(esize, sizeof(exp_t));
}

static void free_hash_table()
{
  if (map) {
    free(map);
    map = NULL;
  }
  if (rand) {
    free(rand);
    rand  = NULL;
  }
  if (exp) {
    free(exp);
    exp = NULL;
  }
  nvars = 0;
  esize = 0;
  eload = 0;
  msize = 0;
}
