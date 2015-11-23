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
 * \file hash_table.h
 * \brief Implementation of the hash tables used in gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_HASH_TABLE_H
#define GB_HASH_TABLE_H

#include "gb_config.h"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
#include "types.h"

// global variables used as random seeds, initialized to max unsigned values
// depending on available wordsize of the machine
#if __GB_WORDSIZE==64
uint64_t random_seed  = 0xFFFFFFFFFFFFFFFF;
#elif __GB_WORDSIZE==32
uint32_t random_seed  = 0xFFFFFFFF;
#endif

/**
 * \brief Generates pseudo random numbers using xorshifts using global defined
 * random_seed variable
 *
 * \return some pseudo random number
 */
static inline void pseudo_random_generator()
{
  random_seed ^=  (random_seed << 13);
  random_seed ^=  (random_seed << 7);
  random_seed ^=  (random_seed << 17);
  random_seed ^=  (random_seed << 5);
}

/********************************************************************************
 * FOLLOWING HASH TABLE IMPLEMENTATION IS COPIED FROM COMPACT F4 IMPLEMENTATION 
 * BY MONAGAN AND PIERCE (see PASCO 2015)
 *******************************************************************************/
/**
 * \brief Generates hash table as defined in compact F4 implementation by
 * Monagan and Pearce (see PASCO 2015)
 *
 * \param hash_table_size given size of hash table
 *
 * \param number_variables number of variables in given polynomial ring
 *
 * \return hash table
 */
static inline mp_cf4_ht_t *init_hash_table(const ht_size_t hash_table_size, const nvars_t number_variables)
{
  nvars_t i;

  mp_cf4_ht_t *ht = (mp_cf4_ht_t *)malloc(sizeof(mp_cf4_ht_t));
  
  // global table data
  ht->size  = hash_table_size;
  ht->load  = 0;
  ht->value = (hash_t *)calloc(ht->size, sizeof(hash_t));
  ht->exp   = (hash_t *)calloc(ht->size*number_variables, sizeof(hash_t)); 
  ht->rand  = (hash_t *)malloc(number_variables * sizeof(hash_t));
  ht->deg   = (deg_t *)calloc(ht->size, sizeof(deg_t));
  ht->div   = (nelts_t *)calloc(ht->size, sizeof(nelts_t));
  ht->idx   = (hash_t *)calloc(ht->size, sizeof(hash_t));
  // use random_seed, no zero values are allowed
  for (i=0; i<number_variables; ++i) {
    pseudo_random_generator();
    ht->rand[i] = random_seed | 1;
  }

  return ht;
}
#endif /* GB_HASH_TABLE_H */
