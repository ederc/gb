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
 * \brief Inserts random seeds in hash table random array
 *
 * \param hash table hash_table
 */
static inline void set_random_seed(mp_cf4_ht_t *hash_table)
{
  nvars_t i;

  // use random_seed, no zero values are allowed
  for (i=0; i<hash_table->nvars; ++i) {
    pseudo_random_generator();
    hash_table->rand[i] = random_seed | 1;
  }
}

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
  ht->nvars = number_variables;
  ht->size  = hash_table_size;
  ht->load  = 0;
  ht->lut   = (hash_t *)calloc(ht->size, sizeof(hash_t));
  ht->val   = (hash_t *)calloc(ht->size, sizeof(hash_t));
  ht->exp   = (exp_t **)malloc(ht->size * sizeof(exp_t *));
  ht->deg   = (deg_t *)calloc(ht->size, sizeof(deg_t));
  ht->div   = (nelts_t *)calloc(ht->size, sizeof(nelts_t));
  ht->idx   = (hash_t *)calloc(ht->size, sizeof(hash_t));
  ht->rand  = (hash_t *)malloc(ht->nvars * sizeof(hash_t));
  // use random_seed, no zero values are allowed
  set_random_seed(ht);
  // get memory for each exponent
  for (i=0; i<ht->size; ++i) {
    ht->exp[i]  = (exp_t *)calloc(ht->nvars, sizeof(exp_t));
  }

  return ht;
}

/**
 * \brief Frees dynamically allocated data in hash table
 *
 * \param hash table hash_table
 */
static inline void free_hash_table(mp_cf4_ht_t *hash_table)
{
  if (hash_table) {

    hash_t i;

    free(hash_table->lut);
    free(hash_table->val);
    free(hash_table->rand);
    free(hash_table->deg);
    free(hash_table->div);
    free(hash_table->idx);
    for (i=0; i<hash_table->size; ++i)
      free(hash_table->exp[i]);
    free(hash_table->exp);
  }
}

/**
 * \brief Inserts in hash table using quadratic probing
 *
 * \param hash value to be inserted hash
 *
 * \param hash table hash_table
 */
static inline void insert_with_quadratic_probing(hash_t hash, mp_cf4_ht_t *ht)
{
  hash_t i, tmp;

  tmp = hash;
  for (i=0; i<ht->size; ++i) {
    tmp = (tmp+i) & (ht->size-1);
    if (ht->lut[tmp])
      continue;
    ht->lut[tmp]  = i;
    break;
  }
}

/**
 * \brief Frees dynamically allocated data in hash table
 *
 * \param hash table hash_table
 *
 * \param new size of hash table new_size
 */
static inline void enlarge_hash_table(mp_cf4_ht_t *hash_table, hash_t new_size)
{
  hash_t i, hash;

  hash_table->lut   = realloc(hash_table->lut, new_size);
  hash_table->val   = realloc(hash_table->val, new_size);
  hash_table->deg   = realloc(hash_table->deg, new_size);
  hash_table->idx   = realloc(hash_table->idx, new_size);
  hash_table->div   = realloc(hash_table->div, new_size);
  hash_table->exp   = realloc(hash_table->exp, new_size);
  for (i=hash_table->size; i<new_size; ++i) {
    hash_table->exp[i]  = (exp_t *)calloc(hash_table->nvars, sizeof(exp_t));
  }
  hash_table->size  = new_size;
  // re-insert all elements in block
  memset(hash_table->lut, 0, new_size * sizeof(hash_t));
  for (i=0; i<hash_table->load; ++i) {
    hash  = hash_table->val[i];

    // insert using quadratic probing
    insert_with_quadratic_probing(hash, hash_table);
  }
}
#endif /* GB_HASH_TABLE_H */
