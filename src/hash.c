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

#include "hash.h"

// extern declaration in src/hash.h
mp_cf4_ht_t *ht;

/********************************************************************************
 * FOLLOWING HASH TABLE IMPLEMENTATION IS COPIED FROM COMPACT F4 IMPLEMENTATION 
 * BY MONAGAN AND PIERCE (see PASCO 2015)
 *******************************************************************************/

// global variables used as random seeds, initialized to max unsigned values
// depending on available wordsize of the machine
#if __GB_WORDSIZE==64
uint64_t random_seed  = 0xFFFFFFFFFFFFFFFF;
inline void pseudo_random_generator()
{
  random_seed ^=  (random_seed << 13);
  random_seed ^=  (random_seed << 7);
  random_seed ^=  (random_seed << 17);
}

#elif __GB_WORDSIZE==32
uint32_t random_seed  = 0xFFFFFFFF;
inline void pseudo_random_generator()
{
  random_seed ^=  (random_seed << 13);
  random_seed ^=  (random_seed << 17);
  random_seed ^=  (random_seed << 5);
}

#endif

inline void set_random_seed(mp_cf4_ht_t *hash_table)
{
  hash_t i;

  // use random_seed, no zero values are allowed
  for (i=0; i<hash_table->nvars; ++i) {
    pseudo_random_generator();
    hash_table->rand[i] = random_seed | 1;
  }
}

inline mp_cf4_ht_t *init_hash_table(const ht_size_t hash_table_size,
    const nvars_t number_variables)
{
  hash_t i;

  mp_cf4_ht_t *ht = (mp_cf4_ht_t *)malloc(sizeof(mp_cf4_ht_t));
  
  // global table data
  ht->nvars = number_variables;
  ht->size  = hash_table_size;
  // for easier divisibility checks we start at index 1. If the divisibility
  // check routines return 0, there is no division.
  ht->load  = 1;
  ht->exp   = (exp_t **)malloc(ht->size * sizeof(exp_t *));
  ht->lut   = (hash_t *)calloc(ht->size, sizeof(hash_t));
  ht->val   = (hash_t *)calloc(ht->size, sizeof(hash_t));
  ht->deg   = (deg_t *)calloc(ht->size, sizeof(deg_t));
  ht->div   = (nelts_t *)calloc(ht->size, sizeof(nelts_t));
  ht->idx   = (hash_t *)calloc(ht->size, sizeof(hash_t));
  ht->rand  = (hash_t *)malloc(ht->nvars * sizeof(hash_t));
  // use random_seed, no zero values are allowed
  set_random_seed(ht);
  // get memory for each exponent
  for (i=0; i<ht->size; ++i) {
    ht->exp[i]  = (exp_t *)malloc(ht->nvars * sizeof(exp_t));
  }

  return ht;
}

inline void free_hash_table(mp_cf4_ht_t *hash_table)
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

  free(ht);
  ht  = NULL;
}

inline hash_t get_hash(const exp_t *exp, mp_cf4_ht_t *ht)
{
  hash_t i;
  hash_t hash = 0;

  for (i=0; i<ht->nvars; ++i)
    hash  +=  ht->rand[i] * (hash_t)exp[i];

  return hash;
}

inline hash_t insert_with_quadratic_probing(const exp_t *exp,
    const hash_t hash, mp_cf4_ht_t *ht)
{
  hash_t i, j, tmp;

  tmp = hash;
  for (i=0; i<ht->size; ++i) {
    tmp = (tmp+i) & (ht->size-1);
    if (!ht->lut[tmp])
      break;
    if (ht->val[ht->lut[tmp]] != hash)
      continue;
    for (j=0; j<ht->nvars; ++j) {
      if (ht->exp[ht->lut[tmp]][j] != exp[j])
        break;
    }
    if (j == ht->nvars)
      return ht->lut[tmp];
  }
  return ht->lut[tmp];
}

inline void insert_while_enlarging(const hash_t hash, mp_cf4_ht_t *ht)
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

inline hash_t insert_in_hash_table(const hash_t hash,
    const hash_t pos,  mp_cf4_ht_t *ht)
{
  nvars_t i;
  // the new exponent is already stored in ht->exp[ht->load]
  // deg is possibly not zero due to some earlier check in the hash table
  ht->deg[ht->load] = 0;
  for (i=0; i<ht->nvars; ++i)
    ht->deg[ht->load] +=  ht->exp[ht->load][i];
  // ht->div and ht->idx are already initialized with 0, so nothing to do there
  ht->val[ht->load] = hash;
  ht->lut[pos]      = ht->load;
  ht->load++;

  // we need to keep one place open in ht->exp since the next element to be
  // checked against the hash table will be intermediately stored there
  if (ht->load >= ht->size-1)
    enlarge_hash_table(ht, 2*ht->size);

  return (ht->load-1);
}

inline hash_t insert_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    const hash_t hash, const hash_t pos,  mp_cf4_ht_t *ht)
{
  hash_t i;

  hash_t last_pos = ht->load;

  for (i=0; i<ht->nvars; ++i)
    ht->exp[last_pos][i] = ht->exp[mon_1][i] + ht->exp[mon_2][i];

  // ht->div and ht->idx are already initialized with 0, so nothing to do there
  ht->deg[last_pos] = ht->deg[mon_1] + ht->deg[mon_2];
  ht->val[last_pos] = hash;
  ht->lut[pos]      = last_pos;
  ht->load++;

  if (ht->load >= ht->size)
    enlarge_hash_table(ht, 2*ht->size);

  return last_pos;
}

inline hash_t check_in_hash_table(mp_cf4_ht_t *ht)
{
  hash_t i,j;
  // element to be checked, intermediately stored in the first free position of
  // ht->exp
  exp_t *exp  = ht->exp[ht->load];

  hash_t hash   = get_hash(exp, ht);
  hash_t tmp_h  = hash; // temporary hash values for quadratic probing
  hash_t tmp_l;         // temporary lookup table value

#if HASH_DEBUG
  for (i=0; i<ht->nvars; ++i)
    printf("%u ",exp[i]);
  printf("\nhash = %u\n",hash);
#endif

  for (i=0; i<ht->size; ++i) {
    tmp_h = (tmp_h + i) & (ht->size - 1);
    tmp_l = ht->lut[tmp_h];
    if (tmp_l == 0)
      break;
    if (ht->val[tmp_l] != hash)
      continue;
    for (j=0; j<ht->nvars; ++j)
      if (exp[j] != ht->exp[tmp_l][j])
        break;
    if (j == ht->nvars)
      return tmp_l;
  }
  
  // at this point we know that we do not have the hash value of exp in the
  // table, so we have to insert it
  return insert_in_hash_table(hash, tmp_h, ht);
}

inline hash_t check_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    mp_cf4_ht_t *ht)
{
  hash_t i,j;

  // hash value of the product is the sum of the hash values in our setting
  hash_t hash   = ht->val[mon_1] + ht->val[mon_2];
  hash_t tmp_h  = hash; // temporary hash values for quadratic probing
  hash_t tmp_l;         // temporary lookup table value

  for (i=0; i<ht->size; ++i) {
    tmp_h = (tmp_h + i) & (ht->size - 1);
    tmp_l = ht->lut[tmp_h];
    if (tmp_l == 0)
      break;
    if (ht->val[tmp_l] != hash)
      continue;
    for (j=0; j<ht->nvars; ++j)
      if (ht->exp[tmp_l][j] != ht->exp[mon_1][j] + ht->exp[mon_2][j])
        break;
    if (j == ht->nvars)
      return tmp_l;
  }
  
  // at this point we know that we do not have the hash value of exp in the
  // table, so we have to insert it
  return insert_in_hash_table_product(mon_1, mon_2, hash, tmp_h, ht);
}

inline void enlarge_hash_table(mp_cf4_ht_t *hash_table, const hash_t new_size)
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
    insert_while_enlarging(hash, hash_table);
  }
}

inline hash_t get_lcm(hash_t h1, hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
  exp_t *lcm, *e1, *e2;

  // use first free entry in hash table ht to store possible new lcm monomial
  lcm = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  for (i=0; i<ht->nvars; ++i) {
    lcm[i]  = e1[i] < e2[i] ? e2[i] : e1[i];
  }
  return check_in_hash_table(ht);
}

inline hash_t monomial_division(hash_t h1, hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
  exp_t *e, *e1, *e2;

  e   = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  for (i=0; i<ht->nvars; ++i) {
    if (e1[i] < e2[i])
      return 0;
    e[i]  = e1[i] - e2[i];
  }
  return check_in_hash_table(ht);
}

inline hash_t get_multiplier(hash_t h1, hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
  exp_t *e, *e1, *e2;

  e   = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  // we know that exp e2 divides exp e1, so no check for e1[i] < e2[i]
  for (i=0; i<ht->nvars; ++i)
    e[i]  = e1[i] - e2[i];
  return check_in_hash_table(ht);
}

inline void clear_hash_table_idx(mp_cf4_ht_t *ht)
{
  memset(ht->idx, 0, ht->size * sizeof(hash_t));
}
