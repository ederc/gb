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

/****************************************************************************************
 * AVX/SSE stuff to try
 *
 * _mm256_add_epi8
 * _mm_add_epi8
 * _mm_adds_epu8 (unsigned with saturation)
 * _mm256_adds_epu8
 * _mm_subs_epu8 (unsigned with saturation)
 * _mm256_subs_epu8
 * _mm_cmplt_epi8 (-1 if true, 0 if false)
 * _mm256_cmpgt_epi8 (-1 if true, 0 if false)
 ***************************************************************************************/
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

void set_random_seed(mp_cf4_ht_t *ht)
{
  hash_t i;

  // use random_seed, no zero values are allowed
  for (i=0; i<ht->nv; ++i) {
    pseudo_random_generator();
    ht->rand[i] = random_seed | 1;
  }
}

mp_cf4_ht_t *init_hash_table(const ht_size_t ht_si,
    const nvars_t nv)
{
  hash_t i;

  mp_cf4_ht_t *ht = (mp_cf4_ht_t *)malloc(sizeof(mp_cf4_ht_t));
  
  // global table data

  // sets possible non-Mersenne primes for hash table size
  // we assume minimal <2^18 and a maximal of <2^32 as feasible range
  ht->primes  = (ht_size_t *)malloc(15 * sizeof(ht_size_t));
  ht->primes[0]   = 262139;     // < 2^18
  ht->primes[1]   = 524269;     // < 2^19
  ht->primes[2]   = 1048573;    // < 2^20
  ht->primes[3]   = 2097143;    // < 2^21
  ht->primes[4]   = 4194301;    // < 2^22
  ht->primes[5]   = 8388593;    // < 2^23
  ht->primes[6]   = 16777213;   // < 2^24
  ht->primes[7]   = 33554393;   // < 2^25
  ht->primes[8]   = 67108859;   // < 2^26
  ht->primes[9]   = 134217689;  // < 2^27
  ht->primes[10]  = 268435399;  // < 2^28
  ht->primes[11]  = 536870909;  // < 2^29
  ht->primes[12]  = 1073741789; // < 2^30
  ht->primes[13]  = 2147483629; // < 2^31
  ht->primes[14]  = 4294967295; // < 2^32

  ht->nv    = nv;
  // sets size index in primes array for hash table size
  ht->si    = ht_si;
  // for easier divisibility checks we start at index 1. If the divisibility
  // check routines return 0, there is no division.
  ht->load  = 1;
  ht->exp   = (exp_t **)malloc(ht->primes[ht->si] * sizeof(exp_t *));
  ht->lut   = (hash_t *)calloc(ht->primes[ht->si], sizeof(hash_t));
  ht->val   = (hash_t *)calloc(ht->primes[ht->si], sizeof(hash_t));
  ht->deg   = (deg_t *)calloc(ht->primes[ht->si], sizeof(deg_t));
  ht->div   = (nelts_t *)calloc(ht->primes[ht->si], sizeof(nelts_t));
  ht->idx   = (hash_t *)calloc(ht->primes[ht->si], sizeof(hash_t));
  ht->rand  = (hash_t *)malloc(ht->nv * sizeof(hash_t));
  // use random_seed, no zero values are allowed
  set_random_seed(ht);
  // get memory for each exponent
  for (i=0; i<ht->primes[ht->si]; ++i) {
    ht->exp[i]  = (exp_t *)malloc(ht->nv * sizeof(exp_t));
  }
#if HAVE_SSE2
  ht->ev    = (exp_v *)malloc(ht->primes[ht->si] * sizeof(exp_v));
#endif

  return ht;
}

inline void enlarge_hash_table(mp_cf4_ht_t *ht)
{
  hash_t i, hash;
  const ht_size_t old_si  = ht->si;
  ht->si++;

  ht->lut   = realloc(ht->lut, ht->primes[ht->si] * sizeof(hash_t));
  ht->val   = realloc(ht->val, ht->primes[ht->si] * sizeof(hash_t));
  ht->deg   = realloc(ht->deg, ht->primes[ht->si] * sizeof(deg_t));
  ht->idx   = realloc(ht->idx, ht->primes[ht->si] * sizeof(hash_t));
  ht->div   = realloc(ht->div, ht->primes[ht->si] * sizeof(nelts_t));
  ht->exp   = realloc(ht->exp, ht->primes[ht->si] * sizeof(exp_t *));
#if HAVE_SSE2
  ht->ev    = realloc(ht->ev, ht->primes[ht->si] * sizeof(exp_v));
#endif
  for (i=old_si; i<ht->primes[ht->si]; ++i) {
    ht->exp[i]  = (exp_t *)calloc(ht->nv, sizeof(exp_t));
  }
  // re-insert all elements in block
  memset(ht->lut, 0, ht->primes[ht->si] * sizeof(hash_t));
  for (i=0; i<ht->load; ++i) {
    hash  = ht->val[i];

    // insert using linear probing
    insert_while_enlarging(hash, ht);
  }
}
