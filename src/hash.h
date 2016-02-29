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
 * \file hash.h
 * \brief Implementation of the hash tables used in gb.
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#ifndef GB_HASH_H
#define GB_HASH_H

#include "gb_config.h"
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#if defined(_MSC_VER)
     /* Microsoft C/C++-compatible compiler */
     #include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
     /* GCC-compatible compiler, targeting x86/x86-64 */
     #include <x86intrin.h>
#elif defined(__GNUC__) && defined(__ARM_NEON__)
     /* GCC-compatible compiler, targeting ARM with NEON */
     #include <arm_neon.h>
#elif defined(__GNUC__) && defined(__IWMMXT__)
     /* GCC-compatible compiler, targeting ARM with WMMX */
     #include <mmintrin.h>
#elif (defined(__GNUC__) || defined(__xlC__)) && (defined(__VEC__) || defined(__ALTIVEC__))
     /* XLC or GCC-compatible compiler, targeting PowerPC with VMX/VSX */
     #include <altivec.h>
#elif defined(__GNUC__) && defined(__SPE__)
     /* GCC-compatible compiler, targeting PowerPC with SPE */
     #include <spe.h>
#endif
#include <omp.h>
#include "types.h"

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#ifndef HASH_DEBUG
#define HASH_DEBUG  0
#endif

#ifndef HASH_QUADRATIC_PROBING
#define HASH_QUADRATIC_PROBING  0
#endif

/***************************
 * OUR HASH TABLE IS GLOBAL
 **************************/
extern mp_cf4_ht_t *ht;

// global variables used as random seeds, initialized to max unsigned values
// depending on available wordsize of the machine
#if __GB_WORDSIZE==64
/**
 * \brief Generates pseudo random numbers using xorshifts using global defined
 * random_seed variable
 *
 * \return some pseudo random number
 */
static inline uint64_t pseudo_random_generator(uint64_t random_seed)
{
  random_seed ^=  (random_seed << 13);
  random_seed ^=  (random_seed << 7);
  random_seed ^=  (random_seed << 17);

  return random_seed;
}

#elif __GB_WORDSIZE==32
/**
 * \brief Generates pseudo random numbers using xorshifts using global defined
 * random_seed variable
 *
 * \return some pseudo random number
 */
static inline uint32_t pseudo_random_generator(uint32_t random_seed)
{
  random_seed ^=  (random_seed << 13);
  random_seed ^=  (random_seed << 17);
  random_seed ^=  (random_seed << 5);

  return random_seed;
}

#endif

/********************************************************************************
 * FOLLOWING HASH TABLE IMPLEMENTATION IS COPIED FROM COMPACT F4 IMPLEMENTATION 
 * BY MONAGAN AND PIERCE (see PASCO 2015)
 *******************************************************************************/
/**
 * \brief Inserts random seeds in hash table random array
 *
 * \param hash table ht
 */
static inline void set_random_seed(mp_cf4_ht_t *ht)
{
  hash_t i;

#if __GB_WORDSIZE==64
uint64_t random_seed  = 0xFFFFFFFFFFFFFFFF;
#elif __GB_WORDSIZE==32
uint32_t random_seed  = 0xFFFFFFFF;
#endif
  // use random_seed, no zero values are allowed
  for (i=0; i<ht->nv; ++i) {
    random_seed = pseudo_random_generator(random_seed);
    ht->rand[i] = random_seed | 1;
  }
}

/**
 * \brief Generates hash table as defined in compact F4 implementation by
 * Monagan and Pearce (see PASCO 2015)
 *
 * \param hash table size index ht_si
 *
 * \param number of variables in given polynomial ring nv
 *
 * \return hash table
 */
static inline mp_cf4_ht_t *init_hash_table(const ht_size_t ht_si,
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
  ht->lut   = (ht_size_t *)calloc(ht->primes[ht->si], sizeof(ht_size_t));
  ht->val   = (hash_t *)calloc(ht->primes[ht->si], sizeof(hash_t));
  ht->deg   = (deg_t *)calloc(ht->primes[ht->si], sizeof(deg_t));
  ht->div   = (nelts_t *)calloc(ht->primes[ht->si], sizeof(nelts_t));
  ht->idx   = (ht_size_t *)calloc(ht->primes[ht->si], sizeof(ht_size_t));
  ht->rand  = (hash_t *)malloc(ht->nv * sizeof(hash_t));
#if __GB_HAVE_SSE2
  ht->nev   = (ht->nv * sizeof(exp_t) / sizeof(exp_v)) + 1;
  ht->vl    = sizeof(exp_v) / sizeof(exp_t);
  ht->ev    = (exp_v **)malloc(ht->primes[ht->si] * sizeof(exp_v *));
  // get memory for each exponent vector
  for (i=0; i<ht->primes[ht->si]; ++i) {
    ht->ev[i]  = (exp_v *)calloc(ht->nev, sizeof(exp_v));
  }
#else
  ht->exp   = (exp_t **)malloc(ht->primes[ht->si] * sizeof(exp_t *));
  // get memory for each exponent
  for (i=0; i<ht->primes[ht->si]; ++i) {
    ht->exp[i]  = (exp_t *)malloc(ht->nv * sizeof(exp_t));
  }
#endif
  // use random_seed, no zero values are allowed
  set_random_seed(ht);

  return ht;
}

/**
 * \brief Inserts in hash table using quadratic probing when enlarging table
 *
 * \param hash value to be inserted hash
 *
 * \param hash table hash_table
 */
static inline void insert_while_enlarging(const hash_t hash, const ht_size_t pos, mp_cf4_ht_t *ht)
{
  ht_size_t i;
  ht_size_t tmp = (ht_size_t)MODP(hash,ht->primes[ht->si]);

  for (i=0; i<ht->primes[ht->si]; ++i) {
#if HASH_QUADRATIC_PROBING
    tmp = MODP(tmp+(i*i),ht->primes[ht->si]);
#else
    tmp = MODP(tmp+(i),ht->primes[ht->si]);
#endif
    if (ht->lut[tmp] != 0) {
      continue;
    } else {
      ht->lut[tmp]  = pos;
    }
    return;
  }
}

/**
 * \brief Enlarges hash table to next prime number size
 *
 * \param hash table ht
 */
static inline void enlarge_hash_table(mp_cf4_ht_t *ht)
{
  ht_size_t i;
  hash_t hash;
  const ht_size_t old_si  = ht->si;
  ht->si++;

#if HASH_DEBUG
  printf("enlarging hash table: load = %u || primes[ %u] = %u --> primes[%u] = %u\n",ht->load,old_si,ht->primes[old_si], ht->si,ht->primes[ht->si]);
#endif
  ht->lut   = realloc(ht->lut, ht->primes[ht->si] * sizeof(ht_size_t));
  ht->val   = realloc(ht->val, ht->primes[ht->si] * sizeof(hash_t));
  ht->deg   = realloc(ht->deg, ht->primes[ht->si] * sizeof(deg_t));
  ht->idx   = realloc(ht->idx, ht->primes[ht->si] * sizeof(ht_size_t));
  ht->div   = realloc(ht->div, ht->primes[ht->si] * sizeof(nelts_t));
#if __GB_HAVE_SSE2
  ht->ev    = realloc(ht->ev, ht->primes[ht->si] * sizeof(exp_v *));
  for (i=ht->primes[old_si]; i<ht->primes[ht->si]; ++i) {
    ht->ev[i] = (exp_v *)calloc(ht->nev, sizeof(exp_v));
  }
#else
  ht->exp   = realloc(ht->exp, ht->primes[ht->si] * sizeof(exp_t *));
  for (i=ht->primes[old_si]; i<ht->primes[ht->si]; ++i) {
    ht->exp[i]  = (exp_t *)calloc(ht->nv, sizeof(exp_t));
  }
#endif
  // re-insert all elements in block
  memset(ht->lut, 0, ht->primes[ht->si] * sizeof(ht_size_t));
  for (i=0; i<ht->load; ++i) {
    hash  = ht->val[i];
    insert_while_enlarging(hash, i, ht);
  }
}

/**
 * \brief Frees dynamically allocated data in hash table
 *
 * \param hash table ht
 */
static inline void free_hash_table(mp_cf4_ht_t *ht)
{
  if (ht) {

    hash_t i;

    free(ht->lut);
    free(ht->val);
    free(ht->rand);
    free(ht->deg);
    free(ht->div);
    free(ht->idx);
#if __GB_HAVE_SSE2
    for (i=0; i<ht->primes[ht->si]; ++i)
      free(ht->ev[i]);
    free(ht->ev);
#else
    for (i=0; i<ht->primes[ht->si]; ++i)
      free(ht->exp[i]);
    free(ht->exp);
#endif
    free(ht->primes);
  }

  free(ht);
  ht  = NULL;
}

/**
 * \brief Get hash value
 *
 * \param exponent vector exp
 *
 * \param hash table hash_table
 *
 * \return hash value
 */
#if __GB_HAVE_SSE2
static inline hash_t get_hash(const exp_v *ev, mp_cf4_ht_t *ht)
{
  nvars_t i;
  hash_t hash = 0;

  exp_t *exp  = (exp_t *)calloc(ht->nev + ht->vl, sizeof(exp_t));

  for (i=0; i<ht->nev; ++i)
    _mm_storeu_si128((exp_v *)exp+(i*ht->vl), ev[i]);

  for (i=0; i<ht->nv; ++i)
    hash  +=  ht->rand[i] * exp[i];

  free(exp);
  return hash;
}
#else
static inline hash_t get_hash(const exp_t *exp, mp_cf4_ht_t *ht)
{
  hash_t i;
  hash_t hash = 0;

  for (i=0; i<ht->nv; ++i)
    hash  +=  ht->rand[i] * exp[i];

  return hash;
}
#endif

/**
 * \brief Inserts a new element to the hash table
 *
 * \note The exponent to be checked is already stored in ht->exp[ht->load]. If
 * it has to be inserted, we have it already at the correct place. Otherwise,
 * ht->load will not be increased and ht->exp[ht->load] will be overwritten
 * later on.
 *
 * \param hash value of exp hash
 *
 * \param position in lookup table pos
 *
 * \param hash table ht
 *
 * \return position of hash of exp in table
 */
static inline hash_t insert_in_hash_table(const hash_t hash,
    const ht_size_t pos,  mp_cf4_ht_t *ht)
{
  // the new exponent is already stored in ht->exp[ht->load] and also
  // ht->deg[ht->load] is already set

  // ht->div and ht->idx are already initialized with 0, so nothing to do there
  ht->val[ht->load] = hash;
  ht->lut[pos]      = ht->load;
  ht->load++;

  // we need to keep one place open in ht->exp since the next element to be
  // checked against the hash table will be intermediately stored there
#if HASH_QUADRATIC_PROBING
  if (ht->load >= ht->primes[ht->si]/2-1)
#else
  if (ht->load >= ht->primes[ht->si]-1)
#endif
    enlarge_hash_table(ht);

  return (ht->load-1);
}

/**
 * \brief Inserts a new element to the hash table coming from a product of two
 * monomials.
 *
 * \note We use the sum of the hash values of mon_1 and mon_2 as hash value for
 * the product of mon_1 and mon_2.
 *
 * \note The exponent to be checked is already stored in ht->exp[ht->load]. If
 * it has to be inserted, we have it already at the correct place. Otherwise,
 * ht->load will not be increased and ht->exp[ht->load] will be overwritten
 * later on.
 *
 * \param monomial 1 mon_1
 *
 * \param monomial 2 mon_2
 *
 * \param hash value of exp hash
 *
 * \param position in lookup table pos
 *
 * \param hash table ht
 *
 * \return position of hash of exp in table
 */
static inline hash_t insert_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    const hash_t hash, const ht_size_t pos,  mp_cf4_ht_t *ht)
{

  ht_size_t last_pos = ht->load;

#if !__GB_HAVE_SSE2
  nvars_t i;
  for (i=0; i<ht->nv; ++i)
    ht->exp[last_pos][i] = ht->exp[mon_1][i] + ht->exp[mon_2][i];
#endif
  // ht->div and ht->idx are already initialized with 0, so nothing to do there
  ht->deg[last_pos] = ht->deg[mon_1] + ht->deg[mon_2];
  ht->val[last_pos] = hash;
  ht->lut[pos]      = last_pos;

  // we do not need this anymore since it is already computed and stored in
  // check_in_hash_table_product()
  /*
#if __GB_HAVE_SSE2
  ht->ev[last_pos]  = _mm_adds_epu8(ht->ev[mon_1], ht->ev[mon_2]);
#endif
  */
  ht->load++;

#if HASH_QUADRATIC_PROBING
  if (ht->load >= ht->primes[ht->si]/2-1)
#else
  if (ht->load >= ht->primes[ht->si]-1)
#endif
    enlarge_hash_table(ht);

  return last_pos;
}

/**
 * \brief Checks if the given monomial exponent is already in the hash table. If
 * not, it is added to the table
 *
 * \note The exponent to be checked is already stored in ht->exp[ht->load]. If
 * it has to be inserted, we have it already at the correct place. Otherwise,
 * ht->load will not be increased and ht->exp[ht->load] will be overwritten
 * later on.
 *
 * \param hash table ht
 *
 * \return position of hash of exp in table
 */
static inline hash_t check_in_hash_table(mp_cf4_ht_t *ht)
{
  nvars_t i, j;
  // element to be checked, intermediately stored in the first free position of
  // ht->exp

#if __GB_HAVE_SSE2
  hash_t hash = get_hash(ht->ev[ht->load], ht);
#else
  exp_t *exp  = ht->exp[ht->load];
  hash_t hash = get_hash(exp, ht);
#endif
  ht_size_t tmp_h = (ht_size_t)MODP(hash, ht->primes[ht->si]); // temporary hash values for linear probing
  ht_size_t tmp_l;         // temporary lookup table value

  // first check directly
  tmp_l = ht->lut[tmp_h];
  if (tmp_l == 0)
    return insert_in_hash_table(hash, tmp_h, ht);
  if (ht->val[tmp_l] == hash) {
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi8(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    return tmp_l;
#else
    for (j=0; j<ht->nv; ++j)
      if (exp[j] != ht->exp[tmp_l][j])
        break;
    if (j == ht->nv) {
      return tmp_l;
    }
#endif
  }
#if HASH_DEBUG
  for (i=0; i<ht->nv; ++i)
    printf("%u ",exp[i]);
  printf("\nhash = %u\n",hash);
#endif

  // remaining checks with probing
  for (i=1; i<ht->primes[ht->si]; ++i) {
#if HASH_QUADRATIC_PROBING
    tmp_h = MODP(tmp_h+(i*i),ht->primes[ht->si]);
#else
    tmp_h = MODP(tmp_h+(i),ht->primes[ht->si]);
#endif
    tmp_l = ht->lut[tmp_h];
    if (tmp_l == 0)
      break;
    if (ht->val[tmp_l] != hash)
      continue;
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi8(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    return tmp_l;
#else
    nvars_t j;
    for (j=0; j<ht->nv; ++j)
      if (exp[j] != ht->exp[tmp_l][j])
        break;
    if (j == ht->nv) {
      return tmp_l;
    }
#endif
  }
  // at this point we know that we do not have the hash value of exp in the
  // table, so we have to insert it
  return insert_in_hash_table(hash, tmp_h, ht);
}

/**
 * \brief Finds the product of the given two monomial exponents is already
 * in the hash table.
 *
 * \note This function is used when constructing the gbla matrix, it is assumed
 * that the product is in the hash table, thus the hash table is used const.
 *
 * \param monomial 1 mon_1
 *
 * \param monomial 2 mon_2
 *
 * \param hash table ht
 *
 * \return position of hash of exp in table
 */
static inline hash_t find_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    const mp_cf4_ht_t *ht)
{
  ht_size_t i, j;

  // hash value of the product is the sum of the hash values in our setting
  hash_t hash     = ht->val[mon_1] + ht->val[mon_2];
  ht_size_t tmp_h = (ht_size_t)MODP(hash, ht->primes[ht->si]); // temporary hash values for linear probing
  ht_size_t tmp_l;         // temporary lookup table value
  
  // first check directly
  tmp_l = ht->lut[tmp_h];
#if __GB_HAVE_SSE2
  for (i=0; i<ht->nev; ++i)
    ht->ev[ht->load][i] = _mm_adds_epu8(ht->ev[mon_1][i], ht->ev[mon_2][i]);
#endif
  if (tmp_l == 0)
    return 0;
  if (ht->val[tmp_l] == hash) {
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi8(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    return tmp_l;
#else
    for (j=0; j<ht->nv; ++j)
      if (ht->exp[tmp_l][j] != ht->exp[mon_1][j] + ht->exp[mon_2][j])
        break;
    if (j == ht->nv)
      return tmp_l;
#endif
  }

  // remaining checks with probing
  for (i=1; i<ht->primes[ht->si]; ++i) {
#if HASH_QUADRATIC_PROBING
    tmp_h = MODP(tmp_h+(i*i),ht->primes[ht->si]);
#else
    tmp_h = MODP(tmp_h+(i),ht->primes[ht->si]);
#endif
    tmp_l = ht->lut[tmp_h];
    if (tmp_l == 0)
      break;
    if (ht->val[tmp_l] != hash)
      continue;
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi8(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    return tmp_l;
#else
    for (j=0; j<ht->nv; ++j)
      if (ht->exp[tmp_l][j] != ht->exp[mon_1][j] + ht->exp[mon_2][j])
        break;
    if (j == ht->nv)
      return tmp_l;
#endif
  }
  return 0;
}

/**
 * \brief Checks if the product of the given two monomial exponents is already
 * in the hash table. If not, it is added to the table
 *
 * \note We use the sum of the hash values of mon_1 and mon_2 as hash value for
 * the product of mon_1 and mon_2.
 *
 * \note The exponent to be checked is already stored in ht->exp[ht->load]. If
 * it has to be inserted, we have it already at the correct place. Otherwise,
 * ht->load will not be increased and ht->exp[ht->load] will be overwritten
 * later on.
 *
 * \param monomial 1 mon_1
 *
 * \param monomial 2 mon_2
 *
 * \param hash table ht
 *
 * \return position of hash of exp in table
 */
static inline hash_t check_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    mp_cf4_ht_t *ht)
{
  ht_size_t i, j;

  // hash value of the product is the sum of the hash values in our setting
  hash_t hash   = ht->val[mon_1] + ht->val[mon_2];
  ht_size_t tmp_h = (ht_size_t)MODP(hash, ht->primes[ht->si]); // temporary hash values for linear probing
  ht_size_t tmp_l;         // temporary lookup table value

  // first check directly
  tmp_l = ht->lut[tmp_h];
#if __GB_HAVE_SSE2
  for (i=0; i<ht->nev; ++i)
    ht->ev[ht->load][i] = _mm_adds_epu8(ht->ev[mon_1][i], ht->ev[mon_2][i]);
#endif
  if (tmp_l == 0)
    return insert_in_hash_table_product(mon_1, mon_2, hash, tmp_h, ht);
  if (ht->val[tmp_l] == hash) {
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi8(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    return tmp_l;
#else
    for (j=0; j<ht->nv; ++j)
      if (ht->exp[tmp_l][j] != ht->exp[mon_1][j] + ht->exp[mon_2][j])
        break;
    if (j == ht->nv)
      return tmp_l;
#endif
  }

  // remaining checks with probing
  for (i=1; i<ht->primes[ht->si]; ++i) {
#if HASH_QUADRATIC_PROBING
    tmp_h = MODP(tmp_h+(i*i),ht->primes[ht->si]);
#else
    tmp_h = MODP(tmp_h+(i),ht->primes[ht->si]);
#endif
    tmp_l = ht->lut[tmp_h];
    if (tmp_l == 0)
      break;
    if (ht->val[tmp_l] != hash)
      continue;
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi8(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    return tmp_l;
#else
    nvars_t j;
    for (j=0; j<ht->nv; ++j)
      if (ht->exp[tmp_l][j] != ht->exp[mon_1][j] + ht->exp[mon_2][j])
        break;
    if (j == ht->nv)
      return tmp_l;
#endif
  }
  // at this point we know that we do not have the hash value of exp in the
  // table, so we have to insert it
  return insert_in_hash_table_product(mon_1, mon_2, hash, tmp_h, ht);
}

/**
 * \brief Returns position of lcm of the exponents stored at position h1 and h2
 * in the hash table.
 *
 * \param position of first generator monomial h1
 *
 * \param position of second generator monomial h2
 *
 * \param hash table ht
 *
 * \return position of lcm of generators h1 and h2 in hash table ht
 */
static inline hash_t get_lcm(hash_t h1, hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
#if __GB_HAVE_SSE2
  exp_t exp[ht->nev*ht->vl];
  for (i=0; i<ht->nev; ++i) {
    ht->ev[ht->load][i] = _mm_max_epu8(ht->ev[h1][i], ht->ev[h2][i]);
    _mm_storeu_si128((exp_v *)exp + i*ht->vl, ht->ev[ht->load][i]);
  }
  ht->deg[ht->load] = 0;
  for (i=0; i<ht->nv; ++i)
    ht->deg[ht->load] += exp[i];
#else
  exp_t *lcm, *e1, *e2;
  deg_t deg = 0;

  // use first free entry in hash table ht to store possible new lcm monomial
  lcm = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  for (i=0; i<ht->nv; ++i) {
    deg +=  lcm[i]  = e1[i] < e2[i] ? e2[i] : e1[i];
  }
  ht->deg[ht->load] = deg;
#endif
  return check_in_hash_table(ht);
}

/**
 * \brief Tests if exp of h1 is divisible by exp of h2. If divisibility is
 * fulfilled we add the multiplier to the hash table and return its hash
 * position. Else we return 0.
 *
 * \param hash position h1
 *
 * \param hash position h2
 *
 * \param hash table ht
 *
 * \return hash position of multiplier or 0
 */
static inline hash_t monomial_division(hash_t h1, hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
#if __GB_HAVE_SSE2
  exp_v cmpv;
  for (i=0; i<ht->nev; ++i) {
    cmpv  = _mm_cmplt_epi8(ht->ev[h1][i], ht->ev[h2][i]);
    if (_mm_movemask_epi8(cmpv) != 0)
      return 0;
  }
  for (i=0; i<ht->nev; ++i)
    ht->ev[ht->load][i] = _mm_subs_epi8(ht->ev[h1][i], ht->ev[h2][i]);
#else
  exp_t *e, *e1, *e2;

  e   = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  for (i=0; i<ht->nv; ++i) {
    if (e1[i] < e2[i])
      return 0;
    e[i]  = e1[i] - e2[i];
  }
#endif
  ht->deg[ht->load] = ht->deg[h1] - ht->deg[h2];
  return check_in_hash_table(ht);
}

/**
 * \brief Tests if exp of h1 is divisible by exp of h2. If divisibility is
 * fulfilled only 1 is returned, else 0.
 *
 * \note This procedure only tests "if" divisible, but not "by which" it is
 * divisible. This is enough for detecting redundant elements when new elements
 * are added to the intermediate groebner basis.
 *
 * \param hash position h1
 *
 * \param hash position h2
 *
 * \param hash table ht
 *
 * \return 0 if not divisible, 1 is divisible
 */
static inline int check_monomial_division(hash_t h1, hash_t h2, const mp_cf4_ht_t *ht)
{
  nvars_t i;
#if __GB_HAVE_SSE2
  exp_v cmpv;
  for (i=0; i<ht->nev; ++i) {
    cmpv  = _mm_cmplt_epi8(ht->ev[h1][i], ht->ev[h2][i]);
    if (_mm_movemask_epi8(cmpv) != 0)
      return 0;
  }
  return 1;
#else
  exp_t *e1, *e2;

  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  for (i=0; i<ht->nv; ++i) {
    if (e1[i] < e2[i])
      return 0;
  }
  return 1;
#endif
}

/**
 * \brief Returns the multiplier needed to multiply h2 with in order to get
 * h1
 *
 * \note It is nearly the same function as monomial_division, but here we know
 * that h2 divides h1, so we do not have to check this.
 *
 * \param hash position h1
 *
 * \param hash position h2
 *
 * \param hash table ht
 *
 * \return hash position of multiplier 
 */
static inline hash_t get_multiplier(hash_t h1, hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
#if __GB_HAVE_SSE2
  for (i=0; i<ht->nev; ++i)
    ht->ev[ht->load][i] = _mm_subs_epi8(ht->ev[h1][i], ht->ev[h2][i]);
#else
  exp_t *e, *e1, *e2;

  e   = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  // we know that exp e2 divides exp e1, so no check for e1[i] < e2[i]
  for (i=0; i<ht->nv; ++i)
    e[i]  = e1[i] - e2[i];
#endif
  ht->deg[ht->load] = ht->deg[h1] - ht->deg[h2];
  return check_in_hash_table(ht);
}

/**
 * \brief Resets all idx entries of the hash table to zero
 *
 * \param hash table ht
 */
static inline void clear_hash_table_idx(mp_cf4_ht_t *ht)
{
  memset(ht->idx, 0, ht->primes[ht->si] * sizeof(ht_size_t));
}
#endif /* GB_HASH_TABLE_H */
