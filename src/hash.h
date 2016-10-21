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
  /*
  random_seed ^=  (random_seed << 3);
  random_seed ^=  (random_seed << 11);
  random_seed ^=  (random_seed << 7);
  */
#if 1
  random_seed = (uint64_t) ((1103515245 * ((uint64_t)random_seed) + 12345) & 0x7fffffffUL);
  return random_seed;
#else
  uint64_t r = random_seed;
  for (int i=0; i<64; i += 30) {
    r = r*(RAND_MAX + (uint64_t)1) + rand();
  }
  return r;
#endif
}

/*
static inline uint32_t pseudo_random_generator(uint32_t random_seed)
{
  random_seed = (uint32_t) ((1103515245 * ((unsigned int)random_seed) + 12345) & 0x7fffffffUL);
  return random_seed;
  random_seed = 36969*(random_seed & 65535) + (random_seed >> 16);
  random_seed = 18000*(random_seed & 65535) ^ (random_seed >> 16);
  return (random_seed << 16) + random_seed;
}
*/

#elif __GB_WORDSIZE==32
/**
 * \brief Generates pseudo random numbers using xorshifts using global defined
 * random_seed variable
 *
 * \return some pseudo random number
 */
static inline uint32_t pseudo_random_generator(uint32_t random_seed)
{
  /*
  random_seed ^=  (random_seed << 11);
  random_seed ^=  (random_seed << 13);
  random_seed ^=  (random_seed << 18);

  return random_seed;
  */
  random_seed = 36969*(random_seed & 65535) * (random_seed >> 16);
  random_seed = 18000*(random_seed & 65535) ^ (random_seed >> 16);
  return (random_seed << 16) * random_seed;
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
//uint64_t random_seed  = 88172645463325252LL;
//uint64_t random_seed  = 2463534242;
//uint64_t random_seed  = 0xF1FF3FAFFFCFF6F7;
uint32_t random_seed  = 2463534243;
#elif __GB_WORDSIZE==32
uint32_t random_seed  = 2463534242;
//uint32_t random_seed  = 0xFFFFFFFF;
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
  ht->sz    = pow(2,ht_si);
  // we add one extra variable in case we have to homogenize the system
  ht->nv    = nv+1;
  // for easier divisibility checks we start at index 1. If the divisibility
  // check routines return 0, there is no division.
  ht->load  = 1;
  ht->lut   = (ht_size_t *)calloc(ht->sz, sizeof(ht_size_t));
  ht->val   = (hash_t *)calloc(ht->sz, sizeof(hash_t));
  ht->deg   = (deg_t *)calloc(ht->sz, sizeof(deg_t));
  ht->div   = (nelts_t *)calloc(ht->sz, sizeof(nelts_t));
  ht->idx   = (ht_size_t *)calloc(ht->sz, sizeof(ht_size_t));
  ht->rand  = (hash_t *)malloc(ht->nv * sizeof(hash_t));
#if __GB_HAVE_SSE2
  ht->nev   = (ht->nv * sizeof(exp_t) / sizeof(exp_v)) + ((ht->nv * sizeof(exp_t)) % sizeof(exp_v) > 0);
  ht->vl    = sizeof(exp_v) / sizeof(exp_t);
  ht->ev    = (exp_v **)malloc(ht->sz * sizeof(exp_v *));
  // get memory for each exponent vector
  for (i=0; i<ht->sz; ++i) {
    ht->ev[i]  = (exp_v *)calloc(ht->nev, sizeof(exp_v));
  }
#else
  ht->exp   = (exp_t **)malloc(ht->sz * sizeof(exp_t *));
  // get memory for each exponent
  for (i=0; i<ht->sz; ++i) {
    ht->exp[i]  = (exp_t *)calloc(ht->nv, sizeof(exp_t));
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
  ht_size_t tmp_h;
  tmp_h = hash & (ht->sz-1);
  //printf("hash %u --> %u\n", hash, tmp_h);

  for (i=0; i<ht->sz; ++i) {
    tmp_h = (tmp_h+i) & (ht->sz-1);
    if (ht->lut[tmp_h] != 0) {
      continue;
    } else {
      ht->lut[tmp_h]  = pos;
#if HASH_DEBUG
      for (int i=0; i<ht->nv; ++i)
        printf("%u ",ht->exp[pos][i]);
      printf(" ||| ");
      printf("%11u | %11u | %5u\n",hash, pos, ht->deg[pos]);
#endif
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

  const ht_size_t old_sz  = ht->sz;
  ht->sz  = 2*ht->sz;
//#if HASH_DEBUG
  printf("enlarging hash table: %10u --> %10u\n", old_sz, ht->sz);
//#endif
  ht->lut   = realloc(ht->lut, ht->sz * sizeof(ht_size_t));
  ht->val   = realloc(ht->val, ht->sz * sizeof(hash_t));
  ht->deg   = realloc(ht->deg, ht->sz * sizeof(deg_t));
  ht->idx   = realloc(ht->idx, ht->sz * sizeof(ht_size_t));
  ht->div   = realloc(ht->div, ht->sz * sizeof(nelts_t));
  // set mew values for divisors and index to zero
  memset(ht->idx+old_sz, 0, (ht->sz-old_sz) * sizeof(ht_size_t));
  memset(ht->lut+old_sz, 0, (ht->sz-old_sz) * sizeof(ht_size_t));
  memset(ht->div+old_sz, 0, (ht->sz-old_sz) * sizeof(nelts_t));
#if __GB_HAVE_SSE2
  ht->ev    = realloc(ht->ev, ht->sz * sizeof(exp_v *));
  for (i=old_sz; i<ht->sz; ++i) {
    ht->ev[i] = (exp_v *)calloc(ht->nev, sizeof(exp_v));
  }
#else
  ht->exp   = realloc(ht->exp, ht->sz * sizeof(exp_t *));
  for (i=old_sz; i<ht->sz; ++i) {
    ht->exp[i]  = (exp_t *)calloc(ht->nv, sizeof(exp_t));
  }
#endif
  // re-insert all elements in block
  memset(ht->lut+1, 0, (ht->sz-1) * sizeof(ht_size_t));
  for (i=1; i<ht->load; ++i) {
    hash  = ht->val[i];
    //printf("coming from position %u ---> ",i);
    insert_while_enlarging(hash, i, ht);
  }
}

/**
 * \brief Frees dynamically allocated data in hash table
 *
 * \param hash table ht
 */
static inline void free_hash_table(mp_cf4_ht_t **ht_in)
{
  mp_cf4_ht_t *ht = *ht_in;
  if (ht) {

    hash_t i;

    free(ht->lut);
    free(ht->val);
    free(ht->rand);
    free(ht->deg);
    free(ht->div);
    free(ht->idx);
#if __GB_HAVE_SSE2
    for (i=0; i<ht->sz; ++i)
      free(ht->ev[i]);
    free(ht->ev);
#else
    for (i=0; i<ht->sz; ++i) {
      free(ht->exp[i]);
    }
    free(ht->exp);
#endif
  }

  free(ht);
  ht      = NULL;
  *ht_in  = ht;
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
static inline hash_t get_hash(const exp_v *ev, const mp_cf4_ht_t *ht)
{

  exp_t exp[ht->nev * ht->vl] __attribute__ ((aligned (16)));
  exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
  for (int i=0; i<ht->nev; ++i) {
    _mm_store_si128((exp_v *)tmp, ev[i]);
    memcpy(exp+(i*ht->vl), tmp, ht->vl*sizeof(exp_t));
  }
  hash_t hash  = 0;
  for (int i=0; i<ht->nv; ++i)
    hash  +=  ht->rand[i] * exp[i];
#if HASH_DEBUG
  for (nelts_t i=0; i<ht->nv; ++i)
    printf("%u ", exp[i]);
  printf(" --> %lu\n", hash);
#endif
  return hash;
}
#else
static inline hash_t get_hash(const exp_t *exp, const mp_cf4_ht_t *ht)
{
  nvars_t i;
  hash_t hash  = 0;
  for (i=0; i<ht->nv; ++i)
    hash  +=  ht->rand[i] * exp[i];
#if HASH_DEBUG
  for (nelts_t i=0; i<ht->nv; ++i)
    printf("%u ", exp[i]);
  printf(" --> %lu\n", hash);
#endif
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
#if HASH_DEBUG
  for (int i=0; i<ht->nv; ++i)
    printf("%u ",ht->exp[ht->load][i]);
  printf(" ||| ");
  printf("%11u | %11u | %5u\n",hash, ht->load, ht->deg[ht->load]);
#endif
  ht->load++;

  // we need to keep one place open in ht->exp since the next element to be
  // checked against the hash table will be intermediately stored there
  if (ht->load >= ht->sz)
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

  /*
#if !__GB_HAVE_SSE2
  nvars_t i;
  for (i=0; i<ht->nv; ++i)
    ht->exp[last_pos][i] = ht->exp[mon_1][i] + ht->exp[mon_2][i];
#endif
*/
  // ht->div and ht->idx are already initialized with 0, so nothing to do there
  ht->deg[last_pos] = ht->deg[mon_1] + ht->deg[mon_2];
  ht->val[last_pos] = hash;
  ht->lut[pos]      = last_pos;
#if HASH_DEBUG
  int i;
  for (unsigned long j=0; j<ht->load; ++j) {
    for (i=0; i<ht->nv; ++i) {
      if (ht->exp[ht->load][i] != ht->exp[j][i]) {
        break;
      }
    }
    if (i==ht->nv) {
      for (i=0; i<ht->nv; ++i)
        printf("%u ",ht->exp[ht->load][i]);
      printf(" ||| ");
      printf("%11u | %11u\n",hash, last_pos);
      printf("------------------------------\n");
      for (i=0; i<ht->nv; ++i)
        printf("%u ",ht->exp[j][i]);
      printf(" ||| ");
      printf("%11u | %11u\n", ht->val[j], j);
      printf("==============================\n");
    }
  }
#endif
  // we do not need this anymore since it is already computed and stored in
  // check_in_hash_table_product()
  /*
#if __GB_HAVE_SSE2
  ht->ev[last_pos]  = _mm_adds_epu8(ht->ev[mon_1], ht->ev[mon_2]);
#endif
  */
  ht->load++;

  if (ht->load >= ht->sz)
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
  nvars_t i;
  // element to be checked, intermediately stored in the first free position of
  // ht->exp

#if __GB_HAVE_SSE2
  nvars_t j;
  hash_t hash = get_hash(ht->ev[ht->load], ht);
#else
  exp_t *exp  = ht->exp[ht->load];
  hash_t hash = get_hash(exp, ht);
#endif
  ht_size_t tmp_h;
  ht_size_t tmp_l;         // temporary lookup table value
  tmp_h = hash & (ht->sz-1);

  // first check directly
  tmp_l = ht->lut[tmp_h];
  if (tmp_l == 0)
    return insert_in_hash_table(hash, tmp_h, ht);
  if (ht->val[tmp_l] == hash) {
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi32(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    if (j == ht->nev)
      return tmp_l;
#else
#if __GB_USE_64_EXP_VEC
    if (exp[0] == ht->exp[tmp_l][0])
      return tmp_l;
#else
    if (memcmp(exp, ht->exp[tmp_l], ht->nv*sizeof(exp_t)) == 0)
      return tmp_l;
    /*
    for (j=0; j<ht->nv; ++j)
      if (exp[j] != ht->exp[tmp_l][j])
        break;
    if (j == ht->nv) {
      return tmp_l;
    }
    */
#endif
#endif
  }
#if HASH_DEBUG
  for (i=0; i<ht->nv; ++i)
    printf("%u ",exp[i]);
  printf("\nhash = %u\n",hash);
#endif

  // remaining checks with probing
  for (i=1; i<ht->sz; ++i) {
    tmp_h = (tmp_h+i) & (ht->sz-1);
    tmp_l = ht->lut[tmp_h];
    if (tmp_l == 0)
      break;
    if (ht->val[tmp_l] != hash)
      continue;
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi32(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    if (j == ht->nev)
      return tmp_l;
#else
#if __GB_USE_64_EXP_VEC
    if (exp[0] == ht->exp[tmp_l][0])
      return tmp_l;
#else
    if (memcmp(exp, ht->exp[tmp_l], ht->nv*sizeof(exp_t)) == 0)
      return tmp_l;
    /*
    nvars_t j;
    for (j=0; j<ht->nv; ++j)
      if (exp[j] != ht->exp[tmp_l][j])
        break;
    if (j == ht->nv) {
      return tmp_l;
    }
    */
#endif
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
 * \note We use this also in the parallel construction of gbla matrices, thus we
 * have to make this procedure threadsafe, i.e. we are NOT allowed to store the
 * product of the hashes intermediately in ht->exp[ht->load] resp.
 * ht->ev[ht->load]. Otherwise, other threads will overwrite this value!
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
  ht_size_t i;
  hash_t hash;
#if __GB_HAVE_SSE2
  ht_size_t j;
  exp_v ev[ht->nev];
  for (i=0; i<ht->nev; ++i)
    ev[i] = _mm_adds_epu8(ht->ev[mon_1][i], ht->ev[mon_2][i]);
  //hash = get_hash(ht->ev[ht->load], ht);
#else
  exp_t exp[ht->nv];
  for (i=0; i<ht->nv; ++i)
    exp[i] = ht->exp[mon_1][i] + ht->exp[mon_2][i];
  //hash = get_hash(ht->exp[ht->load], ht);
#endif
  // hash value of the product is the sum of the hash values in our setting
  hash   = ht->val[mon_1] + ht->val[mon_2];
  ht_size_t tmp_h;
  ht_size_t tmp_l;         // temporary lookup table value
  tmp_h = hash & (ht->sz-1);
  
  // first check directly
  tmp_l = ht->lut[tmp_h];
  if (tmp_l == 0)
    return 0;
  if (ht->val[tmp_l] == hash) {
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi32(ht->ev[tmp_l][j], ev[j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    if (j == ht->nev)
      return tmp_l;
#else
#if __GB_USE_64_EXP_VEC
    if (exp[0] == ht->exp[tmp_l][0])
      return tmp_l;
#else
    if (memcmp(exp, ht->exp[tmp_l], ht->nv*sizeof(exp_t)) == 0)
      return tmp_l;
    /*
    for (j=0; j<ht->nv; ++j)
      if (exp[j] != ht->exp[tmp_l][j])
        break;
    if (j == ht->nv)
      return tmp_l;
      */
#endif
#endif
  }

  // remaining checks with probing
  for (i=1; i<ht->sz; ++i) {
    tmp_h = (tmp_h+i) & (ht->sz-1);
    tmp_l = ht->lut[tmp_h];
    if (tmp_l == 0)
      break;
    if (ht->val[tmp_l] != hash)
      continue;
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi32(ht->ev[tmp_l][j], ev[j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    if (j == ht->nev)
      return tmp_l;
#else
#if __GB_USE_64_EXP_VEC
    if (exp[0] == ht->exp[tmp_l][0])
      return tmp_l;
#else
    if (memcmp(exp, ht->exp[tmp_l], ht->nv*sizeof(exp_t)) == 0)
      return tmp_l;
    /*
    for (j=0; j<ht->nv; ++j)
      if (exp[j] != ht->exp[tmp_l][j])
        break;
    if (j == ht->nv)
      return tmp_l;
    */
#endif
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
  ht_size_t i;
  hash_t hash;
#if __GB_HAVE_SSE2
  ht_size_t j;
  for (i=0; i<ht->nev; ++i)
    ht->ev[ht->load][i] = _mm_adds_epu16(ht->ev[mon_1][i], ht->ev[mon_2][i]);
  //hash = get_hash(ht->ev[ht->load], ht);
#else
  for (i=0; i<ht->nv; ++i)
    ht->exp[ht->load][i] = ht->exp[mon_1][i] + ht->exp[mon_2][i];
  //hash = get_hash(ht->exp[ht->load], ht);
#endif
  // hash value of the product is the sum of the hash values in our setting
  hash   = ht->val[mon_1] + ht->val[mon_2];
  ht_size_t tmp_h;
  ht_size_t tmp_l;         // temporary lookup table value
  tmp_h = hash & (ht->sz-1);

  // first check directly
  tmp_l = ht->lut[tmp_h];
  if (tmp_l == 0)
    return insert_in_hash_table_product(mon_1, mon_2, hash, tmp_h, ht);
  if (ht->val[tmp_l] == hash) {
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi32(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    if (j == ht->nev)
      return tmp_l;
#else
#if __GB_USE_64_EXP_VEC
    if (ht->exp[ht->load][0] == ht->exp[tmp_l][0])
      return tmp_l;
#else
    if (memcmp(ht->exp[ht->load], ht->exp[tmp_l], ht->nv*sizeof(exp_t)) == 0)
      return tmp_l;
    /*
    for (j=0; j<ht->nv; ++j)
      if (ht->exp[tmp_l][j] != ht->exp[ht->load][j])
        break;
    if (j == ht->nv)
      return tmp_l;
    */
#endif
#endif
  }

  // remaining checks with probing
  for (i=1; i<ht->sz; ++i) {
    tmp_h = (tmp_h+i) & (ht->sz-1);
    tmp_l = ht->lut[tmp_h];
    if (tmp_l == 0)
      break;
    if (ht->val[tmp_l] != hash)
      continue;
#if __GB_HAVE_SSE2
    exp_v cmpv;
    for (j=0; j<ht->nev; ++j) {
      cmpv  = _mm_cmpeq_epi32(ht->ev[tmp_l][j], ht->ev[ht->load][j]);
      if (_mm_movemask_epi8(cmpv) == 0)
        break;
    }
    if (j == ht->nev)
      return tmp_l;
#else
#if __GB_USE_64_EXP_VEC
    if (ht->exp[ht->load][0] == ht->exp[tmp_l][0])
      return tmp_l;
#else
    if (memcmp(ht->exp[ht->load], ht->exp[tmp_l], ht->nv*sizeof(exp_t)) == 0)
      return tmp_l;
    /*
    nvars_t j;
    for (j=0; j<ht->nv; ++j)
      if (ht->exp[tmp_l][j] != ht->exp[ht->load][j])
        break;
    if (j == ht->nv)
      return tmp_l;
    */
#endif
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
  exp_t exp[ht->nev * ht->vl] __attribute__ ((aligned (16)));
  exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
  for (i=0; i<ht->nev; ++i) {
    ht->ev[ht->load][i] = _mm_max_epu16(ht->ev[h1][i], ht->ev[h2][i]);
    _mm_store_si128((exp_v *)tmp, ht->ev[ht->load][i]);
    memcpy(exp+(i*ht->vl), tmp, ht->vl*sizeof(exp_t));
  }
  ht->deg[ht->load] = 0;
  for (i=0; i<ht->nv; ++i)
    ht->deg[ht->load] += exp[i];
  /*
  printf("LCM ");
  for (int ii=0; ii<ht->nv; ++ii)
    printf("%u ", exp[ii]);
  printf("\n");
  printf("deg %u\n", deg);
  */
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
  /*
  printf("LCM ");
  for (int ii=0; ii<ht->nv; ++ii)
    printf("%u ", lcm[ii]);
  printf("\n");
  printf("deg %u\n", deg);
  */
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
  if (ht->deg[h1] < ht->deg[h2])
    return 0;
  nvars_t i;
#if __GB_HAVE_SSE2
  exp_v cmpv;
  for (i=0; i<ht->nev; ++i) {
    cmpv  = _mm_cmplt_epi16(ht->ev[h1][i], ht->ev[h2][i]);
    if (_mm_movemask_epi8(cmpv) != 0)
      return 0;
  }
  for (i=0; i<ht->nev; ++i)
    ht->ev[ht->load][i] = _mm_subs_epu16(ht->ev[h1][i], ht->ev[h2][i]);
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
 * \return 0 if not divisible, 1 if divisible
 */
static inline int check_monomial_division(const hash_t h1, const hash_t h2, const mp_cf4_ht_t *ht)
{
  if (ht->deg[h1] < ht->deg[h2])
    return 0;
  nvars_t i;
#if __GB_HAVE_SSE2
  exp_v cmpv;
  for (i=0; i<ht->nev; ++i) {
    cmpv  = _mm_cmplt_epi16(ht->ev[h1][i], ht->ev[h2][i]);
    if (_mm_movemask_epi8(cmpv) != 0)
      return 0;
  }
  return 1;
#else
  const exp_t * const e1  = ht->exp[h1];
  const exp_t * const e2  = ht->exp[h2];

  for (i=0; i<ht->nv; ++i) {
    if (e1[i] < e2[i])
      return 0;
  }
  return 1;
#endif
}

/**
 * \brief Tests if exp of h1 is divisible by exp of h2. If divisibility is
 * fulfilled only 1 is returned, else 0. This is a special version where the
 * homogenization variable (last variable) is not taking care of. It is needed
 * for redundancy checks when saturating polynomials in homogenized
 * computations.
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
 * \return 0 if not divisible, 1 if divisible (w.r.t to all variables besides
 * the last one)
 */
static inline int check_monomial_division_saturated(const hash_t h1, const hash_t h2, const mp_cf4_ht_t *ht)
{
  // do not do the degree check: for saturated polynomials we have not computed
  // the correct degree!
  /*
  if (ht->deg[h1] < ht->deg[h2])
    return 0;
  */
  nvars_t i;
#if __GB_HAVE_SSE2
  exp_t exp1[ht->nev * ht->vl] __attribute__ ((aligned (16)));
  exp_t exp2[ht->nev * ht->vl] __attribute__ ((aligned (16)));
  exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
  for (i=0; i<ht->nev; ++i) {
    _mm_store_si128((exp_v *)tmp, ht->ev[h1][i]);
    memcpy(exp1+(i*ht->vl), tmp, ht->vl*sizeof(exp_t));
    _mm_store_si128((exp_v *)tmp, ht->ev[h2][i]);
    memcpy(exp2+(i*ht->vl), tmp, ht->vl*sizeof(exp_t));
  }
#else
  const exp_t * const exp1  = ht->exp[h1];
  const exp_t * const exp2  = ht->exp[h2];
#endif

  // note that we explicitly do not check w.r.t. the last variable!
  for (i=0; i<ht->nv-1; ++i) {
    if (exp1[i] < exp2[i])
      return 0;
  }
  return 1;
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
static inline hash_t get_multiplier(const hash_t h1, const hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
#if __GB_HAVE_SSE2
  for (i=0; i<ht->nev; ++i)
    ht->ev[ht->load][i] = _mm_subs_epu16(ht->ev[h1][i], ht->ev[h2][i]);
#else
  exp_t *e  = ht->exp[ht->load];
  const exp_t * const e1  = ht->exp[h1];
  const exp_t * const e2  = ht->exp[h2];

  // we know that exp e2 divides exp e1, so no check for e1[i] < e2[i]
  for (i=0; i<ht->nv; ++i) {
    e[i]  = e1[i] - e2[i];
  }
  /*
  for (i=0; i<ht->nv; i=i+2) {
    e[i]  = e1[i] - e2[i];
    e[i+1]  = e1[i+1] - e2[i+1];
  }
  */
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
  memset(ht->idx, 0, ht->sz * sizeof(ht_size_t));
}
#endif /* GB_HASH_TABLE_H */
