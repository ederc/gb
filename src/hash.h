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
#define HASH_QUADRATIC_PROBING  1
#endif

/***************************
 * OUR HASH TABLE IS GLOBAL
 **************************/
extern mp_cf4_ht_t *ht;

/**
 * \brief Generates pseudo random numbers using xorshifts using global defined
 * random_seed variable
 *
 * \return some pseudo random number
 */
void pseudo_random_generator();

/********************************************************************************
 * FOLLOWING HASH TABLE IMPLEMENTATION IS COPIED FROM COMPACT F4 IMPLEMENTATION 
 * BY MONAGAN AND PIERCE (see PASCO 2015)
 *******************************************************************************/
/**
 * \brief Inserts random seeds in hash table random array
 *
 * \param hash table ht
 */
void set_random_seed(mp_cf4_ht_t *ht);

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
mp_cf4_ht_t *init_hash_table(const ht_size_t ht_si,
    const nvars_t nv);

/**
 * \brief Enlarges hash table to next prime number size
 *
 * \param hash table ht
 */
void enlarge_hash_table(mp_cf4_ht_t *ht);

/**
 * \brief Frees dynamically allocated data in hash table
 *
 * \param hash table ht
 */
inline void free_hash_table(mp_cf4_ht_t *ht)
{
  if (ht) {

    hash_t i;

    free(ht->lut);
    free(ht->val);
    free(ht->rand);
    free(ht->deg);
    free(ht->div);
    free(ht->idx);
    free(ht->primes);
    for (i=0; i<ht->primes[ht->si]; ++i)
      free(ht->exp[i]);
    free(ht->exp);
#if HAVE_SSE2
    free(ht->ev);
#endif
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
inline hash_t get_hash(const exp_t *exp, mp_cf4_ht_t *ht)
{
  hash_t i;
  hash_t hash = 0;

  for (i=0; i<ht->nv; ++i)
    hash  +=  ht->rand[i] * (hash_t)exp[i];

  return hash;
}

/**
 * \brief Inserts in hash table using quadratic probing in general
 *
 * \param exponent vector to be inserted exp
 *
 * \param hash value to be inserted hash
 *
 * \param hash table hash_table
 *
 * \return hash value
 */
inline hash_t insert_with_quadratic_probing(const exp_t *exp,
    const hash_t hash, mp_cf4_ht_t *ht)
{
  hash_t i, j, tmp;

  tmp = hash;
  for (i=0; i<ht->primes[ht->si]; ++i) {
    tmp = MODP(tmp+(i*i),ht->primes[ht->si]);
    if (!ht->lut[tmp])
      break;
    if (ht->val[ht->lut[tmp]] != hash)
      continue;
    for (j=0; j<ht->nv; ++j) {
      if (ht->exp[ht->lut[tmp]][j] != exp[j])
        break;
    }
    if (j == ht->nv)
      return ht->lut[tmp];
  }
  return ht->lut[tmp];
}

/**
 * \brief Inserts in hash table using linear probing in general
 *
 * \param exponent vector to be inserted exp
 *
 * \param hash value to be inserted hash
 *
 * \param hash table hash_table
 *
 * \return hash value
 */
inline hash_t insert_with_linear_probing(const exp_t *exp,
    const hash_t hash, mp_cf4_ht_t *ht)
{
  hash_t i, j, tmp;

  tmp = hash;
  for (i=0; i<ht->primes[ht->si]; ++i) {
    tmp = MODP(tmp+(i),ht->primes[ht->si]);
    if (!ht->lut[tmp])
      break;
    if (ht->val[ht->lut[tmp]] != hash)
      continue;
    for (j=0; j<ht->nv; ++j) {
      if (ht->exp[ht->lut[tmp]][j] != exp[j])
        break;
    }
    if (j == ht->nv)
      return ht->lut[tmp];
  }
  return ht->lut[tmp];
}

/**
 * \brief Inserts in hash table using quadratic probing when enlarging table
 *
 * \param hash value to be inserted hash
 *
 * \param hash table hash_table
 */
inline void insert_while_enlarging(const hash_t hash, mp_cf4_ht_t *ht)
{
  hash_t i, tmp;

  tmp = hash;
  for (i=0; i<ht->primes[ht->si]; ++i) {
#if HASH_QUADRATIC_PROBING
    tmp = MODP(tmp+(i*i),ht->primes[ht->si]);
#else
    tmp = MODP(tmp+(i),ht->primes[ht->si]);
#endif
    if (ht->lut[tmp])
      continue;
    ht->lut[tmp]  = i;
    break;
  }
}

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
inline hash_t insert_in_hash_table(const hash_t hash,
    const hash_t pos,  mp_cf4_ht_t *ht)
{
  nvars_t i;
  // the new exponent is already stored in ht->exp[ht->load]
  // deg is possibly not zero due to some earlier check in the hash table
  ht->deg[ht->load] = 0;
  for (i=0; i<ht->nv; ++i)
    ht->deg[ht->load] +=  ht->exp[ht->load][i];
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
inline hash_t insert_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    const hash_t hash, const hash_t pos,  mp_cf4_ht_t *ht)
{
  hash_t i;

  hash_t last_pos = ht->load;

  for (i=0; i<ht->nv; ++i)
    ht->exp[last_pos][i] = ht->exp[mon_1][i] + ht->exp[mon_2][i];

  // ht->div and ht->idx are already initialized with 0, so nothing to do there
  ht->deg[last_pos] = ht->deg[mon_1] + ht->deg[mon_2];
  ht->val[last_pos] = hash;
  ht->lut[pos]      = last_pos;

#if HAVE_SSE2
  exp_t *exp  = (exp_t *)calloc(16, sizeof(exp_t));
  for (i=0; i<ht->nv; ++i)
    exp[i]  = ht->exp[last_pos][i];
  ht->ev[last_pos] = _mm_loadu_si128((exp_v *)exp);
  free(exp);
#endif
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
inline hash_t check_in_hash_table(mp_cf4_ht_t *ht)
{
  hash_t i,j;
  // element to be checked, intermediately stored in the first free position of
  // ht->exp
  exp_t *exp  = ht->exp[ht->load];

  hash_t hash     = get_hash(exp, ht);
  ht_size_t tmp_h = (ht_size_t)MODP(hash, ht->primes[ht->si]); // temporary hash values for linear probing
  ht_size_t tmp_l;         // temporary lookup table value

#if HASH_DEBUG
  for (i=0; i<ht->nv; ++i)
    printf("%u ",exp[i]);
  printf("\nhash = %u\n",hash);
#endif

  for (i=0; i<ht->primes[ht->si]; ++i) {
#if HASH_QUADRATIC_PROBING
    tmp_h = MODP(tmp_h+(i*i),ht->primes[ht->si]);
#else
    tmp_H = MODP(tmp_h+(i),ht->primes[ht->si]);
#endif
    tmp_l = ht->lut[tmp_h];
    if (tmp_l == 0)
      break;
    if (ht->val[tmp_l] != hash)
      continue;
#if HAVE_SSE2
    exp_v cmpv  = _mm_cmpeq_epi64(ht->ev[tmp_l], ht->ev[ht->load]);
    if (_mm_movemask_epi8(cmpv) == 0xFFFF) {
      return tmp_l;
    }
#else
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
inline hash_t find_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    const mp_cf4_ht_t *ht)
{
  hash_t i,j;

  // hash value of the product is the sum of the hash values in our setting
  hash_t hash   = ht->val[mon_1] + ht->val[mon_2];
  ht_size_t tmp_h = (ht_size_t)MODP(hash, ht->primes[ht->si]); // temporary hash values for linear probing
  ht_size_t tmp_l;         // temporary lookup table value
#if HAVE_SSE2
  exp_v prod  = _mm_adds_epu8(ht->ev[mon_1], ht->ev[mon_2]);
#endif

  for (i=0; i<ht->primes[ht->si]; ++i) {
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
#if HAVE_SSE2
    exp_v cmpv  = _mm_cmpeq_epi64(ht->ev[tmp_l], prod);
    if (_mm_movemask_epi8(cmpv) == 0xFFFF) {
      return tmp_l;
    }
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
inline hash_t check_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    mp_cf4_ht_t *ht)
{
  hash_t i,j;

  // hash value of the product is the sum of the hash values in our setting
  hash_t hash   = ht->val[mon_1] + ht->val[mon_2];
  ht_size_t tmp_h = (ht_size_t)MODP(hash, ht->primes[ht->si]); // temporary hash values for linear probing
  ht_size_t tmp_l;         // temporary lookup table value

#if HAVE_SSE2
  exp_v prod  = _mm_adds_epu8(ht->ev[mon_1], ht->ev[mon_2]);
#endif
  for (i=0; i<ht->primes[ht->si]; ++i) {
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
#if HAVE_SSE2
    exp_v cmpv  = _mm_cmpeq_epi64(ht->ev[tmp_l], prod);
    if (_mm_movemask_epi8(cmpv) == 0xFFFF) {
      return tmp_l;
    }
#else
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
inline hash_t get_lcm(hash_t h1, hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
  exp_t *lcm, *e1, *e2;

  // use first free entry in hash table ht to store possible new lcm monomial
  lcm = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  for (i=0; i<ht->nv; ++i) {
    lcm[i]  = e1[i] < e2[i] ? e2[i] : e1[i];
  }
#if HAVE_SSE2
  exp_t *exp  = (exp_t *)calloc(16, sizeof(exp_t));
  for (i=0; i<ht->nv; ++i)
    exp[i]  = ht->exp[ht->load][i];
  ht->ev[ht->load] = _mm_loadu_si128((exp_v *)exp);
  free(exp);
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
inline hash_t monomial_division(hash_t h1, hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
  exp_t *e, *e1, *e2;

  e   = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  for (i=0; i<ht->nv; ++i) {
    if (e1[i] < e2[i])
      return 0;
    e[i]  = e1[i] - e2[i];
  }
#if HAVE_SSE2
  exp_t *exp  = (exp_t *)calloc(16, sizeof(exp_t));
  for (i=0; i<ht->nv; ++i)
    exp[i]  = ht->exp[ht->load][i];
  ht->ev[ht->load] = _mm_loadu_si128((exp_v *)exp);
  free(exp);
#endif
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
inline hash_t check_monomial_division(hash_t h1, hash_t h2, const mp_cf4_ht_t *ht)
{
  nvars_t i;
  exp_t *e, *e1, *e2;

  e   = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  for (i=0; i<ht->nv; ++i) {
    if (e1[i] < e2[i])
      return 0;
    e[i]  = e1[i] - e2[i];
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
inline hash_t get_multiplier(hash_t h1, hash_t h2, mp_cf4_ht_t *ht)
{
  nvars_t i;
  exp_t *e, *e1, *e2;

  e   = ht->exp[ht->load];
  e1  = ht->exp[h1];
  e2  = ht->exp[h2];

  // we know that exp e2 divides exp e1, so no check for e1[i] < e2[i]
  for (i=0; i<ht->nv; ++i)
    e[i]  = e1[i] - e2[i];
#if HAVE_SSE2
  exp_t *exp  = (exp_t *)calloc(16, sizeof(exp_t));
  for (i=0; i<ht->nv; ++i)
    exp[i]  = ht->exp[ht->load][i];
  ht->ev[ht->load] = _mm_loadu_si128((exp_v *)exp);
  free(exp);
#endif
  return check_in_hash_table(ht);
}

/**
 * \brief Resets all idx entries of the hash table to zero
 *
 * \param hash table ht
 */
inline void clear_hash_table_idx(mp_cf4_ht_t *ht)
{
  memset(ht->idx, 0, ht->primes[ht->si] * sizeof(hash_t));
}
#endif /* GB_HASH_TABLE_H */
