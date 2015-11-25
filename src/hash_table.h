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
void pseudo_random_generator()
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
void set_random_seed(mp_cf4_ht_t *hash_table);

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
mp_cf4_ht_t *init_hash_table(const ht_size_t hash_table_size,
    const nvars_t number_variables);

/**
 * \brief Frees dynamically allocated data in hash table
 *
 * \param hash table hash_table
 */
void free_hash_table(mp_cf4_ht_t *hash_table);

/**
 * \brief Get hash value
 *
 * \param exponent vector exp
 *
 * \param hash table hash_table
 *
 * \return hash value
 */
hash_t get_hash(const exp_t *exp, mp_cf4_ht_t *ht);

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
hash_t insert_with_quadratic_probing(const exp_t *exp,
    const hash_t hash, mp_cf4_ht_t *ht);

/**
 * \brief Inserts in hash table using quadratic probing when enlarging table
 *
 * \param hash value to be inserted hash
 *
 * \param hash table hash_table
 */
void insert_while_enlarging(const hash_t hash, mp_cf4_ht_t *ht);

/**
 * \brief Inserts a new element to the hash table
 *
 * \param monomial exponent to be inserted exp
 *
 * \param hash value of exp hash
 *
 * \param position in lookup table pos
 *
 * \param hash table ht
 *
 * \return position of hash of exp in table
 */
hash_t insert_in_hash_table(const exp_t *exp, const hash_t hash,
    const hash_t pos,  mp_cf4_ht_t *ht);

/**
 * \brief Inserts a new element to the hash table coming from a product of two
 * monomials.
 *
 * \note We use the sum of the hash values of mon_1 and mon_2 as hash value for
 * the product of mon_1 and mon_2.
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
hash_t insert_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    const hash_t hash, const hash_t pos,  mp_cf4_ht_t *ht);

/**
 * \brief Checks if the given monomial exponent is already in the hash table. If
 * not, it is added to the table
 *
 * \param monomial exponent to be inserted exp
 *
 * \param hash table ht
 *
 * \return position of hash of exp in table
 */
hash_t check_in_hash_table(const exp_t *exp, mp_cf4_ht_t *ht);

/**
 * \brief Checks if the product of the given two monomial exponents is already
 * in the hash table. If not, it is added to the table
 *
 * \note We use the sum of the hash values of mon_1 and mon_2 as hash value for
 * the product of mon_1 and mon_2.
 *
 * \param monomial 1 mon_1
 *
 * \param monomial 2 mon_2
 *
 * \param hash table ht
 *
 * \return position of hash of exp in table
 */
hash_t check_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    mp_cf4_ht_t *ht);

/**
 * \brief Inserts elements in hash table during the elargement of the table
 *
 * \param hash table hash_table
 *
 * \param new size of hash table new_size
 */
void enlarge_hash_table(mp_cf4_ht_t *hash_table, const hash_t new_size);
#endif /* GB_HASH_TABLE_H */
