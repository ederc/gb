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
#include <omp.h>
#include "types.h"

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#define HASH_DEBUG  0

/***************************
 * OUR HASH TABLE IS GLOBAL
 **************************/
mp_cf4_ht_t *ht;

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
void free_hash_table_dynamic_data(mp_cf4_ht_t *hash_table);

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
hash_t insert_in_hash_table(const hash_t hash,
    const hash_t pos,  mp_cf4_ht_t *ht);

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
hash_t insert_in_hash_table_product(const hash_t mon_1, const hash_t mon_2,
    const hash_t hash, const hash_t pos,  mp_cf4_ht_t *ht);

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
hash_t check_in_hash_table(mp_cf4_ht_t *ht);

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
hash_t get_lcm(hash_t h1, hash_t h2, mp_cf4_ht_t *ht);

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
hash_t monomial_division(hash_t h1, hash_t h2, mp_cf4_ht_t *ht);

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
hash_t get_multiplier(hash_t h1, hash_t h2, mp_cf4_ht_t *ht);
#endif /* GB_HASH_TABLE_H */
