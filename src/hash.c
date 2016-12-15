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
/** extern declaration in src/hash.h */
mp_cf4_ht_t *ht;

/********************************************************************************
 * FOLLOWING HASH TABLE IMPLEMENTATION IS COPIED FROM COMPACT F4 IMPLEMENTATION 
 * BY MONAGAN AND PIERCE (see PASCO 2015)
 *******************************************************************************/
