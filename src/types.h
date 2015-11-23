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
 * \file types.h
 * \brief General typedefs
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_TYPES_H
#define GB_TYPES_H

#if 0
#define GB_USE_FLOAT XXX
/* #define GB_USE_INT16 XXX */
#else
#define GB_USE_UINT16 OK
//#define GB_USE_UINT32 OK
#endif
/* #define GB_USE_UINT32 OK */
/* #define GB_USE_INT32 */
/* #define GB_USE_AVX */

/* number of variables */
typedef uint16_t nvars_t;

/* number of elements */
typedef uint32_t nelts_t;

/* hash table table size */
typedef uint32_t ht_size_t;

/* hash table entry size */
typedef uint32_t hash_t;

/* degree size */
typedef uint8_t deg_t;

#ifdef GB_USE_FLOAT
/** coefficient storage type */
typedef float coeff_t;
/** type of field characteristic */
typedef float mod_t;
#endif

#ifdef GB_USE_DOUBLE
/** coefficient storage type */
typedef double coeff_t;
/** type of field characteristic */
typedef double mod_t;
#endif
#ifdef GB_USE_UINT16
/** coefficient storage type */
typedef uint16_t coeff_t;
/** type of field characteristic */
typedef uint16_t mod_t;
#endif
#ifdef GB_USE_INT16
/** coefficient storage type */
typedef int16_t coeff_t;
/** type of field characteristic */
typedef int16_t mod_t;
#endif
#ifdef GB_USE_UINT32
/** coefficient storage type */
typedef uint32_t coeff_t;
/** type of field characteristic */
typedef uint32_t mod_t;
#endif
#ifdef GB_USE_INT32
/** coefficient storage type */
typedef int32_t coeff_t;
/** type of field characteristic */
typedef int32_t mod_t;
#endif


#define ALIGNT 32

#ifdef GB_USE_DOUBLE
#define MODP(a,b) \
	fmod((a),(b))
	/* static double MODP(double a, double b) { assert(a>=0) ; assert(b>0) ; double c = fmod(a,b) ; assert(c >=0) ; return c; } */
#define CAST(a) \
	(double) (a)
#endif

#ifdef GB_USE_FLOAT
#define MODP(a,b) \
	fmod((a),(b))
	/* static double MODP(double a, double b) { assert(a>=0) ; assert(b>0) ; double c = fmod(a,b) ; assert(c >=0) ; return c; } */
#define CAST(a) \
	(double) (a)
#endif

#ifdef GB_USE_UINT16
#define MODP(a,b) \
	(a) % (b)
#define CAST(a) \
	(a) & 0x000000000000ffff
#endif

#ifdef GB_USE_INT16
#define MODP(a,b) \
	(a) % (b)
#define CAST(a) \
	(a)
#endif


#ifdef GB_USE_UINT32
#define MODP(a,b) \
	(a) % (b)
#define CAST(a) \
	(a) & 0x00000000ffffffff
#endif

#ifdef GB_USE_INT32
#define MODP(a,b) \
	(a) % (b)
#define CAST(a) \
	(a) & 0x00000000ffffffff
#endif


/******************************************
 * general data structures needed in gb
 *****************************************/

/**
 * \brief Verbosity information
 */
typedef struct info_t
{
  nelts_t n_reduced;                          /*!<  number of elements reduced*/
  nelts_t n_pairs_removed;                    /*!<  number of pairs removed by criteria*/
  nelts_t n_zero_reductions;                  /*!<  number of zero reductions*/
  struct timeval read_input_time;             /*!<  time for reading input data*/
  struct timeval criteria_check_time;         /*!<  overall time for criteria checks*/
  struct timeval pair_selection_time;         /*!<  overall time for pair selection*/
  struct timeval symbolic_preprocessing_time; /*!<  overall time for symbolic
                                                    preprocessing*/
  struct timeval poly_to_matrix_time;         /*!<  overall time for matrix construction
                                                    from poly data*/
  struct timeval matrix_to_poly_time;         /*!<  overall time for poly construction
                                                    from matrix data*/
  struct timeval linear_algebra_time;         /*!<  overall time for reduction in
                                                    linear algebra*/
} info_t;


/**
 * \brief Input generators for gb computation
 */
typedef struct gen_t
{
  coeff_t **coeffs; /*!<  coefficients of input elements*/
  nelts_t size;     /*!<  number of elements */
  nvars_t n_vars;   /*!<  number of variables */
  char* var_names;  /*!<  variable names */
  mod_t modulus;    /*!<  modulo/field characteristic */
  float fs;         /*!<  file size of input matrix */
  char *fsu;        /*!<  file size unit of input matrix, e.g. GB */
} gen_t;

/**
 * \brief S-pairs resp. S-polynomials list
 */
typedef struct spair_t
{
  nelts_t gen1; /*!<  index to first generator*/
  nelts_t gen2; /*!<  index to second generator*/
  deg_t degree; /*!<  degree of S-pair*/
  hash_t lcm;   /*!<  hash of lcm of the lead terms of gen1 and gen2*/
} spair_t;

/**
 * \brief Groebner basis data
 */
typedef struct basis_t
{
  coeff_t **coeffs; /*!<  coefficients of gb elements*/
  nelts_t size;     /*!<  size of gb*/
  nelts_t load;     /*!<  current load of gb, i.e. number of elements in gb*/
} basis_t;

/**
 * \brief Hash table as defined by Monagan and Pearce in compact F4
 * implementation (see PASCO 2015)
 */
typedef struct mp_cf4_ht_t
{
  // size and load counters
  ht_size_t size; /*!<  size of hash table*/
  ht_size_t load; /*!<  load of hash table*/
  // data and arrays for storage
  hash_t *value;  /*!<  array of hash values*/
  hash_t *exp;    /*!<  array of exponents*/
  hash_t *rand;   /*!<  array of random values for each variable
                        to generate hash values out of exponents*/
  deg_t *deg;     /*!<  degree of monmial, for faster sorting and searching*/
  nelts_t *div;   /*!<  latest element from gb checked for its leading
                        term dividing corresponding monomial in hash table*/
  hash_t *idx;    /*!<  index used for matrix generation and marking
                        monomials already taken care of in symbolic
                        preprocessing*/
} mp_cf4_ht_t;

#endif /* GB_TYPES_H */
