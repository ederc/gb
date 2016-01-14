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

#include <gbla/matrix.h>

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

/* exponent size */
typedef uint8_t exp_t;

#ifdef GB_USE_FLOAT
/** coefficient storage type */
typedef float coeff_t;
#endif

#ifdef GB_USE_DOUBLE
/** coefficient storage type */
typedef double coeff_t;
#endif
#ifdef GB_USE_UINT16
/** coefficient storage type */
typedef uint16_t coeff_t;
#endif
#ifdef GB_USE_INT16
/** coefficient storage type */
typedef int16_t coeff_t;
#endif
#ifdef GB_USE_UINT32
/** coefficient storage type */
typedef uint32_t coeff_t;
#endif
#ifdef GB_USE_INT32
/** coefficient storage type */
typedef int32_t coeff_t;
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
  nelts_t nred_last;                          /*!<  number of elements reduced in last step*/
  nelts_t nred_total;                         /*!<  number of elements reduced in last step*/
  nelts_t ncrit_last;                         /*!<  number of pairs removed by criteria in the last step*/
  nelts_t ncrit_total;                        /*!<  number of pairs removed by criteria in total*/
  nelts_t nzerored_last;                      /*!<  number of zero reductions in the last step*/
  nelts_t nzerored_total;                     /*!<  number of zero reductions in total*/
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
 * \brief Groebner basis, starts with input data
 * 
 * \note For easier divisibility checks in symbolic preprocessing we put at
 * the first position of basis, i.e. index 0 a NULL element.
 * Thus, basis->load is always one bigger than the actual number of elements
 * in the basis.
 */
typedef struct gb_t
{
  // global data
  nelts_t size;     /*!<  memory allocated */
  nelts_t load;     /*!<  number of elements in basis*/
  nvars_t nvars;    /*!<  number of variables */
  coeff_t modulus;  /*!<  modulo/field characteristic */
  // element data
  nelts_t *nt;      /*!<  number of terms in each element resp. polynomial*/
  deg_t *deg;       /*!<  degree of each element resp. polynomial*/
  coeff_t **cf;     /*!<  coefficients of input elements*/
  hash_t **eh;      /*!<  monomial exponent hash*/
  // meta data
  char **vnames;    /*!<  variable names */
  double fs;         /*!<  file size of input matrix */
  char *fsu;        /*!<  file size unit of input matrix, e.g. GB */
} gb_t;

/**
 * \brief Typedef of enum for criteria recognition: no criterion applies or the
 * product criterion applies or the chain criterion applies.
 */
typedef enum {NO_CRIT, CHAIN_CRIT, PROD_CRIT} criteria_t;

/**
 * \brief S-pairs resp. S-polynomials list
 */
typedef struct spair_t
{
  nelts_t gen1;     /*!<  index to first generator*/
  nelts_t gen2;     /*!<  index to second generator*/
  deg_t deg;        /*!<  degree of S-pair*/
  hash_t lcm;       /*!<  hash of lcm of the lead terms of gen1 and gen2*/
  criteria_t crit;  /*!<  tracker if product and chain criterion applies*/
} spair_t;

/**
 * \brief Pair set
 */
typedef struct ps_t
{
  // global data
  nelts_t size;     /*!<  memory allocated */
  nelts_t load;     /*!<  number of elements in basis*/
  // element data
  spair_t **pairs;   /*!<  pointers of spairs */
} ps_t;

/**
 * \brief Multiplier-polynomial pair for symbolic preprocessing
 */
typedef struct mpp_t
{
  nelts_t idx;  /*!<  index of polynomial in basis */
  hash_t mul;   /*!<  hash of multiplier */
  hash_t mlm;   /*!<  hash of multiplied leading monomial, needed for faster
                      sorting of rows when generating the gbla matrix */
} mpp_t;

/**
 * \brief Selection of multiplied elements for next matrix: Those are generators
 * of spairs and corresponding reducers.
 */
typedef struct sel_t
{
  // global data
  deg_t deg;        /*!<  maximal degree of all elements in selection set*/
  nelts_t size;     /*!<  memory allocated */
  nelts_t load;     /*!<  number of elements in selection*/
  // element data
  mpp_t *mpp;       /*!<  multiplier-polynomial-pair for symbolic preprocessing */
} sel_t;

/**
 * \brief Data structure storing the list of positions in hash table
 * corresponding to monomials. Used for symbolic preprocessing to get all
 * possible reducers.
 */
typedef struct pre_t
{
  hash_t *hpos;   /*!<  position of monomials in hash table*/
  nelts_t size;   /*!<  size of list*/
  nelts_t load;   /*!<  number of elements already stored in list*/
  nelts_t nlm;    /*!<  number of leading monomials*/
} pre_t;

/**
 * \brief Structure storing the full information of the symbolic preprocessing,
 * that means, it consists of the list of selected polynomials from basis
 * including the corresponding multipliers (sel) and the list of appearing
 * monomials including the number of their appearances (mon). Both of these data
 * is essential for the generation of the corresponding GBLA matrix later on.
 */
typedef struct spd_t
{
  sel_t *selu;  /*!<  selected polynomials and their  multipliers for upper
                      part of gbla matrix */
  sel_t *sell;  /*!<  selected polynomials and their  multipliers for lower
                      part of gbla matrix */
  pre_t *col;   /*!<  list of monomials appearing in selection, i.e. columns
                      of matrix */
} spd_t;

/**
 * \brief Structure storing all 4 matrix parts for gbla matrix reduction. Uses
 * gbla internal types:
 *
 * A | B
 * -----
 * C | D
 *
 * After reduction rows of D will correspond to new polynomials to be added to
 * the groebner basis.
 */
typedef struct mat_t
{
  sb_fl_t *A;   /*!<  upper left sparse matrix part*/
  dbm_fl_t *B;  /*!<  upper right dense block matrix part*/
  sb_fl_t *C;   /*!<  lower right sparse matrix part*/
  dbm_fl_t *D;  /*!<  lower right dense block matrix part*/
} mat_t;

/**
 * \brief Hash table as defined by Monagan and Pearce in compact F4
 * implementation (see PASCO 2015)
 *
 * \note We start with the first element at index 1, not at index 0. Thus we can
 * optimize divisibility checks (if index 0 is returned we have not found any
 * divisor). See also procedure init_hash_table().
 */
typedef struct mp_cf4_ht_t
{
  // size and load counters
  ht_size_t nvars;  /*!<  number of variables*/
  ht_size_t size;   /*!<  size of hash table*/
  ht_size_t load;   /*!<  load of hash table*/
  // data and arrays for storage
  hash_t *lut;      /*!<  lookup table between hash value and position in
                          exponent array*/
  hash_t *val;      /*!<  array of hash values*/
  exp_t **exp;      /*!<  array of exponents*/
  hash_t *rand;     /*!<  array of random values for each variable
                          to generate hash values out of exponents*/
  deg_t *deg;       /*!<  degree of monmial, for faster sorting and searching*/
  nelts_t *div;     /*!<  latest element from gb checked for its leading
                          term dividing corresponding monomial in hash table*/
  hash_t *idx;      /*!<  index used for matrix generation and marking
                          monomials already taken care of in symbolic
                          preprocessing*/
} mp_cf4_ht_t;



// global meta_data
extern info_t *meta_data;
#endif /* GB_TYPES_H */
