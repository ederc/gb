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
 * \file data.h
 * \brief General and global data
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_DATA_H
#define GB_DATA_H

#include <stdint.h>

/* meta data about computation */


/* computational data */
typedef int32_t hv_t;   /* hash values */
typedef int32_t deg_t;  /* (total) degree of polynomial */
typedef int32_t bi_t;   /* basis index of element */
typedef int32_t bl_t;   /* basis load */
typedef int32_t pl_t;   /* pair set load */
typedef int32_t *poly_t;  /* polynomials are just arrays of coeffs
                             and exponent hashes. They start with the
                             length of the array and the second entry
                             is used for distinguishing between known
                             basis elements and new basis elements when
                             applying the linear algebra. That means
                             p = [len, known, cf1, eh1, cf2, eh2, ...] */


/* S-pair definitions */

/* S-pair types */
typedef enum {S_PAIR, GCD_PAIR, GEN_PAIR} spt_t;
tyoedef struct spair_t spair_t;
struct spair_T
{
  deg_t deg;   /* if criteria apply, information is stored here */
  spt_t type;
  hash_t lcm;
  bi_t gen1;
  bi_t gen2;
};

/* pair set data */
static spair *ps  = NULL;
static pl_t pload = 0;
static pl_t psize = 0;

/* basis data */
static poly_t *gb   = NULL;
static bl_t bload   = 0;
static bl_t bsize   = 0;

static void initialize();

/******************************************
 * general data structures needed in gb
 *****************************************/

/**
 * \brief Verbosity information
 */
struct info_t
{
  unsigned long long non_div;
  unsigned long long non_div_found;
  nelts_t sel_pairs;                          /*!<  number of selected pairs in last step*/
  nelts_t curr_deg;                           /*!<  degree of pairs in current step*/
  nelts_t mat_cols;                           /*!<  number of columns in last gbla matrix*/
  nelts_t mat_rows;                           /*!<  number of rows in last gbla matrix*/
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
};
typedef struct info_t info_t;


/**
 * \brief Typedef of enum for redundancy: if a new element is added to the
 * groebner basis which lead monomial divides the lead monomial of an existing
 * element, this element is redundant.
 */
typedef hash_t red_t;

/**
 * \brief Data structure storing temporarily a  list of indices of basis
 * elements that already appear in the symbolic preprocessing as a generator of
 * an spair with the given lcm. Thus we can remove duplicate rows from the gbla
 * matrices resp. not enter them to the matrix construction at all.
 */
typedef struct dup_t
{
  nelts_t *idx;
  nelts_t load;
  nelts_t size;
  hash_t lcm;
} dup_t;

typedef uint32_t poly_t;

/**
 * \brief Groebner basis, starts with input data
 * 
 * \note For easier divisibility checks in symbolic preprocessing we put at
 * the first position of basis, i.e. index 0 a NULL element.
 * Thus, basis->load is always one bigger than the actual number of elements
 * in the basis.
 *
 * \note The first elements are just the input data. Those elements are not a
 * priori part of the basis. The basis starts at index st.
 */
typedef struct gb_t
{
  /* global data */
  nelts_t size;     /*!<  memory allocated */
  nelts_t load;     /*!<  number of elements in basis*/
  nelts_t load_ls;  /*!<  number of elements in basis after last step*/
  nelts_t fl;       /*!<  final length of basis without redundant elements*/
  nelts_t st;       /*!<  start of the real basis, everything before is input data */
  nvars_t rnv;      /*!<  real number of variables from input */
  nvars_t nv;       /*!<  number of variables, possibly including homogenization variable */
  hom_t init_hom;   /*!<  homogeneous input? 1=yes, 0=no */
  hom_t hom;        /*!<  homogeneous computation? 1=yes, 0=no (we might homogenize) */
  ord_t ord;        /*!<  monomial ordering */
  nelts_t nred;     /*!<  number of redundant elements in basis */
  cf_t mod;         /*!<  modulo/field characteristic */
  int has_unit;     /*!<  is set to 1 if we have found a unit in the basis */
  nelts_t max_sel;  /*!<  maximal number of spairs handled at once */
  /* element data */
  deg_t *deg;       /*!<  degree of each element resp. polynomial*/
  red_t *red;       /*!<  stores if the element is redundant or not*/
  poly_t **p;       /*!<  polynomial data, i.e. exponent hashes and coefficients */
  /* meta data */
  char **vnames;    /*!<  variable names */
  size_t mtl;       /*!<  maximal length of term (needed for
                          estimating buffer size when writing result)*/
  double fs;        /*!<  file size of input matrix */
  const char *fsu;  /*!<  file size unit of input matrix, e.g. GB */
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
  hash_t lcm;       /*!<  hash of lcm of the lead terms of gen1 and gen2 */
  deg_t deg;        /*!<  degree of S-pair */
  nelts_t gen1;     /*!<  index to first generator */
  nelts_t gen2;     /*!<  index to second generator */
  criteria_t crit;  /*!<  tracker if product and chain criterion applies */
} spair_t;

/**
 * \brief Pair set
 */
typedef struct ps_t
{
  /* global data */
  nelts_t size;     /*!<  memory allocated */
  nelts_t load;     /*!<  number of elements in basis*/
  /* element data */
  spair_t *pairs;   /*!<  pointers of spairs */
} ps_t;

/**
 * \brief Selection of multiplied elements for next matrix: Those are generators
 * of spairs and corresponding reducers.
 */
typedef struct sel_t
{
  /* global data */
  deg_t deg;        /*!<  maximal degree of all elements in selection set*/
  nelts_t size;     /*!<  memory allocated */
  nelts_t load;     /*!<  number of elements in selection*/
} sel_t;

/**
 * \brief Data structure storing the list of positions in hash table
 * corresponding to monomials. Used for symbolic preprocessing to get all
 * possible reducers.
 */
typedef struct pre_t
{
  hash_t *hash;   /*!<  position of monomials in hash table*/
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
  sb_fl_t *A;   /*!<  upper left sparse block matrix part */
  sm_fl_t *AR;  /*!<  upper left sparse row matrix part */
  dbm_fl_t *B;  /*!<  upper right dense block matrix part */
  dm_t *BR;     /*!<  reduced B in dense row matrix format */
  sb_fl_t *C;   /*!<  lower right sparse block matrix part */
  sm_fl_t *CR;  /*!<  lower right sparse row matrix part */
  dbm_fl_t *D;  /*!<  lower right dense block matrix part */
  dm_t *DR;     /*!<  reduced D in dense row matrix format */
  cf_t mod;     /*!<  modulo/field characteristic */
  bi_t bs;      /*!<  block size given by gbla */
  ci_t ncl;     /*!<  number of columns lefthand side, i.e. cols of A resp. C */
  ci_t ncr;     /*!<  number of columns righthand side, i.e. cols of B resp. D */
  ri_t nru;     /*!<  number of rows upper par, i.e. rows of A resp. B*/
  ri_t nrl;     /*!<  number of rows lower par, i.e. rows of C resp. D*/
  ri_t rbu;     /*!<  number of row blocks for upper part of gbla matrix */
  ri_t rbl;     /*!<  number of row blocks for lower part of gbla matrix */
  ci_t cbl;     /*!<  number of column blocks for left part of gbla matrix */
  ci_t cbr;     /*!<  number of column blocks for right part of gbla matrix */
} mat_t;

/**
 * \brief sparse row structure
 */
typedef struct sr_t
{
  cf_t *val;
  nelts_t *pos;
  nelts_t sz;
} sr_t;

typedef uint32_t src_t;
/**
 * \brief sparse matrix structure
 *
 * \note the rows r have the following structure:
 *  r[0] = number of nonzero entries in row
 *  ----
 *  r[1] = column position of r[2]
 *  r[2] = entry
 *  r[3] = column position of r[4]
 *  r[4] = entry
 *  ...
 */
typedef struct smc_t
{
  src_t **row;  /*!<  rows of matrix, stores not only entries, but also positions */
  ci_t ncl;     /*!<  number of columns lefthand side, i.e. cols of A resp. C */
  ci_t ncr;     /*!<  number of columns righthand side, i.e. cols of B resp. D */
  ri_t nr;      /*!<  number of rows lower par, i.e. rows of C resp. D*/
  ri_t rk;      /*<!  rank of the matrix */
  cf_t mod;     /*!<  modulo/field characteristic */
} smc_t;

/**
 * \brief sparse matrix structure
 */
typedef struct smat_t
{
  sr_t **row;     /*!<  rows of the matrix */
  ci_t ncl;       /*!<  number of columns lefthand side, i.e. cols of A resp. C */
  ci_t ncr;       /*!<  number of columns righthand side, i.e. cols of B resp. D */
  ri_t nr;        /*!<  number of rows*/
  ri_t rk;        /*<!  rank of the matrix */
  cf_t mod;       /*!<  modulo/field characteristic */
} smat_t;

/**
 * len stores the lengths of each row: it starts with 0 and adds up how many
 * entries are in the row, i.e. len[i] = starting position in pos and val for
 * row i+1. len[bs] = # elements in block
 */
typedef struct mat_gb_block_t
{
  nelts_t *len;
  bs_t *pos;
  cf_t *val;
  nelts_t nr;
} mat_gb_block_t;

typedef struct mat_gb_meta_data_t
{
  cf_t mod;
  nelts_t bs;
  nelts_t rk;
  nelts_t nc_AC;
  nelts_t nc_BD;
  nelts_t nc;
  nelts_t ncb;
  nelts_t ncb_AC;
  nelts_t ncb_BD;
  nelts_t nr_AB;
  nelts_t nr_CD;
  nelts_t nr;
  nelts_t nrb_AB;
  nelts_t nrb_CD;
} mat_gb_meta_data_t;

/**
 * \brief Struct keeping all function pointers of functions depending on the
 * chosen monomial order.
 */
typedef struct sort_t
{
  nelts_t (*get_pairs_by_minimal_degree)(ps_t *ps);
  void (*sort_presorted_columns_invert_left_side)(pre_t *mon, const int nthreads);
  void (*sort_presorted_columns)(pre_t *mon, const int nthreads);
  void (*sort_rows_by_decreasing_lm)(smc_t *M);
  void (*sort_rows_by_increasing_lm)(smc_t *M);
  void (*sort_columns)(pre_t *mon);
  int (*compare_spairs)(const void *a, const void *b);
  int (*compare_monomials)(const void *a, const void *b);
  int (*compare_monomials_inverse)(const void *a, const void *b);
  int (*compare_polynomials)(const void *a, const void *b);
  int (*compare_polynomials_inverse)(const void *a, const void *b);
} sort_t;

/**
 * \brief Hash table using linear probing combined with Robin Hood hashing
 * (cf. https://cs.uwaterloo.ca/research/tr/1986/CS-86-14.pdf, Robin Hood
 * Hashing by Pedro Celis, 1986
 *
 * \note We start with the first element at index 1, not at index 0. Thus we can
 * optimize divisibility checks (if index 0 is returned we have not found any
 * divisor). See also procedure init_hash_table().
 * 
 * \note The lookup table part "lut" is important: If we store the hashes
 * directly, all pairs, all symbolic preprocessing data, etc. needs to be
 * recomputed when we enlarge the hash table during the computation.
 */
typedef struct ht_t
{
  /* size and load counters */
  nvars_t nv;         /*!<  number of variables*/
  nvars_t offset;     /*!<  offset for loop starts: ht->nv == even => 0 else 1 */
  ht_size_t load_ls;  /*!<  load of hash table after last step */
  ht_size_t load;     /*!<  load of hash table*/
  ht_size_t sz;       /*!<  hash table size*/
  /* data and arrays for storage */
  ht_size_t *lut;     /*!<  lookup table between hash value and position in
                            exponent array*/
  hash_t *val;        /*!<  array of hash values*/
  exp_t *exp;        /*!<  array of exponents, note that exp_t is possibly
                            SSE/AVX vector if available*/
  hash_t *rand;       /*!<  array of random values for each variable
                            to generate hash values out of exponents*/
  deg_t *deg;         /*!<  degree of monmial, for faster sorting and searching*/
  nelts_t *div;       /*!<  latest element from gb checked for its leading
                            term dividing corresponding monomial in hash table*/
  nelts_t bpv;        /*!<  bits per variable in divmask */
  nelts_t ndv;        /*!<  number of divmask variables */
  deg_t *divmap;      /*!<  divmap for each divmask caluclations */ 
  divm_t *dm;         /*!<  divmask for monomial */
  uint32_t rcdm;      /*!<  counter when to recalculate divmaps and divmasks */
  uint32_t muldm;     /*!<  multiplier for divmap recalculation range */
#if HASH_CHECK
  ht_size_t *ctr;     [>!<  counts how often a hash appears <]
#endif
  ht_size_t *idx;     /*!<  index used for matrix generation and marking
                            monomials already taken care of in symbolic
                            preprocessing*/
  sort_t sort;        /*!<  sort functions pointers */
} ht_t;

/**
 * \brief Dense row structure for buffering nonzero elements when generating
 * gbla matrix from symbolic preprocessing data.
 *
 * \note Data type of ctr needs to be nelts_t > bi_t since otherwise we would
 * increment the counter in a dense row to 256 = 0 mod sizeof(bi_t) (=256).
 */
typedef struct dbr_t
{
  nelts_t *ctr; /*!< array of counts how many nonzero elements are stored */
  re_t **cf; /*!< array of coefficients in the given row and block */
  re_t **bl; /*!< array of blocks of coefficients for the dense part of
                     the gbla matrix*/
} dbr_t;

typedef struct src_tmp_t
{
  uint32_t pos;
  uint32_t cf;
} src_tmp_t;


/* global meta_data */
extern info_t *meta_data;

#endif /* GB_TYPES_H */
