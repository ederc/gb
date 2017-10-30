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

#ifndef HASH_CHECK
#define HASH_CHECK 0
#endif

#include <sys/time.h>
#include <stdint.h>

/* 
 * GBLA TYPES
 *  */

#if 0
#define GB_USE_FLOAT XXX
/* #define GB_USE_INT16 XXX */
#else
#define GBLA_USE_UINT16 OK
/* #define GB_USE_UINT32 OK */
#endif
/* #define GB_USE_UINT32 OK */
/* #define GB_USE_INT32 */
/* #define GB_USE_AVX */

/** index type */
typedef uint32_t  mli_t;
/** block index type */
typedef uint32_t  bi_t;

/** storage type for entries */
typedef uint16_t  re_s;
/** storage type for mod */
typedef uint32_t mod_s;

#ifdef GBLA_USE_FLOAT
/** matrix row entry type */
typedef float re_t;
/** matrix row entry type enlarged for delayed modulus */
typedef double re_l_t;
/** matrix row entry type enlarged (half) for delayed modulus */
typedef double re_m_t;
/** type of field characteristic */
typedef float mod_t;
#endif

#ifdef GBLA_USE_DOUBLE
/** matrix row entry type */
typedef double re_t;
/** matrix row entry type enlarged for delayed modulus */
typedef double re_l_t;
/** matrix row entry type enlarged (half) for delayed modulus */
typedef double re_m_t;
/** type of field characteristic */
typedef double mod_t;
#endif
#ifdef GBLA_USE_UINT16
/** matrix row entry type */
typedef uint16_t  re_t;
/// matrix row entry type enlarged (half) for delayed modulus
typedef uint32_t  re_m_t;
/// matrix row entry type enlarged (large) for delayed modulus
typedef uint64_t  re_l_t;
/// type of field characteristic
typedef uint32_t  mod_t;
#endif
#ifdef GBLA_USE_INT16
/** matrix row entry type */
typedef int16_t  re_t;
/** matrix row entry type enlarged for delayed modulus */
typedef int64_t  re_l_t;
/** matrix row entry type enlarged (half) for delayed modulus */
typedef int32_t  re_m_t;
/** type of field characteristic */
typedef int32_t  mod_t;
#endif
#ifdef GBLA_USE_UINT32
/** matrix row entry type */
typedef uint32_t re_t;
/** matrix row entry type enlarged for delayed modulus */
typedef __uint128_t re_l_t;
/** matrix row entry type enlarged (half) for delayed modulus */
typedef uint64_t re_m_t;
/** type of field characteristic */
typedef uint64_t mod_t;
#endif
#ifdef GBLA_USE_INT32
/** matrix row entry type */
typedef int32_t re_t;
/** matrix row entry type enlarged for delayed modulus */
typedef int64_t re_l_t;
/** matrix row entry type enlarged (half) for delayed modulus */
typedef int64_t re_m_t;
/** type of field characteristic */
typedef int64_t mod_t;
#endif

/** row and column index types */
typedef uint32_t  ci_t;
typedef uint32_t  ri_t;
/** number of nonzero elements type */
typedef uint64_t  nnz_t;


#define ALIGNT 32

/* field_ops.h */
#ifdef GBLA_USE_DOUBLE
#define MODP(a,b) \
	fmod((a),(b))
	/* static double MODP(double a, double b) { assert(a>=0) ; assert(b>0) ; double c = fmod(a,b) ; assert(c >=0) ; return c; } */
#define CAST(a) \
	(double) (a)
#endif

#ifdef GBLA_USE_FLOAT
#define MODP(a,b) \
	fmod((a),(b))
	/* static double MODP(double a, double b) { assert(a>=0) ; assert(b>0) ; double c = fmod(a,b) ; assert(c >=0) ; return c; } */
#define CAST(a) \
	(double) (a)
#endif

#ifdef GBLA_USE_UINT16
#define MODP(a,b) \
	(a) % (b)
#define CAST(a) \
	(a) & 0x000000000000ffff
#endif

#ifdef GBLA_USE_INT16
#define MODP(a,b) \
	(a) % (b)
#define CAST(a) \
	(a)
#endif


#ifdef GBLA_USE_UINT32
#define MODP(a,b) \
	(a) % (b)
#define CAST(a) \
	(a) & 0x00000000ffffffff
#endif

#ifdef GBLA_USE_INT32
#define MODP(a,b) \
	(a) % (b)
#define CAST(a) \
	(a) & 0x00000000ffffffff
#endif

/**
 * \brief Sparse matrix structure for reading jcf matrices
 */

typedef struct sm_t {
  ri_t nrows;     /*!<  number of rows */
  ci_t ncols;     /*!<  number of columns */
  nnz_t nnz;      /*!<  number of nonzero entries */
  float density;  /*!<  density used for adjusting memory allocation during
                        splicing and generation of ABCD blocks */
  mod_t mod;      /*!<  modulo/field characteristic */
  float fs;       /*!<  file size of input matrix */
  char *fsu;      /*!<  file size unit of input matrix, e.g. GB */
  re_t **rows;    /*!<  address of row: M->rows[i] gives first
                        address of ith row */
  ci_t **pos;     /*!<  position of entry in row: M->pos[i] gives first
                        address of first position of nonzero entry in row i */
  ci_t *rwidth;   /*!<  width of row: M->rwidth[i] gives number of nonzero
                        entries in row i */
  ci_t *buf;      /*!<  stores buffer of memory allocated for the given row*/
} sm_t;

/**
 * \brief Data type to store rows of dense row matrices: We keep the entries
 * resp. values in an array val of fixed size ncols. Moreover, we store the
 * current column position of the first nonzero entry in the corresponding row
 * in lead. The idea is the following:
 *
 * 1. The initial re_t values are stored in init_val with lead entry lead.
 * 2. During the reducion of D we copy these initial rows to val, which is of
 *    type re_l_t, so we can keep on reducing without doing modulo operations
 *    until the end. They always update lead as lead entry of val.
 * 3. Once we have finished a row, we have a new pivot row. This one is
 *    normalized modulo the field characteristic and then stored in piv_val with
 *    lead entry piv_lead.
 */
typedef struct dr_t {
  re_l_t *val;  /*!< entries in row of dense row matrix */
  re_t *init_val; /*!< entries in row of dense row matrix, once pivots are fixed */
  re_t *piv_val;  /*!< entries in row of dense row matrix, once pivots are fixed */
  ci_t lead;      /*!< lead entry of the corresponding row */
  ci_t piv_lead;  /*!< lead entry of the corresponding pivot row */
} dr_t;

/**
 * \brief Dense row matrix structure for reducing D. We directly store data in
 * size re_l_t in order to not copy to wide arrays all the time. Usually D is
 * very small compared to the whole matrix, so memory should not be a problem in
 * general.
 */
typedef struct dm_t {
  ri_t nrows;     /*!<  number of rows */
  ci_t ncols;     /*!<  number of columns */
  ri_t rank;      /*!<  rank of matrix */
  nnz_t nnz;      /*!<  number of nonzero entries */
  float density;  /*!<  density used for adjusting memory allocation during
                        splicing and generation of ABCD blocks */
  mod_t mod;      /*!<  modulo/field characteristic */
  dr_t **row;     /*!<  rows resp. elements of matrix */
} dm_t;

/**
 * \brief A sparse matrix block
 */
typedef struct sbl_t {
  re_t **val; /*!< row or column entries */
  bi_t **pos; /*!< position in row */
  bi_t *sz;   /*!< size of row */
  bi_t *buf;  /*!< memory buffer already allocated */
} sbl_t;

/**
 * \brief Sparse block matrix structure for Faugère-Lachartre decompositions.
 * For non multiline implementation using small sparse blocks for the A part.
 */

typedef struct sb_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  sbl_t **blocks;   /*!<  address of blocks: M->blocks[i][j] gives address of
                          block. */
} sb_fl_t;

/**
 * \brief Sparse matrix structure for Faugère-Lachartre decompositions.
 * For non multiline implementation using small sparse blocks for the A part.
 */

typedef struct sm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  re_t **row;       /*!< row entries */
  ci_t **pos;       /*!< position in row */
  ci_t *sz;         /*!< size of row */
  ci_t *buf;        /*!< memory buffer already allocated */
} sm_fl_t;

/**
 * \brief An intermediate block representation that can be sparse or dense.
 */
typedef struct ibl_t {
  re_t *val;  /*!< rows, NULL if sparse */
  re_t **row; /*!< entries in sparse row, NULL if dense */
  bi_t **pos; /*!< positions, NULL if dense */
  bi_t *sz;   /*!< length of sparse rows, NULL if dense*/
} ibl_t;

/**
 * \brief A dense block is a block of size  __GBLA_SIMD_BLOCK_SIZE^2 of
 * matrix entries.
 */
typedef struct dbl_t {
  re_t *val;  /*!< rows */
} dbl_t;

/**
 * \brief A multiline block is a vector of an vector of multilines. It consists of
 * an index for the column of the value entries in rows. For each index rows stores
 * __GBLA_NROWS_MULTILINE elements, i.e. this many rows are taken care of at once.
 */
typedef struct mbl_t {
  bi_t *idx;  /*!< column index in the multiline vector */
  re_t *val;  /*!< multiline row, must be __GBLA_NROWS_MULTILINE * length(idx) */
  bi_t sz;    /*!< current length of the block row */
  bi_t dense; /*!< if 1 the multiline row is in dense representation */
} mbl_t;

/**
 * \brief A multiline is a vector of an vector of multilines. It consists of
 * an index for the column of the value entries in rows. For each index rows stores
 * __GBLA_NROWS_MULTILINE elements, i.e. this many rows are taken care of at once.
 *
 * \note In spite of type mbl_t the index idx as well as the size sz might
 * excced 2^16, thus bi_t is not enough and we need them to be of type ci_t
 * resp. uint32_t at least.
 *
 */
typedef struct ml_t {
  mli_t *idx; /*!< column index in the multiline vector */
  re_t  *val; /*!< multiline row, must be __GBLA_NROWS_MULTILINE * length(idx) */
  ci_t  sz;    /*!< current length of the multiline row */
  bi_t  dense; /*!< if 1 the multiline row is in dense representation */
} ml_t;


/**
 * \brief Enum of block alignments: Entries in the block submatrices are stored
 * w.r.t. different alignments, some from top to down, then from left to right,
 * other in combinations and inversions of these.
 */
enum ba_t {
  tdlr, /*!<  top-to-down, left-to-right */
  tdrl, /*!<  top-to-down, right-to-left */
  dtlr, /*!<  down-to-top, left-to-right */
  dtrl  /*!<  down-to-top, right-to-left */
};


/**
 * \brief Mixed block matrix structure for Faugère-Lachartre decompositions.
 * For non multiline implementation using small dense and sparse blocks for
 * exploiting SIMD inctructions.
 */

typedef struct ibm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  ibl_t **blocks;   /*!<  address of blocks: M->blocks[i][j] gives address of
                          block. */
} ibm_fl_t;

/**
 * \brief Dense block matrix structure for Faugère-Lachartre decompositions.
 * For non multiline implementation using small dense blocks for exploiting SIMD
 * inctructions.
 */

typedef struct dbm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  dbl_t **blocks;   /*!<  address of blocks: M->blocks[i][j] gives address of
                          block. */
} dbm_fl_t;

/**
 * \brief Hybrid simd block matrix structure for Faugère-Lachartre decompositions.
 * For non multiline implementation using small inner dense blocks for exploiting
 * SIMD inctructions.
 *
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * |       |       |       |       |
 * ---------------------------------
 * 
 * The above picture represents one block of size __GBLA_SIMD_BLOCK_SIZE_RECT.
 * Inside we store small sub blocks of height 1 and width __GBLA_SIMD_INNER_SIZE.
 * Such an inner sub block is either NULL (if all entries are zero) or
 * represented in dense fashion.
 * All in all this leads to the fact that blocks is of type dbl_t ****blocks:
 * blocks[outer_row][outer_col][inner_row][inner_col].
 */
typedef struct hbm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  dbl_t ****blocks; /*!<  address of blocks: M->blocks[i][j][k][l] as explained
                          above */
} hbm_fl_t;

/**
 * \brief Sparse block matrix structure for Faugère-Lachartre decompositions.
 * Can be used for usual line implementations and multi line implementations,
 * the corresponding multi line functions are labeled with an "_ml".
 *
 */

typedef struct sbm_fl_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  ri_t bheight;     /*!<  number of rows per block */
  ci_t bwidth;      /*!<  number of columns per block */
  enum ba_t ba;     /*!<  memory alignment in the block: depending on
                          the block parts A,B,C and D the entries are
                          stored top-to-down, left-to-right, or any
                          combination and inversion of these two */
  int fe;           /*!<  if 1 then empty blocks are used to fill
                          block matrix
                          if 0 then no fill is done */
  int hr;           /*!<  if 0 then those rows are not accepted
                          if 1 then hybrid (sparse/dense) rows are accepted */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  mbl_t ***blocks;  /*!<  address of blocks: M->blocks[i][j] gives address of
                          block. There are nrows/bheight * ncols/bwidth blocks. */
} sbm_fl_t;

/**
 * \brief Sparse matrix with multiline row structure for Faugère-Lachartre
 * decompositions. Can be used for multiline implementations, the corresponding
 * multiline functions are labeled with an "_ml".
 */

typedef struct sm_fl_ml_t {
  ri_t nrows;       /*!<  number of rows */
  ci_t ncols;       /*!<  number of columns */
  enum ba_t ba;     /*!<  memory alignment in the block: depending on
                          the block parts A,B,C and D the entries are
                          stored top-to-down, left-to-right, or any
                          combination and inversion of these two */
  int fe;           /*!<  if 1 then empty blocks are used to fill
                          block matrix
                          if 0 then no fill is done */
  int hr;           /*!<  if 1 then hybrid (sparse/dense) rows are accepted
                          if 0 then those rows are not accepted */
  nnz_t nnz;        /*!<  number of nonzero elements */
  double density;   /*!<  density of this submatrix */
  ml_t *ml;         /*!<  address of multilines: M->ml[i] gives address of
                          multiline i. There are nrows/__GBLA_NROWS_MULTILINE
                          multilines. */
} sm_fl_ml_t;

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



/* number of variables */
/* typedef uint16_t nvars_t; */
typedef uint32_t nvars_t;

/* number of elements */
typedef uint32_t nelts_t;

/* hash table table size */
typedef uint32_t ht_size_t;

/* hash table entry size */
typedef int32_t hash_t;
/* typedef uint64_t hash_t; */
/* typedef unsigned long hash_t; */

/* degree size */
/* typedef uint16_t deg_t; */
/* typedef unsigned int deg_t; */
typedef int32_t deg_t;

/* homogeneity */
typedef int hom_t;

/* monomial ordering */
typedef int ord_t;

/* divmask */
typedef int divm_t;

/* exponent size */
/* typedef uint8_t exp_s; */
typedef uint16_t exp_s;
typedef exp_s exp_t;

/* typedef re_t cf_t; */
typedef uint32_t cf_t;
typedef re_l_t bf_t;

/* block size */
typedef uint16_t bs_t;

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
