/* gb: Gröbner Basis
 * Copyright (C) 2018 Christian Eder <ederc@mathematik.uni-kl.de>
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
#include <stdio.h>
#include <string.h> /* for memset et al. */
#include <limits.h>
#include <math.h>

/* check if OpenMP is available */
/* #if defined(_OPENMP) */
#include <omp.h>
/* #else
 * typedef int omp_int_t;
 * inline omp_int_t omp_get_thread_num(void) { return 0;}
 * inline omp_int_t omp_get_max_threads(void) { return 1;}
 * #endif */

#define ORDER_COLUMNS 1

/* computational data */
typedef int32_t sdm_t;    /* short divmask for faster divisibility checks */
typedef int32_t len_t;    /* length type for different structures */
typedef int32_t exp_t;    /* exponent type */
typedef int32_t val_t;    /* core values like hashes */
typedef int32_t cf_q_t;   /* coefficient type for rationals */
typedef int32_t cf_32_t;  /* coefficient type 32 bit primes */
typedef int16_t cf_16_t;  /* coefficient type 16 bit primes */
typedef int32_t deg_t;    /* (total) degree of polynomial */
typedef int32_t bi_t;     /* basis index of element */
typedef int32_t bl_t;     /* basis load */
typedef int32_t pl_t;     /* pair set load */

/* S-pair types */
typedef enum {S_PAIR, GCD_PAIR, GEN_PAIR} spt_t;
typedef struct spair_t spair_t;
struct spair_t
{
  deg_t deg;   /* if criteria apply, information is stored here */
  val_t lcm;
  bi_t gen1;
  bi_t gen2;
  spt_t type;
};

/* pair set data */
static spair_t *ps    = NULL;
static pl_t pload     = 0;
static pl_t psize     = 0;
static int32_t mnsel  = 0; /* maximal number of pairs to be selected */

/* monomial order - until we have general orders:
 * 0 = DRL
 * 1 = LEX */
static val_t mo = 0;

/* field characteristic */
static val_t fc = 0;

typedef struct md_t md_t;
struct md_t
{
  int32_t nv;   /* number of variables */
  int32_t os;   /* offset when looping over number of variables */
  int32_t mo;   /* monomial order: 0 = DRL, 1 = LEX */
  int32_t fc;   /* field characteristic */
  int32_t la;   /* linear algebra option */
  int32_t ng;   /* number of generators */
  int32_t mp;   /* maximal number of pairs selected at once */
  int32_t nt;   /* number of threads */
  int32_t rght; /* reset global hash table */
  double density;
  /* statistics */
  double select_ctime;
  double select_rtime;
  double symbol_ctime;
  double symbol_rtime;
  double update_ctime;
  double update_rtime;
  double convert_ctime;
  double convert_rtime;
  double reduce_ctime;
  double reduce_rtime;
  double rght_ctime;
  double rght_rtime;
  double la_ctime;
  double la_rtime;
  double psort_rtime;
  int64_t num_pairsred;
  int64_t num_gb_crit;
  int64_t num_redundant;
  int64_t num_duplicates;
  int64_t num_rowsred;
  int64_t num_zerored;
  int64_t num_ht_enlarge;
  int64_t num_sdm_found;
  int64_t num_not_sdm_found;
};

/* structured rows resp. polynomials (depending on the context)
 * we do not know what type our coefficients have a priori: we support
 * 16-bit prime fields, 32-bit prime fields and the rationals */
typedef struct row_t row_t;
struct row_t
{
  len_t sz;   /* size of row */
  len_t os;   /* offset for loop unrolling */
  int32_t rd; /* redundant? */
  void *cf;   /* coefficients */
  len_t *ch;  /* column positions resp. hash positions */
};

/* structured matrices */
typedef struct mat_t mat_t;
struct mat_t
{
  len_t nr;   /* number of rows */
  len_t nc;   /* number of columns */
  len_t nru;  /* number of upper rows (ABCD split) */
  len_t nrl;  /* number of lower rows (ABCD split) */
  len_t ncl;  /* number of left columns (ABCD split) */
  len_t ncr;  /* number of right columns (ABCD split) */
  len_t np;   /* number of pivots */
  len_t na;   /* number of rows allocated */
  row_t **r;  /* rows */
};

/* structured basis */
typedef struct bs_t bs_t;
struct bs_t
{
  len_t ld;   /* current load of basis */
  len_t sz;   /* size of basis */
  len_t ol;   /* old load of basis before entering new elements */
  row_t **p;  /* polynomials */
  len_t *lm;  /* lead monomials of basis elements */
};

/* basis data */
static val_t **bs = NULL;
static bl_t blold = 0;
static bl_t bload = 0;
static bl_t bsize = 0;

/* lead monomials of all basis elements */
static val_t *lms = NULL;

static const long bred  = (long)1;  /* maRking redundant elements */
static const long bmask = ~(long)1; /* maSking redundant elements */

/* matrix data */
static int32_t nrall  = 0; /* allocated rows for matrix */
static int32_t npivs  = 0; /* new pivots in the current round */
static int32_t nrows  = 0; /* rows used in the current round */
static int32_t ncols  = 0; /* columns used in the current round */
static int32_t nru    = 0; /* number of upper rows (in ABCD splicing) */
static int32_t nrl    = 0; /* number of lower rows (in ABCD splicing) */
static int32_t ncl    = 0; /* number of lefthand columns(in ABCD splicing) */
static int32_t ncr    = 0; /* number of righthand columns(in ABCD splicing) */

/* loop unrolling in sparse linear algebra:
 * we store the offset of the first elements not unrolled
 * in the second entry of the sparse row resp. sparse polynomial. */
#define UNROLL  4

/* threads data */
static int32_t nthrds = 1; /* number of CPU threads */

/* linear algebra options */
static int32_t laopt  = 0;

/* function pointers */
mat_t *(*import_julia_data)(
    const int32_t * const lens,
    const int32_t * const cfs,
    const int32_t * const exps,
    const md_t * const md
    );

int64_t (*export_julia_data)(
    int32_t **bp,
    const bs_t * const bs,
    const md_t * const md
    );

row_t *(*multiplied_polynomial_to_matrix_row)(
    const len_t hm,
    const deg_t deg,
    const exp_t * const em,
    const row_t * const poly
    );

void (*normalize_matrix_row)(
    row_t *row,
    const md_t * const md
    );

int (*matrix_row_initial_input_cmp)(
    const void *a,
    const void *b
    );

int (*monomial_cmp)(
    const exp_t * const ea,
    const exp_t * const eb
    );

int (*spair_cmp)(
    const void *a,
    const void *b
    );

int (*hcm_cmp)(
    const void *a,
    const void *b
    );

/* linear algebra routines */
void (*linear_algebra)(
    mat_t *mat,
    md_t *md
    );

void (*sparse_linear_algebra)(
    mat_t *mat,
    md_t *md
    );

void (*probabilistic_sparse_linear_algebra)(
    mat_t *mat,
    md_t *md
    );

/* -----------------------------------
 * non-static functions and procedures
 * ----------------------------------- */
int64_t f4_julia_ff(
    int32_t **jl_basis,
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const int32_t field_char,
    const int32_t mon_order,
    const int32_t nr_vars,
    const int32_t nr_gens,
    const int32_t ht_size,
    const int32_t nr_threads,
    const int32_t max_nr_pairs,
    const int32_t reset_hash_table,
    const int32_t la_option
    );
#endif
