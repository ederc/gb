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
typedef int32_t cf_t;     /* coefficient type */
typedef int32_t val_t;    /* core values like hashes */
typedef val_t hl_t;       /* length of hash table */
typedef hl_t dt_t;        /* data type for other polynomial informatio */
/* like exponent hashes, etc. */
typedef int32_t ind_t;    /* index in hash table structure */
typedef int32_t sdm_t;    /* short divmask for faster divisibility checks */
typedef int32_t len_t;    /* length type for different structures */
typedef int16_t exp_t;    /* exponent type */
typedef int32_t deg_t;    /* (total) degree of polynomial */
typedef len_t bi_t;     /* basis index of element */
typedef len_t bl_t;     /* basis load */
typedef len_t pl_t;     /* pair set load */

/* hash data structure */
typedef struct hd_t hd_t;
struct hd_t
{
    sdm_t sdm;
    deg_t deg;
    len_t div;
    ind_t idx;
    val_t val;
};
/* basis hash table data */
static hl_t *hmap   = NULL; /* global hash map */
static hl_t hsz     = 0;
static exp_t **ev   = NULL; /* exponents from global hash table */
static hd_t *hd     = NULL;
static hl_t eld     = 0;
static hl_t esz     = 0;
static int32_t rht  = 0;  /* basis hash table resetting */

/* update hash table data */
static hl_t *humap  = NULL; /* local hash map */
static hl_t husz    = 0;
static exp_t **evu  = NULL; /* exponents from local hash table */
static hd_t *hdu    = NULL;
static hl_t euld    = 0;
static hl_t eusz    = 0;

/* symbolic hash table data */
static hl_t *hmaps  = NULL; /* local hash map */
static hl_t hssz    = 0;
static exp_t **evs  = NULL; /* exponents from local hash table */
static hd_t *hds    = NULL;
static hl_t esld    = 0;
static hl_t essz    = 0;

static len_t htes   = 0;  /* hash table exponent at start */
static len_t nvars  = 0; /* number of variables */
static len_t bpv    = 0; /* bits per variable in divmask */
static len_t ndvars = 0; /* number of variables for divmask */

/* statistic stuff */
static int il       = 0; /* info level for printed output */
typedef struct stat_t stat_t;
struct stat_t
{
    double round_ctime;
    double select_ctime;
    double symbol_ctime;
    double la_ctime;
    double update_ctime;
    double convert_ctime;
    double overall_ctime;
    double rht_ctime;

    double round_rtime;
    double select_rtime;
    double symbol_rtime;
    double la_rtime;
    double update_rtime;
    double convert_rtime;
    double overall_rtime;
    double rht_rtime;

    int64_t num_pairsred;
    int64_t num_gb_crit;
    int64_t num_redundant;
    int64_t num_rowsred;
    int64_t num_zerored;

    int64_t max_ht_size;
    int64_t len_output;
    int32_t size_basis;
};

/* random values for generating hash values */
static val_t *rv  = NULL;

/* divisor map for short divisibility tests */
static sdm_t *dm  = NULL;

/* pseudo random number generator for hash value
 * generation */
uint32_t rseed  = 2463534242;

/* temporary exponent vector for diffrerent situations */
exp_t *etmp     = NULL;

/* S-pair types */
typedef enum {S_PAIR, GCD_PAIR, GEN_PAIR} spt_t;
typedef struct spair_t spair_t;
struct spair_t
{
    hl_t lcm;
    bi_t gen1;
    bi_t gen2;
    spt_t type;
};

typedef struct ps_t ps_t;
struct ps_t
{
    len_t ld;
    len_t sz;
    len_t mnsel;  /* maximal number of pairs selected */
    spair_t *p;
};

/* monomial order - until we have general orders:
 * 0 = DRL
 * 1 = LEX */
static val_t mo = 0;

/* field characteristic */
static cf_t fc = 0;

/* number generators */
static int32_t ngens  = 0;

/* maximum number of spair selection */
static int32_t mnsel = 0;
/* basis data */

/* we need to store:
 * idx of coefficient array
 * length of array
 * offset for loop unrolling
 * => all in all length = real length + 3
 *
 * how to achieve the same offset, i.e. 3 in the coefficient array?
 * redundant yes/no ?
 * length ?
 * offset ? */

static cf_t **gbcf  = NULL;
static cf_t **tmpcf = NULL;
static dt_t **gbdt  = NULL;
static int8_t *red  = NULL;
static bl_t blold   = 0;
static bl_t bload   = 0;
static bl_t bsize   = 0;

/* lead monomials of all basis elements */
static sdm_t *lms = NULL;

static const long bred  = (long)1;  /* maRking redundant elements */
static const long bmask = ~(long)1; /* maSking redundant elements */

/* matrix data */
static len_t nrall  = 0; /* allocated rows for matrix */
static len_t npivs  = 0; /* new pivots in the current round */
static len_t nrows  = 0; /* rows used in the current round */
static len_t ncols  = 0; /* columns used in the current round */
static len_t nru    = 0; /* number of upper rows (in ABCD splicing) */
static len_t nrl    = 0; /* number of lower rows (in ABCD splicing) */
static len_t ncl    = 0; /* number of lefthand columns(in ABCD splicing) */
static len_t ncr    = 0; /* number of righthand columns(in ABCD splicing) */

/* loop unrolling in sparse linear algebra:
 * we store the offset of the first elements not unrolled
 * in the second entry of the sparse row resp. sparse polynomial. */
#define UNROLL  4

/* threads data */
static int32_t nthrds = 1; /* number of CPU threads */

/* linear algebra options */
static int32_t laopt  = 0;

/* function pointers */
int (*matrix_row_initial_input_cmp)(
        const void *a,
        const void *b
        );

int (*monomial_cmp)(
        const hl_t a,
        const hl_t b
        );

int (*monomial_update_cmp)(
        const hl_t a,
        const hl_t b
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
dt_t **(*linear_algebra)(
        dt_t **mat,
        stat_t *st
        );

cf_t *(*reduce_dense_row_by_known_pivots)(
        int64_t *dr,
        dt_t *const *pivs,
        const hl_t dpiv
        );

dt_t *(*reduce_dense_row_by_known_pivots_sparse)(
        int64_t *dr,
        dt_t *const *pivs,
        const hl_t dpiv,
        const dt_t tmp_pos
        );

cf_t *(*reduce_dense_row_by_all_pivots)(
        int64_t *dr,
        len_t *pc,
        dt_t *const *pivs,
        cf_t *const *dpivs
        );

cf_t *(*reduce_dense_row_by_dense_new_pivots)(
        int64_t *dr,
        len_t *pc,
        cf_t *const *pivs
        );

/* -----------------------------------
 * non-static functions and procedures
 * ----------------------------------- */
int64_t f4_julia(
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
        const int32_t la_option,
        const int32_t info_level
        );
#endif
