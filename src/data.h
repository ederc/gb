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
typedef int16_t cf16_t;   /* coefficient type */
typedef int32_t cf32_t;   /* coefficient type */
typedef int32_t val_t;    /* core values like hashes */
typedef val_t hl_t;       /* length of hash table */
typedef hl_t dt_t;        /* data type for other polynomial informatio */
typedef int8_t red_t;     /* redundancy type */
/* like exponent hashes, etc. */
typedef int32_t ind_t;    /* index in hash table structure */
typedef int32_t sdm_t;    /* short divmask for faster divisibility checks */
typedef int32_t len_t;    /* length type for different structures */
typedef int16_t exp_t;    /* exponent type */
typedef int32_t deg_t;    /* (total) degree of polynomial */
typedef len_t bi_t;     /* basis index of element */
typedef len_t bl_t;     /* basis load */
typedef len_t pl_t;     /* pair set load */

/* hash table data */

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

typedef struct ht_t ht_t;
struct ht_t
{
    hl_t hsz;         /* size of hash data array */
    hl_t eld;         /* load of exponent vector */
    hl_t esz;         /* allocated exponent vector size */
    len_t nv;         /* number of variables */
    len_t ndv;        /* number of variables used for divmask */
    len_t bpv;        /* bits per variable for divmask */
    uint32_t rseed;   /* random seed */
    hl_t *map;        /* map between hash data and exponent vector */
    hd_t *hd;         /* hash data array */
    val_t *rv;        /* randomizing array for hashing */
    sdm_t *dm;        /* short divisor mask */
    exp_t **ev;       /* exponent vector array */
};

/* statistic stuff */
typedef struct stat_t stat_t;
struct stat_t
{
    int32_t info_level;   /* printout level */

    /* global computational data */
    int32_t mon_order;    /* monomial order */
    int32_t field_char;   /* field characteristic */
    int32_t la_variant;   /* linear algebra variant */
    int32_t regen_ht;     /* regenerating the global hash table */
    int32_t nthrds;       /* number of threads */
    int32_t nr_vars;      /* number of variables */
    int32_t max_nr_pairs; /* maximal number of pairs per matrix */
    int32_t nr_gens;      /* number of generators of input */
    int32_t init_ht_sz;   /* initial hash table size */

    /* CPU timings for different parts */
    double round_ctime;
    double rght_ctime;
    double select_ctime;  
    double symbol_ctime;
    double la_ctime;
    double update_ctime;
    double convert_ctime;
    double overall_ctime;

    /* real timings for different parts */
    double round_rtime;
    double rght_rtime;
    double select_rtime;
    double symbol_rtime;
    double la_rtime;
    double update_rtime;
    double convert_rtime;
    double overall_rtime;

    /* meta data */
    int64_t num_pairsred; 
    int64_t num_gm_crit;
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

/* basis definition */
typedef struct bs_t bs_t;
struct bs_t
{
    bl_t lo;    /* old load before last update */
    bl_t ld;    /* current load of basis */
    bl_t sz;    /* size of basis allocated */
    void **cf;  /* coefficient arrays of eleuments, type depends on */
                /* underlying ring, e.g. finite field or rationals or ... */
    dt_t **hd;  /* hash data arrays of elements */
    sdm_t *lm;  /* lead monomials of elements */
    red_t *red; /* is the element redundant? */
    void **tcf; /* temporary coefficient arrays for storing coefficents */
                /* during linear algebra computation, type depends on */
                /* underlying ring, see above */
};

/* matrix data */
typedef struct mat_t mat_t;
struct mat_t
{
    len_t na;   /* number of rows allocated */
    len_t nr;   /* number of rows */
    len_t nc;   /* number of columns */
    len_t np;   /* number of pivots */
    len_t nru;  /* number of upper rows (ABCD splicing) */
    len_t nrl;  /* number of lower rows (ABCD splicing) */
    len_t ncl;  /* number of left columns (ABCD splicing) */
    len_t ncr;  /* number of right columns (ABCD splicing) */
    dt_t **r;   /* rows of matrix */
};

/* loop unrolling in sparse linear algebra:
 * we store the offset of the first elements not unrolled
 * in the second entry of the sparse row resp. sparse polynomial. */
#define UNROLL  4

/* function pointers */
int (*matrix_row_initial_input_cmp)(
        const void *a,
        const void *b
        );

int (*monomial_cmp)(
        const hl_t a,
        const hl_t b
        );

int (*monomial_local_cmp)(
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

cf_t *(*reduce_dense_row_by_known_pivots_sparse)(
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
