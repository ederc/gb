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

#define _GNU_SOURCE
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h> /* for memset et al. */
#include <limits.h>
#include <math.h>
/* #include <threads.h> */

/* check if OpenMP is available */
/* #if defined(_OPENMP) */
#include <omp.h>
/* #else
 * typedef int omp_int_t;
 * inline omp_int_t omp_get_thread_num(void) { return 0;}
 * inline omp_int_t omp_get_max_threads(void) { return 1;}
 * #endif */

#define ORDER_COLUMNS 1
/* loop unrolling in sparse linear algebra:
 * we store the offset of the first elements not unrolled
 * in the second entry of the sparse row resp. sparse polynomial. */
#define UNROLL  4


/* computational data */
typedef int32_t cf32_t;     /* coefficient type */
typedef int32_t val_t;    /* core values like hashes */
typedef val_t hl_t;       /* length of hash table */
typedef hl_t hm_t;        /* hashed monomials for polynomial entries */
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
    val_t val;
    sdm_t sdm;
    deg_t deg;
    ind_t idx;
};

/* hash table data structure */
typedef struct ht_t ht_t;
struct ht_t
{
    exp_t **ev;   /* exponent vector */
    hd_t *hd;     /* hash data */
    hl_t *hmap;   /* hash map */
    hl_t eld;     /* load of exponent vector */
    hl_t esz;     /* size of exponent vector */
    hl_t hsz;     /* size of hash map */
    len_t nv;     /* number of variables */
    sdm_t *dm;    /* divisor map for divisibility checks */
    len_t ndv;    /* number of variables for divmask */
    len_t bpv;    /* bits per variable in divmask */
    val_t *rn;    /* random numbers for hash generation */
    uint32_t rsd; /* seed for random number generator */
};

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

typedef struct bs_t bs_t;
struct bs_t
{
    bl_t ld;        /* load of basis */
    bl_t sz;        /* size allocated for basis */
    bl_t lo;        /* load before current update */
    deg_t mltdeg;   /* maximal appearing degree in lead term in basis */
    bl_t *lmps;     /* lead monomials as short divmask */
    sdm_t *lm;      /* lead monomials as short divmask */
    bl_t lml;       /* number of lead monomials of non redundant
                       elements in basis */
    int8_t *red;    /* tracks redundancy of basis elements */
    hm_t **hm;      /* hashed monomials representing exponents */
    cf32_t **cf_ff; /* coefficients for finite fields (32bits) */
    mpq_t **cf_qq;  /* coefficients for rationals */
};

typedef struct mat_t mat_t;
struct mat_t
{
    hm_t **r;       /* rows of the matrix, only column entries, coefficients */
                    /* are handled via linking to coefficient arrays */
    cf32_t **cf_ff; /* coefficients for finite fields (32bits) */
    mpq_t **cf_qq;  /* coefficients for rationals */
    len_t sz;       /* number of rows allocated resp. size */
    len_t np;       /* number of new pivots */
    len_t nr;       /* number of rows set */
    len_t nc;       /* number of columns */
    len_t nru;      /* number of upper rows (in ABCD splicing) */
    len_t nrl;      /* number of lower rows (in ABCD splicing) */
    len_t ncl;      /* number of left columns (in ABCD splicing) */
    len_t ncr;      /* number of right columns (in ABCD splicing) */
};

/* statistic stuff */
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
    int64_t num_redundant_old;
    int64_t num_redundant;
    int64_t num_rht;
    int64_t num_rowsred;
    int64_t num_zerored;

    int32_t ngens;
    int32_t nvars;
    int32_t mnsel;
    int32_t homogeneous;
    int32_t fc;
    int32_t mo;
    int32_t laopt;
    int32_t init_hts;
    int32_t nthrds;
    int32_t reset_ht;
    int32_t current_rd;
    int32_t current_deg;
    int64_t max_bht_size;
    int64_t max_sht_size;
    int64_t max_uht_size;
    int64_t nterms_basis;
    int32_t size_basis;

    int32_t info_level;
    int32_t gen_pbm_file;
};

/* function pointers */
bs_t *(*initialize_basis)(
        const int32_t ngens
        );
void (*check_enlarge_basis)(
        bs_t *bs,
        const len_t added
        );

void (*normalize_initial_basis)(
        bs_t *bs,
        const int32_t fc
        );

int (*initial_input_cmp)(
        const void *a,
        const void *b,
        void *ht
        );

int (*monomial_cmp)(
        const hl_t a,
        const hl_t b,
        const ht_t *ht
        );

int (*spair_cmp)(
        const void *a,
        const void *b,
        void *htp
        );

int (*hcm_cmp)(
        const void *a,
        const void *b,
        void *htp
        );

void (*import_julia_data)(
        bs_t *bs,
        ht_t *ht,
        stat_t *st,
        const int32_t *lens,
        const int32_t *exps,
        const void *vcfs
        );

int64_t (*export_julia_data)(
        int32_t *bload,
        int32_t **blen,
        int32_t **bexp,
        void **bcf,
        const bs_t * const bs,
        const ht_t * const ht
        );

/* linear algebra routines */
void (*linear_algebra)(
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st
        );

cf32_t *(*reduce_dense_row_by_old_pivots)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t * const * const pivs,
        const hl_t dpiv,
        const int32_t fc
        );

hm_t *(*reduce_dense_row_by_known_pivots_sparse)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        hm_t *const *pivs,
        const hl_t dpiv,
        const hm_t tmp_pos,
        const int32_t fc
        );

cf32_t *(*reduce_dense_row_by_all_pivots)(
        int64_t *dr,
        mat_t *mat,
        const bs_t * const bs,
        len_t *pc,
        hm_t *const *pivs,
        cf32_t *const *dpivs,
        const int32_t fc
        );


cf32_t *(*reduce_dense_row_by_dense_new_pivots)(
        int64_t *dr,
        len_t *pc,
        cf32_t * const * const pivs,
        const len_t ncr,
        const int32_t fc
        );

/* -----------------------------------
 * non-static functions and procedures
 * ----------------------------------- */
int64_t f4_julia(
        int32_t *bld,   /* basis load */
        int32_t **blen, /* length of each poly in basis */
        int32_t **bexp, /* basis exponent vectors */
        void **bcf,     /* coefficients of basis elements */
        const int32_t *lens,
        const int32_t *exps,
        const void *cfs,
        const int32_t field_char,
        const int32_t mon_order,
        const int32_t nr_vars,
        const int32_t nr_gens,
        const int32_t ht_size,
        const int32_t nr_threads,
        const int32_t max_nr_pairs,
        const int32_t reset_hash_table,
        const int32_t la_option,
        const int32_t pbm_file,
        const int32_t info_level
        );
#endif
