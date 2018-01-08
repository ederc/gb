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
#if defined(ENABLE_OPENMP)
#include <omp.h>
#else
typedef int omp_int_t;
inline omp_int_t omp_get_thread_num(void) { return 0;}
inline omp_int_t omp_get_max_threads(void) { return 1;}
#endif

/* computational data */
typedef int32_t sdm_t;    /* short divmask for faster divisibility checks */
typedef int32_t len_t;    /* length type for different structures */
typedef int32_t exp_t;    /* exponent type */
typedef int32_t val_t;    /* core values like hashes, coefficients, etc. */
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
  spt_t type;
  val_t lcm;
  bi_t gen1;
  bi_t gen2;
};

/* pair set data */
static spair_t *ps  = NULL;
static pl_t pload   = 0;
static pl_t psize   = 0;

/* field characteristic */
static val_t fc = 0;

/* basis data */
static val_t **bs = NULL;
static bl_t bload = 0;
static bl_t bsize = 0;

static const long bred  = (long)1;  /* maRking redundant elements */
static const long bmask = ~(long)1; /* maSking redundant elements */

/* matrix data */
static int32_t nrall  = 0; /* allocated rows for matrix */
static int32_t nrows  = 0; /* rows used in the current round */
static int32_t ncols  = 0; /* columns used in the current round */
static int32_t npivs  = 0; /* new pivots in the current round */
static int32_t nru    = 0; /* number of upper rows (in ABCD splicing) */
static int32_t nrl    = 0; /* number of lower rows (in ABCD splicing) */
static int32_t ncl    = 0; /* number of left columns(in ABCD splicing) */
static int32_t ncr    = 0; /* number of right columns(in ABCD splicing) */

/* loop unrolling in sparse linear algebra:
 * we store the offset of the first elements not unrolled
 * in the second entry of the sparse row resp. sparse polynomial. */
#define UNROLL  4

/* threads data */
static int32_t nthrds = 1; /* number of CPU threads */

/* linear algebra options */
static int32_t laopt  = 0;

/* function pointers */
int (*monomial_cmp)(
    const len_t a,
    const len_t b
    );

int (*hcm_cmp)(
    const void *a,
    const void *b
    );

val_t **(*select_spairs)(
    void
    );

/* linear algebra routines */
val_t **(*sparse_linear_algebra)(
    val_t **mat
    );

/* -----------------------------------
 * non-static functions and procedures
 * ----------------------------------- */
int32_t *f4_julia(
    const int32_t *lens,
    const int32_t *cfs,
    const int32_t *exps,
    const int32_t field_char,
    const int32_t nr_vars,
    const int32_t nr_gens,
    const int32_t ht_size
    );
#endif
