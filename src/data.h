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
#include <stdio.h>
#include <string.h> /* for memset et al. */
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
typedef int32_t len_t;    /* length type for different structures */
typedef int32_t exp_t;    /* exponent type */
typedef int32_t val_t;    /* core values like hashes, coefficients, etc. */
typedef int32_t deg_t;    /* (total) degree of polynomial */
typedef int32_t bi_t;     /* basis index of element */
typedef int32_t bl_t;     /* basis load */
typedef int32_t pl_t;     /* pair set load */
typedef int32_t *poly_t;  /* polynomials are just arrays of coeffs
                             and exponent hashes. They start with the
                             length of the array and the second entry
                             is used for distinguishing between known
                             basis elements and new basis elements when
                             applying the linear algebra. That means
                             p = [len, known, cf1, eh1, cf2, eh2, ...] */


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

/* basis data */
static poly_t *bs = NULL;
static bl_t bload = 0;
static bl_t bsize = 0;

/* -----------------------------------
 * non-static functions and procedures
 * ----------------------------------- */
int32_t *f4_julia(
    int32_t *lens,
    int32_t *cfs,
    int32_t *exps,
    int32_t nr_vars,
    int32_t nr_gens,
    int32_t ht_size
    );
#endif
