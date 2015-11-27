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
 * \file io.h
 * \brief Input/output routines for matrices
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_IO_H
#define GB_IO_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <gb_config.h>
#include <omp.h>
#include "types.h"

/**
 * \brief Gets number of variables, needs to be done before reading file
 * completely as we need to initialize the hash table beforehand
 *
 * \param file name fn
 *
 * \return number of variables
 */
nvars_t get_nvars(const char *fn);

/**
 * \brief Loads input file and intializes basis
 *
 * \param file name fn
 *
 * \param number of variables nvars
 *
 * \param level of verbosity vb
 *
 * \param number of threads nthrds
 *
 * \return initial state of groebner basis
 */
gb_t *load_input(const char *fn, nvars_t nvars, int vb, int nthrds);

/*  ========== TIMINGS and MEMORY PRINTING ========== */

/**
 * \brief Returns walltime since time t_start.
 *
 * \param Start time stamp t_start
 *
 * \return Difference in walltime between current time stamp and t_start
 */
double walltime(struct timeval t_start);

/**
 * \brief Prints current memory usage.
 *
 * \note This clearly depends on the operating system and is only tested for
 * UNIX(>=3.1.4) and OS X (>=10.11).
 */
void print_mem_usage();

/**
 * \brief Initializes meta data information
 */
static inline info_t *init_meta_data()
{
  info_t *meta_data = (info_t *)malloc(sizeof(meta_data));
  meta_data->n_reduced          = 0;
  meta_data->n_pairs_removed    = 0;
  meta_data->n_zero_reductions  = 0;
   
  return meta_data;
}
#endif
