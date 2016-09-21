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
#include <src/gb_config.h>
#include <omp.h>
#include <src/hash.h>
#include <src/poly.h>
#include <src/matrix.h>

#ifndef META_DATA_DEBUG
#define META_DATA_DEBUG 0
#endif

#ifndef IO_DEBUG
#define IO_DEBUG  0
#endif

/**
 * Checks is a line of the input file is just empty resp. consists only
 * of whitespaces
 *
 * \param line line
 *
 * \return 1 if the line is empty, else 0
 */
static inline int is_line_empty(const char *line)
{
  while (*line != '\0') {
    if (!isspace(*line))
      return 0;
    line++;
  }
  return 1;
}

/**
 * Returns a sorted and possibly saturated (i.e. elements in
 * the basis that are redundant after saturation are marked) groebner basis.
 *
 * \param resulting groebner basis basis
 *
 * \return final groebner basis, sorted and possibly saturated
 */
static inline poly_t *final_basis_for_output(gb_t *basis)
{
  // if we have a unit in the basis we add it "by hand"
  if (basis->has_unit == 1) {
    basis->fl = 1;
    poly_t *fb  = (poly_t *)malloc(1 * sizeof(poly_t));
    fb[0].cf    = (coeff_t *)malloc(1 * sizeof(coeff_t));
    fb[0].eh    = (hash_t *)malloc(1 * sizeof(hash_t));
    fb[0].cf[0] = 1;
    fb[0].eh[0] = 0;
    fb[0].nt    = 1;
    fb[0].red   = 0;

    return fb;
  }

  // else we sort the basis
  nelts_t bs  = basis->load - basis->st - basis->nred;
  poly_t *fb  = (poly_t *)malloc(bs * sizeof(poly_t));
  nelts_t np = 0;
  for (int i=basis->st; i<basis->load; ++i) {
    if (basis->red[i] == 0) {
      fb[np].cf   = basis->cf[i];
      fb[np].eh   = basis->eh[i];
      fb[np].nt   = basis->nt[i];
      fb[np].red  = 0;
      np++;
    }
  }
  // sort final basis
  qsort(fb, np, sizeof(poly_t), ht->sort.compare_polynomials_inverse);

  fb  = realloc(fb, np*sizeof(poly_t));

  // if the given monomial order is not degree compatible and the input system
  // is not homogeneous the computation was homogenized, thus we have to
  // saturate, i.e. recheck for redundancy when setting the homogenization
  // variable to 1.
  // note that we have already removed redundant elements, but when saturating
  // more redundancy might appear.
  if (basis->init_hom == 0 && basis->hom == 1) {
    for (int i=0; i<bs; ++i) {
      for (int j=0; j<i; ++j) {
        if (check_monomial_division_saturated(fb[i].eh[0], fb[j].eh[0], ht) == 1) {
          fb[i].red = 1; // any value =/= 0 is OK to mark it for deletion
          np--;
          break;
        }
      }
    }
  }
  // final size of basis
  basis->fl = np;

  return fb;
}

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
 * \brief Sorts input polynomials w.r.t. the given monomial order (stored in the
 * hash table ht). All terms of all polynomials are sorted. Afterwards the
 * polynomials themselves are sorted by increasing lead term.
 *
 * \param input elements resp. intermediate basis basis
 *
 * \param hash table ht
 */
void sort_input_polynomials(gb_t *basis, const mp_cf4_ht_t *ht);

/**
 * \brief Homogenizes input polynomials with an extra variable.
 *
 * \note When initializing the hash table we always create a "+1" offset for
 * such a homogenizing variable, thus we do not have to recreate the hash table.
 *
 * \note This function is mostly used when we have a monomial order that is not
 * degree compatible and the input system is not homogeneous
 *
 * \param input elements resp. intermediate basis basis
 *
 * \param hash table ht
 */
void homogenize_input_polynomials(gb_t *basis, mp_cf4_ht_t *ht);

/**
 * \brief Loads input file and intializes input data structure
 *
 * \note The input files need to have the following format at the moment:
 * First line: comma-separated list of variables, e.g. x1,x2,x3
 * Second line: modulus of the field, e.g. 32003
 * Following lines: in each line one polynomial, no shorthand notation, i,e,
 * between each coefficient and each power of a variable there is a "*" and for
 * each power of a variable there is "^" (whereas you can write "x1" instead of
 * "x1^1").
 *
 * \note The input polynomials are directly normalized when read in.
 *
 * \param file name fn
 *
 * \param number of variables nvars
 *
 * \param chosen monomial ordering ordering
 *
 * \param hash table to store the exponents ht
 *
 * \param should simplifications be used? simplify
 *
 * \param maximal number of spairs handled at once
 *
 * \param level of verbosity vb
 *
 * \param number of threads nthrds
 *
 * \return initial state of input
 */
gb_t *load_input(const char *fn, const nvars_t nvars, const int ordering,
    mp_cf4_ht_t *ht, const int simplify, const long max_spairs,
    const int vb, const int nthrds);

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
 * \brief Gets variable name w.r.t prev_pos in line
 *
 * \param line in file representing polynomial line
 *
 * \param pointer to last position checked in line prev_pos
 *
 * \return variable name as char *
 */
char *get_variable_name(const char *line, char **prev_pos);

/**
 * \brief Gets next term out of line resp. polynomial
 *
 * \param line in file representing polynomial line
 *
 * \param pointer to previous char * for last term prev_pos
 *
 * \param pointer to previous index of last term prev_idx
 *
 * \param pointer to term in which we store the new term term
 */
void get_term(const char *line, char **prev_pos,
    char **term);

/**
 * \brief Stores exponent vector of term in last entry of exp in hash table
 *
 * \param term term
 *
 * \param intermediate groebner basis basis
 *
 * \param hash table ht
 */
void store_exponent(const char *term, const gb_t *basis, mp_cf4_ht_t *ht);

/**
 * \brief Returns number of terms in polynomial represented by a line in a text
 * file
 *
 * \param text line from file line
 *
 * \return number of terms in polynomial
 */
int get_number_of_terms(const char *line);

/**
 * \brief Writes matrix to pbm file for printing.
 *
 * \param matrix mat
 *
 * \param file name fn
 */
void write_matrix_to_pbm(mat_t *mat, const char *fn);

/**
 * \brief Writes one row from the sparse block matrix A and the dense block
 * matrix B to the buffer for generating the pbm file of the matrix.
 *
 * \param row buffer buffer
 *
 * \param row index idx
 *
 * \param sparse block matrix A
 *
 * \param dense block matrix B
 *
 * \param number of column blocks for A (i.e. lefthand side) cbl
 *
 * \param number of column blocks for B (i.e. righthand side) cbr
 *
 * \param block size bs
 */
void write_sparse_dense_block_row_to_buffer(char *buffer, const nelts_t idx,
    const sb_fl_t *A, const dbm_fl_t *B, const nelts_t cbl, const nelts_t cbr,
    const bi_t bs);

/**
 * \brief Writes reduced gbla matrix to pbm file for printing.
 *
 * \param matrix mat
 *
 * \param file name fn
 */
void write_reduced_matrix_to_pbm(mat_t *mat, const char *fn);

/**
 * \brief Writes one row from the upper AB part of the reduced gbla matrix.
 * 
 * \note A is the identity.
 *
 * \param row buffer buffer
 *
 * \param row index idx
 *
 * \param sparse block matrix A
 *
 * \param dense block matrix B
 *
 * \param number of column blocks for A (i.e. lefthand side) cbl
 *
 * \param number of column blocks for B (i.e. righthand side) cbr
 *
 * \param block size bs
 */
void write_upper_part_row_to_buffer(char *buffer, const nelts_t idx,
    const mat_t *mat);

/**
 * \brief Writes one row from the lower CD part of the reduced gbla matrix.
 * 
 * \note C is zero.
 *
 * \note D is not a dense block matrix, but only a dense row matrix.
 *
 * \param row buffer buffer
 *
 * \param row index idx
 *
 * \param sparse block matrix C
 *
 * \param dense row matrix D
 *
 * \param number of column blocks for C (i.e. lefthand side) cbl
 *
 * \param number of column blocks for D (i.e. righthand side) cbr
 *
 * \param block size bs
 */
void write_lower_part_row_to_buffer(char *buffer, const nelts_t idx,
    const mat_t *mat);

/**
 * \brief Initializes meta data information
 */
static inline info_t *init_meta_data()
{
  info_t *meta_data = (info_t *)malloc(sizeof(info_t));
  meta_data->nred_last      = 0;
  meta_data->nred_total     = 0;
  meta_data->ncrit_last     = 0;
  meta_data->ncrit_total    = 0;
  meta_data->nzerored_last  = 0;
  meta_data->nzerored_total = 0;

  return meta_data;
}

/**
 * \brief Inverts coefficient x w.r.t. to modulus. Stores the inverted value in
 * x in place.
 *
 * \param pointer to coefficient x
 *
 * \param modulus w.r.t. which the inverse is computed modulus
 */
void inverse_coefficient(coeff_t *x, const coeff_t modulus);

/**
 * \brief Prints resulting groebner basis to stdout, sorted w.r.t. the given
 * monomial order.
 *
 * \param groebner basis basis
 *
 * \param already sorted and possibly saturated final basis fb
 */
void print_basis(const gb_t *basis, const poly_t *fb);

/**
 * \brief Prints resulting groebner basis in Singular style, sorted w.r.t. the
 * given monomial order.
 *
 * \param groebner basis basis
 *
 * \param already sorted and possibly saturated final basis fb
 */
void print_basis_in_singular_format(const gb_t *basis, const poly_t *fb);
#endif
