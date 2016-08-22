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




#include "gb.h"

// extern declaration in src/types.h
info_t *meta_data;

void print_help()
{
  printf("\n");
  printf("NAME\n");
  printf("    gb - Groebner basis computations\n");
  printf("\n");
  printf("SYNOPSIS\n");
  printf("    gb [options] file\n");
  printf("\n");
  printf("DESCRIPTION\n");
  printf("    Computes a Groebner basis of the given input ideal.\n");
  printf("\n");
  printf("OPTIONS\n");
  printf("    -c HTC      Hash table cache resp. size in log_2: The size is then set\n");
  printf("                to the biggest non-Mersenne prime smaller then 2^(given value).\n");
  printf("                Default: 18.\n");
  printf("    -h HELP     Print help.\n");
  printf("    -m MAT      Generates .pbm files of gbla matrices.\n");
  printf("                Considers as argument a folder to write into.\n");
  printf("    -n NREDMAT  If option is set the gbla matrices are not fully reduced.\n");
  printf("                Otherwise the A and thus B part of the gbla matrices are\n");
  printf("                reduced and can thus be used for simplification of further\n");
  printf("                used polynomials.\n");
  printf("    -o ORDER    Order w.r.t. which the Groebner basis is computed.\n");
  printf("                0: Graded/Degree reverse lexicographical order (DRL)\n");
  printf("                1: Lexicographical order (LEX)\n");
  printf("                Default: 0.\n");
  printf("    -p PRINT    Prints resulting groebner basis.\n");
  printf("                0 -> no printing (default if option is not set at all).\n");
  printf("                1 -> default print out of basis.\n");
  printf("                2 -> Singular format print out of basis.\n");
  printf("    -r REDGB    Compute the reduced Groebner basis.\n");
  printf("                Default: 0 (not reduced).\n");
  printf("    -s SIMP     Use simplify in F4 algorithm.\n");
  printf("                Default: 0 (not simplified).\n");
  printf("                         1 (simplified but B not fully reduced).\n");
  printf("                         2 (simplified and B fully reduced).\n");
  printf("    -t THRDS    Number of threads used.\n");
  printf("                Default: 1.\n");
  printf("    -v LEVEL    Level of verbosity:\n");
  printf("                1 -> only error messages printed\n");
  printf("                2 -> some meta information printed\n");
  printf("                3 -> additionally meta information printed\n");
  printf("                Note: Everything >2 is time consuming and might\n");
  printf("                      slow down the overall computations.\n");
  printf("\n");

  return;
}


int main(int argc, char *argv[])
{
  const char *fn    = NULL;
  int verbose       = 0;
  int nthreads      = 1;
  int reduce_gb     = 0;
  int simplify      = 0;
  int generate_pbm  = 0;
  int print_gb      = 0;
  int order      = 0;
  int keep_A        = 0;
  int htc           = 15;
  // generate file name holder if pbms are generated
  char *pbm_dir = NULL;
  char pbm_fn[400];

  // monomial order names storage, for printing purpose only
  char orders[10][10];
  snprintf(orders[0], 10, "DRL");
  snprintf(orders[1], 10, "LEX");

  int index;
  int opt;
  int done;

  int steps = 0;

  /*  timing structs */
  struct timeval t_load_start;
  struct timeval t_complete;
  double t_linear_algebra         = 0;
  double t_symbolic_preprocessing = 0;
  double t_sorting_columns        = 0;
  double t_generating_gbla_matrix = 0;
  double t_update_pairs           = 0;
  unsigned long n_zero_reductions = 0;

  // keep track of meta data, meta_data is global to be used wherever we need it
  // to keep track of data
  meta_data = init_meta_data();

	opterr  = 0;

  while ((opt = getopt(argc, argv, "c:hm:no:p:r:s:t:v:")) != -1) {
    switch (opt) {
      case 'c':
        htc = (int)strtol(optarg, NULL, 10);
        break;
      case 'h':
        print_help();
        return 0;
      case 'm':
        pbm_dir       = optarg;
        generate_pbm  = 1;
        break;
      case 'n':
        keep_A  = 1;
        break;
      case 'o':
        order  = (int)strtol(optarg, NULL, 10);
        break;
      case 'p':
        print_gb = (int)strtol(optarg, NULL, 10);
        if (print_gb > 2)
          print_gb = 0;
        break;
      case 'r':
        reduce_gb = (int)strtol(optarg, NULL, 10);
        break;
      case 's':
        simplify = (int)strtol(optarg, NULL, 10);
        if (simplify > 2)
          simplify = 1;
        break;
      case 't':
        nthreads  = (int)strtol(optarg, NULL, 10);
        if (nthreads == 0)
          nthreads  = 1;
        break;
      case 'v':
        verbose = (int)strtol(optarg, NULL, 10);
        if (verbose > 4)
          verbose = 4;
        break;
      case '?':
        if (optopt == 'f')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
              "Unknown option character `\\x%x'.\n",
              optopt);
        return 1;
      default:
        abort ();
    }
  }
  for (index = optind; index < argc; index++)
    fn = argv[index];

  if (fn == NULL) {
    fprintf(stderr, "File name is required.\nSee help using '-h' option.\n");
    return 1;
  }
  /*
  if (order != 0) {
    fprintf(stderr, "At the moment only computations w.r.t. the degree reverse lexicographical\norder are possible.\nSee help using '-h' option.\n");
    return 1;
  }
  */
  // we start with hash table sizes being non-Mersenne primes < 2^18, nothing
  // smaller. So we shift htc to the corresponding index of ht->primes.
  if (htc < 15)
    htc =   15;

  // first get number of variables in order to initialize hash table
  nvars_t nvars = get_nvars(fn);
  // initialize hash table
  ht = init_hash_table(htc, nvars);
  set_sort_functions_depending_on_monomial_order(ht, order);

  if (verbose > 0) {
    printf("---------------------------------------------------------------------------\n");
    printf("----------------------------- Computing Groebner --------------------------\n");
    printf("--------------------- with the following options set ----------------------\n");
    printf("---------------------------------------------------------------------------\n");
    printf("number of threads           %15d\n", nthreads);
    printf("hash table size             %15u (2^%u)\n", ht->sz, htc);
    printf("compute reduced basis?      %15d\n", reduce_gb);
    printf("do not reduce A|B in gbla   %15d\n", keep_A);
    printf("use simplify?               %15d\n", simplify);
    printf("generate pbm files?         %15d\n", generate_pbm);
    printf("print resulting basis?      %15d\n", print_gb);
    printf("---------------------------------------------------------------------------\n");
  }
  // input stores input data
  gb_t *basis = load_input(fn, nvars, order, ht, simplify, verbose, nthreads);

  // global simplifier list
  // generate simplifier list if simplification is enabled
  gb_t *sf  = NULL;
  if (basis->sl > 0)
    sf = initialize_simplifier_list(basis);

  if (verbose > 0) {
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Loading input data ...");
    fflush(stdout);
  }
  if (verbose > 0)
    printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

  if (verbose > 1)
    print_mem_usage();

  if (verbose > 0) {
    printf("---------------------------------------------------------------------------\n");
    printf("Data for %s\n", fn);
    printf("---------------------------------------------------------------------------\n");
    printf("field characteristic        %15d\n", basis->mod);
    printf("monomial order              %15s\n", orders[basis->ord]);
    printf("number of variables         %15d\n", basis->rnv);
    // See note on gb_t in src/types.h why we decrement basis->load here.
    printf("number of generators        %15d\n", basis->load-1);
    printf("homogeneous input?          %15d\n", basis->init_hom);
    printf("homogeneous computation?    %15d\n", basis->hom);
    printf("input file size             %18.2f %s\n", basis->fs, basis->fsu);
    printf("---------------------------------------------------------------------------\n");
  }
#if __GB_HAVE_SSE2
  if (verbose > 2)
    printf("ht->nev: %u | %u | %u\n",ht->nv, ht->nev, ht->vl);
#endif

  // track time for the complete reduction process (excluding load)
  if (verbose > 0) {
    gettimeofday(&t_complete, NULL);
    gettimeofday(&t_load_start, NULL);
  }
  // initialize spair set
  ps_t *ps = initialize_pair_set(basis, ht);

  if (verbose > 0)
    t_update_pairs  +=  walltime(t_load_start);

  if (verbose > 2) {
    printf("---------------------------------------------------------------------------\n");
    printf("ps->size                          %9u\n",ps->size);
    printf("ps->load                          %9u\n",ps->load);
    printf("basis->size                       %9u\n",basis->size);
    // See note on gb_t in src/types.h why we decrement basis->load here.
    printf("basis->load                       %9u\n",basis->load-1); 
    printf("---------------------------------------------------------------------------\n");
    printf("criteria applications (last step) %9u\n",meta_data->ncrit_last);
    printf("criteria applications (total)     %9u\n",meta_data->ncrit_total);
  }


  // run while there exist spairs to be handled
  while (ps->load > 0)
  {
    steps++;
    
    // select next bunch of spairs
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);
    spd_t *spd  = symbolic_preprocessing(ps, basis, sf);
    if (verbose > 0)
      t_symbolic_preprocessing +=  walltime(t_load_start);

    if (verbose > 2) {
      printf("---------------------------------------------------------------------------\n");
      printf("sel->deg                          %9u\n",spd->selu->deg);
      printf("selu->load                        %9u\n",spd->selu->load);
      printf("sell->load                        %9u\n",spd->sell->load);
      printf("mon->load                         %9u\n",spd->col->load);
      printf("mon->nlm                          %9u\n",spd->col->nlm);
    }

    // next we have to store arrays for the connection between lead monomials
    // and matrix column indices resp. non-lead monomials and matrix column
    // indices. Note that we already know the sizes of the arrays due to
    // symbolic preprocessing:
    // We first sort spd->col via lead and non lead monomials, i.e. ht->idx[i] =
    // 1 or = 2
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);
    sort_columns_by_lead(spd);

    // next we can sort both parts (we know the number of lead monomials =
    // spd->sel->load - spd->sel->nsp, i.e. all polynomials considered in
    // symbolic preprocessing minus the number of spairs)
    
    // if we do not keep A we use for A and C a row structure and thus do
    // not invert the order of the left side columns.
    // if we keep A we use for A and C a block structure and we invert
    // the order of the left side columns since gbla expects blocks be
    // inverted for a faster reduction of A in the first step of the matrix
    // reduction..
    if (keep_A == 1)
      ht->sort.sort_presorted_columns(spd, nthreads);
    else
      ht->sort.sort_presorted_columns_invert_left_side(spd, nthreads);

    // connect monomial hash positions with columns in to be constructed gbla
    // matrix
    set_column_index_in_hash_table(ht, spd);

    // now sort upper and lower polynomial selections by lead monomials. this
    // corresponds to sorting by columns (=rows as this is one-to-one for lead
    // monomials).
    // if we do not keep A we use for A and C a row structure and thus do not
    // invert the order of the left side columns.
    // if we keep A we use for A and C a block structure and we invert
    // the order of the left side columns since gbla expects blocks be
    // inverted for a faster reduction of A in the first step of the matrix
    // reduction..
    if (keep_A == 1)
      sort_selection_by_inverted_column_index(spd, ht, nthreads);
    else
      sort_selection_by_column_index(spd, ht, nthreads);
    if (verbose > 0)
      t_sorting_columns +=  walltime(t_load_start);

#if GB_DEBUG
    printf("# lead terms: %u\n", spd->col->nlm);
    for (int ii=0; ii<spd->col->load; ++ii) {
      for (int jj=0; jj<ht->nv; ++jj) {
        printf("%u ", ht->exp[spd->col->hpos[ii]][jj]);
      }
      printf("|| %lu | %u\n", spd->col->hpos[ii], ht->idx[spd->col->hpos[ii]]);
    }
#endif

    // generate gbla matrix out of data from symbolic preprocessing
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);
    mat_t *mat  = NULL;
    if (keep_A == 1) {
      mat = generate_gbla_matrix_keep_A(basis, sf, spd, nthreads);
      meta_data->mat_rows = mat->AR->nrows+mat->CR->nrows;
      meta_data->mat_cols = mat->AR->ncols+mat->B->ncols;
      if (verbose > 1) {
        printf("matrix rows %6u + %6u = %6u\n", mat->AR->nrows, mat->CR->nrows, mat->AR->nrows + mat->CR->nrows);
        printf("matrix cols %6u + %6u = %6u\n", mat->AR->ncols, mat->B->ncols, mat->AR->ncols + mat->B->ncols);
      }
    } else {
      mat = generate_gbla_matrix(basis, sf, spd, nthreads);
      meta_data->mat_rows = mat->A->nrows+mat->C->nrows;
      meta_data->mat_cols = mat->A->ncols+mat->B->ncols;
      if (verbose > 1) {
        printf("matrix rows %6u + %6u = %6u\n", mat->A->nrows, mat->C->nrows, mat->A->nrows + mat->C->nrows);
        printf("matrix cols %6u + %6u = %6u\n", mat->A->ncols, mat->B->ncols, mat->A->ncols + mat->B->ncols);
      }
    }
    if (verbose == 1) {
      printf("step %3d : %5u spairs  --->  %7u x %7u matrix\n", steps, meta_data->sel_pairs, meta_data->mat_rows, meta_data->mat_cols);
      //printf("%d:%u|%u,%u\n", steps, meta_data->sel_pairs, meta_data->mat_rows, meta_data->mat_cols);
    }
    //mat_t *mat  = generate_gbla_matrix(basis, sf, spd, nthreads);
    if (verbose > 0)
      t_generating_gbla_matrix  +=  walltime(t_load_start);
    // generate pbm files of gbla matrix
    if (generate_pbm) {
      int pos = 0;
      pos = snprintf(pbm_fn+pos, 200, "%s/",pbm_dir);
      pos = snprintf(pbm_fn+pos, 200, "%s-mat%u.pbm", fn, steps);
      printf("%s\n", pbm_fn);
      write_matrix_to_pbm(mat, pbm_fn);
    }
    // reduce matrix using gbla
    if (verbose == 2) {
      gettimeofday(&t_load_start, NULL);
      printf("[%2u] ", steps);
      printf("%-33s","GBLA matrix reduction ...");
      fflush(stdout);
    }
    int rankDR  = 0;
    if (keep_A == 1)
      rankDR  = reduce_gbla_matrix_keep_A(mat, verbose, nthreads);
    else
      rankDR  = reduce_gbla_matrix(mat, verbose, nthreads);
    if (verbose == 2) {
      printf("%9.3f sec %7d %7d %7d\n",
          walltime(t_load_start) / (1000000), rankDR, mat->DR->nrows - rankDR, mat->DR->nrows);
    }
    if (verbose > 0)
      t_linear_algebra  +=  walltime(t_load_start);

    // generate pbm files of gbla matrix
    if (generate_pbm) {
      int pos = 0;
      pos = snprintf(pbm_fn+pos, 200, "%s/",pbm_dir);
      pos = snprintf(pbm_fn+pos, 200, "%s-mat%u-red.pbm", fn, steps);
      printf("%s\n", pbm_fn);
      write_reduced_matrix_to_pbm(mat, pbm_fn);
    }

    // add new elements to basis and update pair set
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);
    
    if (basis->sl == 0)
      done  = update_basis(basis, ps, spd, mat, ht, rankDR);
    else
      done  = update_basis_and_add_simplifier(basis, sf, ps,
          spd, mat, ht, rankDR, nthreads);
    if (verbose > 2) {
      if (basis->sf != NULL)
        printf("basis->load %u | sf->load %u (%u)\n",basis->load, sf->load, spd->col->nlm);
      else
        printf("basis->load %u (%u)\n",basis->load, spd->col->nlm);
    }

    if (verbose > 0) {
      n_zero_reductions +=  (mat->DR->nrows - mat->DR->rank);
    }
    free_gbla_matrix(&mat);
    free_symbolic_preprocessing_data(&spd);
    clear_hash_table_idx(ht);

    if (verbose > 0)
      t_update_pairs  +=  walltime(t_load_start);

    if (verbose > 2) {
      printf("---------------------------------------------------------------------------\n");
      printf("ps->size                          %9u\n",ps->size);
      printf("ps->load                          %9u\n",ps->load);
      printf("basis->size                       %9u\n",basis->size);
      // See note on gb_t in src/types.h why we decrement basis->load here.
      printf("basis->load                       %9u\n",basis->load-1); 
      printf("---------------------------------------------------------------------------\n");
      printf("criteria applications (last step) %9u\n",meta_data->ncrit_last);
      printf("criteria applications (total)     %9u\n",meta_data->ncrit_total);
    }

    // if we are done then we have found the constant 1 as element in the basis
    if (done)
      break;
  }

  // final basis for possible output data
  poly_t *fb  = NULL;
  // generate final basis for output data
  if (verbose > 0 || print_gb > 0) {
    fb  = final_basis_for_output(basis);
  }
  if (verbose > 0) {
    printf("---------------------------------------------------------------------------\n");
    printf("%-38s","Time for updating pairs ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        t_update_pairs / (1000000));
    printf("%-38s","Time for symbolic preprocessing ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        t_symbolic_preprocessing / (1000000));
    printf("%-38s","Time for sorting of columns ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        t_sorting_columns / (1000000));
    printf("%-38s","Time for constructing matrices ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        t_generating_gbla_matrix / (1000000));
    printf("%-38s","Time for linear algebra ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        t_linear_algebra / (1000000));
    printf("---------------------------------------------------------------------------\n");
    printf("%-38s","Computation completed ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        walltime(t_complete) / (1000000));
    printf("---------------------------------------------------------------------------\n");
    printf("Size of basis                     %9u\n", basis->fl);
    printf("Number of zero reductions         %9lu\n", n_zero_reductions);
    printf("---------------------------------------------------------------------------\n");
    if (verbose > 2)
      print_mem_usage();
  }

  // printing of output
  switch (print_gb) {
    case 0:
      break;
    case 1:
      print_basis(basis, fb);
      break;
    case 2:
      print_basis_in_singular_format(basis, fb);
      break;
    default:
      break;
  }

  // free allocated memory
  free(meta_data);
  free_pair_set(&ps);
  if (basis->sl > 0)
    free_simplifier_list(&sf);
  free_basis(&basis);
  free_hash_table(&ht);
  free(fb);

  return 0;
}

void add_simplifier(gb_t *basis, gb_t *sf, mat_t *mat, const spd_t *spd,
    const mp_cf4_ht_t *ht)
{
  if (spd->col->nlm != 0) {
    ri_t i;
    // we add the polys to sf, we know that there is one coefficient at col pos i
    // for row i.
    // note: we only add these simplifiers if they are not too dense, i.e. at
    // most twice the size of the original polynomial in the basis. this is a
    // heuristic and might be bad in some examples, but in most of our tests it
    // is the fastest choice.
    if (mat->sl == 1) {
      for (i=0; i<mat->BR->nrows; ++i) {
        if (1+mat->DR->ncols < 2 * basis->nt[spd->selu->mpp[i].bi]) {
#if 0
          if (mat->BR->row[i]->init_val != NULL) {
            copy_to_val(mat->BR, i);
            reduce_B_by_D(mat->BR, mat->DR, i);
            add_new_element_to_simplifier_list(basis, sf, mat->BR, i, spd, ht);
          }
#else
          add_new_element_to_simplifier_list(basis, sf, mat->BR, i, spd, ht);
#endif
        }
      }
    }
    if (mat->sl > 1) {
      for (i=0; i<mat->BR->nrows; ++i) {
        if (1+mat->DR->ncols < 2*basis->nt[spd->selu->mpp[i].bi]) {
          if (mat->BR->row[i]->init_val != NULL) {
            copy_to_val(mat->BR, i);
            reduce_B_by_D(mat->BR, mat->DR, i);
            add_new_element_to_simplifier_list(basis, sf, mat->BR, i, spd, ht);
          }
        }
      }
    }
    free_dense_row_submatrix(&(mat->BR), 1);
  }
}

int update_basis_and_add_simplifier(gb_t *basis, gb_t *sf, ps_t *ps,
    spd_t *spd, mat_t *mat, const mp_cf4_ht_t *ht,  const ri_t rankDR,
    const int nthreads)
{
  int done;
  const int t = 2<nthreads ? 2 : nthreads;
#pragma omp parallel num_threads(t)
  {
    // update basis and pair set, mark redundant elements in basis
#pragma omp single nowait
    {
      // add simplifier, i.e. polynomials corresponding to the rows in AB,
      // for further computation
#pragma omp task
      {
        add_simplifier(basis, sf, mat, spd, ht);
      }
#pragma omp task
      {
        done  = update_basis(basis, ps, spd, mat, ht, rankDR);
      }
    }
    #pragma omp taskwait
  }
  return done;
}
