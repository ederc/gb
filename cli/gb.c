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
  printf("    -d ORDER    Ordering w.r.t. which the Groebner basis is computed.\n");
  printf("                Graded reverse lexicographical: 0\n");
  printf("                Default: 0.\n");
  printf("    -h HELP     Print help.\n");
  printf("    -o OUTPUT   Prints resulting groebner basis.\n");
  printf("                0 -> no printing (default if option is not set at all).\n");
  printf("                1 -> default print out of basis.\n");
  printf("                2 -> Singular format print out of basis.\n");
  printf("    -p PBM      Generates .pbm files of gbla matrices.\n");
  printf("                Considers as argument a folder to write into.\n");
  printf("    -r REDGB    Compute the reduced Groebner basis.\n");
  printf("                Default: 0 (not reduced).\n");
  printf("    -s SIMP     Use simplify in F4 algorithm.\n");
  printf("                Default: 0 (not simplified).\n");
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
  int ordering      = 0;
  int htc           = 18;
  // generate file name holder if pbms are generated
  char *pbm_dir = NULL;
  char pbm_fn[400];

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

  // keep track of meta data, meta_data is global to be used wherever we need it
  // to keep track of data
  meta_data = init_meta_data();

	opterr  = 0;

  while ((opt = getopt(argc, argv, "c:d:ho:p:r:s:t:v:")) != -1) {
    switch (opt) {
      case 'c':
        htc = (int)strtol(optarg, NULL, 10);
        break;
      case 'd':
        ordering  = (int)strtol(optarg, NULL, 10);
        break;
      case 'h':
        print_help();
        return 0;
      case 'o':
        print_gb = (int)strtol(optarg, NULL, 10);
        if (print_gb > 2)
          print_gb = 0;
        break;
      case 'p':
        pbm_dir       = optarg;
        generate_pbm  = 1;
        break;
      case 'r':
        reduce_gb = (int)strtol(optarg, NULL, 10);
        break;
      case 's':
        simplify = (int)strtol(optarg, NULL, 10);
        if (simplify > 0)
          simplify = 1;
        break;
      case 't':
        nthreads  = (int)strtol(optarg, NULL, 10);
        break;
      case 'v':
        verbose = (int)strtol(optarg, NULL, 10);
        if (verbose > 3)
          verbose = 2;
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
  if (ordering != 0) {
    fprintf(stderr, "At the moment only computations w.r.t. the degree reverse lexicographical\nordering are possible.\nSee help using '-h' option.\n");
    return 1;
  }
  // we start with hash table sizes being non-Mersenne primes < 2^18, nothing
  // smaller. So we shift htc to the corresponding index of ht->primes.
  if (htc < 18)
    htc =   0;
  else
    htc -=  18;

  // first get number of variables in order to initialize hash table
  nvars_t nvars = get_nvars(fn);
  // initialize hash table
  ht = init_hash_table(htc, nvars);

  if (verbose > 0) {
    printf("---------------------------------------------------------------------\n");
    printf("-------------------------- Computing Groebner -----------------------\n");
    printf("------------------ with the following options set -------------------\n");
    printf("---------------------------------------------------------------------\n");
    printf("number of threads           %12d\n", nthreads);
    printf("hash table size             %12u (non-Mersenne prime < 2^%d)\n", ht->primes[ht->si], ht->si+18);
    printf("compute reduced basis?      %12d\n", reduce_gb);
    printf("use simplify?               %12d\n", simplify);
    printf("generate pbm files?         %12d\n", generate_pbm);
    printf("print resulting basis?      %12d\n", print_gb);
    printf("---------------------------------------------------------------------\n");
  }
  // input stores input data
  gb_t *basis = load_input(fn, nvars, ordering, ht, simplify, verbose, nthreads);

  // global simplifier list
  // generate simplifier list if simplification is enabled
  gb_t *sf  = NULL;
  if (simplify == 1)
    sf = initialize_simplifier_list(basis);

  if (verbose > 0) {
    printf("---------------------------------------------------------------------\n");
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Loading input data ...");
    fflush(stdout);
  }
  if (verbose > 0) {
    printf("%9.3f sec\n", walltime(t_load_start) / (1000000));
  }  
  if (verbose > 0) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("Data for %s\n", fn);
    printf("---------------------------------------------------------------------\n");
    printf("field characteristic        %15d\n", basis->mod);
    printf("monomial ordering           %15d\n", basis->ord);
    printf("number of variables         %15d\n", basis->nv);
    // See note on gb_t in src/types.h why we decrement basis->load here.
    printf("number of generators        %15d\n", basis->load-1);
    printf("homogeneous input?          %15d\n", basis->hom);
    printf("input file size             %18.2f %s\n", basis->fs, basis->fsu);
    printf("---------------------------------------------------------------------\n");
  }

  /*  track time for the complete reduction process (excluding load) */
  if (verbose > 0) {
    gettimeofday(&t_complete, NULL);
    gettimeofday(&t_load_start, NULL);
  }
  // initialize spair set
  ps_t *ps = initialize_pair_set(basis, ht);

  if (verbose > 0)
    t_update_pairs  +=  walltime(t_load_start);

  if (verbose > 1) {
    printf("---------------------------------------------------------------------\n");
    printf("ps->size                          %9u\n",ps->size);
    printf("ps->load                          %9u\n",ps->load);
    printf("basis->size                       %9u\n",basis->size);
    // See note on gb_t in src/types.h why we decrement basis->load here.
    printf("basis->load                       %9u\n",basis->load-1); 
    printf("---------------------------------------------------------------------\n");
    printf("criteria applications (last step) %9u\n",meta_data->ncrit_last);
    printf("criteria applications (total)     %9u\n",meta_data->ncrit_total);
  }


  // run while there exist spairs to be handled
  while (ps->load > 0)
  {
    steps++;
    if (verbose > 1)
      printf(">>> step %u\n", steps);
    
    // select next bunch of spairs
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);
    spd_t *spd  = symbolic_preprocessing(ps, basis, sf);
    if (verbose > 0)
      t_symbolic_preprocessing +=  walltime(t_load_start);

    if (verbose > 1) {
      printf("---------------------------------------------------------------------\n");
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
    sort_presorted_columns_by_grevlex(spd, nthreads);

#if GB_DEBUG
    /*
    for (int ii=0; ii<spd->col->load; ++ii) {
      for (int jj=0; jj<ht->nv; ++jj) {
        printf("%u ", ht->exp[spd->col->hpos[ii]][jj]);
      }
      printf("|| %lu | %u\n", spd->col->hpos[ii], ht->idx[spd->col->hpos[ii]]);
    }
    */
#endif

    // connect monomial hash positions with columns in to be constructed gbla
    // matrix
    set_column_index_in_hash_table(ht, spd);

    // now sort upper and lower polynomial selections by lead monomials. this
    // corresponds to sorting by columns (=rows as this is one-to-one for lead
    // monomials).
    sort_selection_by_column_index(spd, ht, nthreads);
    if (verbose > 0)
      t_sorting_columns +=  walltime(t_load_start);

#if GB_DEBUG
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
    mat_t *mat  = generate_gbla_matrix(basis, spd, nthreads);
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
    if (verbose == 1) {
      gettimeofday(&t_load_start, NULL);
      printf("[%2u] ", steps);
      printf("%-33s","GBLA matrix reduction ...");
      fflush(stdout);
    }
    int rankDR  = reduce_gbla_matrix(mat, verbose, nthreads);
    if (verbose == 1) {
      printf("%9.3f sec (rank DR: %d)\n",
          walltime(t_load_start) / (1000000), rankDR);
      t_linear_algebra  +=  walltime(t_load_start);
    }

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
    
    if (simplify == 0)
      done  = update_basis(basis, ps, spd, mat, ht, rankDR);
    else
      done  = update_basis_and_add_simplifier(basis, sf, ps,
          spd, mat, ht, rankDR, nthreads);
    if (verbose > 1)
      printf("basis->load %u | sf->load %u (%u)\n",basis->load, sf->load, spd->col->nlm);
    free_gbla_matrix(mat);
    free_symbolic_preprocessing_data(spd);
    clear_hash_table_idx(ht);

    if (verbose > 0)
      t_update_pairs  +=  walltime(t_load_start);

    if (verbose > 1)
      printf("<<< step %u\n", steps);

    // if we are done then we have found the constant 1 as element in the basis
    if (done)
      break;
  }
  if (verbose > 0) {
    printf("---------------------------------------------------------------------\n");
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
    printf("---------------------------------------------------------------------\n");
    printf("%-38s","Computation completed ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        walltime(t_complete) / (1000000));
    printf("---------------------------------------------------------------------\n");
    printf("size of basis                     %9u\n",basis->load - basis->st - basis->nred);
    if (verbose > 1)
      print_mem_usage();
  }

  // printing of output
  switch (print_gb) {
    case 0:
      break;
    case 1:
      print_basis(basis);
      break;
    case 2:
      print_basis_in_singular_format(basis);
      break;
    default:
      break;
  }

  // free allocated memory
  free(meta_data);
  free_pair_set(ps);
  free_basis(basis);
  if (simplify == 1)
    free_simplifier_list(sf);
  free_hash_table(ht);
  return 0;

}

int update_basis(gb_t *basis, ps_t *ps, spd_t *spd, const mat_t *mat,
    const mp_cf4_ht_t *ht,  const ri_t rankDR)
{
  ri_t i;
  hash_t hash;
  for (i=0; i<rankDR; ++i) {
    // add lowest row first, it has the smallest new lead monomial
    hash = add_new_element_to_basis_grevlex(basis, mat, rankDR-1-i, spd, ht);
    // if hash value 0 is new lead monomial we are done, since we have found a
    // unit in the basis, i.e. basis = { 1 }
    if (hash == 0)
      return 1;
    update_pair_set(ps, basis, basis->load-1);
    track_redundant_elements_in_basis(basis);
  }
  return 0;
}

void add_simplifier_grevlex(gb_t *basis, gb_t *sf, mat_t *mat, const spd_t *spd,
    const mp_cf4_ht_t *ht)
{
  if (spd->col->nlm != 0) {
    ri_t i;
    // store B in dense non-block matrix
    dm_t *B = copy_block_to_dense_matrix(&(mat->B), 1);
    // we add the polys to sf, we know that there is one coefficient at col pos i
    // for row i.
    for (i=0; i<B->nrows; ++i) {
      printf("B->nrows %u / %u\n", i, B->nrows);
      add_new_element_to_simplifier_list_grevlex(basis, sf, B, i, spd, ht);
    }
  }
}

int update_basis_and_add_simplifier(gb_t *basis, gb_t *sf, ps_t *ps,
    spd_t *spd, mat_t *mat, const mp_cf4_ht_t *ht,  const ri_t rankDR,
    const int nthreads)
{
  int done;
#pragma omp parallel num_threads(nthreads)
  {
    // update basis and pair set, mark redundant elements in basis
#pragma omp single nowait
    {
#pragma omp task
      {
        done  = update_basis(basis, ps, spd, mat, ht, rankDR);
      }
      // add simplifier, i.e. polynomials corresponding to the rows in AB,
      // for further computation
#pragma omp task
      {
        add_simplifier_grevlex(basis, sf, mat, spd, ht);
      }
    }
    #pragma omp taskwait
  }
  return done;
}
