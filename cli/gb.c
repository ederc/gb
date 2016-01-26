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
#include <getopt.h>

// extern declaration in src/types.h
info_t *meta_data;

void print_help()
{
  printf("\n");
  printf("NAME\n");
  printf("    gb - Groebner basis computations\n");
  printf("\n");
  printf("SYNOPSIS\n");
  printf("    gb [options] [file..]\n");
  printf("\n");
  printf("DESCRIPTION\n");
  printf("    Computes a Groebner basis of the given input ideal.\n");
  printf("\n");
  printf("OPTIONS\n");
  printf("    -c HTC      Hash table cache resp. size.\n");
  printf("                Default: 2^(16).\n");
  printf("    -h HELP     Print help.\n");
  printf("    -p PBM      Generates .pbm files of gbla matrices.\n");
  printf("                Note: These files can become huge, handle with care.\n");
  printf("    -r REDGB    Compute the reduced Groebner basis.\n");
  printf("                Default: 0.\n");
  printf("    -s SIMP     Use simplify in F4 algorithm.\n");
  printf("                Default: 0 (off).\n");
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
  const char *fn        = NULL;
  int verbose           = 0;
  int nthreads          = 1;
  int reduce_gb         = 0;
  ht_size_t ht_size     = pow(2,16);
  int simplify          = 0;
  int pbm               = 0;

  // generate file name holder if pbms are generated
  char pbm_fn[400];

  int index;
  int opt;

  int steps = 0;

  /*  timing structs */
  struct timeval t_load_start;
  struct timeval t_complete;

  // keep track of meta data, meta_data is global to be used wherever we need it
  // to keep track of data
  meta_data = init_meta_data();

	opterr  = 0;

  while ((opt = getopt(argc, argv, "c:hpr:s:t:v:")) != -1) {
    switch (opt) {
      case 'c':
        ht_size = atoi(optarg);
        break;
      case 'h':
        print_help();
        return 0;
      case 'p':
        pbm = 1;
        break;
      case 'r':
        reduce_gb = atoi(optarg);
        break;
      case 's':
        simplify = atoi(optarg);
        if (simplify > 0)
          simplify = 1;
        break;
      case 't':
        nthreads  = atoi(optarg);
        break;
      case 'v':
        verbose = atoi(optarg);
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
  if (verbose > 1) {
    double log_size = log2(ht_size);
    printf("---------------------------------------------------------------------\n");
    printf("-------------------------- Computing Groebner -----------------------\n");
    printf("------------------ with the following options set -------------------\n");
    printf("---------------------------------------------------------------------\n");
    printf("number of threads           %15d\n", nthreads);
    printf("hash table size             %15d (+/- 2^%.1f)\n", ht_size, log_size);
    printf("compute reduced basis?      %15d\n", reduce_gb);
    printf("use simplify?               %15d\n", simplify);
    printf("generate pbm files?         %15d\n", pbm);
    printf("---------------------------------------------------------------------\n");
  }

  // first get number of variables in order to initialize hash table
  nvars_t nvars = get_nvars(fn);
  // initialize hash table
  ht = init_hash_table(ht_size, nvars);
  // input stores input data
  gb_t *basis = load_input(fn, nvars, ht, verbose, nthreads);

  if (verbose > 0) {
    printf("---------------------------------------------------------------------\n");
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Loading input data ...");
    fflush(stdout);
  }
  if (verbose > 0) {
    printf("%9.3f sec\n", walltime(t_load_start) / (1000000));
  }  
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("Data for %s\n", fn);
    printf("---------------------------------------------------------------------\n");
    printf("field characteristic        %15d\n", basis->mod);
    printf("number of variables         %15d\n", basis->nv);
    // See note on gb_t in src/types.h why we decrement basis->load here.
    printf("number of generators        %15d\n", basis->load-1);
    printf("basis file size             %18.2f %s\n", basis->fs, basis->fsu);
  }

  /*  track time for the complete reduction process (excluding load) */
  if (verbose > 0)
    gettimeofday(&t_complete, NULL);
  // initialize spair set
  ps_t *ps = initialize_pair_set(basis, ht);

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


  nelts_t i;
  hash_t hv;

  // run while there exist spairs to be handled
  while (ps->load > 0)
  {
    steps++;
    if (verbose > 1)
      printf(">>> step %u\n", steps);
    
    // select next bunch of spairs
    spd_t *spd  = symbolic_preprocessing(ps, basis);

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
    sort_columns_by_lead(spd);

    // next we can sort both parts (we know the number of lead monomials =
    // spd->sel->load - spd->sel->nsp, i.e. all polynomials considered in
    // symbolic preprocessing minus the number of spairs)
    sort_presorted_columns_by_grevlex(spd, nthreads);

    // connect monomial hash positions with columns in to be constructed gbla
    // matrix
    set_column_index_in_hash_table(ht, spd);

    // now sort upper and lower polynomial selections by lead monomials. this
    // corresponds to sorting by columns (=rows as this is one-to-one for lead
    // monomials).
    sort_selection_by_column_index(spd, ht, nthreads);

    for (int ii=0; ii<spd->col->load; ++ii) {
      printf("%u %u | ", ht->idx[spd->col->hpos[ii]], spd->col->hpos[ii]);
      for (int jj=0; jj<ht->nv; ++jj)
        printf("%u ", ht->exp[spd->col->hpos[ii]][jj]);
      printf("\n");
    }

    // generate gbla matrix out of data from symbolic preprocessing
    mat_t *mat  = generate_gbla_matrix(basis, spd, nthreads);
    // generate pbm files of gbla matrix
    if (pbm) {
      snprintf(pbm_fn, 300, "%s-mat%u.pbm", fn, steps);
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
    }
    // generate pbm files of gbla matrix
    if (pbm) {
      snprintf(pbm_fn, 300, "%s-mat%u-red.pbm", fn, steps);
      write_reduced_matrix_to_pbm(mat, pbm_fn);
    }

    // add new elements to basis and update pair set
    for (i=0; i<rankDR; ++i) {
      hv  = add_new_element_to_basis_grevlex(basis, mat, i, spd, ht);
      // if hash value 0 is new lead monomial we are done, since we have found a
      // unit in the basis, i.e. basis = { 1 }
      if (hv == 0)
        goto done;
      update_pair_set(ps, basis, basis->load-1);
    }
    free_gbla_matrix(mat);
    free_symbolic_preprocessing_data(spd);
    clear_hash_table_idx(ht);
    if (verbose > 1)
      printf("<<< step %u\n", steps);
  }
  done:
  // free allocated memory
  free(meta_data);
  free_pair_set(ps);
  free_basis(basis);
  free_hash_table(ht);
  if (verbose > 0) {
    printf("---------------------------------------------------------------------\n");
    printf("%-38s","Computation completed ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        walltime(t_complete) / (1000000));
    if (verbose > 1) 
      print_mem_usage();
  }
  printf("---------------------------------------------------------------------\n");
  return 0;

}
