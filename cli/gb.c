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
  printf("    -h HELP     Print help.\n");
  printf("    -r REDGB    Compute the reduced Groebner basis.\n");
  printf("                Default: 0.\n");
  printf("    -t THRDS    Number of threads used.\n");
  printf("                Default: 1.\n");
  printf("    -s HTS      Hash table size.\n");
  printf("                Default: UINT_MAX.\n");
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


  int index;
  int opt;

  int i, steps;

  /*  timing structs */
  struct timeval t_load_start;
  struct timeval t_complete;

  // keep track of meta data
  info_t *meta_data = init_meta_data();

	opterr  = 0;

  while ((opt = getopt(argc, argv, "hr:s:t:v:")) != -1) {
    switch (opt) {
      case 'h':
        print_help();
        return 0;
      case 'r':
        reduce_gb = atoi(optarg);
        break;
      case 's':
        ht_size = atoi(optarg);
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
    printf("number of threads           %20d\n", nthreads);
    printf("hash table size             %20d (+/- 2^%.1f)\n", ht_size, log_size);
    printf("compute reduced basis?      %20d\n", reduce_gb);
    printf("---------------------------------------------------------------------\n");
  }

  // first get number of variables in order to initialize hash table
  nvars_t nvars = get_nvars(fn);
  // initialize hash table
  ht = init_hash_table(ht_size, nvars);
  // basis stores input data
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
    printf("field characteristic        %14d\n", basis->modulus);
    printf("number of variables         %14d\n", basis->nvars);
    printf("number of generators        %14d\n", basis->load);
    printf("input file size             %14.2f %s\n", basis->fs, basis->fsu);
    printf("---------------------------------------------------------------------\n");
  }

  /*  track time for the complete reduction process (excluding load) */
  if (verbose > 0)
    gettimeofday(&t_complete, NULL);
  // initialize spair set
  ps_t *ps = init_pair_set(basis, ht);

  printf("ps->size %u\n",ps->size);
  printf("ps->load %u\n",ps->load);
  printf("basis->size %u\n",basis->size);
  printf("basis->load %u\n",basis->load);

  for (i=1; i<basis->load; ++i)
    //update_basis(ps, basis, i);

  // run while there exist spairs to be handled
  for (steps=1; ps->load>0; ++steps)
  {
    
  }

  // free allocated memory
  free(meta_data);
  free_pair_set_dynamic_data(ps);
  free(ps);
  free_basis_dynamic_data(basis);
  free(basis);
  free_hash_table_dynamic_data(ht);
  free(ht);
  if (verbose > 0) {
    printf("-------------------------------------------------------------------\n");
    printf("%-38s","Computation completed ...");
    fflush(stdout);
    printf("%9.3f sec\n",
        walltime(t_complete) / (1000000));
    if (verbose > 1) 
      print_mem_usage();
  }
  printf("-------------------------------------------------------------------\n");
  return 0;

}
