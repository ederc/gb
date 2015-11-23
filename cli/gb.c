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
  printf("    -r REDGB    Compute the reduced Groebner basis.\n");
  printf("                Default: 0.\n");
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
  int nbr_threads       = 1;
  int reduce_gb         = 0;

  int index;
  int opt;

  /*  timing structs */
  struct timeval t_load_start;
  struct timeval t_complete;

	opterr  = 0;

  while ((opt = getopt(argc, argv, "hr:t:v:")) != -1) {
    switch (opt) {
      case 'h':
        print_help();
        return 0;
      case 'r':
        reduce_gb  = atoi(optarg);
        break;
      case 't':
        nbr_threads = atoi(optarg);
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
    printf("---------------------------------------------------------------------\n");
    printf("-------------------------- Computing Groebner -----------------------\n");
    printf("------------------ with the following options set -------------------\n");
    printf("---------------------------------------------------------------------\n");
    printf("number of threads           %4d\n", nbr_threads);
    printf("compute reduced basis?      %4d\n", reduce_gb);
    printf("---------------------------------------------------------------------\n");
  }

  // generators stores input data
  gen_t *generators  = NULL;

  if (verbose > 0) {
    printf("---------------------------------------------------------------------\n");
    gettimeofday(&t_load_start, NULL);
    printf("%-38s","Loading input data ...");
    fflush(stdout);
  }
  if (verbose > 0) {
    printf("%9.3f sec (%.3f %s/sec)\n",
        walltime(t_load_start) / (1000000),
        generators->fs / (walltime(t_load_start) / (1000000)), generators->fsu);
  }  
  if (verbose > 1) {
    print_mem_usage();
    printf("---------------------------------------------------------------------\n");
    printf("Data for %s\n", fn);
    printf("---------------------------------------------------------------------\n");
    printf("field characteristic        %14d\n", generators->modulus);
    printf("number of generators        %14d\n", generators->nbr_vars);
    printf("---------------------------------------------------------------------\n");
  }
}
