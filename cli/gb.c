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
#define COL_CHECK 0

/* extern declaration in src/types.h */
info_t *meta_data;

void print_help(void)
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
  printf("    -b BS       Block size of matrix tiles\n");
  printf("                NOTE: This setting does NOT influence the GBLA block size,\n");
  printf("                      i.e. if using GBLA for linear algebra the block size\n");
  printf("                      is set via __GBLA_SIMD_BLOCK_SIZE (see config.h).\n");
  printf("                      It has also NO influence in the sparse row linear\n");
  printf("                      algebra implementations as we do not use blocks there\n");
  printf("                      at all.\n");
  printf("                      For the block linear algebra implementation we use the\n");
  printf("                      given setting as base for adaptively choosing a best\n");
  printf("                      fitting block size depending on the matrix size.\n");
  printf("                DEFAULT: 4096 = 2^12.\n");
  printf("\n");
  printf("    -c HTC      Hash table cache resp. size in log_2: The size is then set\n");
  printf("                to the biggest non-Mersenne prime smaller then 2^(given value).\n");
  printf("                DEFAULT: 18.\n");
  printf("\n");
  printf("    -g GIT      Outputs git commit hash if verbosity level is >0.\n");
  printf("\n");
  printf("    -h HELP     Print help.\n");
  printf("\n");
  printf("    -l LIMIT    Maximal number of spairs handled at once.\n");
  printf("                NOTE: We try keep all spairs with the same lcm together,\n");
  printf("                      so it might happen that the limit is exceeded slightly.\n");
  printf("\n");
  printf("    -m MAT      Generates .pbm files of gbla matrices.\n");
  printf("                Considers as argument a folder to write into.\n");
  printf("                NOTE: This option is BROKEN at the moment\n");
  printf("\n");
  printf("    -o ORDER    Order w.r.t. which the Groebner basis is computed.\n");
  printf("                0 -> Graded/Degree reverse lexicographical order (DRL)\n");
  printf("                1 -> Lexicographical order (LEX)\n");
  printf("                DEFAULT: 0.\n");
  printf("\n");
  printf("    -p PRINT    Prints resulting groebner basis.\n");
  printf("                0 -> No printing.\n");
  printf("                1 -> Simple printout of basis.\n");
  printf("                2 -> Singular format print out of basis.\n");
  printf("                3 -> Print lead ideal.\n");
  printf("                DEFAULT: 0.\n");
  printf("\n");
  printf("    -r REDLA    Variant of linear algebra to be used (usually the matrix M\n");
  printf("                is splied to AB (upper part) and CD (lower part)):\n");
  printf("                  1 -> GBLA matrices, complete reduction of matrix.\n");
  printf("                  2 -> GBLA matrices, only reducing CD. (BROKEN)\n");
  printf("                  3 -> Block implementation, only reducing CD.\n");
  printf("                  4 -> Block implementation, only reducing CD, constructing\n");
  printf("                       AB only blockwise depending on block size.\n");
  printf("                  5 -> Block implementation, complete reduction of matrix.\n");
  printf("                  6 -> Sparse row implementation, ABCD mapping,\n");
  printf("                       directly reducing CD.\n");
  printf("                  7 -> Sparse row implementation, ABCD mapping, first\n");
  printf("                       reducing CD on its own, then reducing via AB\n");
  printf("                  8 -> Sparse row implementation, ABCD mapping, directly\n");
  printf("                       reducing CD, using multiline structure for AB.\n");
  printf("                  9 -> Sparse row implementation, complete reduction\n");
  printf("                       of matrix.\n");
  printf("                 10 -> Sparse row implementation, no column remapping,\n");
  printf("                       directly reducing lower matrix part.\n");
  /* printf("                 11 -> Sparse row implementation, ABCD mapping,\n");
   * printf("                       directly reducing lower matrix part, but storing\n");
   * printf("                       all pivots in a big sparse matrix\n"); */
  printf("                 42 -> Probabilistic linear algebra, error probability is\n");
  printf("                       lower than 1/(characteristic of field). See also the\n");
  printf("                       paper \"An Algorithm For Splitting Polynomial Systems\n");
  printf("                       Based on F4\" by Michael Monagan & Roman Pearce.\n");
  printf("                666 -> A combination of (4) and (6) depending on the density:\n");
  printf("                       If the density of the matrix is < 0.01%% or CD has very\n");
  printf("                       few rows option (6) is used; otherwise (4) is used.\n");
  printf("                DEFAULT: 666.\n");
  printf("\n");
  /* printf("    -s SIMP     Use simplify in F4 algorithm.\n");
   * printf("                0 -> not simplified.\n");
   * printf("                1 -> simplified but B not fully reduced.\n");
   * printf("                2 -> simplified and B fully reduced.\n");
   * printf("                NOTE: If simplification is enabled the GBLA linear algebra\n");
   * printf("                      variant (i.e. option \"-r1\") is used.\n");
   * printf("                DEFAULT: 0.\n");
   * printf("\n"); */
  printf("    -t THRDS    Number of threads used.\n");
  printf("                DEFAULT: 1.\n");
  printf("\n");
  printf("    -v VERBOSE  Verbose output during computations.\n");
  printf("\n");

  return;
}

/*  timing structs */
struct timeval t_load_start;
struct timeval t_complete;
double t_linear_algebra         = 0;
double t_linear_algebra_local   = 0;
double t_symbolic_preprocessing = 0;
double t_sorting_columns        = 0;
double t_generating_gbla_matrix = 0;
double t_update_pairs           = 0;
unsigned long n_zero_reductions = 0;
int steps = 0;

int main(int argc, char *argv[])
{
  const char *fn      = NULL;
  int verbose         = 0;
  int nthreads        = 1;
  int linear_algebra  = 666;
  /* int simplify        = 0; */
  /* int generate_pbm    = 1; */
  int print_gb        = 0;
  int order           = 0;
  long max_spairs     = 0;
  nelts_t block_size  = 4096;
  ht_size_t htc       = 18;
  int git_hash        = 0;
  /* generate file name holder if pbms are generated */
  /* char *pbm_dir = NULL;
   * char pbm_fn[400]; */

  /* monomial order names storage, for printing purpose only */
  char orders[10][10];
  snprintf(orders[0], 10, "DRL");
  snprintf(orders[1], 10, "LEX");

  int index;
  int opt;

  /* keep track of meta data, meta_data is global to be used wherever we need it
   * to keep track of data */
  meta_data = init_meta_data();

	opterr  = 0;

  while ((opt = getopt(argc, argv, "b:c:ghl:m:o:p:r:s:t:v")) != -1) {
    switch (opt) {
      case 'b':
        block_size  = (nelts_t)strtol(optarg, NULL, 10);
        break;
      case 'c':
        htc = (ht_size_t)strtol(optarg, NULL, 10);
        break;
      case 'g':
        git_hash  = 1;
        break;
      case 'h':
        print_help();
        return 0;
      case 'l':
        max_spairs  = (int)strtol(optarg, NULL, 10);
        break;
      case 'm':
        /* pbm_dir       = optarg;
         * generate_pbm  = 1; */
        break;
      case 'o':
        order  = (int)strtol(optarg, NULL, 10);
        break;
      case 'p':
        print_gb = (int)strtol(optarg, NULL, 10);
        if (print_gb > 3)
          print_gb = 0;
        break;
      case 'r':
        linear_algebra = (int)strtol(optarg, NULL, 10);
        break;
      /* case 's':
       *   simplify = (int)strtol(optarg, NULL, 10);
       *   if (simplify > 2)
       *     simplify = 1;
       *   break; */
      case 't':
        nthreads  = (int)strtol(optarg, NULL, 10);
        if (nthreads == 0)
          nthreads  = 1;
        break;
      case 'v':
        verbose = 1;
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

  /* if GBLA linear algebra is used the block size is given by a macro */
  if (linear_algebra == 1 || linear_algebra == 2)
    block_size  = __GBLA_SIMD_BLOCK_SIZE;

  /* set block_size to zero if a sparse row implementation is used */
  if (linear_algebra > 5 && linear_algebra != 666)
    block_size  = 0;

  for (index = optind; index < argc; index++)
    fn = argv[index];

  if (fn == NULL) {
    fprintf(stderr, "File name is required.\nSee help using '-h' option.\n");
    return 1;
  }
  /* we start with hash table sizes being non-Mersenne primes < 2^18, nothing
   * smaller. So we shift htc to the corresponding index of ht->primes. */
  if (htc < 15)
    htc =   15;

  /* first get number of variables in order to initialize hash table */
  nvars_t nvars = get_nvars(fn);

  /* initialize hash table */
  ht = init_hash_table(htc, nvars);
  set_sort_functions_depending_on_monomial_order(ht, order);

  if (verbose > 0) {
    if (git_hash == 1) {
      printf("git commit: ");
      FILE *fp;
      char hash_string[200];
      fp  = popen("git rev-parse HEAD", "r");
      if (fp == NULL) {
        printf("Failed to run command\n");
        exit(1);
      }
      while (fgets(hash_string, sizeof(hash_string)-1, fp) != NULL)
        printf("%s", hash_string);
      pclose(fp);
    }
    printf("---------------------------------------------------------------------------\n");
    printf("----------------------------- Computing Groebner --------------------------\n");
    printf("--------------------- with the following options set ----------------------\n");
    printf("---------------------------------------------------------------------------\n");
    printf("number of threads           %15d\n", nthreads);
    printf("hash table size             %15u (2^%u)\n", ht->sz, htc);
    printf("linear algebra variant      %15d\n", linear_algebra);
    printf("limit for handling spairs?  %15ld\n", max_spairs);
    printf("block size of matrix tiles  %15d\n", block_size);
    /* printf("use simplify?               %15d\n", simplify); */
    /* printf("generate pbm files?         %15d\n", generate_pbm); */
    printf("print resulting basis?      %15d\n", print_gb);
    printf("---------------------------------------------------------------------------\n");
  }
  /* input stores input data */
  gb_t *basis = load_input(fn, nvars, order, ht, max_spairs, verbose);
  lms = (hash_t *)malloc(basis->load * sizeof(hash_t));
  for (size_t i = 1; i < basis->st; ++i) {
    lms[i]  = basis->p[i][2];
  }

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
    /* See note on gb_t in src/types.h why we decrement basis->load here. */
    printf("number of generators        %15d\n", basis->load-1);
    printf("homogeneous input?          %15d\n", basis->init_hom);
    printf("homogeneous computation?    %15d\n", basis->hom);
    printf("input file size             %18.2f %s\n", basis->fs, basis->fsu);
    printf("---------------------------------------------------------------------------\n");
  }

  /* track time for the complete reduction process (excluding load) */
  if (verbose > 0) {
    gettimeofday(&t_complete, NULL);
    gettimeofday(&t_load_start, NULL);
  }
  /* initialize spair set */
  ps_t *ps = initialize_pair_set(basis);

  if (verbose > 0)
    t_update_pairs  +=  walltime(t_load_start);

  /* matrix handlers for linear algebra */
  pre_t *mon  = (pre_t *)malloc(sizeof(pre_t));
  mon->hash   = NULL;
  smc_t *AB   = (smc_t *)malloc(sizeof(smc_t));
  AB->row     = NULL;
  smc_t *CD   = (smc_t *)malloc(sizeof(smc_t));
  CD->row     = NULL;
  AB->mod = CD->mod = basis->mod;

  /* run while there exist spairs to be handled */
  while (ps->load > 0)
  {
    steps++;

    mon->nlm  = 0;
    mon->load = 0;
    mon->size = 2 * ps->size;
    mon->hash = realloc(mon->hash, 2 * ps->size * sizeof(hash_t));
    
    /* select next bunch of spairs */
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);

    symbolic_preprocessing(ps, AB, CD, mon, basis);

    /* set matrix meta data */
    AB->row = realloc(AB->row, AB->rk * sizeof(src_t *));
    AB->nr  = AB->rk;
    CD->row = realloc(CD->row, CD->rk * sizeof(src_t *));
    CD->nr  = CD->rk;
    AB->ncl = CD->ncl  = mon->nlm;
    AB->ncr = CD->ncr  = mon->load - mon->nlm;

    mon->hash = realloc(mon->hash, mon->load * sizeof(hash_t));

    if (verbose > 0) {
      t_symbolic_preprocessing +=  walltime(t_load_start);
      gettimeofday(&t_load_start, NULL);
    }

#if GB_DEBUG
    printf("# lead terms before sorting: %u\n", mon->nlm);
    for (int i = 0; i < mon->load; ++i) {
      if (i == mon->nlm)
        printf("-------------------------\n");
      for (int j = 0; j < ht->nv; ++j) {
        printf("%u ", ht->exp[mon->hashs[i]][j]);
      }
      printf("|| %lu | %u\n", mon->hash[i], ht->idx[mon->hash[i]]);
    }
#endif


    /* next we have to store arrays for the connection between lead monomials
     * and matrix column indices resp. non-lead monomials and matrix column
     * indices. Note that we already know the sizes of the arrays due to
     * symbolic preprocessing:
     * We first sort spd->col via lead and non lead monomials, i.e. ht->idx[i] =
     * 1 or = 2 */
    if (linear_algebra != 10) {
      sort_columns_by_lead(mon);

      /* next we can sort both parts (we know the number of lead monomials =
       * spd->sel->load - spd->sel->nsp, i.e. all polynomials considered in
       * symbolic preprocessing minus the number of spairs)
       *
       * if we do not keep A we use for A and C a row structure and thus do
       * not invert the order of the left side columns.
       * if we keep A we use for A and C a block structure and we invert
       * the order of the left side columns since gbla expects blocks be
       * inverted for a faster reduction of A in the first step of the matrix
       * reduction. */
      if (linear_algebra > 1)
        ht->sort.sort_presorted_columns(mon, nthreads);
      else
        ht->sort.sort_presorted_columns_invert_left_side(mon, nthreads);
    } else { 
        ht->sort.sort_columns(mon);
    }
    /* for (size_t i = 0; i < mon->nlm; ++i) {
     *   printf("%u - %u -- ", i, mon->hash[i]);
     *   for (size_t j = 0; j < ht->nv; ++j)
     *     printf("%u ", ht->exp[mon->hash[i]][j]);
     *   printf("\n");
     * }
     * printf("-------------------------------\n");
     * for (size_t i = mon->nlm; i < mon->load; ++i) {
     *   printf("%u - %u -- ", i, mon->hash[i]);
     *   for (size_t j = 0; j < ht->nv; ++j)
     *     printf("%u ", ht->exp[mon->hash[i]][j]);
     *   printf("\n");
     * } */


    set_column_indices_in_ht_idx(ht, mon);
    /* exchange hash values with column indices in preconstructed matrices */
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp single
      {
        #pragma omp task
        {
          if (AB->rk != 0) {
            set_column_indices(AB, ht);
            sort_rows_by_column_indices(AB);
          }
        }
        #pragma omp task
        {
          set_column_indices(CD, ht);
          sort_rows_by_column_indices(CD);
        }
        #pragma omp taskwait
      }
    }

    /* sort rows in AB and CD */
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp single
      {
        #pragma omp task
        qsort(AB->row, AB->rk, sizeof(src_t *), cmp_rows_by_decreasing_lm);
        /* ht->sort.sort_rows_by_decreasing_lm(AB); */
        #pragma omp task
        qsort(CD->row, CD->rk, sizeof(src_t *), cmp_rows_by_increasing_lm);
        /* ht->sort.sort_rows_by_increasing_lm(CD); */
        #pragma omp taskwait
      }
    }
    /* for (size_t l = 0; l < AB->rk; ++l) {
     *   printf("AB->row[%u] %p   ", l, AB->row[l]);
     *   for (size_t k = 0; k < AB->row[l][1]; k = k+2) {
     *     printf("%u %u  ", AB->row[l][k], AB->row[l][k+1]);
     *   }
     *   printf("\n");
     * }
     * printf("\n");
     * for (size_t l = 0; l < CD->rk; ++l) {
     *   printf("CD->row[%u] %p   ", l, CD->row[l]);
     *   printf("test\n");
     *   for (size_t k = 0; k < CD->row[l][1]; k = k+2) {
     *     printf("%u %u  ", CD->row[l][k], CD->row[l][k+1]);
     *   }
     *   printf("\n");
     * } */

    uint64_t terms  = 0;
    for (nelts_t i = 0; i < AB->nr; ++i) {
      terms +=  (AB->row[i][1]-2)/2;
    }
    for (nelts_t i = 0; i < CD->nr; ++i) {
      terms +=  (CD->row[i][1]-2)/2;
    }
    uint64_t dimension  =
      (uint64_t)((AB->nr+CD->nr) * (CD->ncl + CD->ncr));
    double density      = (double)terms / (double)dimension;

    /* printf("ht->load %u | AB->nr %u\n", ht->load, AB->nr); */
    /* ht =  clear_hash_table_after_symbolic_preprocessing_new(ht, mon, ps, basis, AB->nr); */
    /* clear_hash_table_after_symbolic_preprocessing(ht, mon, ps, basis, AB->nr); */

    /* find corresponding linear algebra implementation */
    switch (linear_algebra) {
   
/*       case 1:
 *         linear_algebra_gbla(basis, AB, CD, density, ps,
 *             0, verbose, nthreads);
 *         break;
 *
 *       case 2:
 *         linear_algebra_gbla(basis, AB. CD, density, ps,
 *             1, verbose, nthreads);
 *         break;
 *
 *       case 3:
 *         linear_algebra_block_ABCD_reduce_CD_directly(
 *             basis, ps, AB, CD, density, block_size, verbose, nthreads);
 *         break;
 *
 *       case 4:
 *         linear_algebra_block_ABCD_reduce_CD_directly_blockwise_AB_construction(
 *             basis, spd, density, ps, block_size, verbose, nthreads);
 *         break;
 *
 *       case 5:
 *         linear_algebra_block_ABCD_reduce_AB_first(
 *             basis, spd, density, ps, block_size, verbose, nthreads);
 *         break;
 *
 *       case 6:
 *         linear_algebra_sparse_rows_ABCD(
 *             basis, spd, density, ps, verbose, nthreads);
 *         break;
 *
 *       case 7:
 *         linear_algebra_sparse_rows_ABCD_reduce_CD_first(
 *             basis, spd, density, ps, verbose, nthreads);
 *         break;
 *
 *       case 8:
 *         linear_algebra_sparse_rows_ABCD_multiline_AB(
 *             basis, spd, density, ps, verbose, nthreads);
 *         break;
 *
 *       case 9:
 *         linear_algebra_sparse_rows_ABCD_reduce_AB_first(
 *             basis, spd, density, ps, verbose, nthreads);
 *         break;
 *
 *       case 10:
 *         linear_algebra_sparse_rows_no_column_mapping(
 *             basis, spd, density, ps, verbose, nthreads);
 *         break; */

      case 11:
        linear_algebra_all_pivs_16_bit(
            basis, AB, CD, mon, density, ps, verbose, nthreads);
        break;

      case 42:
        linear_algebra_probabilistic(
            basis, AB, CD, mon, density, ps, verbose, nthreads);
        break;

      /* case 666:
       *   if (density < 0.01 || spd->sell->load < 40) {
       *     linear_algebra_sparse_rows_ABCD(
       *         basis, spd, density, ps, verbose, nthreads);
       *   } else {
       *     linear_algebra_block_ABCD_reduce_CD_directly_blockwise_AB_construction(
       *         basis, spd, density, ps, block_size, verbose, nthreads);
       *   }
       *   break; */


      default:
        fprintf (stderr,
            "Unknown linear algebra option.\n");
        return 1;
    }
    /* ht =  clear_hash_table_after_symbolic_preprocessing_new(ht, mon, ps, basis, AB->nr); */
  /* free(mon->hash); */
  }

  free(mon->hash);
  free(mon);
  free(AB->row);
  free(AB);
  free(CD->row);
  free(CD);

  /* generate final basis for output data */
  if (verbose > 0 || print_gb > 0) {
    final_basis_for_output(&basis);
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
#if COUNT_DIV_HITS
    printf("Div Hits                          %9lu/%9lu (%4.2f)\n", meta_data->non_div_found,
        meta_data->non_div, (double)(meta_data->non_div_found)/(double)(meta_data->non_div));
#endif

    uint64_t nterms = 0;
    for (size_t l = basis->st; l < basis->fl; ++l) {
      nterms  +=  basis->p[l][1];
    }

    printf("Size of basis                     %9u\n", basis->fl-basis->st);
    printf("Number of terms in basis          %9lu (<= 2^%u)\n", nterms,
        (unsigned int)(ceil(log(nterms) / log(2))));
    printf("criteria applications (total)     %9u\n", meta_data->ncrit_total);
    printf("Number of zero reductions         %9lu\n", n_zero_reductions);
    printf("Number of hashed elements         %9u (<= 2^%u) of 2^%u\n", ht->load,
        (unsigned int)(ceil(log(ht->load) / log(2))), (unsigned int)(log(ht->sz) / log(2)));
    printf("---------------------------------------------------------------------------\n");
    if (verbose > 2)
      print_mem_usage();
  }

  /* printing of output */
  switch (print_gb) {
    case 0:
      break;
    case 1:
      print_basis(basis);
      break;
    case 2:
      print_basis_in_singular_format(basis);
      break;
    case 3:
      print_lead_ideal(basis);
      break;
    default:
      break;
  }

  /* free allocated memory */
  free(meta_data);
  free_pair_set(&ps);
  free_basis(&basis);
  free_hash_table(&ht);

  return 0;
}

#if 0
void linear_algebra_gbla(gb_t *basis, smc_t *AB, smc_t *CD,
    const double density, ps_t *ps, const int keep_A,
    const int verbose, const int nthreads)
{
  int done;
  meta_data->mat_rows = AB->nr + CD->nr;
  meta_data->mat_cols = AB->ncl + AB->ncr;
  if (verbose > 0) {
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }

  mat_t *mat  = NULL;
  if (keep_A == 1)
    mat = generate_gbla_matrix(AB, CD, nthreads);
  else
    mat = generate_gbla_matrix(AB, CD, nthreads);

  /* delete sparse matrix structures */
  for (size_t i = 0; i < AB->nr; ++i) {
    free(AB->row[i]);
    AB->row[i]  = 0;
  }
  for (size_t i = 0; i < CD->nr; ++i) {
    free(CD->row[i]);
    CD->row[i]  = 0;
  }

  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }

  /* generate pbm files of gbla matrix
   * this feature is broken at the moment */
  /* if (generate_pbm) {
   *   int pos = 0;
   *   pos = snprintf(pbm_fn+pos, 200, "%s/",pbm_dir);
   *   pos = snprintf(pbm_fn+pos, 200, "%s-mat%u.pbm", fn, steps);
   *   printf("%s\n", pbm_fn);
   *   write_matrix_to_pbm(mat, pbm_fn);
   * } */

  ri_t rankDR  = 0;
  if (keep_A == 1)
    rankDR  = reduce_gbla_matrix_keep_A(mat, verbose, nthreads);
  else
    rankDR  = reduce_gbla_matrix(mat, verbose, nthreads);
  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    n_zero_reductions +=  (CD->nr - rankDR);
    printf("%6u new %6u zero ", rankDR, CD->nr - rankDR);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d ", density);
    printf("%6u bs\n", __GBLA_SIMD_BLOCK_SIZE);
    gettimeofday(&t_load_start, NULL);
  }

  /* generate pbm files of gbla matrix
   * this feature is broken at the moment */
  /* if (generate_pbm) {
   *   int pos = 0;
   *   pos = snprintf(pbm_fn+pos, 200, "%s/",pbm_dir);
   *   pos = snprintf(pbm_fn+pos, 200, "%s-mat%u-red.pbm", fn, steps);
   *   printf("%s\n", pbm_fn);
   *   write_reduced_matrix_to_pbm(mat, pbm_fn);
   * } */

  /* add new elements to basis and update pair set */
  done  = update_basis(basis, ps, mat, ht, rankDR);

  if (verbose > 0) {
    t_update_pairs  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }

  free_gbla_matrix(&mat);
  clear_hash_table_idx(ht);

  // if we are done then we have found the constant 1 as element in the basis
  if (done) {
    basis->has_unit = 1;
  }
}

void linear_algebra_sparse_rows_no_column_mapping(
    gb_t *basis, const spd_t *spd, const double density,
    ps_t *ps, const int verbose, const int nthreads)
{
  int done;
  meta_data->mat_rows = spd->selu->load + spd->sell->load;
  meta_data->mat_cols = spd->col->load;

  nelts_t init_rk_CD  = 0;
  smc_t *AB = NULL, *CD = NULL;
  /* generate lower matrix, i.e. unkown pivots */
  CD = generate_sparse_compact_matrix_true_columns(basis, spd->sell,
      spd->col->load, nthreads);
  init_rk_CD  = CD->nr;
  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }
  /* check appearing columns in CD */
  if (spd->selu->load > 0) {
    AB = generate_sparse_compact_matrix_offset_true_columns(basis, spd->selu,
        spd->col->load, nthreads);
    if (verbose > 0)
      t_generating_gbla_matrix  +=  walltime(t_load_start);
#pragma omp parallel for num_threads(nthreads)
    for (nelts_t i=0; i<CD->nr; ++i) {
      CD->row[i]  = reduce_lower_by_upper_rows_offset_true_columns(
          CD->row[i], AB, basis);
    }
    for (nelts_t k=0; k<AB->nr; ++k) {
      free(AB->row[k]);
    }
    free(AB->row);
    free(AB);
    AB  = NULL;
  }
#if COL_CHECK
  free(columns);
#endif

  nelts_t ctr = 0;
  for (nelts_t i=0; i<CD->nr; ++i) {
    if (CD->row[i] != NULL) {
      CD->row[ctr] = CD->row[i];
      ctr++;
    }
  }
  CD->nr  = ctr;
  CD->rk  = ctr;
  if (CD->rk > 1)
    reduce_lower_rows_c(CD, 0, nthreads);
  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    n_zero_reductions +=  (spd->sell->load - CD->rk);
    printf("%6u new %6u zero ", CD->rk, spd->sell->load - CD->rk);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d\n", density);
    gettimeofday(&t_load_start, NULL);
  }

  done  = update_basis_new(basis, ps, spd, CD, ht);
  if (verbose > 0) {
    n_zero_reductions +=  (init_rk_CD - CD->rk);
  }
  clear_hash_table_idx(ht);
  for (nelts_t k=0; k<CD->rk; ++k) {
    free(CD->row[k]);
  }
  free(CD->row);
  free(CD);
  CD  = NULL;
  if (verbose > 0)
    t_update_pairs  +=  walltime(t_load_start);
  if (done) {
    basis->has_unit = 1;
  }
}

void linear_algebra_block_ABCD_reduce_CD_directly_blockwise_AB_construction(
    gb_t *basis, const spd_t *spd, const double density,
    ps_t *ps, const nelts_t block_size, const int verbose,
    const int nthreads)
{
  int done;
  meta_data->mat_rows = spd->selu->load + spd->sell->load;
  meta_data->mat_cols = spd->col->load;

  mat_gb_block_t *AB, *CD   = NULL;
  mat_gb_meta_data_t *meta  = NULL;
  if (verbose > 0)
    t_linear_algebra_local  = 0;

  /* generate meta data */
  meta  = generate_matrix_meta_data(block_size, basis->mod, spd);

  /* generate CD part */
  /* printf("gen CD block %u\n", spd->sell->load); */
  CD  = generate_mat_gb_lower(meta, basis, spd, ht, nthreads);

  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }

  if (spd->selu->load > 0) {
    nelts_t i, j;
    for (i=0; i<meta->nrb_AB; ++i) {
      /* printf("gen AB block\n"); */
      AB  = generate_mat_gb_row_block(i, meta, basis, spd, ht);
      if (verbose > 0) {
        t_generating_gbla_matrix  +=  walltime(t_load_start);
        gettimeofday(&t_load_start, NULL);
      }

      /* possibly the first block of A is already the unit matrix, then it is
       * just an empty block and we do not need to update AB at all */
      update_lower_by_upper_row_block_not_prereduced(CD, AB, i, meta, nthreads);
      if (verbose > 0) {
        t_linear_algebra_local  +=  walltime(t_load_start);
        gettimeofday(&t_load_start, NULL);
      }

      for (j=0; j<meta->ncb; ++j)
        free_mat_gb_block(AB+j);
      free(AB);
    }
  }
  smc_t *D  = convert_mat_gb_to_smc_format(CD, meta, nthreads);
  nelts_t j, k;
  for (j=0; j<meta->nrb_CD; ++j)
    for (k=0; k<meta->ncb; ++k)
      free_mat_gb_block(CD+j*meta->ncb+k);
  free(CD);
  clear_hash_table_idx(ht);

  /* get rank of D */
  nelts_t ctr = 0;
  for (nelts_t i=0; i<D->nr; ++i) {
    if (D->row[i] != NULL) {
      D->row[ctr] = D->row[i];
      ctr++;
    }
  }
  D->nr  = ctr;
  D->rk  = ctr;
  if (D->rk > 1)
    reduce_lower_rows_c(D, D->ncl, nthreads);
  if (verbose > 0) {
    t_linear_algebra_local  +=  walltime(t_load_start);
    t_linear_algebra        +=  t_linear_algebra_local;
    n_zero_reductions +=  (spd->sell->load - D->rk);
    printf("%6u new %6u zero ", D->rk, spd->sell->load-D->rk);
    printf("%9.3f sec ", t_linear_algebra_local / (1000000));
    printf("%7.3f%% d ", density);
    printf("%6u bs\n", meta->bs);
    gettimeofday(&t_load_start, NULL);
  }
  free(meta);
  done  = update_basis_new(basis, ps, spd, D, ht);
  for (nelts_t k=0; k<D->rk; ++k) {
    free(D->row[k]);
  }
  free(D->row);
  free(D);
  CD  = NULL;
  if (verbose > 0)
    t_update_pairs  +=  walltime(t_load_start);
  if (done) {
    basis->has_unit = 1;
  }
}

void linear_algebra_block_ABCD_reduce_AB_first(
    gb_t *basis, const spd_t *spd, const double density,
    ps_t *ps, const nelts_t block_size, const int verbose,
    const int nthreads)
{
  int done;
  meta_data->mat_rows = spd->selu->load + spd->sell->load;
  meta_data->mat_cols = spd->col->load;

  nelts_t i, j, k;
  mat_gb_block_t **AB       = NULL;
  mat_gb_block_t *CD        = NULL;
  mat_gb_meta_data_t *meta  = NULL;
  if (verbose > 0)
    t_linear_algebra_local  = 0;

  /* generate meta data */
  meta  = generate_matrix_meta_data(block_size, basis->mod, spd);

  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }

  if (spd->selu->load > 0) {
    /* generate AB */
    AB  = (mat_gb_block_t **)malloc(meta->nrb_AB * sizeof(mat_gb_block_t *));
    for (i=0; i<meta->nrb_AB; ++i) {
      AB[i]  = generate_mat_gb_upper_row_block(i, meta, basis, spd, ht);
    }
    if (verbose > 0) {
      t_generating_gbla_matrix  +=  walltime(t_load_start);
      gettimeofday(&t_load_start, NULL);
    }
    /* compute A^-1 B */
    reduce_upper_part(AB, meta, nthreads);
  }

  /* generate CD part */
  /* printf("gen CD block %u\n", spd->sell->load); */
  CD  = generate_mat_gb_lower(meta, basis, spd, ht, nthreads);

  /* possibly the first block of A is already the unit matrix, then it is
    * just an empty block and we do not need to update AB at all */
  /* update_lower_by_upper_row_block(CD, AB[i], i, meta, nthreads); */
  if (verbose > 0) {
    t_linear_algebra_local  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }

  if (AB != NULL) {
    update_lower_by_reduced_upper_row_block(CD, AB, meta, nthreads);
    for (i=0; j<meta->nrb_AB; ++i) {
      for (j=0; j<meta->ncb; ++j) {
        free_mat_gb_block(AB[i]+j);
      }
    }
    free(AB);
    AB  = NULL;
  }

  smc_t *D  = convert_mat_gb_to_smc_format(CD, meta, nthreads);
  for (j=0; j<meta->nrb_CD; ++j)
    for (k=0; k<meta->ncb; ++k)
      free_mat_gb_block(CD+j*meta->ncb+k);
  free(CD);

  /* get rank of D */
  nelts_t ctr = 0;
  for (nelts_t i=0; i<D->nr; ++i) {
    if (D->row[i] != NULL) {
      D->row[ctr] = D->row[i];
      ctr++;
    }
  }
  D->nr  = ctr;
  D->rk  = ctr;
  if (D->rk > 1)
    reduce_lower_rows_c(D, D->ncl, nthreads);
  if (verbose > 0) {
    t_linear_algebra_local  +=  walltime(t_load_start);
    t_linear_algebra        +=  t_linear_algebra_local;
    n_zero_reductions +=  (spd->sell->load - D->rk);
    printf("%6u new %6u zero ", D->rk, spd->sell->load - D->rk);
    printf("%9.3f sec ", t_linear_algebra_local / (1000000));
    printf("%7.3f%% d ", density);
    printf("%6u bs\n", meta->bs);
    gettimeofday(&t_load_start, NULL);
  }
  free(meta);
  done  = update_basis_new(basis, ps, spd, D, ht);
  clear_hash_table_idx(ht);
  for (nelts_t k=0; k<D->rk; ++k) {
    free(D->row[k]);
  }
  free(D->row);
  free(D);
  CD  = NULL;
  if (verbose > 0)
    t_update_pairs  +=  walltime(t_load_start);
  if (done) {
    basis->has_unit = 1;
  }
}

void linear_algebra_block_ABCD_reduce_CD_directly(
    gb_t *basis, ps_t *ps, smc_t *sAB, 
    smc_t *sCD, const double density,
    const nelts_t block_size, const int verbose,
    const int nthreads)
{
  int done;
  meta_data->mat_rows = sAB->nr + sCD->nr;
  meta_data->mat_cols = sAB->ncl + sAB->ncr;

  nelts_t j, k;
  mat_gb_block_t *AB        = NULL;
  mat_gb_block_t *CD        = NULL;
  mat_gb_meta_data_t *meta  = NULL;
  if (verbose > 0)
    t_linear_algebra_local  = 0;

  /* generate meta data */
  meta  = generate_matrix_meta_dataas(sAB, sCD, block_size);

  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }

  /* generate AB if upper rows are available */
  if (sAB->nr > 0) {
    AB  = generate_mat_gb_upper(meta, basis, spd, ht, nthreads);

    if (verbose > 0) {
      t_generating_gbla_matrix  +=  walltime(t_load_start);
      gettimeofday(&t_load_start, NULL);
    }
  }

  /* generate CD part */
  /* printf("gen CD block %u\n", spd->sell->load); */
  CD  = generate_mat_gb_lower(meta, basis, spd, ht, nthreads);

  /* possibly the first block of A is already the unit matrix, then it is
    * just an empty block and we do not need to update AB at all */
  /* update_lower_by_upper_row_block(CD, AB[i], i, meta, nthreads); */
  if (verbose > 0) {
    t_linear_algebra_local  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }

  if (AB != NULL) {
    update_lower_by_unreduced_upper_row_block(CD, AB, meta, nthreads);
  }
  for (j=0; j<meta->nrb_AB; ++j)
    for (k=0; k<meta->ncb; ++k)
      free_mat_gb_block(AB+j*meta->ncb+k);
  free(AB);
  AB  = NULL;

  smc_t *D  = convert_mat_gb_to_smc_format(CD, meta, nthreads);
  for (j=0; j<meta->nrb_CD; ++j)
    for (k=0; k<meta->ncb; ++k)
      free_mat_gb_block(CD+j*meta->ncb+k);
  free(CD);

  /* get rank of D */
  nelts_t ctr = 0;
  for (nelts_t i=0; i<D->nr; ++i) {
    if (D->row[i] != NULL) {
      D->row[ctr] = D->row[i];
      ctr++;
    }
  }
  D->nr  = ctr;
  D->rk  = ctr;
  if (D->rk > 1)
    reduce_lower_rows_c(D, D->ncl, nthreads);
  if (verbose > 0) {
    t_linear_algebra_local  +=  walltime(t_load_start);
    t_linear_algebra        +=  t_linear_algebra_local;
    n_zero_reductions +=  (spd->sell->load - D->rk);
    printf("%6u new %6u zero ", D->rk, spd->sell->load - D->rk);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d ", density);
    printf("%6u bs\n", meta->bs);
    gettimeofday(&t_load_start, NULL);
  }
  free(meta);
  done  = update_basis_new(basis, ps, spd, D, ht);
  clear_hash_table_idx(ht);
  for (nelts_t k=0; k<D->rk; ++k) {
    free(D->row[k]);
  }
  free(D->row);
  free(D);
  CD  = NULL;
  if (verbose > 0)
    t_update_pairs  +=  walltime(t_load_start);
  if (done) {
    basis->has_unit = 1;
  }
}

void linear_algebra_sparse_rows_ABCD_reduce_CD_first(
    gb_t *basis, const spd_t *spd, const double density,
    ps_t *ps, const int verbose, const int nthreads)
{
  int done;
  meta_data->mat_rows = spd->selu->load + spd->sell->load;
  meta_data->mat_cols = spd->col->load;
  if (verbose > 0) {
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }
  /* generate upper matrix, i.e. already known pivots */
  smat_t *AB = NULL, *CD = NULL;
  /* generate lower matrix, i.e. unkown pivots */
  CD = generate_sparse_matrix(basis, spd->sell,
      spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
      nthreads);
  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }
  reduce_lower_rows(CD, 0, nthreads);
  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }
  nelts_t ctr = 0;
  for (nelts_t i=0; i<CD->nr; ++i) {
    if (CD->row[i] != NULL) {
      CD->row[ctr] = CD->row[i];
      ctr++;
    }
  }
  CD->nr  = ctr;
  CD->rk  = ctr;
  if (CD->rk > 0) {
    nelts_t l;
    for (l=0; l<CD->rk; ++l) {
      if (CD->row[l]->pos[0] < CD->ncl)
        break;
    }
    if (l<CD->rk) {
      if (spd->selu->load > 0)
        AB = generate_sparse_matrix(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
      if (verbose > 0) {
        t_generating_gbla_matrix  +=  walltime(t_load_start);
        gettimeofday(&t_load_start, NULL);
      }
    }
    if (AB != NULL) {
#pragma omp parallel for num_threads(nthreads)
      for (nelts_t i=0; i<CD->nr; ++i) {
        CD->row[i]  = reduce_lower_by_upper_rows(CD->row[i], AB);
      }
      for (nelts_t k=0; k<AB->nr; ++k) {
        free(AB->row[k]->pos);
        free(AB->row[k]->val);
      }
      free(AB->row);
      free(AB);
      AB  = NULL;
    }
    ctr = 0;
    for (nelts_t i=0; i<CD->nr; ++i) {
      if (CD->row[i] != NULL) {
        CD->row[ctr] = CD->row[i];
        ctr++;
      }
    }
    CD->nr  = ctr;
    CD->rk  = ctr;
    reduce_lower_rows(CD, CD->ncl, nthreads);
  }
  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    n_zero_reductions +=  spd->sell->load - CD->rk;
    printf("%6u new %6u zero ", CD->rk, spd->sell->load - CD->rk);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d\n", density);
    gettimeofday(&t_load_start, NULL);
  }

  done  = update_basis_new_new(basis, ps, spd, CD, ht);

  for (nelts_t k=0; k<CD->rk; ++k) {
    free(CD->row[k]->pos);
    free(CD->row[k]->val);
  }
  free(CD->row);
  free(CD);
  CD  = NULL;
  clear_hash_table_idx(ht);
  if (verbose > 0) {
    t_update_pairs  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }
  if (done) {
    basis->has_unit = 1;
  }
}

void linear_algebra_sparse_rows_ABCD_unoptimized(
    gb_t *basis, const spd_t *spd, const double density,
    ps_t *ps, const int verbose, const int nthreads)
{
  /* generate upper matrix, i.e. already known pivots */
  int done;
  smat_t *AB = NULL, *CD = NULL;
  meta_data->mat_rows = spd->selu->load + spd->sell->load;
  meta_data->mat_cols = spd->col->load;
  /* generate lower matrix, i.e. unkown pivots */
  CD = generate_sparse_matrix(basis, spd->sell,
      spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
      nthreads);
  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }
  /* check appearing columns in CD */
  if (spd->selu->load > 0) {
    AB = generate_sparse_matrix(basis, spd->selu,
        spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
        nthreads);
    if (verbose > 0) {
      t_generating_gbla_matrix  +=  walltime(t_load_start);
      gettimeofday(&t_load_start, NULL);
    }
#pragma omp parallel for num_threads(nthreads)
    for (nelts_t i=0; i<CD->nr; ++i) {
      CD->row[i]  = reduce_lower_by_upper_rows(CD->row[i], AB);
    }
    for (nelts_t k=0; k<AB->nr; ++k) {
      free(AB->row[k]);
    }
    free(AB->row);
    free(AB);
    AB  = NULL;
  }
  clear_hash_table_idx(ht);

  nelts_t ctr = 0;
  for (nelts_t i=0; i<CD->nr; ++i) {
    if (CD->row[i] != NULL) {
      CD->row[ctr] = CD->row[i];
      ctr++;
    }
  }
  CD->nr  = ctr;
  CD->rk  = ctr;
  if (CD->rk > 1)
    reduce_lower_rows(CD, CD->ncl, nthreads);
  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    n_zero_reductions +=  (spd->sell->load - CD->rk);
    printf("%6u new %6u zero ", CD->rk, spd->sell->load - CD->rk);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d\n", density);
    gettimeofday(&t_load_start, NULL);
  }

  done  = update_basis_new_new(basis, ps, spd, CD, ht);
  for (nelts_t k=0; k<CD->rk; ++k) {
    free(CD->row[k]);
  }
  free(CD->row);
  free(CD);
  CD  = NULL;
  if (verbose > 0) {
    t_update_pairs  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }

  if (done) {
    basis->has_unit = 1;
  }
}

void linear_algebra_sparse_rows_ABCD(gb_t *basis, const spd_t *spd,
    const double density, ps_t *ps, const int verbose,
    const int nthreads)
{
  /* generate upper matrix, i.e. already known pivots */
  int done;
  smc_t *AB = NULL, *CD = NULL;
  meta_data->mat_rows = spd->selu->load + spd->sell->load;
  meta_data->mat_cols = spd->col->load;
  /* generate lower matrix, i.e. unkown pivots */
  CD = generate_sparse_compact_matrix(basis, spd->sell,
      spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
      nthreads);
  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }
  /* check appearing columns in CD */
  if (spd->selu->load > 0) {
    AB = generate_sparse_compact_matrix_offset(basis, spd->selu,
        spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
        nthreads);
    if (verbose > 0) {
      t_generating_gbla_matrix  +=  walltime(t_load_start);
      gettimeofday(&t_load_start, NULL);
    }
#pragma omp parallel for num_threads(nthreads)
    for (nelts_t i=0; i<CD->nr; ++i) {
      CD->row[i]  = reduce_lower_by_upper_rows_offset_c(CD->row[i], AB);
    }
    for (nelts_t k=0; k<AB->nr; ++k) {
      free(AB->row[k]);
    }
    free(AB->row);
    free(AB);
    AB  = NULL;
  }
  clear_hash_table_idx(ht);

  nelts_t ctr = 0;
  for (nelts_t i=0; i<CD->nr; ++i) {
    if (CD->row[i] != NULL) {
      CD->row[ctr] = CD->row[i];
      ctr++;
    }
  }
  CD->nr  = ctr;
  CD->rk  = ctr;
  if (CD->rk > 1)
    reduce_lower_rows_c(CD, CD->ncl, nthreads);
  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    n_zero_reductions +=  (spd->sell->load- CD->rk);
    printf("%6u new %6u zero ", CD->rk, spd->sell->load - CD->rk);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d\n", density);
    gettimeofday(&t_load_start, NULL);
  }

  done  = update_basis_new(basis, ps, spd, CD, ht);
  for (nelts_t k=0; k<CD->rk; ++k) {
    free(CD->row[k]);
  }
  free(CD->row);
  free(CD);
  CD  = NULL;
  if (verbose > 0) {
    t_update_pairs  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }

  if (done) {
    basis->has_unit = 1;
  }
}

void linear_algebra_sparse_rows_ABCD_multiline_AB(
    gb_t *basis, const spd_t *spd, const double density,
    ps_t *ps, const int verbose, const int nthreads)
{
  int done;
  meta_data->mat_rows = spd->selu->load + spd->sell->load;
  meta_data->mat_cols = spd->col->load;
  /* generate upper matrix, i.e. already known pivots */
  smc_t *AB = NULL, *CD = NULL;
  /* generate lower matrix, i.e. unkown pivots */
  CD = generate_sparse_compact_matrix(basis, spd->sell,
      spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
      nthreads);
  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
    printf("%3d %5u/%5u spairs of deg %3u %7u x %7u matrix ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }
  if (spd->selu->load > 0) {
    AB = generate_sparse_compact_multiline_matrix(basis, spd->selu,
        spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
        nthreads);
    if (verbose > 0) {
      t_generating_gbla_matrix  +=  walltime(t_load_start);
      gettimeofday(&t_load_start, NULL);
    }
#pragma omp parallel for num_threads(nthreads)
    for (nelts_t i=0; i<CD->nr; ++i) {
      CD->row[i]  = reduce_lower_by_upper_multiline_rows_c(CD->row[i], AB);
    }
    for (nelts_t k=0; k<AB->nr; ++k) {
      free(AB->row[k]);
    }
    free(AB->row);
    free(AB);
    AB  = NULL;
  }
  nelts_t ctr = 0;
  for (nelts_t i=0; i<CD->nr; ++i) {
    if (CD->row[i] != NULL) {
      CD->row[ctr] = CD->row[i];
      ctr++;
    }
  }
  CD->nr  = ctr;
  CD->rk  = ctr;
  if (CD->rk > 1)
    reduce_lower_rows_c(CD, CD->ncl, nthreads);
  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    n_zero_reductions +=  (spd->sell->load - CD->rk);
    printf("%6u new %6u zero ", CD->rk, spd->sell->load - CD->rk);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d\n", density);
    gettimeofday(&t_load_start, NULL);
  }
  done  = update_basis_new(basis, ps, spd, CD, ht);
  clear_hash_table_idx(ht);
  for (nelts_t k=0; k<CD->rk; ++k) {
    free(CD->row[k]);
  }
  free(CD->row);
  free(CD);
  CD  = NULL;
  if (verbose > 0) {
    t_update_pairs  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }
  if (done) {
    basis->has_unit = 1;
  }
}
#endif

void linear_algebra_probabilistic(gb_t *basis, smc_t *AB, smc_t *CD,
    const pre_t *mon, const double density, ps_t *ps, const int verbose,
    const int nthreads)
{
  srand((unsigned int)time(NULL));   // should only be called once
  int done;
  const nelts_t nr  = AB->nr + CD->nr;
  const nelts_t nc  = AB->ncl + AB->ncr;
  src_t **pivs      = (src_t **)malloc(nc * sizeof(src_t *));
  for (size_t i = 0; i < nc; ++i)
    pivs[i] = NULL;
  for (size_t i = 0; i < AB->rk; ++i) {
    pivs[AB->row[i][2]]     = AB->row[i];
    pivs[AB->row[i][2]][0]  = 1;
  }

  src_t *np; /* possible new pivot row */

  /* meta data information for printing */
  meta_data->mat_rows = nr;
  meta_data->mat_cols = nc;
  if (verbose > 0) {
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }

  /* global dense random row for checking after all blocks are done */
  bf_t *drg = (bf_t *)calloc(nc, sizeof(bf_t));

  /* number of blocks */
  const nelts_t nrl = CD->nr;
  const nelts_t nb  = (nelts_t)(floor(sqrt(nrl/2))) > 0 ?
    (nelts_t)(floor(sqrt(nrl/2))) :
    (nelts_t)(floor(sqrt(nrl))) ;
  nelts_t rem       = (nrl % nb == 0) ? 0 : 1;
  /* rows per block */
  const nelts_t rpb = (nrl / nb) + rem;

again:
#pragma omp parallel num_threads(nthreads) shared(drg)
  {
    src_t *mul        = (src_t *)malloc(rpb * sizeof(src_t));
    bf_t *dr          = (bf_t *)malloc(nc * sizeof(bf_t));
#pragma omp parallel for num_threads(nthreads)
    for (size_t i = 0; i < nb; ++i) {
      nelts_t nbl   = (nelts_t) (nrl > (i+1)*rpb ? (i+1)*rpb : nrl);
      nelts_t nrbl  = (nelts_t) (nbl - i*rpb); 

      /* printf("rpb %u | nrbl %u\n", rpb, nrbl); */
      if (nrbl != 0) {

        nelts_t bctr  = 0;
        while (bctr < nrbl) {
          /* fill random value array */
          for (size_t j = 0; j < nrbl; ++j)
            mul[j]  = (src_t) rand() % basis->mod;

          /* generate one dense row as random linear combination
           * of the rows of the block */
          memset(dr, 0, nc * sizeof(bf_t));
          nelts_t ctr = 0;
          for (size_t j = i*rpb; j < nbl; ++j) {
            /* printf("ctr %u\n", ctr); */
            for (size_t k = 2; k < CD->row[j][1]; k = k+2) {
              dr[CD->row[j][k]]  +=  (bf_t) mul[ctr] * CD->row[j][k+1];
            }
            ctr++;
          }

          /* reduce the dense random row w.r.t. to the already known pivots  */
          done = 0;
          while (!done) {
            np  = reduce_dense_row_by_known_pivots(
                dr, pivs, mon->load, basis->mod);
            if (!np)
              goto block_done;
            done  = __sync_bool_compare_and_swap(&pivs[np[2]], NULL, np);
          }
          bctr++;
        }
block_done:
        /* fill global dense row for final check at the end */
        for (size_t j = 0; j < nrbl; ++j)
          mul[j]  = (src_t) rand() % basis->mod;
        for (size_t j = i*rpb; j < nbl; ++j)
          for (size_t k = 2; k < CD->row[j][1]; k = k+2)
            drg[CD->row[j][k]]  +=  (bf_t) mul[j] * CD->row[j][k+1];

      }
    }
    free(mul);
    free(dr);
  }

  /* do final check, go back if check fails */
  src_t *fc  = reduce_dense_row_by_known_pivots(
      drg, pivs, mon->load, basis->mod);

  if (fc != NULL)
    goto again;

  /* interreduce new pivs */
  nelts_t nnr = 0; /* number of new rows */
  for (size_t i = nc; i > mon->nlm; --i) {
    if (pivs[i-1] != NULL) {
      ++nnr;
      memset(drg, 0, nc * sizeof(bf_t));
      for (size_t k = 2; k < pivs[i-1][1]; k = k+2)
        drg[pivs[i-1][k]]  +=  (bf_t) pivs[i-1][k+1];
      free(pivs[i-1]);
      pivs[i-1] = NULL;
      np  = reduce_dense_row_by_known_pivots(
          drg, pivs, nc, basis->mod);
      pivs[i-1] = np;
    }
  }

  free(drg);
  clear_hash_table_idx(ht);

  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    n_zero_reductions += nrl - nnr;
    printf("%5u nb ", nb);
    printf("%6u new %6u zero ", nnr, nrl - nnr);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d\n", density);
    gettimeofday(&t_load_start, NULL);
  }

  /* for (size_t i = 0; i < nc; ++i)
   *   printf("%u %p\n", i, pivs[i]); */
  done  = update_basis_all_pivs(basis, ps, mon, pivs, nc, ht);
    /* printf("1 last (%u) basis element: ", basis->load-1);
     * for (size_t i = 0; i < basis->p[basis->load-1][1]; i = i+2)
     *   printf("%u %u  ", basis->p[basis->load-1][i], basis->p[basis->load-1][i+1]);
     * printf("\n"); */

  for (size_t i = 0; i < CD->nr; ++i)
    free(CD->row[i]);
  CD->nr  = CD->rk  = 0;
    /* printf("2 last (%u) basis element: ", basis->load-1);
     * for (size_t i = 0; i < basis->p[basis->load-1][1]; i = i+2)
     *   printf("%u %u  ", basis->p[basis->load-1][i], basis->p[basis->load-1][i+1]);
     * printf("\n"); */

  for (size_t i = 0; i < nc; ++i)
    free(pivs[i]);
  free(pivs);
  pivs  = NULL;
  AB->nr  = AB->rk  = 0;
    /* printf("3 last (%u) basis element: ", basis->load-1);
     * for (size_t i = 0; i < basis->p[basis->load-1][1]; i = i+2)
     *   printf("%u %u  ", basis->p[basis->load-1][i], basis->p[basis->load-1][i+1]);
     * printf("\n"); */

  if (verbose > 0) {
    t_update_pairs  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }

  if (done) {
    basis->has_unit = 1;
  }
}

void linear_algebra_all_pivs_16_bit(gb_t *basis, smc_t *AB, smc_t *CD,
    const pre_t *mon, const double density, ps_t *ps, const int verbose,
    const int nthreads)
{
  int done;
  const nelts_t nr  = AB->nr + CD->nr;
  const nelts_t nc  = AB->ncl + AB->ncr;
  const nelts_t nrl = CD->nr;

  /* meta data information for printing */
  meta_data->mat_rows = nr;
  meta_data->mat_cols = nc;
  if (verbose > 0) {
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }

  src_t **pivs      = (src_t **)malloc(nc * sizeof(src_t *));
  for (size_t i = 0; i < nc; ++i)
    pivs[i] = NULL;
  for (size_t i = 0; i < AB->rk; ++i) {
    pivs[AB->row[i][2]]     = AB->row[i];
    pivs[AB->row[i][2]][0]  = 1;
  }

  src_t *np; /* possible new pivot row */

#pragma omp parallel num_threads(nthreads)
  {
    bf_t *dr          = (bf_t *)malloc(nc * sizeof(bf_t));
#pragma omp parallel for num_threads(nthreads)
    for (size_t i = 0; i < CD->nr; ++i) {
      memset(dr, 0, nc * sizeof(bf_t));
      for (size_t j = 2; j < CD->row[i][1]; j = j+2) {
        dr[CD->row[i][j]] = (bf_t)CD->row[i][j+1];
      }
      /* reduce the dense random row w.r.t. to the already known pivots  */
      done = 0;
      while (!done) {
        np  = reduce_dense_row_by_known_pivots(
            dr, pivs, mon->load, basis->mod);
        if (!np) {
          break;
        }
        done  = __sync_bool_compare_and_swap(&pivs[np[2]], NULL, np);
      }
    }
    free(dr);
  }

  /* interreduce new pivs */
  bf_t *dr  = (bf_t *)malloc(nc * sizeof(bf_t));
  nelts_t nnr = 0; /* number of new rows */
  for (size_t i = nc; i > mon->nlm; --i) {
    if (pivs[i-1] != NULL) {
      ++nnr;
      memset(dr, 0, nc * sizeof(bf_t));
      for (size_t k = 2; k < pivs[i-1][1]; k = k+2)
        dr[pivs[i-1][k]]  +=  (bf_t) pivs[i-1][k+1];
      free(pivs[i-1]);
      pivs[i-1] = NULL;
      np  = reduce_dense_row_by_known_pivots(
          dr, pivs, nc, basis->mod);
      pivs[i-1] = np;
    }
  }
  free(dr);
  clear_hash_table_idx(ht);

  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    n_zero_reductions += nrl - nnr;
    printf("%6u new %6u zero ", nnr, nrl - nnr);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d\n", density);
    gettimeofday(&t_load_start, NULL);
  }

  done  = update_basis_all_pivs(basis, ps, mon, pivs, nc, ht);

  for (size_t i = 0; i < CD->nr; ++i)
    free(CD->row[i]);
  CD->nr  = CD->rk  = 0;

  for (size_t i = 0; i < nc; ++i)
    free(pivs[i]);
  free(pivs);
  pivs  = NULL;
  AB->nr  = AB->rk  = 0;

  if (verbose > 0) {
    t_update_pairs  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }

  if (done) {
    basis->has_unit = 1;
  }
}

#if 0
void linear_algebra_sparse_rows_ABCD_reduce_AB_first(
    gb_t *basis, const spd_t *spd, const double density,
    ps_t *ps, const int verbose, const int nthreads)
{
  int done;
  meta_data->mat_rows = spd->selu->load + spd->sell->load;
  meta_data->mat_cols = spd->col->load;
  if (verbose > 0) {
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load,
        meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }

  smc_t *AB = NULL, *CD = NULL;
  /* generate lower matrix, i.e. unkown pivots */
  CD = generate_sparse_compact_matrix_test(basis, spd->sell,
      spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
      nthreads);
  if (verbose > 0) {
    t_generating_gbla_matrix  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }
  if (spd->selu->load > 0) {
    AB = generate_sparse_compact_matrix_offset_test(basis, spd->selu,
        spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
        nthreads);
    if (verbose > 0) {
      t_generating_gbla_matrix  +=  walltime(t_load_start);
      gettimeofday(&t_load_start, NULL);
    }
    AB  = reduce_upper_rows_c(AB);
    if (verbose > 0) {
      t_linear_algebra  +=  walltime(t_load_start);
      gettimeofday(&t_load_start, NULL);
    }
#pragma omp parallel for num_threads(nthreads)
    for (nelts_t i=0; i<CD->nr; ++i) {
      CD->row[i]  = reduce_lower_by_upper_rows_offset_c(CD->row[i], AB);
    }
    for (nelts_t k=0; k<AB->nr; ++k) {
      free(AB->row[k]);
    }
    free(AB->row);
    free(AB);
    AB  = NULL;
  }
  nelts_t ctr = 0;
  for (nelts_t i=0; i<CD->nr; ++i) {
    if (CD->row[i] != NULL) {
      CD->row[ctr] = CD->row[i];
      ctr++;
    }
  }
  CD->nr  = ctr;
  CD->rk  = ctr;
  if (CD->rk > 1)
    reduce_lower_rows_c(CD, CD->ncl, nthreads);
  if (verbose > 0) {
    t_linear_algebra  +=  walltime(t_load_start);
    n_zero_reductions +=  (spd->sell->load - CD->rk);
    printf("%6u new %6u zero ", CD->rk, spd->sell->load - CD->rk);
    printf("%9.3f sec ", walltime(t_load_start) / (1000000));
    printf("%7.3f%% d\n", density);
    gettimeofday(&t_load_start, NULL);
  }
  done  = update_basis_new(basis, ps, spd, CD, ht);
  if (verbose > 0) {
  }
  clear_hash_table_idx(ht);
  for (nelts_t k=0; k<CD->rk; ++k) {
    free(CD->row[k]);
  }
  free(CD->row);
  free(CD);
  CD  = NULL;
  if (verbose > 0) {
    t_update_pairs  +=  walltime(t_load_start);
    gettimeofday(&t_load_start, NULL);
  }
  if (done) {
    basis->has_unit = 1;
  }
}
#endif
