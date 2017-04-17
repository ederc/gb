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
  printf("                Default: 256 = 2^8.\n");
  printf("    -c HTC      Hash table cache resp. size in log_2: The size is then set\n");
  printf("                to the biggest non-Mersenne prime smaller then 2^(given value).\n");
  printf("                Default: 18.\n");
  printf("    -g GIT      Outputs git commit hash if verbosity level is >0.\n");
  printf("    -h HELP     Print help.\n");
  printf("    -l LIMIT    Maximal number of spairs handled at once.\n");
  printf("                Note: We try keep all spairs with the same lcm together,\n");
  printf("                      so it might happen that the limit is exceeded slightly.\n");
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
  const char *fn    = NULL;
  int verbose       = 0;
  int nthreads      = 1;
  int reduce_gb     = 0;
  int simplify      = 0;
  int generate_pbm  = 0;
  int print_gb      = 0;
  int order         = 0;
  long max_spairs   = 0;
  int keep_A        = 0;
  int block_size    = 256;
  ht_size_t htc     = 15;
  int git_hash      = 0;
  /* generate file name holder if pbms are generated */
  char *pbm_dir = NULL;
  char pbm_fn[400];

  /* monomial order names storage, for printing purpose only */
  char orders[10][10];
  snprintf(orders[0], 10, "DRL");
  snprintf(orders[1], 10, "LEX");

  int index;
  int opt;
  int done;



  /* keep track of meta data, meta_data is global to be used wherever we need it
   * to keep track of data */
  meta_data = init_meta_data();

	opterr  = 0;

  while ((opt = getopt(argc, argv, "b:c:ghl:m:no:p:r:s:t:v:")) != -1) {
    switch (opt) {
      case 'b':
        block_size  = (int)strtol(optarg, NULL, 10);
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
    printf("compute reduced basis?      %15d\n", reduce_gb);
    printf("do not reduce A|B in gbla   %15d\n", keep_A);
    printf("limit for handling spairs?  %15ld\n", max_spairs);
    printf("block size of matrix tiles  %15d\n", block_size);
    printf("use simplify?               %15d\n", simplify);
    printf("generate pbm files?         %15d\n", generate_pbm);
    printf("print resulting basis?      %15d\n", print_gb);
    printf("---------------------------------------------------------------------------\n");
  }
  /* input stores input data */
  gb_t *basis = load_input(fn, nvars, order, ht, simplify, max_spairs, verbose);
  set_simplify_functions(ht, basis);

  /* global simplifier list */
  /* generate simplifier list if simplification is enabled */
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

  if (verbose > 2) {
    printf("---------------------------------------------------------------------------\n");
    printf("ps->size                          %9u\n",ps->size);
    printf("ps->load                          %9u\n",ps->load);
    printf("basis->size                       %9u\n",basis->size);
    /* See note on gb_t in src/types.h why we decrement basis->load here. */
    printf("basis->load                       %9u\n",basis->load-1); 
    printf("---------------------------------------------------------------------------\n");
    printf("criteria applications (last step) %9u\n",meta_data->ncrit_last);
    printf("criteria applications (total)     %9u\n",meta_data->ncrit_total);
  }


  /* run while there exist spairs to be handled */
  while (ps->load > 0)
  {
    steps++;
    
    /* select next bunch of spairs */
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);
    spd_t *spd  = symbolic_preprocessing(ps, basis, sf);
#if HASH_CHECK
    printf("HOW OFTEN DO HASHES APPEAR?\n");
    nelts_t counter = 0;
    for (nelts_t kk=0; kk<ht->load; ++kk) {
      if (ht->ctr[kk] == 1) {
        printf("hash %lu appears %4u times\n", ht->val[kk], ht->ctr[kk]);
        counter++;
      }
    }
    printf("%u hashes appear only once!\n",counter);
#endif
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
#if GB_DEBUG
    printf("# lead terms before sorting: %u\n", spd->col->nlm);
    for (int ii=0; ii<spd->col->load; ++ii) {
      if (ii==spd->col->nlm)
        printf("-------------------------\n");
      for (int jj=0; jj<ht->nv; ++jj) {
        printf("%u ", ht->exp[spd->col->hpos[ii]][jj]);
      }
      printf("|| %lu | %u\n", spd->col->hpos[ii], ht->idx[spd->col->hpos[ii]]);
    }
#endif

    /* next we have to store arrays for the connection between lead monomials
     * and matrix column indices resp. non-lead monomials and matrix column
     * indices. Note that we already know the sizes of the arrays due to
     * symbolic preprocessing:
     * We first sort spd->col via lead and non lead monomials, i.e. ht->idx[i] =
     * 1 or = 2 */
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);

    if (reduce_gb != 111) {
      sort_columns_by_lead(spd);

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
      if (keep_A == 1 || reduce_gb > 1)
        ht->sort.sort_presorted_columns(spd, nthreads);
      else
        ht->sort.sort_presorted_columns_invert_left_side(spd, nthreads);
    } else { 
      /* reduce_gb == 111 */
      if (basis->ord == 0)
        qsort(spd->col->hpos, spd->col->load,
            sizeof(hash_t), cmp_symbolic_preprocessing_monomials_by_grevlex);
      if (basis->ord == 1)
        qsort(spd->col->hpos, spd->col->load,
            sizeof(hash_t), cmp_symbolic_preprocessing_monomials_by_lex);
    }
    /* connect monomial hash positions with columns in to be constructed gbla
     * matrix */
    set_column_index_in_hash_table(ht, spd);

    /* now sort upper and lower polynomial selections by lead monomials. this
     * corresponds to sorting by columns (=rows as this is one-to-one for lead
     * monomials).
     * if we do not keep A we use for A and C a row structure and thus do not
     * invert the order of the left side columns.
     * if we keep A we use for A and C a block structure and we invert
     * the order of the left side columns since gbla expects blocks be
     * inverted for a faster reduction of A in the first step of the matrix
     * reduction. */
    if (keep_A == 1)
      sort_selection_by_inverted_column_index(spd, nthreads);
    else
      sort_selection_by_column_index(spd, nthreads);
    if (verbose > 0)
      t_sorting_columns +=  walltime(t_load_start);

#if GB_DEBUG
    printf("# lead terms: %u\n", spd->col->nlm);
    for (int ii=0; ii<spd->col->load; ++ii) {
      if (ii==spd->col->nlm)
        printf("-------------------------\n");
      for (int jj=0; jj<ht->nv; ++jj) {
        printf("%u ", ht->exp[spd->col->hpos[ii]][jj]);
      }
      printf("|| %lu | %u\n", spd->col->hpos[ii], ht->idx[spd->col->hpos[ii]]);
    }
#endif
  
    /* generate gbla matrix out of data from symbolic preprocessing */
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);

    /*
     * FOR TESTING ONLY
     */

    /* new la implementation */
    /* only sparse rows, no block structure */
#if 1
    /* here we only check the density of the matrix, if it is too low we
     * use only sparse rows (i.e. reduce_gb = 101). else we use a block
     * version (reduce_gb = 23). */
    uint64_t terms  = 0;
    for (nelts_t ii=0; ii<spd->selu->load; ++ii) {
      /* printf("bi %u\n",spd->selu[ii].mpp->bi); */
      terms +=  basis->nt[spd->selu->mpp[ii].bi];
    }
    for (nelts_t ii=0; ii<spd->sell->load; ++ii) {
      terms +=  basis->nt[spd->sell->mpp[ii].bi];
    }
    uint64_t dimension  = (uint64_t)(spd->selu->load+spd->sell->load)*spd->col->load;
    double density      = (double)terms / (double)dimension;
    /* printf("%lu\n", dimension); */

    if (reduce_gb == 666) {
      if (density < 0.01 || spd->sell->load < 40)
        reduce_gb_101(basis, spd, density, ps, block_size, verbose, nthreads);
      else
        reduce_gb_23(basis, spd, density, ps, block_size, verbose, nthreads);
    }

    if (reduce_gb == 24) {
      mat_gb_meta_data_t *meta  = NULL;
      mat_gb_block_t **bAB      = NULL;

      nelts_t init_rk_CD  = 0;
      smc_t *CD, *AB      = NULL;

      /* generate meta data */
      meta  = generate_matrix_meta_data(block_size, basis->mod, spd);
      meta_data->mat_rows = spd->selu->load + spd->sell->load;
      meta_data->mat_cols = meta->nc_AC + meta->nc_BD;
      if (verbose > 1) {
        printf("matrix rows %6u \n", meta_data->mat_rows);
        printf("matrix cols %6u \n", meta_data->mat_cols);
      }
      if (verbose == 1) {
        printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
            steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
        fflush(stdout);
      }
      if (spd->selu->load > 0) {
        bAB =
          (mat_gb_block_t **)malloc(meta->nrb_AB * sizeof(mat_gb_block_t *));

        /* generate AB completely */
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<meta->nrb_AB; ++i) {
          bAB[i]  = generate_mat_gb_upper_row_block(i, meta, basis, spd, ht);
        }
        if (verbose > 0) {
          t_generating_gbla_matrix  +=  walltime(t_load_start);
          gettimeofday(&t_load_start, NULL);
        }
        /* interreduce AB blockwise: The first nonempty block in each
         * block row is the unit matrix after this step */
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<meta->nrb_AB; ++i) {
          if (bAB[i][i].len != NULL)
            update_upper_row_block(bAB[i], i, meta, nthreads);
        }
        if (verbose > 0) {
          t_linear_algebra_local  +=  walltime(t_load_start);
          gettimeofday(&t_load_start, NULL);
        }

        AB  = convert_mat_gb_to_smc_offset_format(bAB, meta, nthreads);
        /* printf("-- AB --\n");
         * for (int ii=0; ii<AB->nr; ++ii) {
         *   printf("row[%u] ",ii);
         *   if (AB->row[ii] == NULL)
         *     printf("NULL\n");
         *   else {
         *     printf("len %u offset %u ||| ", AB->row[ii][0], AB->row[ii][1]);
         *     for (int jj=2; jj<AB->row[ii][0]; jj += 2) {
         *       printf("%u at %u | ",AB->row[ii][jj+1],AB->row[ii][jj]);
         *     }
         *     printf("\n");
         *   }
         * } */

        for (nelts_t i=0; i<meta->nrb_AB; ++i) {
          for (nelts_t j=0; j<meta->ncb; ++j) {
            free_mat_gb_block(bAB[i]+j);
          }
          free(bAB[i]);
        }
        free(bAB);
      }
      CD = generate_sparse_compact_matrix_test(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      init_rk_CD  = CD->rk;

      if (AB != NULL) {
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
        if (verbose == 1 && steps > 1)
          printf("%9.3f sec ", walltime(t_load_start) / (1000000));
      }
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
#if 0
        printf("--CD 2--\n");
        for (int ii=0; ii<CD->nr; ++ii) {
          printf("row[%u] ",ii);
          if (CD->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
              printf("%u at %u | ",CD->row[ii][jj+1],CD->row[ii][jj]);
            }
            printf("\n");
          }
        }
#endif
      CD->nr  = ctr;
      CD->rk  = ctr;
      /* printf("rank of CD %u | %u\n", CD->rk, CD->nr); */
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%4u new %4u zero ", CD->rk, init_rk_CD-CD->rk);
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (init_rk_CD - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
#endif
    /* new la implementation */
    if (reduce_gb == 27) {
      mat_gb_block_t *AB, *CD   = NULL;
      mat_gb_meta_data_t *meta  = NULL;
      if (verbose > 0)
        t_linear_algebra_local  = 0;

      /* generate meta data */
      meta  = generate_matrix_meta_data(block_size, basis->mod, spd);

      /* generate CD part */
      CD  = generate_mat_gb_lower(meta, basis, spd, ht, nthreads);

      meta_data->mat_rows = spd->selu->load + spd->sell->load;
      meta_data->mat_cols = meta->nc_AC + meta->nc_BD;
      if (verbose > 0)
        t_generating_gbla_matrix  +=  walltime(t_load_start);
      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      if (verbose > 1) {
        printf("matrix rows %6u \n", meta_data->mat_rows);
        printf("matrix cols %6u \n", meta_data->mat_cols);
      }
      if (verbose == 1) {
        printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
            steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
        fflush(stdout);
      }

      if (spd->selu->load > 0) {
        nelts_t i, j;
        for (i=0; i<meta->nrb_AB; ++i) {
          AB  = generate_mat_gb_upper_row_block(i, meta, basis, spd, ht);
          if (verbose > 0) {
            t_generating_gbla_matrix  +=  walltime(t_load_start);
            gettimeofday(&t_load_start, NULL);
          }

          /* we do not update AB in this variant. so we directly reduce CD with the AB block */
          /* if (AB[i].len != NULL)
           *   update_upper_row_block(AB, i, meta, nthreads); */

          update_lower_by_upper_row_block(CD, AB, i, meta, nthreads);
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
#if newred
      printf("--D BEGINNING--\n");
      for (int ii=0; ii<D->nr; ++ii) {
        printf("row[%u] ",ii);
        if (D->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<D->row[ii][0]; jj += 2) {
            printf("%u at %u | ",D->row[ii][jj+1],D->row[ii][jj]);
          }
          printf("\n");
        }
      }
#endif
      nelts_t j, k;
      for (j=0; j<meta->nrb_CD; ++j)
        for (k=0; k<meta->ncb; ++k)
          free_mat_gb_block(CD+j*meta->ncb+k);
      free(CD);
      free(meta);

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
      }
      if (verbose == 1 && steps > 0) {
        printf("%4u new %4u zero ", D->rk, spd->sell->load-D->rk);
        printf("%9.3f sec\n", t_linear_algebra_local / (1000000));
      }
      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, D, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (spd->sell->load - D->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    if (reduce_gb == 23) {
      reduce_gb_23(basis, spd, density, ps, block_size, verbose, nthreads);
    }

    if (reduce_gb == 11) {
      /* in this variant we first interreduce AB, so a good comparison to gbla */


      /* we generate the matrix in the following shape:
      * AB = already known lead terms resp. pivots
      * --
      * CD = new data, new lead terms to be computed
      *
      * genearte upper matrix, i.e. already known pivots */
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix_test(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      if (spd->selu->load > 0) {
        AB = generate_sparse_compact_matrix_offset_test(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
        AB  = reduce_upper_rows_c(AB);
        /* printf("-- AB --\n");
         * for (int ii=0; ii<AB->nr; ++ii) {
         *   printf("row[%u] ",ii);
         *   if (AB->row[ii] == NULL)
         *     printf("NULL\n");
         *   else {
         *     printf("len %u offset %u ||| ", AB->row[ii][0], AB->row[ii][1]);
         *     for (int jj=2; jj<AB->row[ii][0]; jj += 2) {
         *       printf("%u at %u | ",AB->row[ii][jj+1],AB->row[ii][jj]);
         *     }
         *     printf("\n");
         *   }
         * } */

#if newred
        printf("--CD BEGINNING--\n");
        for (int ii=0; ii<CD->nr; ++ii) {
          printf("row[%u] ",ii);
          if (CD->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
              printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
            }
            printf("\n");
          }
        }
#endif
        meta_data->mat_rows = spd->selu->load + spd->sell->load;
        meta_data->mat_cols = CD->ncl + CD->ncr;
        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u pairs of deg %3u %7u x %7u matrix ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<CD->nr; ++i) {
          /* printf("CD row idx %u of %u | length %u\n", i, CD->nr, CD->row[i][0]); */
          /* printf("CD row idx %u of %u | length %u | last colum idx %u\n", i, CD->nr, CD->row[i][0], CD->row[i][CD->row[i][0]-2]); */
          CD->row[i]  = reduce_lower_by_upper_rows_offset_c(CD->row[i], AB);
        }
        for (nelts_t k=0; k<AB->nr; ++k) {
          free(AB->row[k]);
        }
        free(AB->row);
        free(AB);
        AB  = NULL;
      }
#if newred
      printf("--CD BEFORE--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
#if newred
        printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
#if newred
      printf("rank of CD %u\n", CD->rk);
      printf("--CD AFTER--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    /* here we do not remap the columns. in this way we can store A via column
     * positions only and map to the corresponding polynomials in basis for the
     * coefficients */
    if (reduce_gb == 111) {
      /* we generate the matrix in the following shape:
      * AB = already known lead terms resp. pivots
      * --
      * CD = new data, new lead terms to be computed
      *
      * genearte upper matrix, i.e. already known pivots */
      nelts_t init_rk_CD  = 0;
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix_true_columns(basis, spd->sell,
          spd->col->load, nthreads);
      init_rk_CD  = CD->nr;
      /* check appearing columns in CD */
      if (spd->selu->load > 0) {
        AB = generate_sparse_compact_matrix_offset_true_columns(basis, spd->selu,
            spd->col->load, nthreads);
        meta_data->mat_rows = spd->selu->load + spd->sell->load;
        meta_data->mat_cols = CD->ncl + CD->ncr;
        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<CD->nr; ++i) {
          CD->row[i]  = reduce_lower_by_upper_rows_offset_true_columns(CD->row[i], AB, basis);
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
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%4u new %4u zero ", CD->rk, init_rk_CD-CD->rk);
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (init_rk_CD - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    if (reduce_gb == 121) {
      /* we generate the matrix in the following shape:
       * AB = already known lead terms resp. pivots
       * --
       * CD = new data, new lead terms to be computed
       *
       * genearte upper matrix, i.e. already known pivots */
      nelts_t init_rk_CD  = 0;
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix_test(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
#if 0
      printf("--CD 1--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj+1],CD->row[ii][jj]);
          }
          printf("\n");
        }
      }
#endif
      init_rk_CD  = CD->nr;
      /* check appearing columns in CD */
      if (spd->selu->load > 0) {
        AB = generate_sparse_compact_matrix_offset_test(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
#if 0
        printf("--AB 1--\n");
        for (int ii=0; ii<AB->nr; ++ii) {
          printf("row[%u] ",ii);
          if (AB->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=2; jj<AB->row[ii][0]; jj += 2) {
              printf("%u at %u | ",AB->row[ii][jj+1],AB->row[ii][jj]);
            }
            printf("\n");
          }
        }
#endif
        meta_data->mat_rows = spd->selu->load + spd->sell->load;
        meta_data->mat_cols = CD->ncl + CD->ncr;
        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
        const nelts_t chunk = CD->nr / nthreads;
        const nelts_t leftover  = CD->nr % nthreads;
#pragma omp parallel num_threads(nthreads)
        {
#pragma omp single nowait
          {
            for (nelts_t i=0; i<CD->nr-leftover; i+=chunk) {
#pragma omp task
              {
                reduce_many_lower_by_upper_rows_offset_c(CD->row+i, chunk, AB);
              }
            }
#pragma omp task
            {
              reduce_many_lower_by_upper_rows_offset_c(CD->row+(CD->nr-leftover), leftover, AB);
            }
          }
#pragma omp taskwait
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
#if 0
      printf("--CD 2--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj+1],CD->row[ii][jj]);
          }
          printf("\n");
        }
      }
#endif
      CD->nr  = ctr;
      CD->rk  = ctr;
      /* printf("rank of CD %u | %u\n", CD->rk, CD->nr); */
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%4u new %4u zero ", CD->rk, init_rk_CD-CD->rk);
      printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (init_rk_CD - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    if (reduce_gb == 123) {
      /* we generate the matrix in the following shape:
      * AB = already known lead terms resp. pivots
      * --
      * CD = new data, new lead terms to be computed
      *
      * genearte upper matrix, i.e. already known pivots */
      nelts_t init_rk_CD  = 0;
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix_test(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
#if 0
        printf("--CD 1--\n");
        for (int ii=0; ii<CD->nr; ++ii) {
          printf("row[%u] ",ii);
          if (CD->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
              printf("%u at %u | ",CD->row[ii][jj+1],CD->row[ii][jj]);
            }
            printf("\n");
          }
        }
#endif
      meta_data->mat_rows = spd->selu->load + spd->sell->load;
      meta_data->mat_cols = CD->ncl + CD->ncr;

      init_rk_CD  = CD->nr;
      if (spd->selu->load > 0) {
        for (nelts_t k=0; k<spd->selu->load; k=k+block_size) {
          const nelts_t nr = block_size < spd->selu->load - k ? 
            block_size : spd->selu->load-k;
          /* printf("sl %u | k %u | nr %u\n", spd->selu->load, k, nr); */
          AB = generate_sparse_compact_matrix_offset_block(basis, spd->selu,
              k, nr, spd->col->nlm, spd->col->load-spd->col->nlm);
        /* printf("--AB %u--\n", k);
         * for (int ii=0; ii<AB->nr; ++ii) {
         *   printf("row[%u] ",ii);
         *   if (AB->row[ii] == NULL)
         *     printf("NULL\n");
         *   else {
         *     for (int jj=2; jj<AB->row[ii][0]; jj += 2) {
         *       printf("%u at %u | ",AB->row[ii][jj+1],AB->row[ii][jj]);
         *     }
         *     printf("\n");
         *   }
         * } */
        interreduce_upper_rows_offset_c(AB, k,  nthreads);
        /* printf("--AB %u--\n", k);
         * for (int ii=0; ii<AB->nr; ++ii) {
         *   printf("row[%u] ",ii);
         *   if (AB->row[ii] == NULL)
         *     printf("NULL\n");
         *   else {
         *     for (int jj=2; jj<AB->row[ii][0]; jj += 2) {
         *       printf("%u at %u | ",AB->row[ii][jj+1],AB->row[ii][jj]);
         *     }
         *     printf("\n");
         *   }
         * } */

#pragma omp parallel for num_threads(nthreads)
          for (nelts_t i=0; i<CD->nr; ++i) {
            /* printf("CD->row[%u] = %p\n", i, CD->row[i]); */
            if (CD->row[i]  != NULL)
              CD->row[i]  = reduce_lower_by_upper_rows_offset_c_block(CD->row[i], AB, k);
          }
          for (nelts_t k=0; k<AB->nr; ++k) {
            free(AB->row[k]);
          }
          free(AB->row);
          free(AB);
          AB  = NULL;
        }

        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
      }
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          if (CD->row[ctr][2] != 1)
            CD->row[ctr] = normalize_row_c(CD->row[ctr], CD->mod);
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
      /* printf("rank of CD %u | %u\n", CD->rk, CD->nr); */
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%4u new %4u zero ", CD->rk, init_rk_CD-CD->rk);
        printf("%9.3f sec\n", t_linear_algebra / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (init_rk_CD - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }

    if (reduce_gb == 101) {
      reduce_gb_101(basis, spd, density, ps, block_size, verbose, nthreads);
    }
    if (reduce_gb == 10) {
      /* we generate the matrix in the following shape:
      * AB = already known lead terms resp. pivots
      * --
      * CD = new data, new lead terms to be computed
      *
      * genearte upper matrix, i.e. already known pivots */
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      /* check appearing columns in CD */
#define COL_CHECK 0
#if COL_CHECK
      uint32_t *columns  = (uint32_t *)calloc(spd->col->nlm, sizeof(uint32_t));
      for (int ii=0; ii<CD->nr; ++ii) {
        for (int jj=1; jj<CD->row[ii][0]; jj=jj+2) {
          if (CD->row[ii][jj] < spd->col->nlm)
            columns[CD->row[ii][jj]] = 1;
        }
      }
      uint32_t colctr  = 0;
      for (uint32_t ii=0; ii<spd->col->nlm; ++ii) {
        if (columns[ii] != 0)
          colctr++;
      }
      printf("\n\nCD %6u / %6u cols | %4.2f %\n", colctr, spd->col->nlm,
          (float)colctr / (float)spd->col->nlm);
#endif
      if (spd->selu->load > 0) {
        AB = generate_sparse_compact_matrix_offset(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
#if COL_CHECK
        uint32_t *appears  = (uint32_t *)calloc(spd->col->nlm, sizeof(uint32_t));
        for (int ii=0; ii<AB->nr; ++ii) {
          // starting at 2 due to offsets!
          for (int jj=2; jj<AB->row[ii][0]; jj=jj+2) {
            if (AB->row[ii][jj] < spd->col->nlm && columns[AB->row[ii][jj]] == 0)
              appears[AB->row[ii][jj]]++;
          }
        }
        /* uint32_t cc=0;
         * for (uint32_t ii=0; ii<spd->col->nlm; ++ii) {
         *   if (appears[ii] != 0) {
         *     cc++;
         *     printf("[%u] column %6u appears %6u times\n", cc, ii, appears[ii]);
         *   }
         * } */
        uint32_t ac = 0;
        for (int ii=0; ii<AB->nr; ++ii) {
          // starting at 2 due to offsets!
          for (int jj=2; jj<AB->row[ii][0]; jj=jj+2) {
            if (AB->row[ii][jj] < spd->col->nlm && columns[AB->row[ii][jj]] == 0) {
              ac++;
              break;
            }
          }
        }
        printf("AB %6u / %6u rows | %4.2f %\n\n",ac, AB->nr, (float)ac/(float)AB->nr);
        free(appears);
#endif
#if newred
        printf("--CD BEGINNING--\n");
        for (int ii=0; ii<CD->nr; ++ii) {
          printf("row[%u] ",ii);
          if (CD->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
              printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
            }
            printf("\n");
          }
        }
#endif
        meta_data->mat_rows = spd->selu->load + spd->sell->load;
        meta_data->mat_cols = CD->ncl + CD->ncr;
        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u spairs of deg %3u %7u x %7u matrix ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<CD->nr; ++i) {
          /* printf("CD row idx %u of %u | length %u\n", i, CD->nr, CD->row[i][0]); */
          /* printf("CD row idx %u of %u | length %u | last colum idx %u\n", i, CD->nr, CD->row[i][0], CD->row[i][CD->row[i][0]-2]); */
          CD->row[i]  = reduce_lower_by_upper_rows_offset_c(CD->row[i], AB);
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
#if newred
      printf("--CD BEFORE--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
#if newred
        printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
#if newred
      printf("rank of CD %u\n", CD->rk);
      printf("--CD AFTER--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    if (reduce_gb == 9) {
      /* we generate the matrix in the following shape:
      * AB = already known lead terms resp. pivots
      * --
      * CD = new data, new lead terms to be computed
      *
      * generate upper matrix, i.e. already known pivots */
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      if (spd->selu->load > 0) {
        AB = generate_sparse_compact_matrix_pos_val(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
#if newred
        printf("--CD BEGINNING--\n");
        for (int ii=0; ii<CD->nr; ++ii) {
          printf("row[%u] ",ii);
          if (CD->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
              printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
            }
            printf("\n");
          }
        }
#endif
        meta_data->mat_rows = spd->selu->load + spd->sell->load;
        meta_data->mat_cols = CD->ncl + CD->ncr;
        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u spairs of deg %3u %7u x %7u matrix ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<CD->nr; ++i) {
          /* printf("CD row idx %u of %u | length %u\n", i, CD->nr, CD->row[i][0]); */
          /* printf("CD row idx %u of %u | length %u | last colum idx %u\n", i, CD->nr, CD->row[i][0], CD->row[i][CD->row[i][0]-2]); */
          CD->row[i]  = reduce_lower_by_upper_rows_pos_val_2_c(CD->row[i], AB);
        }
        for (nelts_t k=0; k<AB->nr; ++k) {
          free(AB->row[k]);
        }
        free(AB->row);
        free(AB);
        AB  = NULL;
      }
#if newred
      printf("--CD BEFORE--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
#if newred
        printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
#if newred
      printf("rank of CD %u\n", CD->rk);
      printf("--CD AFTER--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    if (reduce_gb == 8) {
      /* we generate the matrix in the following shape:
      * AB = already known lead terms resp. pivots
      * --
      * CD = new data, new lead terms to be computed
      *
      * generate upper matrix, i.e. already known pivots */
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      if (spd->selu->load > 0) {
        AB = generate_sparse_compact_matrix_new(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
#if newred
        printf("--CD BEGINNING--\n");
        for (int ii=0; ii<CD->nr; ++ii) {
          printf("row[%u] ",ii);
          if (CD->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
              printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
            }
            printf("\n");
          }
        }
#endif
        meta_data->mat_rows = spd->selu->load + spd->sell->load;
        meta_data->mat_cols = CD->ncl + CD->ncr;
        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u spairs of deg %3u %7u x %7u matrix ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<CD->nr; ++i) {
          /* printf("CD row idx %u of %u | length %u\n", i, CD->nr, CD->row[i][0]); */
          /* printf("CD row idx %u of %u | length %u | last colum idx %u\n", i, CD->nr, CD->row[i][0], CD->row[i][CD->row[i][0]-2]); */
          CD->row[i]  = reduce_lower_by_upper_rows_new_c(CD->row[i], AB);
        }
        for (nelts_t k=0; k<AB->nr; ++k) {
          free(AB->row[k]);
        }
        free(AB->row);
        free(AB);
        AB  = NULL;
      }
#if newred
      printf("--CD BEFORE--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
#if newred
        printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
#if newred
      printf("rank of CD %u\n", CD->rk);
      printf("--CD AFTER--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    if (reduce_gb == 7) {
      /* we generate the matrix in the following shape:
      * AB = already known lead terms resp. pivots
      * --
      * CD = new data, new lead terms to be computed
      *
      * generate upper matrix, i.e. already known pivots */
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      if (spd->selu->load > 0) {
        AB = generate_sparse_compact_multiline_matrix(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
#if newred
        printf("--CD BEGINNING--\n");
        for (int ii=0; ii<CD->nr; ++ii) {
          printf("row[%u] ",ii);
          if (CD->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
              printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
            }
            printf("\n");
          }
        }
#endif
        meta_data->mat_rows = spd->selu->load + spd->sell->load;
        meta_data->mat_cols = CD->ncl + CD->ncr;
        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u spairs of deg %3u %7u x %7u matrix ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<CD->nr; ++i) {
          /* printf("CD row idx %u of %u | length %u\n", i, CD->nr, CD->row[i][0]); */
          /* printf("CD row idx %u of %u | length %u | last colum idx %u\n", i, CD->nr, CD->row[i][0], CD->row[i][CD->row[i][0]-2]); */
          CD->row[i]  = reduce_lower_by_upper_multiline_rows_c(CD->row[i], AB);
        }
        for (nelts_t k=0; k<AB->nr; ++k) {
          free(AB->row[k]);
        }
        free(AB->row);
        free(AB);
        AB  = NULL;
      }
#if newred
      printf("--CD BEFORE--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
#if newred
        printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
#if newred
      printf("rank of CD %u\n", CD->rk);
      printf("--CD AFTER--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    if (reduce_gb == 6) {
      /* we generate the matrix in the following shape:
      * AB = already known lead terms resp. pivots
      * --
      * CD = new data, new lead terms to be computed
      *
      * generate upper matrix, i.e. already known pivots */
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      if (spd->selu->load > 0) {
        AB = generate_sparse_compact_matrix(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
#if newred
        printf("--CD BEGINNING--\n");
        for (int ii=0; ii<CD->nr; ++ii) {
          printf("row[%u] ",ii);
          if (CD->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
              printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
            }
            printf("\n");
          }
        }
#endif
        meta_data->mat_rows = spd->selu->load + spd->sell->load;
        meta_data->mat_cols = CD->ncl + CD->ncr;
        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u spairs of deg %3u %7u x %7u matrix ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
        nelts_t i = 0;
#pragma omp parallel for num_threads(nthreads)
        for (i=0; i<CD->nr-1; i=i+2) {
          /* printf("CD row idx %u of %u | length %u\n", i, CD->nr, CD->row[i][0]); */
          /* printf("CD row idx %u of %u | length %u | last colum idx %u\n", i, CD->nr, CD->row[i][0], CD->row[i][CD->row[i][0]-2]); */
          reduce_lower_by_upper_rows_double_c(&(CD->row[i]), &(CD->row[i+1]), AB);
        }
        if (CD->nr%2 == 1) {
          CD->row[CD->nr-1]  = reduce_lower_by_upper_rows_c(CD->row[CD->nr-1], AB);
        }

        for (nelts_t k=0; k<AB->nr; ++k) {
          free(AB->row[k]);
        }
        free(AB->row);
        free(AB);
        AB  = NULL;
      }
#if newred
      printf("--CD BEFORE--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
#if newred
        printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
#if newred
      printf("rank of CD %u\n", CD->rk);
      printf("--CD AFTER--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    if (reduce_gb == 2) {
      /* we generate the matrix in the following shape:
      * AB = already known lead terms resp. pivots
      * --
      * CD = new data, new lead terms to be computed
      *
      * generate upper matrix, i.e. already known pivots */
      smc_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_compact_matrix(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      if (spd->selu->load > 0) {
        AB = generate_sparse_compact_matrix(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
        if (simplify > 0) {
          interreduce_pivots_c(AB->row, AB->nr, AB->ncl, AB->ncr, AB->mod);
          update_simplifiers_new(basis, sf, spd, AB, ht);
        }
#if newred
        printf("--CD BEGINNING--\n");
        for (int ii=0; ii<CD->nr; ++ii) {
          printf("row[%u] ",ii);
          if (CD->row[ii] == NULL)
            printf("NULL\n");
          else {
            for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
              printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
            }
            printf("\n");
          }
        }
#endif
        meta_data->mat_rows = spd->selu->load + spd->sell->load;
        meta_data->mat_cols = CD->ncl + CD->ncr;
        if (verbose > 0)
          t_generating_gbla_matrix  +=  walltime(t_load_start);
        if (verbose > 0)
          gettimeofday(&t_load_start, NULL);
        if (verbose > 1) {
          printf("matrix rows %6u \n", meta_data->mat_rows);
          printf("matrix cols %6u \n", meta_data->mat_cols);
        }
        if (verbose == 1) {
          printf("%3d %5u/%5u spairs of deg %3u %7u x %7u matrix ",
              steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
          fflush(stdout);
        }
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<CD->nr; ++i) {
          /* printf("CD row idx %u of %u | length %u\n", i, CD->nr, CD->row[i][0]); */
          /* printf("CD row idx %u of %u | length %u | last colum idx %u\n", i, CD->nr, CD->row[i][0], CD->row[i][CD->row[i][0]-2]); */
          CD->row[i]  = reduce_lower_by_upper_rows_c(CD->row[i], AB);
        }
        for (nelts_t k=0; k<AB->nr; ++k) {
          free(AB->row[k]);
        }
        free(AB->row);
        free(AB);
        AB  = NULL;
      }
#if newred
      printf("--CD BEFORE--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
#if newred
        printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
      if (CD->rk > 1)
        reduce_lower_rows_c(CD, CD->ncl, nthreads);
#if newred
      printf("rank of CD %u\n", CD->rk);
      printf("--CD AFTER--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] ",ii);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
            printf("%u at %u | ",CD->row[ii][jj],CD->row[ii][jj+1]);
          }
          printf("\n");
        }
      }
#endif
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1 && steps > 1)
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      done  = update_basis_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      free_symbolic_preprocessing_data(&spd);
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
        break;
      }
    }
    if (reduce_gb == 3) {
      /* generate upper matrix, i.e. already known pivots */
      smat_t *AB = NULL, *CD = NULL;
      if (spd->selu->load > 0)
        AB = generate_sparse_matrix(basis, spd->selu,
            spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
            nthreads);
      /* generate lower matrix, i.e. unkown pivots */
      CD = generate_sparse_matrix(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      meta_data->mat_rows = spd->selu->load + spd->sell->load;
      meta_data->mat_cols = CD->ncl + CD->ncr;
      if (verbose > 0)
        t_generating_gbla_matrix  +=  walltime(t_load_start);
      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      if (verbose > 1) {
        printf("matrix rows %6u \n", meta_data->mat_rows);
        printf("matrix cols %6u \n", meta_data->mat_cols);
      }
      if (verbose == 1) {
        printf("step %3d : %5u spairs of degree %3u  --->  %7u x %7u matrix\n",
            steps, meta_data->sel_pairs, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
      }
#if newred
      printf("--CD BEGINNING--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] (%p)",ii, CD->row[ii]);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=0; jj<CD->row[ii]->sz; ++jj) {
            printf("%u at %u | ",CD->row[ii]->val[jj],CD->row[ii]->pos[jj]);
          }
          printf("\n");
        }
      }
#endif
      if (AB != NULL) {
#pragma omp parallel for num_threads(nthreads)
        for (nelts_t i=0; i<CD->nr; ++i) {
          CD->row[i]  = reduce_lower_by_upper_rows(CD->row[i], AB);
#if newred
          printf("row[%u] after = %p\n", i, CD->row[i]);
#endif
        }
        for (nelts_t k=0; k<AB->nr; ++k) {
          free(AB->row[k]->pos);
          free(AB->row[k]->val);
        }
        free(AB->row);
        free(AB);
        AB  = NULL;
      }
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
#if newred
        printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
        if (CD->row[i] != NULL) {
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
#if newred
      printf("--CD BEFORE--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] (%p | %p)",ii, CD->row[ii], CD->row[ii]->pos);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=0; jj<CD->row[ii]->sz; ++jj) {
            printf("%u at %u | ",CD->row[ii]->val[jj],CD->row[ii]->pos[jj]);
          }
          printf("\n");
        }
      }
#endif
      //if (CD->nr > 1)
        reduce_lower_rows(CD, CD->ncl, nthreads);
#if newred
      printf("--CD AFTER--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] (%p)",ii, CD->row[ii]);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=0; jj<CD->row[ii]->sz; ++jj) {
            printf("%u at %u | ",CD->row[ii]->val[jj],CD->row[ii]->pos[jj]);
          }
          printf("\n");
        }
      }
#endif
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);

      done  = update_basis_new_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      for (nelts_t k=0; k<CD->rk; ++k) {
        free(CD->row[k]->pos);
        free(CD->row[k]->val);
      }
      free(CD->row);
      free(CD);
      CD  = NULL;
      free_symbolic_preprocessing_data(&spd);
      clear_hash_table_idx(ht);
      if (verbose > 0)
        t_update_pairs  +=  walltime(t_load_start);
      if (done) {
        basis->has_unit = 1;
        break;
      }
    }
    if (reduce_gb == 5) {
      /* generate upper matrix, i.e. already known pivots */
      smat_t *CD = NULL;
      /* generate lower matrix, i.e. unkown pivots */
      CD = generate_sparse_matrix(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      meta_data->mat_rows = spd->selu->load + spd->sell->load;
      meta_data->mat_cols = CD->ncl + CD->ncr;
      if (verbose > 0)
        t_generating_gbla_matrix  +=  walltime(t_load_start);
      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      if (verbose > 1) {
        printf("matrix rows %6u \n", meta_data->mat_rows);
        printf("matrix cols %6u \n", meta_data->mat_cols);
      }
      if (verbose == 1) {
        printf("step %3d : %5u spairs of degree %3u  --->  %7u x %7u matrix\n",
            steps, meta_data->sel_pairs, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
      }
#if newred
      printf("--CD BEGINNING--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] (%p)",ii, CD->row[ii]);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=0; jj<CD->row[ii]->sz; ++jj) {
            printf("%u at %u | ",CD->row[ii]->val[jj],CD->row[ii]->pos[jj]);
          }
          printf("\n");
        }
      }
#endif
      if (spd->selu->load > 0) {
        for (nelts_t l=0; l<spd->selu->load; ++l) {
          sr_t *reducer = poly_to_sparse_matrix_row(spd->selu->mpp+l, CD->ncl+CD->ncr, basis);
#pragma omp parallel for num_threads(nthreads)
          for (nelts_t i=0; i<CD->nr; ++i) {
            if (CD->row[i] != NULL) {
              CD->row[i]  = reduce_lower_by_one_upper_row(CD->row[i], reducer, CD->mod, CD->ncl+CD->ncr);
            }
          }
          free(reducer->val);
          free(reducer->pos);
          free(reducer);
        }
        /* free(row); */
      }
      nelts_t ctr = 0;
      for (nelts_t i=0; i<CD->nr; ++i) {
#if newred
        printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
        if (CD->row[i] != NULL) {
          CD->row[i] = normalize_row(CD->row[i], CD->mod);
          CD->row[ctr] = CD->row[i];
          ctr++;
        }
      }
      CD->nr  = ctr;
      CD->rk  = ctr;
#if newred
      printf("--CD BEFORE--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] (%p | %p)",ii, CD->row[ii], CD->row[ii]->pos);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=0; jj<CD->row[ii]->sz; ++jj) {
            printf("%u at %u | ",CD->row[ii]->val[jj],CD->row[ii]->pos[jj]);
          }
          printf("\n");
        }
      }
#endif
      /* if (CD->nr > 1) */
      reduce_lower_rows(CD, CD->ncl, nthreads);
#if newred
      printf("--CD AFTER--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] (%p)",ii, CD->row[ii]);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=0; jj<CD->row[ii]->sz; ++jj) {
            printf("%u at %u | ",CD->row[ii]->val[jj],CD->row[ii]->pos[jj]);
          }
          printf("\n");
        }
      }
#endif
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);

      done  = update_basis_new_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      for (nelts_t k=0; k<CD->rk; ++k) {
        free(CD->row[k]->pos);
        free(CD->row[k]->val);
      }
      free(CD->row);
      free(CD);
      CD  = NULL;
      free_symbolic_preprocessing_data(&spd);
      clear_hash_table_idx(ht);
      if (verbose > 0)
        t_update_pairs  +=  walltime(t_load_start);
      if (done) {
        basis->has_unit = 1;
        break;
      }
    }
    if (reduce_gb == 4) {
      /* genearte upper matrix, i.e. already known pivots */
      smat_t *AB = NULL, *CD = NULL;
      /* genearte lower matrix, i.e. unkown pivots */
      CD = generate_sparse_matrix(basis, spd->sell,
          spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
          nthreads);
      meta_data->mat_rows = spd->selu->load + spd->sell->load;
      meta_data->mat_cols = CD->ncl + CD->ncr;
      if (verbose > 0)
        t_generating_gbla_matrix  +=  walltime(t_load_start);
      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);
      if (verbose > 1) {
        printf("matrix rows %6u \n", meta_data->mat_rows);
        printf("matrix cols %6u \n", meta_data->mat_cols);
      }
      if (verbose == 1) {
        printf("step %3d : %5u spairs of degree %3u  --->  %7u x %7u matrix\n",
            steps, meta_data->sel_pairs, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
      }
#if newred
      printf("--CD BEGINNING-1--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] (%p)",ii, CD->row[ii]);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=0; jj<CD->row[ii]->sz; ++jj) {
            printf("%u at %u | ",CD->row[ii]->val[jj],CD->row[ii]->pos[jj]);
          }
          printf("\n");
        }
      }
#endif
      reduce_lower_rows(CD, 0, nthreads);
#if newred
      printf("--CD BEGINNING-2--\n");
      for (int ii=0; ii<CD->nr; ++ii) {
        printf("row[%u] (%p)",ii, CD->row[ii]);
        if (CD->row[ii] == NULL)
          printf("NULL\n");
        else {
          for (int jj=0; jj<CD->row[ii]->sz; ++jj) {
            printf("%u at %u | ",CD->row[ii]->val[jj],CD->row[ii]->pos[jj]);
          }
          printf("\n");
        }
      }
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
        }
        if (AB != NULL) {
#pragma omp parallel for num_threads(nthreads)
          for (nelts_t i=0; i<CD->nr; ++i) {
            CD->row[i]  = reduce_lower_by_upper_rows(CD->row[i], AB);
#if newred
            printf("row[%u] after = %p\n", i, CD->row[i]);
#endif
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
#if newred
          printf("test CD->row[%u] = %p\n", i, CD->row[i]);
#endif
          if (CD->row[i] != NULL) {
            CD->row[ctr] = CD->row[i];
            ctr++;
          }
        }
        CD->nr  = ctr;
        CD->rk  = ctr;
        reduce_lower_rows(CD, CD->ncl, nthreads);
      }
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);

      if (verbose > 0)
        gettimeofday(&t_load_start, NULL);

      done  = update_basis_new_new(basis, ps, spd, CD, ht);
      if (verbose > 0) {
        n_zero_reductions +=  (CD->nr - CD->rk);
      }
      for (nelts_t k=0; k<CD->rk; ++k) {
        free(CD->row[k]->pos);
        free(CD->row[k]->val);
      }
      free(CD->row);
      free(CD);
      CD  = NULL;
      free_symbolic_preprocessing_data(&spd);
      clear_hash_table_idx(ht);
      if (verbose > 0)
        t_update_pairs  +=  walltime(t_load_start);
      if (done) {
        basis->has_unit = 1;
        break;
      }
    }
    if (reduce_gb == 0) {
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

        /* if (mat->B->blocks != NULL) {
         *   const uint32_t clB  = (uint32_t) ceil((float)mat->B->ncols / __GBLA_SIMD_BLOCK_SIZE);
         *   // row loops
         *   const uint32_t rlB  = (uint32_t) ceil((float)mat->B->nrows / __GBLA_SIMD_BLOCK_SIZE);
         *   uint64_t neltsnz  = 0;
         *   for (int ii=0; ii<rlB; ++ii) {
         *     for (int jj=0; jj<clB; ++jj) {
         *       if (mat->B->blocks[ii][jj].val != NULL) {
         *         neltsnz = 0;
         *         for (int kk=0; kk<__GBLA_SIMD_BLOCK_SIZE; ++kk) {
         *           for (int ll=0; ll<__GBLA_SIMD_BLOCK_SIZE; ++ll) {
         *             if (mat->B->blocks[ii][jj].val[kk*__GBLA_SIMD_BLOCK_SIZE+ll] != 0)
         *               neltsnz++;
         *           }
         *         }
         *         printf("B %6lu elements of block [%3u, %3u] are nonzero (%6.3f %%)", neltsnz, ii, jj,(float)neltsnz/(float)65536);
         *         if ((float)neltsnz/(float)65536 > 0.4)
         *           printf(" <---\n");
         *         else
         *           printf("\n");
         *       }
         *     }
         *   }
         * } */
        meta_data->mat_rows = mat->A->nrows+mat->C->nrows;
        meta_data->mat_cols = mat->A->ncols+mat->B->ncols;
        if (verbose > 1) {
          printf("matrix rows %6u + %6u = %6u\n", mat->A->nrows, mat->C->nrows, mat->A->nrows + mat->C->nrows);
          printf("matrix cols %6u + %6u = %6u\n", mat->A->ncols, mat->B->ncols, mat->A->ncols + mat->B->ncols);
        }
      }
      if (verbose == 1) {
        printf("%3d %5u/%5u spairs of deg %3u %7u x %7u matrix ",
            steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
        fflush(stdout);
      }
      /* mat_t *mat  = generate_gbla_matrix(basis, sf, spd, nthreads); */
      if (verbose > 0)
        t_generating_gbla_matrix  +=  walltime(t_load_start);
      /* generate pbm files of gbla matrix */
      if (generate_pbm) {
        int pos = 0;
        pos = snprintf(pbm_fn+pos, 200, "%s/",pbm_dir);
        pos = snprintf(pbm_fn+pos, 200, "%s-mat%u.pbm", fn, steps);
        printf("%s\n", pbm_fn);
        write_matrix_to_pbm(mat, pbm_fn);
      }
      /* reduce matrix using gbla */
      if (verbose == 2) {
        gettimeofday(&t_load_start, NULL);
        printf("[%2u] ", steps);
        printf("%-33s","GBLA matrix reduction ...");
        fflush(stdout);
      }
      ri_t rankDR  = 0;
      if (keep_A == 1)
        rankDR  = reduce_gbla_matrix_keep_A(mat, verbose, nthreads);
      else
        rankDR  = reduce_gbla_matrix(mat, verbose, nthreads);
      if (verbose == 2) {
        printf("%9.3f sec %7d %7d %7d\n",
            walltime(t_load_start) / (1000000), rankDR, (int)mat->DR->nrows - (int)rankDR, mat->DR->nrows);
      }
      if (verbose > 0)
        t_linear_algebra  +=  walltime(t_load_start);
      if (verbose == 1)
        printf("%9.3f sec\n", walltime(t_load_start) / (1000000));

      /* generate pbm files of gbla matrix */
      if (generate_pbm) {
        int pos = 0;
        pos = snprintf(pbm_fn+pos, 200, "%s/",pbm_dir);
        pos = snprintf(pbm_fn+pos, 200, "%s-mat%u-red.pbm", fn, steps);
        printf("%s\n", pbm_fn);
        write_reduced_matrix_to_pbm(mat, pbm_fn);
      }

      /* add new elements to basis and update pair set */
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
        /* See note on gb_t in src/types.h why we decrement basis->load here. */
        printf("basis->load                       %9u\n",basis->load-1); 
        printf("---------------------------------------------------------------------------\n");
        printf("criteria applications (last step) %9u\n",meta_data->ncrit_last);
        printf("criteria applications (total)     %9u\n",meta_data->ncrit_total);
      }

      /*
         if (ps->load == 13936 || ps->load == 13950) {
         for (int ii=0; ii<ps->load; ++ii) {
         printf("pair[%u]   gen1 %6u | gen2 %6u || deg %4u || lcm ", ii, ps->pairs[ii]->gen1, ps->pairs[ii]->gen2, ps->pairs[ii]->deg);
         exp_t expa[ht->nev * ht->vl] __attribute__ ((aligned (16)));
         exp_t tmp[ht->vl] __attribute__ ((aligned (16)));
         for (int jj=0; jj<ht->nev; ++jj) {
         _mm_store_si128((exp_v *)tmp, ht->ev[ps->pairs[ii]->lcm][jj]);
         memcpy(expa+(jj*ht->vl), tmp, ht->vl*sizeof(exp_t));
         }
         for (int kk=0; kk<ht->nv; ++kk) {
         printf("%u ", expa[kk]);
         }
         printf("\n");
         }
         }
         */

      // if we are done then we have found the constant 1 as element in the basis
      if (done) {
        basis->has_unit = 1;
        break;
      }
    }
  }

  /* final basis for possible output data */
  poly_t *fb  = NULL;
  /* generate final basis for output data */
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
#if COUNT_DIV_HITS
    printf("Div Hits                          %9lu/%9lu (%4.2f)\n", meta_data->non_div_found,
        meta_data->non_div, (double)(meta_data->non_div_found)/(double)(meta_data->non_div));
#endif
    printf("Size of basis                     %9u\n", basis->fl);
    printf("criteria applications (total)     %9u\n", meta_data->ncrit_total);
    printf("Number of zero reductions         %9lu\n", n_zero_reductions);
    printf("---------------------------------------------------------------------------\n");
    if (verbose > 2)
      print_mem_usage();
  }

  /* printing of output */
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

  /* free allocated memory */
  free(meta_data);
  free_pair_set(&ps);
  if (basis->sl > 0)
    free_simplifier_list(&sf);
  /* if we have found a unit we have allocated memory for the unit */
  if (basis->has_unit == 1) {
    free(fb[0].cf);
    free(fb[0].eh);
  }
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
    /* we add the polys to sf, we know that there is one coefficient at col pos i
     * for row i.
     * note: we only add these simplifiers if they are not too dense, i.e. at
     * most twice the size of the original polynomial in the basis. this is a
     * heuristic and might be bad in some examples, but in most of our tests it
     * is the fastest choice. */
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
    /* update basis and pair set, mark redundant elements in basis */
#pragma omp single nowait
    {
      /* add simplifier, i.e. polynomials corresponding to the rows in AB,
       * for further computation */
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

void reduce_gb_23(gb_t *basis, const spd_t *spd,
    const double density, const ps_t *ps, const nelts_t block_size,
    const int verbose, const int nthreads)
{
  mat_gb_block_t *AB, *CD   = NULL;
  mat_gb_meta_data_t *meta  = NULL;
  if (verbose > 0)
    t_linear_algebra_local  = 0;

  /* generate meta data */
  meta  = generate_matrix_meta_data(block_size, basis->mod, spd);

  /* generate CD part */
  /* printf("gen CD block %u\n", spd->sell->load); */
  CD  = generate_mat_gb_lower(meta, basis, spd, ht, nthreads);

  meta_data->mat_rows = spd->selu->load + spd->sell->load;
  meta_data->mat_cols = meta->nc_AC + meta->nc_BD;
  if (verbose > 0)
    t_generating_gbla_matrix  +=  walltime(t_load_start);
  if (verbose > 0)
    gettimeofday(&t_load_start, NULL);
  if (verbose > 1) {
    printf("matrix rows %6u \n", meta_data->mat_rows);
    printf("matrix cols %6u \n", meta_data->mat_cols);
  }
  if (verbose == 1) {
    printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
        steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
    fflush(stdout);
  }

  if (spd->selu->load > 0) {
    nelts_t i, j;
    for (i=0; i<meta->nrb_AB; ++i) {
#define OLD 0
#if OLD
      AB  = generate_mat_gb_upper_row_block(i, meta, basis, spd, ht);
#else
      /* printf("gen AB block\n"); */
      AB  = generate_mat_gb_row_block(i, meta, basis, spd, ht);
#endif
      if (verbose > 0) {
        t_generating_gbla_matrix  +=  walltime(t_load_start);
        gettimeofday(&t_load_start, NULL);
      }

      /* possibly the first block of A is already the unit matrix, then it is
       * just an empty block and we do not need to update AB at all */
#if OLD
      if (AB[i].len != NULL)
        update_upper_row_block(AB, i, meta, nthreads);

      update_lower_by_upper_row_block(CD, AB, i, meta, nthreads);
#else
      update_lower_by_upper_row_block_not_prereduced(CD, AB, i, meta, nthreads);
#endif
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
#if newred
  printf("--D BEGINNING--\n");
  for (int ii=0; ii<D->nr; ++ii) {
    printf("row[%u] ",ii);
    if (D->row[ii] == NULL)
      printf("NULL\n");
    else {
      for (int jj=1; jj<D->row[ii][0]; jj += 2) {
        printf("%u at %u | ",D->row[ii][jj+1],D->row[ii][jj]);
      }
      printf("\n");
    }
  }
#endif
  nelts_t j, k;
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
  }
  if (verbose == 1 && steps > 0) {
    printf("%4u new %4u zero ", D->rk, spd->sell->load-D->rk);
    printf("%9.3f sec ", t_linear_algebra_local / (1000000));
    printf("%6.3f%% d ", density);
    printf("%6u bs\n", meta->bs);
  }
  free(meta);
  if (verbose > 0)
    gettimeofday(&t_load_start, NULL);
  int done  = update_basis_new(basis, ps, spd, D, ht);
  if (verbose > 0) {
    n_zero_reductions +=  (spd->sell->load - D->rk);
  }
  free_symbolic_preprocessing_data(&spd);
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

void reduce_gb_101(gb_t *basis, const spd_t *spd,
    const double density, const ps_t *ps, const nelts_t block_size,
    const int verbose, const int nthreads)
{
  /* we generate the matrix in the following shape:
   * AB = already known lead terms resp. pivots
   * --
   * CD = new data, new lead terms to be computed
   *
   * genearte upper matrix, i.e. already known pivots */
  nelts_t init_rk_CD  = 0;
  smc_t *AB = NULL, *CD = NULL;
  /* printf("selu->load %u | sell->load %u | nc %u\n", spd->selu->load, spd->sell->load, spd->col->load); */
  /* genearte lower matrix, i.e. unkown pivots */
  CD = generate_sparse_compact_matrix_test(basis, spd->sell,
      spd->sell->load, spd->col->nlm, spd->col->load-spd->col->nlm,
      nthreads);
#if 0
  printf("--CD 1--\n");
  for (int ii=0; ii<CD->nr; ++ii) {
    printf("row[%u] ",ii);
    if (CD->row[ii] == NULL)
      printf("NULL\n");
    else {
      for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
        printf("%u at %u | ",CD->row[ii][jj+1],CD->row[ii][jj]);
      }
      printf("\n");
    }
  }
#endif
  init_rk_CD  = CD->nr;
  /* check appearing columns in CD */
#define COL_CHECK 0
#if COL_CHECK
  uint32_t *columns  = (uint32_t *)calloc(spd->col->nlm, sizeof(uint32_t));
  for (int ii=0; ii<CD->nr; ++ii) {
    for (int jj=1; jj<CD->row[ii][0]; jj=jj+2) {
      if (CD->row[ii][jj] < spd->col->nlm)
        columns[CD->row[ii][jj]] = 1;
    }
  }
  uint32_t colctr  = 0;
  for (uint32_t ii=0; ii<spd->col->nlm; ++ii) {
    if (columns[ii] != 0)
      colctr++;
  }
  /* printf("\n\nCD %6u / %6u cols | %4.2f %\n", colctr, spd->col->nlm, */
  (float)colctr / (float)spd->col->nlm);
#endif
  if (spd->selu->load > 0) {
    AB = generate_sparse_compact_matrix_offset_test(basis, spd->selu,
        spd->selu->load, spd->col->nlm, spd->col->load-spd->col->nlm,
        nthreads);
#if 0
    printf("--AB 1--\n");
    for (int ii=0; ii<AB->nr; ++ii) {
      printf("row[%u] ",ii);
      if (AB->row[ii] == NULL)
        printf("NULL\n");
      else {
        for (int jj=2; jj<AB->row[ii][0]; jj += 2) {
          printf("%u at %u | ",AB->row[ii][jj+1],AB->row[ii][jj]);
        }
        printf("\n");
      }
    }
#endif
#if COL_CHECK
    uint32_t *appears  = (uint32_t *)calloc(spd->col->nlm, sizeof(uint32_t));
    for (int ii=0; ii<AB->nr; ++ii) {
      // starting at 2 due to offsets!
      for (int jj=2; jj<AB->row[ii][0]; jj=jj+2) {
        if (AB->row[ii][jj] < spd->col->nlm && columns[AB->row[ii][jj]] == 0)
          appears[AB->row[ii][jj]]++;
      }
    }
    uint32_t ac = 0;
    for (int ii=0; ii<AB->nr; ++ii) {
      // starting at 2 due to offsets!
      for (int jj=2; jj<AB->row[ii][0]; jj=jj+2) {
        if (AB->row[ii][jj] < spd->col->nlm && columns[AB->row[ii][jj]] == 0) {
          ac++;
          break;
        }
      }
    }
    /* printf("AB %6u / %6u rows | %4.2f %\n\n",ac, AB->nr, (float)ac/(float)AB->nr); */
    free(appears);
#endif
    meta_data->mat_rows = spd->selu->load + spd->sell->load;
    meta_data->mat_cols = CD->ncl + CD->ncr;
    if (verbose > 0)
      t_generating_gbla_matrix  +=  walltime(t_load_start);
    if (verbose > 0)
      gettimeofday(&t_load_start, NULL);
    if (verbose > 1) {
      printf("matrix rows %6u \n", meta_data->mat_rows);
      printf("matrix cols %6u \n", meta_data->mat_cols);
    }
    if (verbose == 1) {
      printf("%3d %5u/%5u pairs deg %3u %7u x %7u mat ",
          steps-1, meta_data->sel_pairs, meta_data->sel_pairs+ps->load, meta_data->curr_deg, meta_data->mat_rows, meta_data->mat_cols);
      fflush(stdout);
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
  /* if (verbose == 1 && steps > 0) {
   *   printf("%9.3f sec ", walltime(t_load_start) / (1000000));
   * } */
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
#if 0
  printf("--CD 2--\n");
  for (int ii=0; ii<CD->nr; ++ii) {
    printf("row[%u] ",ii);
    if (CD->row[ii] == NULL)
      printf("NULL\n");
    else {
      for (int jj=1; jj<CD->row[ii][0]; jj += 2) {
        printf("%u at %u | ",CD->row[ii][jj+1],CD->row[ii][jj]);
      }
      printf("\n");
    }
  }
#endif
  CD->nr  = ctr;
  CD->rk  = ctr;
  /* printf("rank of CD %u | %u\n", CD->rk, CD->nr); */
  if (CD->rk > 1)
    reduce_lower_rows_c(CD, CD->ncl, nthreads);
  if (verbose > 0)
    t_linear_algebra  +=  walltime(t_load_start);
  if (verbose == 1 && steps > 1)
    printf("%4u new %4u zero ", CD->rk, init_rk_CD-CD->rk);
  printf("%9.3f sec ", walltime(t_load_start) / (1000000));
  printf("%6.3f%% d\n", density);

  if (verbose > 0)
    gettimeofday(&t_load_start, NULL);
  int done  = update_basis_new(basis, ps, spd, CD, ht);
  if (verbose > 0) {
    n_zero_reductions +=  (init_rk_CD - CD->rk);
  }
  free_symbolic_preprocessing_data(&spd);
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
