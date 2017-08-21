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

#include "io.h"

/*  ========== TIMINGS and MEMORY PRINTING ========== */
double walltime(struct timeval t_start)
{
	struct timeval t_end;
	gettimeofday(&t_end, NULL);
	return (double)((t_end.tv_sec - t_start.tv_sec) * 1000000 + t_end.tv_usec - t_start.tv_usec);
}

void print_mem_usage()
{
  const char *units[4] = {"KB", "MB", "GB", "TB"};
  int ctr = 0;
  double  vms   = 0.0; /*  virtual memory size */
  double  rss   = 0.0; /*  resident set size */
  /*  possibly x86-64 is configured to use 2MB pages */
  long    page_size_kb  = sysconf(_SC_PAGE_SIZE) / 1024;

  unsigned long _vms = 0;
  long          _rss = 0;
  /*  get memory usage from 'file' stat which is not perfect, but gives most */
  /*  reliable information. */
  /*  Note: This corresponds to Martani's memory usage printing in his LELA */
  /*  implementation, thus it is used for comparison reasons. It might be changed */
  /*  in later versions of the gb. */
  const char *fn  = "/proc/self/stat";
  FILE *fh        = fopen(fn,"r");

  /*  dummy vars for leading entries in /proc/self/stat we are not interested in */
  char pid[1024] = "\0", comm[1024] ="\0", state[1024] ="\0", ppid[1024] ="\0", pgrp[1024] ="\0", session[1024] ="\0", tty_nr[1024] ="\0";
  char tpgid[1024] ="\0", flags[1024] ="\0", minflt[1024] ="\0", cminflt[1024] ="\0", majflt[1024] ="\0", cmajflt[1024] ="\0";
  char utime[1024] ="\0", stime[1024] ="\0", cutime[1024] ="\0", cstime[1024] ="\0", priority[1024] ="\0", nice[1024] ="\0";
  char nthrds[1024] ="\0", itrealvalue[1024] ="\0", starttime[1024] ="\0";

  /*  dummy reading of useless information */
  fscanf(fh, "%1023s", &pid[0]);
  fscanf(fh, "%1023s", &comm[0]);
  fscanf(fh, "%1023s", &state[0]);
  fscanf(fh, "%1023s", &ppid[0]);
  fscanf(fh, "%1023s", &pgrp[0]);
  fscanf(fh, "%1023s", &session[0]);
  fscanf(fh, "%1023s", &tty_nr[0]);
  fscanf(fh, "%1023s", &tpgid[0]);
  fscanf(fh, "%1023s", &flags[0]);
  fscanf(fh, "%1023s", &minflt[0]);
  fscanf(fh, "%1023s", &cminflt[0]);
  fscanf(fh, "%1023s", &majflt[0]);
  fscanf(fh, "%1023s", &cmajflt[0]);
  fscanf(fh, "%1023s", &utime[0]);
  fscanf(fh, "%1023s", &stime[0]);
  fscanf(fh, "%1023s", &cutime[0]);
  fscanf(fh, "%1023s", &cstime[0]);
  fscanf(fh, "%1023s", &priority[0]);
  fscanf(fh, "%1023s", &nice[0]);
  fscanf(fh, "%1023s", &nthrds[0]);
  fscanf(fh, "%1023s", &itrealvalue[0]);
  fscanf(fh, "%1023s", &starttime[0]);

  /*  get real memory information */
  fscanf(fh, "%lu", &_vms);
  fscanf(fh, "%ld", &_rss);

  /*  close file */
  fclose(fh);

  /*  TODO: How to read /proc/self/stat ??? */

  vms = (double)_vms / 1024.0;
  rss = (double)(_rss * page_size_kb);

  /*  MB ? */
  if (vms > 1024) {
    vms   = vms/1024.0;
    rss   = rss/1024.0;
    ctr++;
  }
  /*  GB ? */
  if (vms > 1024) {
    vms   = vms/1024.0;
    rss   = rss/1024.0;
    ctr++;
  }
  /*  TB ? Just joking! */
  if (vms > 1024) {
    vms   = vms/1024.0;
    rss   = rss/1024.0;
    ctr++;
  }
  printf("MMRY\tRSS - %.3f %s | VMS - %.3f %s\n", rss, units[ctr], vms, units[ctr]);
}

inline char *get_variable_name(const char *line, char **prev_pos)
{
  const char comma_splicer  = ',';

  char *tmp_var   = (char *)malloc(50 * sizeof(char));
  char *curr_pos  = strchr(*prev_pos, comma_splicer);
  if (curr_pos != NULL) {
    size_t var_diff   = (size_t)(curr_pos - *prev_pos);
    memcpy(tmp_var, *prev_pos, var_diff);
    tmp_var[var_diff] = '\0';
    *prev_pos = curr_pos+1;
  } else { /** we are at the last variable */
    int prev_idx    = (int)(*prev_pos - line);
    int curr_idx    = (int)(strlen(line)+1);
    size_t var_diff = (size_t)(curr_idx - prev_idx);
    memcpy(tmp_var, *prev_pos, var_diff);
    tmp_var[var_diff] = '\0';
  }

  /** trim variable, remove blank spaces */
  char *tmp_var_begin, *tmp_var_end;
  tmp_var_begin = tmp_var_end =  tmp_var;
  while (isspace(*tmp_var_begin))
    tmp_var_begin++;
  if (*tmp_var_begin != 0) {
    while (*tmp_var_end)
      tmp_var_end++;
    do {
      tmp_var_end--;
    } while (isspace(*tmp_var_end));
  }
  size_t final_length  = (size_t)(tmp_var_end - tmp_var_begin + 1);
  char *final_var   = (char *)malloc((final_length+1) * sizeof(char));
  memcpy(final_var, tmp_var_begin, final_length);
  final_var[final_length] = '\0';
  
  free(tmp_var);

  return final_var;
}

inline nelts_t get_number_of_terms(const char *line)
{
  const char add_splicer        = '+';
  const char minus_splicer      = '-';
  const char whitespace_splicer = ' ';
  char *tmp;
  nelts_t nterms  = 1;
  /** remove useless whitespaces at the beginning */
  int i = 0;
  while (strncmp(&whitespace_splicer, line+i, 1) == 0)
    i++;
  /** check if first non-whitespace char is "-", set term counter -1 in this case */
  if (strncmp(&minus_splicer, line+i, 1) == 0)
    nterms--;
  /** now count terms */
  tmp = strchr(line, add_splicer);
  while (tmp != NULL) {
    nterms++;
    tmp = strchr(tmp+1, add_splicer);
  }
  tmp = strchr(line, minus_splicer);
  while (tmp != NULL) {
    nterms++;
    tmp = strchr(tmp+1, minus_splicer);
  }

  return nterms;
}
inline void store_exponent(const char *term, const gb_t *basis, ht_t *ht)
{
  nvars_t k;
  /** first we have to fill the buffers with zeroes */
  memset(ht->exp[ht->load], 0, ht->nv * sizeof(exp_t));
  const char mult_splicer = '*';
  const char exp_splicer  = '^';
  exp_t exp = 0;
  deg_t deg = 0;

  for (k=0; k<basis->rnv; ++k) {
    exp = 0;
    char *var = strstr(term, basis->vnames[k]);
    if (var != NULL) {
#if IO_DEBUG
      printf("var1 %s\n", var);
#endif
      var   = strtok(var, "\n");
      var   = strtok(var, ",");
#if IO_DEBUG
      printf("var2 %s\n", var);
#endif
      /** if the next variable follows directly => exp = 1 */
      if (strncmp(&mult_splicer, var+strlen(basis->vnames[k]), 1) == 0) {
        exp = 1;
      } else {
        /** if there follows an exp symbol "^" */
        if (strncmp(&exp_splicer, var+strlen(basis->vnames[k]), 1) == 0) {
          char exp_str[1000];
          char *mult_pos;
          mult_pos  = strchr(var, mult_splicer);
          if (mult_pos != NULL) {
            size_t exp_len = (size_t)(mult_pos - (var+strlen(basis->vnames[k])) - 1);
            memcpy(exp_str, var+strlen(basis->vnames[k])+1, exp_len);
            exp_str[exp_len] = '\0';
            exp = (exp_s)strtol(exp_str, NULL, 10);
          } else { /** no further variables in this term */
            size_t exp_len = (size_t)((var+strlen(var)) + 1 - (var+strlen(basis->vnames[k])) - 1);
            memcpy(exp_str, var+strlen(basis->vnames[k])+1, exp_len);
            exp = (exp_s)strtol(exp_str, NULL, 10);
            exp_str[exp_len] = '\0';
          }
        } else { /** we are at the last variable with exp = 1 */
          if (strcmp(basis->vnames[k], var) == 0)
            exp = 1;
          else
            continue;
        }
      }
    }
    /** if we use graded reverse lexicographical order (basis->ord = 0) or
      * lexicographic order (basis->ord = 1)  we store
      * the exponents in reverse order so that we can use memcmp to sort the terms
      * efficiently later on */
    if (basis->ord == 0)
      deg +=  ht->exp[ht->load][ht->nv-1-k] = exp;
    else
      deg +=  ht->exp[ht->load][k] = exp;
  }
  ht->deg[ht->load] = deg;
}

inline void get_term(const char *line, char **prev_pos,
    char **term)
{
  /** note that maximal term length we handle */
  const char add_splicer    = '+';
  const char minus_splicer  = '-';

  char *start_pos;
  char *curr_pos_add    = strchr(*prev_pos, add_splicer);
  char *curr_pos_minus  = strchr(*prev_pos, minus_splicer);
  if (curr_pos_minus == line)
    curr_pos_minus  = strchr(*prev_pos+1, minus_splicer);

  if (*prev_pos != line)
    start_pos = *prev_pos - 1;
  else
    start_pos = *prev_pos;

  if (curr_pos_add != NULL && curr_pos_minus != NULL) {
    size_t term_diff_add   = (size_t)(curr_pos_add - start_pos);
    size_t term_diff_minus = (size_t)(curr_pos_minus - start_pos);
    /** if "-" is the first char in the line, we have to adjust the
      * if minus is nearer */
    if (term_diff_add > term_diff_minus) {
      memcpy(*term, start_pos, term_diff_minus);
      (*term)[term_diff_minus]  = '\0';
      *prev_pos                 = curr_pos_minus+1;
      return;
    /** if plus is nearer */
    } else {
      memcpy(*term, start_pos, term_diff_add);
      (*term)[term_diff_add]  = '\0';
      *prev_pos               = curr_pos_add+1;
      return;
    }
  } else {
    if (curr_pos_add != NULL) {
      size_t term_diff_add   = (size_t)(curr_pos_add - start_pos);
      memcpy(*term, start_pos, term_diff_add);
      (*term)[term_diff_add]  = '\0';
      *prev_pos               = curr_pos_add+1;
      return;
    }
    if (curr_pos_minus != NULL) {
      size_t term_diff_minus = (size_t)(curr_pos_minus - start_pos);
      memcpy(*term, start_pos, term_diff_minus);
      (*term)[term_diff_minus]  = '\0';
      *prev_pos                 = curr_pos_minus+1;
      return;
    }
    if (curr_pos_add == NULL && curr_pos_minus == NULL) {
      size_t prev_idx  = (size_t)(start_pos - line);
      size_t term_diff = strlen(line) + 1 - prev_idx;
      memcpy(*term, start_pos, term_diff);
      (*term)[term_diff]  = '\0';
      return;
    }
  }
}

nvars_t get_nvars(const char *fn)
{
  FILE *fh  = fopen(fn,"r");
  /** load lines and store data */
  const int max_line_size  = 1000;
  char *line  = (char *)malloc((nelts_t)max_line_size * sizeof(char));

  /** get first line (variables) */
  const char comma_splicer  = ',';

  /** get number of variables */
  nvars_t nvars = 1; /** number of variables is number of commata + 1 in first line */
  if (fgets(line, max_line_size, fh) != NULL) {
    char *tmp = strchr(line, comma_splicer);
    while (tmp != NULL) {
      /** if there is a comma at the end of the line, i.e. strlen(line)-2 (since
        * we have "\0" at the end of the string line, then we do not get another
        * variable */
      if ((uint32_t)(tmp-line) < strlen(line)-2) {
      nvars++;
      tmp = strchr(tmp+1, comma_splicer);
      } else {
        break;
      }
    }
  } else {
    printf("Bad file format.\n");
    nvars = 0;
  }
  free(line);
  fclose(fh);

  return nvars;
}

void check_for_same_exponents(gb_t *basis)
{
  cf_t *tmp_cf;
  hash_t *tmp_eh;
  cf_t *sort_cf   = (cf_t *)malloc(sizeof(cf_t));
  hash_t *sort_eh = (hash_t *)malloc(sizeof(hash_t));
  nelts_t i, j;
  nelts_t ctr;
  bf_t coeff;

  /** first we sort the terms of each input polynomial w.r.t. the monomial order */
  for (i=1; i<basis->load; ++i) {
    sort_eh  = realloc(sort_eh, basis->nt[i] * sizeof(hash_t));
    sort_cf  = realloc(sort_cf, basis->nt[i] * sizeof(cf_t));
    ctr = 0;
    /** lead term is clear */
    sort_eh[0]  = basis->eh[i][0];
    sort_cf[0]  = basis->cf[i][0];
    ctr++;
    /** check all other exponents with previous exponents, probably add
      * coefficients */
    for (j=1; j<basis->nt[i]; ++j) {
      if (basis->eh[i][j] == sort_eh[ctr-1]) {
        coeff = (bf_t)sort_cf[ctr-1] + basis->cf[i][j];
        coeff %=  basis->mod;
        sort_cf[ctr-1]  = (cf_t)coeff;
      } else {
        sort_eh[ctr]  = basis->eh[i][j];
        sort_cf[ctr]  = basis->cf[i][j];
        ctr++;
      }
    }
    /** reallocate memory */
    sort_eh  = realloc(sort_eh, ctr * sizeof(hash_t));
    sort_cf  = realloc(sort_cf, ctr * sizeof(cf_t));

    /** swap arrays */
    tmp_cf  = basis->cf[i];
    tmp_eh  = basis->eh[i];

    basis->cf[i]  = sort_cf;
    basis->eh[i]  = sort_eh;
    basis->nt[i]  = ctr;
    
    sort_cf = tmp_cf;
    sort_eh = tmp_eh;
  }
  /** free temporary allocated memory */
  free(sort_cf);
  free(sort_eh);
}

void sort_input_polynomials(gb_t *basis, const ht_t *ht)
{
  cf_t *tmp_cf;
  hash_t *tmp_eh;
  cf_t *sort_cf  = (cf_t *)malloc(sizeof(cf_t));
  hash_t *sort_eh   = (hash_t *)malloc(sizeof(hash_t));
  nelts_t i, j, k;

  /** first we sort the terms of each input polynomial w.r.t. the monomial order */
  for (i=1; i<basis->load; ++i) {
    /** sort exponent hashes */
    sort_eh  = realloc(sort_eh, basis->nt[i] * sizeof(hash_t));
    sort_cf  = realloc(sort_cf, basis->nt[i] * sizeof(cf_t));
    memcpy(sort_eh, basis->eh[i], basis->nt[i] * sizeof(hash_t));
    qsort(sort_eh, basis->nt[i], sizeof(hash_t), ht->sort.compare_monomials);
    /** sort coefficients like the exponent hashes */
    for (j=0; j<basis->nt[i]; ++j) {
      for (k=0; k<basis->nt[i]; ++k) {
        if (basis->eh[i][k] == sort_eh[j]) {
          sort_cf[j]  = basis->cf[i][k];
        }
      }
    }
    /** swap arrays */
    tmp_cf  = basis->cf[i];
    tmp_eh  = basis->eh[i];

    basis->cf[i]  = sort_cf;
    basis->eh[i]  = sort_eh;

    /* normalize polynomials */
    if (basis->cf[i][0] != 1) {
      cf_t iv = basis->cf[i][0];
      inverse_coefficient(&iv, basis->mod);
      basis->cf[i][0]  = 1;
      for (j=1; j<basis->nt[i]; ++j)
        basis->cf[i][j] = MODP(basis->cf[i][j]*iv,basis->mod);
    }
    
    sort_cf = tmp_cf;
    sort_eh = tmp_eh;
  }
  /** sort the polynomials by increasing lead term w.r.t. the monomial order
    * note that basis elements start at position 1, not zero! thus we have
    * basis->load and not basis->load-1 elements stored at the moment. */
  sort_eh = realloc(sort_eh, basis->load * sizeof(hash_t));
  for (i=1; i<basis->load; ++i)
    sort_eh[i]  = basis->eh[i][0];
  qsort(sort_eh+1, basis->load-1, sizeof(hash_t), ht->sort.compare_monomials_inverse);

  /** stores if a position is already set */
  uint8_t *pos_set = (uint8_t *)calloc(basis->load, sizeof(uint8_t));

  for (i=1; i<basis->load; ++i) {
    basis->deg[i] = basis->deg[i];
    basis->nt[i]  = basis->nt[i];
    basis->eh[i]  = basis->eh[i];
    basis->cf[i]  = basis->cf[i];
  }
  /** free temporary allocated memory */
  free(sort_cf);
  free(sort_eh);
  free(pos_set);
}

void homogenize_input_polynomials(gb_t *basis, ht_t *ht) {
  nelts_t i, j;
  /** make sure that hash table is big enough */
  if (2*basis->load > ht->sz)
    enlarge_hash_table(ht);
  /** use extra variable in exponent vector representation in hash table to
    * homogenize the terms correspondingly */
  for (i=1; i<basis->load; ++i) {
    for (j=0; j<basis->nt[i]; ++j) {
      memcpy(ht->exp[ht->load], ht->exp[basis->eh[i][j]], ht->nv * sizeof(exp_t));
      /** add homogenizing variable entry in exponent */
      ht->exp[ht->load][ht->nv-1] = (exp_t)(basis->deg[i] -  ht->deg[basis->eh[i][j]]);
      /** add new exponent hash to table */
      ht->deg[ht->load] = basis->deg[i];
      ht->val[ht->load] = get_hash(ht->exp[ht->load], ht);
      basis->eh[i][j] = check_in_hash_table(ht);
    }
  }
  basis->hom  = 1;
}

gb_t *load_input(const char *fn, const nvars_t nvars, const int order,
    ht_t *ht, const long max_spairs,
    const int vb)
{
  int64_t fl;

  hash_t i, j;

  struct timeval t_load_start;
  if (vb > 1) {
    gettimeofday(&t_load_start, NULL);
  }

  /** open file in binary mode and get file size */
  FILE	*fh  = fopen(fn,"rb");
  if (fh == NULL) {
    if (vb > 0)
      printf("File not found!\n");
    return NULL;
  } else {
    fseek(fh, 0L, SEEK_END);
    fl  = ftell(fh);
    fclose(fh);
  }

  /** get number of lines in file:
    * number of lines - 2 is number of generators in input system */
  int nlines  = 0;
  char buf[10000];
  fh  = fopen(fn,"r");
  while(fgets(buf, sizeof(buf), fh) != NULL) {
    /** check if there are empty lines in the input file */
    if (is_line_empty(buf) == 0)
      nlines++;
  }
  fclose(fh);

  fh  = fopen(fn,"r");
  /** load lines and store data */
  const int max_line_size  = 100000;
  char *line  = (char *)malloc((nelts_t)max_line_size * sizeof(char));

  /** we already know the number of variables
    * basis->rnv  = nvars; */

  char *tmp;
  /** allocate memory for storing variable names */
  char **vnames = (char **)malloc(nvars * sizeof(char *));
  if (fgets(line, max_line_size, fh) != NULL) {
    tmp = line;
  } else {
    printf("Bad file format.\n");
    return NULL;
  }
  for (i=0; i<nvars; ++i) {
    vnames[i]  = get_variable_name(line, &tmp);
  }

  /** get second line (modulus) */
  mod_t mod = 0;
  if (fgets(line, max_line_size, fh) != NULL) {
    int64_t tmp_mod = atol(line);
    if (tmp_mod > 0) {
      mod  = (mod_t)tmp_mod;
    } else {
      printf("Bad file format.\n");
      return NULL;
    }
  } else {
    printf("Bad file format.\n");
    return NULL;
  }

  /** initialize basis with information from above */
  gb_t *basis = initialize_basis(order, nlines, nvars, vnames, mod,
      max_spairs, fl);

  char *prev_pos;
  char *term  = (char *)malloc(200 * sizeof(char));
  nelts_t nterms;
  deg_t max_deg;

  /** NOTE: For easier divisibility checks in symbolic preprocessing we put at
    * the first position of basis, i.e. index 0 a NULL element.
    * Thus, basis->load is always one bigger than the actual number of elements
    * in the basis. */
  basis->p[0]  = NULL;

  /** get all remaining lines, i.e. generators */
  int cf_tmp  = 0; /** temp for coefficient value, possibly coeff is negative. */
  int iv_tmp  = 0; /** temp for inverse value, possibly coeff is negative. */
  cf_t iv  = 0; /** inverse value of lead coeff in order to normalize input */
  for (i=1; i<basis->load; ++i) {
    if (fgets(line, max_line_size, fh) != NULL && is_line_empty(line) != 1) {
      /** get number of terms first */
      nterms          = get_number_of_terms(line);
      basis->p[i]     = (poly_t *)malloc((2*nterms+2) * sizeof(poly_t));
      basis->p[i][0]  = 2 * nterms + 2; 
      basis->p[i][1]  = nterms % 4; 


#if IO_DEBUG
      printf("nterms %d\n",nterms);
#endif
      prev_pos  = line;
      max_deg   = 0;
      /** next: go through line, term by term
        * we do first term differently since we normalize polynomial with lead
        * coefficient */
      get_term(line, &prev_pos, &term);
#if IO_DEBUG
      printf("%u : %s ",i, term);
#endif
      /** get coefficient first */
      if (term != NULL) {
        iv_tmp  = (int)strtol(term, NULL, 10);
        /** if shortcut notation is used coeff 1 is not written down and strtol
          * boils down to 0. so we adjust this value to 1 again */
        if (iv_tmp == 0) {
          if (term[0] == '-') {
            iv_tmp = -1;
          } else {
            iv_tmp = 1;
          }
        }
        while (iv_tmp < 0) {
          iv_tmp  +=  (int)basis->mod;
        }
        iv  = (cf_t)iv_tmp;
        basis->p[i][3] = iv;
        /* normalization is done at the end when sorting the polynomials
         * printf("iv %u\n", iv);
         * inverse_coefficient(&iv, basis->mod);
         * printf("iv %u\n", iv);
         * basis->cf[i][0] = 1; */
      }
      store_exponent(term, basis, ht);
      /** hash exponent and store degree */
      max_deg         = max_deg > ht->deg[ht->load] ? max_deg : ht->deg[ht->load];
      ht->val[ht->load] = get_hash(ht->exp[ht->load], ht);
      basis->p[i][2] = check_in_hash_table(ht);
#if IO_DEBUG
      printf("cf[%lu] = %u | eh[%lu][%u] = %lu --> %lu\n",i,basis->cf[i][0],i,0,basis->eh[i][0], ht->val[basis->eh[i][0]]);
#endif
      for (j = 4; j < basis->p[i][0]; j = j+2) {
        get_term(line, &prev_pos, &term);
#if IO_DEBUG
        printf("%s ",term);
#endif
        /** get coefficient first */
        if (term != NULL) {
          cf_tmp  = (int)strtol(term, NULL, 10);
          if (cf_tmp == 0) {
            if (term[0] == '-') {
              cf_tmp = -1;
            } else {
              cf_tmp = 1;
            }
          }
          while (cf_tmp < 0) {
            cf_tmp  += (int)basis->mod;
          }
          basis->p[i][j+1] = (cf_t)cf_tmp;
          /* normalization is done at the end when sorting the polynomials
           * basis->cf[i][j] = MODP(basis->cf[i][j]*iv,basis->mod); */
        }
        store_exponent(term, basis, ht);
        /** hash exponent and store degree */
        max_deg           = max_deg > ht->deg[ht->load] ? max_deg : ht->deg[ht->load];
        ht->val[ht->load] = get_hash(ht->exp[ht->load], ht);
        basis->p[i][j]    = check_in_hash_table(ht);
#if IO_DEBUG
        printf("cf[%lu] = %u | eh[%lu][%lu] = %lu --> %lu\n",i,basis->cf[i][j],i,j,basis->eh[i][j], ht->val[basis->eh[i][j]]);
#endif
      }
      /** if basis->init_hom is 0 then we have already found an inhomogeneous
        * polynomial and the system of polynomials is not homogeneous */
      if (basis->init_hom == 1) {
        if (ht->deg[basis->p[i][2]] == ht->deg[basis->p[i][basis->p[i][0]-2]])
          basis->init_hom  = 1;
        else
          basis->init_hom  = 0;
      }
      basis->deg[i] = max_deg;
    } else {
      /** the line is empty, thus we have to reset i by -1 and continue the loop */
      i--;
      continue;
    }
#if IO_DEBUG
    printf("\n");
#endif
  }
  basis->hom  = basis->init_hom;
  /** for orderings not degree compatible it is at the moment better to just
    * homogenize the input and dehomogenize at the end. */
  if (basis->ord != 0 && basis->hom == 0)
    homogenize_input_polynomials(basis, ht);
  /** sort polynomial w.r.t. the chosen monomial order */
  sort_input_polynomials(basis, ht);
  check_for_same_exponents(basis);
  recalculate_divmaps(ht);
  ht->rcdm  = 10000;
  free(term);
  free(line);
  fclose(fh);

  return basis;
}

void write_reduced_matrix_to_pbm(mat_t *mat, const char *fn)
{
  ri_t i;
	ri_t m        = mat->nru + mat->nrl;
	ci_t n        = mat->ncl + mat->ncr;
  ri_t min      = n > 512 ? n : 512;
  /** min+2 since we need to end line with '\n\0' */
	char *buffer  = malloc(min+2 * sizeof(char));

	FILE *fh  = fopen(fn, "wb");

	/*  magic PBM header */
#ifdef __LP64__ /*  64bit machine */
	sprintf(buffer, "P1\n# matrix size(%u, %u)\n%u %u\n", m, n, n, m);
#else /*  32bit machine */
	sprintf(buffer, "P1\n# matrix size(%u, %u)\n%u %u\n", m, n, n, m);
#endif

	fwrite(buffer, sizeof(char), strlen(buffer), fh);
  
  /** write top block AB */
  for (i=0; i<mat->nru; ++i) {
    write_upper_part_row_to_buffer(buffer, i, mat);
	  fwrite(buffer, sizeof(char), strlen(buffer), fh);
    fflush(fh);
  }
  /** write bottom block CD */
  for (i=0; i<mat->nrl; ++i) {
    write_lower_part_row_to_buffer(buffer, i, mat);
	  fwrite(buffer, sizeof(char), strlen(buffer), fh);
    fflush(fh);
  }
  fclose(fh);
  free(buffer);
}

void write_matrix_to_pbm(mat_t *mat, const char *fn)
{
  ri_t i;
	ri_t m        = mat->nru + mat->nrl;
	ci_t n        = mat->ncl + mat->ncr;
  ri_t min      = n > 512 ? n : 512;
  /** min+2 since we need to end line with '\n\0' */
	char *buffer  = malloc(min+2 * sizeof(char));
	FILE *fh  = fopen(fn, "wb");

	/*  magic PBM header */
#ifdef __LP64__ /*  64bit machine */
	sprintf(buffer, "P1\n# matrix size(%u, %u)\n%u %u\n", m, n, n, m);
#else /*  32bit machine */
	sprintf(buffer, "P1\n# matrix size(%u, %u)\n%u %u\n", m, n, n, m);
#endif

	fwrite(buffer, sizeof(char), strlen(buffer), fh);

  /** write top block AB
    * for (i=0; i<mat->nru; ++i) { */
  for (i=mat->nru; i>0; --i) {
    write_sparse_dense_block_row_to_buffer(buffer, i-1, mat->A, mat->B, mat->cbl,
        mat->cbr, mat->bs);
	  fwrite(buffer, sizeof(char), strlen(buffer), fh);
    fflush(fh);
  }
  /** write bottom block CD
    * for (i=0; i<mat->nrl; ++i) { */
  for (i=mat->nrl; i>0; --i) {
    write_sparse_dense_block_row_to_buffer(buffer, i-1, mat->C, mat->D, mat->cbl,
        mat->cbr, mat->bs);
	  fwrite(buffer, sizeof(char), strlen(buffer), fh);
    fflush(fh);
  }
  fclose(fh);
  free(buffer);
}

void write_sparse_dense_block_row_to_buffer(char *buffer, const nelts_t idx,
    const sb_fl_t *A, const dbm_fl_t *B, const nelts_t cbl, const nelts_t cbr,
    const bi_t bs)
{
  nelts_t i, j;

  nelts_t nc  = A->ncols + B->ncols;

  memset(buffer, '0', nc);
  buffer[nc]    = '\n';
  buffer[nc+1]  = '\0';
  
  ri_t rbi = idx / bs;
  ri_t rib = idx % bs;

  /** row in A, lefthand side */
  for (i=0; i<cbl; ++i) {
    if (A->blocks[rbi][i].val != NULL) {
      for (j=0; j<A->blocks[rbi][i].sz[rib]; ++j) {
        buffer[A->ncols-1 - (i*bs +  A->blocks[rbi][i].pos[rib][j])]  = '1';
      }
    }
  }
   /** row in B, righthand side */
  for (i=0; i<cbr; ++i) {
    if (B->blocks[rbi][i].val != NULL) {
      for (j=0; j<bs; ++j) {
        if (B->blocks[rbi][i].val[(rib*bs) + j] != 0) {
          buffer[A->ncols + i*bs + j]  = '1';
        }
      }
    }
  }
}

void write_upper_part_row_to_buffer(char *buffer, const nelts_t idx,
    const mat_t *mat)
{
  nelts_t i, j;

  nelts_t nc  = mat->ncl + mat->ncr;

  /** need inverse index since rows of B has to be written in inverse order */
  nelts_t iidx  = mat->nru-1-idx;

  memset(buffer, '0', nc);
  buffer[nc]    = '\n';
  buffer[nc+1]  = '\0';
  
  ri_t rbi = iidx / mat->bs;
  ri_t rib = iidx % mat->bs;

  /** row in A, lefthand side
    * has only 1 at diagonal at position idx not iidx! */
  buffer[idx] = '1';
  /** row in B, righthand side */
  for (i=0; i<mat->cbr; ++i) {
    if (mat->B->blocks[rbi][i].val != NULL) {
      for (j=0; j<mat->bs; ++j) {
        if (mat->B->blocks[rbi][i].val[rib*mat->bs + j] != 0) {
          buffer[mat->ncl + i*mat->bs + j]  = '1';
        }
      }
    }
  }
}

void write_lower_part_row_to_buffer(char *buffer, const nelts_t idx,
    const mat_t *mat)
{
  nelts_t i;

  nelts_t nc  = mat->ncl + mat->ncr;

  memset(buffer, '0', nc);
  buffer[nc]    = '\n';
  buffer[nc+1]  = '\0';
  
  /** row in C, lefthand side is zero! */

  /** row in B, righthand side.
    * if idx > D->rank just keep buffer ={0} */
  if (idx < mat->DR->rank) {
    if (mat->DR->row[idx]->piv_val != NULL) {
      for (i=0; i<mat->ncr; ++i) {
        if (mat->DR->row[idx]->piv_val[i] != 0) {
          buffer[mat->ncl + i] ='1';
        }
      }
    }
  }
}

void print_basis(const gb_t *basis, const poly_t *fb)
{
  nelts_t i, j;
  nvars_t k;
  exp_t *exp = NULL;

  /** depending on the chosen order we have different start and end points in the
    * exponent vectors:
    * for DRL (ord==0) we have stored the exponents in reverse order,
    * for LEX (prd==1) we kept the usual order, but have possibly homogenized */
  const nvars_t ev_start  = basis->ord == 0 ? basis->rnv : 0;
  const nvars_t ev_end    = basis->ord == 0 ? 0 : basis->rnv;

  for (k=0; k<basis->rnv-1; ++k) {
    printf("%s, ", basis->vnames[k]);
  }
  printf("%s\n", basis->vnames[k]);
  printf("%u\n", basis->mod);
  /** prints groebner basis */
  nelts_t bs  = basis->load - basis->st - basis->nred;

  for (i=0; i<bs; ++i) {
    if (fb[i].red == 0) {
      /** we do the first term differently, since we do not have a "+" in front of
        * it */
      printf("%u", fb[i].cf[0]);
      exp = ht->exp[fb[i].eh[0]];
      switch (basis->ord) {
        case 0:
          for (k=ev_start; k>ev_end; --k) {
            if (exp[k] != 0) {
              printf("*%s^%u", basis->vnames[basis->rnv-k], exp[k]);
            }
          }
          break;
        case 1:
          for (k=ev_start; k<ev_end; ++k) {
            if (exp[k] != 0) {
              printf("*%s^%u", basis->vnames[k], exp[k]);
            }
          }
          break;
        default:
          abort();
      }
      for (j=1; j<fb[i].nt; ++j) {
        printf("+%u", fb[i].cf[j]);
        exp = ht->exp[fb[i].eh[j]];
        switch (basis->ord) {
          case 0:
            for (k=ev_start; k>ev_end; --k) {
              if (exp[k] != 0) {
                printf("*%s^%u", basis->vnames[basis->rnv-k], exp[k]);
              }
            }
            break;
          case 1:
            for (k=ev_start; k<ev_end; ++k) {
              if (exp[k] != 0) {
                printf("*%s^%u", basis->vnames[k], exp[k]);
              }
            }
            break;
          default:
            abort();
        }
      }
      printf("\n");
    }
  }
}

void print_basis_in_singular_format(const gb_t *basis, const poly_t *fb)
{
  nelts_t i, j;
  nvars_t k;
  exp_t *exp = NULL;
  /** prints ring */
  printf("ring r = %u, (%s", basis->mod, basis->vnames[0]);
  for (k=1; k<basis->rnv; ++k)
    printf(",%s", basis->vnames[k]);
  switch (basis->ord) {
    case 0:
      printf("), dp;\n");
      break;
    case 1:
      printf("), lp;\n");
      break;
    default:
      abort();
  }

  /** depending on the chosen order we have different start and end points in the
    * exponent vectors:
    * for DRL (ord==0) we have stored the exponents in reverse order,
    * for LEX (prd==1) we kept the usual order, but have possibly homogenized */
  const nvars_t ev_start  = basis->ord == 0 ? basis->rnv : 0;
  const nvars_t ev_end    = basis->ord == 0 ? 0 : basis->rnv;

  /** prints input ideal */
  printf("ideal i;\n");
  for (i=1; i<basis->st; ++i) {
    printf("i[%u]=", i);
    /** we do the first term differently, since we do not have a "+" in front of
      * it */
    printf("%u", basis->cf[i][0]);
    exp = ht->exp[basis->eh[i][0]];
    switch (basis->ord) {
      case 0:
        for (k=ev_start; k>ev_end; --k) {
          if (exp[k] != 0) {
            printf("*%s^%u", basis->vnames[basis->rnv-k], exp[k]);
          }
        }
        break;
      case 1:
        for (k=ev_start; k<ev_end; ++k) {
          if (exp[k] != 0) {
            printf("*%s^%u", basis->vnames[k], exp[k]);
          }
        }
        break;
      default:
        abort();
    }
    for (j=1; j<basis->nt[i]; ++j) {
      printf("+%u", basis->cf[i][j]);
      exp = ht->exp[basis->eh[i][j]];
      switch (basis->ord) {
        case 0:
          for (k=ev_start; k>ev_end; --k) {
            if (exp[k] != 0) {
              printf("*%s^%u", basis->vnames[basis->rnv-k], exp[k]);
            }
          }
          break;
        case 1:
          for (k=ev_start; k<ev_end; ++k) {
            if (exp[k] != 0) {
              printf("*%s^%u", basis->vnames[k], exp[k]);
            }
          }
          break;
        default:
          abort();
      }
    }
    printf(";\n");
  }

  /** prints groebner basis sorted w.r.t. the given monomial order */
  printf("ideal g;\n");
  /** we need to get the basis size without the probably removed elements due to
    * saturation */
  nelts_t bs  = basis->load - basis->st - basis->nred;
  int ctr = 0;
  for (i=0; i<bs; ++i) {
    if (fb[i].red == 0) {
      printf("g[%u]=", ctr+1);
      ctr++;
      /** we do the first term differently, since we do not have a "+" in front of
        * it */
      printf("%u", fb[i].cf[0]);
      exp = ht->exp[fb[i].eh[0]];
      switch (basis->ord) {
        case 0:
          for (k=ev_start; k>ev_end; --k) {
            if (exp[k] != 0) {
              printf("*%s^%u", basis->vnames[basis->rnv-k], exp[k]);
            }
          }
          break;
        case 1:
          for (k=ev_start; k<ev_end; ++k) {
            if (exp[k] != 0) {
              printf("*%s^%u", basis->vnames[k], exp[k]);
            }
          }
          break;
        default:
          abort();
      }
      for (j=1; j<fb[i].nt; ++j) {
        printf("+%u", fb[i].cf[j]);
        exp = ht->exp[fb[i].eh[j]];
        switch (basis->ord) {
          case 0:
            for (k=ev_start; k>ev_end; --k) {
              if (exp[k] != 0) {
                printf("*%s^%u", basis->vnames[basis->rnv-k], exp[k]);
              }
            }
            break;
          case 1:
            for (k=ev_start; k<ev_end; ++k) {
              if (exp[k] != 0) {
                printf("*%s^%u", basis->vnames[k], exp[k]);
              }
            }
            break;
          default:
            abort();
        }
      }
      printf(";\n");
    }
  }
}

void inverse_coefficient(cf_t *x, const cf_t modulus)
{
  assert(*x);
  if ( *x == 1 ) return ;
  assert((int32_t)modulus > 0);
  int32_t u1 = 1, u2 = 0;
  int32_t v1 = 0, v3 = (int32_t)modulus;
  int32_t u3 = (int32_t)*x, v2 = 1;
  while (v3 != 0) {
    int32_t q  = u3 / v3;
    int32_t t1 = u1 - v1 * q;
    u1  = v1; v1  = t1;

    int32_t t3 = u3 - v3 * q;
    u3  = v3; v3  = t3;

    int32_t t2 = u2 - v2 * q;
    u2  = v2; v2  = t2;
  }
  if (u1 < 0) {
    u1  +=  (int32_t)modulus;
    /* check_inverse(*x,u1,modulus); */
    *x  =   (cf_t)u1;
    return;
  }
  if (u1 > (int32_t)modulus) {
    u1  -=  (int32_t)modulus;
    /* check_inverse(*x,u1,modulus); */
    *x  =   (cf_t) u1;
    return;
  }
  /* check_inverse(*x,u1,modulus); */
  *x  = (cf_t)u1;
  return;
}
