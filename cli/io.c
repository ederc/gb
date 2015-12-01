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
  char    *unit = "KB";
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

  int ret;
  /*  dummy reading of useless information */
  ret = fscanf(fh, "%1023s", &pid[0]);
  ret = fscanf(fh, "%1023s", &comm[0]);
  ret = fscanf(fh, "%1023s", &state[0]);
  ret = fscanf(fh, "%1023s", &ppid[0]);
  ret = fscanf(fh, "%1023s", &pgrp[0]);
  ret = fscanf(fh, "%1023s", &session[0]);
  ret = fscanf(fh, "%1023s", &tty_nr[0]);
  ret = fscanf(fh, "%1023s", &tpgid[0]);
  ret = fscanf(fh, "%1023s", &flags[0]);
  ret = fscanf(fh, "%1023s", &minflt[0]);
  ret = fscanf(fh, "%1023s", &cminflt[0]);
  ret = fscanf(fh, "%1023s", &majflt[0]);
  ret = fscanf(fh, "%1023s", &cmajflt[0]);
  ret = fscanf(fh, "%1023s", &utime[0]);
  ret = fscanf(fh, "%1023s", &stime[0]);
  ret = fscanf(fh, "%1023s", &cutime[0]);
  ret = fscanf(fh, "%1023s", &cstime[0]);
  ret = fscanf(fh, "%1023s", &priority[0]);
  ret = fscanf(fh, "%1023s", &nice[0]);
  ret = fscanf(fh, "%1023s", &nthrds[0]);
  ret = fscanf(fh, "%1023s", &itrealvalue[0]);
  ret = fscanf(fh, "%1023s", &starttime[0]);

  /*  get real memory information */
  ret = fscanf(fh, "%lu", &_vms);
  ret = fscanf(fh, "%ld", &_rss);

  /*  close file */
  fclose(fh);

  /*  TODO: How to read /proc/self/stat ??? */

  vms = _vms / 1024.0;
  rss = _rss * page_size_kb;

  /*  MB ? */
  if (vms > 1024) {
    vms   = vms/1024.0;
    rss   = rss/1024.0;
    unit  = "MB";
  }
  /*  GB ? */
  if (vms > 1024) {
    vms   = vms/1024.0;
    rss   = rss/1024.0;
    unit  = "GB";
  }
  /*  TB ? Just joking! */
  if (vms > 1024) {
    vms   = vms/1024.0;
    rss   = rss/1024.0;
    unit  = "TB";
  }
  printf("MMRY\tRSS - %.3f %s | VMS - %.3f %s\n", rss, unit, vms, unit);
}

inline int get_number_of_terms(const char *line)
{
  const char add_splicer  = '+';
  char *tmp;
  int nterms  = 1;
  tmp = strchr(line, add_splicer);
  while (tmp != NULL) {
    nterms++;
    tmp = strchr(tmp+1, add_splicer);
  }

  return nterms;
}

inline exp_t get_exponent(const char *term, const char *var_name)
{
  const char mult_splicer = '*';
  const char exp_splicer  = '^';
  exp_t exp = 0;

  char *var = strstr(term, var_name);
  if (var != NULL) {
    // if the next variable follows directly => exp = 1
    if (strncmp(&mult_splicer, var+strlen(var_name), 1) == 0) {
      exp = 1;
    } else {
      // if there follows an exp symbol "^"
      if (strncmp(&exp_splicer, var+strlen(var_name), 1) == 0) {
        char exp_str[100];
        char *mult_pos;
        mult_pos  = strchr(var, mult_splicer); 
        if (mult_pos != NULL) {
          int exp_len = (int)(mult_pos - (var+strlen(var_name)) - 1);
          memcpy(exp_str, var+strlen(var_name)+1, exp_len);
          exp_str[exp_len] = '\0';
          exp = (exp_t)atoi(exp_str);
        } else { // no further variables in this term
          int exp_len = (int)((var+strlen(var)) + 1 - (var+strlen(var_name)) - 1);
          memcpy(exp_str, var+strlen(var_name)+1, exp_len);
          exp = (exp_t)atoi(exp_str);
          exp_str[exp_len] = '\0';
        }
      } else { // we are at the last variable with exp = 1
        exp = 1;
      }
    }
  }
  return exp;
}

inline void get_term(const char *line, char **prev_pos,
    char **term)
{
  // note that maximal term length we handle
  const char add_splicer  = '+';

  char *curr_pos  = strchr(*prev_pos, add_splicer);
  if (curr_pos != NULL) {
    int prev_idx  = (int)(*prev_pos - line);
    int curr_idx  = (int)(curr_pos - line);
    int term_diff = curr_idx - prev_idx;
    memcpy(*term, *prev_pos, term_diff);
    (*term)[term_diff]  = '\0';
    *prev_pos           = curr_pos+1;
  } else { // we are at the last term
    int prev_idx  = (int)(*prev_pos - line);
    int term_diff = strlen(line) + 1 - prev_idx;
    memcpy(*term, *prev_pos, strlen(line)+1);
    (*term)[term_diff]  = '\0';
  }
}

nvars_t get_nvars(const char *fn)
{
  FILE *fh  = fopen(fn,"r");
  // load lines and store data
  const size_t max_line_size  = 1000;
  char *line  = (char *)malloc(max_line_size * sizeof(char));

  // get first line (variables)
  const char comma_splicer  = ',';

  // get number of variables
  nvars_t nvars = 1; // number of variables is number of commata + 1 in first line
  if (fgets(line, max_line_size, fh) != NULL) {
    char *tmp = strchr(line, comma_splicer);
    while (tmp != NULL) {
      nvars++;
      tmp = strchr(tmp+1, comma_splicer);
    }
  } else {
    printf("Bad file format.\n");
    return NULL;
  }

  return nvars;
}

gb_t *load_input(const char *fn, nvars_t nvars, mp_cf4_ht_t *ht, int vb, int nthrds)
{
  uint64_t fl;

  hash_t i, j, k;

  gb_t *basis = (gb_t *)malloc(sizeof(gb_t));

  struct timeval t_load_start;
  if (vb > 1) {
    gettimeofday(&t_load_start, NULL);
  }

  // open file in binary mode and get file size
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

  basis->fs   = (double) fl / 1024;
  basis->fsu  = "KB";
  if (basis->fs > 1000) {
    basis->fs   = basis->fs / 1024;
    basis->fsu  = "MB";
  }
  if (basis->fs > 1000) {
    basis->fs   = basis->fs / 1024;
    basis->fsu  = "GB";
  }

  // get number of lines in file:
  // number of lines - 2 is number of generators in input system
  int nlines  = 0;
  char buf[1000];
  fh  = fopen(fn,"r");
  while(fgets(buf, sizeof(buf), fh) != NULL)
    nlines++;
  fclose(fh);

  // #generators of the input system = nlines - 2 since the first line has the
  // variable names and second line is the field modulus
  basis->load = nlines-2;

  printf("bload = %u\n",basis->load);

  // now initialize the size of the basis 3 times as big as the input system
  basis->size = 3 * basis->load;
  basis->nt   = (nelts_t *)malloc(basis->size * sizeof(nelts_t));
  basis->deg  = (deg_t *)malloc(basis->size * sizeof(deg_t));
  basis->cf   = (coeff_t **)malloc(basis->size * sizeof(coeff_t *));
  basis->eh   = (hash_t **)malloc(basis->size * sizeof(hash_t *));

  fh  = fopen(fn,"r");
  // load lines and store data
  const size_t max_line_size  = 1000;
  char *line  = (char *)malloc(max_line_size * sizeof(char));

  // get first line (variables)
  const char comma_splicer  = ',';

  // we already know the number of variables
  basis->nvars  = nvars;

  char *tmp_vname, *tmp_vname_old;
  int index, index_old;
  // allocate memory for storing variable names
  basis->vnames = (char **)malloc(basis->nvars * sizeof(char *));
  if (fgets(line, max_line_size, fh) != NULL) {
    index_old  = 0;
    tmp_vname_old = line;
    for (i=0; i<basis->nvars-1; ++i) {
      tmp_vname = strchr(tmp_vname_old, comma_splicer);
      index = (int)(tmp_vname - line);
      basis->vnames[i]  = (char *)malloc(50*sizeof(char));
      memcpy(basis->vnames[i], line+index_old, index-index_old);
      basis->vnames[i][index-index_old] = '\0';
      basis->vnames[i]  = realloc(basis->vnames[i], index-index_old+1);
      tmp_vname_old = tmp_vname + 1;
      index_old     = index + 1;
    }
    // now do the last variable separately
    basis->vnames[i]  = (char *)malloc(50*sizeof(char));
    memcpy(basis->vnames[i], tmp_vname_old, strlen(line)+1);
    basis->vnames[i][strlen(line)-index] = '\0';
    basis->vnames[i]  = realloc(basis->vnames[i], strlen(line)-index+1);
    char *tmp_end;
    // trim variable names
    for (i=0; i<basis->nvars; ++i) {
      while (isspace(*basis->vnames[i]))
        basis->vnames[i]++;
      if (*basis->vnames[i] != 0) {
        tmp_end = basis->vnames[i];
        while (*tmp_end)
          tmp_end++;
        do {
          tmp_end--;
        } while (isspace(*tmp_end));
        tmp_end++;
        *tmp_end = '\0';
      }
    }

    for (i=0; i<basis->nvars; ++i) {
      printf(basis->vnames[i]);
    }
  } else {
    printf("Bad file format.\n");
    return NULL;
  }
  // get second line (modulus)
  if (fgets(line, max_line_size, fh) != NULL) {
    int64_t tmp_mod = atol(line);
    if (tmp_mod > 0) {
      basis->modulus  = (mod_t)tmp_mod;
    } else {
      printf("Bad file format.\n");
      return NULL;
    }
  } else {
    printf("Bad file format.\n");
    return NULL;
  }

  const char add_splicer  = '+';
  const char mult_splicer = '*';
  const char exp_splicer  = '^';

  char *prev_pos;
  char *term  = (char *)malloc(200 * sizeof(char));;
  int nterms;
  exp_t *exp  = (exp_t *)malloc(basis->nvars * sizeof(exp_t));
  deg_t deg, max_deg;

  // get all remaining lines, i.e. generators
  for (i=0; i<basis->load; ++i) {
    if (fgets(line, max_line_size, fh) != NULL) {
      // get number of terms first
      nterms  = get_number_of_terms(line);
      printf("nterms %d\n",nterms);

      // allocate memory for all terms
      basis->cf[i]  = (coeff_t *)malloc(nterms * sizeof(coeff_t));
      basis->eh[i]  = (hash_t *)malloc(nterms * sizeof(hash_t));
      printf(line);
      prev_pos  = line;
      max_deg   = 0;
      // next: go through line, term by term
      for (j=0; j<nterms; ++j) {
        memset(exp, 0, basis->nvars * sizeof(exp_t));
        deg = 0;
        get_term(line, &prev_pos, &term);
        printf("get term ---> ");
        printf(term);
        // get coefficient first
        if (term != NULL)
          basis->cf[i][j] = (coeff_t)atoi(term);
        printf("term: ");
        printf(term);
        printf(" || coeff: ");
        printf("%u\n",basis->cf[i][j]);

        // now loop over variables of term
        for (k=0; k<basis->nvars; ++k) {
          exp[k]  =   get_exponent(term, basis->vnames[k]);
          deg     +=  exp[k];
        }
        // hash exponent and store degree
        basis->eh[i][j] = check_in_hash_table(exp, ht);
        max_deg = max_deg > deg ? max_deg : deg;
        printf("eh[%u][%u] = %u --> %lu\n",i,j,basis->eh[i][j], ht->lut[basis->eh[i][j]]);
      }
      basis->deg[i] = max_deg;
    }
  }
  free(line);

  return basis;
}
