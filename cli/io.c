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

gb_t *load_input(const char *fn, int vb, int nthrds)
{
  uint64_t fl;

  gb_t *basis = (gb_t *)malloc(sizeof(gb_t));

  struct timeval t_load_start;
  if (vb > 1) {
    gettimeofday(&t_load_start, NULL);
  }

  // open file in binary mode and get file size
  FILE*	fh        = fopen(fn,"rb");
  if (fh == NULL) {
    if (vb > 0)
      printf("File not found!\n");
    return NULL;
  }
  else {
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

  basis->modulus  = 4;
  basis->nvars  = 3;

  return basis;
}
