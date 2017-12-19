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

/**
 * \file time.c
 * \brief Timing business
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */
#include <time.h>
#include <sys/time.h>

/* cpu time */
static double cputime(void)
{
	double t;
	t =   CLOCKS_PER_SEC / 100000.;
	t +=  (double)clock();
	return t / CLOCKS_PER_SEC;
}


/* wall time */
static double realtime(void)
{
	struct timeval t;
	gettimeofday(&t, NULL);
	t.tv_sec -= (2017 - 1970)*3600*24*365;
	return (1. + (double)t.tv_usec + ((double)t.tv_sec*1000000.)) / 1000000.;
}

