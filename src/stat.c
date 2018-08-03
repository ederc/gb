/* gb: Gröbner Basis
 * Copyright (C) 2018 Christian Eder <ederc@mathematik.uni-kl.de>
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
 * \file statistics.c
 * \brief Global data for covering computational statistics
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"


static stat_t *initialize_statistics(
    void
    )
{
    stat_t *st  = (stat_t *)malloc(sizeof(stat_t));

    st->round_ctime   = 0;
    st->rght_ctime    = 0;
	  st->select_ctime  = 0;
    st->symbol_ctime  = 0;
    st->la_ctime      = 0;
    st->update_ctime  = 0;
    st->convert_ctime = 0;
    st->overall_ctime = 0;

    st->round_rtime   = 0;
    st->rght_rtime    = 0;
    st->select_rtime  = 0;
    st->symbol_rtime  = 0;
    st->la_rtime      = 0;
    st->update_rtime  = 0;
    st->convert_rtime = 0;
    st->overall_rtime = 0;

    st->num_pairsred    = 0;
    st->num_gb_crit     = 0;
    st->num_redundant   = 0;
    st->num_duplicates  = 0;
    st->num_rowsred     = 0;
    st->num_zerored     = 0;

    st->max_ht_size       = 0;
    st->num_sdm_found     = 0;
    st->num_not_sdm_found = 0;
    st->len_output        = 0;
    st->size_basis        = 0;

    return st;
}

static void print_final_statistics(
        const stat_t *st
        )
{
    printf("\n---------------- TIMINGS ---------------\n");
    printf("overall      %15.3f sec\n", st->overall_rtime);
    printf("overall(cpu) %15.3f sec\n", st->overall_ctime);
    printf("select       %15.3f sec %5.1f%%\n",
            st->select_rtime,
            (double)100*(double)st->select_rtime
            / (double)(st->overall_rtime));
    printf("symbol       %15.3f sec %5.1f%%\n",
            st->symbol_rtime,
            (double)100*(double)st->symbol_rtime
            / (double)(st->overall_rtime));
    printf("update       %15.3f sec %5.1f%%\n",
            st->update_rtime,
            (double)100*(double)st->update_rtime
            / (double)(st->overall_rtime));
    printf("convert      %15.3f sec %5.1f%%\n",
            st->convert_rtime,
            (double)100*(double)st->convert_rtime
            / (double)(st->overall_rtime));
    printf("rght         %15.3f sec %5.1f%%\n",
            st->rght_rtime,
            (double)100*(double)st->rght_rtime
            / (double)(st->overall_rtime));
    printf("la           %15.3f sec %5.1f%%\n",
            st->la_rtime,
            (double)100*(double)st->la_rtime
            / (double)(st->overall_rtime));
    printf("-----------------------------------------\n");
    printf("\n---------- COMPUTATIONAL DATA -----------\n");
    printf("size of basis      %9d\n", st->size_basis);
    printf("#terms in basis    %9ld\n",
            (st->len_output-st->size_basis-1)/(1+nvars));
    printf("#pairs reduced     %9ld\n", st->num_pairsred);
    printf("#GM criterion      %9ld\n", st->num_gb_crit);
    printf("#redundant         %9ld\n", st->num_redundant);
    printf("#rows reduced      %9ld\n", st->num_rowsred);
    printf("#zero reductions   %9ld\n", st->num_zerored);
    printf("#global hash table %9d <= 2^%d\n",
            eld, (int32_t)((ceil(log(eld)/log(2)))));
    printf("#local hash table  %9d <= 2^%d\n",
            elld, (int32_t)(ceil(log(elld)/log(2))));
    printf("maximal ht size         2^%d\n",
            (int32_t)(ceil(log(st->max_ht_size)/log(2))));
    printf("#no reducer found  %9ld\n",
            st->num_sdm_found + st->num_not_sdm_found);
    printf("sdm findings       %8.3f%% \n",
            (double)100*(double)st->num_sdm_found/
            (double)(st->num_sdm_found + st->num_not_sdm_found));
    printf("-----------------------------------------\n\n");
}
