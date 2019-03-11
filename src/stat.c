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
	  st->select_ctime  = 0;
    st->symbol_ctime  = 0;
    st->la_ctime      = 0;
    st->update_ctime  = 0;
    st->convert_ctime = 0;
    st->overall_ctime = 0;
    st->rht_ctime     = 0;

    st->round_rtime   = 0;
    st->select_rtime  = 0;
    st->symbol_rtime  = 0;
    st->la_rtime      = 0;
    st->update_rtime  = 0;
    st->convert_rtime = 0;
    st->overall_rtime = 0;
    st->rht_rtime     = 0;

    st->num_pairsred  = 0;
    st->num_gb_crit   = 0;
    st->num_redundant = 0;
    st->num_rowsred   = 0;
    st->num_zerored   = 0;

    st->reset_ht      = 0;
    st->current_rd    = 0;
    st->current_deg   = 0;
    st->max_ht_size   = 0;
    st->len_output    = 0;
    st->size_basis    = 0;

    st->info_level    = 0;
    st->gen_pbm_file  = 0;

    return st;
}

static void print_initial_statistics(
        const stat_t *st
        )
{

    printf("\n--------------- INPUT DATA ---------------\n");
    printf("#variables             %11d\n", nvars);
    printf("#equations             %11d\n", ngens);
    printf("field characteristic   %11d\n", fc);
    if (mo == 0) {
        printf("monomial order                 DRL\n");
    }
    if (mo == 1) {
        printf("monomial order                 LEX\n");
    }
    if ((mo != 0) && (mo != 1)) {
        printf("monomial order           DONT KNOW\n");
    }
    if (st->reset_ht == 2147483647) {
        printf("basis hash table resetting     OFF\n");
    } else {
        printf("basis hash table resetting  %6d\n", st->reset_ht);
    }
    printf("linear algebra option  %11d\n", laopt);
    printf("intial hash table size %11d (2^%d)\n",
            (int32_t)pow(2,htes), htes);
    if (mnsel == 2147483647) {
        printf("max pair selection             ALL\n");
    } else {
        printf("max pair selection     %11d\n", mnsel);
    }
    printf("#threads               %11d\n", nthrds);
    printf("info level             %11d\n", st->info_level);
    printf("generate pbm files     %11d\n", st->gen_pbm_file);
    printf("------------------------------------------\n");
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
    printf("la           %15.3f sec %5.1f%%\n",
            st->la_rtime,
            (double)100*(double)st->la_rtime
            / (double)(st->overall_rtime));
    if (rht != 2147483647) {
    printf("rht          %15.3f sec %5.1f%%\n",
            st->rht_rtime,
            (double)100*(double)st->rht_rtime
            / (double)(st->overall_rtime));
    }
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
    printf("#update hash table %9d <= 2^%d\n",
            euld, (int32_t)(ceil(log(euld)/log(2))));
    printf("maximal ht size         2^%d\n",
            (int32_t)(ceil(log((double)st->max_ht_size)/log(2))));
    printf("-----------------------------------------\n\n");
}
