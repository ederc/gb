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

#include "stat.h"


stat_t *initialize_statistics(
    void
    )
{
    stat_t *st  = (stat_t *)malloc(sizeof(stat_t));

    st->info_level    = 0;
    st->mon_order     = 0;
    st->field_char    = 0;
    st->la_variant    = 0;
    st->regen_ht      = 0;
    st->nthrds        = 0;
    st->nr_vars       = 0;
    st->max_nr_pairs  = 0;
    st->nr_gens       = 0;
    st->init_ht_sz    = 0; 
    st->cf_sz         = 0;

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

    st->num_pairsred  = 0;
    st->num_gm_crit   = 0;
    st->num_redundant = 0;
    st->num_rowsred   = 0;
    st->num_zerored   = 0;

    st->max_ht_size   = 0;
    st->len_output    = 0;
    st->size_basis    = 0;
    st->num_matrices  = 0;

    return st;
}

void print_round_statistics_header(
        void
        )
{
        printf("\ndeg     sel   pairs        mat          density \
          new data             time(rd)\n");
        printf("-------------------------------------------------\
----------------------------------------\n");
}

void print_round_statistics_footer(
        void
        )
{
        printf("-------------------------------------------------\
----------------------------------------\n");
}


void print_initial_statistics(
        const stat_t *st
        )
{
    printf("\n--------------- INPUT DATA ---------------\n");
    printf("#variables             %11d\n", st->nr_vars);
    printf("#equations             %11d\n", st->nr_gens);
    printf("field characteristic   %11d\n", st->field_char);
    if (st->mon_order == 0) {
        printf("monomial order                 DRL\n");
    }
    if (st->mon_order == 1) {
        printf("monomial order                 LEX\n");
    }
    if ((st->mon_order != 0) && (st->mon_order != 1)) {
        printf("monomial order           DONT KNOW\n");
    }
    printf("linear algebra option  %11d\n", st->la_variant);
    printf("intial hash table size %11d (2^%d)\n",
            (int32_t)pow(2,st->init_ht_sz), st->init_ht_sz);
    printf("regenerate hash table after %6d step(s)\n", st->regen_ht);
    printf("max pair selection     %11d\n", st->max_nr_pairs);
    printf("#threads               %11d\n", st->nthrds);
    printf("info level             %11d\n", st->info_level);
    printf("------------------------------------------\n");
}

void print_final_statistics(
        const stat_t *st,
        const hl_t eld,
        const hl_t elld
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
            (st->len_output-st->size_basis-1)/(1+st->nr_vars));
    printf("#pairs reduced     %9ld\n", st->num_pairsred);
    printf("#GM criterion      %9ld\n", st->num_gm_crit);
    printf("#redundant         %9ld\n", st->num_redundant);
    printf("#matrices reduced  %9d\n", st->num_matrices);
    printf("#rows reduced      %9ld\n", st->num_rowsred);
    printf("#zero reductions   %9ld\n", st->num_zerored);
    printf("#global hash table %9d <= 2^%d\n",
            eld, (int32_t)((ceil(log(eld)/log(2)))));
    printf("#local hash table  %9d <= 2^%d\n",
            elld, (int32_t)(ceil(log(elld)/log(2))));
    printf("maximal ht size         2^%d\n",
            (int32_t)(ceil(log(st->max_ht_size)/log(2))));
    printf("-----------------------------------------\n\n");
}
