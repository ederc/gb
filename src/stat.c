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
    stat_t *st  = (stat_t *)calloc(1, sizeof(stat_t));

    return st;
}

static void print_initial_statistics(
        const stat_t *st
        )
{

    printf("\n--------------- INPUT DATA ---------------\n");
    printf("#variables             %11d\n", st->nvars);
    printf("#equations             %11d\n", st->ngens);
    printf("field characteristic   %11d\n", st->fc);
    printf("homogeneous input?     %11d\n", st->homogeneous);
    if (st->mo == 0) {
        printf("monomial order                 DRL\n");
    }
    if (st->mo == 1) {
        printf("monomial order                 LEX\n");
    }
    if ((st->mo != 0) && (st->mo != 1)) {
        printf("monomial order           DONT KNOW\n");
    }
    if (st->reset_ht == 2147483647) {
        printf("basis hash table resetting     OFF\n");
    } else {
        printf("basis hash table resetting  %6d\n", st->reset_ht);
    }
    printf("linear algebra option  %11d\n", st->laopt);
    printf("intial hash table size %11d (2^%d)\n",
            (int32_t)pow(2,st->init_hts), st->init_hts);
    if (st->mnsel == 2147483647) {
        printf("max pair selection             ALL\n");
    } else {
        printf("max pair selection     %11d\n", st->mnsel);
    }
    printf("#threads               %11d\n", st->nthrds);
    printf("info level             %11d\n", st->info_level);
    printf("generate pbm files     %11d\n", st->gen_pbm_file);
    printf("------------------------------------------\n");
}

static void print_final_statistics(
        const stat_t * const st
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
    if (st->reduce_gb == 1) {
        printf("reduce gb    %15.3f sec %5.1f%%\n",
                st->reduce_gb_rtime,
                (double)100*(double)st->reduce_gb_rtime
                / (double)(st->overall_rtime));
    }
    if (st->reset_ht != 2147483647) {
        printf("rht          %15.3f sec %5.1f%%\n",
                st->rht_rtime,
                (double)100*(double)st->rht_rtime
                / (double)(st->overall_rtime));
    }
    printf("-----------------------------------------\n");
    printf("\n---------- COMPUTATIONAL DATA -----------\n");
    printf("size of basis      %9d\n", st->size_basis);
    printf("#terms in basis    %9ld\n", st->nterms_basis);
    printf("#pairs reduced     %9ld\n", st->num_pairsred);
    printf("#GM criterion      %9ld\n", st->num_gb_crit);
    printf("#redundant         %9ld\n", st->num_redundant);
    printf("#reset bht         %9ld\n", st->num_rht);
    printf("#rows reduced      %9ld\n", st->num_rowsred);
    printf("#zero reductions   %9ld\n", st->num_zerored);
    printf("maximal uht size           2^%d\n",
            (int32_t)(ceil(log((double)st->max_uht_size)/log(2))));
    printf("maximal sht size           2^%d\n",
            (int32_t)(ceil(log((double)st->max_sht_size)/log(2))));
    printf("maximal bht size           2^%d\n",
            (int32_t)(ceil(log((double)st->max_bht_size)/log(2))));
    printf("-----------------------------------------\n\n");
}
