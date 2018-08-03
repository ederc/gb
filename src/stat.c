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

    st->round_rtime   = 0;
    st->rght_rtime    = 0;
    st->select_rtime  = 0;
    st->symbol_rtime  = 0;
    st->la_rtime      = 0;
    st->update_rtime  = 0;
    st->convert_rtime = 0;

    st->num_pairsred    = 0;
    st->num_gb_crit     = 0;
    st->num_redundant   = 0;
    st->num_duplicates  = 0;
    st->num_rowsred     = 0;
    st->num_zerored     = 0;

    st->max_ht_size       = 0;
    st->num_sdm_found     = 0;
    st->num_not_sdm_found = 0;

    return st;
}
