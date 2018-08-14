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
 * \file stat.h
 * \brief Global data for covering computational statistics
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_STAT_H
#define GB_STAT_H

#include "data.h"

stat_t *initialize_statistics(void);

void print_round_statistics_header(void);

void print_initial_statistics(const stat_t *st);

void print_final_statistics(const stat_t *st);
#endif
