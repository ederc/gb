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
 * \file basis.h
 * \brief General and global data
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_BASIS_H
#define GB_BASIS_H

#include "data.h"

inline bs_t *initialize_basis(
        const stat_t *st
        )
{
    bs_t *bs  = (bs_t *)malloc(sizeof(bs_t));

    bs->lo  = 0;
    bs->ld  = 0;
    bs->sz  = 2 * st->nr_gens;

    bs->cf  = (void **)malloc((unsigned long)bs->sz * sizeof(void *));
    bs->hd  = (dt_t **)malloc((unsigned long)bs->sz * sizeof(dt_t *));
    bs->lm  = (sdm_t *)malloc((unsigned long)bs->sz * sizeof(sdm_t));
    bs->red = (red_t *)malloc((unsigned long)bs->sz * sizeof(red_t));

    bs->tcf = NULL;

    return bs;
}

inline void free_basis(
        bs_t **bsp
        )
{
    bs_t *bs  = *bsp;
    len_t i;
    if (bs->hd) {
        for (i = 0; i < bs->ld; ++i) {
            free(bs->cf[i]);
            free(bs->hd[i]);
        }
        free(bs->cf);
        bs->cf  = NULL;
        free(bs->hd);
        bs->hd  = NULL;
        free(bs->lm);
        bs->lm  = NULL;
        free(bs->red);
        bs->red = NULL;
    }

    bs    = NULL;
    *bsp  = bs;
}

inline void check_enlarge_basis(
        bs_t *bs,
        len_t added
        )
{
    if (bs->ld+added >= bs->sz) {
        bs->sz  = bs->sz*2 > bs->ld+added ? bs->sz*2 : bs->ld+added;
        bs->cf  = realloc(bs->cf, (unsigned long)bs->sz * sizeof(void *));
        bs->hd  = realloc(bs->hd, (unsigned long)bs->sz * sizeof(dt_t *));
        bs->lm  = realloc(bs->lm, (unsigned long)bs->sz * sizeof(sdm_t));
        bs->red = realloc(bs->red, (unsigned long)bs->sz * sizeof(red_t));
    }
}

inline void check_old_elements_for_redundancy(
        bs_t *bs,
        stat_t *st,
        const ht_t *const ht
        )
{
    len_t i;
    const len_t bl  = bs->ld;
    const hl_t  nh  = bs->hd[bs->ld][3];

    for (i = 0; i < bl; ++i) {
        if (!bs->red[i]
            && check_monomial_division(bs->hd[i][3], nh, ht)) {
            bs->red[i]  = 1;
            st->num_redundant++;
        }
    }
}

#endif
