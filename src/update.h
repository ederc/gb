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
 * \file update.h
 * \brief Update process and pairset handling
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#ifndef GB_UPDATE_H
#define GB_UPDATE_H

#include "data.h" 
#include "time.h"
#include "stat.h"
#include "basis.h"
#include "hash.h"

void update_basis(
        ps_t *ps,
        bs_t *bs,
        ht_t *ght,
        ht_t *lht,
        stat_t *st,
        const len_t ne
        );

void insert_and_update_spairs(
        ps_t *psl,
        bs_t *bs,
        ht_t *ght,
        ht_t *lht,
        stat_t *st
        );

static inline ps_t *initialize_pairset(
        const stat_t *st
        )
{
    ps_t *ps  = (ps_t *)malloc(sizeof(ps_t));
    ps->ld    = 0;
    ps->sz    = 10 * st->nr_gens;
    ps->mnsel = st->max_nr_pairs;
    ps->p     = (spair_t *)malloc((unsigned long)ps->sz * sizeof(spair_t));

    return ps;
}

static inline void check_enlarge_pairset(
        ps_t *ps,
        const len_t bsl,
        const len_t np
        )
{
    len_t i;
    /* compute number of new pairs we need to handle at most */
    len_t sum  = bsl * np;
    for (i = 1; i < np; ++i) {
        sum = sum + i;
    }
    if (ps->ld+sum >= ps->sz) {
        ps->sz  = ps->sz*2 > ps->ld+sum ? ps->sz*2 : ps->ld+sum;
        ps->p   = realloc(ps->p, (unsigned long)ps->sz * sizeof(spair_t));
        memset(ps->p+ps->ld, 0,
                (unsigned long)(ps->sz-ps->ld) * sizeof(spair_t));
    }
}

static inline void free_pairset(
        ps_t **psp
        )
{
    ps_t *ps  = *psp;
    if (ps->p) {
        free(ps->p);
        ps->p   = NULL;
        ps->ld  = 0;
        ps->sz  = 0;
    }
    free(ps);
    ps  = NULL;
    *psp  = ps;
}

static inline int is_new_generator_redundant(
        ps_t *psl,
        bs_t *bs,
        ht_t *ht,
        stat_t *st
        )
{
    len_t i;
    const len_t bl  = bs->ld;
    const hd_t *nh  = bs->m[bl].h[0];
    
    for (i = bs->lo; i < bl; ++i) {
        if (!bs->red[i] 
            && check_monomial_division(nh, bs->m[i].h[0])) {
            /* printf("Mark polynomial %d unnecessary for new pairs\n", bload); */
            spair_t p = psl->p[psl->ld];
            p.gen1 = i;
            p.gen2 = bl;
            if (ht->eld >= ht->esz) {
                enlarge_hash_table(ht);
            }
            p.lcm  = get_lcm(bs->m[i].h[0], nh, ht);
            bs->red[bl] = 1;
            st->num_redundant++;
            bs->ld++;
            psl->ld++;
            return 1;
        }
    }
    return 0;
}

static inline hd_t **generate_new_pairs(
        ps_t *psl,
        ht_t *lht,
        const bs_t *const bs
        )
{
    len_t i, k;

    const len_t bl  = bs->ld;
    const hd_t *nh  = bs->m[bl].h[0];

    spair_t *ps = psl->p;
    hd_t **plcm = (hd_t **)malloc((unsigned long)(bl+1) * sizeof(hd_t *));

    /* create all possible new pairs */
    for (i = 0, k = psl->ld; i < bl; ++i, ++k) {
        /* b = (hl_t *)((long)bs[i] & bmask); */
        ps[k].gen1  = i;
        ps[k].gen2  = bl;

        ps[k].lcm = get_lcm(bs->m[i].h[0], nh, lht);

        if (bs->red[i]) {
            ps[k].lcm = NULL; /* redundant pair */
        } else {
            if (lcm_equals_multiplication(bs->m[i].h[0], nh, ps[k].lcm)) {
                ps[k].lcm = NULL; /* criterion */
            }
        }
        plcm[i] = ps[k].lcm;
    }

    return plcm;
}

#endif