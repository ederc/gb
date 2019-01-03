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
 * \file update.c
 * \brief Update process and pairset handling
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "update.h"

void insert_and_update_spairs(
        ps_t *psl,
        bs_t *bs,
        ht_t *ght,
        ht_t *lht,
        stat_t *st
        )
{
    len_t i, j, l;

    spair_t *ps = psl->p;

    const len_t pl        = psl->ld;
    const len_t bl        = bs->ld;
    const hd_t * const nh = bs->m[bl].h[0];

    reset_hash_table(lht, bl);

    if (is_new_generator_redundant(psl, bs, ght, st)) {
        return;
    }

    hd_t **plcm = generate_new_pairs(psl, lht, bs);
    len_t nl    = pl + bl;

    /* Gebauer-Moeller: check old pairs first */
    /* note: old pairs are sorted by the given spair order */
    for (i = 0; i < pl; ++i) {
        j = ps[i].gen1;
        l = ps[i].gen2;
        if (check_monomial_division(ps[i].lcm, nh)
                && ps[i].lcm->val != plcm[j]->val
                && ps[i].lcm->val != plcm[l]->val
           ) {
            ps[i].lcm = NULL;
        }
    }

    /* sort new pairs by increasing lcm, earlier polys coming first */
    spair_t *pp = ps+pl;
    j = 0;
    for (i = 0; i < bl; ++i) {
        if (pp[i].lcm != NULL) {
            pp[j++] = pp[i];
        }
    }
    qsort(pp, (unsigned long)j, sizeof(spair_t), spair_cmp);
    for (i = 0; i < j; ++i) {
        plcm[i] = pp[i].lcm;
    }
    plcm[j]  = 0;
    const len_t pc  = j;

    j = 0;

    for (; j < pc; ++j) {
        if (plcm[j] == NULL) {
            continue;
        }
        const hd_t * const plcmj = plcm[j];
        i = j+1;
        while (plcm[i] == plcmj) {
            plcm[i] = NULL;
            ++i;
        }
        j = i-1;
        while (i < pc) {
            if (plcm[i] != NULL &&
                    check_monomial_division(plcm[i], plcmj) != 0) {
                plcm[i]  = NULL;
            }
            ++i;
        }
    }

    /* remove useless pairs from pairset */
    j = 0;
    /* old pairs */
    for (i = 0; i < psl->ld; ++i) {
        if (ps[i].lcm == NULL) {
            continue;
        }
        ps[j++] = ps[i];
    }
    if (ght->esz - ght->eld <= nl-psl->ld) {
        enlarge_hash_table(ght);
    }
    /* new pairs, wee need to add the lcm to the global hash table */
    for (i = 0; i < pc; ++i) {
        if (plcm[i] == NULL) {
            continue;
        }
        pp[i].lcm = insert_in_hash_table(
                plcm[i]->exp, ght);
        ps[j++]   = pp[i];
    }
    free(plcm);
    psl->ld = j;
    st->num_gm_crit +=  nl - j;

    /* mark redundant elements in basis */
    check_old_elements_for_redundancy(bs, st);
    bs->ld++;
}

void update_basis(
        ps_t *ps,
        bs_t *bs,
        ht_t *ght,
        ht_t *lht,
        stat_t *st,
        const len_t ne
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    for (i = 0; i < ne; ++i) {
        insert_and_update_spairs(ps, bs, ght, lht, st);
    }

    bs->lo  = bs->ld;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->update_ctime  +=  ct1 - ct0;
    st->update_rtime  +=  rt1 - rt0;
}
