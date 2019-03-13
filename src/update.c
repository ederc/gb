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

#include "data.h" 

static ps_t *initialize_pairset(
        void
        )
{
    ps_t *ps  = (ps_t *)malloc(sizeof(ps_t));
    ps->ld  = 0;
    ps->sz  = 192;
    ps->p = (spair_t *)malloc((unsigned long)ps->sz * sizeof(spair_t));

    return ps;
}

static inline void check_enlarge_pairset(
        ps_t *ps,
        len_t added
        )
{
    if (ps->ld+added >= ps->sz) {
        ps->sz  = ps->sz*2 > ps->ld+added ? ps->sz*2 : ps->ld+added;
        ps->p   = realloc(ps->p, (unsigned long)ps->sz * sizeof(spair_t));
        memset(ps->p+ps->ld, 0,
                (unsigned long)(ps->sz-ps->ld) * sizeof(spair_t));
    }
}

static void free_pairset(
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

static void insert_and_update_spairs(
        ps_t *psl,
        bs_t *bs,
        ht_t *bht,
        ht_t *uht,
        stat_t *st
        )
{
    len_t i, j, k, l;

    spair_t *ps = psl->p;

    const len_t pl  = psl->ld;
    const len_t bl  = bs->ld;

    const hm_t nch = bs->hm[bl][3];

    reinitialize_hash_table(uht, bl);

    lms[bl] = hd[nch].sdm;

    for (i = bs->lo; i < bl; ++i) {
        if (red[i]) {
            continue;
        }
        if (check_monomial_division(nch, bs->hm[i][3], bht)) {
            /* printf("Mark polynomial %d unnecessary for new pairs\n", bload); */
            ps[pl].gen1 = i;
            ps[pl].gen2 = bl;
            ps[pl].lcm  = get_lcm(bs->hm[i][3], nch, bht, bht);
            bs->red[bl] = 1;
            st->num_redundant++;
            bs->ld++;
            psl->ld++;
            return;
        }
    }

    hl_t *plcm  = (hl_t *)malloc((unsigned long)(bl+1) * sizeof(hl_t));

    /* create all possible new pairs */
    for (i = 0, k = pl; i < bl; ++i, ++k) {
        ps[k].gen1  = i;
        ps[k].gen2  = bl;
        ps[k].lcm   = get_lcm(gbdt[i][3], nch, bht, uht);
        plcm[i]     = ps[k].lcm;
    }

    len_t nl  = k;
    /* Gebauer-Moeller: check old pairs first */
    /* note: old pairs are sorted by the given spair order */
    for (i = 0; i < pl; ++i) {
        j = ps[i].gen1;
        l = ps[i].gen2;
        if (check_monomial_division(ps[i].lcm, nch, bht)
                && hd[ps[i].lcm].val != hdu[plcm[j]].val
                && hd[ps[i].lcm].val != hdu[plcm[l]].val
           ) {
            ps[i].lcm = -1;
        }
    }
    /* check new pairs for redundancy */
    spair_t *pp = ps+pl;
    j = 0;
    for (i = 0; i < bl; ++i) {
        if (red[i] == 0) {
            pp[j++] = pp[i];
        }
    }
    /* sort new pairs by increasing lcm, earlier polys coming first */
    sort_r(pp, (unsigned long)j, sizeof(spair_t), spair_cmp, uht);
    for (i = 0; i < j; ++i) {
        plcm[i] = pp[i].lcm;
    }
    plcm[j]  = 0;
    const len_t pc  = j;

    j = 0;
    for (; j < pc; ++j) {
        if (plcm[j] < 0) {
            continue;
        }
        const hl_t plcmj = plcm[j];
        check_monomial_division_update(plcm, j, pc, plcmj);
    }

    /* remove useless pairs from pairset */
    j = 0;
    /* old pairs */
    for (i = 0; i < psl->ld; ++i) {
        if (ps[i].lcm < 0) {
            continue;
        }
        ps[j++] = ps[i];
    }
    if (bht->esz - bht->eld <= pc) {
        enlarge_basis_hash_table();
    }
    /* new pairs, wee need to add the lcm to the basis hash table */
    insert_in_basis_hash_table_plcms(psl, pp, j, pc, plcm);
    free(plcm);
    st->num_gb_crit +=  nl - psl->ld;

    /* mark redundant elements in basis */
    for (i = 0; i < bl; ++i) {
        if (red[i]) {
            continue;
        }
        if (check_monomial_division(gbdt[i][3], nch, bht)) {
            bs->red[i]  = 1;
            st->num_redundant++;
        }
    }
    bs->ld++;
}

static void update_basis(
        ps_t *ps,
        bs_t *bs,
        ht_t *bht,
        ht_t *uht,
        stat_t *st,
        const len_t npivs
        )
{
    len_t i;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* compute number of new pairs we need to handle at most */
    len_t np  = bload * npivs;
    for (i = 1; i < npivs; ++i) {
        np  = np + i;
    }
    check_enlarge_pairset(ps, np);

    for (i = 0; i < npivs; ++i) {
        insert_and_update_spairs(ps, bs, bht, uht, st);
    }

    bs->lo  = bs->ld;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->update_ctime  +=  ct1 - ct0;
    st->update_rtime  +=  rt1 - rt0;
}
