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
 * \file symbol.c
 * \brief Symbolic preprocessing routines
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

static void select_spairs_by_minimal_degree(
        mat_t *mat,
        const bs_t * const bs,
        ps_t *psl,
        stat_t *st,
        ht_t *sht,
        ht_t *bht
        )
{
    len_t i, j, k, l, md, nps, npd;
    hm_t *b;
    deg_t d = 0;
    len_t load = 0;
    hl_t lcm;
    len_t *gens;
    exp_t *elcm, *eb;
    exp_t *etmp = bht->ev[0];
    hm_t **rows = mat->r;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    spair_t *ps     = psl->p;
    const len_t nv  = bht->nv;

    /* sort pair set */
    sort_r(ps, (unsigned long)psl->ld, sizeof(spair_t), spair_degree_cmp, bht);
    /* get minimal degree */
    md  = bht->hd[ps[0].lcm].deg;

    /* select pairs of this degree respecting maximal selection size mnsel */
    for (i = 0; i < psl->ld; ++i) {
        if (bht->hd[ps[i].lcm].deg > md) {
            break;
        }
    }
    npd  = i;
    sort_r(ps, (unsigned long)npd, sizeof(spair_t), spair_cmp, bht);
    /* now do maximal selection if it applies */
    
    /* if we stopped due to maximal selection size we still get the following
     * pairs of the same lcm in this matrix */
    if (npd > psl->mnsel) {
        nps = psl->mnsel;
        lcm = ps[nps].lcm;
        while (nps < npd && ps[nps+1].lcm == lcm) {
            nps++;
        }
    } else {
        nps = npd;
    }
    if (st->info_level > 1) {
        printf("%3d  %6d %7d", md, nps, psl->ld);
        fflush(stdout);
    }
    /* statistics */
    st->num_pairsred  +=  nps;

    /* list for generators */
    gens  = (len_t *)malloc(2 * (unsigned long)nps * sizeof(len_t));
    /* preset matrix meta data */
    rows    = (hm_t **)malloc(2 * (unsigned long)nps * sizeof(hm_t *));
    mat->sz = 2 * nps;
    mat->nc = mat->ncl = mat->ncr = 0;
    mat->nr = 0;

    i = 0;
    while (i < nps) {
        /* ncols initially counts number of different lcms */
        mat->nc++;
        load  = 0;
        lcm   = ps[i].lcm;
        j = i;
        while (j < nps && ps[j].lcm == lcm) {
            gens[load++] = ps[j].gen1;
            gens[load++] = ps[j].gen2;
            ++j;
        }
        /* sort gens set */
        qsort(gens, (unsigned long)load, sizeof(len_t), gens_cmp);

        len_t prev  = -1;
        for (k = 0; k < load; ++k) {
            /* check sorted list for doubles */
            if (gens[k] ==  prev) {
                continue;
            }
            prev  = gens[k];
            /* ev might change when enlarging the hash table during insertion of a new
             * row in the matrix, thus we have to reset elcm inside the for loop */
            elcm  = bht->ev[lcm];
            d     = 0;
            b     = bs->hm[prev];
            eb    = bht->ev[b[3]];
            for (l = 0; l < nv; ++l) {
                etmp[l] = (exp_t)(elcm[l] - eb[l]);
                d     +=  etmp[l];
            }
            const hl_t h  = bht->hd[lcm].val - bht->hd[b[3]].val;
            rows[mat->nr] = multiplied_poly_to_matrix_row(sht, bht, h, d, etmp, b);
            /* mark lcm column as lead term column */
            sht->hd[rows[mat->nr][3]].idx = 2;
            mat->nr++;
        }

        i = j;
    }
    st->num_rowsred +=  mat->nr - mat->nc;
    st->current_deg =   md;

    free(gens);

    /* remove selected spairs from pairset */
    memmove(ps, ps+nps, (unsigned long)(psl->ld-nps) * sizeof(spair_t));
    psl->ld -=  nps;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->select_ctime  +=  ct1 - ct0;
    st->select_rtime  +=  rt1 - rt0;
}

static inline hm_t *find_multiplied_reducer(
        const bs_t * const bs,
        const hm_t m,
        const ht_t * const bht,
        ht_t *sht
        )
{
    len_t i, k;
    const exp_t * const e  = sht->ev[m];

    const hd_t hdm  = sht->hd[m];
    const len_t bl  = bs->ld;
    const sdm_t ns  = ~hdm.sdm;
    const deg_t hdd = hdm.deg;

    const len_t nv  = bht->nv;

    const sdm_t * const lms = bs->lm;

    exp_t *etmp = bht->ev[0];
    const hd_t * const hdb  = bht->hd;
    exp_t * const * const evb = bht->ev;

    i = 0;
start:
    while (i < bl && lms[i] & ns) {
        i++;
    }
    if (i < bl) {
        const hm_t *b = bs->hm[i];
        const deg_t d = hdd - hdb[b[3]].deg;
        if (d < 0) {
            i++;
            goto start;
        }
        const exp_t * const f = evb[b[3]];
        for (k=nv-1; k >= 0; --k) {
            etmp[k] = (exp_t)(e[k]-f[k]);
            if (etmp[k] < 0) {
                i++;
                goto start;
            }
        }
        const hl_t h  = hdm.val - hdb[b[3]].val;
        return multiplied_poly_to_matrix_row(sht, bht, h, d, etmp, b);
    } else {
        return NULL;
    }
}

static void symbolic_preprocessing(
        mat_t *mat,
        const bs_t * const bs,
        stat_t *st,
        ht_t *sht,
        const ht_t * const bht
        )
{
    len_t i;
    hm_t *red;
    hm_t **rows = mat->r;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* note that we have already counted the different lcms, i.e.
     * ncols until this step. moreover, we have also already marked
     * the corresponding hash indices to represent lead terms. so
     * we only have to do the bookkeeping for newly added reducers
     * in the following. */

    const len_t oesld = sht->eld;
    hd_t *hds = sht->hd; 
    i = 1;
    /* we only have to check if idx is set for the elements already set
     * when selecting spairs, afterwards (second for loop) we do not
     * have to do this check */
    while (mat->sz <= mat->nr + oesld) {
        mat->sz *=  2;
        rows    =   realloc(rows, (unsigned long)mat->sz * sizeof(hm_t *));
    }
    for (; i < oesld; ++i) {
        if (!hds[i].idx) {
            hds[i].idx = 1;
            mat->nc++;
            red = find_multiplied_reducer(bs, i, bht, sht);
            if (red) {
                hds[i].idx = 2;
                /* add new reducer to matrix */
                rows[mat->nr++] = red;
            }
        }
    }
    for (; i < sht->eld; ++i) {
        if (mat->sz == mat->nr) {
            mat->sz *=  2;
            rows    =   realloc(rows, (unsigned long)mat->sz * sizeof(hm_t *));
        }
        hds[i].idx = 1;
        mat->nc++;
        red = find_multiplied_reducer(bs, i, bht, sht);
        if (red) {
            hds[i].idx = 2;
            /* add new reducer to matrix */
            rows[mat->nr++]  = red;
        }
    }
    /* realloc to real size */
    rows    = realloc(rows, (unsigned long)mat->nr * sizeof(hm_t *));
    mat->sz = mat->nr;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->symbol_ctime  +=  ct1 - ct0;
    st->symbol_rtime  +=  rt1 - rt0;
}
