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

#include "symbol.h"

mat_t *select_spairs_by_minimal_degree(
        mat_t *mat,
        ps_t *psl,
        ht_t *ht,
        const bs_t * const bs,
        stat_t *st
        )
{
    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    len_t i, j, k, l, md, npairs;
    mon_t b;
    hd_t *lcm;
    len_t *gens;
    exp_t *elcm, *eb;
    deg_t d = 0;
    len_t load = 0, load_all = 0;
    exp_t *etmp = ht->hd[0].exp;
    spair_t *ps = psl->p;

    const len_t nv  = gbnv;

    /* sort pair set */
    qsort(ps, (unsigned long)psl->ld, sizeof(spair_t), &spair_cmp);
    /* get minimal degree */
    md  = ps[0].lcm->deg;

    /* select pairs of this degree respecting maximal selection size mnsel */
    for (i = 0; i < psl->ld; ++i) {
        if (ps[i].lcm->deg > md || i >= psl->mnsel) {
            break;
        }
    }
    npairs  = i;
    /* if we stopped due to maximal selection size we still get the following
     * pairs of the same lcm in this matrix */
    if (i > psl->mnsel) {
        j = i+1;
        while (ps[j].lcm->val== ps[i].lcm->val) {
            ++j;
        }
        npairs = j;
    }
    if (st->info_level > 1) {
        printf("%3d  %6d %7d", md, npairs, psl->ld);
        fflush(stdout);
    }
    /* statistics */
    st->num_pairsred  +=  npairs;
    /* printf("npairs %d\n", npairs); */

    gens  = (bi_t *)malloc(2 * (unsigned long)npairs * sizeof(bi_t));

    /* preset matrix meta data */
    mat->mp = realloc(mat->mp, 2 * (unsigned long)npairs * sizeof(mon_t));
    mat->na = 2 * npairs;
    mat->nc = mat->ncl = mat->ncr = 0;
    mat->nr = mat->nru = mat->nrl = 0;

    i = 0;
    while (i < npairs) {
        /* ncols initially counts number of different lcms */
        mat->nc++;
        load_all  += load;
        load  = 0;
        lcm   = ps[i].lcm;
        gens[load++] = ps[i].gen1;
        gens[load++] = ps[i].gen2;
        j = i+1;
        while (j < npairs && ps[j].lcm == lcm) {
            for (k = 0; k < load; ++k) {
                if (ps[j].gen1 == gens[k]) {
                    break;
                }
            }
            if (k == load) {
                gens[load++]  = ps[j].gen1;
            }
            for (k = 0; k < load; ++k) {
                if (ps[j].gen2 == gens[k]) {
                    break;
                }
            }
            if (k == load) {
            k    gens[load++]  = ps[j].gen2;
            }
            j++;
        }
        for (k = 0; k < load; ++k) {
            /* ev might change when enlarging the hash table during
             * insertion of a new row in the matrix, thus we have to
             * reset elcm inside the for loop */
            elcm  = lcm->exp;
            d = 0;
            b   = bs->m[gens[k]];
            eb  = b.h[0]->exp;
            for (l = 0; l < nv; ++l) {
                etmp[l] =   elcm[l] - eb[l];
                d       +=  etmp[l];
            }
            const hl_t h    = lcm->val - b.h[0]->val;
            multiplied_polynomial_to_matrix_row(mat, ht, h, d, etmp, b);
            /* mark lcm column as lead term column */
            x = 2;
            mat->nr++;
        }
        i = j;
    }

    /* load_all are all rows from the chosen spairs. for each column, i.e.
     * each lcm we have one reducer, the other rows have to be reduced.
     * thus in the end we have to subtract ncols from the number of rows
     * considered. */
    st->num_rowsred     +=  load_all - mat->nc;

    free(gens);

    /* remove selected spairs from pairset */
    memmove(ps, ps+npairs, (unsigned long)(psl->ld-npairs) * sizeof(spair_t));
    psl->ld -=  npairs;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->select_ctime  +=  ct1 - ct0;
    st->select_rtime  +=  rt1 - rt0;

    return mat;
}

mat_t *symbolic_preprocessing(
        mat_t *mat,
        ht_t *ht,
        const bs_t *const bs,
        stat_t *st
        )
{
    len_t i, j;
    int red;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* note that we have already counted the different lcms, i.e.
     * ncols until this step. moreover, we have also already marked
     * the corresponding hash indices to represent lead terms. so
     * we only have to do the bookkeeping for newly added reducers
     * in the following. */

    for (i = 1; i < ht->eld; ++i) {
        if (ht->hd[i].idx == 0) {
            ht->hd[i].idx  = 1;
            mat->nc++;
            red = find_multiplied_reducer(mat, ht, &(ht->hd[i]), bs);
            if (red == 1) {
                ht->hd[i].idx  = 2;
            }
        }
    }

    /* realloc to real size */
    mat->mp = realloc(mat->mp, (unsigned long)mat->nr * sizeof(mon_t));
    mat->na = mat->nr;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->symbol_ctime  +=  ct1 - ct0;
    st->symbol_rtime  +=  rt1 - rt0;

    return mat;
}

int find_multiplied_reducer(
        mat_t *mat,
        ht_t *ht,
        hd_t *m,
        const bs_t *const bs
        )
{
    len_t i, k;
    deg_t d = 0;
    mon_t b;
    exp_t *f;

    const exp_t * const e   = m->exp;
    const sdm_t *const lms  = bs->lm;

    const len_t nv  = gbnv;
    const len_t bl  = bs->ld;
    const len_t os  = nv & 1 ? 1 : 0;
    const sdm_t ns  = ~(m->sdm);

    exp_t *etmp = ht->hd[0].exp;

    i = 0;
start1:
    while (i < bl-3) {
        if (lms[i] & ns &&
                lms[i+1] & ns &&
                lms[i+2] & ns &&
                lms[i+3] & ns) {
            i +=  4;
            continue;
        }
        while (lms[i] & ns) {
            i++;
        }
        b = bs->m[i];
        f = b.h[0]->exp;
        if ((e[0]-f[0]) < 0) {
            i++;
            goto start1;
        }
        for (k = os; k < nv; k += 2) {
            if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
                i++;
                goto start1;
            }
        }
        for (k = 0; k < nv; ++k) {
            etmp[k] = e[k] - f[k];
        }
        const val_t h = m->val - b.h[0]->val;
        for (k = 0; k < nv; ++k) {
            d += etmp[k];
        }
        multiplied_polynomial_to_matrix_row(h, d, etmp, b, ht);
        return 1;
    }
start2:
    while (i < bl) {
        if (lms[i] & ns) {
            i++;
            continue;
        }
        b = bs->m[i];
        f = b.h[0]->exp;
        if ((e[0]-f[0]) < 0) {
            i++;
            goto start2;
        }
        for (k = os; k < nv; k += 2) {
            if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
                i++;
                goto start2;
            }
        }
        for (k = 0; k < nv; ++k) {
            etmp[k] = e[k] - f[k];
        }
        const hl_t h  = m->val - b.h[0]->val;
        for (k = 0; k < nv; ++k) {
            d += etmp[k];
        }
        multiplied_polynomial_to_matrix_row(h, d, etmp, b, ht);
        return 1;
    }
    return 0;
}
