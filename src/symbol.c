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
        stat_t *st
        )
{
    len_t i, j, k, l, md, npairs;
    dt_t *b;
    deg_t d = 0;
    len_t load = 0, load_all = 0;
    hl_t lcm;
    len_t *gens;
    exp_t *elcm, *eb;

    exp_t etmp[ht->nv];

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    spair_t *ps = psl->p;

    /* sort pair set */
    qsort(ps, (unsigned long)psl->ld, sizeof(spair_t), spair_cmp);
    /* get minimal degree */
    md  = ht->hd[ps[0].lcm].deg;

    /* select pairs of this degree respecting maximal selection size mnsel */
    for (i = 0; i < psl->ld; ++i) {
        if (ht->hd[ps[i].lcm].deg > md || i >= psl->mnsel) {
            break;
        }
    }
    npairs  = i;
    /* if we stopped due to maximal selection size we still get the following
     * pairs of the same lcm in this matrix */
    if (i > psl->mnsel) {
        j = i+1;
        while (ps[j].lcm == ps[i].lcm) {
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

    gens  = (len_t *)malloc(2 * (unsigned long)npairs * sizeof(len_t));

    /* preset matrix meta data */
    mat->r  = realloc(mat->r, 2 * (unsigned long)npairs * sizeof(hl_t *));
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
                gens[load++]  = ps[j].gen2;
            }
            j++;
        }
        for (k = 0; k < load; ++k) {
            /* ev might change when enlarging the hash table during
             * insertion of a new row in the matrix, thus we have to
             * reset elcm inside the for loop */
            elcm  = ht->ev[lcm];
            d = 0;
            b   = ht->hd[gens[k]];
            eb  = ht->ev[b[3]];
            for (l = 0; l < ht->nv; ++l) {
                etmp[l] =   elcm[l] - eb[l];
                d       +=  etmp[l];
            }
            const hl_t h    = ht->hd[lcm].val - ht->hd[b[3]].val;
            mat->r[mat->nr] = multiplied_polynomial_to_matrix_row(
                    h, d, etmp, b, ht);
            /* mark lcm column as lead term column */
            ht->hd[mat->r[mat->nr][3]].idx = 2;
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
    dt_t *red;
    dt_t m;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* note that we have already counted the different lcms, i.e.
     * ncols until this step. moreover, we have also already marked
     * the corresponding hash indices to represent lead terms. so
     * we only have to do the bookkeeping for newly added reducers
     * in the following. */

    /* get reducers from basis */
    for (i = 0; i < mat->nr; ++i) {
        /* printf("%p | %d / %d (i / nrows)\n",mat[i], i, nrows); */
        const hl_t len = mat->r[i][2];
        /* check row reallocation only once per polynomial */
        if ((mat->na - mat->nr) < (len-3)) {
            mat->na +=  mat->na > len ? mat->na : (len_t)len;
            mat-r   =   realloc(mat-r,
                    (unsigned long)mat->na * sizeof(dt_t *));
        }
        for (j = 4; j < len; ++j) {
            m = mat->r[i][j];
            /* printf("hd[%d].idx = %d\n", m, hd[m].idx); */
            if (!ht->hd[m].idx) {
                ht->hd[m].idx = 1;
                mat->nc++;
                red = find_multiplied_reducer(m, ht, bs);
                if (red) {
                    ht->hd[m].idx = 2;
                    /* add new reducer to matrix */
                    mat->r[mat->nr++]  = red;
                }
            }
        }
    }

    /* realloc to real size */
    mat-r   = realloc(mat->r, (unsigned long)mat->nr * sizeof(dt_t *));
    mat->na = mat->nr;;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->symbol_ctime  +=  ct1 - ct0;
    st->symbol_rtime  +=  rt1 - rt0;

    return mat;
}
