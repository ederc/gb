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

static dt_t **select_spairs_by_minimal_degree(
        ps_t *psl,
        dt_t **mat,
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

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    spair_t *ps = psl->p;

    /* sort pair set */
    qsort(ps, (unsigned long)psl->ld, sizeof(spair_t), spair_cmp);
    /* get minimal degree */
    md  = hd[ps[0].lcm].deg;

    /* select pairs of this degree respecting maximal selection size mnsel */
    for (i = 0; i < psl->ld; ++i) {
        if (hd[ps[i].lcm].deg > md || i >= psl->mnsel) {
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
    mat   = (dt_t **)malloc(2 * (unsigned long)npairs * sizeof(dt_t *));
    nrall = 2 * npairs;
    ncols = ncl = ncr = 0;
    nrows = 0;

    i = 0;
    while (i < npairs) {
        /* ncols initially counts number of different lcms */
        ncols++;
        load_all  += load;
        load  = 0;
        lcm   = ps[i].lcm;
        j = i;
        while (j < npairs && ps[j].lcm == lcm) {
            gens[load++] = ps[i].gen1;
            gens[load++] = ps[i].gen2;
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
            elcm  = ev[lcm];
            d = 0;
            b = gbdt[prev];
            /* b = (hl_t *)((long)bs[gens[k]] & bmask); */
            eb  = ev[b[3]];
            /* m = monomial_division_no_check(lcm, b[2]); */
            for (l = 0; l < nvars; ++l) {
                etmp[l] = elcm[l] - eb[l];
                d     +=  etmp[l];
            }
            const hl_t h  = hd[lcm].val - hd[b[3]].val;
            mat[nrows]    = multiplied_polynomial_to_matrix_row(h, d, etmp, b);
            /* mark lcm column as lead term column */
            /* for (int ii=3; ii < mat[nrows][2]; ++ii) {
             *     printf("%d | ", mat[nrows][ii]);
             *     for (int jj=0; jj < nvars; ++jj) {
             *         printf("%d ", evs[mat[nrows][ii]][jj]);
             *     }
             *     printf("\n");
             * }
             * printf("\n");
             * printf("%d || idx %d\n", nrows, hds[mat[nrows][3]].idx); */
            hds[mat[nrows][3]].idx = 2;
            nrows++;
        }

        i = j;
    }

    /* load_all are all rows from the chosen spairs. for each column, i.e.
     * each lcm we have one reducer, the other rows have to be reduced.
     * thus in the end we have to subtract ncols from the number of rows
     * considered. */
    st->num_rowsred +=  load_all - ncols;
    st->current_deg =   md;

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

static inline dt_t *find_multiplied_reducer(
        const dt_t m
        )
{
    len_t i, k;
    deg_t d = 0;
    const exp_t * const e  = evs[m];
    exp_t *f;

    const len_t bl  = bload;
    const len_t nv  = nvars;
    const len_t os  = nv & 1 ? 1 : 0;
    const sdm_t ns  = ~hds[m].sdm;

    i = 0;
start:
    while (i < bl-3) {
        if ((lms[i] & ns) &&
            (lms[i+1] & ns) &&
            (lms[i+2] & ns) &&
            (lms[i+3] & ns)) {
            i +=  4;
            continue;
        }
        while (lms[i] & ns) {
            i++;
        }
        const dt_t *b = gbdt[i];
        f = ev[b[3]];
        if ((e[0]-f[0]) < 0) {
            i++;
            goto start;
        }
        for (k = os; k < nv; k += 2) {
            if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
                i++;
                goto start;
            }
        }
        for (k = 0; k < nv; ++k) {
            etmp[k] = e[k] - f[k];
        }
        const hl_t h  = hds[m].val - hd[b[3]].val;
        for (k = 0; k < nv; ++k) {
            d += etmp[k];
        }
        return multiplied_polynomial_to_matrix_row(h, d, etmp, b);
    }
start2:
    while (i < bl) {
        if (lms[i] & ns) {
            i++;
            continue;
        }
        const dt_t *b = gbdt[i];
        f = ev[b[3]];
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
        const hl_t h  = hds[m].val - hd[b[3]].val;
        for (k = 0; k < nv; ++k) {
            d += etmp[k];
        }
        return multiplied_polynomial_to_matrix_row(h, d, etmp, b);
    }
    return NULL;
}

static dt_t **symbolic_preprocessing(
        dt_t **mat,
        stat_t *st
        )
{
    len_t i;
    dt_t *red;

    /* timings */
    double ct0, ct1, rt0, rt1;
    ct0 = cputime();
    rt0 = realtime();

    /* note that we have already counted the different lcms, i.e.
     * ncols until this step. moreover, we have also already marked
     * the corresponding hash indices to represent lead terms. so
     * we only have to do the bookkeeping for newly added reducers
     * in the following. */

    for (i = 1; i < esld; ++i) {
        if (esld >= essz) {
            enlarge_symbolic_hash_table();
        }
        if (nrall == nrows) {
            nrall *=  2;
            mat   =   realloc(mat, (unsigned long)nrall * sizeof(dt_t *));
        }
        if (!hds[i].idx) {
            hds[i].idx = 1;
            ncols++;
            red = find_multiplied_reducer(i);
            if (red) {
                hds[i].idx = 2;
                /* add new reducer to matrix */
                mat[nrows++]  = red;
            }
        }
    }
    /* realloc to real size */
    mat   = realloc(mat, (unsigned long)nrows * sizeof(dt_t *));
    nrall = nrows;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    st->symbol_ctime  +=  ct1 - ct0;
    st->symbol_rtime  +=  rt1 - rt0;

    return mat;
}
