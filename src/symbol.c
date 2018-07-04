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
        dt_t **mat
        )
{
    len_t i, j, k, l, md, npairs;
    dt_t *b;
    deg_t d = 0;
    len_t load = 0, load_old = 0;
    hl_t lcm;
    len_t *gens;
    exp_t *elcm, *eb;

    /* timings */
    double ct0, ct1, rt0, rt1, rrt0, rrt1;
    ct0 = cputime();
    rt0 = realtime();

    /* sort pair set */
    rrt0 = realtime();
    qsort(ps, (unsigned long)pload, sizeof(spair_t), spair_cmp);
    rrt1 = realtime();
    pair_sort_rtime +=  rrt1 - rrt0;
    /* get minimal degree */
    md  = hd[ps[0].lcm].deg;

    /* select pairs of this degree respecting maximal selection size mnsel */
    for (i = 0; i < pload; ++i) {
        if (hd[ps[i].lcm].deg > md || i >= mnsel) {
            break;
        }
    }
    npairs  = i;
    /* if we stopped due to maximal selection size we still get the following
     * pairs of the same lcm in this matrix */
    if (i > mnsel) {
        j = i+1;
        while (ps[j].lcm == ps[i].lcm) {
            ++j;
        }
        npairs = j;
    }
    GB_DEBUG(SELDBG, " %6d/%6d pairs - deg %2d", npairs, pload, md);
    /* statistics */
    num_pairsred  +=  npairs;
    printf("npairs %d\n", npairs);

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
        load_old  = load;
        lcm   = ps[i].lcm;
        gens[load++] = ps[i].gen1;
        gens[load++] = ps[i].gen2;
        j = i+1;
        while (j < npairs && ps[j].lcm == lcm) {
            for (k = load_old; k < load; ++k) {
                if (ps[j].gen1 == gens[k]) {
                    break;
                }
            }
            if (k == load) {
                gens[load++]  = ps[j].gen1;
            }
            for (k = load_old; k < load; ++k) {
                if (ps[j].gen2 == gens[k]) {
                    break;
                }
            }
            if (k == load) {
                gens[load++]  = ps[j].gen2;
            }
            j++;
        }
        for (k = load_old; k < load; ++k) {
            /* ev might change when enlarging the hash table during insertion of a new
             * row in the matrix, thus we have to reset elcm inside the for loop */
            elcm  = ev[lcm];
            d = 0;
            b = gbdt[gens[k]];
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
            hd[mat[nrows][3]].idx = 2;
            nrows++;
        }

        i = j;
    }

    num_duplicates  +=  0;
    num_rowsred     +=  load;

    free(gens);

    /* remove selected spairs from pairset */
    for (j = npairs; j < pload; ++j) {
        ps[j-npairs]  = ps[j];
    }
    pload = pload - npairs;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    select_ctime  +=  ct1 - ct0;
    select_rtime  +=  rt1 - rt0;

    return mat;
}

static inline dt_t *find_multiplied_reducer(
        const dt_t m
        )
{
    len_t i, k;
    deg_t d = 0;
    dt_t *b;
    const exp_t * const e  = ev[m];
    exp_t *f;

    const len_t bl  = bload;
    const len_t os  = nvars & 1 ? 1 : 0;
    i = hd[m].div;

    const sdm_t ns  = ~hd[m].sdm;
start:
    while (i < bl-3) {
        if (lms[i] & ns &&
                lms[i+1] & ns &&
                lms[i+2] & ns &&
                lms[i+3] & ns) {
            num_sdm_found +=  4;
            i +=  4;
            continue;
        }
        while (lms[i] & ns) {
            i++;
        }
        b = gbdt[i];
        f = ev[b[3]];
        if ((e[0]-f[0]) < 0) {
            i++;
            goto start;
        }
        for (k = os; k < nvars; k += 2) {
            if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
                i++;
                goto start;
            }
        }
        for (k = 0; k < nvars; ++k) {
            etmp[k] = e[k] - f[k];
        }
        const hl_t h  = hd[m].val - hd[b[3]].val;
        for (k = 0; k < nvars; ++k) {
            d += etmp[k];
        }
        b = multiplied_polynomial_to_matrix_row(h, d, etmp, b);
        hd[m].div = i;
        return b;
    }
start2:
    while (i < bl) {
        if (lms[i] & ns) {
            num_sdm_found++;
            i++;
            continue;
        }
        b = gbdt[i];
        f = ev[b[3]];
        if ((e[0]-f[0]) < 0) {
            i++;
            goto start2;
        }
        for (k = os; k < nvars; k += 2) {
            if ((e[k]-f[k]) < 0 || (e[k+1]-f[k+1]) < 0) {
                i++;
                goto start2;
            }
        }
        for (k = 0; k < nvars; ++k) {
            etmp[k] = e[k] - f[k];
        }
        const hl_t h  = hd[m].val - hd[b[3]].val;
        for (k = 0; k < nvars; ++k) {
            d += etmp[k];
        }
        b = multiplied_polynomial_to_matrix_row(h, d, etmp, b);
        hd[m].div = i;
        return b;
    }
    hd[m].div = i;
    num_not_sdm_found++;
    return NULL;
}

static dt_t **symbolic_preprocessing(
        dt_t **mat
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
    for (i = 0; i < nrows; ++i) {
        const hl_t len = mat[i][2];
        /* check row reallocation only once per polynomial */
        if ((nrall - nrows) < (mat[i][2]-3)) {
            nrall = 2*nrall > mat[i][2] ? 2*nrall : (len_t)mat[i][2];
            mat   = realloc(mat, (unsigned long)nrall * sizeof(dt_t *));
        }
        for (j = 4; j < len; ++j) {
            m = mat[i][j];
            if (!hd[m].idx) {
                hd[m].idx = 1;
                ncols++;
                red = find_multiplied_reducer(m);
                if (red) {
                    hd[m].idx = 2;
                    /* add new reducer to matrix */
                    mat[nrows++]  = red;
                }
            }
        }
    }

    /* realloc to real size */
    mat   = realloc(mat, (unsigned long)nrows * sizeof(dt_t *));
    nrall = nrows;

    /* timings */
    ct1 = cputime();
    rt1 = realtime();
    symbol_ctime  +=  ct1 - ct0;
    symbol_rtime  +=  rt1 - rt0;

    return mat;
}
