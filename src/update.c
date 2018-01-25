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

static void initialize_pairset(
    void
    )
{
  pload = 0;
  psize = 192;

  ps  = malloc((unsigned long)psize * sizeof(spair_t));
}

static inline void check_enlarge_pairset(
    int32_t added
    )
{
  if (pload+added >= psize) {
    psize = psize*2 > pload+added ? psize*2 : pload+added;
    ps    = realloc(ps, (unsigned long)psize * sizeof(spair_t));
    memset(ps+pload, 0, (unsigned long)(psize-pload) * sizeof(spair_t));
  }
}

static void free_pairset(
    void
    )
{
  if (ps) {
    free(ps);
    ps    = NULL;
    pload = 0;
    psize = 0;
  }
}

static void insert_and_update_spairs(
    val_t *nelt
    )
{
  int32_t i, j, k, l;
  val_t *b;

  check_enlarge_pairset(bload);

  /* move exponents to global hash table */
  for (i = 2; i < nelt[0]; i = i+2) {
    nelt[i] = insert_in_global_hash_table(evl + nelt[i]);
  }
  bs[bload] = nelt;
  /* printf("element added to basis: ");
   * for (int32_t o = 2; o < nelt[0]; o += 2) {
   *   printf("%d ", nelt[o+1]);
   *   for (int32_t p = 0; p < nvars; ++p) {
   *     printf("%d",(ev+nelt[o])[p]);
   *   }
   *   printf(" | ");
   * }
   * printf("\n"); */

  /* create all possible new pairs */
  for (i = 0, k = pload; i < bload; ++i, ++k) {
    b = (val_t *)((long)bs[i] & bmask);
    ps[k].gen1  = i;
    ps[k].gen2  = bload;
    ps[k].lcm   = get_lcm(b[2], nelt[2]);

    if (b != bs[i]) {
      ps[k].deg = -1; /* redundant pair */
    } else {
      if (lcm_equals_multiplication(b[2], nelt[2], ps[k].lcm)) {
        ps[k].deg = -2; /* criterion */
      } else {
        ps[k].deg = (ev + ps[k].lcm)[HASH_DEG];
      }
    }
  }
  
  /* Gebauer-Moeller: check old pairs first */
  /* note: old pairs are sorted by the given spair order */
  for (i = 0; i < pload; ++i) {
    j = ps[i].gen1;
    l = ps[i].gen2;
    /* if (ps[i].lcm != ps[pload+j].lcm */
    if (check_monomial_division(ps[i].lcm, nelt[2])
        && ps[i].lcm != ps[pload+l].lcm
        && ps[i].lcm != ps[pload+j].lcm) {
      ps[i].deg = -1;
    }
  }

  /* sort new pairs by increasing lcm, earlier polys coming first */
  qsort(ps+pload, (unsigned long)bload, sizeof(spair_t), &spair_cmp);

  /* check with earlier new pairs */
  for (j = pload; j < k; ++j) {
    l = j;
    if (ps[j].deg != -1) {
      i = j+1;
      while (i < k && ps[i].lcm == ps[j].lcm)
        ++i;
      l = i-1;
      while (i < k) {
        /* if (check_monomial_division(sp[i].lcm, sp[j].lcm, ht) != 0) { */
        if (ps[i].deg >= 0 &&
            check_monomial_division(ps[i].lcm, ps[j].lcm) != 0) {
          ps[i].deg  = -1;
        }
        ++i;
      }
    }
    j = l;
  }
  /* note that k = pload + bload from very first loop */
  for (i = pload; i < k; ++i) {
    if (ps[i].deg == -1) {
      continue;
    }
    j = i+1;
    while (j < k && ps[j].lcm == ps[i].lcm) {
      ps[j++].deg = -1;
    }
    i = j-1;
  } 

  /* remove useless pairs from pairset */
  j = 0;
  for (i = 0; i < k; ++i) {
    if (ps[i].deg < 0) {
      continue;
    }
    ps[j++] = ps[i];
  }
  num_gb_crit = k - j;
  pload       = j;

  /* mark redundant elements in basis */
  for (i = 0; i < bload; ++i) {
    if ((long)bs[i] & bred) {
      continue;
    }
    if (check_monomial_division(bs[i][2], nelt[2])) {
      bs[i] = (val_t *)((long)bs[i] | bred);
      num_redundant++;
    }
  }
  bload++;
}

static void update_basis(
    val_t **mat
    )
{
  int32_t i;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  check_enlarge_basis(npivs);

  for (i = 1; i <= npivs; ++i) {
    insert_and_update_spairs(mat[nrows-i]);
  }

  clear_local_hash_table();

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  update_ctime  +=  ct1 - ct0;
  update_rtime  +=  rt1 - rt0;
}
