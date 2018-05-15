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

#define INSERT_SMALL_FIRST 1

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
    bs_t *bs,
    row_t *row,
    md_t *md
    )
{
  int32_t i, j, k, l;
  exp_t *ej;
  row_t *p;

  const len_t pl  = pload;
  const len_t bl  = bs->ld;

  reset_local_hash_table(bl);
  bs->p[bl]     = row;
  bs->lm[bl]    = (ev+row->ch[0])[HASH_SDM];
  bs->p[bl]->rd = 0;
  /* printf("new element added to basis: (%d) ", bl);
   * for (j = 0; j < row->sz; ++j) {
   *   for (i = 0; i < nvars; ++i) {
   *     printf("%d ", (ev+row->ch[j])[i]);
   *   }
   *   printf("  +  ");
   * }
   * printf("\n"); */

#if INSERT_SMALL_FIRST
  for (i = bs->ol; i < bl; ++i) {
    p = bs->p[i];
    if (p->rd) {
      continue;
    }
    if (check_monomial_division(ev+row->ch[0], ev+p->ch[0])) {
      ps[pl].gen1 = i;
      ps[pl].gen2 = bl;
      ps[pl].lcm  = get_lcm(p->ch[0], row->ch[0]);
      ps[pl].deg  = (evl + ps[pl].lcm)[HASH_DEG];
      ps[pl].lcm  = insert_in_global_hash_table(evl+ps[pl].lcm);
      row->rd     = 1;
      md->num_redundant++;
      bs->ld++;
      pload++;
      return;
    }
  }
#endif

  /* create all possible new pairs */
  for (i = 0, k = pl; i < bl; ++i, ++k) {
    p = bs->p[i];
    ps[k].gen1  = i;
    ps[k].gen2  = bl;
    ps[k].lcm   = get_lcm(p->ch[0], row->ch[0]);

    if (p->rd) {
      ps[k].deg = -1; /* redundant pair */
    } else {
      if (lcm_equals_multiplication(p->ch[0], row->ch[0], ps[k].lcm)) {
        ps[k].deg = -2; /* criterion */
      } else {
        ps[k].deg = (evl + ps[k].lcm)[HASH_DEG];
      }
    }
  }
  
  len_t nl  = k;
  /* Gebauer-Moeller: check old pairs first */
  /* note: old pairs are sorted by the given spair order */
  for (i = 0; i < pl; ++i) {
    j = ps[i].gen1;
    l = ps[i].gen2;
    if (check_monomial_division(ev+ps[i].lcm, ev+row->ch[0])
        && (ev+ps[i].lcm)[HASH_VAL] != (evl+ps[pl+j].lcm)[HASH_VAL]
        && (ev+ps[i].lcm)[HASH_VAL] != (evl+ps[pl+l].lcm)[HASH_VAL]
        ) {
      ps[i].deg = -1;
    }
  }

  /* sort new pairs by increasing lcm, earlier polys coming first */
  qsort(ps+pl, (unsigned long)bl, sizeof(spair_t), &spair_local_cmp);


  for (j = pl; j < nl; ++j) {
    if (ps[j].deg < 0) {
      continue;
    }
    ej  = evl+ps[j].lcm;
    l = j;
    i = j+1;
    while (i < nl && ps[i].lcm == ps[j].lcm) {
      ++i;
    }
    l = i-1;
    while (i < nl) {
      if (ps[i].deg >= 0 &&
          check_monomial_division(evl+ps[i].lcm, ej) != 0) {
        ps[i].deg  = -1;
      }
      ++i;
    }
    j = l;
  }

  /* remove deg == -1 pairs from list */
  j = pl;
  for (i = pl; i < nl; ++i) {
    if (ps[i].deg == -1) {
      continue;
    }
    ps[j++] = ps[i];
  }
  nl = j;

  for (i = pl; i < nl; ++i) {
    j = i+1;
    while (j < nl && ps[j].lcm == ps[i].lcm) {
      ps[j++].deg = -1;
    }
    i = j-1;
  } 

  /* remove useless pairs from pairset */
  j = 0;
  /* old pairs */
  for (i = 0; i < pload; ++i) {
    if (ps[i].deg < 0) {
      continue;
    }
    ps[j++] = ps[i];
  }
  /* new pairs, wee need to add the lcm to the global hash table */
  for (; i < nl; ++i) {
    if (ps[i].deg < 0) {
      continue;
    }
    ps[i].lcm = insert_in_global_hash_table(evl+ps[i].lcm);
    ps[j++]   = ps[i];
  }
  md->num_gb_crit +=  nl - j;
  pload       =   j;

  /* mark redundant elements in basis */
  for (i = 0; i < bl; ++i) {
    p = bs->p[i];
    if (p->rd) {
      continue;
    }
    if (check_monomial_division(ev+p->ch[0], ev+row->ch[0])) {
      p->rd = 1;
      md->num_redundant++;
    }
  }
  bs->ld++;
}

static void update_basis(
    bs_t *bs,
    md_t * md,
    const mat_t * const mat
    )
{
  int32_t i;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  check_enlarge_basis(bs, mat->np);

  /* compute number of new pairs we need to handle at most */
  len_t np  = bs->ld * mat->np;
  for (i = mat->np-1; i > 0; --i) {
    np  = np + i;
  }
  check_enlarge_pairset(np);

#if INSERT_SMALL_FIRST
  for (i = 0; i < mat->np; ++i) {
    insert_and_update_spairs(bs, mat->r[i], md);
  }
#else
  for (i = 1; i <= mat->np; ++i) {
    insert_and_update_spairs(bs, mat->r[mat->nr-i], md);
  }
#endif
  bs->ol  = bs->ld;

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  md->update_ctime  +=  ct1 - ct0;
  md->update_rtime  +=  rt1 - rt0;
}
