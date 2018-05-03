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
    val_t *nelt
    )
{
  int32_t i, j, k, l;
  val_t *b;
  exp_t *ej;

  const len_t pl  = pload;
  const len_t bl  = bload;

  const int64_t hl  = HASH_LEN;

  reset_local_hash_table(bl);
  /* move exponents to global hash table */
  /* for (i = 2; i < nelt[0]; i = i+2) {
   *   nelt[i] = insert_in_global_hash_table(evl + nelt[i]);
   * } */
  bs[bl]  = nelt;
  /* i = bl * LM_LEN; */
  lms[bl] = (ev+bs[bl][2]*hl)[HASH_SDM];
  /* j = 0;
   * while (j < nvars) {
   *   lms[i++]  = (ev+bs[bl][2])[j];
   *   j++;
   * } */
  /* printf("element added to basis %d (%d terms) at %p: ", bload, (nelt[0]-2)/2, nelt);
   * [> for (int32_t o = 2; o < 3; o += 2) { <]
   * for (int32_t o = 2; o < nelt[0]; o += 2) {
   *   printf("%d ", nelt[o+1]);
   *   for (int32_t p = 0; p < nvars; ++p) {
   *     printf("%d",(ev+nelt[o]*hl)[p]);
   *   }
   *   printf(" | ");
   * }
   * printf("\n"); */

#if INSERT_SMALL_FIRST
  for (i = blold; i < bl; ++i) {
    if ((long)bs[i] & bred) {
      continue;
    }
    if (check_monomial_division(ev+nelt[2]*hl, ev+bs[i][2]*hl)) {
      /* printf("Mark polynomial %d unnecessary for new pairs\n", bload); */
      ps[pl].gen1 = i;
      ps[pl].gen2 = bl;
      ps[pl].lcm  = get_lcm(bs[i][2]*hl, nelt[2]*hl);
      ps[pl].deg  = (evl + ps[pl].lcm*hl)[HASH_DEG];
      ps[pl].lcm  = insert_in_global_hash_table(evl+ps[pl].lcm*hl);
      bs[bl]  = (val_t *)((long)bs[bl] | bred);
      num_redundant++;
      bload++;
      pload++;
      return;
    }
  }
#endif

  /* create all possible new pairs */
  for (i = 0, k = pl; i < bl; ++i, ++k) {
    b = (val_t *)((long)bs[i] & bmask);
    ps[k].gen1  = i;
    ps[k].gen2  = bl;
    ps[k].lcm   = get_lcm(b[2]*hl, nelt[2]*hl);

    if (b != bs[i]) {
      ps[k].deg = -1; /* redundant pair */
    } else {
      if (lcm_equals_multiplication(b[2]*hl, nelt[2]*hl, ps[k].lcm*hl)) {
        ps[k].deg = -2; /* criterion */
      } else {
        ps[k].deg = (evl + ps[k].lcm*hl)[HASH_DEG];
      }
    }
  }
  
  len_t nl  = k;
  /* Gebauer-Moeller: check old pairs first */
  /* note: old pairs are sorted by the given spair order */
  for (i = 0; i < pl; ++i) {
    j = ps[i].gen1;
    l = ps[i].gen2;
    /* if (ps[i].lcm != ps[pload+j].lcm */
    if (check_monomial_division(ev+ps[i].lcm*hl, ev+nelt[2]*hl)
        && (ev+ps[i].lcm*hl)[HASH_VAL] != (evl+ps[pl+j].lcm*hl)[HASH_VAL]
        && (ev+ps[i].lcm*hl)[HASH_VAL] != (evl+ps[pl+l].lcm*hl)[HASH_VAL]
        ) {
        /* && ps[i].lcm != ps[pload+j].lcm) { */
      ps[i].deg = -1;
    }
  }

  /* sort new pairs by increasing lcm, earlier polys coming first */
  double rrt0, rrt1;
  rrt0 = realtime();
  qsort(ps+pl, (unsigned long)bl, sizeof(spair_t), &spair_local_cmp);
  rrt1 = realtime();
  update1_rtime  +=  rrt1 - rrt0;

/*   for (i = pl; i < nl; ++i) {
 *     printf("%d || lcm %d | deg %d\n", i, ps[i].lcm, ps[i].deg);
 *   }
 *   printf("ctr %d\n", ctr);
 *   for (i = 0; i < ctr; ++i) {
 *     printf("lcm %d | fp %d\n", lcms[i], fp[i]);
 *   }
 *  */

  /* check with earlier new pairs */
  /* for (j = 0; j < ctr; ++j) {
   *   if (ps[fp[j]].deg != -1) {
   *     ej  = evl+ps[fp[j]].lcm;
   *     for (i = j+1; i < ctr; ++i) {
   *       if ([> ps[fp[i]].deg >= 0 && <]
   *           check_monomial_division(evl+ps[fp[i]].lcm, ej) != 0) {
   *         for (k = fp[i]; k < fp[i+1]; ++k) {
   *           if (ps[k].deg >= 0) {
   *             ps[k].deg = -1;
   *           }
   *         }
   *       }
   *     }
   *   }
   * } */

  for (j = pl; j < nl; ++j) {
    if (ps[j].deg < 0) {
      continue;
    }
    ej  = evl+ps[j].lcm*hl;
    l = j;
    i = j+1;
    while (i < nl && ps[i].lcm == ps[j].lcm) {
      ++i;
    }
    l = i-1;
    while (i < nl) {
      if (ps[i].deg >= 0 &&
          check_monomial_division(evl+ps[i].lcm*hl, ej) != 0) {
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

  /* for (i = 0; i < pload; ++i) {
   *   printf("pair[%d] %d | %d | %d | %d || ", i, ps[i].gen1, ps[i].gen2, ps[i].deg,
   *       ps[i].lcm);
   *   for (j = 0; j < nvars; ++j) {
   *     printf("%d ", (ev+ps[i].lcm*hl)[j]);
   *   }
   *   printf("\n");
   * }
   * for (; i < nl; ++i) {
   *   printf("pair[%d] %d | %d | %d | %d || ", i, ps[i].gen1, ps[i].gen2, ps[i].deg,
   *       ps[i].lcm);
   *   for (j = 0; j < nvars; ++j) {
   *     printf("%d ", (evl+ps[i].lcm*hl)[j]);
   *   }
   *   printf("\n");
   * } */

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
    ps[i].lcm = insert_in_global_hash_table(evl+ps[i].lcm*hl);
    ps[j++]   = ps[i];
  }
  num_gb_crit +=  nl - j;
  pload       =   j;

  /* mark redundant elements in basis */
  for (i = 0; i < bl; ++i) {
    if ((long)bs[i] & bred) {
      continue;
    }
    if (check_monomial_division(ev+bs[i][2]*hl, ev+nelt[2]*hl)) {
      /* printf("Mark polynomial %d unnecessary for new pairs\n", i); */
      bs[i] = (val_t *)((long)bs[i] | bred);
      num_redundant++;
    }
  }
  bload++;
  /* printf("lms:\n");
   * for (j = 0; j < bload; ++j) {
   *   for (i = 0; i < LM_LEN; ++i) {
   *     printf("%d ", lms[j*LM_LEN+i]);
   *   }
   * }
   * printf("\n"); */
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

  /* compute number of new pairs we need to handle at most */
  len_t np  = bload * npivs;
  for (i = 1; i < npivs; ++i) {
    np  = np + i;
  }
  check_enlarge_pairset(np);

#if INSERT_SMALL_FIRST
  for (i = 0; i < npivs; ++i) {
    insert_and_update_spairs(mat[i]);
  }
#else
  for (i = 1; i <= npivs; ++i) {
    insert_and_update_spairs(mat[nrows-i]);
  }
#endif
  blold = bload;

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  update_ctime  +=  ct1 - ct0;
  update_rtime  +=  rt1 - rt0;
}
