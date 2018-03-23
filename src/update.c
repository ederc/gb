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
    basis_t *bs,
    const pr_t *nelt
    )
{
  int32_t i, j, k, l;
  val_t *b;

  const len_t pl  = pload;
  const len_t bl  = bs->ld;

  reset_local_hash_table(bl);
  bs->pr+bl = nelt;
  /* i = bl * LM_LEN; */
  bs->lm[bl] = (ev+bs->pr[bl].mn[0])[HASH_SDM];

  /* printf("element added to basis %d (%d terms) at %p: ", bload, (nelt[0]-2)/2, nelt);
   * [> for (int32_t o = 2; o < 3; o += 2) { <]
   * for (int32_t o = 2; o < nelt[0]; o += 2) {
   *   printf("%d ", nelt[o+1]);
   *   for (int32_t p = 0; p < nvars; ++p) {
   *     printf("%d",(ev+nelt[o])[p]);
   *   }
   *   printf(" | ");
   * }
   * printf("\n"); */

#if INSERT_SMALL_FIRST
  for (i = bs->ol; i < bl; ++i) {
    if ((long)bs->pr[i].cf & bred) {
      continue;
    }
    if (check_monomial_division(ev+nelt[2], ev+bs->pr[i].mn[0])) {
      /* printf("Mark polynomial %d unnecessary for new pairs\n", bload); */
      ps[pl].gen1 = i;
      ps[pl].gen2 = bs->ld;
      ps[pl].lcm  = get_lcm(bs->pr[i].mn[0], nelt->mn[0]);
      ps[pl].deg  = (evl + ps[pl].lcm)[HASH_DEG];
      ps[pl].lcm  = insert_in_global_hash_table(evl+ps[pl].lcm);
      bs->pr[bl].cf = (val_t *)((long)bs->pr[bl].cf | bred);
      num_redundant++;
      bs->ld++;
      pload++;
      return;
    }
  }
#endif

  /* create all possible new pairs */
  for (i = 0, k = pl; i < bl; ++i, ++k) {
    ps[k].gen1  = i;
    ps[k].gen2  = bl;
    ps[k].lcm   = get_lcm(bs->pr[i].mn[0], nelt->mn[0]);

    b = (val_t *)((long)bs->pr[i].cf & bmask);
    if (b != bs->pr[i].cf) {
      ps[k].deg = -1; /* redundant pair */
    } else {
      if (lcm_equals_multiplication(bs->pr[i].mn[0], nelt->mn[0], ps[k].lcm)) {
        ps[k].deg = -2; /* criterion */
      } else {
        ps[k].deg = (evl + ps[k].lcm)[HASH_DEG];
      }
    }
  }
  
  const len_t nl  = k;
  /* Gebauer-Moeller: check old pairs first */
  /* note: old pairs are sorted by the given spair order */
  for (i = 0; i < pl; ++i) {
    j = ps[i].gen1;
    l = ps[i].gen2;
    /* if (ps[i].lcm != ps[pload+j].lcm */
    if (check_monomial_division(ev+ps[i].lcm, ev+nelt[2])
        && (ev+ps[i].lcm)[HASH_VAL] != (evl+ps[pl+j].lcm)[HASH_VAL]
        && (ev+ps[i].lcm)[HASH_VAL] != (evl+ps[pl+l].lcm)[HASH_VAL]
        ) {
        /* && ps[i].lcm != ps[pload+j].lcm) { */
      ps[i].deg = -1;
    }
  }

  /* sort new pairs by increasing lcm, earlier polys coming first */
  qsort(ps+pl, (unsigned long)bl, sizeof(spair_t), &spair_local_cmp);

  /* check with earlier new pairs */
  for (j = pl; j < nl; ++j) {
    l = j;
    if (ps[j].deg != -1) {
      i = j+1;
      while (i < nl && ps[i].lcm == ps[j].lcm)
        ++i;
      l = i-1;
      while (i < nl) {
        /* if (check_monomial_division(sp[i].lcm, sp[j].lcm, ht) != 0) { */
        if (ps[i].deg >= 0 &&
            check_monomial_division(evl+ps[i].lcm, evl+ps[j].lcm) != 0) {
          ps[i].deg  = -1;
        }
        ++i;
      }
    }
    j = l;
  }
  /* note that k = pload + bload from very first loop */
  for (i = pl; i < nl; ++i) {
    if (ps[i].deg == -1) {
      continue;
    }
    j = i+1;
    while (j < k && ps[j].lcm == ps[i].lcm) {
      ps[j++].deg = -1;
    }
    i = j-1;
  } 

  /* for (i = 0; i < pload; ++i) {
   *   printf("pair %d | g1 %d | g2 %d | red %d | ",
   *       i, ps[i].gen1, ps[i].gen2, ps[i].deg);
   *   for (j = 0; j < nvars; ++j) {
   *     printf("%d ", (ev+ps[i].lcm)[j]);
   *   }
   *   printf("\n");
   * }
   * for (; i < nl; ++i) {
   *   printf("pair %d | g1 %d | g2 %d | red %d | ",
   *       i, ps[i].gen1, ps[i].gen2, ps[i].deg);
   *   for (j = 0; j < nvars; ++j) {
   *     printf("%d ", (evl+ps[i].lcm)[j]);
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
    ps[i].lcm = insert_in_global_hash_table(evl+ps[i].lcm);
    ps[j++]   = ps[i];
  }
  num_gb_crit +=  nl - j;
  pload       =   j;

  /* mark redundant elements in basis */
  for (i = 0; i < bl; ++i) {
    if ((long)bs[i] & bred) {
      continue;
    }
    if (check_monomial_division(ev+bs->pr[i].mn[0], ev+nelt[2])) {
      /* printf("Mark polynomial %d unnecessary for new pairs\n", i); */
      bs->pr[i].cf  = (val_t *)((long)bs->pr[i].cf | bred);
      num_redundant++;
    }
  }
  bs->ld++;
  /* printf("lms:\n");
   * for (j = 0; j < bload; ++j) {
   *   for (i = 0; i < LM_LEN; ++i) {
   *     printf("%d ", lms[j*LM_LEN+i]);
   *   }
   * }
   * printf("\n"); */
}

static void update_basis(
    basis_t *bs,
    const mat_t *mat
    )
{
  int32_t i;

  /* timings */
  double ct0, ct1, rt0, rt1;
  ct0 = cputime();
  rt0 = realtime();

  check_enlarge_basis(bs, mat->np);

  /* compute number of new pairs we need to handle at most */
  len_t np  = bs->ld + mat->np;
  for (i = 1; i < mat->np; ++i) {
    np  = np + i;
  }
  check_enlarge_pairset(np);

#if INSERT_SMALL_FIRST
  for (i = 0; i < npivs; ++i) {
    insert_and_update_spairs(bs, mat[i]);
  }
#else
  for (i = 1; i <= npivs; ++i) {
    insert_and_update_spairs(bs, mat[nrows-i]);
  }
#endif
  bs->ol  = bs->ld;;

  /* timings */
  ct1 = cputime();
  rt1 = realtime();
  update_ctime  +=  ct1 - ct0;
  update_rtime  +=  rt1 - rt0;
}
