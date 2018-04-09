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
 * \file basis.c
 * \brief General and global data
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

static bs_t *initialize_basis(
    int32_t ngens
    )
{
  bs_t *bs = (bs_t *)malloc(sizeof(bs_t));
  bs->ld  = 0;
  bs->sz  = 2*ngens;
  bs->ol  = bs->ld;

  bs->p   = (void **)malloc((unsigned long)bs->sz * sizeof(void *));
  bs->lm  = (len_t *)malloc((unsigned long)bs->sz * sizeof(len_t)); 

  return bs;
}

static inline void check_enlarge_basis(
    bs_t *bs,
    const int32_t added
    )
{
  if ((bs->ld + added) >= bs->sz) {
    bs->sz = (bs->sz * 2) > (bs->ld + added) ?
      (bs->sz * 2) : (bs->ld + added);
    bs->p   = realloc(bs->p, (unsigned long)bs->sz * sizeof(void *));
    bs->lm  = realloc(bs->lm, (unsigned long)bs->sz * sizeof(val_t));
  }
}

static void free_basis(
    bs_t **bsp
    )
{
  int32_t i;
  bs_t *bs = *bsp;

  if (bs) {
    for (i = 0; i < bs->ld; ++i) {
      /* reset pointer of possible redundant elements */
      bs->p[i]->cf = (void *)((long)bs->p[i]->cf & bmask);
      free(bs->p[i]->cf);
      free(bs->p[i]->ch);
    }
    free(bs->p);
    free(bs->lm);
    bs->lm  = NULL;
    free(bs);
    bs      = NULL;
  }
  *bsp  = bs;
}
