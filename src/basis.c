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

static basis_t *initialize_basis(
    const int32_t ngens
    )
{
  basis_t *bs = (basis_t *)malloc(sizeof(basis_t));;
  bs->ld  = 0;
  bs->sz  = 2*ngens;

  bs->pr  = (pr_t *)malloc((unsigned long)bs->sz * sizeof(pr_t));
  bs->lm  = (sdm_t *)malloc((unsigned long)bs->sz * sizeof(sdm_t));

  return bs;
}

static inline void check_enlarge_basis(
    basis_t *bs,
    const int32_t added
    )
{
  if (bs->ld + added >= bs->sz) {
    bs->sz  = bs->sz*2 > bs->ld+added ? bs->sz*2 : bs->ld+added;
    bs->pr  = realloc(bs->pr, (unsigned long)bs->sz * sizeof(pr_t));
    bs->lm  = realloc(bs->lm, (unsigned long)bs->sz * sizeof(sdm_t));
  }
}

static void free_basis(
    basis_t **bsp
    )
{
  int32_t i;
  basis_t *bs = *bsp;
  if (bs) {
    for (i = 0; i < bs->ld; ++i) {
      /* reset pointer of possible redundant elements */
      bs->pr[i].cf = 
        (val_t *)((long)bs->pr[i].cf & bmask);
      free(bs->pr[i].cf);
      free(bs->pr[i].mn);
    }
    free(bs->pr);
    free(bs->lm);
    free(bs);
    bs    = NULL;
    *bsp  = bs;
  }
}
