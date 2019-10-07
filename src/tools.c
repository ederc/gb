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
 * \file tools.c
 * \brief Implementation of smaller tools
 *
 * \author Christian Eder <ederc@mathematik.uni-kl.de>
 */

#include "data.h"

static inline int32_t mod_p_inverse_32(
        const int32_t val,
        const int32_t p
        )
{
    int32_t a, b, c, d, e, f;
    a =   p;
    b =   val % p;
    /* if b < 0 we shift correspondingly */
    b +=  (b >> 31) & p;
    c =   1;
    d =   0;

    while (b != 0) {
        f = b;
        e = a/f;
        b = a - e*f;
        a = f;
        f = c;
        c = d - e*f;
        d = f;
    }

    /* if d < 0 we shift correspondingly */
    d +=  (d >> 31) & p;

    return d;
}

static inline int16_t mod_p_inverse_16(
        const int16_t val,
        const int16_t p
        )
{
    int16_t a, b, c, d, e, f;
    a =   p;
    b =   val % p;
    /* if b < 0 we shift correspondingly */
    b +=  (b >> 15) & p;
    c =   1;
    d =   0;

    while (b != 0) {
        f = b;
        e = a/f;
        b = a - e*f;
        a = f;
        f = c;
        c = d - e*f;
        d = f;
    }

    /* if d < 0 we shift correspondingly */
    d +=  (d >> 15) & p;

    return d;
}

static inline val_t compare_and_swap(
        long *ptr,
        long old,
        long new
        )
{
    val_t prev;

#if 0
    __asm__ __volatile__(
            "lock; cmpxchgl %2, %1" : "=a"(prev),
            "+m"(*ptr) : "r"(new), "0"(old) : "memory");
    /* on which systems do we need "cmpxchgq" instead of "cmpxchgl" ? */
#else
    __asm__ __volatile__(
            "lock; cmpxchgq %2, %1" : "=a"(prev),
            "+m"(*ptr) : "r"(new), "0"(old) : "memory");
#endif
    return prev;
}
