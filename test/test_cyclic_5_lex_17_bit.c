#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    int32_t i, j, k;
    len_t *hcm;

    int32_t round = 0;

    const int32_t lens[]  = {5, 5, 5, 5, 2}; 
    const int32_t cfs[]   = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1};
    const int32_t exps[]  = {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0};

    const int32_t nr_vars           = 5;
    const int32_t nr_gens           = 5;
    const int32_t ht_size           = 12;
    const int32_t field_char        = 65521;
    const int32_t mon_order         = 1;
    const int32_t nr_threads        = 1;
    const int32_t info_level				=	2;
		const int32_t la_option         = 1;
    const int32_t pbm_file          = 0;
    const int32_t max_nr_pairs      = 0;
    const int32_t reset_hash_table  = 0;

    int32_t failure = 0;

    int32_t **basis = (int32_t **)malloc(sizeof(int32_t *));
    int64_t len     = f4_julia(
            basis, lens, cfs, exps, field_char, mon_order, nr_vars,
            nr_gens, ht_size, nr_threads, max_nr_pairs, reset_hash_table,
            la_option, pbm_file, info_level);

    if (len != 690) {
        failure = 1;
        free(*basis);
        free(basis);
        basis = NULL;
        return failure;
    }
    int32_t val[690]  = {11, 30, 36, 108, 78, 72, 90, 60, 60, 72, 48, 24, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 2, 0, 0, 0, 1, 0, 1, 0, 1, 0, 2, 0, 1, 0, 0, 1, 65520, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 5, 65520, 0, 1, 0, 0, 0, 32760, 0, 0, 2, 4, 0, 32760, 0, 0, 2, 2, 2, 32761, 0, 0, 2, 1, 3, 32761, 0, 0, 2, 0, 4, 32759, 0, 0, 1, 4, 1, 65518, 0, 0, 1, 3, 2, 32763, 0, 0, 1, 1, 4, 32762, 0, 0, 1, 0, 5, 32761, 0, 0, 1, 0, 0, 32760, 0, 0, 0, 4, 2, 65519, 0, 0, 0, 3, 3, 65519, 0, 0, 0, 2, 4, 32761, 0, 0, 0, 1, 5, 32762, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 6, 32762, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 65520, 0, 1, 0, 0, 1, 30573, 0, 0, 1, 3, 3, 39314, 0, 0, 1, 2, 4, 8736, 0, 0, 1, 1, 5, 61154, 0, 0, 1, 1, 0, 56786, 0, 0, 1, 0, 1, 13103, 0, 0, 0, 5, 2, 26206, 0, 0, 0, 4, 3, 61153, 0, 0, 0, 3, 4, 17474, 0, 0, 0, 2, 5, 56782, 0, 0, 0, 2, 0, 21845, 0, 0, 0, 0, 2, 1, 0, 1, 1, 0, 0, 65520, 0, 1, 0, 0, 1, 1, 0, 0, 2, 0, 0, 65519, 0, 0, 1, 3, 3, 1, 0, 0, 1, 2, 4, 1, 0, 0, 1, 1, 0, 2, 0, 0, 1, 0, 1, 65520, 0, 0, 0, 5, 2, 65519, 0, 0, 0, 4, 3, 1, 0, 0, 0, 2, 5, 65518, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 3, 0, 0, 2, 0, 0, 2, 0, 1, 52943, 0, 0, 1, 2, 0, 60620, 0, 0, 1, 1, 6, 37739, 0, 0, 1, 1, 1, 45259, 0, 0, 1, 0, 2, 24344, 0, 0, 0, 6, 2, 56854, 0, 0, 0, 5, 3, 12576, 0, 0, 0, 4, 4, 11938, 0, 0, 0, 3, 5, 20745, 0, 0, 0, 3, 0, 12912, 0, 0, 0, 2, 6, 7350, 0, 0, 0, 2, 1, 41176, 0, 0, 0, 1, 2, 8667, 0, 0, 0, 0, 3, 1, 0, 0, 1, 0, 5, 65520, 0, 0, 1, 0, 0, 26211, 0, 0, 0, 3, 8, 39310, 0, 0, 0, 3, 3, 52418, 0, 0, 0, 2, 9, 13103, 0, 0, 0, 2, 4, 13105, 0, 0, 0, 1, 5, 52416, 0, 0, 0, 1, 0, 39317, 0, 0, 0, 0, 6, 26204, 0, 0, 0, 0, 1, 1, 0, 0, 0, 7, 0, 3, 0, 0, 0, 6, 1, 1, 0, 0, 0, 5, 2, 10, 0, 0, 0, 2, 10, 65377, 0, 0, 0, 2, 5, 133, 0, 0, 0, 2, 0, 65199, 0, 0, 0, 1, 6, 319, 0, 0, 0, 1, 1, 76, 0, 0, 0, 0, 7, 65444, 0, 0, 0, 0, 2, 1, 0, 0, 1, 1, 0, 65520, 0, 0, 1, 0, 1, 26210, 0, 0, 0, 6, 1, 39317, 0, 0, 0, 5, 2, 65520, 0, 0, 0, 3, 4, 4291, 0, 0, 0, 2, 10, 3900, 0, 0, 0, 2, 5, 57331, 0, 0, 0, 2, 0, 6477, 0, 0, 0, 1, 6, 32834, 0, 0, 0, 1, 1, 19674, 0, 0, 0, 0, 7, 6530, 0, 0, 0, 0, 2, 1, 0, 0, 0, 2, 5, 65520, 0, 0, 0, 2, 0, 8339, 0, 0, 0, 1, 11, 13100, 0, 0, 0, 1, 6, 44082, 0, 0, 0, 1, 1, 33356, 0, 0, 0, 0, 12, 52399, 0, 0, 0, 0, 7, 45287, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 15, 122, 0, 0, 0, 0, 10, 65399, 0, 0, 0, 0, 5, 65520, 0, 0, 0, 0, 0};

    for (i = 0; i < len; ++i) {
        if (val[i] != (*basis)[i]) {
            printf("difference at position %d: %d -- %d\n", i, val[i], (*basis)[i]);
            failure = 1;
            break;
        }
    }

    free(*basis);
    free(basis);
    basis = NULL;

    return failure;
}
