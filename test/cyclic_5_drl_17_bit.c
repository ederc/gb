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
    const int32_t exps[]  = {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0};

    const int32_t nr_vars           = 5;
    const int32_t nr_gens           = 5;
    const int32_t ht_size           = 12;
    const int32_t field_char        = 101;
    const int32_t mon_order         = 0;
    const int32_t nr_threads        = 1;
    const int32_t info_level				=	2;
		const int32_t la_option         = 2;
    const int32_t reduce_gb         = 0;
    const int32_t pbm_file          = 0;
    const int32_t max_nr_pairs      = 0;
    const int32_t reset_hash_table  = 0;

    int32_t failure = 0;

    /* returned basis data as pointers for interfaces */
    int32_t *bld    = (int32_t *)malloc(sizeof(int32_t));
    int32_t **blen  = (int32_t **)malloc(sizeof(int32_t *));
    int32_t **bexp  = (int32_t **)malloc(sizeof(int32_t *));
    void **bcf      = (void **)malloc(sizeof(void *));

    /* f4 computation:
     * ret = 0 if computation correct
     * ret = 1 if something went wrong */
    int ret = f4_julia(
            bld, blen, bexp, bcf, lens, exps, cfs, field_char, mon_order, nr_vars,
            nr_gens, ht_size, nr_threads, max_nr_pairs, reset_hash_table,
            la_option, reduce_gb, pbm_file, info_level);

    int32_t nterms  = 0;
    for (i = 0; i < (*bld); ++i) {
        nterms += (*blen)[i];
    }
    if ((*bld) != 20 || nterms != 242) {
        failure = 1;
    }

    int32_t tlen[20] = {5, 6, 11, 13, 12, 17, 13, 15, 17, 15, 4, 14, 11, 10, 6, 15, 15, 15, 15, 13};
    int32_t tcf[242] = {1, 1, 1, 1, 1, 1, 1, 100, 2, 1, 1, 1, 100, 1, 100, 1, 1, 1, 100, 99, 1, 100, 1, 1, 1, 100, 1, 100, 100, 100, 99, 1, 1, 100, 1, 1, 51, 1, 50, 100, 50, 50, 51, 51, 100, 51, 50, 1, 30, 61, 100, 12, 60, 90, 91, 60, 99, 69, 42, 92, 41, 8, 10, 43, 1, 1, 99, 100, 100, 3, 99, 99, 98, 3, 3, 99, 2, 1, 28, 14, 80, 89, 89, 90, 3, 9, 17, 71, 96, 29, 73, 18, 1, 3, 1, 3, 1, 99, 3, 3, 97, 99, 99, 94, 1, 4, 100, 97, 2, 1, 5, 1, 5, 1, 97, 4, 4, 95, 98, 99, 90, 7, 95, 4, 1, 100, 100, 1, 1, 68, 66, 66, 34, 32, 2, 86, 20, 94, 1, 13, 81, 42, 1, 34, 34, 34, 67, 68, 33, 67, 100, 68, 100, 1, 51, 4, 51, 52, 100, 50, 97, 50, 49, 1, 3, 1, 100, 98, 100, 1, 42, 21, 37, 42, 46, 25, 46, 13, 71, 80, 85, 21, 59, 17, 1, 92, 46, 52, 60, 39, 29, 75, 67, 100, 19, 63, 82, 8, 75, 1, 85, 93, 63, 85, 21, 29, 21, 96, 50, 8, 30, 92, 16, 17, 1, 14, 6, 74, 2, 86, 100, 7, 91, 92, 68, 24, 33, 87, 22, 1, 99, 100, 5, 4, 100, 99, 2, 7, 97, 94, 2, 97};
    int32_t texp[1210] = {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 0, 1, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 1, 1, 2, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 2, 1, 0, 0, 0, 3, 1, 0, 1, 1, 0, 2, 0, 1, 0, 1, 2, 0, 0, 0, 2, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 1, 1, 1, 2, 0, 0, 2, 1, 2, 0, 0, 1, 2, 2, 0, 1, 1, 0, 3, 0, 1, 0, 1, 3, 0, 0, 1, 1, 3, 0, 0, 0, 2, 3, 0, 1, 0, 0, 4, 0, 0, 1, 0, 4, 0, 0, 0, 1, 4, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 3, 4, 0, 0, 2, 0, 5, 0, 1, 0, 1, 5, 0, 0, 1, 1, 5, 0, 0, 0, 2, 5, 0, 1, 0, 0, 6, 0, 0, 1, 0, 6, 0, 0, 0, 1, 6, 0, 0, 0, 0, 7, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 0, 2, 6, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 2, 1, 0, 0, 0, 3, 1, 0, 1, 1, 0, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1, 2, 0, 0, 1, 1, 2, 0, 0, 0, 2, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 0, 2, 2, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 2, 1, 0, 0, 0, 3, 1, 0, 1, 1, 0, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1, 2, 0, 0, 1, 1, 2, 0, 1, 0, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 1, 0, 0, 5, 0, 0, 1, 0, 5, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 3, 3, 0, 1, 1, 0, 4, 0, 0, 2, 0, 4, 0, 1, 0, 1, 4, 0, 0, 1, 1, 4, 0, 0, 0, 2, 4, 0, 0, 1, 0, 5, 0, 0, 0, 1, 5, 0, 0, 0, 0, 6, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 3, 0, 0, 0, 3, 3, 0, 1, 1, 0, 4, 0, 0, 2, 0, 4, 0, 1, 0, 1, 4, 0, 0, 1, 1, 4, 0, 0, 1, 0, 5, 0, 0, 0, 0, 6, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 5, 0, 0, 0, 2, 5, 0, 0, 1, 0, 6, 0, 0, 0, 1, 6, 0, 0, 0, 0, 7, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 2, 0, 5, 0, 0, 1, 0, 6, 0, 0, 0, 0, 7, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 8, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 0, 1, 7, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 0, 7, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 0, 4, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 2, 1, 0, 0, 0, 3, 1, 0, 1, 1, 0, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1, 2, 0, 0, 1, 1, 2, 0, 0, 0, 2, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 0, 3, 1, 0, 1, 1, 0, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1, 2, 0, 0, 0, 2, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4};

    for (i = 0; i < (*bld); ++i) {
        if (tlen[i] != (*blen)[i]) {
            failure = 1;
            break;
        }
    }
    int32_t **cf  = (int32_t **)bcf;
    for (i = 0; i < nterms; ++i) {
        if (tcf[i] != (*cf)[i]) {
            failure = 1;
            break;
        }
    }
    for (i = 0; i < nterms * nr_vars; ++i) {
        if (texp[i] != (*bexp)[i]) {
            failure = 1;
            break;
        }
    }
    free(*blen);
    free(*bexp);
    free(*cf);
    free(blen);
    free(bexp);
    free(bcf);
    free(bld);
    return failure;
}
