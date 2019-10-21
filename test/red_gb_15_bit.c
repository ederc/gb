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
    const int32_t field_char        = 32003;
    const int32_t mon_order         = 0;
    const int32_t nr_threads        = 1;
    const int32_t info_level				=	2;
		const int32_t la_option         = 2;
    const int32_t reduce_gb         = 1;
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
    if ((*bld) != 20 || nterms != 232) {
        failure = 1;
    }

    int32_t tlen[20] = {5, 6, 13, 11, 15, 13, 14, 12, 13, 12, 4, 14, 11, 10, 6, 13, 15, 15, 15, 15};
    int32_t tcf[232] = {1, 1, 1, 1, 1, 1, 1, 32002, 2, 1, 1, 1, 1, 32001, 32002, 32002, 3, 32001, 32001, 32000, 3, 3, 32001, 2, 1, 32002, 1, 32002, 1, 1, 1, 32002, 32001, 1, 32002, 1, 14, 6, 31976, 2, 31988, 32002, 7, 31993, 31994, 31970, 24, 33, 31989, 22, 1, 32001, 32002, 5, 4, 32002, 32001, 2, 7, 31999, 31996, 2, 31999, 1, 31998, 32001, 10, 1, 6, 32000, 2, 2, 13, 31995, 31990, 4, 31995, 1, 1, 32001, 2, 32002, 32001, 32001, 32001, 3, 2, 32001, 2, 1, 1, 1, 32002, 1, 32002, 32002, 32002, 32001, 1, 1, 32002, 1, 1, 16002, 1, 16001, 32002, 16001, 16001, 16002, 16002, 32002, 16002, 16001, 1, 32002, 32002, 1, 1, 21336, 21334, 21334, 10668, 10666, 2, 2132, 12801, 17068, 1, 29869, 19202, 25604, 1, 10668, 10668, 10668, 21335, 21336, 10667, 21335, 32002, 21336, 32002, 1, 16002, 4, 16002, 16003, 32002, 16001, 31999, 16001, 16000, 1, 3, 1, 32002, 32000, 32002, 1, 19204, 12803, 6401, 25603, 32001, 32001, 2, 32002, 12802, 19197, 25602, 6402, 1, 42, 21, 31838, 42, 31948, 31927, 31948, 13, 31872, 31982, 186, 21, 31961, 219, 1, 31893, 31948, 52, 60, 39, 29, 31977, 31969, 31901, 120, 63, 31883, 109, 31977, 1, 31987, 31995, 63, 31987, 21, 29, 21, 31998, 50, 8, 31932, 31994, 16, 31919, 1, 28, 14, 31982, 31991, 31991, 31992, 3, 9, 17, 31973, 31998, 29, 31975, 18};
    int32_t texp[1210] = {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 0, 0, 3, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 1, 2, 0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 0, 4, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 2, 1, 0, 0, 0, 3, 1, 0, 1, 1, 0, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1, 2, 0, 0, 1, 1, 2, 0, 0, 0, 2, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 0, 3, 1, 0, 1, 1, 0, 2, 0, 0, 2, 0, 2, 0, 1, 0, 1, 2, 0, 0, 0, 2, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 1, 0, 3, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 2, 1, 0, 0, 0, 3, 1, 0, 0, 2, 0, 2, 0, 1, 0, 1, 2, 0, 0, 1, 1, 2, 0, 0, 0, 2, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 0, 2, 2, 0, 0, 1, 1, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 2, 1, 0, 0, 0, 3, 1, 0, 1, 0, 1, 2, 0, 0, 1, 1, 2, 0, 0, 0, 2, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 1, 1, 2, 0, 0, 1, 1, 1, 1, 0, 0, 2, 1, 1, 0, 1, 0, 2, 1, 0, 0, 1, 2, 1, 0, 0, 0, 3, 1, 0, 1, 1, 0, 2, 0, 1, 0, 1, 2, 0, 0, 0, 2, 2, 0, 1, 0, 0, 3, 0, 0, 1, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 4, 0, 1, 1, 1, 2, 0, 0, 2, 1, 2, 0, 0, 1, 2, 2, 0, 1, 1, 0, 3, 0, 1, 0, 1, 3, 0, 0, 1, 1, 3, 0, 0, 0, 2, 3, 0, 1, 0, 0, 4, 0, 0, 1, 0, 4, 0, 0, 0, 1, 4, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 1, 0, 0, 5, 0, 0, 1, 0, 5, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 2, 3, 0, 0, 0, 3, 3, 0, 1, 1, 0, 4, 0, 0, 2, 0, 4, 0, 1, 0, 1, 4, 0, 0, 1, 1, 4, 0, 0, 0, 2, 4, 0, 0, 1, 0, 5, 0, 0, 0, 1, 5, 0, 0, 0, 0, 6, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 3, 0, 0, 0, 3, 3, 0, 1, 1, 0, 4, 0, 0, 2, 0, 4, 0, 1, 0, 1, 4, 0, 0, 1, 1, 4, 0, 0, 1, 0, 5, 0, 0, 0, 0, 6, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 5, 0, 0, 0, 2, 5, 0, 0, 1, 0, 6, 0, 0, 0, 1, 6, 0, 0, 0, 0, 7, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 2, 0, 5, 0, 0, 1, 0, 6, 0, 0, 0, 0, 7, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 3, 4, 0, 0, 0, 2, 5, 0, 0, 1, 0, 6, 0, 0, 0, 1, 6, 0, 0, 0, 0, 7, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 8, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 0, 1, 7, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 0, 7, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 0, 2, 6, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 1, 0, 0, 2, 0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 2, 1, 0, 1, 0, 0, 2, 0, 0, 1, 0, 2, 0, 0, 0, 1, 2, 0, 0, 0, 0, 3};

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
