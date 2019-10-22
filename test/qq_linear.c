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

    const int32_t nterms  = 6;
    mpz_t **cfs = (mpz_t **)malloc((unsigned long)(2*nterms)*sizeof(mpz_t *));
    for (i = 0; i < 2*nterms; ++i) {
        cfs[i]  = (mpz_t *)malloc(sizeof(mpz_t));
        mpz_init(*(cfs[i]));
    }
    mpz_set_si(*(cfs[0]), 1);
    mpz_set_si(*(cfs[1]), 1);
    mpz_set_si(*(cfs[2]), 1);
    mpz_set_si(*(cfs[3]), 1);
    mpz_set_si(*(cfs[4]), 2);
    mpz_set_si(*(cfs[5]), 1);

    mpz_set_si(*(cfs[6]), 1);
    mpz_set_si(*(cfs[7]), 1);
    mpz_set_si(*(cfs[8]), 3);
    mpz_set_si(*(cfs[9]), 1);
    mpz_set_si(*(cfs[10]), 1);
    mpz_set_si(*(cfs[11]), 1);

    const int32_t lens[]  = {3,3,3}; 
    const int32_t exps[]  = {1, 0, 0, 1, 0, 0, 
                             1, 0, 0, 1, 0, 0};

    const int32_t nr_vars           = 2;
    const int32_t nr_gens           = 2;
    const int32_t ht_size           = 6;
    const int32_t field_char        = 0;
    const int32_t mon_order         = 0;
    const int32_t nr_threads        = 1;
    const int32_t info_level				=	2;
		const int32_t la_option         = 1;
    const int32_t reduce_gb         = 0;
    const int32_t max_nr_pairs      = 0;
    const int32_t pbm_file          = 0;
    const int32_t reset_hash_table  = 0;

    int32_t failure = 0;

    /* returned basis data as pointers for interfaces */
    int32_t *bld    = (int32_t *)malloc(sizeof(int32_t));
    int32_t **blen  = (int32_t **)malloc(sizeof(int32_t *));
    int32_t **bexp  = (int32_t **)malloc(sizeof(int32_t *));
    void **bcf      = (void **)malloc(sizeof(void *));

    int ret = f4_julia(
            bld, blen, bexp, bcf, lens, exps, cfs, field_char, mon_order, nr_vars,
            nr_gens, ht_size, nr_threads, max_nr_pairs, reset_hash_table,
            la_option, reduce_gb, pbm_file, info_level);

    for (i = 0; i < 2*nterms; ++i) {
        mpz_clear(*(cfs[i]));
        free(cfs[i]);
    }
    free(cfs);


    free_julia_data(blen, bexp, bcf, *bld, field_char);

    free(blen);
    free(bexp);
    free(bld);

    return failure;
}
