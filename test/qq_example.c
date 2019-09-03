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

    const int32_t nterms  = 22;
    mpz_t **cfs = (mpz_t **)malloc((unsigned long)(2*nterms)*sizeof(mpz_t *));
    for (i = 0; i < 2*nterms; ++i) {
        cfs[i]  = (mpz_t *)malloc(sizeof(mpz_t));
        mpz_init(*(cfs[i]));
    }
    mpz_set_si(*(cfs[0]), -3);
    mpz_set_si(*(cfs[1]), 1);
    mpz_set_si(*(cfs[2]), 1);
    mpz_set_si(*(cfs[3]), 1);
    mpz_set_si(*(cfs[4]), 1);
    mpz_set_si(*(cfs[5]), 1);
    mpz_set_si(*(cfs[6]), 1);
    mpz_set_si(*(cfs[7]), 1);
    mpz_set_si(*(cfs[8]), 1);
    mpz_set_si(*(cfs[9]), 1);
    mpz_set_si(*(cfs[10]), -1);
    mpz_set_si(*(cfs[11]), 1);
    mpz_set_si(*(cfs[12]), 1);
    mpz_set_si(*(cfs[13]), 1);
    mpz_set_si(*(cfs[14]), 1);
    mpz_set_si(*(cfs[15]), 1);
    mpz_set_si(*(cfs[16]), 1);
    mpz_set_si(*(cfs[17]), 1);
    mpz_set_si(*(cfs[18]), 1);
    mpz_set_si(*(cfs[19]), 1);
    mpz_set_si(*(cfs[20]), 1);
    mpz_set_si(*(cfs[21]), 1);
    mpz_set_si(*(cfs[22]), 1);
    mpz_set_si(*(cfs[23]), 1);
    mpz_set_si(*(cfs[24]), 1);
    mpz_set_si(*(cfs[25]), 1);
    mpz_set_si(*(cfs[26]), 1);
    mpz_set_si(*(cfs[27]), 1);
    mpz_set_si(*(cfs[28]), 1);
    mpz_set_si(*(cfs[29]), 1);
    mpz_set_si(*(cfs[30]), 3);
    mpz_set_si(*(cfs[31]), 4);
    mpz_set_si(*(cfs[32]), 2);
    mpz_set_si(*(cfs[33]), 1);
    mpz_set_si(*(cfs[34]), 1);
    mpz_set_si(*(cfs[35]), 1);
    mpz_set_si(*(cfs[36]), 1);
    mpz_set_si(*(cfs[37]), 1);
    mpz_set_si(*(cfs[38]), 1);
    mpz_set_si(*(cfs[39]), 1);
    mpz_set_si(*(cfs[40]), 14);
    mpz_set_si(*(cfs[41]), 1);
    mpz_set_si(*(cfs[42]), -2);
    mpz_set_si(*(cfs[43]), 1);

    const int32_t lens[]  = {5, 5, 5, 5, 2}; 
    const int32_t exps[]  = {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0};

    const int32_t nr_vars           = 5;
    const int32_t nr_gens           = 5;
    const int32_t ht_size           = 6;
    const int32_t field_char        = 0;
    const int32_t mon_order         = 0;
    const int32_t nr_threads        = 1;
    const int32_t info_level				=	2;
		const int32_t la_option         = 2;
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
            la_option, pbm_file, info_level);

    mpz_clears(
        *(cfs[0]),
        *(cfs[1]),
        *(cfs[2]),
        *(cfs[3]),
        *(cfs[4]),
        *(cfs[5]),
        *(cfs[6]),
        *(cfs[7]),
        *(cfs[8]),
        *(cfs[9]),
        *(cfs[10]),
        *(cfs[11]),
        *(cfs[12]),
        *(cfs[13]),
        *(cfs[14]),
        *(cfs[15]),
        *(cfs[16]),
        *(cfs[17]),
        *(cfs[18]),
        *(cfs[19]),
        *(cfs[20]),
        *(cfs[21]),
        *(cfs[22]),
        *(cfs[23]),
        *(cfs[24]),
        *(cfs[25]),
        *(cfs[26]),
        *(cfs[27]),
        *(cfs[28]),
        *(cfs[29]),
        *(cfs[30]),
        *(cfs[31]),
        *(cfs[32]),
        *(cfs[33]),
        *(cfs[34]),
        *(cfs[35]),
        *(cfs[36]),
        *(cfs[37]),
        *(cfs[38]),
        *(cfs[39]),
        *(cfs[40]),
        *(cfs[41]),
        *(cfs[42]),
        *(cfs[43]),
        NULL);

    for (i = 0; i < 2*nterms; ++i) {
        free(cfs[i]);
    }
    free(cfs);

    mpz_t **bcfs = (mpz_t **)bcf;
    printf("number terms %d\n", ret);
    for (i = 0; i < ret; ++i) {
        gmp_printf("%d ---- %Zd\n", i, (*bcfs)[i]);
    }

    free_julia_data(blen, bexp, bcf, *bld, field_char);

    free(blen);
    free(bexp);
    free(bld);

    return failure;
}
