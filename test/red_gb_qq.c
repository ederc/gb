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
    const int32_t nr_threads        = 4;
    const int32_t info_level				=	2;
		const int32_t la_option         = 1;
    const int32_t reduce_gb         = 1;
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

    if ((*bld) != 23 || ret != 318) {
        failure = 1;
    }

    mpz_t *tcfs = (mpz_t *)malloc((unsigned long)ret * sizeof(mpz_t));
    for (i = 0; i < ret; ++i) {
        mpz_init(tcfs[i]);
    }

    mpz_set_si(tcfs[0], 3);
    mpz_set_si(tcfs[1], -1);
    mpz_set_si(tcfs[2], -1);
    mpz_set_si(tcfs[3], -1);
    mpz_set_si(tcfs[4], -1);
    mpz_set_si(tcfs[5], 1);
    mpz_set_si(tcfs[6], -2);
    mpz_set_si(tcfs[7], 1);
    mpz_set_si(tcfs[8], -3);
    mpz_set_si(tcfs[9], -1);
    mpz_set_si(tcfs[10], -4);
    mpz_set_si(tcfs[11], -1);
    mpz_set_si(tcfs[12], 12);
    mpz_set_si(tcfs[13], 12);
    mpz_set_si(tcfs[14], 27);
    mpz_set_si(tcfs[15], 12);
    mpz_set_si(tcfs[16], -9);
    mpz_set_si(tcfs[17], -136);
    mpz_set_si(tcfs[18], -76);
    mpz_set_si(tcfs[19], -73);
    mpz_set_si(tcfs[20], -410);
    mpz_set_si(tcfs[21], -73);
    mpz_set_si(tcfs[22], -64);
    mpz_set_si(tcfs[23], -44);
    mpz_set_si(tcfs[24], -185);
    mpz_set_si(tcfs[25], -28);
    mpz_set_si(tcfs[26], 3);
    mpz_set_si(tcfs[27], 3);
    mpz_set_si(tcfs[28], 3);
    mpz_set_si(tcfs[29], 4);
    mpz_set_si(tcfs[30], 1);
    mpz_set_si(tcfs[31], 1);
    mpz_set_si(tcfs[32], 11);
    mpz_set_si(tcfs[33], 1);
    mpz_set_si(tcfs[34], 1);
    mpz_set_si(tcfs[35], 2);
    mpz_set_si(tcfs[36], 5);
    mpz_set_si(tcfs[37], 1);
    mpz_set_si(tcfs[38], 108);
    mpz_set_si(tcfs[39], 3079);
    mpz_set_si(tcfs[40], 3376);
    mpz_set_si(tcfs[41], -416);
    mpz_set_si(tcfs[42], 13628);
    mpz_set_si(tcfs[43], 2296);
    mpz_set_si(tcfs[44], 2712);
    mpz_set_si(tcfs[45], 1332);
    mpz_set_si(tcfs[46], 7681);
    mpz_set_si(tcfs[47], 4892);
    mpz_set_si(tcfs[48], 8537);
    mpz_set_si(tcfs[49], 2376);
    mpz_set_si(tcfs[50], 1224);
    mpz_set_si(tcfs[51], 12157);
    mpz_set_si(tcfs[52], 2508);
    mpz_set_si(tcfs[53], 9);
    mpz_set_si(tcfs[54], 5);
    mpz_set_si(tcfs[55], 20);
    mpz_set_si(tcfs[56], -14);
    mpz_set_si(tcfs[57], 66);
    mpz_set_si(tcfs[58], 22);
    mpz_set_si(tcfs[59], 32);
    mpz_set_si(tcfs[60], 8);
    mpz_set_si(tcfs[61], 62);
    mpz_set_si(tcfs[62], 33);
    mpz_set_si(tcfs[63], 42);
    mpz_set_si(tcfs[64], 24);
    mpz_set_si(tcfs[65], 8);
    mpz_set_si(tcfs[66], 110);
    mpz_set_si(tcfs[67], 24);
    mpz_set_si(tcfs[68], -36);
    mpz_set_si(tcfs[69], -111);
    mpz_set_si(tcfs[70], -240);
    mpz_set_si(tcfs[71], 144);
    mpz_set_si(tcfs[72], -852);
    mpz_set_si(tcfs[73], -120);
    mpz_set_si(tcfs[74], -280);
    mpz_set_si(tcfs[75], -100);
    mpz_set_si(tcfs[76], -661);
    mpz_set_si(tcfs[77], -452);
    mpz_set_si(tcfs[78], -673);
    mpz_set_si(tcfs[79], -232);
    mpz_set_si(tcfs[80], -104);
    mpz_set_si(tcfs[81], -1073);
    mpz_set_si(tcfs[82], -220);
    mpz_set_si(tcfs[83], 36);
    mpz_set_si(tcfs[84], 12);
    mpz_set_si(tcfs[85], -36);
    mpz_set_si(tcfs[86], 84);
    mpz_set_si(tcfs[87], 12);
    mpz_set_si(tcfs[88], 160);
    mpz_set_si(tcfs[89], 16);
    mpz_set_si(tcfs[90], 76);
    mpz_set_si(tcfs[91], 392);
    mpz_set_si(tcfs[92], 124);
    mpz_set_si(tcfs[93], 64);
    mpz_set_si(tcfs[94], 32);
    mpz_set_si(tcfs[95], 140);
    mpz_set_si(tcfs[96], 16);
    mpz_set_si(tcfs[97], 9);
    mpz_set_si(tcfs[98], 5);
    mpz_set_si(tcfs[99], -4);
    mpz_set_si(tcfs[100], 3);
    mpz_set_si(tcfs[101], 5);
    mpz_set_si(tcfs[102], 3);
    mpz_set_si(tcfs[103], 24);
    mpz_set_si(tcfs[104], 7);
    mpz_set_si(tcfs[105], 51);
    mpz_set_si(tcfs[106], 7);
    mpz_set_si(tcfs[107], 8);
    mpz_set_si(tcfs[108], 8);
    mpz_set_si(tcfs[109], 39);
    mpz_set_si(tcfs[110], 8);
    mpz_set_si(tcfs[111], 28);
    mpz_set_si(tcfs[112], 28);
    mpz_set_si(tcfs[113], 112);
    mpz_set_si(tcfs[114], -168);
    mpz_set_si(tcfs[115], -28);
    mpz_set_si(tcfs[116], -336);
    mpz_set_si(tcfs[117], 56);
    mpz_set_si(tcfs[118], -56);
    mpz_set_si(tcfs[119], -56);
    mpz_set_si(tcfs[120], -252);
    mpz_set_si(tcfs[121], -56);
    mpz_set_si(tcfs[122], 9);
    mpz_set_si(tcfs[123], 2632);
    mpz_set_si(tcfs[124], 5516);
    mpz_set_si(tcfs[125], 70168);
    mpz_set_si(tcfs[126], 3472);
    mpz_set_si(tcfs[127], 23744);
    mpz_set_si(tcfs[128], 120596);
    mpz_set_si(tcfs[129], 9968);
    mpz_set_si(tcfs[130], 26152);
    mpz_set_si(tcfs[131], 13272);
    mpz_set_si(tcfs[132], 113120);
    mpz_set_si(tcfs[133], 23016);
    mpz_set_si(tcfs[134], 954);
    mpz_set_si(tcfs[135], -636);
    mpz_set_si(tcfs[136], -279);
    mpz_set_si(tcfs[137], -1389);
    mpz_set_si(tcfs[138], 658);
    mpz_set_si(tcfs[139], 1190);
    mpz_set_si(tcfs[140], 18872);
    mpz_set_si(tcfs[141], 896);
    mpz_set_si(tcfs[142], 4578);
    mpz_set_si(tcfs[143], 31164);
    mpz_set_si(tcfs[144], 2254);
    mpz_set_si(tcfs[145], 6664);
    mpz_set_si(tcfs[146], 2576);
    mpz_set_si(tcfs[147], 28322);
    mpz_set_si(tcfs[148], 5600);
    mpz_set_si(tcfs[149], 369);
    mpz_set_si(tcfs[150], -246);
    mpz_set_si(tcfs[151], -72);
    mpz_set_si(tcfs[152], -222);
    mpz_set_si(tcfs[153], -7896);
    mpz_set_si(tcfs[154], 21504);
    mpz_set_si(tcfs[155], 288512);
    mpz_set_si(tcfs[156], 11144);
    mpz_set_si(tcfs[157], 95144);
    mpz_set_si(tcfs[158], 494116);
    mpz_set_si(tcfs[159], 34244);
    mpz_set_si(tcfs[160], 106736);
    mpz_set_si(tcfs[161], 56056);
    mpz_set_si(tcfs[162], 465472);
    mpz_set_si(tcfs[163], 95312);
    mpz_set_si(tcfs[164], 3717);
    mpz_set_si(tcfs[165], -2196);
    mpz_set_si(tcfs[166], -1107);
    mpz_set_si(tcfs[167], -6057);
    mpz_set_si(tcfs[168], 1001107128);
    mpz_set_si(tcfs[169], 4447371432);
    mpz_set_si(tcfs[170], 1610336);
    mpz_set_si(tcfs[171], -683518864);
    mpz_set_si(tcfs[172], -2286164608);
    mpz_set_si(tcfs[173], -920068240);
    mpz_set_si(tcfs[174], 317918085);
    mpz_set_si(tcfs[175], 41847252);
    mpz_set_si(tcfs[176], 189595305);
    mpz_set_si(tcfs[177], 320362143);
    mpz_set_si(tcfs[178], -68194323);
    mpz_set_si(tcfs[179], 161068047);
    mpz_set_si(tcfs[180], -238973091);
    mpz_set_si(tcfs[181], -209071785);
    mpz_set_si(tcfs[182], 10115403);
    mpz_set_si(tcfs[183], 1001107128);
    mpz_set_si(tcfs[184], 66836448);
    mpz_set_si(tcfs[185], 204715168);
    mpz_set_si(tcfs[186], 558992560);
    mpz_set_si(tcfs[187], 2317456792);
    mpz_set_si(tcfs[188], 626802736);
    mpz_set_si(tcfs[189], 28845165);
    mpz_set_si(tcfs[190], 12203004);
    mpz_set_si(tcfs[191], 3257109);
    mpz_set_si(tcfs[192], 34101075);
    mpz_set_si(tcfs[193], -6498639);
    mpz_set_si(tcfs[194], -29010345);
    mpz_set_si(tcfs[195], -8968359);
    mpz_set_si(tcfs[196], -70468665);
    mpz_set_si(tcfs[197], -87821457);
    mpz_set_si(tcfs[198], 667404752);
    mpz_set_si(tcfs[199], -11261651324);
    mpz_set_si(tcfs[200], 8558816);
    mpz_set_si(tcfs[201], 2203992112);
    mpz_set_si(tcfs[202], 4849319776);
    mpz_set_si(tcfs[203], 2074811424);
    mpz_set_si(tcfs[204], -954498078);
    mpz_set_si(tcfs[205], -79242024);
    mpz_set_si(tcfs[206], -613981590);
    mpz_set_si(tcfs[207], -766277541);
    mpz_set_si(tcfs[208], 199887966);
    mpz_set_si(tcfs[209], -399966876);
    mpz_set_si(tcfs[210], 710575101);
    mpz_set_si(tcfs[211], 831563652);
    mpz_set_si(tcfs[212], 112438875);
    mpz_set_si(tcfs[213], 4004428512);
    mpz_set_si(tcfs[214], -22886502660);
    mpz_set_si(tcfs[215], 1221118528);
    mpz_set_si(tcfs[216], 3254258224);
    mpz_set_si(tcfs[217], 16020592336);
    mpz_set_si(tcfs[218], 5465581744);
    mpz_set_si(tcfs[219], -1714625187);
    mpz_set_si(tcfs[220], -243006540);
    mpz_set_si(tcfs[221], -1000792071);
    mpz_set_si(tcfs[222], -1785712614);
    mpz_set_si(tcfs[223], 366510897);
    mpz_set_si(tcfs[224], -762530397);
    mpz_set_si(tcfs[225], 1185626490);
    mpz_set_si(tcfs[226], 1018036875);
    mpz_set_si(tcfs[227], 1774128);
    mpz_set_si(tcfs[228], 54059784912);
    mpz_set_si(tcfs[229], -775188286236);
    mpz_set_si(tcfs[230], -4181508576);
    mpz_set_si(tcfs[231], 163530786384);
    mpz_set_si(tcfs[232], 328173016752);
    mpz_set_si(tcfs[233], 145386641232);
    mpz_set_si(tcfs[234], -44514047907);
    mpz_set_si(tcfs[235], -7186528332);
    mpz_set_si(tcfs[236], -24199277391);
    mpz_set_si(tcfs[237], -46760893704);
    mpz_set_si(tcfs[238], 10050436017);
    mpz_set_si(tcfs[239], -27087982821);
    mpz_set_si(tcfs[240], 42429647772);
    mpz_set_si(tcfs[241], 46406970279);
    mpz_set_si(tcfs[242], -14966046810);
    mpz_set_si(tcfs[243], 1344);
    mpz_set_si(tcfs[244], 421514);
    mpz_set_si(tcfs[245], -64063);
    mpz_set_si(tcfs[246], -112722);
    mpz_set_si(tcfs[247], -221635);
    mpz_set_si(tcfs[248], 40242);
    mpz_set_si(tcfs[249], 812646);
    mpz_set_si(tcfs[250], -179553);
    mpz_set_si(tcfs[251], 869875);
    mpz_set_si(tcfs[252], 922221);
    mpz_set_si(tcfs[253], 48067);
    mpz_set_si(tcfs[254], -740602);
    mpz_set_si(tcfs[255], 1077353);
    mpz_set_si(tcfs[256], 138828);
    mpz_set_si(tcfs[257], 3962);
    mpz_set_si(tcfs[258], 336);
    mpz_set_si(tcfs[259], -24527);
    mpz_set_si(tcfs[260], 3802);
    mpz_set_si(tcfs[261], 6165);
    mpz_set_si(tcfs[262], 13165);
    mpz_set_si(tcfs[263], -2271);
    mpz_set_si(tcfs[264], -46268);
    mpz_set_si(tcfs[265], 10570);
    mpz_set_si(tcfs[266], -49677);
    mpz_set_si(tcfs[267], -52129);
    mpz_set_si(tcfs[268], -1941);
    mpz_set_si(tcfs[269], 42836);
    mpz_set_si(tcfs[270], -62790);
    mpz_set_si(tcfs[271], -8671);
    mpz_set_si(tcfs[272], -208);
    mpz_set_si(tcfs[273], 1344);
    mpz_set_si(tcfs[274], -66632);
    mpz_set_si(tcfs[275], 7669);
    mpz_set_si(tcfs[276], 26514);
    mpz_set_si(tcfs[277], 27841);
    mpz_set_si(tcfs[278], -7770);
    mpz_set_si(tcfs[279], -146164);
    mpz_set_si(tcfs[280], 26147);
    mpz_set_si(tcfs[281], -155315);
    mpz_set_si(tcfs[282], -174437);
    mpz_set_si(tcfs[283], -26627);
    mpz_set_si(tcfs[284], 125844);
    mpz_set_si(tcfs[285], -171133);
    mpz_set_si(tcfs[286], -8114);
    mpz_set_si(tcfs[287], -840);
    mpz_set_si(tcfs[288], 1344);
    mpz_set_si(tcfs[289], -3493);
    mpz_set_si(tcfs[290], -8329);
    mpz_set_si(tcfs[291], 39120);
    mpz_set_si(tcfs[292], -31213);
    mpz_set_si(tcfs[293], -6468);
    mpz_set_si(tcfs[294], -83487);
    mpz_set_si(tcfs[295], -9387);
    mpz_set_si(tcfs[296], -84812);
    mpz_set_si(tcfs[297], -132024);
    mpz_set_si(tcfs[298], -85988);
    mpz_set_si(tcfs[299], 39305);
    mpz_set_si(tcfs[300], -214);
    mpz_set_si(tcfs[301], 63273);
    mpz_set_si(tcfs[302], -6835);
    mpz_set_si(tcfs[303], 84);
    mpz_set_si(tcfs[304], 1445);
    mpz_set_si(tcfs[305], -241);
    mpz_set_si(tcfs[306], -327);
    mpz_set_si(tcfs[307], -814);
    mpz_set_si(tcfs[308], 129);
    mpz_set_si(tcfs[309], 2694);
    mpz_set_si(tcfs[310], -627);
    mpz_set_si(tcfs[311], 2878);
    mpz_set_si(tcfs[312], 2997);
    mpz_set_si(tcfs[313], 34);
    mpz_set_si(tcfs[314], -2458);
    mpz_set_si(tcfs[315], 3662);
    mpz_set_si(tcfs[316], 567);
    mpz_set_si(tcfs[317], 17);

    mpz_t **bcfs = (mpz_t **)bcf;

    for (i = 0; i < ret; ++i) {
        if (mpz_cmp(tcfs[i], (*bcfs)[i]) != 0) {
            failure = 1;
        }
    }

    for (i = 0; i < ret; ++i) {
        mpz_clear(tcfs[i]);
    }
    free(tcfs);

    free_julia_data(blen, bexp, bcf, *bld, field_char);

    free(blen);
    free(bexp);
    free(bld);

    printf("failure %d\n", failure);
    return failure;
}
