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

    if ((*bld) != 23 || ret != 325) {
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
    mpz_set_si(tcfs[12], 3);
    mpz_set_si(tcfs[13], 3);
    mpz_set_si(tcfs[14], 3);
    mpz_set_si(tcfs[15], 4);
    mpz_set_si(tcfs[16], 1);
    mpz_set_si(tcfs[17], 1);
    mpz_set_si(tcfs[18], 11);
    mpz_set_si(tcfs[19], 1);
    mpz_set_si(tcfs[20], 1);
    mpz_set_si(tcfs[21], 2);
    mpz_set_si(tcfs[22], 5);
    mpz_set_si(tcfs[23], 1);
    mpz_set_si(tcfs[24], 9);
    mpz_set_si(tcfs[25], 5);
    mpz_set_si(tcfs[26], -4);
    mpz_set_si(tcfs[27], 3);
    mpz_set_si(tcfs[28], 5);
    mpz_set_si(tcfs[29], 3);
    mpz_set_si(tcfs[30], 24);
    mpz_set_si(tcfs[31], 7);
    mpz_set_si(tcfs[32], 51);
    mpz_set_si(tcfs[33], 7);
    mpz_set_si(tcfs[34], 8);
    mpz_set_si(tcfs[35], 8);
    mpz_set_si(tcfs[36], 39);
    mpz_set_si(tcfs[37], 8);
    mpz_set_si(tcfs[38], 28);
    mpz_set_si(tcfs[39], 28);
    mpz_set_si(tcfs[40], 112);
    mpz_set_si(tcfs[41], -168);
    mpz_set_si(tcfs[42], -28);
    mpz_set_si(tcfs[43], -336);
    mpz_set_si(tcfs[44], 56);
    mpz_set_si(tcfs[45], -56);
    mpz_set_si(tcfs[46], -56);
    mpz_set_si(tcfs[47], -252);
    mpz_set_si(tcfs[48], -56);
    mpz_set_si(tcfs[49], 9);
    mpz_set_si(tcfs[50], 168);
    mpz_set_si(tcfs[51], -56);
    mpz_set_si(tcfs[52], 280);
    mpz_set_si(tcfs[53], 28);
    mpz_set_si(tcfs[54], -280);
    mpz_set_si(tcfs[55], 56);
    mpz_set_si(tcfs[56], 112);
    mpz_set_si(tcfs[57], -336);
    mpz_set_si(tcfs[58], 140);
    mpz_set_si(tcfs[59], -56);
    mpz_set_si(tcfs[60], -280);
    mpz_set_si(tcfs[61], -56);
    mpz_set_si(tcfs[62], -9);
    mpz_set_si(tcfs[63], 756);
    mpz_set_si(tcfs[64], 3080);
    mpz_set_si(tcfs[65], 392);
    mpz_set_si(tcfs[66], 1904);
    mpz_set_si(tcfs[67], 8008);
    mpz_set_si(tcfs[68], 644);
    mpz_set_si(tcfs[69], 1288);
    mpz_set_si(tcfs[70], 1680);
    mpz_set_si(tcfs[71], 5880);
    mpz_set_si(tcfs[72], 1288);
    mpz_set_si(tcfs[73], 96);
    mpz_set_si(tcfs[74], 24);
    mpz_set_si(tcfs[75], 54);
    mpz_set_si(tcfs[76], 150);
    mpz_set_si(tcfs[77], -18);
    mpz_set_si(tcfs[78], 33);
    mpz_set_si(tcfs[79], -6);
    mpz_set_si(tcfs[80], 114);
    mpz_set_si(tcfs[81], -228);
    mpz_set_si(tcfs[82], 12);
    mpz_set_si(tcfs[83], 12);
    mpz_set_si(tcfs[84], 27);
    mpz_set_si(tcfs[85], 12);
    mpz_set_si(tcfs[86], -9);
    mpz_set_si(tcfs[87], -136);
    mpz_set_si(tcfs[88], -76);
    mpz_set_si(tcfs[89], -73);
    mpz_set_si(tcfs[90], -410);
    mpz_set_si(tcfs[91], -73);
    mpz_set_si(tcfs[92], -64);
    mpz_set_si(tcfs[93], -44);
    mpz_set_si(tcfs[94], -185);
    mpz_set_si(tcfs[95], -28);
    mpz_set_si(tcfs[96], 84);
    mpz_set_si(tcfs[97], 1445);
    mpz_set_si(tcfs[98], -241);
    mpz_set_si(tcfs[99], -327);
    mpz_set_si(tcfs[100], -814);
    mpz_set_si(tcfs[101], 129);
    mpz_set_si(tcfs[102], 2694);
    mpz_set_si(tcfs[103], -627);
    mpz_set_si(tcfs[104], 2878);
    mpz_set_si(tcfs[105], 2997);
    mpz_set_si(tcfs[106], 34);
    mpz_set_si(tcfs[107], -2458);
    mpz_set_si(tcfs[108], 3662);
    mpz_set_si(tcfs[109], 567);
    mpz_set_si(tcfs[110], 17);
    mpz_set_si(tcfs[111], 9);
    mpz_set_si(tcfs[112], 198);
    mpz_set_si(tcfs[113], -27);
    mpz_set_si(tcfs[114], -632);
    mpz_set_si(tcfs[115], -344);
    mpz_set_si(tcfs[116], -240);
    mpz_set_si(tcfs[117], -1742);
    mpz_set_si(tcfs[118], -60);
    mpz_set_si(tcfs[119], 96);
    mpz_set_si(tcfs[120], -132);
    mpz_set_si(tcfs[121], -391);
    mpz_set_si(tcfs[122], -384);
    mpz_set_si(tcfs[123], -1042);
    mpz_set_si(tcfs[124], -8);
    mpz_set_si(tcfs[125], -104);
    mpz_set_si(tcfs[126], -351);
    mpz_set_si(tcfs[127], -44);
    mpz_set_si(tcfs[128], 36);
    mpz_set_si(tcfs[129], 873);
    mpz_set_si(tcfs[130], -108);
    mpz_set_si(tcfs[131], -2594);
    mpz_set_si(tcfs[132], -1424);
    mpz_set_si(tcfs[133], -978);
    mpz_set_si(tcfs[134], -7142);
    mpz_set_si(tcfs[135], -150);
    mpz_set_si(tcfs[136], 552);
    mpz_set_si(tcfs[137], -540);
    mpz_set_si(tcfs[138], -1591);
    mpz_set_si(tcfs[139], -1299);
    mpz_set_si(tcfs[140], -4339);
    mpz_set_si(tcfs[141], 16);
    mpz_set_si(tcfs[142], -416);
    mpz_set_si(tcfs[143], -1347);
    mpz_set_si(tcfs[144], -164);
    mpz_set_si(tcfs[145], 2632);
    mpz_set_si(tcfs[146], 5516);
    mpz_set_si(tcfs[147], 70168);
    mpz_set_si(tcfs[148], 3472);
    mpz_set_si(tcfs[149], 23744);
    mpz_set_si(tcfs[150], 120596);
    mpz_set_si(tcfs[151], 9968);
    mpz_set_si(tcfs[152], 26152);
    mpz_set_si(tcfs[153], 13272);
    mpz_set_si(tcfs[154], 113120);
    mpz_set_si(tcfs[155], 23016);
    mpz_set_si(tcfs[156], 954);
    mpz_set_si(tcfs[157], -636);
    mpz_set_si(tcfs[158], -279);
    mpz_set_si(tcfs[159], -1389);
    mpz_set_si(tcfs[160], 658);
    mpz_set_si(tcfs[161], 1190);
    mpz_set_si(tcfs[162], 18872);
    mpz_set_si(tcfs[163], 896);
    mpz_set_si(tcfs[164], 4578);
    mpz_set_si(tcfs[165], 31164);
    mpz_set_si(tcfs[166], 2254);
    mpz_set_si(tcfs[167], 6664);
    mpz_set_si(tcfs[168], 2576);
    mpz_set_si(tcfs[169], 28322);
    mpz_set_si(tcfs[170], 5600);
    mpz_set_si(tcfs[171], 369);
    mpz_set_si(tcfs[172], -246);
    mpz_set_si(tcfs[173], -72);
    mpz_set_si(tcfs[174], -222);
    mpz_set_si(tcfs[175], 1001107128);
    mpz_set_si(tcfs[176], 4447371432);
    mpz_set_si(tcfs[177], 1610336);
    mpz_set_si(tcfs[178], -683518864);
    mpz_set_si(tcfs[179], -2286164608);
    mpz_set_si(tcfs[180], -920068240);
    mpz_set_si(tcfs[181], 317918085);
    mpz_set_si(tcfs[182], 41847252);
    mpz_set_si(tcfs[183], 189595305);
    mpz_set_si(tcfs[184], 320362143);
    mpz_set_si(tcfs[185], -68194323);
    mpz_set_si(tcfs[186], 161068047);
    mpz_set_si(tcfs[187], -238973091);
    mpz_set_si(tcfs[188], -209071785);
    mpz_set_si(tcfs[189], 10115403);
    mpz_set_si(tcfs[190], 1001107128);
    mpz_set_si(tcfs[191], 66836448);
    mpz_set_si(tcfs[192], 204715168);
    mpz_set_si(tcfs[193], 558992560);
    mpz_set_si(tcfs[194], 2317456792);
    mpz_set_si(tcfs[195], 626802736);
    mpz_set_si(tcfs[196], 28845165);
    mpz_set_si(tcfs[197], 12203004);
    mpz_set_si(tcfs[198], 3257109);
    mpz_set_si(tcfs[199], 34101075);
    mpz_set_si(tcfs[200], -6498639);
    mpz_set_si(tcfs[201], -29010345);
    mpz_set_si(tcfs[202], -8968359);
    mpz_set_si(tcfs[203], -70468665);
    mpz_set_si(tcfs[204], -87821457);
    mpz_set_si(tcfs[205], 667404752);
    mpz_set_si(tcfs[206], -11261651324);
    mpz_set_si(tcfs[207], 8558816);
    mpz_set_si(tcfs[208], 2203992112);
    mpz_set_si(tcfs[209], 4849319776);
    mpz_set_si(tcfs[210], 2074811424);
    mpz_set_si(tcfs[211], -954498078);
    mpz_set_si(tcfs[212], -79242024);
    mpz_set_si(tcfs[213], -613981590);
    mpz_set_si(tcfs[214], -766277541);
    mpz_set_si(tcfs[215], 199887966);
    mpz_set_si(tcfs[216], -399966876);
    mpz_set_si(tcfs[217], 710575101);
    mpz_set_si(tcfs[218], 831563652);
    mpz_set_si(tcfs[219], 112438875);
    mpz_set_si(tcfs[220], 4004428512);
    mpz_set_si(tcfs[221], -22886502660);
    mpz_set_si(tcfs[222], 1221118528);
    mpz_set_si(tcfs[223], 3254258224);
    mpz_set_si(tcfs[224], 16020592336);
    mpz_set_si(tcfs[225], 5465581744);
    mpz_set_si(tcfs[226], -1714625187);
    mpz_set_si(tcfs[227], -243006540);
    mpz_set_si(tcfs[228], -1000792071);
    mpz_set_si(tcfs[229], -1785712614);
    mpz_set_si(tcfs[230], 366510897);
    mpz_set_si(tcfs[231], -762530397);
    mpz_set_si(tcfs[232], 1185626490);
    mpz_set_si(tcfs[233], 1018036875);
    mpz_set_si(tcfs[234], 1774128);
    mpz_set_si(tcfs[235], 1344);
    mpz_set_si(tcfs[236], 421514);
    mpz_set_si(tcfs[237], -64063);
    mpz_set_si(tcfs[238], -112722);
    mpz_set_si(tcfs[239], -221635);
    mpz_set_si(tcfs[240], 40242);
    mpz_set_si(tcfs[241], 812646);
    mpz_set_si(tcfs[242], -179553);
    mpz_set_si(tcfs[243], 869875);
    mpz_set_si(tcfs[244], 922221);
    mpz_set_si(tcfs[245], 48067);
    mpz_set_si(tcfs[246], -740602);
    mpz_set_si(tcfs[247], 1077353);
    mpz_set_si(tcfs[248], 138828);
    mpz_set_si(tcfs[249], 3962);
    mpz_set_si(tcfs[250], 336);
    mpz_set_si(tcfs[251], -24527);
    mpz_set_si(tcfs[252], 3802);
    mpz_set_si(tcfs[253], 6165);
    mpz_set_si(tcfs[254], 13165);
    mpz_set_si(tcfs[255], -2271);
    mpz_set_si(tcfs[256], -46268);
    mpz_set_si(tcfs[257], 10570);
    mpz_set_si(tcfs[258], -49677);
    mpz_set_si(tcfs[259], -52129);
    mpz_set_si(tcfs[260], -1941);
    mpz_set_si(tcfs[261], 42836);
    mpz_set_si(tcfs[262], -62790);
    mpz_set_si(tcfs[263], -8671);
    mpz_set_si(tcfs[264], -208);
    mpz_set_si(tcfs[265], 1344);
    mpz_set_si(tcfs[266], -66632);
    mpz_set_si(tcfs[267], 7669);
    mpz_set_si(tcfs[268], 26514);
    mpz_set_si(tcfs[269], 27841);
    mpz_set_si(tcfs[270], -7770);
    mpz_set_si(tcfs[271], -146164);
    mpz_set_si(tcfs[272], 26147);
    mpz_set_si(tcfs[273], -155315);
    mpz_set_si(tcfs[274], -174437);
    mpz_set_si(tcfs[275], -26627);
    mpz_set_si(tcfs[276], 125844);
    mpz_set_si(tcfs[277], -171133);
    mpz_set_si(tcfs[278], -8114);
    mpz_set_si(tcfs[279], -840);
    mpz_set_si(tcfs[280], 1344);
    mpz_set_si(tcfs[281], -3493);
    mpz_set_si(tcfs[282], -8329);
    mpz_set_si(tcfs[283], 39120);
    mpz_set_si(tcfs[284], -31213);
    mpz_set_si(tcfs[285], -6468);
    mpz_set_si(tcfs[286], -83487);
    mpz_set_si(tcfs[287], -9387);
    mpz_set_si(tcfs[288], -84812);
    mpz_set_si(tcfs[289], -132024);
    mpz_set_si(tcfs[290], -85988);
    mpz_set_si(tcfs[291], 39305);
    mpz_set_si(tcfs[292], -214);
    mpz_set_si(tcfs[293], 63273);
    mpz_set_si(tcfs[294], -6835);
    mpz_set_si(tcfs[295], 108);
    mpz_set_si(tcfs[296], 3079);
    mpz_set_si(tcfs[297], 3376);
    mpz_set_si(tcfs[298], -416);
    mpz_set_si(tcfs[299], 13628);
    mpz_set_si(tcfs[300], 2296);
    mpz_set_si(tcfs[301], 2712);
    mpz_set_si(tcfs[302], 1332);
    mpz_set_si(tcfs[303], 7681);
    mpz_set_si(tcfs[304], 4892);
    mpz_set_si(tcfs[305], 8537);
    mpz_set_si(tcfs[306], 2376);
    mpz_set_si(tcfs[307], 1224);
    mpz_set_si(tcfs[308], 12157);
    mpz_set_si(tcfs[309], 2508);
    mpz_set_si(tcfs[310], 9);
    mpz_set_si(tcfs[311], 5);
    mpz_set_si(tcfs[312], 20);
    mpz_set_si(tcfs[313], -14);
    mpz_set_si(tcfs[314], 66);
    mpz_set_si(tcfs[315], 22);
    mpz_set_si(tcfs[316], 32);
    mpz_set_si(tcfs[317], 8);
    mpz_set_si(tcfs[318], 62);
    mpz_set_si(tcfs[319], 33);
    mpz_set_si(tcfs[320], 42);
    mpz_set_si(tcfs[321], 24);
    mpz_set_si(tcfs[322], 8);
    mpz_set_si(tcfs[323], 110);
    mpz_set_si(tcfs[324], 24);

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
