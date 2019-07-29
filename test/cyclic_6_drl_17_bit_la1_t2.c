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

    const int32_t lens[]  = {6, 6, 6, 6, 6, 2}; 
    const int32_t cfs[]   = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1};
    const int32_t exps[]  = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0};

    const int32_t nr_vars           = 6;
    const int32_t nr_gens           = 6;
    const int32_t ht_size           = 12;
    const int32_t field_char        = 65521;
    const int32_t mon_order         = 0;
    /* const int32_t nr_threads        = 1; */
    const int32_t nr_threads        = 2;
    const int32_t info_level				=	2;
    const int32_t pbm_file          = 0;
		const int32_t la_option         = 1;
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
            la_option, pbm_file, info_level);

    int32_t nterms  = 0;
    for (i = 0; i < (*bld); ++i) {
        nterms += (*blen)[i];
    }
    if ((*bld) != 45 || nterms != 1161) {
        failure = 1;
    }

    int32_t tlen[45] = {6, 9, 13, 23, 15, 19, 24, 31, 31, 31, 31, 31, 31, 31, 31, 28, 25, 27, 27, 27, 27, 26, 26, 26, 27, 27, 27, 27, 27, 27, 26, 27, 26, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27, 27};
    int32_t tcf[1161] = {1, 1, 1, 1, 1, 1, 1, 1, 65520, 1, 65520, 2, 1, 1, 1, 1, 65520, 1, 65520, 1, 1, 1, 1, 65520, 65519, 65520, 1, 65520, 1, 28081, 28081, 4680, 14040, 9360, 4680, 37441, 37439, 14040, 56161, 28080, 18720, 9360, 14041, 56162, 9360, 42120, 46800, 37441, 46800, 28081, 18720, 1, 2, 65519, 1, 65519, 65519, 65520, 2, 3, 2, 65518, 65520, 65520, 2, 65520, 1, 1, 32758, 32759, 1, 32760, 2, 32759, 1, 1, 1, 32761, 65520, 65520, 65520, 32760, 1, 1, 65520, 1, 37438, 37440, 28083, 18722, 56160, 28083, 28077, 28077, 18722, 9361, 37440, 46799, 56163, 18722, 9361, 56160, 56160, 18716, 28083, 18719, 37444, 46799, 1, 1, 54908, 18365, 55036, 10168, 50777, 18565, 18439, 65386, 46684, 11228, 625, 15495, 3361, 15228, 5610, 14754, 37521, 56223, 11160, 46838, 13557, 23178, 45734, 45411, 61073, 26364, 61173, 22126, 37521, 24785, 1, 20972, 42959, 27391, 54844, 41689, 15607, 9188, 63324, 45521, 19562, 36311, 1043, 10517, 62173, 11891, 19611, 46247, 64050, 53337, 39277, 52694, 33229, 25917, 46856, 37432, 58429, 56516, 64286, 46247, 6736, 1, 63233, 62026, 57769, 21415, 44295, 62797, 2742, 62018, 61236, 22539, 29224, 44983, 42762, 37980, 39548, 30196, 52646, 53427, 16308, 17159, 25568, 32964, 45671, 50367, 60248, 36342, 27871, 40810, 52646, 48108, 1, 56306, 38701, 52493, 22561, 49188, 60792, 64604, 41524, 22260, 56, 18646, 31169, 2190, 48696, 18882, 44141, 6910, 26175, 39672, 36352, 42650, 34754, 13147, 60394, 2414, 32755, 46221, 19545, 6910, 42706, 1, 50838, 48554, 26691, 33677, 49442, 19963, 44114, 16512, 13067, 9702, 2654, 60111, 7622, 27619, 1023, 23709, 40642, 44087, 14136, 11826, 64206, 32907, 5267, 58003, 58971, 26488, 11921, 64511, 40643, 8387, 1, 43216, 15061, 22689, 61390, 22370, 38965, 59493, 64116, 23741, 6841, 6473, 29094, 56554, 26287, 57981, 65088, 38281, 13135, 27942, 3495, 52706, 42939, 35813, 32443, 63236, 10350, 42426, 53903, 38280, 59548, 1, 10733, 7319, 40032, 50523, 8645, 43984, 14692, 43493, 10232, 19262, 22743, 10621, 62862, 64807, 29644, 16831, 37553, 5293, 11857, 50473, 58864, 54023, 12592, 44863, 55690, 25715, 8507, 45299, 37553, 12588, 1, 9220, 26815, 13031, 42965, 16333, 4730, 916, 24000, 43259, 32700, 46875, 34346, 30577, 16831, 13893, 54152, 58601, 39342, 58625, 29166, 22873, 63521, 52346, 37878, 30357, 32766, 19302, 45978, 58600, 22816, 1, 60324, 15663, 5673, 23844, 1266, 22083, 31495, 23177, 32737, 28964, 33430, 5209, 22760, 5109, 49813, 14999, 59618, 55367, 46747, 27195, 12325, 33109, 2462, 1074, 38347, 5108, 62832, 1, 50961, 50959, 32761, 61877, 3634, 10918, 43682, 32761, 36394, 50961, 40044, 58242, 54599, 32764, 21840, 10924, 36404, 7281, 65520, 36403, 21839, 40044, 25480, 25480, 1, 46177, 4877, 62712, 45778, 25862, 42896, 15208, 41142, 41141, 30991, 44283, 10827, 23600, 28454, 59047, 64678, 40257, 21564, 43808, 63091, 9414, 34139, 63907, 22077, 57589, 39295, 1, 55091, 847, 45963, 29016, 24218, 29208, 50657, 49524, 49525, 24559, 36495, 7375, 51158, 35244, 30363, 42577, 50184, 3169, 41371, 49679, 63413, 2551, 48602, 60939, 31052, 4513, 1, 17036, 41568, 43562, 4670, 52469, 26689, 33362, 41694, 41694, 35798, 19655, 41491, 20274, 46050, 33151, 33135, 1515, 9018, 51568, 39984, 3846, 10055, 21843, 49211, 23598, 43315, 1, 44683, 59906, 1072, 5247, 14199, 3074, 1298, 53987, 53987, 8010, 15223, 20070, 52754, 11120, 51087, 52235, 10525, 44557, 33023, 17925, 26513, 21357, 16956, 604, 29047, 6750, 1, 40651, 40651, 20319, 60970, 36100, 60970, 41571, 13036, 16701, 60970, 63095, 38219, 38219, 21840, 13349, 40651, 60344, 35474, 10604, 35474, 14269, 16698, 54920, 35474, 46724, 1, 6149, 6149, 36787, 42936, 49085, 42936, 41766, 33647, 47915, 42936, 20478, 27575, 27574, 43836, 33724, 6149, 42779, 48928, 55077, 48928, 3820, 48383, 9969, 48928, 35318, 1, 55515, 55514, 3801, 59315, 49308, 59315, 18414, 64870, 8407, 59315, 39261, 38060, 38061, 12389, 28053, 55514, 1904, 57418, 47411, 57418, 56473, 45575, 46466, 57418, 33140, 1, 38018, 22150, 35659, 29079, 39405, 5091, 21203, 30823, 3388, 18041, 51774, 53819, 5064, 5064, 36795, 17824, 64772, 20898, 24675, 34928, 32440, 8383, 38792, 4810, 28139, 49696, 1, 32568, 57947, 8120, 21973, 64497, 51841, 60005, 11103, 26748, 1694, 5199, 8446, 16778, 16778, 20756, 31610, 64659, 65504, 47099, 62419, 13658, 60449, 18995, 64849, 11772, 6305, 1, 30096, 10440, 33601, 12791, 53784, 26065, 40272, 44232, 26278, 11855, 34933, 13152, 49343, 49343, 48970, 13560, 6577, 48597, 47240, 26445, 62347, 58150, 62248, 23377, 36790, 46807, 1, 9335, 21257, 3093, 6666, 64990, 23948, 25519, 31246, 48129, 39366, 47707, 17284, 54127, 54127, 50641, 25636, 55193, 24695, 61341, 51312, 19455, 29188, 6549, 20988, 21301, 38679, 1, 58907, 16502, 60872, 27006, 13614, 1180, 2522, 30576, 28129, 49209, 20265, 61249, 12802, 12802, 515, 30340, 38690, 23379, 23386, 56459, 1235, 36764, 57476, 47479, 14469, 60424, 1, 49734, 60053, 36124, 60816, 28675, 62779, 18762, 12178, 997, 26429, 43017, 58542, 57802, 57802, 16188, 50173, 16300, 60906, 14526, 47082, 56836, 54193, 18982, 12512, 24832, 36574, 1, 45703, 38690, 45484, 35415, 314, 55924, 46038, 35577, 35415, 20509, 38690, 18711, 28780, 14088, 28424, 4846, 23686, 49367, 4182, 38993, 58703, 43273, 63549, 18819, 58592, 1, 43387, 20667, 5327, 49263, 59582, 39231, 4781, 52025, 56135, 54577, 49700, 33751, 52626, 42307, 30132, 46531, 65088, 12820, 24071, 61570, 5447, 21877, 23998, 49746, 27994, 50181, 1, 50848, 9027, 16548, 41692, 22400, 35634, 29564, 54162, 41692, 59463, 9027, 15915, 56292, 27759, 55739, 41754, 13965, 59517, 21692, 31121, 45814, 50521, 22047, 15895, 23684, 1, 22134, 49868, 3371, 47709, 19029, 2506, 21233, 60242, 61574, 24034, 28745, 40468, 5739, 34419, 9774, 53989, 21513, 1072, 26282, 64594, 24083, 2667, 50459, 61399, 30344, 19004, 1, 55638, 40156, 59604, 13846, 27598, 9158, 47366, 21891, 57227, 33147, 869, 12502, 59782, 46030, 43325, 24573, 64346, 58405, 63522, 65212, 13443, 20341, 12120, 31846, 15862, 19484, 1, 51894, 31612, 28812, 17601, 28588, 59118, 41450, 920, 11619, 2369, 65101, 49117, 47554, 36567, 12541, 11919, 39821, 22221, 23121, 37385, 12325, 25909, 19473, 29734, 16743, 62737, 1, 5754, 40332, 57516, 55342, 8881, 36632, 23479, 27529, 50629, 54909, 49327, 37298, 12784, 59245, 27380, 3677, 13339, 24517, 31378, 63501, 9393, 54800, 50602, 64408, 37457, 17184, 1, 17757, 53537, 41305, 42142, 60700, 23373, 38905, 3344, 21132, 35342, 41178, 22799, 4809, 51772, 64846, 46985, 16484, 44632, 42404, 15511, 30639, 61149, 54722, 39159, 54915, 53273, 1, 38294, 62754, 1016, 50746, 18905, 42921, 37390, 17554, 24619, 44727, 549, 51240, 10271, 42112, 59137, 12893, 37513, 306, 62955, 26866, 56029, 8563, 27853, 30810, 16438, 3790, 1, 15224, 18759, 43449, 8764, 6701, 46506, 60380, 2791, 21511, 52266, 43369, 57032, 27991, 30054, 36909, 63369, 36747, 65417, 1, 59440, 21807, 54986, 46531, 10698, 61901, 24690, 1, 9305, 8470, 45164, 60560, 50229, 21163, 10446, 13555, 49091, 38850, 59954, 21262, 32682, 43013, 45510, 15317, 62928, 15615, 46063, 38495, 22003, 37463, 29538, 10336, 28086, 36674, 1, 54855, 48977, 61167, 37369, 5442, 22423, 17563, 7910, 60220, 6964, 10596, 59643, 23000, 54927, 57759, 18791, 3990, 12769, 34611, 65374, 19808, 48544, 54201, 19012, 23000, 22857, 1, 1524, 3699, 11319, 33888, 18613, 59642, 13317, 31851, 39978, 60517, 51083, 2177, 64684, 14438, 26742, 3555, 55147, 12660, 1738, 39686, 20571, 24490, 16326, 42305, 64684, 6096, 1, 7619, 9142, 47236, 14984, 22856, 22095, 43175, 37586, 45457, 35808, 7367, 34284, 507, 58156, 63490, 50536, 60442, 16506, 27427, 29966, 37332, 33523, 11427, 38348, 507, 30475, 1, 53331, 4715, 31127, 364, 10231, 43897, 55728, 16614, 60803, 2614, 53117, 38747, 44115, 34248, 54422, 15236, 42445, 26628, 48469, 16323, 32001, 31855, 6675, 45350, 44115, 38602};
    int32_t texp[6966] = {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 2, 0, 1, 2, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 2, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 2, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 1, 1, 1, 0, 0, 0, 0, 1, 2, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 3, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 1, 0, 1, 2, 0, 0, 0, 0, 2, 2, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 2, 0, 2, 0, 0, 0, 0, 2, 2, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 1, 1, 0, 2, 0, 0, 0, 0, 2, 2, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 3, 1, 0, 0, 0, 0, 2, 2, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 1, 2, 1, 0, 0, 0, 0, 2, 2, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 4, 0, 0, 0, 0, 0, 2, 2, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 1, 3, 0, 0, 0, 0, 0, 2, 2, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 5, 0, 0, 0, 0, 1, 3, 1, 0, 0, 0, 0, 4, 1, 0, 0, 1, 2, 0, 2, 0, 0, 0, 3, 0, 2, 0, 1, 1, 0, 1, 2, 0, 0, 2, 0, 1, 2, 0, 1, 0, 1, 1, 2, 0, 0, 1, 1, 1, 2, 0, 0, 0, 2, 1, 2, 0, 1, 0, 0, 2, 2, 0, 0, 1, 0, 2, 2, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 3, 2, 0, 1, 1, 0, 0, 3, 0, 0, 2, 0, 0, 3, 0, 1, 0, 1, 0, 3, 0, 0, 1, 1, 0, 3, 0, 0, 0, 2, 0, 3, 0, 1, 0, 0, 1, 3, 0, 0, 1, 0, 1, 3, 0, 0, 0, 1, 1, 3, 0, 0, 0, 0, 2, 3, 0, 1, 0, 0, 0, 4, 0, 0, 1, 0, 0, 4, 0, 0, 0, 1, 0, 4, 0, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 5, 0, 0, 1, 2, 0, 3, 0, 0, 0, 3, 0, 3, 0, 1, 1, 0, 1, 3, 0, 0, 2, 0, 1, 3, 0, 1, 0, 1, 1, 3, 0, 0, 1, 1, 1, 3, 0, 0, 0, 2, 1, 3, 0, 1, 0, 0, 2, 3, 0, 0, 1, 0, 2, 3, 0, 0, 0, 1, 2, 3, 0, 0, 0, 0, 3, 3, 0, 1, 1, 0, 0, 4, 0, 0, 2, 0, 0, 4, 0, 1, 0, 1, 0, 4, 0, 0, 1, 1, 0, 4, 0, 0, 0, 2, 0, 4, 0, 1, 0, 0, 1, 4, 0, 0, 1, 0, 1, 4, 0, 0, 0, 1, 1, 4, 0, 1, 0, 0, 0, 5, 0, 0, 1, 0, 0, 5, 0, 0, 0, 1, 0, 5, 0, 0, 0, 0, 1, 5, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 4, 0, 1, 0, 1, 1, 4, 0, 0, 1, 1, 1, 4, 0, 0, 0, 2, 1, 4, 0, 1, 0, 0, 2, 4, 0, 0, 1, 0, 2, 4, 0, 0, 0, 1, 2, 4, 0, 0, 0, 0, 3, 4, 0, 1, 1, 0, 0, 5, 0, 0, 2, 0, 0, 5, 0, 1, 0, 1, 0, 5, 0, 0, 1, 1, 0, 5, 0, 0, 0, 2, 0, 5, 0, 1, 0, 0, 1, 5, 0, 0, 1, 0, 1, 5, 0, 0, 0, 1, 1, 5, 0, 0, 0, 0, 2, 5, 0, 1, 0, 0, 0, 6, 0, 0, 1, 0, 0, 6, 0, 0, 0, 1, 0, 6, 0, 0, 0, 0, 1, 6, 0, 0, 0, 0, 0, 7, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 4, 0, 1, 0, 1, 1, 4, 0, 0, 1, 1, 1, 4, 0, 0, 0, 2, 1, 4, 0, 1, 0, 0, 2, 4, 0, 0, 1, 0, 2, 4, 0, 0, 0, 1, 2, 4, 0, 0, 0, 0, 3, 4, 0, 1, 1, 0, 0, 5, 0, 0, 2, 0, 0, 5, 0, 1, 0, 1, 0, 5, 0, 0, 1, 1, 0, 5, 0, 0, 0, 2, 0, 5, 0, 1, 0, 0, 1, 5, 0, 0, 1, 0, 1, 5, 0, 0, 0, 1, 1, 5, 0, 0, 0, 0, 2, 5, 0, 1, 0, 0, 0, 6, 0, 0, 1, 0, 0, 6, 0, 0, 0, 1, 0, 6, 0, 0, 0, 0, 1, 6, 0, 0, 0, 0, 0, 7, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 4, 0, 1, 0, 1, 1, 4, 0, 0, 1, 1, 1, 4, 0, 0, 0, 2, 1, 4, 0, 1, 0, 0, 2, 4, 0, 0, 1, 0, 2, 4, 0, 0, 0, 1, 2, 4, 0, 0, 0, 0, 3, 4, 0, 1, 1, 0, 0, 5, 0, 0, 2, 0, 0, 5, 0, 1, 0, 1, 0, 5, 0, 0, 1, 1, 0, 5, 0, 0, 0, 2, 0, 5, 0, 1, 0, 0, 1, 5, 0, 0, 1, 0, 1, 5, 0, 0, 0, 1, 1, 5, 0, 0, 0, 0, 2, 5, 0, 1, 0, 0, 0, 6, 0, 0, 1, 0, 0, 6, 0, 0, 0, 1, 0, 6, 0, 0, 0, 0, 1, 6, 0, 0, 0, 0, 0, 7, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 4, 3, 0, 1, 0, 1, 1, 4, 0, 0, 1, 1, 1, 4, 0, 0, 0, 2, 1, 4, 0, 1, 0, 0, 2, 4, 0, 0, 1, 0, 2, 4, 0, 0, 0, 1, 2, 4, 0, 0, 0, 0, 3, 4, 0, 1, 1, 0, 0, 5, 0, 0, 2, 0, 0, 5, 0, 1, 0, 1, 0, 5, 0, 0, 1, 1, 0, 5, 0, 0, 0, 2, 0, 5, 0, 1, 0, 0, 1, 5, 0, 0, 1, 0, 1, 5, 0, 0, 0, 1, 1, 5, 0, 0, 0, 0, 2, 5, 0, 1, 0, 0, 0, 6, 0, 0, 1, 0, 0, 6, 0, 0, 0, 1, 0, 6, 0, 0, 0, 0, 1, 6, 0, 0, 0, 0, 0, 7, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 6, 0, 0, 1, 1, 0, 6, 0, 0, 0, 2, 0, 6, 0, 1, 0, 0, 1, 6, 0, 0, 1, 0, 1, 6, 0, 0, 0, 1, 1, 6, 0, 0, 0, 0, 2, 6, 0, 1, 0, 0, 0, 7, 0, 0, 1, 0, 0, 7, 0, 0, 0, 1, 0, 7, 0, 0, 0, 0, 1, 7, 0, 0, 0, 0, 0, 8, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 2, 0, 0, 6, 0, 0, 1, 1, 0, 6, 0, 0, 0, 2, 0, 6, 0, 1, 0, 0, 1, 6, 0, 0, 1, 0, 1, 6, 0, 0, 0, 1, 1, 6, 0, 0, 0, 0, 2, 6, 0, 1, 0, 0, 0, 7, 0, 0, 1, 0, 0, 7, 0, 0, 0, 1, 0, 7, 0, 0, 0, 0, 1, 7, 0, 0, 0, 0, 0, 8, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 1, 1, 0, 0, 6, 0, 0, 1, 1, 0, 6, 0, 0, 0, 2, 0, 6, 0, 1, 0, 0, 1, 6, 0, 0, 1, 0, 1, 6, 0, 0, 0, 1, 1, 6, 0, 0, 0, 0, 2, 6, 0, 1, 0, 0, 0, 7, 0, 0, 1, 0, 0, 7, 0, 0, 0, 1, 0, 7, 0, 0, 0, 0, 1, 7, 0, 0, 0, 0, 0, 8, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 2, 5, 0, 0, 0, 0, 3, 5, 0, 0, 1, 1, 0, 6, 0, 0, 0, 2, 0, 6, 0, 1, 0, 0, 1, 6, 0, 0, 1, 0, 1, 6, 0, 0, 0, 1, 1, 6, 0, 0, 0, 0, 2, 6, 0, 1, 0, 0, 0, 7, 0, 0, 1, 0, 0, 7, 0, 0, 0, 1, 0, 7, 0, 0, 0, 0, 1, 7, 0, 0, 0, 0, 0, 8, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 2, 5, 0, 0, 0, 0, 3, 5, 0, 0, 1, 1, 0, 6, 0, 0, 0, 2, 0, 6, 0, 1, 0, 0, 1, 6, 0, 0, 1, 0, 1, 6, 0, 0, 0, 1, 1, 6, 0, 0, 0, 0, 2, 6, 0, 1, 0, 0, 0, 7, 0, 0, 1, 0, 0, 7, 0, 0, 0, 1, 0, 7, 0, 0, 0, 0, 1, 7, 0, 0, 0, 0, 0, 8, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 2, 5, 0, 0, 0, 0, 3, 5, 0, 0, 1, 1, 0, 6, 0, 0, 0, 2, 0, 6, 0, 1, 0, 0, 1, 6, 0, 0, 1, 0, 1, 6, 0, 0, 0, 1, 1, 6, 0, 0, 0, 0, 2, 6, 0, 1, 0, 0, 0, 7, 0, 0, 1, 0, 0, 7, 0, 0, 0, 1, 0, 7, 0, 0, 0, 0, 1, 7, 0, 0, 0, 0, 0, 8, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 2, 1, 5, 0, 0, 0, 0, 3, 5, 0, 0, 1, 1, 0, 6, 0, 0, 0, 2, 0, 6, 0, 1, 0, 0, 1, 6, 0, 0, 1, 0, 1, 6, 0, 0, 0, 1, 1, 6, 0, 0, 0, 0, 2, 6, 0, 1, 0, 0, 0, 7, 0, 0, 1, 0, 0, 7, 0, 0, 0, 1, 0, 7, 0, 0, 0, 0, 1, 7, 0, 0, 0, 0, 0, 8, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 1, 5, 0, 0, 0, 0, 3, 5, 0, 0, 1, 1, 0, 6, 0, 0, 0, 2, 0, 6, 0, 1, 0, 0, 1, 6, 0, 0, 1, 0, 1, 6, 0, 0, 0, 1, 1, 6, 0, 0, 0, 0, 2, 6, 0, 1, 0, 0, 0, 7, 0, 0, 1, 0, 0, 7, 0, 0, 0, 1, 0, 7, 0, 0, 0, 0, 1, 7, 0, 0, 0, 0, 0, 8, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 1, 0, 1, 1, 5, 0, 0, 0, 0, 3, 5, 0, 0, 1, 1, 0, 6, 0, 0, 0, 2, 0, 6, 0, 1, 0, 0, 1, 6, 0, 0, 1, 0, 1, 6, 0, 0, 0, 1, 1, 6, 0, 0, 0, 0, 2, 6, 0, 1, 0, 0, 0, 7, 0, 0, 1, 0, 0, 7, 0, 0, 0, 1, 0, 7, 0, 0, 0, 0, 1, 7, 0, 0, 0, 0, 0, 8, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 9, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 8, 0, 0, 0, 0, 1, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 0, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 8, 0, 0, 0, 0, 1, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 2, 7, 0, 0, 0, 0, 1, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 7, 0, 0, 0, 0, 1, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 0, 1, 7, 0, 0, 0, 0, 1, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 1, 7, 0, 0, 0, 0, 1, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2, 0, 7, 0, 0, 0, 0, 1, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 7, 0, 0, 0, 0, 1, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 3, 6, 0, 0, 0, 0, 1, 8, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 1, 1, 0, 0, 1, 0, 0, 2, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 2, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 2, 1, 0, 1, 0, 0, 0, 2, 0, 0, 1, 0, 0, 2, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 1, 0, 0, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 4, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 3, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 2, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 2, 1, 1, 0, 1, 0, 0, 2, 1, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, 2, 1, 0, 0, 0, 0, 3, 1, 0, 1, 1, 0, 0, 2, 0, 0, 2, 0, 0, 2, 0, 1, 0, 1, 0, 2, 0, 0, 1, 1, 0, 2, 0, 0, 0, 2, 0, 2, 0, 1, 0, 0, 1, 2, 0, 0, 1, 0, 1, 2, 0, 0, 0, 1, 1, 2, 0, 0, 0, 0, 2, 2, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0, 3, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 4};

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
