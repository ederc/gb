#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "../src/gb.c"

int main(
        void
        )
{
    int32_t i;

    const int32_t lens[]  = {2,2,2}; 
    const int32_t cfs[]   = {1, 1, 1, 3, 5, 3};
    const int32_t exps[]  = {1, 0, 0, 1, 3, 0, 0, 1, 1, 1, 1, 0};

    const int32_t nr_gens     = 3;

    nvars = 2;
    htes  = 12;
    fc    = 101;

    initialize_basis(nr_gens);
    initialize_global_hash_table();

    if (ndvars != nvars) {
        return 1;
    }
    if (bpv != 16) {
        return 1;
    }

    initialize_local_hash_table();

    import_julia_data(lens, cfs, exps, nr_gens);

    calculate_divmask();
    for (i = 0; i < 16; ++i) {
        if (dm[i] != i) {
            return 1;
        }
    }
    for (i = 16; i < 32; ++i) {
        if (dm[i] + 16 != i) {
            return 1;
        }
    }
    return 0;
}
