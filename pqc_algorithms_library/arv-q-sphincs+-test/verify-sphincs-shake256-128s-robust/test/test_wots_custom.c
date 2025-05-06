#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../hash.h"
#include "../wots.h"
#include "../params.h"
#include "../../apis/custom_inst_api.h"

int test_wots_custom()
{
    /* Make stdout buffer more responsive. */
    setbuf(stdout, NULL);

    unsigned char seed[SPX_N];
    unsigned char pub_seed[SPX_N];
    unsigned char pk1[SPX_WOTS_PK_BYTES];
    unsigned char pk2[SPX_WOTS_PK_BYTES];
    unsigned char sig[SPX_WOTS_BYTES];
    unsigned char m[SPX_N];
    uint32_t addr[8] = {0};

    // randombytes(seed, SPX_N);
    // randombytes(pub_seed, SPX_N);
    // randombytes(m, SPX_N);
    // randombytes((unsigned char *)addr, 8 * sizeof(uint32_t));
    srand((unsigned)time(NULL));
    for(int i_temp = 0; i_temp < SPX_N; i_temp++) {
        seed[i_temp] = rand() % 0xff;
    }
    for(int i_temp = 0; i_temp < SPX_N; i_temp++) {
        pub_seed[i_temp] = rand() % 0xff;
    }
    for(int i_temp = 0; i_temp < SPX_N; i_temp++) {
        m[i_temp] = rand() % 0xff;
    }
    for(int i_temp = 0; i_temp < 8 * sizeof(uint32_t); i_temp++) {
        ((unsigned char *)addr)[i_temp] = rand() % 0xff;
    }

    printf("Testing WOTS signature and PK derivation.. ");

    initialize_hash_function(pub_seed, seed);

    wots_gen_pk_custom(pk1, seed, pub_seed, addr);
    wots_sign_custom(sig, m, seed, pub_seed, addr);
    wots_pk_from_sig_custom(pk2, sig, m, pub_seed, addr);

    if (memcmp(pk1, pk2, SPX_WOTS_PK_BYTES)) {
        printf("test wots_custom failed!\n");
        return -1;
    }
    printf("test wots_custom successful.\n");
    return 0;
}

bool wots_custom_profiling()
{
    size_t start, end;
    
    /* Make stdout buffer more responsive. */
    setbuf(stdout, NULL);

    unsigned char seed[SPX_N];
    unsigned char pub_seed[SPX_N];
    unsigned char pk1[SPX_WOTS_PK_BYTES];
    unsigned char pk2[SPX_WOTS_PK_BYTES];
    unsigned char sig[SPX_WOTS_BYTES];
    unsigned char m[SPX_N];
    uint32_t addr[8] = {0};

    srand((unsigned)time(NULL));
    for(int i_temp = 0; i_temp < SPX_N; i_temp++) {
        seed[i_temp] = rand() % 0xff;
    }
    for(int i_temp = 0; i_temp < SPX_N; i_temp++) {
        pub_seed[i_temp] = rand() % 0xff;
    }
    for(int i_temp = 0; i_temp < SPX_N; i_temp++) {
        m[i_temp] = rand() % 0xff;
    }
    for(int i_temp = 0; i_temp < 8 * sizeof(uint32_t); i_temp++) {
        ((unsigned char *)addr)[i_temp] = rand() % 0xff;
    }

    initialize_hash_function(pub_seed, seed);

    printf("begin wots_gen_pk_custom for warm cache\n");
    wots_gen_pk_custom(pk1, seed, pub_seed, addr);

    printf("begin wots_sign_custom for warm cache\n");
    wots_sign_custom(sig, m, seed, pub_seed, addr);

    printf("begin wots_pk_from_sig_custom for warm cache\n");
    wots_pk_from_sig_custom(pk2, sig, m, pub_seed, addr);

    printf("begin wots_gen_pk_custom for profiling\n");
    start = read_cycle();
    wots_gen_pk_custom(pk1, seed, pub_seed, addr);
    end = read_cycle();
    printf("sphincs-shake256-128s-robust wots_gen_pk_custom takes %ld cycles in HW/SW co-design\n", end - start);

    printf("begin wots_sign_custom for profiling\n");
    start = read_cycle();
    wots_sign_custom(sig, m, seed, pub_seed, addr);
    end = read_cycle();
    printf("sphincs-shake256-128s-robust wots_sign_custom takes %ld cycles in HW/SW co-design\n", end - start);

    printf("begin wots_pk_from_sig_custom for profiling\n");
    start = read_cycle();
    wots_pk_from_sig_custom(pk2, sig, m, pub_seed, addr);
    end = read_cycle();
    printf("sphincs-shake256-128s-robust wots_pk_from_sig_custom takes %ld cycles in HW/SW co-design\n", end - start);

    return true;    
}

int main()
{
    wots_custom_profiling();
    return 0;
}
