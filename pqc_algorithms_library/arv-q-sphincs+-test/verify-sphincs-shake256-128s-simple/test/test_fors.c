#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../hash.h"
#include "../fors.h"
// #include "../rng.h"
#include "../params.h"

int main()
{
    /* Make stdout buffer more responsive. */
    setbuf(stdout, NULL);

    unsigned char sk_seed[SPX_N];
    unsigned char pub_seed[SPX_N];
    unsigned char pk1[SPX_FORS_PK_BYTES];
    unsigned char pk2[SPX_FORS_PK_BYTES];
    unsigned char sig[SPX_FORS_BYTES];
    unsigned char m[SPX_FORS_MSG_BYTES];
    uint32_t addr[8] = {0};

    // randombytes(sk_seed, SPX_N);
    // randombytes(pub_seed, SPX_N);
    // randombytes(m, SPX_FORS_MSG_BYTES);
    // randombytes((unsigned char *)addr, 8 * sizeof(uint32_t));    
    srand((unsigned)time(NULL));
    for(int i_temp = 0; i_temp < SPX_N; i_temp++) {
        sk_seed[i_temp] = rand() % 0xff;
    }
    for(int i_temp = 0; i_temp < SPX_N; i_temp++) {
        pub_seed[i_temp] = rand() % 0xff;
    }
    for(int i_temp = 0; i_temp < SPX_FORS_MSG_BYTES; i_temp++) {
        m[i_temp] = rand() % 0xff;
    }
    for(int i_temp = 0; i_temp < 8 * sizeof(uint32_t); i_temp++) {
        ((unsigned char *)addr)[i_temp] = rand() % 0xff;
    }
    printf("Testing FORS signature and PK derivation.. ");

    initialize_hash_function(pub_seed, sk_seed);

    fors_sign(sig, pk1, m, sk_seed, pub_seed, addr);
    fors_pk_from_sig(pk2, sig, m, pub_seed, addr);

    if (memcmp(pk1, pk2, SPX_FORS_PK_BYTES)) {
        printf("failed!\n");
        return -1;
    }
    printf("successful.\n");
    return 0;
}
