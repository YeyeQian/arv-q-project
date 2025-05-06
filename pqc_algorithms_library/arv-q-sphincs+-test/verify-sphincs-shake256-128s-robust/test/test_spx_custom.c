#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "../../apis/custom_inst_api.h"

#include "../api.h"
#include "../params.h"
// #include "../rng.h"

#define SPX_MLEN 32
#define SPX_SIGNATURES 1

bool test_spx_custom()
{
    int ret = 0;
    int i;
    size_t start, end;

    unsigned char pk[SPX_PK_BYTES];
    unsigned char sk[SPX_SK_BYTES];
    unsigned char *m = malloc(SPX_MLEN);
    unsigned char *sm = malloc(SPX_BYTES + SPX_MLEN);
    unsigned char *mout = malloc(SPX_BYTES + SPX_MLEN);
    unsigned long long smlen;
    unsigned long long mlen;

    // randombytes(m, SPX_MLEN);
    srand((unsigned)time(NULL));
    for(int i_temp = 0; i_temp < SPX_MLEN; i_temp++) {
        m[i_temp] = rand() % 0xff;
    }

    printf("Generating keypair.. ");

    if (crypto_sign_keypair_custom(pk, sk)) {
        printf("failed!\n");
        return -1;
    }
    printf("successful.\n");

    printf("Testing %d signatures.. \n", SPX_SIGNATURES);

    for (i = 0; i < SPX_SIGNATURES; i++) {
        printf("  - iteration #%d:\n", i);

        crypto_sign_custom(sm, &smlen, m, SPX_MLEN, sk);

        if (smlen != SPX_BYTES + SPX_MLEN) {
            printf("  X smlen incorrect [%llu != %u]!\n",
                   smlen, SPX_BYTES);
            ret = -1;
        }
        else {
            printf("    smlen as expected [%llu].\n", smlen);
        }

        /* Test if signature is valid. */
        if (crypto_sign_open_custom(mout, &mlen, sm, smlen, pk)) {
            printf("  X verification failed!\n");
            ret = -1;
        }
        else {
            printf("    verification succeeded.\n");
        }

        /* Test if the correct message was recovered. */
        if (mlen != SPX_MLEN) {
            printf("  X mlen incorrect [%llu != %u]!\n", mlen, SPX_MLEN);
            ret = -1;
        }
        else {
            printf("    mlen as expected [%llu].\n", mlen);
        }
        if (memcmp(m, mout, SPX_MLEN)) {
            printf("  X output message incorrect!\n");
            ret = -1;
        }
        else {
            printf("    output message as expected.\n");
        }

        /* Test if signature is valid when validating in-place. */
        if (crypto_sign_open_custom(sm, &mlen, sm, smlen, pk)) {
            printf("  X in-place verification failed!\n");
            ret = -1;
        }
        else {
            printf("    in-place verification succeeded.\n");
        }

        /* Test if flipping bits invalidates the signature (it should). */

        /* Flip the first bit of the message. Should invalidate. */
        sm[smlen - 1] ^= 1;
        if (!crypto_sign_open_custom(mout, &mlen, sm, smlen, pk)) {
            printf("  X flipping a bit of m DID NOT invalidate signature!\n");
            ret = -1;
        }
        else {
            printf("    flipping a bit of m invalidates signature.\n");
        }
        sm[smlen - 1] ^= 1;
    }

    free(m);
    free(sm);
    free(mout);

    bool result = (ret == 0) ? true : false;

    return result;    
}

bool get_spx_keypair_custom()
{
    int i;
    int ret = 0;

    unsigned char pk[SPX_PK_BYTES];
    unsigned char sk[SPX_SK_BYTES];

    printf("Generating keypair.. \n");
    ret = crypto_sign_keypair_custom(pk, sk);
    printf("Finish generating keypair.. \n");

    printf("pk is:\n");
    for(i = 0; i < SPX_PK_BYTES; i++) {
        if(i < (SPX_PK_BYTES -1)) {
            printf("%d, ", pk[i]);
        }
        else {
            printf("%d\n", pk[i]);
        }
    }

    printf("sk is:\n");
    for(i = 0; i < SPX_SK_BYTES; i++) {
        if(i < (SPX_SK_BYTES -1)) {
            printf("%d, ", sk[i]);
        }
        else {
            printf("%d\n", sk[i]);
        }
    }

    bool result = (ret == 0) ? true : false;

    return result;
}

bool spx_keypair_custom_profiling()
{
    int ret = 0;
    size_t start, end;

    unsigned char pk[SPX_PK_BYTES];
    unsigned char sk[SPX_SK_BYTES];

    printf("begin crypto_sign_keypair_custom for warm cache\n");
    ret = crypto_sign_keypair_custom(pk, sk);
    printf("finish crypto_sign_keypair_custom for warm cache\n");

    printf("begin crypto_sign_keypair_custom for getting performance\n");
    start = read_cycle();
    ret = crypto_sign_keypair_custom(pk, sk);
    end = read_cycle();
    printf("finish crypto_sign_keypair_custom for getting performance\n");

    printf("sphincs-shake256-128s-robust sign_keypair takes %ld cycles in HW/SW co-design\n", end - start);
    printf("sphincs-shake256-128s-robust sign_keypair takes %ld cycles in HW/SW co-design\n", end - start);

    bool result = (ret == 0) ? true : false;

    return result;
}

const unsigned char pk_ref[SPX_PK_BYTES] = {
32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 154, 73, 22, 146, 0, 206, 51, 122, 158, 6, 6, 66, 107, 248, 211, 214
};
const unsigned char sk_ref[SPX_SK_BYTES] = {
0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 154, 73, 22, 146, 0, 206, 51, 122, 158, 6, 6, 66, 107, 248, 211, 214    
};


bool spx_sign_and_sign_open_custom_profiling()
{
    int ret = 0;
    int i;
    size_t start, end;

    unsigned char *m = malloc(SPX_MLEN);
    unsigned char *sm = malloc(SPX_BYTES + SPX_MLEN);
    unsigned char *mout = malloc(SPX_BYTES + SPX_MLEN);
    unsigned long long smlen;
    unsigned long long mlen;

    // randombytes(m, SPX_MLEN);
    srand((unsigned)time(NULL));
    for(int i_temp = 0; i_temp < SPX_MLEN; i_temp++) {
        m[i_temp] = rand() % 0xff;
    }

    printf("begin crypto_sign_custom for warm cache\n");
    ret = crypto_sign_custom(sm, &smlen, m, SPX_MLEN, sk_ref);
    printf("finish crypto_sign_custom for warm cache\n");

    printf("begin crypto_sign_open_custom for warm cache\n");
    ret = crypto_sign_open_custom(mout, &mlen, sm, smlen, pk_ref);
    printf("finish crypto_sign_open_custom for warm cache\n");

    printf("begin crypto_sign_custom for getting performance\n");
    start = read_cycle();
    ret = crypto_sign_custom(sm, &smlen, m, SPX_MLEN, sk_ref);
    end = read_cycle();
    printf("finish crypto_sign_custom for getting performance\n");
    printf("sphincs-shake256-128s-robust sign takes %ld cycles in HW/SW co-design\n", end - start);

    printf("begin crypto_sign_open_custom for getting performance\n");
    start = read_cycle();
    ret = crypto_sign_open_custom(mout, &mlen, sm, smlen, pk_ref);
    end = read_cycle();
    printf("finish crypto_sign_open_custom for getting performance\n");
    printf("sphincs-shake256-128s-robust sign_open takes %ld cycles in HW/SW co-design\n", end - start);
    printf("sphincs-shake256-128s-robust sign_open takes %ld cycles in HW/SW co-design\n", end - start);

    free(m);
    free(sm);
    free(mout);

    bool result = (ret == 0) ? true : false;

    return result;
}


int main()
{
    get_spx_keypair_custom();
}