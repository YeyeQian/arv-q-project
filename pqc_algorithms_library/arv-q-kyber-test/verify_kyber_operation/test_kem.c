
//
//  PQCgenKAT_kem.c
//
//  Created by Bassham, Lawrence E (Fed) on 8/29/17.
//  Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
//
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "api.h"
#include "../apis/custom_inst_api.h"
#include <stdbool.h>

#define	MAX_MARKER_LEN		50
#define KAT_SUCCESS          0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_DATA_ERROR      -3
#define KAT_CRYPTO_FAILURE  -4

bool test_kem()
{
    unsigned char       pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    unsigned char       ct[CRYPTO_CIPHERTEXTBYTES], ss[CRYPTO_BYTES],ss1[CRYPTO_BYTES];
    int                 ret_val;
    uint32_t security_level = 0;
    unsigned long start, end;
    uint32_t i;
    bool flag;

#if   (KYBER_K == 2)
    security_level = 512;
#elif (KYBER_K == 3)
    security_level = 768;
#elif (KYBER_K == 4)
    security_level = 1024;
#else
#error "KYBER_K must be in {2,3,4}"
#endif

    start = read_cycle();
    ret_val = crypto_kem_keypair(pk, sk);
    end = read_cycle();
    printf("kyber%d.CCAKEM.KeyGen() execution took %lu cycles in pure software\n", security_level, end - start);

    start = read_cycle();
    ret_val = crypto_kem_enc(ct, ss, pk);
    end = read_cycle();
    printf("kyber%d.CCAKEM.Enc() execution took %lu cycles in pure software\n", security_level, end - start);

    start = read_cycle();
    ret_val = crypto_kem_dec(ss1, ct, sk);
    end = read_cycle();
    printf("kyber%d.CCAKEM.Dec() execution took %lu cycles in pure software\n", security_level, end - start);

    flag = true;
    for(i = 0; i < CRYPTO_BYTES; i++) {
        if(ss[i] != ss1[i]) {
            flag = false;
        }
    }

    if(flag) {
        printf("kem test pass with flag=%d\n", flag);
    }
    else {
        printf("kem test fail with flag=%d\n", flag);
    }

    return flag;
}

bool test_kem_custom()
{
    unsigned char       pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    unsigned char       ct[CRYPTO_CIPHERTEXTBYTES], ss[CRYPTO_BYTES],ss1[CRYPTO_BYTES];
    int                 ret_val;
    uint32_t security_level = 0;
    unsigned long start, end;
    uint32_t i;
    bool flag;

#if   (KYBER_K == 2)
    security_level = 512;
#elif (KYBER_K == 3)
    security_level = 768;
#elif (KYBER_K == 4)
    security_level = 1024;
#else
#error "KYBER_K must be in {2,3,4}"
#endif

    printf("start kyber%d.CCAKEM cache warming\n", security_level);
    for(int i=0;i<3;i++){
        ret_val = crypto_kem_keypair_custom(pk, sk);
        printf("warmup keypair fnished\n");
        ret_val = crypto_kem_enc_custom(ct, ss, pk);
        printf("warmup enc fnished\n");
        ret_val = crypto_kem_dec_custom(ss1, ct, sk);
        printf("warmup dec fnished\n");
    }
    printf("end kyber%d.CCAKEM cache warming\n", security_level);

    start = read_cycle();
    ret_val = crypto_kem_keypair_custom(pk, sk);
    end = read_cycle();
    printf("kyber%d.CCAKEM.KeyGen() execution took %lu cycles in RVV\n", security_level, end - start);

    start = read_cycle();
    ret_val = crypto_kem_enc_custom(ct, ss, pk);
    end = read_cycle();
    printf("kyber%d.CCAKEM.Enc() execution took %lu cycles in RVV\n", security_level, end - start);

    start = read_cycle();
    ret_val = crypto_kem_dec_custom(ss1, ct, sk);
    end = read_cycle();
    printf("kyber%d.CCAKEM.Dec() execution took %lu cycles in RVV\n", security_level, end - start);

    flag = true;
    for(i = 0; i < CRYPTO_BYTES; i++) {
        if(ss[i] != ss1[i]) {
            flag = false;
        }
    }

    if(flag) {
        printf("kem_custom test pass with flag=%d\n", flag);
    }
    else {
        printf("kem_custom test fail with flag=%d\n", flag);
    }

    return flag;
}

bool test_kem_custom_asm()
{
    unsigned char       pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
    unsigned char       ct[CRYPTO_CIPHERTEXTBYTES], ss[CRYPTO_BYTES],ss1[CRYPTO_BYTES];
    int                 ret_val;
    uint32_t security_level = 0;
    unsigned long start, end;
    uint32_t i;
    bool flag;

#if   (KYBER_K == 2)
    security_level = 512;
#elif (KYBER_K == 3)
    security_level = 768;
#elif (KYBER_K == 4)
    security_level = 1024;
#else
#error "KYBER_K must be in {2,3,4}"
#endif

    printf("start kyber%d.CCAKEM cache warming\n", security_level);
    for(int i=0;i<3;i++){
        ret_val = crypto_kem_keypair_custom_asm(pk, sk);
        ret_val = crypto_kem_enc_custom_asm(ct, ss, pk);
        ret_val = crypto_kem_dec_custom_asm(ss1, ct, sk);
    }
    printf("end kyber%d.CCAKEM cache warming\n", security_level);

    start = read_cycle();
    ret_val = crypto_kem_keypair_custom_asm(pk, sk);
    //ret_val = crypto_kem_keypair(pk, sk);
    end = read_cycle();
    printf("kyber%d.CCAKEM.KeyGen() execution took %lu cycles in RVV asm optimized\n", security_level, end - start);

    start = read_cycle();
    ret_val = crypto_kem_enc_custom_asm(ct, ss, pk);
    //ret_val = crypto_kem_enc(ct, ss, pk);
    end = read_cycle();
    printf("kyber%d.CCAKEM.Enc() execution took %lu cycles in RVV asm optimized\n", security_level, end - start);

    start = read_cycle();
    ret_val = crypto_kem_dec_custom_asm(ss1, ct, sk);
    //ret_val = crypto_kem_dec(ss1, ct, sk);
    end = read_cycle();
    printf("kyber%d.CCAKEM.Dec() execution took %lu cycles in RVV asm optimized\n", security_level, end - start);

    flag = true;
    for(i = 0; i < CRYPTO_BYTES; i++) {
        if(ss[i] != ss1[i]) {
            flag = false;
        }
    }

    if(flag) {
        printf("kem_custom_asm test pass with flag=%d\n", flag);
    }
    else {
        printf("kem_custom_asm test fail with flag=%d\n", flag);
    }

    return flag;
}

int main()
{
    //test_kem_custom_asm();
    test_kem_custom();

    return 0;
}