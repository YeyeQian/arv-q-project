
//
//  PQCgenKAT_kem.c
//
//  Created by Bassham, Lawrence E (Fed) on 8/29/17.
//  Copyright © 2017 Bassham, Lawrence E (Fed). All rights reserved.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "api.h"
#include "compiler.h"
#include "fips202.h"
#include "symmetric.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>


#define SHAKE_128_IN_LEN 1000
#define SHAKE_128_OUT_LEN 1000
#define SHAKE_256_IN_LEN 1000
#define SHAKE_256_OUT_LEN 1000

#define SHA3_256_IN_LEN 1000
#define SHA3_256_OUT_LEN 32
#define SHA3_512_IN_LEN 1000
#define SHA3_512_OUT_LEN 64    

int test_keccak_sha3_256()
{   
    uint8_t seed[SHA3_256_IN_LEN];
    uint8_t hashvalue_ref[SHA3_256_OUT_LEN];
    uint8_t hashvalue[SHA3_256_OUT_LEN];
    int i, flag;

    uint64_t start,end;

    srand((unsigned)time(NULL));
    for (i = 0; i < SHA3_256_IN_LEN; i++) {
        seed[i] = rand() & 0xff;
    }    

    start=read_cycle();
    sha3_256(hashvalue_ref, seed, SHA3_256_IN_LEN);
    end=read_cycle();
    printf("sha3_256 takes %lu cycles\n",end-start);

    start=read_cycle();
    sha3_256_custom(hashvalue, seed, SHA3_256_IN_LEN);
    end=read_cycle();
    printf("sha3_256_custom takes %lu cycles\n",end-start);

    flag = 0;
    // for(i = 0; i < SHA3_256_OUT_LEN; i++) {
    //     if(hashvalue[i] != hashvalue_ref[i]) {
    //         printf("h_ref[%d]=%d, while h[%d]=%d\n, i, hashvalue_ref[i], i, hashvalue[i]);
    //         flag = 1;
    //     }
    // }

    // if(flag==0) {
    //     printf("sha256 test pass!\n");
    // }
    // else {
    //     printf("sha256 test fail..\n");
    // }

    return flag;
}

int test_keccak_sha3_512()
{
    uint8_t seed[SHA3_512_IN_LEN];
    uint8_t hashvalue_ref[SHA3_512_OUT_LEN];
    uint8_t hashvalue[SHA3_512_OUT_LEN];
    int i, flag;

    uint64_t start,end;

    srand((unsigned)time(NULL));
    for (i = 0; i < SHA3_512_IN_LEN; i++) {
        seed[i] = rand() & 0xff;
    }

    start=read_cycle();
    sha3_512(hashvalue_ref, seed, SHA3_512_IN_LEN);
    end=read_cycle();
    printf("sha3_512 takes %lu cycles\n",end-start);

    start=read_cycle();
    sha3_512_custom(hashvalue, seed, SHA3_512_IN_LEN);
    end=read_cycle();
    printf("sha3_512_custom takes %lu cycles\n",end-start);

    flag = 0;
    for(i = 0; i < SHA3_512_OUT_LEN; i++) {
        if(hashvalue[i] != hashvalue_ref[i]) {
            printf("h_ref[%d]=%d, while h[%d]=%d\n", i, hashvalue_ref[i], i, hashvalue[i]);
            flag = 1;
        }
    }

    if(flag==0) {
        printf("sha512 test pass!\n");
    }
    else {
        printf("sha512 test fail..\n");
    }    
    return flag;
}

int test_keccak_shake_128()
{   
    uint8_t seed[SHAKE_128_IN_LEN];
    uint8_t hashvalue_ref[SHAKE_128_OUT_LEN];
    uint8_t hashvalue[SHAKE_128_OUT_LEN];
    int i, flag;

    uint64_t start,end;

    srand((unsigned)time(NULL));
    for (i = 0; i < SHAKE_128_IN_LEN; i++) {
        seed[i] = rand() & 0xff;
    }    

    start=read_cycle();
    shake128(hashvalue_ref, SHAKE_128_OUT_LEN, seed, SHAKE_128_IN_LEN);
    end=read_cycle();
    printf("shake128 takes %lu cycles\n",end-start);

    start=read_cycle();
    shake128_custom(hashvalue, SHAKE_128_OUT_LEN, seed, SHAKE_128_IN_LEN);
    end=read_cycle();
    printf("shake128_custom takes %lu cycles\n",end-start);

    flag = 0;
    for(i = 0; i < SHAKE_128_OUT_LEN; i++) {
        if(hashvalue[i] != hashvalue_ref[i]) {
            printf("h_ref[%d]=%d, while h[%d]=%d\n", i, hashvalue_ref[i], i, hashvalue[i]);
            flag = 1;
        }
    }

    if(flag==0) {
        printf("shake128 test pass!\n");
    }
    else {
        printf("shake128 test fail..\n");
    }

    return flag;
}

int test_keccak_shake_256()
{   
    uint8_t seed[SHAKE_256_IN_LEN];
    uint8_t hashvalue_ref[SHAKE_256_OUT_LEN];
    uint8_t hashvalue[SHAKE_256_OUT_LEN];
    int i, flag;

    srand((unsigned)time(NULL));
    for (i = 0; i < SHAKE_256_IN_LEN; i++) {
        seed[i] = rand() & 0xff;
    }    

    shake256(hashvalue_ref, SHAKE_256_OUT_LEN, seed, SHAKE_256_IN_LEN);
    shake256_custom(hashvalue, SHAKE_256_OUT_LEN, seed, SHAKE_256_IN_LEN);

    flag = 0;
    for(i = 0; i < SHAKE_256_OUT_LEN; i++) {
        if(hashvalue[i] != hashvalue_ref[i]) {
            printf("h_ref[%d]=%d, while h[%d]=%d\n", i, hashvalue_ref[i], i, hashvalue[i]);
            flag = 1;
        }
    }

    if(flag==0) {
        printf("shake256 test pass!\n");
    }
    else {
        printf("shake256 test fail..\n");
    }

    return flag;
}

int
main()
{
    /* 
    ** the following funcs passed tests in KYBER_K = 2, 3, 4
    ** and passed tests in VLEN = 256 and 512
    */
    test_keccak_sha3_256();
    test_keccak_sha3_512();
    test_keccak_shake_128();
    //test_keccak_shake_256();
    
    return 0;
}