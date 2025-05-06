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
#include "api.h"
#include "kem.h"

#define	MAX_MARKER_LEN		50
#define KAT_SUCCESS          0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_DATA_ERROR      -3
#define KAT_CRYPTO_FAILURE  -4

#define CUSTOM

#ifdef CUSTOM
int
main()
{
    uint32_t security_level = 0;

    #ifdef PARAM64
        security_level = 1;
    #elif defined(PARAM96)
        security_level = 3;
    #elif defined(PARAM128)
        security_level = 5;
    #else
    #endif

    printf("PQCgenKAT_kem custom test for security level-%d:\n",security_level);
    unsigned char       ct[CRYPTO_CIPHERTEXTBYTES]={0}, ss[CRYPTO_BYTES]={0}, ss1[CRYPTO_BYTES]={0};
    unsigned char       pk[CRYPTO_PUBLICKEYBYTES]={0}, sk[CRYPTO_SECRETKEYBYTES]={0};
    int                 ret_val;
    uint64_t start,end;

    //cache warm up
    printf("start cache warmup\n");
    for(int i=0;i<1;i++){
        ret_val = crypto_kem_keypair_custom(pk, sk);
        ret_val = crypto_kem_enc_custom(ct, ss, pk);
        //ret_val = crypto_kem_dec_custom(ss1, ct, sk);
    }
    printf("end cache warmup\n");

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    ret_val = crypto_kem_keypair_custom(pk, sk);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_keypair_custom cost %lu cycles\n",end-start);
    if (ret_val != 0) {
        printf("crypto_kem_keypair returned <%d>\n", ret_val);
        //return KAT_CRYPTO_FAILURE;
    }

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    ret_val = crypto_kem_enc_custom(ct, ss, pk);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_enc_custom cost %lu cycles\n",end-start);
    if (ret_val != 0) {
        printf("crypto_kem_enc returned <%d>\n", ret_val);
        //return KAT_CRYPTO_FAILURE;
    }

    // __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    // ret_val = crypto_kem_dec_custom(ss1, ct, sk);
    // __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    // printf("crypto_kem_dec_custom cost %lu cycles\n",end-start);
    // if (ret_val != 0) {
    //     printf("crypto_kem_dec returned <%d>\n", ret_val);
    //     //return KAT_CRYPTO_FAILURE;
    // }

    // if ( memcmp(ss, ss1, CRYPTO_BYTES) ) {
    //     printf("crypto_kem_dec returned bad 'ss' value\n");
    //     //return KAT_CRYPTO_FAILURE;
    // }

    return KAT_SUCCESS;
}
#else
int
main()
{
    uint32_t security_level = 0;

    #ifdef PARAM64
        security_level = 1;
    #elif defined(PARAM96)
        security_level = 3;
    #elif defined(PARAM128)
        security_level = 5;
    #else
    #endif

    printf("PQCgenKAT_kem reference test for security level-%d:\n",security_level);
    unsigned char       ct[CRYPTO_CIPHERTEXTBYTES]={0}, ss[CRYPTO_BYTES]={0}, ss1[CRYPTO_BYTES]={0};
    unsigned char       pk[CRYPTO_PUBLICKEYBYTES]={0}, sk[CRYPTO_SECRETKEYBYTES]={0};
    int                 ret_val;
    uint64_t start,end;

    //cache warm up
    printf("start cache warmup\n");
    for(int i=0;i<3;i++){
        //ret_val = crypto_kem_keypair_custom(pk, sk);
        //ret_val = crypto_kem_enc_custom(ct, ss, pk);
        ret_val = crypto_kem_dec_custom(ss1, ct, sk);
    }
    printf("end cache warmup\n");

    // __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    // ret_val = crypto_kem_keypair(pk, sk);
    // __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    // printf("crypto_kem_keypair cost %lu cycles\n",end-start);
    // if (ret_val != 0) {
    //     printf("crypto_kem_keypair returned <%d>\n", ret_val);
    //     return KAT_CRYPTO_FAILURE;
    // }

    // __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    // ret_val = crypto_kem_enc(ct, ss, pk);
    // __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    // printf("crypto_kem_enc cost %lu cycles\n",end-start);
    // if (ret_val != 0) {
    //     printf("crypto_kem_enc returned <%d>\n", ret_val);
    //     return KAT_CRYPTO_FAILURE;
    // }

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    ret_val = crypto_kem_dec(ss1, ct, sk);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_dec cost %lu cycles\n",end-start);
    if (ret_val != 0) {
        printf("crypto_kem_dec returned <%d>\n", ret_val);
        //return KAT_CRYPTO_FAILURE;
    }

    if ( memcmp(ss, ss1, CRYPTO_BYTES) ) {
        printf("crypto_kem_dec returned bad 'ss' value\n");
        //return KAT_CRYPTO_FAILURE;
    }

    return KAT_SUCCESS;
}
#endif


