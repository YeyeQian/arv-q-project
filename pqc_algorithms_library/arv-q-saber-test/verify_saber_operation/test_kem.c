#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "api.h"
#include "../apis/custom_inst_api.h"
#include <stdbool.h>
#include <ctype.h>
#include "../common/debug.h"

#define	MAX_MARKER_LEN		50
#define KAT_SUCCESS          0
#define KAT_FILE_OPEN_ERROR -1
#define KAT_DATA_ERROR      -3
#define KAT_CRYPTO_FAILURE  -4

bool test_crypto_kem(){
    /********************************
     * CRYPTO_KEM_KEYPAIR
    ********************************/
    unsigned char       pk_ref[CRYPTO_PUBLICKEYBYTES], sk_ref[CRYPTO_SECRETKEYBYTES];
    unsigned char       pk_res[CRYPTO_PUBLICKEYBYTES], sk_res[CRYPTO_SECRETKEYBYTES];
    unsigned char       ct_ref[CRYPTO_CIPHERTEXTBYTES], ss_ref[CRYPTO_BYTES];
    unsigned char       ct_res[CRYPTO_CIPHERTEXTBYTES], ss_res[CRYPTO_BYTES];
    unsigned char       ss1_ref[CRYPTO_BYTES],ss1_res[CRYPTO_BYTES];

    int i;
    bool flag_pk, flag_sk;
    uint64_t start, end;

    //cache warm up
    for(int i=0;i<3;i++){
        crypto_kem_keypair_bf_custom(pk_res,sk_res);
        crypto_kem_enc_bf_custom(ct_res,ss_res,pk_res);
        crypto_kem_dec_bf_custom(ss1_res,ct_res,sk_res);
    }

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    crypto_kem_keypair(pk_ref,sk_ref);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_keypair takes %d cycles\n",end-start);

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    crypto_kem_keypair_bf_custom(pk_res,sk_res);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_keypair_bf_custom takes %d cycles\n",end-start);

    //check pk
    flag_pk = true;
    for(i=0;i<CRYPTO_PUBLICKEYBYTES;i++){
        if(pk_res[i] != pk_ref[i]) {
            log_e("pk_ref[%d]=%d, while pk_res[%d]=%d", i, pk_ref[i], i, pk_res[i]);
            flag_pk = false;
        }
    }
    if(flag_pk) {
        printf("pk of test_crypto_kem_keypair test pass!\n\n");
    }
    else {
        printf("pk of test_crypto_kem_keypair test fail..\n\n");
    }

    //check sk
    flag_sk = true;
    for(i = 0; i < CRYPTO_SECRETKEYBYTES; i++) {
        if(sk_res[i] != sk_ref[i]) {
            log_e("sk_ref[%d]=%d, while sk_res[%d]=%d", i, sk_ref[i], i, sk_res[i]);
            flag_sk = false;
        }
    }  
    if(flag_sk) {
        printf("sk of test_crypto_kem_keypair test pass!\n\n");
    }
    else {
        printf("sk of test_crypto_kem_keypair test fail..\n\n");
    }
    /********************************
     * CRYPTO_KEM_ENC
    ********************************/

   __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    crypto_kem_enc(ct_ref,ss_ref,pk_ref);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_enc takes %d cycles\n",end-start);

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    crypto_kem_enc_bf_custom(ct_res,ss_res,pk_res);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_enc_bf_custom takes %d cycles\n",end-start);

    //check ct
    bool flag_ct=true;
    for(i=0;i<CRYPTO_CIPHERTEXTBYTES;i++){
        if(ct_res[i] != ct_ref[i]) {
            log_e("ct_ref[%d]=%d, while ct_res[%d]=%d", i, ct_ref[i], i, ct_res[i]);
            flag_ct = false;
        }
    }
    if(flag_ct) {
        printf("ct of test_crypto_kem_enc test pass!\n\n");
    }
    else {
        printf("ct of test_crypto_kem_enc test fail..\n\n");
    }
    //check ss
    bool flag_ss=true;
    for(i=0;i<CRYPTO_BYTES;i++){
        if(ss_res[i] != ss_ref[i]) {
            log_e("ss_ref[%d]=%d, while ss_res[%d]=%d", i, ss_ref[i], i, ss_res[i]);
            flag_ss = false;
        }
    }
    if(flag_ss) {
        printf("ss of test_crypto_kem_enc test pass!\n\n");
    }
    else {
        printf("ss of test_crypto_kem_enc test fail..\n\n");
    }
    /********************************
     * CRYPTO_KEM_DEC
    ********************************/

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    crypto_kem_dec(ss1_ref,ct_ref,sk_ref);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_dec takes %d cycles\n",end-start);

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    crypto_kem_dec_bf_custom(ss1_res,ct_res,sk_res);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_dec_bf_custom takes %d cycles\n",end-start);

    //check ss1
    bool flag_ss1=true;
    for(i=0;i<CRYPTO_BYTES;i++){
        if(ss1_res[i] != ss1_ref[i]) {
            log_e("ss1_ref[%d]=%d, while ss1_res[%d]=%d", i, ss1_ref[i], i, ss1_res[i]);
            flag_ss1 = false;
        }
    }
    if(flag_ss1) {
        printf("ss1 of test_crypto_kem_dec test pass!\n\n");
    }
    else {
        printf("ss1 of test_crypto_kem_dec test fail..\n\n");
    }

    if ( memcmp(ss_ref, ss1_ref, CRYPTO_BYTES) ) {
        printf("crypto_kem_dec returned bad 'ss' value\n");
        return false;
    }

    return (flag_pk && flag_sk && flag_ct && flag_ss && flag_ss1);
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

#if   (SABER_L == 2)
    security_level = 512;
#elif (SABER_L == 3)
    security_level = 768;
#elif (SABER_L == 4)
    security_level = 1024;
#else
#endif

    printf("start saber%d.CCAKEM cache warming\n", security_level);
    for(int i=0;i<3;i++){
        ret_val = crypto_kem_keypair_bf_custom(pk, sk);
        printf("warmup keypair fnished\n");
        ret_val = crypto_kem_enc_bf_custom(ct, ss, pk);
        printf("warmup enc fnished\n");
        ret_val = crypto_kem_dec_bf_custom(ss1, ct, sk);
        printf("warmup dec fnished\n");
    }
    printf("end saber%d.CCAKEM cache warming\n", security_level);

    start = read_cycle();
    ret_val = crypto_kem_keypair_bf_custom(pk, sk);
    end = read_cycle();
    printf("saber%d.CCAKEM.KeyGen() execution took %lu cycles in RVV\n", security_level, end - start);

    start = read_cycle();
    ret_val = crypto_kem_enc_bf_custom(ct, ss, pk);
    end = read_cycle();
    printf("saber%d.CCAKEM.Enc() execution took %lu cycles in RVV\n", security_level, end - start);

    start = read_cycle();
    ret_val = crypto_kem_dec_bf_custom(ss1, ct, sk);
    end = read_cycle();
    printf("saber%d.CCAKEM.Dec() execution took %lu cycles in RVV\n", security_level, end - start);

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


int main(){
    //test_crypto_kem();
    test_kem_custom();
    return 0;
}