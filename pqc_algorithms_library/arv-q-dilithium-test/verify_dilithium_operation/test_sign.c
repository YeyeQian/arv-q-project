#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "../common/debug.h"
#include "param.h"
#include "sign.h"
#include "sign_custom_asm.h"
#include "sign_custom_optimal_c908.h"

// bool test_crypto_sign_keypair(){
//     uint8_t pk[CRYPTO_PUBLICKEYBYTES];
//     uint8_t sk[CRYPTO_SECRETKEYBYTES];
//     uint8_t pk_ref[CRYPTO_PUBLICKEYBYTES];
//     uint8_t sk_ref[CRYPTO_SECRETKEYBYTES];
//     int i;
//     bool flag_pk, flag_sk;
//     uint64_t start, end;

//     crypto_sign_keypair(pk_ref,sk_ref);
//     __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
//     crypto_sign_keypair_custom_asm_fortest(pk,sk);
//     __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
//     printf("crypto_sign_keypair_custom_fortest takes %d cycles\n",end-start);

//     //check pk
//     flag_pk = true;
//     for(i=0;i<CRYPTO_PUBLICKEYBYTES;i++){
//         if(pk[i] != pk_ref[i]) {
//             log_e("pk_ref[%d]=%d, while pk[%d]=%d", i, pk_ref[i], i, pk[i]);
//             flag_pk = false;
//         }
//     }
//     if(flag_pk) {
//         printf("pk of test_crypto_sign_keypair test pass!\n\n");
//     }
//     else {
//         printf("pk of test_crypto_sign_keypair test fail..\n\n");
//     }

//     //check sk
//     flag_sk = true;
//     for(i = 0; i < CRYPTO_SECRETKEYBYTES; i++) {
//         if(sk[i] != sk_ref[i]) {
//             log_e("sk_ref[%d]=%d, while sk[%d]=%d", i, sk_ref[i], i, sk[i]);
//             flag_sk = false;
//         }
//     }  
//     if(flag_sk) {
//         printf("sk of test_crypto_sign_keypair test pass!\n\n");
//     }
//     else {
//         printf("sk of test_crypto_sign_keypair test fail..\n\n");
//     }
//     return (flag_pk && flag_sk); 
// }

// bool test_crypto_sign(){
//     uint8_t pk[CRYPTO_PUBLICKEYBYTES];
//     uint8_t sk[CRYPTO_SECRETKEYBYTES];
//     uint8_t pk_ref[CRYPTO_PUBLICKEYBYTES];
//     uint8_t sk_ref[CRYPTO_SECRETKEYBYTES];
//     int i;
//     bool flag_pk, flag_sk;
//     uint64_t start, end;
//     /********************************
//      * CRYPTO_SIGN_KEYPAIR
//     ********************************/

//     crypto_sign_keypair(pk_ref,sk_ref);
//     crypto_sign_keypair_custom_fortest(pk,sk);

//     //check pk
//     flag_pk = true;
//     for(i=0;i<CRYPTO_PUBLICKEYBYTES;i++){
//         if(pk[i] != pk_ref[i]) {
//             log_e("pk_ref[%d]=%d, while pk[%d]=%d", i, pk_ref[i], i, pk[i]);
//             flag_pk = false;
//         }
//     }
//     if(flag_pk) {
//         printf("pk of test_crypto_sign_keypair test pass!\n\n");
//     }
//     else {
//         printf("pk of test_crypto_sign_keypair test fail..\n\n");
//     }

//     //check sk
//     flag_sk = true;
//     for(i = 0; i < CRYPTO_SECRETKEYBYTES; i++) {
//         if(sk[i] != sk_ref[i]) {
//             log_e("sk_ref[%d]=%d, while sk[%d]=%d", i, sk_ref[i], i, sk[i]);
//             flag_sk = false;
//         }
//     }  
//     if(flag_sk) {
//         printf("sk of test_crypto_sign_keypair test pass!\n\n");
//     }
//     else {
//         printf("sk of test_crypto_sign_keypair test fail..\n\n");
//     }

//     /********************************
//      * CRYPTO_SIGN
//     ********************************/
//    size_t mlen=33;
//    uint8_t* m=(uint8_t*)malloc(mlen*sizeof(uint8_t));
//    size_t smlen_ref;
//    size_t smlen_res;
//    uint8_t* sm_ref=(uint8_t*)malloc((mlen+CRYPTO_BYTES)*sizeof(uint8_t));
//    uint8_t* sm_res=(uint8_t*)malloc((mlen+CRYPTO_BYTES)*sizeof(uint8_t));
//    srand((unsigned)time(NULL));
//    for(i=0;i<mlen;i++){
//     m[i]=i%256;
//     // m[i]=rand()%256;
//    }
//    crypto_sign(sm_ref,&smlen_ref,m,mlen,sk_ref);
//    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
//    crypto_sign_custom_asm_fortest(sm_res,&smlen_res,m,mlen,sk);
//    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
//    printf("crypto_sign_custom_fortest takes %d cycles\n",end-start);
//    //check sm and smlen
//    if(smlen_ref!=smlen_res){
//     printf("test_crypto_sign fail!\n");
//     return false;
//    }
//    bool flag_sm=true;
//    for(i=0;i<smlen_ref;i++){
//         if(sm_ref[i] != sm_res[i]) {
//             log_e("sm_ref[%d]=%d, while sm_res[%d]=%d", i, sm_ref[i], i, sm_res[i]);
//             flag_sm = false;
//         }
//    }
//    free(m);
//    free(sm_ref);
//    free(sm_res);
//    if(flag_sm) {
//         printf("test_crypto_sign test pass!\n\n");
//         return true;
//     }
//     else {
//         printf("test_crypto_sign test fail..\n\n");
//         return false;
//     }
// }

// bool test_crypto_sign_open(){
//     uint8_t pk[CRYPTO_PUBLICKEYBYTES];
//     uint8_t sk[CRYPTO_SECRETKEYBYTES];
//     uint8_t pk_ref[CRYPTO_PUBLICKEYBYTES];
//     uint8_t sk_ref[CRYPTO_SECRETKEYBYTES];
//     int i;
//     bool flag_pk, flag_sk;
//     uint64_t start, end;
//     /********************************
//      * CRYPTO_SIGN_KEYPAIR
//     ********************************/

//     crypto_sign_keypair(pk_ref,sk_ref);
//     crypto_sign_keypair_custom_asm_fortest(pk,sk);

//     //check pk
//     flag_pk = true;
//     for(i=0;i<CRYPTO_PUBLICKEYBYTES;i++){
//         if(pk[i] != pk_ref[i]) {
//             log_e("pk_ref[%d]=%d, while pk[%d]=%d", i, pk_ref[i], i, pk[i]);
//             flag_pk = false;
//         }
//     }
//     if(flag_pk) {
//         printf("pk of test_crypto_sign_keypair test pass!\n\n");
//     }
//     else {
//         printf("pk of test_crypto_sign_keypair test fail..\n\n");
//     }

//     //check sk
//     flag_sk = true;
//     for(i = 0; i < CRYPTO_SECRETKEYBYTES; i++) {
//         if(sk[i] != sk_ref[i]) {
//             log_e("sk_ref[%d]=%d, while sk[%d]=%d", i, sk_ref[i], i, sk[i]);
//             flag_sk = false;
//         }
//     }  
//     if(flag_sk) {
//         printf("sk of test_crypto_sign_keypair test pass!\n\n");
//     }
//     else {
//         printf("sk of test_crypto_sign_keypair test fail..\n\n");
//     }

//     /********************************
//      * CRYPTO_SIGN
//     ********************************/
//    size_t mlen=33;
//    uint8_t* m=(uint8_t*)malloc(mlen*sizeof(uint8_t));
//    size_t smlen_ref;
//    size_t smlen_res;
//    uint8_t* sm_ref=(uint8_t*)malloc((mlen+CRYPTO_BYTES)*sizeof(uint8_t));
//    uint8_t* sm_res=(uint8_t*)malloc((mlen+CRYPTO_BYTES)*sizeof(uint8_t));
//    srand((unsigned)time(NULL));
//    for(i=0;i<mlen;i++){
//     m[i]=i%256;
//     // m[i]=rand()%256;
//    }
//    crypto_sign(sm_ref,&smlen_ref,m,mlen,sk_ref);
//    crypto_sign_custom_asm_fortest(sm_res,&smlen_res,m,mlen,sk);
//    //check sm and smlen
//    if(smlen_ref!=smlen_res){
//     printf("test_crypto_sign fail!\n");
//     return false;
//    }
//    bool flag_sm=true;
//    for(i=0;i<smlen_ref;i++){
//         if(sm_ref[i] != sm_res[i]) {
//             log_e("sm_ref[%d]=%d, while sm_res[%d]=%d", i, sm_ref[i], i, sm_res[i]);
//             flag_sm = false;
//         }
//    }
//    if(flag_sm) {
//         printf("test_crypto_sign test pass!\n\n");
//     }
//     else {
//         printf("test_crypto_sign test fail..\n\n");
//     }

//     /********************************
//      * CRYPTO_SIGN_OPEN
//     ********************************/
//    size_t mlen1_ref;
//    size_t mlen1_res;
//    uint8_t* m1_ref=(uint8_t*)malloc((mlen+CRYPTO_BYTES)*sizeof(uint8_t));
//    uint8_t* m1_res=(uint8_t*)malloc((mlen+CRYPTO_BYTES)*sizeof(uint8_t));
//    int ret_val_ref=0;
//    int ret_val_res=0;
//    ret_val_ref=crypto_sign_open(m1_ref,&mlen1_ref,sm_ref,smlen_ref,pk_ref);
//    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
//    ret_val_res=crypto_sign_open_custom_asm_fortest(m1_res,&mlen1_res,sm_res,smlen_res,pk);
//    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
//    printf("crypto_sign_open_custom_fortest takes %d cycles\n",end-start);
//    if(ret_val_ref!=0){
//         printf("test_crypto_sign_open ref fail with bad ret\n");
//         return false;
//    }
//    if(ret_val_res!=0){
//         printf("test_crypto_sign_open res fail with bad ret\n");
//         return false;
//    }
//    if(mlen1_ref!=mlen){
//         printf("test_crypto_sign_open ref fail with bad mlen\n");
//         return false;
//    }
//    if(mlen1_res!=mlen){
//         printf("test_crypto_sign_open res fail with bad mlen\n");
//         return false;
//    }
//    //check m
//    bool flag_m_1=true;
//    for(i=0;i<mlen;i++){
//         if(m1_ref[i] != m1_res[i]) {
//             log_e("m1_ref[%d]=%d, while m1_res[%d]=%d", i, m1_ref[i], i, m1_res[i]);
//             flag_m_1 = false;
//         }
//    }
//     if(!flag_m_1) {
//         printf("m1_ref!=m1_res, test_crypto_sign_open fail..\n\n");
//         return false;
//     }

//     bool flag_m_2=true;
//     for(i=0;i<mlen;i++){
//         if(m1_res[i] != m[i]) {
//             log_e("m1_res[%d]=%d, while m[%d]=%d", i, m1_res[i], i, m[i]);
//             flag_m_2 = false;
//         }
//    }
//    if(flag_m_2) {
//         printf("m1_res=m, test_crypto_sign_open success!\n\n");
//         return true;
//     }
//     else {
//         printf("m1_res!=m, test_crypto_sign_open fail..\n\n");
//         return false;
//     }

//     //free 
//     free(m),free(sm_ref),free(sm_res),free(m1_ref),free(m1_res);
// }
#define CUS_OPT
bool test_crypto_sign_overall(){
    int i;
    uint64_t start, end;

    uint8_t pk[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk[CRYPTO_SECRETKEYBYTES];

    size_t mlen=59;
    uint8_t* m=(uint8_t*)malloc(mlen*sizeof(uint8_t));
    srand((unsigned)time(NULL));
    for(i=0;i<mlen;i++){
        // m[i]=i%256;
        m[i]=rand()%256;
    }
    size_t smlen;
    uint8_t* sm=(uint8_t*)malloc((mlen+CRYPTO_BYTES)*sizeof(uint8_t));

    size_t mlen1;
    uint8_t* m1=(uint8_t*)malloc((mlen+CRYPTO_BYTES)*sizeof(uint8_t));
    int ret_val=0;

    //cache warm up
    printf("start cache warmup\n");
    for(int i=0;i<3;i++){
        crypto_sign_keypair_custom(pk,sk);
        crypto_sign_custom(sm,&smlen,m,mlen,sk);
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
        ret_val=crypto_sign_open_custom(m1,&mlen1,sm,smlen,pk);
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
        printf("crypto_sign_keypair_custom takes %d cycles\n", end-start);
    }
    printf("end cache warmup\n");

    /********************************
     * CRYPTO_SIGN_KEYPAIR
    ********************************/
    #ifdef CUS_OPT
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
        // crypto_sign_keypair_custom_asm(pk,sk);
        crypto_sign_keypair_custom(pk,sk);
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
        // printf("crypto_sign_keypair_custom_asm takes %d cycles\n", end-start);
        printf("crypto_sign_keypair_custom takes %d cycles\n", end-start);
    #else
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
        crypto_sign_keypair(pk,sk);
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
        printf("crypto_sign_keypair takes %d cycles\n", end-start);
    #endif

   /********************************
     * CRYPTO_SIGN
    ********************************/
    #ifdef CUS_OPT
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
        // crypto_sign_custom_asm(sm,&smlen,m,mlen,sk);
        crypto_sign_custom(sm,&smlen,m,mlen,sk);
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
        // printf("crypto_sign_custom_asm takes %d cycles\n", end-start);
        printf("crypto_sign_custom takes %d cycles\n", end-start);
    #else
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
        crypto_sign(sm,&smlen,m,mlen,sk);
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
        printf("crypto_sign takes %d cycles\n", end-start);
    #endif


    /********************************
     * CRYPTO_SIGN_OPEN
    ********************************/
    #ifdef CUS_OPT
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
        // ret_val=crypto_sign_open_custom_asm(m1,&mlen1,sm,smlen,pk);
        ret_val=crypto_sign_open_custom(m1,&mlen1,sm,smlen,pk);
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
        // printf("crypto_sign_open_custom_asm takes %d cycles\n", end-start);
        printf("crypto_sign_open_custom takes %d cycles\n", end-start);
    #else
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
        ret_val=crypto_sign_open(m1,&mlen1,sm,smlen,pk);
        __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
        printf("crypto_sign_open takes %d cycles\n", end-start);
    #endif

    // if(ret_val!=0){
    //     printf("test_crypto_sign_open fail with bad ret\n");
    //     return false;
    // }
    // if(mlen1!=mlen){
    //     printf("test_crypto_sign_open fail with bad mlen\n");
    //     return false;
    // }
    bool flag_m=true;
    // for(i=0;i<mlen;i++){
    //     if(m1[i] != m[i]) {
    //         printf("m1[%d]=%d, while m[%d]=%d\n", i, m1[i], i, m[i]);
    //         flag_m = false;
    //     }
    // }
    //free 
    free(m),free(sm),free(m1);

    // if(flag_m) {
    //     printf("m1_res=m, test_crypto_sign_overall success!\n\n");
    //     return true;
    // }
    // else {
    //     printf("m1_res!=m, test_crypto_sign_overall fail..\n\n");
    //     return false;
    // }
}


int main(){
    printf("security level %d\n",DILITHIUM_MODE);
    // test_crypto_sign_keypair();
    // test_crypto_sign();
    // test_crypto_sign_open();//has already passed all security level tests
    test_crypto_sign_overall();//has already passed all security level tests
}