#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "operations.h"
#include "crypto_kem_mceliece460896f.h"
#include "crypto_kem.h"

bool test_crypto_kem_keypair_custom(){
    uint8_t pk_ref[crypto_kem_PUBLICKEYBYTES];
    uint8_t sk_ref[crypto_kem_SECRETKEYBYTES];
    uint8_t pk_res[crypto_kem_PUBLICKEYBYTES];
    uint8_t sk_res[crypto_kem_SECRETKEYBYTES];

    int                 ret_val_ref;
    int                 ret_val_res;

    uint64_t start, end;

    start=read_cycle();
    ret_val_ref = crypto_kem_keypair(pk_ref, sk_ref);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    ret_val_res = crypto_kem_keypair_custom(pk_res, sk_res);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    if(ret_val_ref!=ret_val_res)flag=false;
    for(int i=0;i<(crypto_kem_PUBLICKEYBYTES);i++){
        if(pk_ref[i]!=pk_res[i]){
            flag=false;
            break;
        }
    }
    for(int i=0;i<(crypto_kem_SECRETKEYBYTES);i++){
        if(sk_ref[i]!=sk_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_crypto_kem_keypair_custom return with flag %d\n",flag);

    return flag;
}

bool test_crypto_kem_enc_custom(){
    uint8_t pk[crypto_kem_PUBLICKEYBYTES];
    srand((unsigned)time(NULL));
    for(int i=0;i<crypto_kem_PUBLICKEYBYTES;i++)pk[i]=rand()&255;
    uint8_t ct_ref[crypto_kem_CIPHERTEXTBYTES];
    uint8_t ct_res[crypto_kem_CIPHERTEXTBYTES];
    uint8_t ss_ref[crypto_kem_BYTES];
    uint8_t ss_res[crypto_kem_BYTES];
    int                 ret_val_ref;
    int                 ret_val_res;

    uint64_t start, end;

    start=read_cycle();
    ret_val_ref = crypto_kem_enc(ct_ref, ss_ref, pk);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    ret_val_res = crypto_kem_enc_custom(ct_res, ss_res, pk);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    if(ret_val_ref!=ret_val_res)flag=false;
    for(int i=0;i<crypto_kem_CIPHERTEXTBYTES;i++){
        if(ct_ref[i]!=ct_res[i]){
            flag=false;
            break;
        }
    }
    for(int i=0;i<crypto_kem_BYTES;i++){
        if(ss_ref[i]!=ss_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_crypto_kem_enc_custom return with flag %d\n",flag);

    return flag;
}

bool test_crypto_kem_dec_custom(){
    uint8_t ct[crypto_kem_CIPHERTEXTBYTES];
    uint8_t sk[crypto_kem_SECRETKEYBYTES];
    srand((unsigned)time(NULL));
    for(int i=0;i<crypto_kem_CIPHERTEXTBYTES;i++)ct[i]=rand()&255;
    for(int i=0;i<crypto_kem_SECRETKEYBYTES;i++)sk[i]=rand()&255;
    uint8_t ss1_ref[crypto_kem_BYTES];
    uint8_t ss1_res[crypto_kem_BYTES];

    int                 ret_val_ref;
    int                 ret_val_res;

    uint64_t start, end;

    start=read_cycle();
    ret_val_ref = crypto_kem_dec(ss1_ref, ct, sk);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    ret_val_res = crypto_kem_dec_custom(ss1_res, ct, sk);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    if(ret_val_ref!=ret_val_res)flag=false;
    for(int i=0;i<crypto_kem_BYTES;i++){
        if(ss1_ref[i]!=ss1_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_crypto_kem_dec_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    // test_crypto_kem_keypair_custom();
    // test_crypto_kem_enc_custom();
    test_crypto_kem_dec_custom();
    return 0;
}