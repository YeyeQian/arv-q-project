#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "parameters.h"
#include "api.h"

bool test_crypto_kem_keypair_custom(){
    uint8_t pk_seed[SEED_BYTES]={0};
    uint8_t sk_seed[SEED_BYTES]={0};

    unsigned char       pk_ref[CRYPTO_PUBLICKEYBYTES], sk_ref[CRYPTO_SECRETKEYBYTES];
    unsigned char       pk_res[CRYPTO_PUBLICKEYBYTES], sk_res[CRYPTO_SECRETKEYBYTES];

    srand((unsigned)time(NULL));
    for(int i=0;i<SEED_BYTES;i++){
        pk_seed[i]=rand()&255;
        sk_seed[i]=rand()&255;
    }

    uint64_t start, end;

    //cache warm up
    printf("start cache warmup\n");
    for(int i=0;i<3;i++){
        crypto_kem_keypair_custom(pk_res, sk_res, sk_seed, pk_seed);
    }
    printf("end cache warmup\n");

    start=read_cycle();
    crypto_kem_keypair(pk_ref, sk_ref, sk_seed, pk_seed);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    crypto_kem_keypair_custom(pk_res, sk_res, sk_seed, pk_seed);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<CRYPTO_PUBLICKEYBYTES;i++){
        if(pk_ref[i]!=pk_res[i]){
            flag=false;
            break;
        }
    }
    for(int i=0;i<CRYPTO_SECRETKEYBYTES;i++){
        if(sk_ref[i]!=sk_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_crypto_kem_keypair_custom return with flag %d\n",flag);

    return flag;
}

bool test_crypto_kem_enc_custom(){
    unsigned char pk[CRYPTO_PUBLICKEYBYTES];
    uint64_t m[VEC_K_SIZE_64];
    uint64_t salt[SALT_SIZE_64];
    srand((unsigned)time(NULL));
    for(int i=0;i<CRYPTO_PUBLICKEYBYTES;i++){
        pk[i]=rand()&255;
    }
    for(int i=0;i<VEC_K_SIZE_64;i++){
        m[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());
    }
    for(int i=0;i<SALT_SIZE_64;i++){
        salt[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());
    }

    unsigned char ct_ref[CRYPTO_CIPHERTEXTBYTES], ss_ref[CRYPTO_BYTES];
    unsigned char ct_res[CRYPTO_CIPHERTEXTBYTES], ss_res[CRYPTO_BYTES];

    uint64_t start, end;

    //cache warm up
    printf("start cache warmup\n");
    for(int i=0;i<3;i++){
        crypto_kem_enc_custom(ct_res, ss_res, pk,m,salt);
    }
    printf("end cache warmup\n");

    start=read_cycle();
    crypto_kem_enc(ct_ref, ss_ref, pk,m,salt);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    crypto_kem_enc_custom(ct_res, ss_res, pk,m,salt);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<CRYPTO_CIPHERTEXTBYTES;i++){
        if(ct_ref[i]!=ct_res[i]){
            flag=false;
            break;
        }
    }
    for(int i=0;i<CRYPTO_BYTES;i++){
        if(ss_ref[i]!=ss_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_crypto_kem_enc_custom return with flag %d\n",flag);

    return flag;
}

bool test_crypto_kem_dec_custom(){
    unsigned char ct[CRYPTO_CIPHERTEXTBYTES];
    unsigned char sk[CRYPTO_SECRETKEYBYTES];

    unsigned char ss1_ref[CRYPTO_BYTES];
    unsigned char ss1_res[CRYPTO_BYTES];

    srand((unsigned)time(NULL));
    for(int i=0;i<CRYPTO_SECRETKEYBYTES;i++){
        sk[i]=rand()&255;
    }
    for(int i=0;i<CRYPTO_CIPHERTEXTBYTES;i++){
        ct[i]=rand()&255;
    }

    uint64_t start, end;

    //cache warm up
    printf("start cache warmup\n");
    for(int i=0;i<3;i++){
        crypto_kem_dec_custom(ss1_res, ct, sk);
    }
    printf("end cache warmup\n");

    start=read_cycle();
    crypto_kem_dec(ss1_ref, ct, sk);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    crypto_kem_dec_custom(ss1_res, ct, sk);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<CRYPTO_BYTES;i++){
        if(ss1_ref[i]!=ss1_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_crypto_kem_dec_custom return with flag %d\n",flag);

    return flag;
}

bool test_crypto_kem_overall(){

    printf("PQCgenKAT_kem custom test for security level-%d:\n",192);
    unsigned char       ct[CRYPTO_CIPHERTEXTBYTES]={0}, ss[CRYPTO_BYTES]={0}, ss1[CRYPTO_BYTES]={0};
    unsigned char       pk[CRYPTO_PUBLICKEYBYTES]={0}, sk[CRYPTO_SECRETKEYBYTES]={0};
    int                 ret_val;
    uint64_t start,end;
    uint64_t m[VEC_K_SIZE_64];
    uint64_t salt[SALT_SIZE_64];
    uint8_t pk_seed[SEED_BYTES]={0};
    uint8_t sk_seed[SEED_BYTES]={0};
    srand((unsigned)time(NULL));
    for(int i=0;i<VEC_K_SIZE_64;i++){
        m[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());
    }
    for(int i=0;i<SALT_SIZE_64;i++){
        salt[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());
    }
    for(int i=0;i<SEED_BYTES;i++){
        pk_seed[i]=rand()&255;
        sk_seed[i]=rand()&255;
    }

    //cache warm up
    printf("start cache warmup\n");
    for(int i=0;i<2;i++){
        ret_val = crypto_kem_keypair_custom(pk, sk, sk_seed, pk_seed);
        ret_val = crypto_kem_enc_custom(ct, ss, pk,m,salt);
        ret_val = crypto_kem_dec_custom(ss1, ct, sk);
    }
    printf("end cache warmup\n");

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    ret_val = crypto_kem_keypair_custom(pk, sk, sk_seed, pk_seed);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_keypair_custom cost %lu cycles\n",end-start);
    if (ret_val != 0) {
        printf("crypto_kem_keypair returned <%d>\n", ret_val);
    }

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    ret_val = crypto_kem_enc_custom(ct, ss, pk,m,salt);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_enc_custom cost %lu cycles\n",end-start);
    if (ret_val != 0) {
        printf("crypto_kem_enc returned <%d>\n", ret_val);
    }

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    ret_val = crypto_kem_dec_custom(ss1, ct, sk);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("crypto_kem_dec_custom cost %lu cycles\n",end-start);

    return true;

}

int main(){
    //test_crypto_kem_keypair_custom();
    //test_crypto_kem_enc_custom();
    //test_crypto_kem_dec_custom();
    test_crypto_kem_overall();

    return 0;
}