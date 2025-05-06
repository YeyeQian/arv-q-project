#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "encrypt.h"
#include "crypto_kem_mceliece6960119f.h"
#include "crypto_kem.h"
#include "encrypt.h"
#include "pk_sk_cons.h"

bool test_syndrome(){
    unsigned char e[ SYS_N/8 ];
    uint8_t ct_ref[crypto_kem_CIPHERTEXTBYTES];
    uint8_t ct_res[crypto_kem_CIPHERTEXTBYTES];
    srand(time(NULL));
    for(int i=0;i<(SYS_N/8);i++)e[i]=i;

    uint64_t start, end;
    start=read_cycle();
    syndrome(ct_ref, cons_pk, e);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    syndrome_custom(ct_res, cons_pk, e);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<crypto_kem_CIPHERTEXTBYTES;i++){
        if(ct_ref[i]!=ct_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_syndrome return with flag %d\n",flag);

    return flag;
}

bool test_gen_e(){
    unsigned char e_ref[ SYS_N/8 ];
    unsigned char e_res[ SYS_N/8 ];

    uint64_t start, end;
    start=read_cycle();
    gen_e(e_ref);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    gen_e_custom(e_res);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(SYS_N/8);i++){
        if(e_ref[i]!=e_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_gen_e return with flag %d\n",flag);

    return flag;
}

bool test_encrypt(){
    uint8_t ct_ref[crypto_kem_CIPHERTEXTBYTES];
    uint8_t ct_res[crypto_kem_CIPHERTEXTBYTES];
    unsigned char e_ref[ SYS_N/8 ];
    unsigned char e_res[ SYS_N/8 ];

    uint64_t start, end;
    start=read_cycle();
    encrypt(ct_ref,cons_pk,e_ref);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    encrypt_custom(ct_res,cons_pk,e_res);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(SYS_N/8);i++){
        if(e_ref[i]!=e_res[i]){
            flag=false;
            printf("e check fail\n");
            break;
        }
    }

    for(int i=0;i<(crypto_kem_CIPHERTEXTBYTES);i++){
        if(ct_ref[i]!=ct_res[i]){
            flag=false;
            printf("ct check fail\n");
            break;
        }
    }

    printf("test_encrypt return with flag %d\n",flag);

    return flag;
}

int main(){
    //test_syndrome();
    //test_gen_e();
    test_encrypt();
    return 0;
}