#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "encrypt.h"
#include "crypto_kem_mceliece348864.h"
#include "crypto_kem.h"

bool test_encrypt_custom(){
    uint8_t pk[crypto_kem_PUBLICKEYBYTES];
    srand((unsigned)time(NULL));
    for(int i=0;i<crypto_kem_PUBLICKEYBYTES;i++)pk[i]=rand()&255;

    unsigned char e_ref[ SYS_N/8 ];
    unsigned char e_res[ SYS_N/8 ];
    uint8_t ct_ref[crypto_kem_CIPHERTEXTBYTES];
    uint8_t ct_res[crypto_kem_CIPHERTEXTBYTES];

    uint64_t start, end;

    start=read_cycle();
    encrypt(ct_ref, pk, e_ref);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    encrypt_custom_asm(ct_res, pk, e_res);
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
    printf("gen_e return with flag %d\n",flag);
    for(int i=0;i<(crypto_kem_CIPHERTEXTBYTES);i++){
        if(ct_ref[i]!=ct_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_encrypt_custom return with flag %d\n",flag);

    return flag;
}

bool test_encrypt_custom_trans(){
    uint8_t pk[crypto_kem_PUBLICKEYBYTES];
    uint8_t pk_trans[crypto_kem_PUBLICKEYBYTES];
    srand((unsigned)time(NULL));
    for(int i=0;i<crypto_kem_PUBLICKEYBYTES;i++)pk[i]=rand()&255;

    //get pk_trans
    for(int i=0;i<PK_ROW_BYTES;i++){
        for(int j=0;j<PK_NROWS;j++){
            pk_trans[i*PK_NROWS+j]=pk[j*PK_ROW_BYTES+i];
        }
    }

    unsigned char e_ref[ SYS_N/8 ];
    unsigned char e_res[ SYS_N/8 ];
    uint8_t ct_ref[crypto_kem_CIPHERTEXTBYTES];
    uint8_t ct_res[crypto_kem_CIPHERTEXTBYTES];

    uint64_t start, end;

    start=read_cycle();
    encrypt(ct_ref, pk, e_ref);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    encrypt_custom_asm_trans(ct_res, pk_trans, e_res);
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
    printf("gen_e return with flag %d\n",flag);
    for(int i=0;i<(crypto_kem_CIPHERTEXTBYTES);i++){
        if(ct_ref[i]!=ct_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_encrypt_custom_trans return with flag %d\n",flag);

    return flag;
}

int main(){
    //test_encrypt_custom();
    test_encrypt_custom_trans();
    return 0;
}