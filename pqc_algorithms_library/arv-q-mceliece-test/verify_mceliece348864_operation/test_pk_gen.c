#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "pk_gen.h"
#include "crypto_kem_mceliece348864.h"
#include "crypto_kem.h"

bool test_pk_gen_custom(){
    uint32_t perm[ 1 << GFBITS ];
    uint8_t sk[IRR_BYTES];
    srand((unsigned)time(NULL));
    for(int i=0;i<(1<<GFBITS);i++){
        perm[i]=((uint32_t)rand()<<16)|(uint32_t)rand();
    }
    for(int i=0;i<IRR_BYTES;i++)sk[i]=rand()&255;

    int16_t pi_ref[ 1 << GFBITS ];
    uint8_t pk_ref[crypto_kem_PUBLICKEYBYTES];
    int16_t pi_res[ 1 << GFBITS ];
    uint8_t pk_res[crypto_kem_PUBLICKEYBYTES];

    uint64_t start, end;

    start=read_cycle();
    pk_gen(pk_ref, sk, perm, pi_ref);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    pk_gen_custom(pk_res, sk, perm, pi_res);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(1 << GFBITS);i++){
        if(pi_ref[i]!=pi_res[i]){
            flag=false;
            break;
        }
    }
    for(int i=0;i<(crypto_kem_PUBLICKEYBYTES);i++){
        if(pk_ref[i]!=pk_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_pk_gen_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    test_pk_gen_custom();

    return 0;
}