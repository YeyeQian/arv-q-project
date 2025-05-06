#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "decrypt.h"
#include "crypto_kem_mceliece460896f.h"
#include "crypto_kem.h"
#include "benes.h"
#include "synd.h"

bool test_decrypt_custom(){
    uint8_t sk[crypto_kem_SECRETKEYBYTES];
    uint8_t ct[crypto_kem_CIPHERTEXTBYTES];
    srand((unsigned)time(NULL));
    for(int i=0;i<crypto_kem_SECRETKEYBYTES;i++)sk[i]=rand()&255;
    for(int i=0;i<crypto_kem_CIPHERTEXTBYTES;i++)ct[i]=rand()&255;

    unsigned char e_ref[ SYS_N/8 ];
    unsigned char e_res[ SYS_N/8 ];
    unsigned char ret_decrypt_ref = 0;
    unsigned char ret_decrypt_res = 0;

    uint64_t start, end;

    start=read_cycle();
    ret_decrypt_ref = decrypt(e_ref, sk+40, ct);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    ret_decrypt_res = decrypt_custom(e_res, sk+40, ct);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    if(ret_decrypt_ref!=ret_decrypt_res)flag=false;
    for(int i=0;i<(SYS_N/8);i++){
        if(e_ref[i]!=e_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_decrypt_custom return with flag %d\n",flag);

    return flag;
}

bool test_support_gen(){
    uint8_t* sk=NULL;
    sk=(uint8_t*)malloc(crypto_kem_SECRETKEYBYTES);
    gf L[ SYS_N ]={0};
    uint64_t start, end;

    start=read_cycle();
    support_gen(L, sk+40);
    end=read_cycle();
    printf("support_gen finished with %lu cycles\n",end-start);

    free(sk);

    return true;
}

bool test_synd(){
    unsigned char r[ SYS_N/8 ];

	gf g[ SYS_T+1 ];
	gf L[ SYS_N ];

	gf s[ SYS_T*2 ];

    uint64_t start, end;

    start=read_cycle();
    synd(s, g, L, r);
    end=read_cycle();
    printf("mceliece460896f synd finished with %lu cycles\n",end-start);

    return true;
}

int main(){
    //test_decrypt_custom();
    //test_support_gen();
    test_synd();

    return 0;
}