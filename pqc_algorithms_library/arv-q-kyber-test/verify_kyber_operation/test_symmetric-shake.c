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


#define SHAKE_256_OUT_LEN KYBER_ETA2*KYBER_N/4

bool test_kyber_shake256_prf_custom()
{
    uint8_t key[KYBER_SYMBYTES];
    uint8_t hashvalue_ref[SHAKE_256_OUT_LEN];
    uint8_t hashvalue[SHAKE_256_OUT_LEN];
    int i;
    bool flag;
    uint8_t nonce;

    srand((unsigned)time(NULL));
    nonce = rand() & 0xff;

    for (i = 0; i < KYBER_SYMBYTES; i++) {
        key[i] = rand() & 0xff;
    }  

    kyber_shake256_prf(hashvalue_ref, SHAKE_256_OUT_LEN, key, nonce);
    kyber_shake256_prf_custom(hashvalue, SHAKE_256_OUT_LEN, key, nonce);
    
    flag = true;
    for(i = 0; i < SHAKE_256_OUT_LEN; i++) {
        if(hashvalue[i] != hashvalue_ref[i]) {
            printf("h_ref[%d]=%d, while h[%d]=%d\n", i, hashvalue_ref[i], i, hashvalue[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("kyber_shake256_prf_custom test pass!\n");
    }
    else {
        printf("kyber_shake256_prf_custom test fail..\n");
    }
}

int main()
{   
    /* 
    ** passed tests in SHAKE_256_OUT_LEN = KYBER_ETA1*KYBER_N/4
    ** and SHAKE_256_OUT_LEN = KYBER_ETA2*KYBER_N/4
    */    
    test_kyber_shake256_prf_custom();
    return 0;
}