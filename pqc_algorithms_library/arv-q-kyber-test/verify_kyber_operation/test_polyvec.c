#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "api.h"
#include "compiler.h"
#include "fips202.h"
#include "symmetric.h"
#include "poly.h"
#include "polyvec.h"
#include "ntt.h"
#include "reduce.h"
#include "params.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>


bool test_polyvec_compress_custom()
{   
    polyvec a;
    uint8_t r_ref[KYBER_POLYVECCOMPRESSEDBYTES];
    uint8_t r[KYBER_POLYVECCOMPRESSEDBYTES];
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for(int i=0;i<KYBER_K;i++){
        for(int j=0;j<KYBER_N;j++){
            a.vec[i].coeffs[j] = rand() % KYBER_Q;
        }
    }

    polyvec_compress(r_ref,&a);

    polyvec_compress_custom(r,&a);

    flag = true;
    for(i = 0; i < KYBER_POLYVECCOMPRESSEDBYTES; i++) {
        if(r[i] != r_ref[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("polyvec_compress_custom test pass!\n");
    }
    else {
        printf("polyvec_compress_custom test fail..\n");
    }

    return flag;
}

bool test_polyvec_decompress_custom()
{
    uint8_t a[KYBER_POLYVECCOMPRESSEDBYTES];
    polyvec r_ref,r;
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for(i=0;i<KYBER_POLYVECCOMPRESSEDBYTES;i++){
        a[i]=rand()&255;
    }

    polyvec_decompress(&r_ref,a);

    polyvec_decompress_custom(&r,a);

    flag=true;
    for(int j=0;j<KYBER_K;j++){
        for(int i=0;i<KYBER_N;i++){
        if(r.vec[j].coeffs[i]!=r_ref.vec[j].coeffs[i]){
            flag = false;
        }
    }
    }
    if(flag) {
        printf("polyvec_decompress_custom test pass!\n");
    }
    else {
        printf("polyvec_decompress_custom test fail..\n");
    }

    return flag;
}

int main()
{
    /* 
    ** the following funcs passed tests in KYBER_K = 2, 3, 4
    ** and passed tests in VLEN = 256 and 512
    */    
    test_polyvec_compress_custom();
    // test_polyvec_decompress_custom();
    return 0;
}