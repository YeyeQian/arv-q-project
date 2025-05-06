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
#include "param.h"
#include "rounding.h"
#include "reduce.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "../common/debug.h"
#include <math.h>

bool test_rejsample(){
    poly a_ref,a;
    int i;
    uint8_t seed[SEEDBYTES]={
        41,35,190,132,225,108,214,174,
        82,144,73,241,241,187,233,235,
        179,166,219,60,135,12,62,153,
        36,94,13,28,6,183,71,222
    };
    uint8_t seed1[CRHBYTES]={
        41,35,190,132,225,108,214,174,
        82,144,73,241,241,187,233,235,
        179,166,219,60,135,12,62,153,
        36,94,13,28,6,183,71,222,
        179,18,77,200,67,187,139,166,
        31,3,90,125,9,56,37,31
    };
    uint16_t nonce=1145;
    poly_uniform_gamma1(&a_ref,seed1,nonce);
    poly_uniform_gamma1_custom(&a,seed1,nonce);
    bool flag = true;
    for(i = 0; i < N; i++) {
        if(a.coeffs[i] != a_ref.coeffs[i]) {
            log_e("a_ref.coeffs[%d]=%d, while a.coeffs[%d]=%d", i, a_ref.coeffs[i], i, a.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("test_rejsample test pass!\n");
    }
    else {
        printf("test_rejsample test fail..\n");
    }

    return flag;
}

bool test_poly_challenge(){
    uint64_t start, end;
    poly a_ref,a;
    int i;
    uint8_t seed[SEEDBYTES]={
        41,35,190,132,225,108,214,174,
        82,144,73,241,241,187,233,235,
        179,166,219,60,135,12,62,153,
        36,94,13,28,6,183,71,222
    };
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    poly_challenge(&a_ref,seed);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("poly_challenge takes %d cycles\n", end-start);

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    poly_challenge_custom(&a,seed);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("poly_challenge_custom takes %d cycles\n", end-start);

    bool flag = true;
    for(i = 0; i < N; i++) {
        if(a.coeffs[i] != a_ref.coeffs[i]) {
            log_e("a_ref.coeffs[%d]=%d, while a.coeffs[%d]=%d", i, a_ref.coeffs[i], i, a.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("test_poly_challenge test pass!\n");
    }
    else {
        printf("test_poly_challenge test fail..\n");
    }

    return flag;
}

bool test_poly_uniform_eta_custom(){
    poly ref,res;
    uint8_t seed[SEEDBYTES];
    uint16_t nonce;

    for(int i=0;i<SEEDBYTES;i++){
        seed[i]=i;
    }
    nonce=11451;

    poly_uniform_eta(&ref,seed,nonce);
    poly_uniform_eta_custom(&res,seed,nonce);

    bool flag = true;
    for(int i = 0; i < N; i++) {
        if(ref.coeffs[i] != res.coeffs[i]) {
            log_e("ref.coeffs[%d]=%d, while res.coeffs[%d]=%d", i, ref.coeffs[i], i, res.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("test_poly_uniform_eta_custom test pass!\n");
    }
    else {
        printf("test_poly_uniform_eta_custom test fail..\n");
    }

    return flag;
}

int main()
{
    // test_rejsample();
    // test_poly_challenge();
    test_poly_uniform_eta_custom();
    return 0;
}