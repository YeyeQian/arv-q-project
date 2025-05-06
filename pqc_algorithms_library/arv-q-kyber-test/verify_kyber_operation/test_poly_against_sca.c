#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "compiler.h"
#include "fips202.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "params.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>

void test_poly_basemul_montgomery_custom_redundant(){
    size_t start, end;

    srand(time(NULL));
    uint8_t seed[32];
    for(int i=0;i<32;i++) {
        seed[i] = i + 60;
    }

    //Generate 128 bytes random array
    uint8_t buffer[128];
    shake128_custom(buffer,128,seed,32);

    poly a,b,ref,res;
    for(int i=0;i<KYBER_N;i++){
        a.coeffs[i]=rand()%KYBER_Q;
        b.coeffs[i]=rand()%KYBER_Q;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    poly_basemul_montgomery_custom(&ref,&a,&b);

    start = read_cycle();
    poly_basemul_montgomery_custom_redundant(&res,&a,&b,buffer);
    // poly_basemul_montgomery_custom_redundant_max(&res,&a,&b,buffer);
    end = read_cycle();

    printf("it takes %ld cycles\n", end - start);

    //check
    bool flag=true;
    for(int i=0;i<KYBER_N;i++){
        if(ref.coeffs[i]!=res.coeffs[i]){
            printf("ref.coeffs[%d]=%d, while res.coeffs[%d]=%d\n", i, ref.coeffs[i], i, res.coeffs[i]);
            flag=false;
        }
    }
    if(flag) {
        printf("test_poly_basemul_montgomery_custom_redundant test pass!\n");
    }
    else {
        printf("test_poly_basemul_montgomery_custom_redundant test fail..\n");
    }
    return;
}

void test_poly_basemul_montgomery_custom_shuffling() {
    size_t start, end;

    srand(time(NULL));
    uint8_t seed[32];
    for(int i=0;i<32;i++) {
        seed[i] = i + 60;
    }

    //Generate 128 bytes random array
    uint8_t buffer[128];
    shake128_custom(buffer,128,seed,32);

    poly a,b,ref,res;
    for(int i=0;i<KYBER_N;i++){
        a.coeffs[i]=rand()%KYBER_Q;
        b.coeffs[i]=rand()%KYBER_Q;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    poly_basemul_montgomery_custom(&ref,&a,&b);

    // start = read_cycle();
    poly_basemul_montgomery_custom_shuffling(&res,&a,&b,buffer);
    // end = read_cycle();

    // printf("it takes %ld cycles\n", end - start);

    //check
    bool flag=true;
    for(int i=0;i<KYBER_N;i++){
        if(ref.coeffs[i]!=res.coeffs[i]){
            printf("ref.coeffs[%d]=%d, while res.coeffs[%d]=%d\n", i, ref.coeffs[i], i, res.coeffs[i]);
            flag=false;
        }
    }
    if(flag) {
        printf("test_poly_basemul_montgomery_custom_shuffling test pass!\n");
    }
    else {
        printf("test_poly_basemul_montgomery_custom_shuffling test fail..\n");
    }
    return;
}

void test_poly_basemul_montgomery_custom_outside_shuffling() {
    size_t start, end;

    srand(time(NULL));
    uint8_t seed[32];
    for(int i=0;i<32;i++) {
        seed[i] = i + 60;
    }

    //Generate 128 bytes random array
    uint8_t buffer[128];
    shake128_custom(buffer,128,seed,32);

    uint16_t poly_idx[256];
    uint16_t zeta_idx[128];

    for(int i=0;i<128;i++) {
        zeta_idx[i]=i;
    }

    for(int i=127; i>=0; i--) {             //shuffling in pairs
        int selection=buffer[i]%(i+1);

        //shuffling zeta index
        uint16_t temp=zeta_idx[i];
        zeta_idx[i]=zeta_idx[selection];
        zeta_idx[selection]=temp;
    }

    // get poly index
    for(int i=0;i<128;i++){
        poly_idx[2*i]=zeta_idx[i]*2;
        poly_idx[2*i+1]=zeta_idx[i]*2+1;
    }

    poly a,b,ref,res;
    for(int i=0;i<KYBER_N;i++){
        a.coeffs[i]=rand()%KYBER_Q;
        b.coeffs[i]=rand()%KYBER_Q;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    poly_basemul_montgomery_custom(&ref,&a,&b);

    start = read_cycle();
    poly_basemul_montgomery_custom_outside_shuffling(&res, &a, &b, zeta_idx, poly_idx);
    end = read_cycle();

    printf("it takes %ld cycles\n", end - start);

    //check
    bool flag=true;
    for(int i=0;i<KYBER_N;i++){
        if(ref.coeffs[i]!=res.coeffs[i]){
            printf("ref.coeffs[%d]=%d, while res.coeffs[%d]=%d\n", i, ref.coeffs[i], i, res.coeffs[i]);
            flag=false;
        }
    }
    if(flag) {
        printf("test_poly_basemul_montgomery_custom_outside_shuffling test pass!\n");
    }
    else {
        printf("test_poly_basemul_montgomery_custom_outside_shuffling test fail..\n");
    }
    return;
}


int main() 
{  
    // test_poly_basemul_montgomery_custom_redundant();
    // test_poly_basemul_montgomery_custom_outside_shuffling();
    
    return 0;
}