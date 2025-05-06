#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"
#include "fips202.h"
#include "symmetric.h"
#include "poly.h"
#include "polyvec.h"
#include "ntt.h"
#include "reduce.h"
#include "params.h"
#include "poly_reorg.h"

bool test_poly_eta1_add_ntt()
{
    poly r, r_ref;
    uint8_t seed[KYBER_SYMBYTES];
    uint8_t nonce = 0;
    int i;
    bool flag;

    srand((unsigned)time(NULL));
    for(i = 0; i < KYBER_SYMBYTES; i++) {
        seed[i] = rand() % 0xff;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    poly_eta1_add_ntt(&r_ref, seed, nonce);
    poly_eta1_add_ntt_asm(&r, seed, nonce);

    flag = true;
    for(int i = 0;i < KYBER_N; i++){
        if(r.coeffs[i] != r_ref.coeffs[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("poly_poly_eta1_add_ntt_custom test pass with flag=%d\n", flag);
    }
    else {
        printf("poly_poly_eta1_add_ntt_custom test fail with flag=%d\n", flag);
    }

    return flag;
}

bool test_poly_eta2_add()
{
    poly r, r_ref;
    uint8_t seed[KYBER_SYMBYTES];
    uint8_t nonce = 0;
    int i;
    bool flag;

    srand((unsigned)time(NULL));
    for(i = 0; i < KYBER_SYMBYTES; i++) {
        seed[i] = rand() % 0xff;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    poly_eta2_add(&r_ref, seed, nonce);
    poly_eta2_add_asm(&r, seed, nonce);

    flag = true;
    for(int i = 0;i < KYBER_N; i++){
        if(r.coeffs[i] != r_ref.coeffs[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("poly_eta2_add test pass with flag=%d\n", flag);
    }
    else {
        printf("poly_eta2_add test fail with flag=%d\n", flag);
    }

    return flag;
}

bool test_polyvec_base_mul_acc_tomont()
{
    poly r, r_ref;
    polyvec a, b;
    int i, j;
    bool flag;

    srand((unsigned)time(NULL));
    for(i = 0; i < KYBER_K; i++) {
        for(j = 0; j < KYBER_N; j++) {
            a.vec[i].coeffs[j] = rand() % KYBER_Q;
            b.vec[i].coeffs[j] = rand() % KYBER_Q;
        }
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    polyvec_base_mul_acc_tomont_custom(&r_ref, &a, &b);

    polyvec_base_mul_acc_tomont_custom_asm(&r, &a, &b);

    flag = true;
    for(int i = 0;i < KYBER_N; i++){
        if(r.coeffs[i] != r_ref.coeffs[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
            flag = false;
        }
    }    

    if(flag) {
        printf("polyvec_base_mul_acc_tomont test pass!\n");
    }
    else {
        printf("polyvec_base_mul_acc_tomont test fail..\n");
    }

    return flag;    
}

bool test_polyvec_base_mul_acc_intt_tomont()
{
    poly r, r_ref;
    polyvec a, b;
    int i, j;
    bool flag;

    srand((unsigned)time(NULL));
    for(i = 0; i < KYBER_K; i++) {
        for(j = 0; j < KYBER_N; j++) {
            a.vec[i].coeffs[j] = rand() % KYBER_Q;
            b.vec[i].coeffs[j] = rand() % KYBER_Q;
        }
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    polyvec_base_mul_acc_intt_tomont_custom(&r_ref, &a, &b);

    polyvec_base_mul_acc_intt_tomont_custom_asm(&r, &a, &b);

    flag = true;
    for(int i = 0;i < KYBER_N; i++){
        if(r.coeffs[i] != r_ref.coeffs[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
            flag = false;
        }
    }    

    if(flag) {
        printf("test_polyvec_base_mul_acc_intt_tomont test pass!\n");
    }
    else {
        printf("test_polyvec_base_mul_acc_intt_tomont test fail..\n");
    }

    return flag;    
}

int main()
{
    test_poly_eta2_add();

    return 0;
}
