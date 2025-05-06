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


bool absolute_justify(int16_t a, int16_t b) {
    bool flag;

    int16_t diff = a - b;

    flag = ((diff % KYBER_Q) == 0);

    return flag;
}

bool test_poly_compress_custom()
{   
    poly a;
    uint8_t r_ref[KYBER_POLYCOMPRESSEDBYTES];
    uint8_t r[KYBER_POLYCOMPRESSEDBYTES];
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for (i = 0; i < KYBER_N; i++) {
        a.coeffs[i] = rand() & KYBER_Q;
    }    

    poly_compress(r_ref, &a);

    poly_compress_custom(r, &a);

    flag = true;
    for(i = 0; i < KYBER_POLYCOMPRESSEDBYTES; i++) {
        if(r[i] != r_ref[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("poly_compress_custom test pass!\n");
    }
    else {
        printf("poly_compress_custom test fail..\n");
    }

    return flag;
}

bool test_poly_decompress_custom()
{
    uint8_t a[KYBER_POLYCOMPRESSEDBYTES];
    poly r_ref,r;
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for(i=0;i<KYBER_POLYCOMPRESSEDBYTES;i++){
        a[i]=rand()&255;
    }

    poly_decompress(&r_ref,a);

    poly_decompress_custom(&r,a);

    flag=true;
    for(int i=0;i<KYBER_N;i++){
        if(r.coeffs[i]!=r_ref.coeffs[i]){
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("poly_decompress_custom test pass!\n");
    }
    else {
        printf("poly_decompress_custom test fail..\n");
    }

    return flag;
}

bool test_poly_tobytes_custom()
{   
    poly a;
    uint8_t r_ref[KYBER_POLYBYTES];
    uint8_t r[KYBER_POLYBYTES];
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for (i = 0; i < KYBER_N; i++) {
        a.coeffs[i] = rand() % KYBER_Q;
    }    

    poly_tobytes(r_ref,&a);

    poly_tobytes_custom(r,&a);

    flag = true;
    for(i = 0; i < KYBER_POLYBYTES; i++) {
        if(r[i] != r_ref[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("poly_tobytes_custom test pass!\n");
    }
    else {
        printf("poly_tobytes_custom test fail..\n");
    }

    return flag;
}

bool test_poly_frombytes_custom()
{
    uint8_t a[KYBER_POLYBYTES];
    poly r_ref,r;
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for(i=0;i<KYBER_POLYBYTES;i++){
        a[i]=rand()&255;
    }

    poly_frombytes(&r_ref,a);

    poly_frombytes_custom(&r,a);

    flag=true;
    for(int i=0;i<KYBER_N;i++){
        if(r.coeffs[i]!=r_ref.coeffs[i]){
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("poly_frombytes_custom test pass!\n");
    }
    else {
        printf("poly_frombytes_custom test fail..\n");
    }

    return flag;
}

bool test_poly_frommsg_custom()
{
    uint8_t a[KYBER_INDCPA_MSGBYTES];
    poly r_ref,r;
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for(i=0;i<KYBER_INDCPA_MSGBYTES;i++){
        a[i]=rand()&255;
    }

    poly_frommsg(&r_ref,a);

    poly_frommsg_custom(&r,a);

    flag=true;
    for(int i=0;i<KYBER_N;i++){
        if(r.coeffs[i]!=r_ref.coeffs[i]){
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("poly_frommsg_custom test pass!\n");
    }
    else {
        printf("poly_frommsg_custom test fail..\n");
    }

    return flag;
}

bool test_poly_tomsg_custom()
{   
    poly a;
    uint8_t msg_ref[KYBER_INDCPA_MSGBYTES];
    uint8_t msg[KYBER_INDCPA_MSGBYTES];
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for (i = 0; i < KYBER_N; i++) {
        a.coeffs[i] = rand() % (KYBER_Q<<1);
    }    
    
    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    poly_tomsg(msg_ref,&a);

    poly_tomsg_custom(msg,&a);

    flag = true;
    for(i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
        if(msg[i] != msg_ref[i]) {
            printf("msg_ref[%d]=%d, while msg[%d]=%d\n", i, msg_ref[i], i, msg[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("poly_tomsg_custom test pass!\n");
    }
    else {
        printf("poly_tomsg_custom test fail..\n");
    }

    return flag;
}

bool test_poly_getnoise_eta1_custom()
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

    poly_getnoise_eta1(&r_ref, seed, nonce);
    poly_getnoise_eta1_custom(&r, seed, nonce);

    flag = true;
    for(int i = 0;i < KYBER_N; i++){
        if(r.coeffs[i] != r_ref.coeffs[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("poly_getnoise_eta1_custom test pass!\n");
    }
    else {
        printf("poly_getnoise_eta1_custom test fail..\n");
    }

    return flag;
}

bool test_poly_getnoise_eta2_custom()
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

    poly_getnoise_eta2(&r_ref, seed, nonce);
    poly_getnoise_eta2_custom(&r, seed, nonce);

    flag = true;
    for(int i = 0;i < KYBER_N; i++){
        if(r.coeffs[i] != r_ref.coeffs[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
            flag = false;
        }
    }
    if(flag) {
        printf("poly_getnoise_eta2_custom test pass!\n");
    }
    else {
        printf("poly_getnoise_eta2_custom test fail..\n");
    }

    return flag;
}

bool test_poly_ntt_custom()
{
    poly r, r_ref;
    int i;

    for(i = 0; i < KYBER_N; i++) {
        r.coeffs[i] = i;
        r_ref.coeffs[i] = i;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    poly_ntt(&r_ref);
    bitreverse_to_standard_all(r_ref.coeffs);

    poly_ntt_custom(&r);

    bool flag = true;
    for(i = 0; i < 256; i++) {
        if(r.coeffs[i] != r_ref.coeffs[i]) {
            if(!absolute_justify(r.coeffs[i], r_ref.coeffs[i])) {
                printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
                flag = false;
            }

        }
    }    

    if(flag) {
        printf("poly_ntt_custom test pass!\n");
    }
    else {
        printf("poly_ntt_custom test fail..\n");
    }

    return flag;    
}

/*
** test_poly_invntt_tomont_custom also test poly_tomont_custom,
** because poly_tomont_custom is called in poly_invntt_tomont_custom
*/
bool test_poly_invntt_tomont_custom()   
{
    poly r, r_ref;
    int i;

    for(i = 0; i < KYBER_N; i++) {
        r.coeffs[i] = i;
        r_ref.coeffs[i] = i;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    poly_ntt(&r_ref);
    poly_invntt_tomont(&r_ref);

    poly_ntt_custom(&r);
    poly_invntt_tomont_custom(&r);

    bool flag = true;
    for(i = 0; i < 256; i++) {
        if(r.coeffs[i] != r_ref.coeffs[i]) {
            if(!absolute_justify(r.coeffs[i], r_ref.coeffs[i])) {
                printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
                flag = false;
            }

        }
    }    

    if(flag) {
        printf("poly_invntt_tomont_custom test pass!\n");
    }
    else {
        printf("poly_invntt_tomont_custom test fail..\n");
    }

    return flag;    
}

bool test_poly_basemul_montgomery_custom()
{
	unsigned int i;
	poly a0_ref, a1_ref, r_ref;
	poly a0, a1, r;
	bool flag;

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

	for (i = 0; i < 256; i++) {
		a0_ref.coeffs[i] = i + 1;
		a0.coeffs[i] = a0_ref.coeffs[i];
	}
	for (i = 0; i < 256; i++) {
		a1_ref.coeffs[i] = 0;
		a1.coeffs[i] = a1_ref.coeffs[i];
		r_ref.coeffs[i] = 0;
		r.coeffs[i] = r_ref.coeffs[i];
	}
	a1_ref.coeffs[12] = 10;
	a1.coeffs[12] = a1_ref.coeffs[12];

	poly_ntt(&a0_ref);
	poly_ntt(&a1_ref);
	poly_basemul_montgomery(&r_ref, &a0_ref, &a1_ref);
	bitreverse_to_standard_all(r_ref.coeffs);

	poly_ntt_custom(&a0);
	poly_ntt_custom(&a1);
	poly_basemul_montgomery_custom(&r, &a0, &a1);

    flag = true;
    for(i = 0; i < 256; i++) {
        if(r.coeffs[i] != r_ref.coeffs[i]) {
            if(!absolute_justify(r.coeffs[i], r_ref.coeffs[i])) {
                printf("r_ref.coeffs[%d]=%d, while r.coeffs[%d]=%d\n", i, r_ref.coeffs[i], i, r.coeffs[i]);
                flag = false;
            }

        }
    }    

    if(flag) {
        printf("poly_basemul_montgomery_custom test pass!\n");
    }
    else {
        printf("poly_basemul_montgomery_custom test fail..\n");
    }

    return flag;    
}



int main() 
{
    /* 
    ** the following funcs passed tests in KYBER_K = 2, 3, 4
    ** and passed tests in VLEN = 256 and 512
    */    
    test_poly_compress_custom();
    test_poly_decompress_custom();
    test_poly_tobytes_custom();
    test_poly_frombytes_custom();
    test_poly_frommsg_custom();
    test_poly_tomsg_custom();
    test_poly_ntt_custom();
    test_poly_invntt_tomont_custom();
    test_poly_basemul_montgomery_custom();
    test_poly_getnoise_eta1_custom();
    test_poly_getnoise_eta2_custom();

    return 0;
}