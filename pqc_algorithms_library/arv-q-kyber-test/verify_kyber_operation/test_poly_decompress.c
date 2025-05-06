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
#include "params.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>


int test_poly_decompress_custom()
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

int test_poly_frombytes_custom()
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

int test_poly_frommsg_custom()
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

int test_polyvec_decompress_custom()
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
    test_polyvec_decompress_custom();
    return 0;
}