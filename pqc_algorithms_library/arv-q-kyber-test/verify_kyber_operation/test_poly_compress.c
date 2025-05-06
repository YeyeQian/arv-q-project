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


int test_poly_compress_rvv()
{   
    poly a;
    uint8_t r_ref[KYBER_POLYCOMPRESSEDBYTES];
    uint8_t r[KYBER_POLYCOMPRESSEDBYTES];
    bool flag;
    int i;
    unsigned long start, end;

    srand((unsigned)time(NULL));
    for (i = 0; i < KYBER_N; i++) {
        a.coeffs[i] = rand() & KYBER_Q;
    }    

    poly_compress(r_ref, &a);

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    poly_compress_rvv(r, &a);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("poly_compress_rvv: %d cycles\n", end-start);

    flag = true;
    for(i = 0; i < KYBER_POLYCOMPRESSEDBYTES; i++) {
        if(r[i] != r_ref[i]) {
            printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("poly_compress_rvv test pass!\n");
    }
    else {
        printf("poly_compress_rvv test fail..\n");
    }

    return flag;
}

int test_poly_compress_custom()
{   
    poly a;
    uint8_t r_ref[KYBER_POLYCOMPRESSEDBYTES];
    uint8_t r[KYBER_POLYCOMPRESSEDBYTES];
    bool flag;
    int i;
    unsigned long start, end;

    srand((unsigned)time(NULL));
    for (i = 0; i < KYBER_N; i++) {
        a.coeffs[i] = rand() & KYBER_Q;
    }    

    poly_compress(r_ref, &a);

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    poly_compress_custom(r, &a);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("poly_compress_custom: %d cycles\n", end-start);

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

bool test_pack_5valid()
{
    uint8_t t[32];
    int i;
    uint8_t packed_r[20];
    uint8_t packed_r_ref[20];
    size_t vl;

    uint8_t *t_ptr = t;
    uint8_t *r_ptr = packed_r_ref;

    for(i = 0; i < 31; i++) {
        t[i] = 0x07;
    }

  for(i=0;i<32/8;i++) {
    r_ptr[0] = (t_ptr[0] >> 0) | (t_ptr[1] << 5);
    r_ptr[1] = (t_ptr[1] >> 3) | (t_ptr[2] << 2) | (t_ptr[3] << 7);
    r_ptr[2] = (t_ptr[3] >> 1) | (t_ptr[4] << 4);
    r_ptr[3] = (t_ptr[4] >> 4) | (t_ptr[5] << 1) | (t_ptr[6] << 6);
    r_ptr[4] = (t_ptr[6] >> 2) | (t_ptr[7] << 3);
    r_ptr += 5;
    t_ptr += 8;
  }

  t_ptr = t;
  r_ptr = packed_r;  

  vuint8m1_t vx;
  vl = vsetvl_e8m1(32);
  vx = vle8_v_u8m1((uint8_t*)t_ptr, vl);
  vx = pack_vx_u8m1(vx, 5);
  vse8_v_u8m1(r_ptr, vx, 20);

  bool flag = true;  
  for(i = 0; i < 20; i++) {
    if(packed_r[i] != packed_r_ref[i]) {
        printf("packed_r_ref[%d]=%x, while packed_r[%d]=%x", i, packed_r_ref[i], i, packed_r[i]);
        flag = false;
    }
  }  

    return flag;
}

int test_poly_tobytes_custom()
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

int test_poly_tomsg_custom()
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

int test_polyvec_compress_custom()
{   
    polyvec a;
    uint8_t r_ref[KYBER_POLYVECCOMPRESSEDBYTES];
    uint8_t r[KYBER_POLYVECCOMPRESSEDBYTES];
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for(int i=0;i<KYBER_K;i++){
        for(int j=0;j<KYBER_N;j++){
            a.vec[i].coeffs[j]=rand()%KYBER_Q;
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

int main()
{
    // test_poly_compress_rvv();
    // test_poly_compress_custom();
    //test_pack_5valid();
    //test_poly_tobytes_custom();
    //test_poly_tomsg_custom();
    test_polyvec_compress_custom();
    return 0;
}