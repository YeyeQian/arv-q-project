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
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>

#include "params.h"
#include "poly.h"
#include "indcpa.h"

bool test_rejsample()
{
    poly a_ref;
    poly a;
    int i;
    uint8_t buf[SHAKE128_RATE];
    unsigned int valid_num, valid_num_ref;
    unsigned long start, end;

    for(i = 0; i < KYBER_N; i++) {
        a_ref.coeffs[i] = 0;
        a.coeffs[i] = 0;
    }

    for(i = 0; i < SHAKE128_RATE; i++) {
        buf[i] = i;
    }

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    valid_num_ref = rej_uniform(a_ref.coeffs, KYBER_N, buf, SHAKE128_RATE);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("rej_uniform: %d cycles\n", end-start);

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    valid_num = rej_uniform_custom(a.coeffs, KYBER_N, buf, SHAKE128_RATE);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("rej_uniform_custom: %d cycles\n", end-start);

    if(valid_num_ref != valid_num) {
      printf("valid_num_ref=%d, while valid_num=%d\n", valid_num_ref, valid_num);
      return false;
    }

    bool flag = true;
    for(i = 0; i < valid_num_ref; i++) {
        if(a.coeffs[i] != a_ref.coeffs[i]) {
            printf("a_ref.coeffs[%d]=%d, while a.coeffs[%d]=%d\n", i, a_ref.coeffs[i], i, a.coeffs[i]);
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

bool test_unpack()
{
    int i;
    uint8_t buf[24];
    uint16_t unpacked_ref[16] = {0};
    uint16_t unpacked[16] = {0};

    for(i = 0; i < 24; i++) {
        buf[i] = i;
    }

  unsigned int ctr, pos, buflen;
  buflen = 24;
  uint16_t val0, val1;

  ctr = pos = 0;
  while(pos + 3 <= buflen) {
    val0 = ((buf[pos+0] >> 0) | ((uint16_t)buf[pos+1] << 8)) & 0xFFF;
    val1 = ((buf[pos+1] >> 4) | ((uint16_t)buf[pos+2] << 4)) & 0xFFF;
    pos += 3;
    unpacked_ref[ctr++] = val0;
    unpacked_ref[ctr++] = val1;
  }    
  
  size_t vl;
  vuint8m1_t vreg_m;
  vint16m1_t vreg_rej;
  vl = vsetvl_e8m1(24);
  vreg_m = vle8_v_u8m1(buf, vl);
  vl = vsetvl_e16m1((vl / 3) << 1);   // in such case, vl is always divisible by 3
  vreg_rej = unpack_vx_i16m1(vreg_m, 12);
  vse16_v_i16m1((int16_t*)unpacked, vreg_rej, vl);

  bool flag = true;
  for(i = 0; i < 16; i++) {
    if(unpacked_ref[i] != unpacked[i]) {
        printf("unpacked_ref[%d]=%x, while unpacked[%d]=%x", i, unpacked_ref[i], i, unpacked[i]);
        flag = false;
    }
  }

    if(flag) {
        printf("test_unpack test pass!\n");
    }
    else {
        printf("test_unpack test fail..\n");
    }

  return flag;
}

int main()
{
    test_rejsample();
    return 0;
}