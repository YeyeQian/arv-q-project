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
#include "polyvec.h"
#include "indcpa.h"

bool test_gen_matrix()
{
  polyvec a_ref[KYBER_K];
  polyvec a[KYBER_K];

  uint8_t seed[KYBER_SYMBYTES];
  int i, j, k;

  for(i = 0; i < KYBER_SYMBYTES; i++) {
    seed[i] = i;
  }

  gen_matrix(a_ref, seed, 0);

  gen_matrix_custom(a, seed, 0);

  bool flag = true;

  for(i = 0; i < KYBER_K; i++) {
    for(j = 0; j < KYBER_K; j++) {
      for(k = 0; k < KYBER_N; k++) {
        if(a_ref[i].vec[j].coeffs[k] != a[i].vec[j].coeffs[k]) {
          printf("a_ref[%d].vec[%d].coeffs[%d]=%d, while a[%d].vec[%d].coeffs[%d]=%d\n", i, j, k, a_ref[i].vec[j].coeffs[k], i, j, k, a[i].vec[j].coeffs[k]);
          flag = false;
        }
      }

    }
  }

  if(flag) {
    printf("test_gen_matrix test pass with flag=%d\n", flag);
  }
  else {
    printf("test_gen_matrix test fail with flag=%d\n", flag);
  }

  return flag;   

}

bool test_indcpa_keypair_custom()
{
  uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES];
  uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES];
  uint8_t pk_ref[KYBER_INDCPA_PUBLICKEYBYTES];
  uint8_t sk_ref[KYBER_INDCPA_SECRETKEYBYTES];  
  int i;
  bool flag_pk, flag_sk;

  indcpa_keypair(pk_ref, sk_ref);
  indcpa_keypair_custom_for_test(pk, sk);  

  // verify pk
  flag_pk = true;
  for(i = 0; i < KYBER_INDCPA_PUBLICKEYBYTES; i++) {
    if(pk[i] != pk_ref[i]) {
      printf("pk_ref[%d]=%d, while pk[%d]=%d\n", i, pk_ref[i], i, pk[i]);
      flag_pk = false;
    }
  }  
  if(flag_pk) {
    printf("pk of test_indcpa_keypair_custom test pass!\n");
  }
  else {
    printf("pk of test_indcpa_keypair_custom test fail..\n");
  }

  // verify sk
  flag_sk = true;
  for(i = 0; i < KYBER_INDCPA_SECRETKEYBYTES; i++) {
    if(sk[i] != sk_ref[i]) {
      printf("sk_ref[%d]=%d, while sk[%d]=%d\n", i, sk_ref[i], i, sk[i]);
      flag_sk = false;
    }
  }  
  if(flag_sk) {
    printf("sk of test_indcpa_keypair_custom test pass!\n\n");
  }
  else {
    printf("sk of test_indcpa_keypair_custom test fail..\n\n");
  }

  return (flag_pk && flag_sk);   

}

bool test_indcpa_enc_custom()
{
  int i;

  /*
  * indcpa_keypair stage
  */
  uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES];
  uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES];
  uint8_t pk_ref[KYBER_INDCPA_PUBLICKEYBYTES];
  uint8_t sk_ref[KYBER_INDCPA_SECRETKEYBYTES];  
  bool flag_pk, flag_sk;

  indcpa_keypair(pk_ref, sk_ref);
  indcpa_keypair_custom_for_test(pk, sk);  

  // verify pk
  flag_pk = true;
  for(i = 0; i < KYBER_INDCPA_PUBLICKEYBYTES; i++) {
    if(pk[i] != pk_ref[i]) {
      printf("pk_ref[%d]=%d, while pk[%d]=%d\n", i, pk_ref[i], i, pk[i]);
      flag_pk = false;
    }
  }  
  if(flag_pk) {
    printf("pk of test_indcpa_keypair_custom test pass!\n");
  }
  else {
    printf("pk of test_indcpa_keypair_custom test fail..\n");
  }

  // verify sk
  flag_sk = true;
  for(i = 0; i < KYBER_INDCPA_SECRETKEYBYTES; i++) {
    if(sk[i] != sk_ref[i]) {
      printf("sk_ref[%d]=%d, while sk[%d]=%d\n", i, sk_ref[i], i, sk[i]);
      flag_sk = false;
    }
  }  
  if(flag_sk) {
    printf("sk of test_indcpa_keypair_custom test pass!\n\n");
  }
  else {
    printf("sk of test_indcpa_keypair_custom test fail..\n\n");
  }

  /*
  * indcpa_enc stage
  */
  uint8_t coins[KYBER_SYMBYTES];
  uint8_t m[KYBER_INDCPA_MSGBYTES];
  uint8_t c_ref[KYBER_INDCPA_BYTES];
  uint8_t c[KYBER_INDCPA_BYTES];
  bool flag_c;

  for(i = 0; i < KYBER_SYMBYTES; i++) {
    coins[i] = i;
  }
  for(i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
    m[i] = i;
  }
  
  indcpa_enc(c_ref, m, pk_ref, coins);
  indcpa_enc_custom_for_test(c, m, pk, coins);

  // verify c
  flag_c = true;
  for(i = 0; i < KYBER_INDCPA_BYTES; i++) {
    if(c_ref[i] != c[i]) {
      printf("c_ref[%d]=%d, while c[%d]=%d\n", i, c_ref[i], i, c[i]);
      flag_c = false;
    }
  }

  if(flag_c) {
    printf("c of test_indcpa_enc_custom test pass!\n\n");
  }
  else {
    printf("c of test_indcpa_enc_custom test fail..\n\n");
  }  

  return flag_c;
}

bool test_indcpa_dec_custom()
{
  int i;

  /*
  * indcpa_keypair stage
  */
  uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES];
  uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES];
  uint8_t pk_ref[KYBER_INDCPA_PUBLICKEYBYTES];
  uint8_t sk_ref[KYBER_INDCPA_SECRETKEYBYTES];  
  bool flag_pk, flag_sk;

  indcpa_keypair(pk_ref, sk_ref);
  indcpa_keypair_custom_for_test(pk, sk);  

  // verify pk
  flag_pk = true;
  for(i = 0; i < KYBER_INDCPA_PUBLICKEYBYTES; i++) {
    if(pk[i] != pk_ref[i]) {
      printf("pk_ref[%d]=%d, while pk[%d]=%d\n", i, pk_ref[i], i, pk[i]);
      flag_pk = false;
    }
  }  
  if(flag_pk) {
    printf("pk of test_indcpa_keypair_custom test pass!\n");
  }
  else {
    printf("pk of test_indcpa_keypair_custom test fail..\n");
  }

  // verify sk
  flag_sk = true;
  for(i = 0; i < KYBER_INDCPA_SECRETKEYBYTES; i++) {
    if(sk[i] != sk_ref[i]) {
      printf("sk_ref[%d]=%d, while sk[%d]=%d\n", i, sk_ref[i], i, sk[i]);
      flag_sk = false;
    }
  }  
  if(flag_sk) {
    printf("sk of test_indcpa_keypair_custom test pass!\n");
  }
  else {
    printf("sk of test_indcpa_keypair_custom test fail..\n");
  }

  /*
  * indcpa_enc stage
  */
  uint8_t coins[KYBER_SYMBYTES];
  uint8_t m[KYBER_INDCPA_MSGBYTES];
  uint8_t c_ref[KYBER_INDCPA_BYTES];
  uint8_t c[KYBER_INDCPA_BYTES];
  bool flag_c;

  for(i = 0; i < KYBER_SYMBYTES; i++) {
    coins[i] = i;
  }
  for(i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
    m[i] = i;
  }
  
  indcpa_enc(c_ref, m, pk_ref, coins);
  indcpa_enc_custom_for_test(c, m, pk, coins);

  // verify c
  flag_c = true;
  for(i = 0; i < KYBER_INDCPA_BYTES; i++) {
    if(c_ref[i] != c[i]) {
      printf("c_ref[%d]=%d, while c[%d]=%d\n", i, c_ref[i], i, c[i]);
      flag_c = false;
    }
  }

  if(flag_c) {
    printf("c of test_indcpa_enc_custom test pass!\n");
  }
  else {
    printf("c of test_indcpa_enc_custom test fail..\n");
  }  

  /*
  * indcpa_dec stage
  */
  uint8_t m_dec[KYBER_INDCPA_MSGBYTES];
  uint8_t m_dec_ref[KYBER_INDCPA_MSGBYTES];
  bool flag_m;

  indcpa_dec(m_dec_ref, c_ref, sk_ref);
  indcpa_dec_custom_for_test(m_dec, c, sk);

  flag_m = true;
  for(i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
    if(m_dec_ref[i] != m_dec[i]) {
      printf("m_dec_ref[%d]=%d, while m_dec[%d]=%d\n", i, m_dec_ref[i], i, m_dec[i]);
      flag_m = false;
    }
    if(m_dec_ref[i] != m[i]) {
      printf("m_dec_ref[%d]=%d, while m[%d]=%d\n", i, m_dec_ref[i], i, m[i]);
      flag_m = false;
    }    
  }  

  if(flag_m) {
    printf("c of test_indcpa_dec_custom test pass!\n");
  }
  else {
    printf("c of test_indcpa_dec_custom test fail..\n");
  } 

  return flag_m;
}


bool test_indcpa_custom()
{
  int i;
  uint64_t start, end;

  /*
  * indcpa_keypair stage
  */
  uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES];
  uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES];

  start = read_cycle();
  indcpa_keypair_custom(pk, sk);
  end = read_cycle();

  printf("indcpa_keypair_custom takes %d cycles\n", end-start);

  /*
  * indcpa_enc stage
  */
  uint8_t coins[KYBER_SYMBYTES];
  uint8_t m[KYBER_INDCPA_MSGBYTES];
  uint8_t c[KYBER_INDCPA_BYTES];

  for(i = 0; i < KYBER_SYMBYTES; i++) {
    coins[i] = i;
  }
  for(i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
    m[i] = i;
  }
  start = read_cycle();
  indcpa_enc_custom(c, m, pk, coins);
  end = read_cycle();

  printf("indcpa_enc_custom takes %d cycles\n", end-start);
  /*
  * indcpa_dec stage
  */
  uint8_t m_dec[KYBER_INDCPA_MSGBYTES];
  bool flag_m;
  start = read_cycle();
  indcpa_dec_custom(m_dec, c, sk);
  end = read_cycle();

  printf("indcpa_dec_custom takes %d cycles\n", end-start);

  flag_m = true;
  for(i = 0; i < KYBER_INDCPA_MSGBYTES; i++) {
    if(m_dec[i] != m[i]) {
      printf("m_dec[%d]=%d, while m[%d]=%d\n", i, m_dec[i], i, m[i]);
      flag_m = false;
    }
  }  

  if(flag_m) {
    printf("test_indcpa_custom test pass!\n");
  }
  else {
    printf("test_indcpa_custom test fail..\n");
  } 

  return flag_m;
}

int main()
{
    /* 
    ** the following funcs passed tests in KYBER_K = 2, 3, 4
    ** and passed tests in VLEN = 256 and 512
    */
    // test_gen_matrix();
    test_indcpa_keypair_custom();
    //test_indcpa_custom();
    return 0;
}