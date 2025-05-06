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
#include "poly_reorg.h"

bool test_indcpa_keypair_custom_asm()
{
  uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES];
  uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES];
  uint8_t pk_ref[KYBER_INDCPA_PUBLICKEYBYTES];
  uint8_t sk_ref[KYBER_INDCPA_SECRETKEYBYTES];  
  int i;
  bool flag_pk, flag_sk;
  uint64_t start, end;

  start = read_cycle();  
  indcpa_keypair_custom(pk_ref, sk_ref);
	end = read_cycle();
	printf("indcpa_keypair_custom takes %ld cycles\n", end-start);

  start = read_cycle();  
  indcpa_keypair_custom_asm(pk, sk);  
	end = read_cycle();
	printf("indcpa_keypair_custom_asm takes %ld cycles\n", end-start);

  // verify pk
  flag_pk = true;
  for(i = 0; i < KYBER_INDCPA_PUBLICKEYBYTES; i++) {
    if(pk[i] != pk_ref[i]) {
      printf("pk_ref[%d]=%d, while pk[%d]=%d\n", i, pk_ref[i], i, pk[i]);
      flag_pk = false;
    }
  }  
  if(flag_pk) {
    printf("pk of test_indcpa_keypair_custom_asm test pass!\n");
  }
  else {
    printf("pk of test_indcpa_keypair_custom_asm test fail..\n");
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
    printf("sk of test_indcpa_keypair_custom_asm test pass!\n\n");
  }
  else {
    printf("sk of test_indcpa_keypair_custom_asm test fail..\n\n");
  }

  return (flag_pk && flag_sk);   

}

bool test_indcpa_enc_custom_asm()
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
  uint64_t start, end;

  start = read_cycle();  
  indcpa_keypair_custom(pk_ref, sk_ref);
	end = read_cycle();
	printf("indcpa_keypair_custom takes %ld cycles\n", end-start);

  start = read_cycle();  
  indcpa_keypair_custom_asm(pk, sk);  
	end = read_cycle();
	printf("indcpa_keypair_custom_asm takes %ld cycles\n", end-start);  

  // verify pk
  flag_pk = true;
  for(i = 0; i < KYBER_INDCPA_PUBLICKEYBYTES; i++) {
    if(pk[i] != pk_ref[i]) {
      printf("pk_ref[%d]=%d, while pk[%d]=%d\n", i, pk_ref[i], i, pk[i]);
      flag_pk = false;
    }
  }  
  if(flag_pk) {
    printf("pk of test_indcpa_keypair_custom_asm test pass with flag_pk=%d\n", flag_pk);
  }
  else {
    printf("pk of test_indcpa_keypair_custom_asm test fail with flag_pk=%d\n", flag_pk);
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
    printf("sk of test_indcpa_keypair_custom_asm test pass with flag_sk=%d\n", flag_sk);
  }
  else {
    printf("sk of test_indcpa_keypair_custom_asm test fail with flag_sk=%d\n", flag_sk);
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
  
  start = read_cycle();  
  indcpa_enc_custom(c_ref, m, pk_ref, coins);
  end = read_cycle();
  printf("indcpa_enc_custom takes %ld cycles\n", end-start);

  start = read_cycle();  
  indcpa_enc_custom_asm(c, m, pk, coins);
  end = read_cycle();
  printf("indcpa_enc_custom_asm takes %ld cycles\n", end-start);

  // verify c
  flag_c = true;
  for(i = 0; i < KYBER_INDCPA_BYTES; i++) {
    if(c_ref[i] != c[i]) {
      printf("c_ref[%d]=%d, while c[%d]=%d\n", i, c_ref[i], i, c[i]);
      flag_c = false;
    }
  }

  if(flag_c) {
    printf("c of test_indcpa_enc_custom test pass with flag_c=%d\n", flag_c);
  }
  else {
    printf("c of test_indcpa_enc_custom test fail with flag_c=%d\n", flag_c);
  }

  return flag_c;
}

bool test_indcpa_dec_custom_asm()
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
  uint64_t start, end;

  start = read_cycle();  
  indcpa_keypair_custom(pk_ref, sk_ref);
	end = read_cycle();
	printf("indcpa_keypair_custom takes %ld cycles\n", end-start);

  start = read_cycle();  
  indcpa_keypair_custom_asm(pk, sk);  
	end = read_cycle();
	printf("indcpa_keypair_custom_asm takes %ld cycles\n", end-start);   

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
  
  start = read_cycle();  
  indcpa_enc_custom(c_ref, m, pk_ref, coins);
  end = read_cycle();
  printf("indcpa_enc_custom takes %ld cycles\n", end-start);

  start = read_cycle();  
  indcpa_enc_custom_asm(c, m, pk, coins);
  end = read_cycle();
  printf("indcpa_enc_custom_asm takes %ld cycles\n", end-start);

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

  start = read_cycle();
  indcpa_dec_custom(m_dec_ref, c_ref, sk_ref);
  end = read_cycle();
  printf("indcpa_dec_custom takes %ld cycles\n", end-start);

  start = read_cycle();
  indcpa_dec_custom_asm(m_dec, c, sk);
  end = read_cycle();
  printf("indcpa_dec_custom_asm takes %ld cycles\n", end-start);

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
    printf("c of test_indcpa_dec_custom_asm test pass with flag_m=%d\n", flag_m);
  }
  else {
    printf("c of test_indcpa_dec_custom_asm test fail with flag_m=%d\n", flag_m);
  } 

  return flag_m;
}

int main()
{
    test_indcpa_dec_custom_asm();
    return 0;
}