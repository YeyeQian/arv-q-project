#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include "param.h"
#include "fips202.h"
#include "symmetric.h"
#include "../apis/custom_inst_api.h"

void dilithium_shake128_stream_init_custom(const uint8_t seed[SEEDBYTES],uint16_t nonce)
{
  uint8_t t[SEEDBYTES+2];
  memcpy(t,seed,SEEDBYTES);
  t[SEEDBYTES]=nonce;
  t[SEEDBYTES+1]=nonce >> 8;
  csr_keccakmode_rw(SHAKE_128_MODE);
  keccak_absorb_custom(t,SEEDBYTES+2,SHAKE128_RATE);
}

void dilithium_shake256_stream_init_custom(const uint8_t seed[CRHBYTES],uint16_t nonce)
{
  uint8_t t[CRHBYTES+2];
  memcpy(t,seed,CRHBYTES);
  t[CRHBYTES]=nonce;
  t[CRHBYTES+1]=nonce>>8;
  csr_keccakmode_rw(SHAKE_256_MODE);
  keccak_absorb_custom(t,CRHBYTES+2,SHAKE256_RATE);
}