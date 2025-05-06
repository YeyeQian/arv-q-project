#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"
#include "hash_wrapper.h"
#define SHAKE256_RATE 136
/*************************************************
* Name:        shake256_prng_custom
*
* Description: SHAKE256 XOF with non-incremental API
*
* Arguments:   - uint8_t *out:      pointer to output
*              - size_t outlen:     requested output length in bytes
*              - const uint8_t *in: pointer to input
*              - size_t inlen:      length of input in bytes
**************************************************/
void shake256_prng_custom(uint8_t *out, int outlen, const uint8_t *in, size_t inlen)
{
  // unsigned int i;
  // size_t nblocks = outlen/SHAKE256_RATE;
  // uint8_t t[SHAKE256_RATE];
  // keccak_state state;

  // shake256_absorb(&state, in, inlen);
  // shake256_squeezeblocks(out, nblocks, &state);

  // out += nblocks*SHAKE256_RATE;
  // outlen -= nblocks*SHAKE256_RATE;

  // if(outlen) {
  //   shake256_squeezeblocks(t, 1, &state);
  //   for(i=0;i<outlen;i++)
  //     out[i] = t[i];
  // }
  uint32_t i;
  size_t nblocks = outlen/SHAKE256_RATE;
  uint8_t t[SHAKE256_RATE];

  csr_keccakmode_rw(SHAKE_256_MODE);
  keccak_absorb_custom(in, inlen, SHAKE256_RATE);
  keccak_squeezeblocks_custom(out, nblocks, SHAKE256_RATE, false);

  out += nblocks*SHAKE256_RATE;
  outlen -= nblocks*SHAKE256_RATE;
  if(outlen > 0) {
    if(nblocks == 0) {
      keccak_squeezeblocks_custom(t, 1, SHAKE256_RATE, false);
    }
    else {
      keccak_squeezeblocks_custom(t, 1, SHAKE256_RATE, true);
    }
    for(i = 0; i < outlen; i++)
      out[i] = t[i];
  }   
}
