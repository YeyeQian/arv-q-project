#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include "crypto_hash.h"

/*************************************************
* Name:        keccak_absorb_custom
*
* Description: Absorb step of Keccak;
*              non-incremental, starts by zeroeing the state.
*
**************************************************/
void keccak_absorb_custom(const uint8_t *m, int mlen, uint32_t rate)
{
  // csr_keccakmode_rw(keccakmode);
  // uint32_t rate;
  // case(keccakmode) {
  //   SHAKE_128_MODE: {
  //     rate = SHAKE128_RATE;
  //     break;
  //   }
  //   SHAKE_256_MODE: {
  //     rate = SHAKE256_RATE;
  //     break;
  //   }
  //   SHA3_256_MODE: {
  //     rate = SHA3_256_RATE;
  //     break;
  //   }
  //   SHA3_512_MODE: {
  //     rate = SHA3_512_RATE;
  //     break;
  //   }
  //   default: rate = 0;
  // }
  size_t vl, avl;  // to be absorbed length in one absorb block
  vuint8m8_t vreg0;

  keccak_init();
  while(mlen >= 0) {
    avl = (mlen > rate)? rate : mlen;
    vl = vsetvl_e8m8(avl);
    vreg0 = vle8_v_u8m8 (m, vl);
    keccak_loadm_v_u8m8(vreg0);
    m += vl;
    mlen = keccak_absorb_one_block(mlen);
  }

}

/*************************************************
* Name:        keccak_squeezeblocks_custom
*
* Description: Squeeze step of Keccak. Squeezes full blocks of r bytes each.
*              Modifies the state. Can be called multiple times to keep
*              squeezing, i.e., is incremental.
*
* Arguments:   - uint8_t *h: pointer to output blocks
*              - size_t nblocks: number of blocks to be squeezed (written to h)
*              - uint64_t *s: pointer to input/output Keccak state
*              - unsigned int r: rate in bytes (e.g., 168 for SHAKE128)
*              - bool squeeze_first: whether squeeze before keccak store.
*                                    if keccak_squeezeblocks_custom closely follows keccak_absorb_custom,
*                                    squeeze_first should be false; else squeeze_first should be true
**************************************************/
void keccak_squeezeblocks_custom(uint8_t *out,
                                 size_t nblocks,
                                 unsigned int rate,
                                 bool squeeze_first)
{
  // unsigned int i;
  // while(nblocks > 0) {
  //   KeccakF1600_StatePermute(s);
  //   for(i=0;i<r/8;i++)
  //     store64(out + 8*i, s[i]);
  //   out += r;
  //   --nblocks;
  // }
  if(nblocks == 0)
    return;

  size_t vl;
  vuint8m8_t vreg0;
  vl = vsetvl_e8m8(rate);

  // 1st block
  if(squeeze_first) {
    keccak_squeeze();
  }
  vreg0 = keccak_stores_v_u8m8();
  vse8_v_u8m8(out, vreg0, vl);
  out += rate;
  nblocks -= 1;

  while(nblocks > 0) {
    keccak_squeeze();
    vreg0 = keccak_stores_v_u8m8();
    vse8_v_u8m8(out, vreg0, vl);
    out += rate;
    nblocks -= 1;
  }
}

/*************************************************
* Name:        shake256_custom
*
* Description: SHAKE256 XOF with non-incremental API
*
* Arguments:   - uint8_t *out:      pointer to output
*              - size_t outlen:     requested output length in bytes
*              - const uint8_t *in: pointer to input
*              - size_t inlen:      length of input in bytes
**************************************************/
void shake256_custom(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen)
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