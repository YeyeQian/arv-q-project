#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "param.h"
#include "../apis/custom_inst_api.h"

#if (VLEN == 256)
    #define UNPACK_REJ_LOAD_BYTE 24
    #define UNPACK_REJ_ETA_LOAD_BYTE 4
#elif (VLEN == 512)
    #define UNPACK_REJ_LOAD_BYTE 48
    #define UNPACK_REJ_ETA_LOAD_BYTE 8
#elif (VLEN == 1024)
    #define UNPACK_REJ_LOAD_BYTE 96
    #define UNPACK_REJ_ETA_LOAD_BYTE 16
# else
#error "VLEN must be 256/512/1024"
#endif   


typedef struct {
  int32_t coeffs[N];
} poly;

#define poly_print DILITHIUM_NAMESPACE(_poly_print)
void poly_print(const poly* a);

#define poly_copy DILITHIUM_NAMESPACE(_poly_copy)
void poly_copy(poly* dst, const poly* src);

#define poly_chorder DILITHIUM_NAMESPACE(_poly_chorder)
void poly_chorder(poly* a);

#define poly_reduce DILITHIUM_NAMESPACE(_poly_reduce)
void poly_reduce(poly *a);
#define poly_reduce_rvv DILITHIUM_NAMESPACE(_poly_reduce_rvv)
void poly_reduce_rvv(poly *a);

#define poly_caddq DILITHIUM_NAMESPACE(_poly_caddq)
void poly_caddq(poly *a);
#define poly_caddq_rvv DILITHIUM_NAMESPACE(_poly_caddq_rvv)
void poly_caddq_rvv(poly *a);
#define poly_caddq_custom DILITHIUM_NAMESPACE(_poly_caddq_custom)
void poly_caddq_custom(poly *a);

#define poly_freeze DILITHIUM_NAMESPACE(_poly_freeze)
void poly_freeze(poly *a);
#define poly_freeze_rvv DILITHIUM_NAMESPACE(_poly_freeze_rvv)
void poly_freeze_rvv(poly *a);

#define poly_add DILITHIUM_NAMESPACE(_poly_add)
void poly_add(poly *c, const poly *a, const poly *b);
#define poly_add_custom DILITHIUM_NAMESPACE(_poly_add_custom)
void poly_add_custom(poly *c, const poly *a, const poly *b);
#define poly_add_rvv DILITHIUM_NAMESPACE(_poly_add_rvv)
void poly_add_rvv(poly *c, const poly *a, const poly *b);

#define poly_sub DILITHIUM_NAMESPACE(_poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);
#define poly_sub_custom DILITHIUM_NAMESPACE(_poly_sub_custom)
void poly_sub_custom(poly *c, const poly *a, const poly *b);

#define poly_mon2nor_custom DILITHIUM_NAMESPACE(_poly_mon2nor_custom)
void poly_mon2nor_custom(poly *a);

#define poly_shiftl DILITHIUM_NAMESPACE(_poly_shiftl)
void poly_shiftl(poly *a);
#define poly_shiftl_rvv DILITHIUM_NAMESPACE(_poly_shiftl_rvv)
void poly_shiftl_rvv(poly *a);

#define poly_ntt DILITHIUM_NAMESPACE(_poly_ntt)
void poly_ntt(poly *a);
#define poly_ntt_custom DILITHIUM_NAMESPACE(_poly_ntt_custom)
void poly_ntt_custom(poly *a);

#define poly_invntt_tomont DILITHIUM_NAMESPACE(_poly_invntt_tomont)
void poly_invntt_tomont(poly *a);
#define poly_invntt_custom DILITHIUM_NAMESPACE(_poly_invntt_custom)
void poly_invntt_custom(poly *a);

#define poly_pointwise_montgomery DILITHIUM_NAMESPACE(_poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);
#define poly_pointwise_montgomery_custom DILITHIUM_NAMESPACE(_poly_pointwise_montgomery_custom)
void poly_pointwise_montgomery_custom(poly *c, const poly *a, const poly *b);

#define poly_power2round DILITHIUM_NAMESPACE(_poly_power2round)
void poly_power2round(poly *a1, poly *a0, const poly *a);
#define poly_power2round_rvv DILITHIUM_NAMESPACE(_poly_power2round_rvv)
void poly_power2round_rvv(poly *a1, poly *a0, const poly *a);

#define poly_decompose DILITHIUM_NAMESPACE(_poly_decompose)
void poly_decompose(poly *a1, poly *a0, const poly *a);
#define poly_decompose_rvv DILITHIUM_NAMESPACE(_poly_decompose_rvv)
void poly_decompose_rvv(poly *a1, poly *a0, const poly *a);

#define poly_make_hint DILITHIUM_NAMESPACE(_poly_make_hint)
unsigned int poly_make_hint(poly *h, const poly *a0, const poly *a1);
#define poly_make_hint_rvv DILITHIUM_NAMESPACE(_poly_make_hint_rvv)
unsigned int poly_make_hint_rvv(poly *h, const poly *a0, const poly *a1);

#define poly_use_hint DILITHIUM_NAMESPACE(_poly_use_hint)
void poly_use_hint(poly *b, const poly *a, const poly *h);
#define poly_use_hint_rvv DILITHIUM_NAMESPACE(_poly_use_hint_rvv)
void poly_use_hint_rvv(poly *b, const poly *a, const poly *h);

#define poly_chknorm DILITHIUM_NAMESPACE(_poly_chknorm)
int poly_chknorm(const poly *a, int32_t B);
#define poly_chknorm_rvv DILITHIUM_NAMESPACE(_poly_chknorm_rvv)
int poly_chknorm_rvv(const poly *a, int32_t B);

#define poly_uniform DILITHIUM_NAMESPACE(_poly_uniform)
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce);

#define poly_uniform_custom DILITHIUM_NAMESPACE(_poly_uniform_custom)
void poly_uniform_custom(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce);

#define poly_uniform_eta DILITHIUM_NAMESPACE(_poly_uniform_eta)
void poly_uniform_eta(poly *a,
                      const uint8_t seed[SEEDBYTES],
                      uint16_t nonce);

#define poly_uniform_eta_custom DILITHIUM_NAMESPACE(_poly_uniform_eta_custom)
void poly_uniform_eta_custom(poly *a,
                      const uint8_t seed[SEEDBYTES],
                      uint16_t nonce);

#define poly_uniform_gamma1 DILITHIUM_NAMESPACE(_poly_uniform_gamma1)
void poly_uniform_gamma1(poly *a,
                         const uint8_t seed[CRHBYTES],
                         uint16_t nonce);

#define poly_uniform_gamma1_custom DILITHIUM_NAMESPACE(_poly_uniform_gamma1_custom)
void poly_uniform_gamma1_custom(poly *a,
                         const uint8_t seed[CRHBYTES],
                         uint16_t nonce);

#define poly_challenge DILITHIUM_NAMESPACE(_poly_challenge)
void poly_challenge(poly *c, const uint8_t seed[SEEDBYTES]);
#define poly_challenge_custom DILITHIUM_NAMESPACE(_poly_challenge_custom)
void poly_challenge_custom(poly *c, const uint8_t seed[SEEDBYTES]);

#define polyeta_pack DILITHIUM_NAMESPACE(_polyeta_pack)
void polyeta_pack(uint8_t *r, const poly *a);
#define polyeta_pack_custom DILITHIUM_NAMESPACE(_polyeta_pack_custom)
void polyeta_pack_custom(uint8_t *r, const poly *a);
#define polyeta_unpack DILITHIUM_NAMESPACE(_polyeta_unpack)
void polyeta_unpack(poly *r, const uint8_t *a);
#define polyeta_unpack_custom DILITHIUM_NAMESPACE(_polyeta_unpack_custom)
void polyeta_unpack_custom(poly *r, const uint8_t *a);

#define polyt1_pack DILITHIUM_NAMESPACE(_polyt1_pack)
void polyt1_pack(uint8_t *r, const poly *a);
#define polyt1_pack_custom DILITHIUM_NAMESPACE(_polyt1_pack_custom)
void polyt1_pack_custom(uint8_t *r, const poly *a);
#define polyt1_unpack DILITHIUM_NAMESPACE(_polyt1_unpack)
void polyt1_unpack(poly *r, const uint8_t *a);
#define polyt1_unpack_custom DILITHIUM_NAMESPACE(_polyt1_unpack_custom)
void polyt1_unpack_custom(poly *r, const uint8_t *a);

#define polyt0_pack DILITHIUM_NAMESPACE(_polyt0_pack)
void polyt0_pack(uint8_t *r, const poly *a);
#define polyt0_pack_custom DILITHIUM_NAMESPACE(_polyt0_pack_custom)
void polyt0_pack_custom(uint8_t *r, const poly *a);
#define polyt0_unpack DILITHIUM_NAMESPACE(_polyt0_unpack)
void polyt0_unpack(poly *r, const uint8_t *a);
#define polyt0_unpack_custom DILITHIUM_NAMESPACE(_polyt0_unpack_custom)
void polyt0_unpack_custom(poly *r, const uint8_t *a);

#define polyz_pack DILITHIUM_NAMESPACE(_polyz_pack)
void polyz_pack(uint8_t *r, const poly *a);
#define polyz_pack_custom DILITHIUM_NAMESPACE(_polyz_pack_custom)
void polyz_pack_custom(uint8_t *r, const poly *a);
#define polyz_unpack DILITHIUM_NAMESPACE(_polyz_unpack)
void polyz_unpack(poly *r, const uint8_t *a);
#define polyz_unpack_custom DILITHIUM_NAMESPACE(_polyz_unpack_custom)
void polyz_unpack_custom(poly *r, const uint8_t *a);

#define polyw1_pack DILITHIUM_NAMESPACE(_polyw1_pack)
void polyw1_pack(uint8_t *r, const poly *a);

#define polyw1_pack_custom DILITHIUM_NAMESPACE(_polyw1_pack_custom)
void polyw1_pack_custom(uint8_t *r, const poly *a);

/******************************************************************
 *              ASM OPTIMIZED VERSIONS 
*****************************************************************/
#define poly_ntt_custom_asm DILITHIUM_NAMESPACE(_poly_ntt_custom_asm)
void poly_ntt_custom_asm(poly *a);

#define poly_invntt_custom_asm DILITHIUM_NAMESPACE(_poly_invntt_custom_asm)
void poly_invntt_custom_asm(poly *a);

#endif
