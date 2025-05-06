#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"
#include "../apis/custom_inst_api.h"

#define zetas KYBER_NAMESPACE(_zetas)
extern const int16_t zetas[128];

#if (VLEN == 256)
#define zetas_cg KYBER_NAMESPACE(_zetas_cg)
extern const int16_t zetas_cg[11][16];
#elif (VLEN == 512)
#define zetas_cg KYBER_NAMESPACE(_zetas_cg)
extern const int16_t zetas_cg[8][32];
#elif (VLEN == 1024)
#define zetas_cg KYBER_NAMESPACE(_zetas_cg)
extern const int16_t zetas_cg[7][64];
# else
#error "VLEN must be 256/512/1024"
#endif

#define GET_ZETA_MODE(SAME_TIMES, REPEAT_TIMES) ((SAME_TIMES&0xff) | (REPEAT_TIMES<<8))

#define zetas_inv KYBER_NAMESPACE(_zetas_inv)
extern const int16_t zetas_inv[128];

#define zetas_inv_in_order KYBER_NAMESPACE(_zetas_inv_in_order)
extern const int16_t zetas_inv_in_order[128];

//for SCA usage
#define zetas_inv_in_order_stage1 KYBER_NAMESPACE(_zetas_inv_in_order_stage1)
extern const int16_t zetas_inv_in_order_stage1[128];
#define zetas_inv_in_order_stage2 KYBER_NAMESPACE(_zetas_inv_in_order_stage2)
extern const int16_t zetas_inv_in_order_stage2[128];
#define zetas_inv_in_order_stage3 KYBER_NAMESPACE(_zetas_inv_in_order_stage3)
extern const int16_t zetas_inv_in_order_stage3[128];
#define zetas_inv_in_order_stage4 KYBER_NAMESPACE(_zetas_inv_in_order_stage4)
extern const int16_t zetas_inv_in_order_stage4[128];
#define zetas_inv_in_order_stage5 KYBER_NAMESPACE(_zetas_inv_in_order_stage5)
extern const int16_t zetas_inv_in_order_stage5[128];
#define zetas_inv_in_order_stage6 KYBER_NAMESPACE(_zetas_inv_in_order_stage6)
extern const int16_t zetas_inv_in_order_stage6[128];
#define zetas_inv_in_order_stage7 KYBER_NAMESPACE(_zetas_inv_in_order_stage7)
extern const int16_t zetas_inv_in_order_stage7[128];

#define zetas_basemul_cg KYBER_NAMESPACE(_zetas_basemul_cg)
extern const int16_t zetas_basemul_cg[128];

#define tree_128_NWC KYBER_NAMESPACE(_tree_128_NWC)
extern const uint16_t tree_128_NWC[128];
#define byteoffset_even KYBER_NAMESPACE(_byteoffset_even)
extern const uint16_t byteoffset_even[128];
#define byteoffset_total KYBER_NAMESPACE(_byteoffset_total)
extern const uint16_t byteoffset_total[256];

#define ntt KYBER_NAMESPACE(_ntt)
void ntt(int16_t poly[256]);

#define invntt KYBER_NAMESPACE(_invntt)
void invntt(int16_t poly[256]);

#define invntt_without_post_process KYBER_NAMESPACE(_invntt_without_post_process)
void invntt_without_post_process(int16_t r[KYBER_N]);

#define basemul KYBER_NAMESPACE(_basemul)
void basemul(int16_t r[2],
             const int16_t a[2],
             const int16_t b[2],
             int16_t zeta);

void copy_array(int16_t a[KYBER_N], int16_t temp[KYBER_N]);

#define ntt_CG KYBER_NAMESPACE(_ntt_CG)
void ntt_CG(int16_t poly[256]);

#define invntt_CG KYBER_NAMESPACE(_invntt_CG)
void invntt_CG(int16_t a[KYBER_N]);

#define bitreverse_to_standard_all KYBER_NAMESPACE(_bitreverse_to_standard_all)
void bitreverse_to_standard_all(int16_t a[KYBER_N]);

#define bitreverse_cg_to_standard_all KYBER_NAMESPACE(_bitreverse_cg_to_standard_all)
void bitreverse_cg_to_standard_all(int16_t a[KYBER_N]);

#define standard_to_bitreverse_all KYBER_NAMESPACE(_standard_to_bitreverse_all)
void standard_to_bitreverse_all(int16_t a[KYBER_N]);

#if (VLEN == 256)
#define ntt_cg_custom_asm KYBER_NAMESPACE(_ntt_cg_custom_asm)
void ntt_cg_custom_asm(int16_t poly[256]);

#define ntt_cg_custom_reorder_asm KYBER_NAMESPACE(_ntt_cg_custom_reorder_asm)
void ntt_cg_custom_reorder_asm(int16_t r[256]);

#define ntt_cg_custom KYBER_NAMESPACE(_ntt_cg_custom)
void ntt_cg_custom(int16_t poly[256]);

#define ntt_cg_custom_zeta_compact KYBER_NAMESPACE(_ntt_cg_custom_zeta_compact)
void ntt_cg_custom_zeta_compact(int16_t poly[256]);

#define ntt_cg_custom_reorder KYBER_NAMESPACE(_ntt_cg_custom_reorder)
void ntt_cg_custom_reorder(int16_t poly[256]);

#define invntt_cg_zeta_compact KYBER_NAMESPACE(_invntt_cg_zeta_compact)
void invntt_cg_zeta_compact(int16_t a[KYBER_N]);

#define invntt_cg_custom_reorder KYBER_NAMESPACE(_invntt_cg_custom_reorder)
void invntt_cg_custom_reorder(int16_t a[KYBER_N]);
#elif (VLEN == 512)
#define ntt_cg_custom_asm KYBER_NAMESPACE(_ntt_cg_custom_asm)
void ntt_cg_custom_asm(int16_t poly[256]);

#define ntt_cg_custom_reorder_asm KYBER_NAMESPACE(_ntt_cg_custom_reorder_asm)
void ntt_cg_custom_reorder_asm(int16_t poly[256]);

#define ntt_cg_custom_zeta_compact KYBER_NAMESPACE(_ntt_cg_custom_zeta_compact)
void ntt_cg_custom_zeta_compact(int16_t poly[256]);

#define ntt_cg_custom_reorder KYBER_NAMESPACE(_ntt_cg_custom_reorder)
void ntt_cg_custom_reorder(int16_t poly[256]);

#define invntt_cg_custom_reorder KYBER_NAMESPACE(_invntt_cg_custom_reorder)
void invntt_cg_custom_reorder(int16_t a[KYBER_N]);

#define invntt_cg_custom_reorder_asm KYBER_NAMESPACE(_invntt_cg_custom_reorder_asm)
void invntt_cg_custom_reorder_asm(int16_t r[KYBER_N]);
#elif (VLEN == 1024)
#define ntt_cg_custom_reorder KYBER_NAMESPACE(_ntt_cg_custom_reorder)
void ntt_cg_custom_reorder(int16_t r[256]);

#define ntt_cg_custom_reorder_asm KYBER_NAMESPACE(_ntt_cg_custom_reorder_asm)
void ntt_cg_custom_reorder_asm(int16_t r[256]);

#define invntt_cg_custom_reorder KYBER_NAMESPACE(_invntt_cg_custom_reorder)
void invntt_cg_custom_reorder(int16_t r[KYBER_N]);
# else
#error "VLEN must be 256/512/1024"
#endif

#endif
