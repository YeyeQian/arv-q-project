#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "param.h"

void chorder(int32_t a[N]);

#define zetas_pos DILITHIUM_NAMESPACE(_zetas_pos)
extern const int32_t zetas_pos[256];

#define zetas_inv_in_order DILITHIUM_NAMESPACE(_zetas_inv_in_order)
extern const int32_t zetas_inv_in_order[256];

#define tree_byteoffset DILITHIUM_NAMESPACE(_tree_byteoffset)
extern const uint32_t tree_byteoffset[256];

#define tree DILITHIUM_NAMESPACE(_tree)
extern const uint32_t tree[256];

#define ntt DILITHIUM_NAMESPACE(_ntt)
void ntt(int32_t a[N]);

#define ntt_cg_custom_zeta_compact DILITHIUM_NAMESPACE(_ntt_cg_custom_zeta_compact)
void ntt_cg_custom_zeta_compact(int32_t a[N]);

#define ntt_cg_custom_reorder DILITHIUM_NAMESPACE(_ntt_cg_custom_reorder)
void ntt_cg_custom_reorder(int32_t a[N]);

#define ntt_cg_custom_reorder_asm DILITHIUM_NAMESPACE(_ntt_cg_custom_reorder_asm)
void ntt_cg_custom_reorder_asm(int32_t a[N]);

#define invntt_tomont DILITHIUM_NAMESPACE(_invntt_tomont)
void invntt_tomont(int32_t a[N]);

#define invntt_cg_custom_zeta_compact DILITHIUM_NAMESPACE(_invntt_cg_custom_zeta_compact)
void invntt_cg_custom_zeta_compact(int32_t a[N]);

#define invntt_cg_custom_reorder DILITHIUM_NAMESPACE(_invntt_cg_custom_reorder)
void invntt_cg_custom_reorder(int32_t a[N]);

#define invntt_cg_custom_reorder_asm DILITHIUM_NAMESPACE(_invntt_cg_custom_reorder_asm)
void invntt_cg_custom_reorder_asm(int32_t a[N]);

#define invntt_CG DILITHIUM_NAMESPACE(_invntt_CG)
void invntt_CG(int32_t a[N]);

#define GET_ZETA_MODE(SAME_TIMES, REPEAT_TIMES) ((SAME_TIMES&0xff) | (REPEAT_TIMES<<8))

#endif
