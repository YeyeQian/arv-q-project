#ifndef NTT_H
#define NTT_H
#include <stdint.h>
#include "SABER_params.h"
extern const int32_t zetas[SABER_N];
extern const int32_t inv_zetas[SABER_N];
extern const uint32_t tree_byteoffset[256];
extern const uint32_t tree[256];
extern const int32_t inv_zetas_inorder[SABER_N];

void ntt(int32_t a[SABER_N]);
void invntt_tomont(int32_t a[SABER_N]);
// int bitreverse(int i, int m);
// void get_invzeta_normorder();

//Customized version
void ntt_cg_custom_reorder(int32_t r[SABER_N]);
void ntt_cg_custom_reorder_asm(int32_t r[SABER_N]);
void invntt_cg_custom_reorder(int32_t r[SABER_N]);
void invntt_cg_custom_reorder_asm(int32_t r[SABER_N]);

#define GET_ZETA_MODE(SAME_TIMES, REPEAT_TIMES) ((SAME_TIMES&0xff) | (REPEAT_TIMES<<8))
#endif