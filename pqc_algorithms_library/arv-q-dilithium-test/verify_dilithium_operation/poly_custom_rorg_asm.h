#include <stdint.h>
#include "poly.h"
#include "polyvec.h"
#include "symmetric.h"
#include "../apis/custom_inst_api.h"

#define poly_eta_pack_caddq_ntt_asm DILITHIUM_NAMESPACE(_poly_eta_pack_caddq_ntt_asm)
void poly_eta_pack_caddq_ntt_asm(poly* r, uint8_t* packed_addr);

#define polyvecl_macc_invntt_mon2nor_asm DILITHIUM_NAMESPACE(_polyvecl_macc_invntt_mon2nor_asm)
void polyvecl_macc_invntt_mon2nor_asm(poly* r, const polyvecl *a, const polyvecl *b);

#define poly_eta_pack_caddq_asm DILITHIUM_NAMESPACE(_poly_eta_pack_caddq_asm)
void poly_eta_pack_caddq_asm(poly* r, uint8_t* packed_addr);

#define poly_eta_unpack_caddq_ntt_asm DILITHIUM_NAMESPACE(_poly_eta_unpack_caddq_ntt_asm)
void poly_eta_unpack_caddq_ntt_asm(poly* r, const uint8_t* packed_addr);

#define poly_t0_unpack_caddq_ntt_asm DILITHIUM_NAMESPACE(_poly_t0_unpack_caddq_ntt_asm)
void poly_t0_unpack_caddq_ntt_asm(poly* r, const uint8_t* packed_addr);

#define poly_gamma1_sample_caddq_ntt_asm DILITHIUM_NAMESPACE(_poly_gamma1_sample_caddq_ntt_asm)
void poly_gamma1_sample_caddq_ntt_asm(poly* y, poly* z, const uint8_t seed[CRHBYTES], uint16_t nonce);

#define poly_caddq_ntt_asm DILITHIUM_NAMESPACE(_poly_caddq_ntt_asm)
void poly_caddq_ntt_asm(poly* r);

#define poly_pointwisemul_invntt_mon2nor_asm DILITHIUM_NAMESPACE(_poly_pointwisemul_invntt_mon2nor_asm)
void poly_pointwisemul_invntt_mon2nor_asm(poly *z, const poly *a, const poly *b);

#define poly_add_reduce_asm DILITHIUM_NAMESPACE(_poly_add_reduce_asm)
void poly_add_reduce_asm(poly *z, const poly *a, const poly *b);

#define poly_caddq_modsub_reduce_asm DILITHIUM_NAMESPACE(_poly_caddq_modsub_reduce_asm)
void poly_caddq_modsub_reduce_asm(poly *a, const poly *b);

#define poly_add_caddq_asm DILITHIUM_NAMESPACE(_poly_add_caddq_asm)
void poly_add_caddq_asm(poly *c, const poly *a, const poly *b);

#define poly_unpackt1_shiftl_ntt_asm DILITHIUM_NAMESPACE(_poly_unpackt1_shiftl_ntt_asm)
void poly_unpackt1_shiftl_ntt_asm(poly *r, const uint8_t* packed_addr);

#define poly_macc_sub_invntt_mon2nor_asm DILITHIUM_NAMESPACE(_poly_macc_sub_invntt_mon2nor_asm)
void poly_macc_sub_invntt_mon2nor_asm(poly *r, const polyvecl *a, const polyvecl *b, const poly *c);