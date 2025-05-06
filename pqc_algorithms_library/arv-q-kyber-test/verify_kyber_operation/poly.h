#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

/*
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1]
 */
typedef struct{
  int16_t coeffs[KYBER_N];
} poly;

#define poly_compress KYBER_NAMESPACE(_poly_compress)
void poly_compress(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], poly *a);
#define poly_decompress KYBER_NAMESPACE(_poly_decompress)
void poly_decompress(poly *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES]);

#define poly_tobytes KYBER_NAMESPACE(_poly_tobytes)
void poly_tobytes(uint8_t r[KYBER_POLYBYTES], poly *a);
#define poly_frombytes KYBER_NAMESPACE(_poly_frombytes)
void poly_frombytes(poly *r, const uint8_t a[KYBER_POLYBYTES]);

#define poly_frommsg KYBER_NAMESPACE(_poly_frommsg)
void poly_frommsg(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES]);
#define poly_tomsg KYBER_NAMESPACE(_poly_tomsg)
void poly_tomsg(uint8_t msg[KYBER_INDCPA_MSGBYTES], poly *r);

#define poly_getnoise_eta1 KYBER_NAMESPACE(_poly_getnoise_eta1)
void poly_getnoise_eta1(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

#define poly_getnoise_eta2 KYBER_NAMESPACE(_poly_getnoise_eta2)
void poly_getnoise_eta2(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

#define poly_ntt KYBER_NAMESPACE(_poly_ntt)
void poly_ntt(poly *r);
#define poly_invntt_tomont KYBER_NAMESPACE(_poly_invntt_tomont)
void poly_invntt_tomont(poly *r);
#define poly_basemul_montgomery KYBER_NAMESPACE(_poly_basemul_montgomery)
void poly_basemul_montgomery(poly *r, const poly *a, const poly *b);
#define poly_tomont KYBER_NAMESPACE(_poly_tomont)
void poly_tomont(poly *r);

#define poly_reduce KYBER_NAMESPACE(_poly_reduce)
void poly_reduce(poly *r);
#define poly_csubq KYBER_NAMESPACE(_poly_csubq)
void poly_csubq(poly *r);

#define poly_add KYBER_NAMESPACE(_poly_add)
void poly_add(poly *r, const poly *a, const poly *b);
#define poly_sub KYBER_NAMESPACE(_poly_sub)
void poly_sub(poly *r, const poly *a, const poly *b);

#define poly_standard_to_bitreverse_all KYBER_NAMESPACE(_poly_standard_to_bitreverse_all)
void poly_standard_to_bitreverse_all(poly *r);

#define poly_bitreverse_to_standard_all KYBER_NAMESPACE(_poly_bitreverse_to_standard_all)
void poly_bitreverse_to_standard_all(poly *r);

/*
 * Poly functions accelerated with custom instructions
 */
#define poly_compress_custom KYBER_NAMESPACE(_poly_compress_custom)
void poly_compress_custom(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], poly *a);

#define poly_compress_rvv KYBER_NAMESPACE(_poly_compress_rvv)
void poly_compress_rvv(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], poly *a);

#define poly_decompress_custom KYBER_NAMESPACE(_poly_decompress_custom)
void poly_decompress_custom(poly *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES]);

#define poly_tobytes_custom KYBER_NAMESPACE(_poly_tobytes_custom)
void poly_tobytes_custom(uint8_t r[KYBER_POLYBYTES], poly *a);

#define poly_frombytes_custom KYBER_NAMESPACE(_poly_frombytes_custom)
void poly_frombytes_custom(poly *r, const uint8_t a[KYBER_POLYBYTES]);

#define poly_frommsg_custom KYBER_NAMESPACE(_poly_frommsg_custom)
void poly_frommsg_custom(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES]);

#define poly_tomsg_custom KYBER_NAMESPACE(_poly_tomsg_custom)
void poly_tomsg_custom(uint8_t msg[KYBER_INDCPA_MSGBYTES], poly *a);

#define poly_getnoise_eta1_custom KYBER_NAMESPACE(_poly_getnoise_eta1_custom)
void poly_getnoise_eta1_custom(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

#define poly_getnoise_eta2_custom KYBER_NAMESPACE(_poly_getnoise_eta2_custom)
void poly_getnoise_eta2_custom(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

#define poly_add_custom KYBER_NAMESPACE(_poly_add_custom)
void poly_add_custom(poly *r, const poly *a, const poly *b);

#define poly_sub_custom KYBER_NAMESPACE(_poly_sub_custom)
void poly_sub_custom(poly *r, const poly *a, const poly *b);

#define poly_tomont_custom KYBER_NAMESPACE(_poly_tomont_custom)
void poly_tomont_custom(poly *r);

#define poly_basemul_montgomery_custom KYBER_NAMESPACE(_poly_basemul_montgomery_custom)
void poly_basemul_montgomery_custom(poly *r, const poly *a, const poly *b);

#define poly_ntt_custom KYBER_NAMESPACE(_poly_ntt_custom)
void poly_ntt_custom(poly *r);

#define poly_invntt_tomont_custom KYBER_NAMESPACE(_poly_invntt_tomont_custom)
void poly_invntt_tomont_custom(poly *r);

#define poly_mod_add_q KYBER_NAMESPACE(_poly_mod_add_q)
void poly_mod_add_q(poly *a);

#define poly_frommsg_add_custom KYBER_NAMESPACE(_poly_frommsg_add_custom)
void poly_frommsg_add_custom(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES]);

#define poly_sub_tomsg_custom KYBER_NAMESPACE(_poly_sub_tomsg_custom)
void poly_sub_tomsg_custom(uint8_t msg[KYBER_INDCPA_MSGBYTES], poly *a, poly *b);

void poly_basemul_montgomery_custom_redundant(poly* r, const poly* a, const poly* b, const uint8_t buffer_r[128]);
void poly_basemul_montgomery_custom_redundant_max(poly* r, const poly* a, const poly* b, const uint8_t buffer_r[128]);

void poly_basemul_montgomery_custom_shuffling(poly* r, const poly* a, const poly* b, const uint8_t buffer_r[128]);
#endif
