#ifndef POLY_REORG_H
#define POLY_REORG_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

void poly_eta1_add_ntt(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);
void poly_eta1_add_ntt_asm(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

void polyvec_base_mul_acc_tomont_custom(poly *r,
                                      const polyvec *a,
                                      const polyvec *b);
void polyvec_base_mul_acc_tomont_custom_asm(poly *r,
                                            const polyvec *a,
                                            const polyvec *b);

void polyvec_base_mul_acc_intt_tomont_custom_asm(poly *r,
                                                const polyvec *a,
                                                const polyvec *b);
void polyvec_base_mul_acc_intt_tomont_custom(poly *r,
                                                const polyvec *a,
                                                const polyvec *b);

void poly_eta2_add(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);
void poly_eta2_add_asm(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);
void poly_eta2_addq_add_asm(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce);

void polyvec_ntt_custom_asm(polyvec *r);
#endif

