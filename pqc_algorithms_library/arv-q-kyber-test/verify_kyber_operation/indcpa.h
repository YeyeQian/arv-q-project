#ifndef INDCPA_H
#define INDCPA_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "../apis/custom_inst_api.h"

#if (VLEN == 256)
    #define UNPACK_REJ_LOAD_BYTE_KYBER 24
#elif (VLEN == 512)
    #define UNPACK_REJ_LOAD_BYTE_KYBER 48
#elif (VLEN == 1024)
    #define UNPACK_REJ_LOAD_BYTE_KYBER 96
# else
#error "VLEN must be 256/512/1024"
#endif

unsigned int rej_uniform(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen);
                                
unsigned int rej_uniform_custom(int16_t *r,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen);

#define gen_matrix KYBER_NAMESPACE(_gen_matrix)
void gen_matrix(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed);
#define indcpa_keypair KYBER_NAMESPACE(_indcpa_keypair)
void indcpa_keypair(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);

#define indcpa_enc KYBER_NAMESPACE(_indcpa_enc)
void indcpa_enc(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES]);

#define indcpa_dec KYBER_NAMESPACE(_indcpa_dec)
void indcpa_dec(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);

void gen_matrix_custom(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed);

#define indcpa_keypair_custom KYBER_NAMESPACE(_indcpa_keypair_custom)
void indcpa_keypair_custom(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);

#define indcpa_enc_custom KYBER_NAMESPACE(_indcpa_enc_custom)
void indcpa_enc_custom(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES]);

#define indcpa_dec_custom KYBER_NAMESPACE(_indcpa_dec_custom)
void indcpa_dec_custom(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);


#define indcpa_keypair_custom_for_test KYBER_NAMESPACE(_indcpa_keypair_custom_for_test)
void indcpa_keypair_custom_for_test(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);

#define indcpa_enc_custom_for_test KYBER_NAMESPACE(_indcpa_enc_custom_for_test)
void indcpa_enc_custom_for_test(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES]);

#define indcpa_dec_custom_for_test KYBER_NAMESPACE(_indcpa_dec_custom_for_test)
void indcpa_dec_custom_for_test(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);

void indcpa_keypair_custom_asm(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                              uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);

void indcpa_enc_custom_asm(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES]);

void indcpa_dec_custom_asm(uint8_t m[KYBER_INDCPA_MSGBYTES],
                          const uint8_t c[KYBER_INDCPA_BYTES],
                          const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES]);
#endif
