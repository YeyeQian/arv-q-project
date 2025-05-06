/**
 * @file api.h
 * @brief NIST KEM API used by the HQC_KEM IND-CCA2 scheme
 */

#ifndef API_H
#define API_H
#include <stdint.h>
#define CRYPTO_ALGNAME                      "HQC-192"

#define CRYPTO_SECRETKEYBYTES               4562
#define CRYPTO_PUBLICKEYBYTES               4522
#define CRYPTO_BYTES                        64
#define CRYPTO_CIPHERTEXTBYTES              9042

// As a technicality, the public key is appended to the secret key in order to respect the NIST API.
// Without this constraint, CRYPTO_SECRETKEYBYTES would be defined as 32

int crypto_kem_keypair(unsigned char* pk, unsigned char* sk, unsigned char* sk_seed, unsigned char* pk_seed);
int crypto_kem_enc(unsigned char* ct, unsigned char* ss, const unsigned char* pk, const uint64_t* m, const uint64_t* salt);
int crypto_kem_dec(unsigned char* ss, const unsigned char* ct, const unsigned char* sk);

//Customized Versions
int crypto_kem_keypair_custom(unsigned char *pk, unsigned char *sk, unsigned char* sk_seed, unsigned char* pk_seed);
int crypto_kem_enc_custom(unsigned char *ct, unsigned char *ss, const unsigned char *pk, const uint64_t* m, const uint64_t* salt);
int crypto_kem_dec_custom(unsigned char *ss, const unsigned char *ct, const unsigned char *sk);
#endif
