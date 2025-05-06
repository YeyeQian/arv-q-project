#ifndef HQC_H
#define HQC_H

/**
 * @file hqc.h
 * @brief Functions of the HQC_PKE IND_CPA scheme
 */

#include <stdint.h>

void hqc_pke_keygen(unsigned char* pk, unsigned char* sk, unsigned char* sk_seed, unsigned char* pk_seed);
void hqc_pke_encrypt(uint64_t* u, uint64_t* v, const uint64_t* m, unsigned char* theta, const unsigned char* pk);
void hqc_pke_decrypt(uint64_t* m, const uint64_t* u, const uint64_t* v, const unsigned char* sk);

//Customized Versions
void hqc_pke_keygen_custom(unsigned char* pk, unsigned char* sk, unsigned char* sk_seed, unsigned char* pk_seed);
void hqc_pke_encrypt_custom(uint64_t *u, uint64_t *v,const uint64_t *m, unsigned char *theta, const unsigned char *pk);
void hqc_pke_decrypt_custom(uint64_t *m, const uint64_t *u, const uint64_t *v, const unsigned char *sk);
#endif