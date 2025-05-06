#ifndef KEM_H
#define KEM_H

#include "params.h"

#define crypto_kem_keypair KYBER_NAMESPACE(_keypair)
int crypto_kem_keypair(unsigned char *pk, unsigned char *sk);

#define crypto_kem_enc KYBER_NAMESPACE(_enc)
int crypto_kem_enc(unsigned char *ct,
                   unsigned char *ss,
                   const unsigned char *pk);

#define crypto_kem_dec KYBER_NAMESPACE(_dec)
int crypto_kem_dec(unsigned char *ss,
                   const unsigned char *ct,
                   const unsigned char *sk);

#define crypto_kem_keypair_custom KYBER_NAMESPACE(_keypair_custom)
int crypto_kem_keypair_custom(unsigned char *pk, unsigned char *sk);

#define crypto_kem_enc_custom KYBER_NAMESPACE(_enc_custom)
int crypto_kem_enc_custom(unsigned char *ct,
                   unsigned char *ss,
                   const unsigned char *pk);

#define crypto_kem_dec_custom KYBER_NAMESPACE(_dec_custom)
int crypto_kem_dec_custom(unsigned char *ss,
                   const unsigned char *ct,
                   const unsigned char *sk);

#define crypto_kem_keypair_custom_asm KYBER_NAMESPACE(_keypair_custom_asm)
int crypto_kem_keypair_custom_asm(unsigned char *pk, unsigned char *sk);

#define crypto_kem_enc_custom_asm KYBER_NAMESPACE(_enc_custom_asm)
int crypto_kem_enc_custom_asm(unsigned char *ct,
                   unsigned char *ss,
                   const unsigned char *pk);

#define crypto_kem_dec_custom_asm KYBER_NAMESPACE(_dec_custom_asm)
int crypto_kem_dec_custom_asm(unsigned char *ss,
                   const unsigned char *ct,
                   const unsigned char *sk);
#endif
