#include "api.h"
#include "hqc.h"
#include "parameters.h"
#include "parsing.h"
#include "shake_ds.h"
#include "fips202.h"
#include "vector.h"
#include <stdint.h>
#include <string.h>
#ifdef VERBOSE
#include <stdio.h>
#endif

int crypto_kem_keypair_custom(unsigned char *pk, unsigned char *sk, unsigned char* sk_seed, unsigned char* pk_seed) {
    #ifdef VERBOSE
        printf("\n\n\n\n### KEYGEN ###");
    #endif

    hqc_pke_keygen_custom(pk, sk, sk_seed, pk_seed);
    return 0;
}

int crypto_kem_enc_custom(unsigned char *ct, unsigned char *ss, const unsigned char *pk, const uint64_t* m, const uint64_t* salt) {
    #ifdef VERBOSE
        printf("\n\n\n\n### ENCAPS ###");
    #endif

    uint8_t theta[SHAKE256_512_BYTES] = {0};
    uint64_t u[VEC_N_SIZE_64] = {0};
    uint64_t v[VEC_N1N2_SIZE_64] = {0};
    uint8_t d[SHAKE256_512_BYTES] = {0};
    uint8_t mc[VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES+1] = {0};
    uint8_t tmp[VEC_K_SIZE_BYTES + SEED_BYTES + SALT_SIZE_BYTES+1] = {0};

    memcpy(tmp, m, VEC_K_SIZE_BYTES);
    memcpy(tmp + VEC_K_SIZE_BYTES, pk, SEED_BYTES);
    memcpy(tmp + VEC_K_SIZE_BYTES + SEED_BYTES, salt, SALT_SIZE_BYTES);
    tmp[VEC_K_SIZE_BYTES + SEED_BYTES + SALT_SIZE_BYTES] = G_FCT_DOMAIN;
    shake256_custom(theta, SHAKE256_512_BYTES, tmp, sizeof(tmp));

    // Encrypting m
    hqc_pke_encrypt_custom(u, v, m, theta, pk);

    // Computing d
    uint8_t m_pad[VEC_K_SIZE_BYTES + 1] = { 0 };
    memcpy(m_pad, (uint8_t*)m, VEC_K_SIZE_BYTES);
    m_pad[VEC_K_SIZE_BYTES] = H_FCT_DOMAIN;
    shake256_custom(d, SHAKE256_512_BYTES, m_pad, sizeof(m_pad));

    // Computing shared secret
    memcpy(mc, m, VEC_K_SIZE_BYTES);
    memcpy(mc + VEC_K_SIZE_BYTES, u, VEC_N_SIZE_BYTES);
    memcpy(mc + VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES, v, VEC_N1N2_SIZE_BYTES);
    mc[VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES] = K_FCT_DOMAIN;
    shake256_custom(ss, CRYPTO_BYTES, mc, sizeof(mc));

    // Computing ciphertext
    hqc_ciphertext_to_string(ct, u, v, d, salt);

    #ifdef VERBOSE
        printf("\n\npk: "); for(int i = 0 ; i < PUBLIC_KEY_BYTES ; ++i) printf("%02x", pk[i]);
        printf("\n\nm: "); vect_print(m, VEC_K_SIZE_BYTES);
        printf("\n\ntheta: "); for(int i = 0 ; i < SHAKE256_512_BYTES ; ++i) printf("%02x", theta[i]);
        printf("\n\nd: "); for(int i = 0 ; i < SHAKE256_512_BYTES ; ++i) printf("%02x", d[i]);
        printf("\n\nciphertext: "); for(int i = 0 ; i < CIPHERTEXT_BYTES ; ++i) printf("%02x", ct[i]);
        printf("\n\nsecret 1: "); for(int i = 0 ; i < SHARED_SECRET_BYTES ; ++i) printf("%02x", ss[i]);
    #endif

    return 0;
}

int crypto_kem_dec_custom(unsigned char *ss, const unsigned char *ct, const unsigned char *sk) {
    #ifdef VERBOSE
        printf("\n\n\n\n### DECAPS ###");
    #endif

    uint8_t result;
    uint64_t u[VEC_N_SIZE_64] = {0};
    uint64_t v[VEC_N1N2_SIZE_64] = {0};
    uint8_t d[SHAKE256_512_BYTES] = {0};
    uint8_t pk[PUBLIC_KEY_BYTES] = {0};
    uint64_t m[VEC_K_SIZE_64] = {0};
    uint8_t theta[SHAKE256_512_BYTES] = {0};
    uint64_t u2[VEC_N_SIZE_64] = {0};
    uint64_t v2[VEC_N1N2_SIZE_64] = {0};
    uint8_t d2[SHAKE256_512_BYTES] = {0};
    uint8_t mc[VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES+1] = {0};
    uint64_t salt[SALT_SIZE_64] = {0};
    uint8_t tmp[VEC_K_SIZE_BYTES + SALT_SIZE_BYTES + SEED_BYTES+1] = {0};

    // Retrieving u, v and d from ciphertext
    hqc_ciphertext_from_string(u, v , d, salt, ct);

    // Retrieving pk from sk
    memcpy(pk, sk + SEED_BYTES, PUBLIC_KEY_BYTES);

    // Decryting
    hqc_pke_decrypt_custom(m, u, v, sk);

    // Computing theta
    memcpy(tmp, m, VEC_K_SIZE_BYTES);
    memcpy(tmp + VEC_K_SIZE_BYTES, pk, SEED_BYTES);
    memcpy(tmp + VEC_K_SIZE_BYTES + SEED_BYTES, salt, SALT_SIZE_BYTES);
    tmp[VEC_K_SIZE_BYTES + SALT_SIZE_BYTES + SEED_BYTES] = G_FCT_DOMAIN;
    shake256_custom(theta,SHAKE256_512_BYTES,tmp,sizeof(tmp));

    // Encrypting m'
    hqc_pke_encrypt_custom(u2, v2, m, theta, pk);

    // Computing d'
    uint8_t m_pad[VEC_K_SIZE_BYTES + 1] = { 0 };
    memcpy(m_pad, m, VEC_K_SIZE_BYTES);
    m_pad[VEC_K_SIZE_BYTES] = H_FCT_DOMAIN;
    shake256_custom(d2,SHAKE256_512_BYTES,m_pad,sizeof(m_pad));

    // Computing shared secret
    memcpy(mc, m, VEC_K_SIZE_BYTES);
    memcpy(mc + VEC_K_SIZE_BYTES, u, VEC_N_SIZE_BYTES);
    memcpy(mc + VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES, v, VEC_N1N2_SIZE_BYTES);
    mc[VEC_K_SIZE_BYTES + VEC_N_SIZE_BYTES + VEC_N1N2_SIZE_BYTES] = K_FCT_DOMAIN;
    shake256_custom(ss,CRYPTO_BYTES,mc,sizeof(mc));

    // Abort if c != c' or d != d'
    result = vect_compare_custom((uint8_t *)u, (uint8_t *)u2, VEC_N_SIZE_BYTES);
    result |= vect_compare_custom((uint8_t *)v, (uint8_t *)v2, VEC_N1N2_SIZE_BYTES);
    result |= vect_compare_custom(d, d2, SHAKE256_512_BYTES);

    result = (uint8_t) (-((int16_t) result) >> 15);

    for (size_t i = 0 ; i < SHARED_SECRET_BYTES ; i++) {
        ss[i] &= ~result;
    }

    #ifdef VERBOSE
        printf("\n\npk: "); for(int i = 0 ; i < PUBLIC_KEY_BYTES ; ++i) printf("%02x", pk[i]);
        printf("\n\nsk: "); for(int i = 0 ; i < SECRET_KEY_BYTES ; ++i) printf("%02x", sk[i]);
        printf("\n\nciphertext: "); for(int i = 0 ; i < CIPHERTEXT_BYTES ; ++i) printf("%02x", ct[i]);
        printf("\n\nm: "); vect_print(m, VEC_K_SIZE_BYTES);
        printf("\n\ntheta: "); for(int i = 0 ; i < SHAKE256_512_BYTES ; ++i) printf("%02x", theta[i]);
        printf("\n\n\n# Checking Ciphertext- Begin #");
        printf("\n\nu2: "); vect_print(u2, VEC_N_SIZE_BYTES);
        printf("\n\nv2: "); vect_print(v2, VEC_N1N2_SIZE_BYTES);
        printf("\n\nd2: "); for(int i = 0 ; i < SHAKE256_512_BYTES ; ++i) printf("%02x", d2[i]);
        printf("\n\n# Checking Ciphertext - End #\n");
    #endif

    return -(result & 1);
}