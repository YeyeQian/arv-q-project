#include "hqc.h"
#include "gf.h"
#include "parameters.h"
#include "parsing.h"
#include "shake_prng.h"
#include "code.h"
#include "vector.h"
#include <stdint.h>
#ifdef VERBOSE
#include <stdio.h>
#endif

void hqc_pke_keygen_custom(unsigned char* pk, unsigned char* sk, unsigned char* sk_seed, unsigned char* pk_seed) {

    uint64_t x[VEC_N_SIZE_64] = {0};
    uint64_t y[VEC_N_SIZE_64] = {0};
    uint64_t h[VEC_N_SIZE_64] = {0};
    uint64_t s[VEC_N_SIZE_64] = {0};
    uint32_t support_y[PARAM_OMEGA_R] = { 0 };//store support for y
    //Unify the following three functions
    //// Create seed_expanders for public key and secret key
    //seedexpander_init(&sk_seedexpander, sk_seed, SEED_BYTES);
    //// Compute secret key
    //vect_set_random_fixed_weight(&sk_seedexpander, x, PARAM_OMEGA);
    //vect_set_random_fixed_weight(&sk_seedexpander, y, PARAM_OMEGA);
    vect_set_random_fixed_weight_xy_custom(x, y, support_y, sk_seed);

    ////Unify the following two functions
    //seedexpander_init(&pk_seedexpander, pk_seed, SEED_BYTES);
    //vect_set_random(&pk_seedexpander, h);
    vect_set_random_h_custom(pk_seed,h);

    vect_mul_cycshift_custom(s,h,support_y,PARAM_OMEGA);
    vect_add_custom(s, x, s, VEC_N_SIZE_64);

    // Parse keys to string
    hqc_public_key_to_string(pk, pk_seed, s);
    hqc_secret_key_to_string(sk, sk_seed, pk);

    #ifdef VERBOSE
        printf("\n\nsk_seed: "); for(int i = 0 ; i < SEED_BYTES ; ++i) printf("%02x", sk_seed[i]);
        printf("\n\nx: "); vect_print(x, VEC_N_SIZE_BYTES);
        printf("\n\ny: "); vect_print(y, VEC_N_SIZE_BYTES);

        printf("\n\npk_seed: "); for(int i = 0 ; i < SEED_BYTES ; ++i) printf("%02x", pk_seed[i]);
        printf("\n\nh: "); vect_print(h, VEC_N_SIZE_BYTES);
        printf("\n\ns: "); vect_print(s, VEC_N_SIZE_BYTES);

        printf("\n\nsk: "); for(int i = 0 ; i < SECRET_KEY_BYTES ; ++i) printf("%02x", sk[i]);
        printf("\n\npk: "); for(int i = 0 ; i < PUBLIC_KEY_BYTES ; ++i) printf("%02x", pk[i]);
    #endif
}

void hqc_pke_encrypt_custom(uint64_t *u, uint64_t *v,const uint64_t *m, unsigned char *theta, const unsigned char *pk) {
    uint64_t h[VEC_N_SIZE_64] = {0};
    uint64_t s[VEC_N_SIZE_64] = {0};
    uint64_t r1[VEC_N_SIZE_64] = {0};
    uint64_t r2[VEC_N_SIZE_64] = {0};
    uint64_t e[VEC_N_SIZE_64] = {0};
    uint64_t tmp1[VEC_N_SIZE_64] = {0};
    uint64_t tmp2[VEC_N_SIZE_64] = {0};
    uint32_t support_r2[PARAM_OMEGA_R] = { 0 };//store support for r2
    // Retrieve h and s from public key
    hqc_public_key_from_string_custom(h, s, pk);

    //Unify the following four functions, results become different with  reference
    //// Create seed_expander from theta
    //seedexpander_init(&seedexpander, theta, SEED_BYTES);
    //// Generate r1, r2 and e
    //vect_set_random_fixed_weight(&seedexpander, r1, PARAM_OMEGA_R);
    //vect_set_random_fixed_weight(&seedexpander, r2, PARAM_OMEGA_R);
    //vect_set_random_fixed_weight(&seedexpander, e, PARAM_OMEGA_E);
    vect_set_random_fixed_weight_r1r2e_custom(r1,r2,support_r2,e,theta);

    // Compute u = r1 + r2.h
    vect_mul_cycshift_custom(u, h, support_r2, PARAM_OMEGA_R);
    vect_add_custom(u, r1, u, VEC_N_SIZE_64);

    // Compute v = m.G by encoding the message
    code_encode_custom(v, m);
    vect_resize(tmp1, PARAM_N, v, PARAM_N1N2);

    // Compute v = m.G + s.r2 + e
    vect_mul_cycshift_custom(tmp2, s, support_r2, PARAM_OMEGA_R);
    vect_add_custom(tmp2, e, tmp2, VEC_N_SIZE_64);
    vect_add_custom(tmp2, tmp1, tmp2, VEC_N_SIZE_64);
    vect_resize(v, PARAM_N1N2, tmp2, PARAM_N);

    #ifdef VERBOSE
        printf("\n\nh: "); vect_print(h, VEC_N_SIZE_BYTES);
        printf("\n\ns: "); vect_print(s, VEC_N_SIZE_BYTES);
        printf("\n\nr1: "); vect_print(r1, VEC_N_SIZE_BYTES);
        printf("\n\nr2: "); vect_print(r2, VEC_N_SIZE_BYTES);
        printf("\n\ne: "); vect_print(e, VEC_N_SIZE_BYTES);
        printf("\n\ntmp2: "); vect_print(tmp2, VEC_N_SIZE_BYTES);

        printf("\n\nu: "); vect_print(u, VEC_N_SIZE_BYTES);
        printf("\n\nv: "); vect_print(v, VEC_N1N2_SIZE_BYTES);
    #endif
}

void hqc_pke_decrypt_custom(uint64_t *m, const uint64_t *u, const uint64_t *v, const unsigned char *sk) {
    uint64_t x[VEC_N_SIZE_64] = {0};
    uint64_t y[VEC_N_SIZE_64] = {0};
    uint64_t tmp1[VEC_N_SIZE_64] = {0};
    uint64_t tmp2[VEC_N_SIZE_64] = {0};
    uint32_t support_y[PARAM_OMEGA_R] = { 0 };
    // Retrieve x, y, pk from secret key
    vect_set_random_fixed_weight_xy_custom(x,y,support_y,sk);

    // Compute v - u.y
    vect_resize(tmp1, PARAM_N, v, PARAM_N1N2);
    vect_mul_cycshift_custom(tmp2, u, support_y, PARAM_OMEGA);
    vect_add_custom(tmp2, tmp1, tmp2, VEC_N_SIZE_64);

    #ifdef VERBOSE
        printf("\n\nu: "); vect_print(u, VEC_N_SIZE_BYTES);
        printf("\n\nv: "); vect_print(v, VEC_N1N2_SIZE_BYTES);
        printf("\n\ny: "); vect_print(y, VEC_N_SIZE_BYTES);
        printf("\n\nv - u.y: "); vect_print(tmp2, VEC_N_SIZE_BYTES);
    #endif

    // Compute m by decoding v - u.y
    code_decode_custom(m, tmp2);
}