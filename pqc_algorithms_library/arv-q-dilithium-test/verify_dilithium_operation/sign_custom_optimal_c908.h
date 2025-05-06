#include <stdint.h>
#include "poly.h"
#include "polyvec.h"
#include "symmetric.h"
#include "../apis/custom_inst_api.h"

#define crypto_sign_keypair_custom_optimal DILITHIUM_NAMESPACE(_crypto_sign_keypair_custom_optimal)
int crypto_sign_keypair_custom_optimal(uint8_t *pk, uint8_t *sk);

#define crypto_sign_signature_custom_optimal DILITHIUM_NAMESPACE(_crypto_sign_signature_custom_optimal)
int crypto_sign_signature_custom_optimal(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk);

#define crypto_sign_custom_optimal DILITHIUM_NAMESPACE(_crypto_sign_custom_optimal)
int crypto_sign_custom_optimal(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk);

#define crypto_sign_verify_custom_optimal DILITHIUM_NAMESPACE(_crypto_sign_verify_custom_optimal)
int crypto_sign_verify_custom_optimal(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk);

#define crypto_sign_open_custom_optimal DILITHIUM_NAMESPACE(_crypto_sign_open_custom_optimal)
int crypto_sign_open_custom_optimal(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *pk);