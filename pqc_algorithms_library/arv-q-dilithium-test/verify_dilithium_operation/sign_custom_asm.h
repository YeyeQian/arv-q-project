#include <stddef.h>
#include <stdint.h>
#include "param.h"
#include "polyvec.h"
#include "poly.h"

#define crypto_sign_keypair_custom_asm DILITHIUM_NAMESPACE(_crypto_sign_keypair_custom_asm)
int crypto_sign_keypair_custom_asm(uint8_t *pk, uint8_t *sk);

#define crypto_sign_signature_custom_asm DILITHIUM_NAMESPACE(_crypto_sign_signature_custom_asm)
int crypto_sign_signature_custom_asm(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk);

#define crypto_sign_custom_asm DILITHIUM_NAMESPACE(_crypto_sign_custom_asm)
int crypto_sign_custom_asm(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk);

#define crypto_sign_verify_custom_asm DILITHIUM_NAMESPACE(_crypto_sign_verify_custom_asm)
int crypto_sign_verify_custom_asm(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk);

#define crypto_sign_open_custom_asm DILITHIUM_NAMESPACE(_crypto_sign_open_custom_asm)
int crypto_sign_open_custom_asm(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *pk);

//////////////////////////////////////////////////////
//Functions below are for consistency comparison only
/////////////////////////////////////////////////////

#define crypto_sign_keypair_custom_asm_fortest DILITHIUM_NAMESPACE(_crypto_sign_keypair_custom_asm_fortest)
int crypto_sign_keypair_custom_asm_fortest(uint8_t *pk, uint8_t *sk);

#define crypto_sign_signature_custom_asm_fortest DILITHIUM_NAMESPACE(_crypto_sign_signature_custom_asm_fortest)
int crypto_sign_signature_custom_asm_fortest(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk);

#define crypto_sign_custom_asm_fortest DILITHIUM_NAMESPACE(_crypto_sign_custom_asm_fortest)
int crypto_sign_custom_asm_fortest(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk);

#define crypto_sign_verify_custom_asm_fortest DILITHIUM_NAMESPACE(_crypto_sign_verify_custom_asm_fortest)
int crypto_sign_verify_custom_asm_fortest(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk);

#define crypto_sign_open_custom_asm_fortest DILITHIUM_NAMESPACE(_crypto_sign_open_custom_asm_fortest)
int crypto_sign_open_custom_asm_fortest(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *pk);