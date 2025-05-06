#ifndef SIGN_H
#define SIGN_H

#include <stddef.h>
#include <stdint.h>
#include "param.h"
#include "polyvec.h"
#include "poly.h"

#define challenge DILITHIUM_NAMESPACE(_challenge)
void challenge(poly *c, const uint8_t seed[SEEDBYTES]);

#define crypto_sign_keypair DILITHIUM_NAMESPACE(_keypair)
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk);
#define crypto_sign_keypair_custom DILITHIUM_NAMESPACE(_keypair_custom)
int crypto_sign_keypair_custom(uint8_t *pk, uint8_t *sk);
#define crypto_sign_keypair_custom_rorg DILITHIUM_NAMESPACE(_keypair_custom_rorg)
int crypto_sign_keypair_custom_rorg(uint8_t *pk, uint8_t *sk);

#define crypto_sign_signature DILITHIUM_NAMESPACE(_signature)
int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);
#define crypto_sign_signature_custom DILITHIUM_NAMESPACE(_signature_custom)
int crypto_sign_signature_custom(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);
#define crypto_sign_signature_custom_rorg DILITHIUM_NAMESPACE(_signature_custom_rorg)
int crypto_sign_signature_custom_rorg(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);

#define crypto_sign DILITHIUM_NAMESPACE(_sign)
int crypto_sign(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);
#define crypto_sign_custom DILITHIUM_NAMESPACE(_sign_custom)
int crypto_sign_custom(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);
#define crypto_sign_custom_rorg DILITHIUM_NAMESPACE(_sign_custom_rorg)
int crypto_sign_custom_rorg(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);

#define crypto_sign_verify DILITHIUM_NAMESPACE(_verify)
int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk);
#define crypto_sign_verify_custom DILITHIUM_NAMESPACE(_verify_custom)
int crypto_sign_verify_custom(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk);
#define crypto_sign_verify_custom_rorg DILITHIUM_NAMESPACE(_verify_custom_rorg)
int crypto_sign_verify_custom_rorg(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk);

#define crypto_sign_open DILITHIUM_NAMESPACE(_open)
int crypto_sign_open(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);
#define crypto_sign_open_custom DILITHIUM_NAMESPACE(_open_custom)
int crypto_sign_open_custom(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);
#define crypto_sign_open_custom_rorg DILITHIUM_NAMESPACE(_open_custom_rorg)
int crypto_sign_open_custom_rorg(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);

//////////////////////////////////////////////////////
//Functions below are for consistency comparison only
/////////////////////////////////////////////////////
#define crypto_sign_keypair_custom_fortest DILITHIUM_NAMESPACE(_keypair_custom_fortest)
int crypto_sign_keypair_custom_fortest(uint8_t *pk, uint8_t *sk);

#define crypto_sign_keypair_custom_rorg_fortest DILITHIUM_NAMESPACE(_keypair_custom_rorg_fortest)
int crypto_sign_keypair_custom_rorg_fortest(uint8_t *pk, uint8_t *sk);

#define crypto_sign_signature_custom_fortest DILITHIUM_NAMESPACE(_signature_custom_fortest)
int crypto_sign_signature_custom_fortest(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);

#define crypto_sign_signature_custom_rorg_fortest DILITHIUM_NAMESPACE(_signature_custom_rorg_fortest)
int crypto_sign_signature_custom_rorg_fortest(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const uint8_t *sk);

#define crypto_sign_custom_fortest DILITHIUM_NAMESPACE(_sign_custom_fortest)
int crypto_sign_custom_fortest(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);

#define crypto_sign_custom_rorg_fortest DILITHIUM_NAMESPACE(_sign_custom_rorg_fortest)
int crypto_sign_custom_rorg_fortest(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const uint8_t *sk);

#define crypto_sign_verify_custom_fortest DILITHIUM_NAMESPACE(_verify_custom_fortest)
int crypto_sign_verify_custom_fortest(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk);

#define crypto_sign_verify_custom_rorg_fortest DILITHIUM_NAMESPACE(_verify_custom_rorg_fortest)
int crypto_sign_verify_custom_rorg_fortest(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const uint8_t *pk);

#define crypto_sign_open_custom_fortest DILITHIUM_NAMESPACE(_open_custom_fortest)
int crypto_sign_open_custom_fortest(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);

#define crypto_sign_open_custom_rorg_fortest DILITHIUM_NAMESPACE(_open_custom_rorg_fortest)
int crypto_sign_open_custom_rorg_fortest(uint8_t *m, size_t *mlen,
                     const uint8_t *sm, size_t smlen,
                     const uint8_t *pk);
#endif
