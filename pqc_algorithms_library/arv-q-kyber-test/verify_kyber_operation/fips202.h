#ifndef FIPS202_H
#define FIPS202_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#define FIPS202_NAMESPACE(s) pqcrystals_fips202_ref##s

#define SHAKE128_RATE 168
#define SHAKE256_RATE 136
#define SHA3_256_RATE 136
#define SHA3_512_RATE 72

typedef struct {
  uint64_t s[25];
} keccak_state;

#define shake128_absorb FIPS202_NAMESPACE(_shake128_absorb)
void shake128_absorb(keccak_state *state, const uint8_t *in, size_t inlen);
#define shake128_squeezeblocks FIPS202_NAMESPACE(_shake128_squeezeblocks)
void shake128_squeezeblocks(uint8_t *out, size_t nblocks, keccak_state *state);

#define shake256_absorb FIPS202_NAMESPACE(_shake256_absorb)
void shake256_absorb(keccak_state *state, const uint8_t *in, size_t inlen);
#define shake256_squeezeblocks FIPS202_NAMESPACE(_shake256_squeezeblocks)
void shake256_squeezeblocks(uint8_t *out, size_t nblocks,  keccak_state *state);
#define shake128 FIPS202_NAMESPACE(_shake128)
void shake128(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);
#define shake256 FIPS202_NAMESPACE(_shake256)
void shake256(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);
#define sha3_256 FIPS202_NAMESPACE(_sha3_256)
void sha3_256(uint8_t h[32], const uint8_t *in, size_t inlen);
#define sha3_512 FIPS202_NAMESPACE(_sha3_512)
void sha3_512(uint8_t h[64], const uint8_t *in, size_t inlen);


#define keccak_absorb_custom FIPS202_NAMESPACE(_keccak_absorb_custom)
void keccak_absorb_custom(const uint8_t *m, int mlen, uint32_t rate);
#define keccak_squeezeblocks_custom FIPS202_NAMESPACE(_keccak_squeezeblocks_custom)
void keccak_squeezeblocks_custom(uint8_t *out, size_t nblocks, unsigned int rate, bool squeeze_first);
#define shake128_custom FIPS202_NAMESPACE(_shake128_custom)
void shake128_custom(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);
#define shake256_custom FIPS202_NAMESPACE(_shake256_custom)
void shake256_custom(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);
#define sha3_256_custom FIPS202_NAMESPACE(_sha3_256_custom)
void sha3_256_custom(uint8_t h[32], const uint8_t *in, size_t inlen);
#define sha3_512_custom FIPS202_NAMESPACE(_sha3_512_custom)
void sha3_512_custom(uint8_t h[64], const uint8_t *in, size_t inlen);

#endif