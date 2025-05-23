#ifndef SPX_FIPS202_H
#define SPX_FIPS202_H

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#define SHAKE128_RATE 168
#define SHAKE256_RATE 136
#define SHA3_256_RATE 136
#define SHA3_512_RATE 72

void shake128_absorb(uint64_t *s, const uint8_t *input, size_t inlen);

void shake128_squeezeblocks(uint8_t *output, size_t nblocks, uint64_t *s);

void shake128_inc_init(uint64_t *s_inc);
void shake128_inc_absorb(uint64_t *s_inc, const uint8_t *input, size_t inlen);
void shake128_inc_finalize(uint64_t *s_inc);
void shake128_inc_squeeze(uint8_t *output, size_t outlen, uint64_t *s_inc);

void shake256_absorb(uint64_t *s, const uint8_t *input, size_t inlen);
void shake256_squeezeblocks(uint8_t *output, size_t nblocks, uint64_t *s);

void shake256_inc_init(uint64_t *s_inc);
void shake256_inc_absorb(uint64_t *s_inc, const uint8_t *input, size_t inlen);
void shake256_inc_finalize(uint64_t *s_inc);
void shake256_inc_squeeze(uint8_t *output, size_t outlen, uint64_t *s_inc);

void shake128(uint8_t *output, size_t outlen,
              const uint8_t *input, size_t inlen);

void shake256(uint8_t *output, size_t outlen,
              const uint8_t *input, size_t inlen);

void sha3_256_inc_init(uint64_t *s_inc);
void sha3_256_inc_absorb(uint64_t *s_inc, const uint8_t *input, size_t inlen);
void sha3_256_inc_finalize(uint8_t *output, uint64_t *s_inc);

void sha3_256(uint8_t *output, const uint8_t *input, size_t inlen);

void sha3_512_inc_init(uint64_t *s_inc);
void sha3_512_inc_absorb(uint64_t *s_inc, const uint8_t *input, size_t inlen);
void sha3_512_inc_finalize(uint8_t *output, uint64_t *s_inc);

void sha3_512(uint8_t *output, const uint8_t *input, size_t inlen);

void keccak_absorb_custom(const uint8_t *m, int mlen, uint32_t rate);
void keccak_squeezeblocks_custom(uint8_t *out, size_t nblocks, unsigned int rate, bool squeeze_first);
void shake128_custom(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);
void shake256_custom(uint8_t *out, int outlen, const uint8_t *in, size_t inlen);
void sha3_256_custom(uint8_t h[32], const uint8_t *in, size_t inlen);
void sha3_512_custom(uint8_t h[64], const uint8_t *in, size_t inlen);

#endif
