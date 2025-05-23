#ifndef FIPS202_H
#define FIPS202_H

#include <stdint.h>
#include <stddef.h>

#define SHAKE128_RATE 168
#define SHAKE256_RATE 136
#define SHA3_256_RATE 136
#define SHA3_512_RATE 72

void shake128(unsigned char *output, unsigned long long outlen, const unsigned char *input, unsigned long long inlen);
void sha3_256(unsigned char *output, const unsigned char *input, unsigned long long inlen);
void sha3_512(unsigned char *output, const unsigned char *input, unsigned long long inlen);

//Customized version
void shake128_custom(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);
void sha3_256_custom(uint8_t h[32], const uint8_t *in, size_t inlen);
void sha3_512_custom(uint8_t *h, const uint8_t *in, size_t inlen);
#endif
