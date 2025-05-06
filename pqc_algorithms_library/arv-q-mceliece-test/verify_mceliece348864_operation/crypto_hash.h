/* This file uses SHAKE256 implemented in the Keccak Code Package */

//#include <libkeccak.a.headers/SimpleFIPS202.h> 

#include <stdint.h>
#include <stddef.h>
#define SHAKE256_RATE 136

typedef struct {
	uint64_t s[25];
} keccak_state;


void shake256(uint8_t* out, size_t outlen, const uint8_t* in, size_t inlen);

#define crypto_hash_32b(out,in,inlen) \
  shake256(out,32,in,inlen)

#define shake(out,outlen,in,inlen) \
  shake256(out,outlen,in,inlen)

//Customzied Versions
void shake256_custom(uint8_t *out, size_t outlen, const uint8_t *in, size_t inlen);

#define crypto_hash_32b_custom(out,in,inlen) \
  shake256_custom(out,32,in,inlen)

#define shake_custom(out,outlen,in,inlen) \
  shake256_custom(out,outlen,in,inlen)