/*
  This file is for public-key generation
*/

#ifndef PK_GEN_H
#define PK_GEN_H

#include "params.h"

#define pk_gen CRYPTO_NAMESPACE(pk_gen)

#include "gf.h"

int pk_gen(unsigned char *, unsigned char *, uint32_t *, int16_t *);

//Customized Versions
#define pk_gen_custom CRYPTO_NAMESPACE(pk_gen_custom)
int pk_gen_custom(unsigned char * pk, unsigned char * sk, uint32_t * perm, int16_t * pi);

#endif

