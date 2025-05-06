/*
  This file is for secret-key generation
*/

#ifndef SK_GEN_H
#define SK_GEN_H

#include "params.h"

#define genpoly_gen CRYPTO_NAMESPACE(genpoly_gen)
#define perm_check CRYPTO_NAMESPACE(perm_check)
#define genpoly_gen_vectorized CRYPTO_NAMESPACE(genpoly_gen_vectorized)

#include "gf.h"

#include <stdint.h>

int genpoly_gen(gf *, gf *);
int perm_check(uint32_t *);
int genpoly_gen_vectorized(gf* out, gf* f);

//Customized Versions
#define genpoly_gen_custom CRYPTO_NAMESPACE(genpoly_gen_custom)
int genpoly_gen_custom(gf *out, gf *f);

#endif

