/*
  This file is for evaluating a polynomial at one or more field elements
*/

#ifndef ROOT_H
#define ROOT_H

#include "params.h"
#include <stdint.h>

#define eval CRYPTO_NAMESPACE(eval)
#define root CRYPTO_NAMESPACE(root)

#include "gf.h"

gf eval(gf *, gf);
void root(gf *, gf *, gf *);

//Customized Versions
#define root_custom CRYPTO_NAMESPACE(root_custom)
void root_custom(gf *out, gf *f, gf *L);

void root_shuffling_custom(gf *out, gf *f, gf *L, uint16_t shuffle_index[SYS_N]);
#endif

