/*
  This file is for syndrome computation
*/

#ifndef SYND_H
#define SYND_H

#include "params.h"

#define synd CRYPTO_NAMESPACE(synd)

#include "gf.h"

void synd(gf *, gf *, gf *, unsigned char *);

//Customized Versions
#define synd_custom CRYPTO_NAMESPACE(synd_custom)
void synd_custom(gf *out, gf *f, gf *L, unsigned char *r);

#endif

