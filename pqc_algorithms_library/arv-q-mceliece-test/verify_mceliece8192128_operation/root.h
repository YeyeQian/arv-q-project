/*
  This file is for evaluating a polynomial at one or more field elements
*/

#ifndef ROOT_H
#define ROOT_H

#include "params.h"

#define eval CRYPTO_NAMESPACE(eval)
#define root CRYPTO_NAMESPACE(root)

#include "gf.h"

gf eval(gf *, gf);
void root(gf *, gf *, gf *);

//Customized Versions
#define root_custom CRYPTO_NAMESPACE(root_custom)
void root_custom(gf *out, gf *f, gf *L);

#endif

