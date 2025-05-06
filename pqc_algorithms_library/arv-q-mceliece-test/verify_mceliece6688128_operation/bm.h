/*
  This file is for the Berlekamp-Massey algorithm
  see http://crypto.stanford.edu/~mironov/cs359/massey.pdf
*/

#ifndef BM_H
#define BM_H

#include "params.h"

#define bm CRYPTO_NAMESPACE(bm)

void bm(gf *, gf *);

//Customized Versions
#define bm_custom CRYPTO_NAMESPACE(bm_custom)
void bm_custom(gf *out, gf *s);

#endif

