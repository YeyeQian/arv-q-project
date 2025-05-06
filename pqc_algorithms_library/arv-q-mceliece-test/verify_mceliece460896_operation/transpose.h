/*
  This file is for matrix transposition
*/

#ifndef TRANSPOSE_H
#define TRANSPOSE_H

#include "params.h"

#define transpose_64x64 CRYPTO_NAMESPACE(transpose_64x64)

#include <stdint.h>

void transpose_64x64(uint64_t *, uint64_t *);

//Customized Versions
#define transpose_64x64_custom CRYPTO_NAMESPACE(transpose_64x64_custom)
void transpose_64x64_custom(uint64_t * out, uint64_t * in);

#endif

