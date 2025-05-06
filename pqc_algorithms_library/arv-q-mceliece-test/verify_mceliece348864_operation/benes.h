/*
  This file is for Benes network related functions
*/

#ifndef BENES_H
#define BENES_H

#include "params.h"
#include <stdint.h>

#define apply_benes CRYPTO_NAMESPACE(apply_benes)
#define support_gen CRYPTO_NAMESPACE(support_gen)

#include "gf.h"

void apply_benes(unsigned char *, const unsigned char *, int);
void support_gen(gf *, const unsigned char *);

//Customzied Versions
#define apply_benes_custom CRYPTO_NAMESPACE(apply_benes_custom)
void apply_benes_custom(unsigned char * r, const unsigned char * bits, int rev);
#define support_gen_custom CRYPTO_NAMESPACE(support_gen_custom)
void support_gen_custom(gf * s, const unsigned char *c);

//m4 optimized versions
void benes_4096(uint32_t * r32, const unsigned char * bits, int rev);
#endif

