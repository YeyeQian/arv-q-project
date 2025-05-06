/*
  This file is for Niederreiter encryption
*/

#ifndef ENCRYPT_H
#define ENCRYPT_H

#include "params.h"

#define encrypt CRYPTO_NAMESPACE(encrypt)

void encrypt(unsigned char *, const unsigned char *, unsigned char *);
void syndrome(unsigned char *s, const unsigned char *pk, unsigned char *e);
void gen_e(unsigned char *e);

//Customized Versions
#define encrypt_custom CRYPTO_NAMESPACE(encrypt_custom)
void encrypt_custom(unsigned char *s, const unsigned char *pk, unsigned char *e);
void syndrome_custom(unsigned char *s, const unsigned char *pk, unsigned char *e);
void gen_e_custom(unsigned char *e);
#endif

