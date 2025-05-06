/*
  This file is for Niederreiter encryption
*/

#ifndef ENCRYPT_H
#define ENCRYPT_H

#include "params.h"

#define encrypt CRYPTO_NAMESPACE(encrypt)

void encrypt(unsigned char *, const unsigned char *, unsigned char *);

//Customized Versions
#define encrypt_custom CRYPTO_NAMESPACE(encrypt_custom)
void encrypt_custom(unsigned char *s, const unsigned char *pk, unsigned char *e);

#define encrypt_custom_asm CRYPTO_NAMESPACE(encrypt_custom_asm)
void encrypt_custom_asm(unsigned char *s, const unsigned char *pk, unsigned char *e);

#define encrypt_custom_asm_trans CRYPTO_NAMESPACE(encrypt_custom_asm_trans)
void encrypt_custom_asm_trans(unsigned char *s, const unsigned char *pk_trans, unsigned char *e);

#define encrypt_custom_fortest CRYPTO_NAMESPACE(encrypt_custom_fortest)
void encrypt_custom_fortest(unsigned char *s, const unsigned char *pk, unsigned char *e);
#endif

