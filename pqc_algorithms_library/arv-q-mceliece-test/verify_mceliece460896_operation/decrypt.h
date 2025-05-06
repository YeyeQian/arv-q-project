/*
  This file is for Nieddereiter decryption
*/

#ifndef DECRYPT_H
#define DECRYPT_H

#include "params.h"

#define decrypt CRYPTO_NAMESPACE(decrypt)

int decrypt(unsigned char *, const unsigned char *, const unsigned char *);

//Customized Versions
#define decrypt_custom CRYPTO_NAMESPACE(decrypt_custom)
int decrypt_custom(unsigned char *e, const unsigned char *sk, const unsigned char *c);

#endif

