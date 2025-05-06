#ifndef REDUCE_H
#define REDUCE_H
#include "SABER_params.h"
unsigned int montgomery_reduce(unsigned int a, unsigned int b);
unsigned int mod_add(unsigned int a, unsigned int b);
unsigned int mod_sub(unsigned int a, unsigned int b);
#endif
