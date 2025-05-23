#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "params.h"

#define MONT 2285 // 2^16 mod q
#define QINV 62209 // q^-1 mod 2^16
#define QINV_HW 3327 // -q^-1 mod 2^16

#define montgomery_reduce KYBER_NAMESPACE(_montgomery_reduce)
int16_t montgomery_reduce(int32_t a);

#define barrett_reduce KYBER_NAMESPACE(_barrett_reduce)
int16_t barrett_reduce(int16_t a);

#define csubq KYBER_NAMESPACE(_csubq)
int16_t csubq(int16_t x);

#define mod_div2 KYBER_NAMESPACE(_mod_div2)
int16_t mod_div2(int16_t a);

#define central_reduce_oneQ KYBER_NAMESPACE(_central_reduce_oneQ)
int16_t central_reduce_oneQ(int16_t a);

#endif
