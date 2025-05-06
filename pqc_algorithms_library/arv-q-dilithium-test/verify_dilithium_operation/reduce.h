#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "param.h"

#define MONT -4186625 // 2^32 % Q
#define QINV 58728449 // q^(-1) mod 2^32
#define QINV_HW 4236238847 // -q^(-1) mod 2^32

#define montgomery_reduce DILITHIUM_NAMESPACE(_montgomery_reduce)
int32_t montgomery_reduce(int64_t a);

#define reduce32 DILITHIUM_NAMESPACE(_reduce32)
int32_t reduce32(int32_t a);

#define caddq DILITHIUM_NAMESPACE(_caddq)
int32_t caddq(int32_t a);

#define freeze DILITHIUM_NAMESPACE(_freeze)
int32_t freeze(int32_t a);

#define central_reduce_oneQ DILITHIUM_NAMESPACE(_central_reduce_oneQ)
int32_t central_reduce_oneQ(int32_t a);

#define mod_div2 DILITHIUM_NAMESPACE(_mod_div2)
int32_t mod_div2(int32_t a);

#endif
