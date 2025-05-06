#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

typedef struct{
  poly vec[KYBER_K];
} polyvec;

#define polyvec_compress KYBER_NAMESPACE(_polyvec_compress)
void polyvec_compress(uint8_t r[KYBER_POLYVECCOMPRESSEDBYTES], polyvec *a);
#define polyvec_decompress KYBER_NAMESPACE(_polyvec_decompress)
void polyvec_decompress(polyvec *r,
                        const uint8_t a[KYBER_POLYVECCOMPRESSEDBYTES]);

#define polyvec_tobytes KYBER_NAMESPACE(_polyvec_tobytes)
void polyvec_tobytes(uint8_t r[KYBER_POLYVECBYTES], polyvec *a);
#define polyvec_frombytes KYBER_NAMESPACE(_polyvec_frombytes)
void polyvec_frombytes(polyvec *r, const uint8_t a[KYBER_POLYVECBYTES]);

#define polyvec_ntt KYBER_NAMESPACE(_polyvec_ntt)
void polyvec_ntt(polyvec *r);
#define polyvec_invntt_tomont KYBER_NAMESPACE(_polyvec_invntt_tomont)
void polyvec_invntt_tomont(polyvec *r);

#define polyvec_pointwise_acc_montgomery \
        KYBER_NAMESPACE(_polyvec_pointwise_acc_montgomery)
void polyvec_pointwise_acc_montgomery(poly *r,
                                      const polyvec *a,
                                      const polyvec *b);

#define polyvec_reduce KYBER_NAMESPACE(_polyvec_reduce)
void polyvec_reduce(polyvec *r);
#define polyvec_csubq KYBER_NAMESPACE(_polyvec_csubq)
void polyvec_csubq(polyvec *r);

#define polyvec_add KYBER_NAMESPACE(_polyvec_add)
void polyvec_add(polyvec *r, const polyvec *a, const polyvec *b);

#define polyvec_standard_to_bitreverse_all KYBER_NAMESPACE(_polyvec_standard_to_bitreverse_all)
void polyvec_standard_to_bitreverse_all(polyvec *r);

#define polyvec_bitreverse_to_standard_all KYBER_NAMESPACE(_polyvec_bitreverse_to_standard_all)
void polyvec_bitreverse_to_standard_all(polyvec *r);

#define polyvec_compress_custom KYBER_NAMESPACE(_polyvec_compress_custom)
void polyvec_compress_custom(uint8_t r[KYBER_POLYVECCOMPRESSEDBYTES], polyvec *a);

#define polyvec_decompress_custom KYBER_NAMESPACE(_polyvec_decompress_custom)
void polyvec_decompress_custom(polyvec *r,const uint8_t a[KYBER_POLYVECCOMPRESSEDBYTES]);

#define polyvec_tobytes_custom KYBER_NAMESPACE(_polyvec_tobytes_custom)
void polyvec_tobytes_custom(uint8_t r[KYBER_POLYVECBYTES], polyvec *a);

#define polyvec_frombytes_custom KYBER_NAMESPACE(_polyvec_frombytes_custom)
void polyvec_frombytes_custom(polyvec *r, const uint8_t a[KYBER_POLYVECBYTES]);

#define polyvec_ntt_custom KYBER_NAMESPACE(_polyvec_ntt_custom)
void polyvec_ntt_custom(polyvec *r);

#define polyvec_invntt_tomont_custom KYBER_NAMESPACE(_polyvec_invntt_tomont_custom)
void polyvec_invntt_tomont_custom(polyvec *r);

#define polyvec_pointwise_acc_montgomery_custom \
        KYBER_NAMESPACE(_polyvec_pointwise_acc_montgomery_custom)
void polyvec_pointwise_acc_montgomery_custom(poly *r,
                                      const polyvec *a,
                                      const polyvec *b);

#define polyvec_add_custom KYBER_NAMESPACE(_polyvec_add_custom)
void polyvec_add_custom(polyvec *r, const polyvec *a, const polyvec *b);

#define polyvec_mod_add_q KYBER_NAMESPACE(_polyvec_mod_add_q)
void polyvec_mod_add_q(polyvec *r);
#endif
