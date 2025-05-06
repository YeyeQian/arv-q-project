/*
  This file is for functions for field arithmetic
*/

#ifndef GF_H
#define GF_H

#include "params.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include <string.h>

#define gf_add CRYPTO_NAMESPACE(gf_add)
#define gf_frac CRYPTO_NAMESPACE(gf_frac)
#define gf_inv CRYPTO_NAMESPACE(gf_inv)
#define gf_iszero CRYPTO_NAMESPACE(gf_iszero)
#define gf_mul CRYPTO_NAMESPACE(gf_mul)
#define GF_mul CRYPTO_NAMESPACE(GF_mul)

#include <stdint.h>
#include "params.h"

typedef uint16_t gf;

gf gf_iszero(gf);
gf gf_add(gf, gf);
gf gf_mul(gf, gf);
gf gf_frac(gf, gf);
gf gf_inv(gf);

void GF_mul(gf *, gf *, gf *);

//Customized Versions
static inline void gfmul_vx_custom_u16(uint16_t* dst_v, const uint16_t* src_v, uint16_t s, uint32_t len){
    csr_primpoly_rw(MC_GF_POLY);

    size_t vl,avl;
    avl=len;

    vuint16m2_t va,vb;

    while(avl>0){
        vl=vsetvl_e16m2(avl);
        va=vle16_v_u16m2(src_v,vl);
        vb=vgfmul_vx_u16m2(va,s);
        vse16_v_u16m2(dst_v,vb,vl);

        dst_v+=vl,src_v+=vl,avl-=vl;
    }
}

#define GF_mul_custom CRYPTO_NAMESPACE(GF_mul_custom)
void GF_mul_custom(gf *out, gf *in0, gf *in1);
gf gf_inv_custom(gf den);
#endif

