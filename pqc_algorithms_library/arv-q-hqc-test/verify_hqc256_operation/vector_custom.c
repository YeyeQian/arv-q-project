#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include "parameters.h"
#include "vector.h"
#include "fips202.h"
#include <stdint.h>
#include <string.h>
#include <stdio.h>

static inline uint32_t compare_u32(const uint32_t v1, const uint32_t v2) {
    return 1 ^ (((v1 - v2)|(v2 - v1)) >> 31);
}

void vect_set_random_fixed_weight_noctx_custom(uint32_t rand_u32[PARAM_OMEGA_R], uint64_t* v, uint16_t weight){
    uint32_t support[PARAM_OMEGA_R] = { 0 };
    uint32_t index_tab[PARAM_OMEGA_R] = { 0 };
    uint64_t bit_tab[PARAM_OMEGA_R] = { 0 };

    size_t avl=weight;
    size_t vl;

    vuint32m4_t va,vb,vc;
    uint32_t* addr0=rand_u32;
    uint32_t* addr1=support;
    while(avl>0){
        size_t offset=weight-avl;
        vl=vsetvl_e32m4(avl);
        va=vid_v_u32m4(vl);
        va=vadd_vx_u32m4(va,offset,vl);//denotes i
        vb=vle32_v_u32m4(addr0,vl);//denotes rand_u32[i]
        vc=vrsub_vx_u32m4(va,PARAM_N,vl);//denotes PARAM_N-i
        vb=vremu_vv_u32m4(vb,vc,vl);
        va=vadd_vv_u32m4(va,vb,vl);
        vse32_v_u32m4(addr1,va,vl);

        addr0+=vl,addr1+=vl,avl-=vl;
    }

    for (int32_t i = (weight - 1); i-- > 0;) {
        uint32_t found = 0;

        for (size_t j = i + 1; j < weight; ++j) {
            found |= compare_u32(support[j], support[i]);
        }

        uint32_t mask = -found;
        support[i] = (mask & i) ^ (~mask & support[i]);
    }

    for(int i=0;i<weight;i++){
        uint32_t position=support[i];
        int index = position>>6;
        int pos = position&0x3f;
        v[index]|=1ULL<<(pos);
    }
}

void vect_set_random_fixed_weight_noctx_support_custom(uint32_t rand_u32[PARAM_OMEGA_R],uint32_t* support, uint64_t* v, uint16_t weight){
    uint32_t index_tab[PARAM_OMEGA_R] = { 0 };
    uint64_t bit_tab[PARAM_OMEGA_R] = { 0 };

    size_t avl=weight;
    size_t vl;

    vuint32m4_t va,vb,vc;
    uint32_t* addr0=rand_u32;
    uint32_t* addr1=support;
    while(avl>0){
        size_t offset=weight-avl;
        vl=vsetvl_e32m4(avl);
        va=vid_v_u32m4(vl);
        va=vadd_vx_u32m4(va,offset,vl);//denotes i
        vb=vle32_v_u32m4(addr0,vl);//denotes rand_u32[i]
        vc=vrsub_vx_u32m4(va,PARAM_N,vl);//denotes PARAM_N-i
        vb=vremu_vv_u32m4(vb,vc,vl);
        va=vadd_vv_u32m4(va,vb,vl);
        vse32_v_u32m4(addr1,va,vl);

        addr0+=vl,addr1+=vl,avl-=vl;
    }

    for (int32_t i = (weight - 1); i-- > 0;) {
        uint32_t found = 0;

        for (size_t j = i + 1; j < weight; ++j) {
            found |= compare_u32(support[j], support[i]);
        }

        uint32_t mask = -found;
        support[i] = (mask & i) ^ (~mask & support[i]);
    }

    for(int i=0;i<weight;i++){
        uint32_t position=support[i];
        int index = position>>6;
        int pos = position&0x3f;
        v[index]|=1ULL<<(pos);
    }
}

void vect_set_random_fixed_weight_xy_custom(uint64_t* x, uint64_t* y, uint32_t* support, const unsigned char* sk_seed) {
    uint8_t domain = SEEDEXPANDER_DOMAIN;
    uint8_t sk_seed_pad[SEED_BYTES + 1];
    memcpy(sk_seed_pad,sk_seed,SEED_BYTES);
    sk_seed_pad[SEED_BYTES] = domain;
    uint8_t buffer[8 * PARAM_OMEGA] = { 0 };
    shake256_custom(buffer,sizeof(buffer),sk_seed_pad,sizeof(sk_seed_pad));

    uint32_t* rand_u32_x = (uint32_t*)buffer;
    uint32_t* rand_u32_y = (uint32_t*)buffer + PARAM_OMEGA;

    vect_set_random_fixed_weight_noctx_custom(rand_u32_x,x,PARAM_OMEGA);
    vect_set_random_fixed_weight_noctx_support_custom(rand_u32_y,support,y, PARAM_OMEGA);
}

void vect_set_random_fixed_weight_r1r2e_custom(uint64_t* r1, uint64_t* r2, uint32_t* support,uint64_t* e, const unsigned char* theta) {
    uint8_t domain = SEEDEXPANDER_DOMAIN;
    uint8_t theta_pad[SEED_BYTES + 1];
    memcpy(theta_pad, theta, SEED_BYTES);
    theta_pad[SEED_BYTES] = domain;

    uint8_t buffer[8*PARAM_OMEGA_R+4*PARAM_OMEGA_E] = { 0 };
    shake256_custom(buffer, sizeof(buffer), theta_pad, sizeof(theta_pad));

    uint32_t* rand_u32_r1 = (uint32_t*)buffer;
    uint32_t* rand_u32_r2 = rand_u32_r1 + PARAM_OMEGA_R;
    uint32_t* rand_u32_e = rand_u32_r2 + PARAM_OMEGA_R;

    vect_set_random_fixed_weight_noctx_custom(rand_u32_r1, r1, PARAM_OMEGA_R);
    vect_set_random_fixed_weight_noctx_support_custom(rand_u32_r2,support,r2, PARAM_OMEGA_R);
    vect_set_random_fixed_weight_noctx_custom(rand_u32_e, e, PARAM_OMEGA_E);
}

void vect_set_random_h_custom(unsigned char* pk_seed, uint64_t* v) {
    uint8_t domain = SEEDEXPANDER_DOMAIN;
    uint8_t pk_seed_pad[SEED_BYTES + 1];
    memcpy(pk_seed_pad, pk_seed, SEED_BYTES);
    pk_seed_pad[SEED_BYTES] = domain;

    shake256_custom((uint8_t*)v, VEC_N_SIZE_BYTES, pk_seed_pad, sizeof(pk_seed_pad));

    v[VEC_N_SIZE_64 - 1] &= BITMASK(PARAM_N, 64);
}

void vect_add_custom(uint64_t *o, const uint64_t *v1, const uint64_t *v2, uint32_t size){
    size_t avl,vl;
    avl=size;
    vuint64m8_t va,vb;
    while(avl>0){
        vl=vsetvl_e64m8(avl);
        va=vle64_v_u64m8(v1,vl);
        vb=vle64_v_u64m8(v2,vl);
        va=vxor_vv_u64m8(va,vb,vl);
        vse64_v_u64m8(o,va,vl);

        v1+=vl,v2+=vl,o+=vl,avl-=vl;
    }
}

uint8_t vect_compare_custom(const uint8_t *v1, const uint8_t *v2, uint32_t size) {
    uint64_t r = 0;

    size_t avl,vl;
    avl=size;
    vuint8m8_t va,vb;
    vuint8m1_t vc;
    vl=vsetvl_e8m1(avl);
    vc=vmv_v_x_u8m1(0,vl);
    while(avl>0){
        vl=vsetvl_e8m8(avl);
        va=vle8_v_u8m8(v1,vl);
        vb=vle8_v_u8m8(v2,vl);
        va=vxor_vv_u8m8(va,vb,vl);
        vc=vredor_vs_u8m8_u8m1(vc,va,vc,vl);

        v1+=vl,v2+=vl,avl-=vl;
    }

    r=vmv_x_s_u8m1_u8(vc);
    r = (~r + 1) >> 63;
    return (uint8_t) r;
}