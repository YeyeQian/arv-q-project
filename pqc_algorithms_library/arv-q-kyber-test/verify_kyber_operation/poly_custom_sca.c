#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include "ntt.h"
#include "fips202.h"
#include "params.h"
#include "poly.h"

void poly_basemul_montgomery_custom_shuffling(poly* r, const poly* a, const poly* b, const uint8_t buffer_r[128]){

    uint16_t poly_idx[256];
    uint16_t zeta_idx[128];

    for(int i=0;i<128;i++) {
        zeta_idx[i]=i;
    }

    for(int i=127; i>=0; i--) {             //shuffling in pairs
        int selection=buffer_r[i]%(i+1);

        //shuffling zeta index
        uint16_t temp=zeta_idx[i];
        zeta_idx[i]=zeta_idx[selection];
        zeta_idx[selection]=temp;
    }

    // for test
    // printf("shuffling order:\n");
    // for(int i=0;i<128;i++)printf("%d\n",zeta_idx[i]);

    // get poly index
    for(int i=0;i<128;i++){
        poly_idx[2*i]=zeta_idx[i]*2;
        poly_idx[2*i+1]=zeta_idx[i]*2+1;
    }

    size_t start, end;
    start = read_cycle();
    // basemul with shuffling
    const int16_t *a_ptr = a->coeffs;
    const int16_t *b_ptr = b->coeffs;
    int16_t *r_ptr = r->coeffs;
    const int16_t *zeta_ptr = zetas_basemul_cg;
    const uint16_t* zeta_idx_ptr=zeta_idx;
    const uint16_t* poly_idx_ptr=poly_idx;

    vint16m1_t v_zetas;
    vint16m1_t v_in0, v_in1;
    vuint16m1_t v_poly_idx, v_zeta_idx;
    vint16m1_t v_out0;
    
    size_t avl = KYBER_N;
    size_t vl, zeta_vl;
    while(avl > 0) {
        vl = vsetvl_e16m1(avl);
        zeta_vl = vl >> 1;
        
        v_zeta_idx=vle16_v_u16m1(zeta_idx_ptr,zeta_vl);
        v_zeta_idx=vsll_vx_u16m1(v_zeta_idx,1,zeta_vl);
        v_zetas = vluxei16_v_i16m1(zeta_ptr, v_zeta_idx, zeta_vl);

        v_poly_idx=vle16_v_u16m1(poly_idx_ptr,vl);
        v_poly_idx=vsll_vx_u16m1(v_poly_idx,1,vl);
        v_in0 = vluxei16_v_i16m1(a_ptr, v_poly_idx, vl);
        v_in1 = vluxei16_v_i16m1(b_ptr, v_poly_idx, vl);
        v_out0 = vmod_basemul_vvm_i16m1(v_in0, v_in1, v_zetas);
        vsuxei16_v_i16m1(r_ptr, v_poly_idx, v_out0, vl);
        
        zeta_idx_ptr+=zeta_vl;
        poly_idx_ptr+=vl;
        avl -= vl;
  }
  end = read_cycle();
  printf("it takes %ld cycles\n", end - start);
}

void poly_basemul_montgomery_custom_outside_shuffling(poly* r, const poly* a, const poly* b, const uint16_t zeta_idx[128], const uint16_t poly_idx[256]){
    // basemul with shuffling
    const int16_t *a_ptr = a->coeffs;
    const int16_t *b_ptr = b->coeffs;
    int16_t *r_ptr = r->coeffs;
    const int16_t *zeta_ptr = zetas_basemul_cg;
    const uint16_t* zeta_idx_ptr=zeta_idx;
    const uint16_t* poly_idx_ptr=poly_idx;

    vint16m1_t v_zetas;
    vint16m1_t v_in0, v_in1;
    vuint16m1_t v_poly_idx, v_zeta_idx;
    vint16m1_t v_out0;
    
    size_t avl = KYBER_N;
    size_t vl, zeta_vl;
    while(avl > 0) {
        vl = vsetvl_e16m1(avl);
        zeta_vl = vl >> 1;
        
        v_zeta_idx=vle16_v_u16m1(zeta_idx_ptr,zeta_vl);
        v_zeta_idx=vsll_vx_u16m1(v_zeta_idx,1,zeta_vl);
        v_zetas = vluxei16_v_i16m1(zeta_ptr, v_zeta_idx, zeta_vl);
        v_poly_idx=vle16_v_u16m1(poly_idx_ptr,vl);
        v_poly_idx=vsll_vx_u16m1(v_poly_idx,1,vl);
        v_in0 = vluxei16_v_i16m1(a_ptr, v_poly_idx, vl);
        v_in1 = vluxei16_v_i16m1(b_ptr, v_poly_idx, vl);
        v_out0 = vmod_basemul_vvm_i16m1(v_in0, v_in1, v_zetas);
        vsuxei16_v_i16m1(r_ptr, v_poly_idx, v_out0, vl);
        
        zeta_idx_ptr+=zeta_vl;
        poly_idx_ptr+=vl;
        avl -= vl;
    }
}

/*************************************************
* Name:        poly_basemul_montgomery_custom_redundant
*
* Description: Multiplication of two polynomials in NTT domain
*                     Cooperate with CG NTT and INTT
*                     CG NTT output is reordered to standard order, to prepare for CG INTT
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul_montgomery_custom_redundant(poly* r, const poly* a, const poly* b, const uint8_t buffer_r[128])
{
    //basemul with redundant operation
    const int16_t *a_ptr = a->coeffs;
    const int16_t *b_ptr = b->coeffs;
    int16_t *r_ptr = r->coeffs;
    const int16_t *zeta_ptr = zetas_basemul_cg;

    vint16m1_t v_zetas;
    vint16m1_t v_in0, v_in1;
    vint16m1_t v_out0;
    
    size_t avl = KYBER_N;
    size_t vl, zeta_vl;

    while(avl > 0) {
        //randomize vl and redundant op num
        uint8_t nonce=buffer_r[(avl>>1)-1];
        uint8_t nonce_low_weigth=(nonce&0x01)+((nonce>>1)&0x01)+((nonce>>2)&0x01);
        uint8_t vl_pre=(nonce_low_weigth+1)*8;
        vl=vl_pre>=avl?avl:vl_pre;
        uint8_t nonce_high=nonce>>5;

        
        vl = vsetvl_e16m1(vl);
        zeta_vl = vl >> 1;

        v_zetas = vle16_v_i16m1(zeta_ptr, zeta_vl);
        v_in0 = vle16_v_i16m1(a_ptr, vl);
        v_in1 = vle16_v_i16m1(b_ptr, vl);

        keccak_squeeze();
        keccak_squeeze();
        v_out0 = vmod_basemul_vvm_i16m1(v_in0, v_in1, v_zetas);
        keccak_squeeze();
        keccak_squeeze();

        vse16_v_i16m1(r_ptr, v_out0, vl);

        zeta_ptr += zeta_vl;
        a_ptr += vl;
        b_ptr += vl;
        r_ptr += vl;
        avl -= vl;
        // round_cnt++;
    } 
    // printf("cost %d rounds\n",round_cnt);
}

void poly_basemul_montgomery_custom_redundant_max(poly* r, const poly* a, const poly* b, const uint8_t buffer_r[128])
{
    //basemul with redundant operation
    const int16_t *a_ptr = a->coeffs;
    const int16_t *b_ptr = b->coeffs;
    int16_t *r_ptr = r->coeffs;
    const int16_t *zeta_ptr = zetas_basemul_cg;

    vint16m1_t v_zetas;
    vint16m1_t v_in0, v_in1;
    vint16m1_t v_out0;
    
    size_t avl = KYBER_N;
    size_t vl, zeta_vl;
    while(avl > 0) {
        uint8_t nonce=buffer_r[(avl>>1)-1];
        uint8_t vl_pre=8;
        vl=vl_pre>=avl?avl:vl_pre;
        uint8_t nonce_high=7;

        
        vl = vsetvl_e16m1(vl);
        zeta_vl = vl >> 1;

        v_zetas = vle16_v_i16m1(zeta_ptr, zeta_vl);
        v_in0 = vle16_v_i16m1(a_ptr, vl);
        v_in1 = vle16_v_i16m1(b_ptr, vl);

        keccak_squeeze();
        keccak_squeeze();
        v_out0 = vmod_basemul_vvm_i16m1(v_in0, v_in1, v_zetas);
        keccak_squeeze();
        keccak_squeeze();

        vse16_v_i16m1(r_ptr, v_out0, vl);
        
        zeta_ptr += zeta_vl;
        a_ptr += vl;
        b_ptr += vl;
        r_ptr += vl;
        avl -= vl;
    } 
}