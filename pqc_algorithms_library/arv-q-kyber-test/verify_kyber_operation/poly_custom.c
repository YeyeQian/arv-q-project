#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include "fips202.h"
#include "params.h"
#include "symmetric.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "cbd.h"


/*************************************************
* Name:        poly_compress_custom
*
* Description: Compression and subsequent serialization of a polynomial
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (of length KYBER_POLYCOMPRESSEDBYTES)
*              - poly *a:    pointer to input polynomial
**************************************************/
void poly_compress_custom(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], poly *a)
{
  int16_t *coeffs_ptr = a->coeffs;
  size_t vl;
  size_t avl = KYBER_N;

#if (KYBER_POLYCOMPRESSEDBYTES == 128)
  vuint16m2_t vx;
  vuint8m1_t vt;

  for(; avl > 0; avl -= vl) {
    vl = vsetvl_e16m2(avl);
    vx = vle16_v_u16m2((uint16_t*)coeffs_ptr, vl);
    coeffs_ptr += vl;

    vx = vsll_vx_u16m2(vx, 4, vl);
    vx = vadd_vx_u16m2(vx, KYBER_Q/2, vl);

#if (USE_MULTIPLY_SHIFT_OPT == 1)
    /*
    **  use multiply and shift to replace divide
    ** ceil(a/q) = (a*m) >> k
    ** a_max = 3328 * 2^4 + KYBER_Q/2
    ** k = floor(log2(a_max*q)) = 28
    ** m = floor(2^k/q) = 80636
    ** however, floor(log2(80636)) = 17, thus vx * 80636 is 33bit
    ** as a result, instead of using 32bit width, 64bit should be used
    ** which needs width extensiona and narrow instructions
    ** so that use k=27 and m=40318
    */
    vuint32m4_t vz;
    vz = vwmulu_vx_u32m4(vx, 40318, vl);
    vx = vnsrl_wx_u16m2(vz, 27, vl);
#else
    /*
    **  directly use divide
    */
    vx = vdivu_vx_u16m2(vx, KYBER_Q, vl);
#endif

    vx = vand_vx_u16m2(vx, 15, vl);
    vl = vsetvl_e8m1(vl);
    vt = vnsrl_wx_u8m1(vx, 0, vl);

    // pack with valid bits number = 4
    vt = pack_vx_u8m1(vt, 4);
    size_t vl_packed = vl >> 1;
    vse8_v_u8m1(r, vt, vl_packed);
    r += vl_packed;
  }

#elif (KYBER_POLYCOMPRESSEDBYTES == 160)
  vuint16m2_t vx;
  vuint32m4_t vy;
  vuint8m1_t vt;

  for(; avl > 0; avl -= vl) {
    vl = vsetvl_e16m2(avl);
    vx = vle16_v_u16m2((uint16_t*)coeffs_ptr, vl);
    coeffs_ptr += vl;

    vy = vzext_vf2_u32m4(vx, vl);
    vy = vsll_vx_u32m4(vy, 5, vl);
    vy = vadd_vx_u32m4(vy, KYBER_Q/2, vl);

#if (USE_MULTIPLY_SHIFT_OPT == 1)
    /*
    **  use multiply and shift to replace divide
    ** ceil(a/q) = (a*m) >> k
    ** a_max = 3328 * 2^5 + KYBER_Q/2
    ** k = floor(log2(a_max*q)) = 29
    ** m = floor(2^k/q) = 161271
    */
   vuint64m8_t vz;
    vz = vwmulu_vx_u64m8(vy, 161271, vl);
    vy = vnsrl_wx_u32m4(vz, 29, vl);
#else
    /*
    **  directly use divide
    */
    vy = vdivu_vx_u32m4(vy, KYBER_Q, vl);
#endif

    vy = vand_vx_u32m4(vy, 31, vl);
    vx = vnsrl_wx_u16m2(vy, 0, vl);
    vt = vnsrl_wx_u8m1(vx, 0, vl);

    // pack with valid bits number = 5
    vt = pack_vx_u8m1(vt, 5);
    size_t vl_packed = (5 * vl) >> 3;
    vse8_v_u8m1(r, vt, vl_packed);
    r += vl_packed;
  }     
#else
#error "KYBER_POLYCOMPRESSEDBYTES needs to be in {128, 160}"
#endif
}



/*************************************************
* Name:        poly_compress_rvv
*
* Description: Compression and subsequent serialization of a polynomial
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (of length KYBER_POLYCOMPRESSEDBYTES)
*              - poly *a:    pointer to input polynomial
**************************************************/
void poly_compress_rvv(uint8_t r[KYBER_POLYCOMPRESSEDBYTES], poly *a)
{
  unsigned int i;
  int16_t *coeffs_ptr = a->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  uint8_t t[KYBER_N];
  uint8_t *t_ptr = t;

#if (KYBER_POLYCOMPRESSEDBYTES == 128)
  vuint16m2_t vx;
  vuint8m1_t vt;

  for(; avl > 0; avl -= vl) {
    vl = vsetvl_e16m2(avl);
    vx = vle16_v_u16m2((uint16_t*)coeffs_ptr, vl);
    coeffs_ptr += vl;

    vx = vsll_vx_u16m2(vx, 4, vl);
    vx = vadd_vx_u16m2(vx, KYBER_Q/2, vl);

#if (USE_MULTIPLY_SHIFT_OPT == 1)
    /*
    **  use multiply and shift to replace divide
    ** ceil(a/q) = (a*m) >> k
    ** a_max = 3328 * 2^4 + KYBER_Q/2
    ** k = floor(log2(a_max*q)) = 28
    ** m = floor(2^k/q) = 80636
    ** however, floor(log2(80636)) = 17, thus vx * 80636 is 33bit
    ** as a result, instead of using 32bit width, 64bit should be used
    ** which needs width extensiona and narrow instructions
    ** so that use k=27 and m=40318
    */
    vuint32m4_t vz;
    vz = vwmulu_vx_u32m4(vx, 40318, vl);
    vx = vnsrl_wx_u16m2(vz, 27, vl);
#else
    /*
    **  directly use divide
    */
    vx = vdivu_vx_u16m2(vx, KYBER_Q, vl);
#endif

    vx = vand_vx_u16m2(vx, 15, vl);
    vl = vsetvl_e8m1(vl);
    vt = vnsrl_wx_u8m1(vx, 0, vl);

    vse8_v_u8m1(t_ptr, vt, vl);
    t_ptr += vl;
  }

  t_ptr = t;
  for(i=0;i<KYBER_N/8;i++) {
    r[0] = t_ptr[0] | (t_ptr[1] << 4);
    r[1] = t_ptr[2] | (t_ptr[3] << 4);
    r[2] = t_ptr[4] | (t_ptr[5] << 4);
    r[3] = t_ptr[6] | (t_ptr[7] << 4);
    r += 4;
    t_ptr += 8;
  }

#elif (KYBER_POLYCOMPRESSEDBYTES == 160)
  vuint16m2_t vx;
  vuint32m4_t vy;
  vuint8m1_t vt;

  for(; avl > 0; avl -= vl) {
    vl = vsetvl_e16m2(avl);
    vx = vle16_v_u16m2((uint16_t*)coeffs_ptr, vl);
    coeffs_ptr += vl;

    vy = vzext_vf2_u32m4(vx, vl);
    vy = vsll_vx_u32m4(vy, 5, vl);
    vy = vadd_vx_u32m4(vy, KYBER_Q/2, vl);

#if (USE_MULTIPLY_SHIFT_OPT == 1)
    /*
    **  use multiply and shift to replace divide
    ** ceil(a/q) = (a*m) >> k
    ** a_max = 3328 * 2^5 + KYBER_Q/2
    ** k = floor(log2(a_max*q)) = 29
    ** m = floor(2^k/q) = 161271
    */
   vuint64m8_t vz;
    vz = vwmulu_vx_u64m8(vy, 161271, vl);
    vy = vnsrl_wx_u32m4(vz, 29, vl);
#else
    /*
    **  directly use divide
    */
    vy = vdivu_vx_u32m4(vy, KYBER_Q, vl);
#endif

    vy = vand_vx_u32m4(vy, 31, vl);
    vx = vnsrl_wx_u16m2(vy, 0, vl);
    vt = vnsrl_wx_u8m1(vx, 0, vl);

    vse8_v_u8m1(t_ptr, vt, vl);
    t_ptr += vl;
  }

  t_ptr = t;
  for(i=0;i<KYBER_N/8;i++) {
    r[0] = (t_ptr[0] >> 0) | (t_ptr[1] << 5);
    r[1] = (t_ptr[1] >> 3) | (t_ptr[2] << 2) | (t_ptr[3] << 7);
    r[2] = (t_ptr[3] >> 1) | (t_ptr[4] << 4);
    r[3] = (t_ptr[4] >> 4) | (t_ptr[5] << 1) | (t_ptr[6] << 6);
    r[4] = (t_ptr[6] >> 2) | (t_ptr[7] << 3);
    r += 5;
    t_ptr += 8;
  }       
#else
#error "KYBER_POLYCOMPRESSEDBYTES needs to be in {128, 160}"
#endif
}

/*************************************************
* Name:        poly_decompress_custom
*
* Description: De-serialization and subsequent decompression of a polynomial;
*              approximate inverse of poly_compress
*
* Arguments:   - poly *r:          pointer to output polynomial
*              - const uint8_t *a: pointer to input byte array
*                                  (of length KYBER_POLYCOMPRESSEDBYTES bytes)
**************************************************/
void poly_decompress_custom(poly *r, const uint8_t a[KYBER_POLYCOMPRESSEDBYTES])
{
  int16_t *coeffs_ptr = r->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
#if (KYBER_POLYCOMPRESSEDBYTES == 128)
  vuint8m1_t vt;
  vuint16m2_t vx;
  for(;avl>0;avl-=vl){
    vl=vsetvl_e8m1(avl);
    size_t vl_unpacked = (4 * vl) >> 3;
    vt=vle8_v_u8m1(a,vl_unpacked);
    a+=vl_unpacked;

    //should use asm volatile otherwise will be optimized out and tirger panic if
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e8, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );

    vt=unpack_vx_u8m1(vt,4);
    vx=vzext_vf2_u16m2(vt,vl);
    vx=vmul_vx_u16m2(vx,KYBER_Q,vl);
    vx=vadd_vx_u16m2(vx,8,vl);
    vx=vsrl_vx_u16m2(vx,4,vl);
    vse16_v_u16m2((uint16_t*)coeffs_ptr,vx,vl);
    coeffs_ptr+=vl;
  }
#elif (KYBER_POLYCOMPRESSEDBYTES == 160)
  vuint8m1_t vt;
  vuint16m2_t vx;
  vuint32m4_t vy;
  for(;avl>0;avl-=vl){
    vl=vsetvl_e8m1(avl);
    size_t vl_unpacked=(5*vl)>>3;
    vt=vle8_v_u8m1(a,vl_unpacked);
    a+=vl_unpacked;

    //should use asm volatile otherwise will be optimized out and tirger panic if
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e8, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );

    vt=unpack_vx_u8m1(vt,5);
    vx=vzext_vf2_u16m2(vt,vl);
    vy=vzext_vf2_u32m4(vx,vl);
    vy=vmul_vx_u32m4(vy,KYBER_Q,vl);
    vy=vadd_vx_u32m4(vy,16,vl);
    vx=vnsrl_wx_u16m2(vy,5,vl);
    vse16_v_u16m2(coeffs_ptr,vx,vl);
    coeffs_ptr+=vl;
  }
#else
#error "KYBER_POLYCOMPRESSEDBYTES needs to be in {128, 160}"
#endif
}

/*************************************************
* Name:        poly_tobytes_custom
*
* Description: Serialization of a polynomial
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (needs space for KYBER_POLYBYTES bytes)
*              - poly *a:    pointer to input polynomial
**************************************************/
void poly_tobytes_custom(uint8_t r[KYBER_POLYBYTES], poly *a)
{
  int16_t *coeffs_ptr = a->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  vuint16m1_t vt;
  vuint8m1_t vx;
  while(avl>0){
    vl=vsetvl_e16m1(avl);
    vt=vle16_v_u16m1((uint16_t*)coeffs_ptr,vl);
    avl-=vl;
    coeffs_ptr+=vl;
    vx=pack_vx_u16m1(vt,12);
    size_t vl_packed=(12*vl)>>3;
    vse8_v_u8m1(r,vx,vl_packed);
    r+=vl_packed;
  }
}

/*************************************************
* Name:        poly_frombytes_custom
*
* Description: De-serialization of a polynomial;
*              inverse of poly_tobytes
*
* Arguments:   - poly *r:          pointer to output polynomial
*              - const uint8_t *a: pointer to input byte array
*                                  (of KYBER_POLYBYTES bytes)
**************************************************/
void poly_frombytes_custom(poly *r, const uint8_t a[KYBER_POLYBYTES])
{
  int16_t *coeffs_ptr = r->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  vuint16m1_t vx;
  vuint8m1_t vt;
  while(avl > 0) {
    vl = vsetvl_e16m1(avl);
    size_t vl_packed = (12*vl) >> 3;
    vt = vle8_v_u8m1(a,vl_packed);
    a += vl_packed;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e16, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    vx = unpack_vx_u16m1(vt,12);
    vse16_v_u16m1((uint16_t*)coeffs_ptr, vx, vl);
    coeffs_ptr += vl;
    avl -= vl;
  }
}

/*************************************************
* Name:        poly_frommsg_custom
*
* Description: Convert 32-byte message to polynomial
*
* Arguments:   - poly *r:            pointer to output polynomial
*              - const uint8_t *msg: pointer to input message
**************************************************/
void poly_frommsg_custom(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
{

#if (KYBER_INDCPA_MSGBYTES != KYBER_N/8)
#error "KYBER_INDCPA_MSGBYTES must be equal to KYBER_N/8 bytes!"
#endif

  int16_t *coeffs_ptr = r->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  vint16m1_t vx;
  vuint8m1_t vt;
  int16_t cons = (KYBER_Q+1) >> 1;
  while(avl > 0) {
    vl = vsetvl_e16m1(avl);
    size_t vl_packed =vl >> 3;
    vt = vle8_v_u8m1(msg, vl_packed);
    msg += vl_packed;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e16, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    vx = unpack_vx_i16m1(vt, 1);
    vx = vrsub_vx_i16m1(vx, 0, vl);
    vx = vand_vx_i16m1(vx, cons, vl);
    vse16_v_i16m1(coeffs_ptr, vx, vl);
    coeffs_ptr += vl;
    avl -= vl;
  }
}

/*************************************************
* Name:        poly_tomsg_custom
*
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - poly *a:      pointer to input polynomial
*
* Requirements: -coefficient of a in {0..2q-1}
**************************************************/
void poly_tomsg_custom(uint8_t msg[KYBER_INDCPA_MSGBYTES], poly *a)
{
  int16_t *coeffs_ptr = a->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  vint16m1_t vx;
  vuint8m1_t vt;
  while(avl>0){
    vl=vsetvl_e16m1(avl);
    vx=vle16_v_i16m1(coeffs_ptr,vl);
    // vx=vmod_add_vx_i16m1(vx,KYBER_Q);
    vx=vsll_vx_i16m1(vx,1,vl);
    vx=vadd_vx_i16m1(vx,KYBER_Q>>1,vl);
#if (USE_MULTIPLY_SHIFT_OPT == 1)
    /*
    **  use multiply and shift to replace divide
    ** ceil(a/q) = (a*m) >> k
    ** a_max = 3328 * 2^1 + KYBER_Q/2
    ** k = floor(log2(a_max*q)) = 25
    ** m = floor(2^k/q) = 10080
    */
    vint32m2_t vz;
    vz = vwmul_vx_i32m2(vx, 10080, vl);
    vx = vnsra_wx_i16m1(vz, 25, vl);
#else
    /*
    **  directly use divide
    */
    vx = vdiv_vx_i16m1(vx, KYBER_Q, vl);
#endif
    vx=vand_vx_i16m1(vx,1,vl);
    coeffs_ptr+=vl;
    avl-=vl;
    size_t vl_packed=vl>>3;
    vt=pack_vx_i16m1(vx,1);
    vse8_v_u8m1(msg,vt,vl_packed);
    msg+=vl_packed;
  }
}

/*************************************************
* Name:        poly_getnoise_eta1_custom
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA1
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length KYBER_SYMBYTES bytes)
*              - uint8_t nonce:       one-byte input nonce
**************************************************/
void poly_getnoise_eta1_custom(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce)
{
  uint8_t buf[KYBER_ETA1*KYBER_N/4];
  kyber_shake256_prf_custom(buf, sizeof(buf), seed, nonce);
  cbd_eta1_custom(r, buf);
}

/*************************************************
* Name:        poly_getnoise_eta2_custom
*
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA2
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length KYBER_SYMBYTES bytes)
*              - uint8_t nonce:       one-byte input nonce
**************************************************/
void poly_getnoise_eta2_custom(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce)
{
  uint8_t buf[KYBER_ETA2*KYBER_N/4];
  kyber_shake256_prf_custom(buf, sizeof(buf), seed, nonce);
  cbd_eta2_custom(r, buf);
}

/*************************************************
* Name:        poly_add_custom
*
* Description: Add two polynomials
*
* Arguments: - poly *r:       pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add_custom(poly *r, const poly *a, const poly *b)
{
  size_t vl;
  size_t avl = KYBER_N;
  const int16_t* a_ptr=a->coeffs;
  const int16_t* b_ptr=b->coeffs;
  int16_t* r_ptr=r->coeffs;
  vint16m8_t va;
  vint16m8_t vb;
  vint16m8_t vr;
  while(avl > 0){
    vl = vsetvl_e16m8(avl);
    va = vle16_v_i16m8(a_ptr, vl);
    vb = vle16_v_i16m8(b_ptr, vl);
    vr = vmod_add_vv_i16m8(va, vb);
    vse16_v_i16m8(r_ptr, vr, vl);
    r_ptr += vl;
    a_ptr += vl;
    b_ptr += vl;
    avl -= vl;
  }
}

/*************************************************
* Name:        poly_sub_custom
*
* Description: Subtract two polynomials
*
* Arguments: - poly *r:       pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub_custom(poly *r, const poly *a, const poly *b)
{
  size_t vl;
  size_t avl = KYBER_N;
  const int16_t* a_ptr=a->coeffs;
  const int16_t* b_ptr=b->coeffs;
  int16_t* r_ptr=r->coeffs;
  vint16m8_t va;
  vint16m8_t vb;
  vint16m8_t vr;
  while(avl>0){
    vl=vsetvl_e16m8(avl);
    va=vle16_v_i16m8(a_ptr,vl);
    vb=vle16_v_i16m8(b_ptr,vl);
    vr=vmod_sub_vv_i16m8(va,vb);
    vse16_v_i16m8(r_ptr,vr,vl);
    r_ptr+=vl;
    a_ptr+=vl;
    b_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_tomont_custom
*
* Description: Inplace conversion of all coefficients of a polynomial
*              from normal domain to Montgomery domain
*
* Arguments:   - poly *r: pointer to input/output polynomial
**************************************************/
void poly_tomont_custom(poly *r)
{
  const int16_t f = (1ULL << 32) % KYBER_Q;
  int16_t *coeffs_ptr = r->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  vint16m8_t vx;
  while(avl>0){
    vl=vsetvl_e16m8(avl);
    vx=vle16_v_i16m8(coeffs_ptr,vl);
    vx=vmod_mul_vx_i16m8(vx,f);
    vse16_v_i16m8(coeffs_ptr,vx,vl);
    coeffs_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_basemul_montgomery_custom
*
* Description: Multiplication of two polynomials in NTT domain
*                     Cooperate with CG NTT and INTT
*                     CG NTT output is reordered to standard order, to prepare for CG INTT
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_basemul_montgomery_custom(poly* r, const poly* a, const poly* b)
{
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
    vl = vsetvl_e16m1(avl);
    zeta_vl = vl >> 1;

    v_zetas = vle16_v_i16m1(zeta_ptr, zeta_vl);
    v_in0 = vle16_v_i16m1(a_ptr, vl);
    v_in1 = vle16_v_i16m1(b_ptr, vl);
    v_out0 = vmod_basemul_vvm_i16m1(v_in0, v_in1, v_zetas);
    vse16_v_i16m1(r_ptr, v_out0, vl);
    
    zeta_ptr += zeta_vl;
    a_ptr += vl;
    b_ptr += vl;
    r_ptr += vl;
    avl -= vl;
  } 
}

/*************************************************
* Name:        poly_ntt_custom
*
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in normal order
*
* Arguments:   - uint16_t *r: pointer to in/output polynomial
**************************************************/
void poly_ntt_custom(poly *r)
{
  ntt_cg_custom_reorder(r->coeffs);
}

/*************************************************
* Name:        poly_invntt_tomont_custom
*
* Description: Computes inverse of negacyclic number-theoretic transform (NTT)
*              of a polynomial in place;
*              inputs assumed to be in normal order, output in normal order
*
* Arguments:   - uint16_t *a: pointer to in/output polynomial
**************************************************/
void poly_invntt_tomont_custom(poly *r)
{
  invntt_cg_custom_reorder(r->coeffs);
  poly_tomont_custom(r);
}

void mod_add_q(int16_t a[KYBER_N])
{
  size_t vl;
  vint16m8_t vx;
  size_t avl = KYBER_N;
  int16_t* coeffs_ptr = a;

  while(avl > 0) {
    vl = vsetvl_e16m8(avl);
    vx = vle16_v_i16m8(coeffs_ptr, vl);
    vx = vmod_add_vx_i16m8(vx, KYBER_Q);
    vse16_v_i16m8(coeffs_ptr, vx, vl);
    coeffs_ptr += vl;
    avl -= vl;
  }
}

/*************************************************
* Name:        poly_mod_add_q
*
* Description: modular add q to each coefficients in a polynomial
*              in this way turn each coefficients to non-negetive number
*
* Arguments:   - int16_t p[N]: input/output coefficient array
**************************************************/
void poly_mod_add_q(poly *a)
{
  mod_add_q(a->coeffs);
}

/*************************************************
* Name:        poly_frommsg_add_custom
*
* Description: Convert 32-byte message to polynomial,
*              and modular add this poly with poly r 
*
* Arguments:   - poly *r:            pointer to output polynomial
*              - const uint8_t *msg: pointer to input message
**************************************************/
void poly_frommsg_add_custom(poly *r, const uint8_t msg[KYBER_INDCPA_MSGBYTES])
{

#if (KYBER_INDCPA_MSGBYTES != KYBER_N/8)
#error "KYBER_INDCPA_MSGBYTES must be equal to KYBER_N/8 bytes!"
#endif

  int16_t *coeffs_ptr = r->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  vint16m1_t vr;
  vint16m1_t vx;
  vuint8m1_t vt;
  int16_t cons = (KYBER_Q+1) >> 1;
  while(avl > 0) {
    vl = vsetvl_e16m1(avl);
    size_t vl_packed =vl >> 3;
    vt = vle8_v_u8m1(msg, vl_packed);
    msg += vl_packed;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e16, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    vx = unpack_vx_i16m1(vt, 1);
    vx = vrsub_vx_i16m1(vx, 0, vl);
    vx = vand_vx_i16m1(vx, cons, vl);
    vr = vle16_v_i16m1(coeffs_ptr, vl);
    vr = vmod_add_vv_i16m1(vr, vx);
    vse16_v_i16m1(coeffs_ptr, vr, vl);
    coeffs_ptr += vl;
    avl -= vl;
  }
}

/*************************************************
* Name:        poly_sub_tomsg_custom
*
* Description: 1. poly r = poly a - poly b
*              2. Convert polynomial r to 32-byte message
*
* Arguments:   - uint8_t *msg: pointer to output message
*              - poly *a:      pointer to input polynomial
*              - poly *b:      pointer to input polynomial
*
* Requirements: -coefficient of a in {0..2q-1}
**************************************************/
void poly_sub_tomsg_custom(uint8_t msg[KYBER_INDCPA_MSGBYTES], poly *a, poly *b)
{
  int16_t *coeffs_ptr_a = a->coeffs;
  int16_t *coeffs_ptr_b = b->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  vint16m1_t vx;
  vint16m1_t vy;
  vuint8m1_t vt;
  while(avl>0){
    vl=vsetvl_e16m1(avl);
    vx=vle16_v_i16m1(coeffs_ptr_a,vl);
    vy=vle16_v_i16m1(coeffs_ptr_b,vl);
    vx=vmod_sub_vv_i16m1(vx, vy);
    // vx=vmod_add_vx_i16m1(vx,KYBER_Q);
    vx=vsll_vx_i16m1(vx,1,vl);
    vx=vadd_vx_i16m1(vx,KYBER_Q>>1,vl);
#if (USE_MULTIPLY_SHIFT_OPT == 1)
    /*
    **  use multiply and shift to replace divide
    ** ceil(a/q) = (a*m) >> k
    ** a_max = 3328 * 2^1 + KYBER_Q/2
    ** k = floor(log2(a_max*q)) = 25
    ** m = floor(2^k/q) = 10080
    */
    vint32m2_t vz;
    vz = vwmul_vx_i32m2(vx, 10080, vl);
    vx = vnsra_wx_i16m1(vz, 25, vl);
#else
    /*
    **  directly use divide
    */
    vx = vdiv_vx_i16m1(vx, KYBER_Q, vl);
#endif
    vx=vand_vx_i16m1(vx,1,vl);
    coeffs_ptr_a+=vl;
    coeffs_ptr_b+=vl;
    avl-=vl;
    size_t vl_packed=vl>>3;
    vt=pack_vx_i16m1(vx,1);
    vse8_v_u8m1(msg,vt,vl_packed);
    msg+=vl_packed;
  }
}