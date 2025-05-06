#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "polyvec.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>

/*************************************************
* Name:        polyvec_compress_custom
*
* Description: Compress and serialize vector of polynomials
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (needs space for KYBER_POLYVECCOMPRESSEDBYTES)
*              - polyvec *a: pointer to input vector of polynomials in {0,q}
**************************************************/
void polyvec_compress_custom(uint8_t r[KYBER_POLYVECCOMPRESSEDBYTES], polyvec *a)
{

  //polyvec_csubq(a);
    size_t vl;
#if (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 352))
  for(int i=0;i<KYBER_K;i++){
    int16_t* coeff_ptr=a->vec[i].coeffs;
    size_t avl = KYBER_N;
    vuint16m1_t vx;
    vuint32m2_t vy;
    vuint8m1_t vt;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        vx=vle16_v_u16m1((uint16_t*)coeff_ptr,vl);
        vy=vzext_vf2_u32m2(vx,vl);
        vy=vsll_vx_u32m2(vy,11,vl);
        vy=vadd_vx_u32m2(vy,KYBER_Q>>1,vl);
    #if (USE_MULTIPLY_SHIFT_OPT == 1)
        /*
        **  use multiply and shift to replace divide
        ** ceil(a/q) = (a*m) >> k
        ** a_max = 3328 * 2^11 + KYBER_Q/2
        ** k = floor(log2(a_max*q)) = 35
        ** m = floor(2^k/q) = 10321340
        */
        vuint64m4_t vz;
        vz = vwmulu_vx_u64m4(vy, 10321340, vl);
        vy = vnsrl_wx_u32m2(vz, 35, vl);
    #else
        /*
        **  directly use divide
        */
        vy = vdivu_vx_u32m2(vy, KYBER_Q, vl);
    #endif
        vy=vand_vx_u32m2(vy,0x7ff,vl);
        vx=vnsrl_wx_u16m1(vy,0,vl);
        size_t vl_packed=(vl*11)>>3;
        vt=pack_vx_u16m1(vx,11);
        avl-=vl;
        coeff_ptr+=vl;
        vse8_v_u8m1(r,vt,vl_packed);
        r+=vl_packed;
    }
  }
#elif (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 320))
  for(int i=0;i<KYBER_K;i++){
    int16_t* coeff_ptr = a->vec[i].coeffs;
    size_t avl = KYBER_N;
    vuint16m1_t vx;
    vuint32m2_t vy;
    vuint8m1_t vt;
    while(avl > 0) {
        vl = vsetvl_e16m1(avl);
        vx = vle16_v_u16m1((uint16_t*)coeff_ptr,vl);
        vy = vzext_vf2_u32m2(vx,vl);
        vy = vsll_vx_u32m2(vy,10,vl);
        vy = vadd_vx_u32m2(vy, KYBER_Q/2, vl);
    #if (USE_MULTIPLY_SHIFT_OPT == 1)
        /*
        **  use multiply and shift to replace divide
        ** ceil(a/q) = (a*m) >> k
        ** a_max = 3328 * 2^10 + KYBER_Q/2
        ** k = floor(log2(a_max*q)) = 34
        ** m = floor(2^k/q) = 5160670
        */
        vuint64m4_t vz;
        vz = vwmulu_vx_u64m4(vy, 5160670, vl);
        vy = vnsrl_wx_u32m2(vz, 34, vl);
    #else
        /*
        **  directly use divide
        */
        vy = vdivu_vx_u32m2(vy, KYBER_Q, vl);
    #endif
        vy = vand_vx_u32m2(vy, 0x3ff, vl);
        vx = vnsrl_wx_u16m1(vy, 0, vl);
        size_t vl_packed = (vl*10) >> 3;
        vt = pack_vx_u16m1(vx, 10);
        avl -= vl;
        coeff_ptr += vl;
        vse8_v_u8m1(r, vt, vl_packed);
        r += vl_packed;
    }
  }
#else
#error "KYBER_POLYVECCOMPRESSEDBYTES needs to be in {320*KYBER_K, 352*KYBER_K}"
#endif
}

/*************************************************
* Name:        polyvec_decompress_custom
*
* Description: De-serialize and decompress vector of polynomials;
*              approximate inverse of polyvec_compress
*
* Arguments:   - polyvec *r:       pointer to output vector of polynomials
*              - const uint8_t *a: pointer to input byte array
*                                  (of length KYBER_POLYVECCOMPRESSEDBYTES)
**************************************************/
void polyvec_decompress_custom(polyvec *r,
                        const uint8_t a[KYBER_POLYVECCOMPRESSEDBYTES])
{
size_t vl;
#if (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 352))
  for(int i=0;i<KYBER_K;i++){
    int16_t* coeff_ptr=r->vec[i].coeffs;
    size_t avl = KYBER_N;
    vuint16m1_t vx;
    vuint32m2_t vy;
    vuint8m1_t vt;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        size_t vl_packed=(vl*11)>>3;
        vt=vle8_v_u8m1(a,vl_packed);
        __asm__ __volatile__ ( "vsetvli %[vl], %[n], e16, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
        vx=unpack_vx_u16m1(vt,11);
        vy=vzext_vf2_u32m2(vx,vl);
        vy=vand_vx_u32m2(vy,0x7ff,vl);
        vy=vmul_vx_u32m2(vy,KYBER_Q,vl);
        vy=vadd_vx_u32m2(vy,1024,vl);
        vx=vnsrl_wx_u16m1(vy,11,vl);
        vse16_v_u16m1((uint16_t*)coeff_ptr,vx,vl);
        coeff_ptr+=vl;
        avl-=vl;
        a+=vl_packed;
    }
  }
#elif (KYBER_POLYVECCOMPRESSEDBYTES == (KYBER_K * 320))
  for(int i=0;i<KYBER_K;i++){
    int16_t* coeff_ptr=r->vec[i].coeffs;
    size_t avl = KYBER_N;
    vuint16m1_t vx;
    vuint32m2_t vy;
    vuint8m1_t vt;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        size_t vl_packed=(vl*10)>>3;
        vt=vle8_v_u8m1(a,vl_packed);
        __asm__ __volatile__ ( "vsetvli %[vl], %[n], e16, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
        vx=unpack_vx_u16m1(vt,10);
        vy=vzext_vf2_u32m2(vx,vl);
        vy=vand_vx_u32m2(vy,0x3ff,vl);
        vy=vmul_vx_u32m2(vy,KYBER_Q,vl);
        vy=vadd_vx_u32m2(vy,512,vl);
        vx=vnsrl_wx_u16m1(vy,10,vl);
        vse16_v_u16m1((uint16_t*)coeff_ptr,vx,vl);
        coeff_ptr+=vl;
        avl-=vl;
        a+=vl_packed;
    }
  }
#else
#error "KYBER_POLYVECCOMPRESSEDBYTES needs to be in {320*KYBER_K, 352*KYBER_K}"
#endif
}

/*************************************************
* Name:        polyvec_tobytes_custom
*
* Description: Serialize vector of polynomials
*
* Arguments:   - uint8_t *r: pointer to output byte array
*                            (needs space for KYBER_POLYVECBYTES)
*              - polyvec *a: pointer to input vector of polynomials
**************************************************/
void polyvec_tobytes_custom(uint8_t r[KYBER_POLYVECBYTES], polyvec *a)
{
  unsigned int i;
  for(i=0;i<KYBER_K;i++)
    poly_tobytes_custom(r+i*KYBER_POLYBYTES, &a->vec[i]);
}

/*************************************************
* Name:        polyvec_frombytes_custom
*
* Description: De-serialize vector of polynomials;
*              inverse of polyvec_tobytes
*
* Arguments:   - uint8_t *r:       pointer to output byte array
*              - const polyvec *a: pointer to input vector of polynomials
*                                  (of length KYBER_POLYVECBYTES)
**************************************************/
void polyvec_frombytes_custom(polyvec *r, const uint8_t a[KYBER_POLYVECBYTES])
{
  unsigned int i;
  for(i=0;i<KYBER_K;i++)
    poly_frombytes_custom(&r->vec[i], a+i*KYBER_POLYBYTES);
}

/*************************************************
* Name:        polyvec_ntt_custom
*
* Description: Apply forward NTT to all elements of a vector of polynomials
*
* Arguments:   - polyvec *r: pointer to in/output vector of polynomials
**************************************************/
void polyvec_ntt_custom(polyvec *r)
{
  unsigned int i;
  for(i=0;i<KYBER_K;i++)
    poly_ntt_custom(&r->vec[i]);
}

/*************************************************
* Name:        polyvec_invntt_tomont_custom
*
* Description: Apply inverse NTT to all elements of a vector of polynomials
*              and multiply by Montgomery factor 2^16
*
* Arguments:   - polyvec *r: pointer to in/output vector of polynomials
**************************************************/
void polyvec_invntt_tomont_custom(polyvec *r)
{
  unsigned int i;
  for(i=0;i<KYBER_K;i++) {
    poly_invntt_tomont_custom(&r->vec[i]);
  }
}

/*************************************************
* Name:        polyvec_pointwise_acc_montgomery_custom
*
* Description: Pointwise multiply elements of a and b, accumulate into r,
*              and multiply by 2^-16.
*
* Arguments: - poly *r:          pointer to output polynomial
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_pointwise_acc_montgomery_custom(poly *r,
                                      const polyvec *a,
                                      const polyvec *b)
{
  unsigned int i;
  poly t;

  poly_basemul_montgomery_custom(r, &a->vec[0], &b->vec[0]);
  for(i=1;i<KYBER_K;i++) {
    poly_basemul_montgomery_custom(&t, &a->vec[i], &b->vec[i]);
    poly_add_custom(r, r, &t);
  }

}

/*************************************************
* Name:        polyvec_add_custom
*
* Description: Add vectors of polynomials
*
* Arguments: - polyvec *r:       pointer to output vector of polynomials
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_add_custom(polyvec *r, const polyvec *a, const polyvec *b)
{
  unsigned int i;
  for(i=0;i<KYBER_K;i++)
    poly_add_custom(&r->vec[i], &a->vec[i], &b->vec[i]);
}

/*************************************************
* Name:        polyvec_mod_add_q
*
* Description: modular add q to each coefficients in every polynomial in polyvec
*              in this way turn each coefficients to non-negetive number
*
* Arguments:   - polyvec *r: pointer to in/output vector of polynomials
**************************************************/
void polyvec_mod_add_q(polyvec *r)
{
  unsigned int i;
  for(i=0; i < KYBER_K; i++)
    poly_mod_add_q(&r->vec[i]);
}