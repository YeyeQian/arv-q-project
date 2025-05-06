#include <stdio.h>
#include <stdint.h>
#include "param.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "rounding.h"
#include "symmetric.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

void poly_print(const poly* a){//for debug usage
  int32_t* a_ptr=a->coeffs;
  for(int i=0;i<N;i++)printf("%d\n",a_ptr[i]);
}

void poly_copy(poly* dst, const poly* src){
  int32_t* dst_ptr=dst->coeffs;
  int32_t* src_ptr=src->coeffs;
  for(int i=0;i<N;i++)dst_ptr[i]=src_ptr[i];
}

//change coefficients order of polynomial.
//Input normal order then output bit-reverse order; Input bit-reverse order then output normal order.
void poly_chorder(poly* a){
  chorder(a->coeffs);
}

/*************************************************
* Name:        poly_reduce_rvv
*
* Description: Inplace reduction of all coefficients of polynomial to
*              representative in [-6283009,6283007].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_reduce_rvv(poly *a) {
  size_t vl;
  size_t avl=N;
  int32_t* coeff_ptr=a->coeffs;
  vint32m8_t vx;
  vint32m8_t vy;
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    vx=vle32_v_i32m8(coeff_ptr,vl);
    vy=vadd_vx_i32m8(vx,1<<22,vl);
    vy=vsra_vx_i32m8(vy,23,vl);
    vy=vmul_vx_i32m8(vy,Q,vl);
    vx=vsub_vv_i32m8(vx,vy,vl);
    vse32_v_i32m8(coeff_ptr,vx,vl);
    coeff_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_caddq_rvv
*
* Description: For all coefficients of in/out polynomial add Q if
*              coefficient is negative.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_caddq_rvv(poly *a) {
  size_t vl;
  size_t avl=N;
  int32_t* coeff_ptr=a->coeffs;
  vint32m8_t vx;
  vint32m8_t vy;
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    vx=vle32_v_i32m8(coeff_ptr,vl);
    vy=vsra_vx_i32m8(vx,31,vl);
    vy=vand_vx_i32m8(vy,Q,vl);
    vx=vadd_vv_i32m8(vy,vx,vl);
    vse32_v_i32m8(coeff_ptr,vx,vl);
    coeff_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_caddq_custom
*
* Description: Coefficients of input polynomial range in (-Q,Q).
*              Conduct modular addition of Q for all coefficients.
*              Coefficients of output polynomial range in (0,Q).
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_caddq_custom(poly *a) {
  size_t vl;
  size_t avl=N;
  int32_t* coeff_ptr=a->coeffs;
  vint32m8_t vx;
  csr_modulusq_rw(Q);
  csr_qinv_rw(QINV_HW);
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    vx=vle32_v_i32m8(coeff_ptr,vl);
    vx=vmod_add_vx_i32m8(vx,Q);
    vse32_v_i32m8(coeff_ptr,vx,vl);
    coeff_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_freeze_rvv
*
* Description: Inplace reduction of all coefficients of polynomial to
*              standard representatives.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_freeze_rvv(poly *a) {
  poly_reduce_rvv(a);
  poly_caddq_rvv(a);
}

/*************************************************
* Name:        poly_add_custom
*
* Description: Add polynomials. Modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first summand
*              - const poly *b: pointer to second summand
*
* Requirments: poly's coefficients must all be positive 
**************************************************/
void poly_add_custom(poly *c, const poly *a, const poly *b)  {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  int32_t* b_ptr=b->coeffs;
  int32_t* c_ptr=c->coeffs;
  vint32m8_t va;
  vint32m8_t vb;
  vint32m8_t vc;
  csr_modulusq_rw(Q);
  csr_qinv_rw(QINV_HW);
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    va=vle32_v_i32m8(a_ptr,vl);
    vb=vle32_v_i32m8(b_ptr,vl);
    vc=vmod_add_vv_i32m8(va,vb);
    vse32_v_i32m8(c_ptr,vc,vl);
    a_ptr+=vl;
    b_ptr+=vl;
    c_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_add_rvv
*
* Description: Add polynomials. No Modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first summand
*              - const poly *b: pointer to second summand
*
* Requirments: poly's coefficients must all be positive 
**************************************************/
void poly_add_rvv(poly *c, const poly *a, const poly *b)  {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  int32_t* b_ptr=b->coeffs;
  int32_t* c_ptr=c->coeffs;
  vint32m8_t va;
  vint32m8_t vb;
  vint32m8_t vc;
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    va=vle32_v_i32m8(a_ptr,vl);
    vb=vle32_v_i32m8(b_ptr,vl);
    vc=vadd_vv_i32m8(va,vb,vl);
    vse32_v_i32m8(c_ptr,vc,vl);
    a_ptr+=vl;
    b_ptr+=vl;
    c_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_sub_custom
*
* Description: Subtract polynomials. Modular reduction is
*              performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial to be
*                               subtraced from first input polynomial
*
* Requirments: poly's coefficients must all be positive
**************************************************/
void poly_sub_custom(poly *c, const poly *a, const poly *b)  {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  int32_t* b_ptr=b->coeffs;
  int32_t* c_ptr=c->coeffs;
  vint32m8_t va;
  vint32m8_t vb;
  vint32m8_t vc;
  csr_modulusq_rw(Q);
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    va=vle32_v_i32m8(a_ptr,vl);
    vb=vle32_v_i32m8(b_ptr,vl);
    vc=vmod_sub_vv_i32m8(va,vb);
    vse32_v_i32m8(c_ptr,vc,vl);
    a_ptr+=vl;
    b_ptr+=vl;
    c_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_tomont_custom
*
* Description: Input polynomial in montgomery domain, output polynomial in normal domain.
*
* Arguments:   - poly *a: pointer to input/output polynomial
*
* Requirments: poly's coefficients must all be positive
**************************************************/
void poly_mon2nor_custom(poly *a)  {
  size_t vl;
  size_t avl=N;
  int32_t cons=2365951;
  int32_t* a_ptr=a->coeffs;
  vint32m8_t va;
  csr_modulusq_rw(Q);
  csr_qinv_rw(QINV_HW);
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    va=vle32_v_i32m8(a_ptr,vl);
    va=vmod_mul_vx_i32m8(va,cons);
    vse32_v_i32m8(a_ptr,va,vl);
    a_ptr+=vl;
    avl-=vl;
  }
}


/*************************************************
* Name:        poly_shiftl_rvv
*
* Description: Multiply polynomial by 2^D without modular reduction. Assumes
*              input coefficients to be less than 2^{31-D} in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_shiftl_rvv(poly *a) {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  vint32m8_t va;
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    va=vle32_v_i32m8(a_ptr,vl);
    va=vsll_vx_i32m8(va,D,vl);
    vse32_v_i32m8(a_ptr,va,vl);
    a_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_ntt_custom
*
* Description: Outof-place NTT(Constant Geometry NTT). Input and output are both in normal order.
*              Coefficients of both input and output polynomial range in (0,Q).
*              Coefficients of both input and output polynomial are in normal domain (not montgomery domain).
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_ntt_custom(poly *a) {
  ntt_cg_custom_reorder(a);
}

/*************************************************
* Name:        poly_invntt_custom
*
* Description: Outof-place inverse NTT. Input and output are both in normal order.
*              Coefficients of both input and output polynomial range in (0,Q).
*              Input and output polynomial are both in montgomery domain.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invntt_custom(poly *a) {
  invntt_cg_custom_reorder(a);
}

/*************************************************
* Name:        poly_pointwise_montgomery_custom
*
* Description: Pointwise multiplication of polynomials in NTT domain
*              representation. Input polynomial in normal domain.
*              Output polynomial in montgomery domain.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_pointwise_montgomery_custom(poly *c, const poly *a, const poly *b) {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  int32_t* b_ptr=b->coeffs;
  int32_t* c_ptr=c->coeffs;
  vint32m8_t va,vb,vc;
  csr_modulusq_rw(Q);
  csr_qinv_rw(QINV_HW);
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    va=vle32_v_i32m8(a_ptr,vl);
    vb=vle32_v_i32m8(b_ptr,vl);
    vc=vmod_mul_vv_i32m8(va,vb);
    vse32_v_i32m8(c_ptr,vc,vl);
    a_ptr+=vl,b_ptr+=vl,c_ptr+=vl,avl-=vl;
  }
}

/*************************************************
* Name:        poly_power2round_rvv
*
* Description: For all coefficients c of the input polynomial,
*              compute c0, c1 such that c mod Q = c1*2^D + c0
*              with -2^{D-1} < c0 <= 2^{D-1}. Assumes coefficients to be
*              standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_power2round_rvv(poly *a1, poly *a0, const poly *a) {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  int32_t* a0_ptr=a0->coeffs;
  int32_t* a1_ptr=a1->coeffs;
  vint32m8_t vx,vy;
  while(avl>0){
    vl=vsetvl_e32m8(avl);
    vx=vle32_v_i32m8(a_ptr,vl);
    vy=vadd_vx_i32m8(vx,(1 << (D-1)) - 1,vl);
    vy=vsra_vx_i32m8(vy,D,vl);
    vse32_v_i32m8(a1_ptr,vy,vl);
    vy=vsll_vx_i32m8(vy,D,vl);
    vx=vsub_vv_i32m8(vx,vy,vl);
    vse32_v_i32m8(a0_ptr,vx,vl);
    a_ptr+=vl;
    a0_ptr+=vl;
    a1_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_decompose_rvv
*
* Description: For all coefficients c of the input polynomial,
*              compute high and low bits c0, c1 such c mod Q = c1*ALPHA + c0
*              with -ALPHA/2 < c0 <= ALPHA/2 except c1 = (Q-1)/ALPHA where we
*              set c1 = 0 and -ALPHA/2 <= c0 = c mod Q - Q < 0.
*              Assumes coefficients to be standard representatives.
*
* Arguments:   - poly *a1: pointer to output polynomial with coefficients c1
*              - poly *a0: pointer to output polynomial with coefficients c0
*              - const poly *a: pointer to input polynomial
**************************************************/
void poly_decompose_rvv(poly *a1, poly *a0, const poly *a) {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  int32_t* a0_ptr=a0->coeffs;
  int32_t* a1_ptr=a1->coeffs;
  vint32m8_t v_a, v_a1, v_a0, v_temp;
  while(avl>0){
    vl = vsetvl_e32m8(N);
    // load poly a
    v_a = vle32_v_i32m8(a_ptr, vl);
    // get a + 127
    v_a1 = vadd_vx_i32m8(v_a, 127, vl);
    // get a1 = (a + 127) >> 7
    v_a1 = vsra_vx_i32m8(v_a1, 7, vl);

  #if GAMMA2 == (Q-1)/32
    // get a1*1025
    v_a1 = vmul_vx_i32m8(v_a1, 1025, vl);
    // get a1*1025 + (1 << 21)
    v_a1 = vadd_vx_i32m8(v_a1, 1<<21, vl);
    // get a1 = (a1*1025 + (1 << 21)) >> 22
    v_a1 = vsra_vx_i32m8(v_a1, 22, vl);
    // get a1 = a1 & 15
    v_a1 = vand_vx_i32m8(v_a1, 15, vl);

    // store poly a1
    vse32_v_i32m8(a1_ptr, v_a1, vl);

  #elif GAMMA2 == (Q-1)/88
    // get a1*11275
    v_a1 = vmul_vx_i32m8(v_a1, 11275, vl);
    // get a1*11275 + (1 << 23)
    v_a1 = vadd_vx_i32m8(v_a1, 1<<23, vl);
    // get a1 = (a1*1025 + (1 << 21)) >> 24
    v_a1 = vsra_vx_i32m8(v_a1, 24, vl);
    // get 43 - a1
    v_temp = vrsub_vx_i32m8(v_a1, 43, vl);
    // get (43 - a1) >> 31
    v_temp = vsra_vx_i32m8(v_temp, 31, vl);
    // get ((43 - a1) >> 31) & a1
    v_temp = vand_vv_i32m8(v_temp, v_a1, vl);
    // get a1 ^= ((43 - a1) >> 31) & a1
    v_a1 = vxor_vv_i32m8(v_a1, v_temp, vl);

    // store poly a1
    vse32_v_i32m8(a1_ptr, v_a1, vl);
  #endif

    // get a1*2*GAMMA2
    v_a1 = vmul_vx_i32m8(v_a1, 2*GAMMA2, vl);
    // get a0 = a - a1*2*GAMMA2
    v_a0 = vsub_vv_i32m8(v_a, v_a1, vl);
    // get (Q-1)/2 - a0
    v_temp = vrsub_vx_i32m8(v_a0, (Q-1)/2, vl);
    // get ((Q-1)/2 - a0) >> 31
    v_temp = vsra_vx_i32m8(v_temp, 31, vl);
    // get (((Q-1)/2 - *a0) >> 31) & Q
    v_temp = vand_vx_i32m8(v_temp, Q, vl);
    // get a0 -= (((Q-1)/2 - *a0) >> 31) & Q
    v_a0 = vsub_vv_i32m8(v_a0, v_temp, vl);

    // store poly a0
    vse32_v_i32m8(a0_ptr, v_a0, vl);

    a_ptr+=vl;
    a0_ptr+=vl;
    a1_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_make_hint_rvv
*
* Description: Compute hint polynomial. The coefficients of which indicate
*              whether the low bits of the corresponding coefficient of
*              the input polynomial overflow into the high bits.
*
* Arguments:   - poly *h: pointer to output hint polynomial
*              - const poly *a0: pointer to low part of input polynomial
*              - const poly *a1: pointer to high part of input polynomial
*
* Returns number of 1 bits.
**************************************************/
unsigned int poly_make_hint_rvv(poly *h, const poly *a0, const poly *a1) {
  unsigned int s;
  size_t vl;
  size_t avl=N;
  int32_t* h_ptr=h->coeffs;
  const int32_t* a0_ptr=a0->coeffs;
  const int32_t* a1_ptr=a1->coeffs;
  vint32m8_t vx, vy;
  vbool4_t vb0, vb1, vb2, vb3;
  while(avl>0){
    vl = vsetvl_e32m8(avl);  
  
    // load poly a0 to vx
    vx = vle32_v_i32m8(a0_ptr, vl);
    // get a0 <= GAMMA2
    vb0 = vmsle_vx_i32m8_b4(vx, GAMMA2, vl);
    // get a0 > Q - GAMMA2
    vb1 = vmsgt_vx_i32m8_b4(vx, Q - GAMMA2, vl);
    // get a0 == Q - GAMMA2
    vb2 = vmseq_vx_i32m8_b4(vx, Q - GAMMA2, vl);

    // load poly a1 to vy
    vy = vle32_v_i32m8(a1_ptr, vl);
    // get a1 == 0
    vb3 = vmseq_vx_i32m8_b4(vy, 0, vl);

    // get a0 <= GAMMA2 || a0 > Q - GAMMA2
    vb0 = vmor_mm_b4(vb0, vb1, vl);
    // get a0 == Q - GAMMA2 && a1 == 0
    vb2 = vmand_mm_b4(vb2, vb3, vl);
    // get a0 <= GAMMA2 || a0 > Q - GAMMA2 || (a0 == Q - GAMMA2 && a1 == 0)
    vb0 = vmor_mm_b4(vb0, vb2, vl);
    vb0=vmnot_m_b4(vb0,vl);
    //store result into h
    vy=vmv_v_x_i32m8(0,vl);
    vy=vmerge_vxm_i32m8(vb0,vy,1,vl);
    vse32_v_i32m8(h_ptr,vy,vl);
    
    // vmor, vmand, vcpop instruction remains unimplemented in PCIT RVV gem5
    s+=(unsigned int)vcpop_m_b4(vb0, vl);

    //update process
    h_ptr+=vl;
    a1_ptr+=vl;
    a0_ptr+=vl;
    avl-=vl;
  }
  return s;
}

/*************************************************
* Name:        poly_use_hint_rvv
*
* Description: Use hint polynomial to correct the high bits of a polynomial.
*
* Arguments:   - poly *b: pointer to output polynomial with corrected high bits
*              - const poly *a: pointer to input polynomial
*              - const poly *h: pointer to input hint polynomial
**************************************************/
void poly_use_hint_rvv(poly *b, const poly *a, const poly *h) {
  poly a1_poly, a0_poly;
  poly *a1 = &a1_poly;
  poly *a0 = &a0_poly;
  poly_decompose_rvv(a1,a0,a);

  vint32m2_t vx,vy,vz;
  vbool16_t va,vb,vc;
  size_t vl;
  size_t avl=N;
  const int32_t* h_ptr=h->coeffs;
  int32_t* b_ptr=b->coeffs;
  int32_t* a0_ptr=a0_poly.coeffs;
  int32_t* a1_ptr=a1_poly.coeffs;
  while(avl>0){
    vl=vsetvl_e32m2(avl);
    vx=vle32_v_i32m2(a0_ptr,vl);
    vy=vle32_v_i32m2(a1_ptr,vl);
    va=vmsgt_vx_i32m2_b16(vx,0,vl);//a0>0
    vx=vmv_v_v_i32m2(vy,vl);//a1
    vy=vadd_vx_i32m2(vx,1,vl); //a1+1
    vz=vsub_vx_i32m2(vx,1,vl);//a1-1
  #if GAMMA2 == (Q-1)/32
    vz=vand_vx_i32m2(vz,15,vl);//(a1-1)&15
    //__asm__ __volatile__ ( "vand.vi %[vd], %[vt], %[imm]" : [vd] "=vr"(vz) : [vt] "vr"(vz),[imm] "i"(15));
    vy=vand_vx_i32m2(vy,15,vl);//(a1+1)&15
    //__asm__ __volatile__ ( "vand.vi %[vd], %[vt], %[imm]" : [vd] "=vr"(vy) : [vt] "vr"(vy),[imm] "i"(15));
    vy=vmerge_vvm_i32m2(va,vz,vy,vl);//now vy holds final result of this branch
  #elif GAMMA2 == (Q-1)/88
    vb=vmseq_vx_i32m2_b16(vx,43,vl);//a1==43
    vc=vmseq_vx_i32m2_b16(vx,0,vl);//a1==0
    vy=vmerge_vxm_i32m2(vb,vy,0,vl);
    vz=vmerge_vxm_i32m2(vc,vz,43,vl);
    vy=vmerge_vvm_i32m2(va,vz,vy,vl);//now vy holds final result of this branch
  #endif
    vz=vle32_v_i32m2(h_ptr,vl);
    va=vmseq_vx_i32m2_b16(vz,0,vl);//hint==0
    vx=vmerge_vvm_i32m2(va,vy,vx,vl);
    vse32_v_i32m2(b_ptr,vx,vl);

    h_ptr+=vl;
    b_ptr+=vl;
    a0_ptr+=vl;
    a1_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        poly_chknorm_rvv
*
* Description: Check infinity norm of polynomial against given bound.
*              Assumes input coefficients were reduced by reduce32().
*
* Arguments:   - const poly *a: pointer to polynomial, all in range (0,Q)
*              - int32_t B: norm bound
*
* Returns 0 if norm is strictly smaller than B <= (Q-1)/8 and 1 otherwise.
**************************************************/
int poly_chknorm_rvv(const poly *a, int32_t B) {
  if(B > (Q-1)/8)
    return 1;

  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  vint32m4_t vx,vy,vz;
  vbool8_t va;
  uint32_t s=0;
  while(avl>0){
    vl=vsetvl_e32m4(avl);
    vx=vle32_v_i32m4(a_ptr,vl);
    vy=vsra_vx_i32m4(vx,31,vl);
    vz=vmul_vx_i32m4(vx,2,vl);
    vy=vand_vv_i32m4(vy,vz,vl);
    vx=vsub_vv_i32m4(vx,vy,vl);
    va=vmslt_vx_i32m4_b8(vx,B,vl);
    va=vmnot_m_b8(va,vl);
    s+=(uint32_t)vcpop_m_b8(va,vl);
    a_ptr+=vl;
    avl-=vl;
  }
  if(s>0){
    return 1;
  }else{
    return 0;
  }
}

/*************************************************
* Name:        rej_uniform_custom
*
* Description: Sample uniformly random coefficients in [0, Q-1] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_uniform_custom(int32_t *a,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  size_t ctr, pos, vl, valid_num;
  vuint8m1_t vreg_m;
  vint32m1_t vreg_rej;

  csr_validnum_rw();
  ctr=pos=0;

  while((ctr<len)&&(pos<buflen)){
    if((pos + UNPACK_REJ_LOAD_BYTE) < buflen) {
      vl = vsetvl_e8m1(UNPACK_REJ_LOAD_BYTE);
    }
    else {
      vl = vsetvl_e8m1(buflen - pos);
    }
    
    vreg_m = vle8_v_u8m1(buf + pos, vl);
    pos += vl;

    size_t unpacked_len=vl/3;

    vsetvl_e32m1_wrapper(unpacked_len);

    vreg_rej=unpack_vx_i32m1(vreg_m,24);
    vreg_rej=vand_vx_i32m1(vreg_rej,0x7FFFFF,unpacked_len);
    vreg_rej=sample_rej_vx_i32m1(vreg_rej,Q);
    valid_num=csr_validnum_rw();

    if(ctr+valid_num>len){
      valid_num=len-ctr;
    }

    vl=vsetvl_e32m1(valid_num);
    vse32_v_i32m1(a+ctr,vreg_rej,valid_num);
    ctr+=valid_num;
  }
  return ctr;
}

/*************************************************
* Name:        poly_uniform_custom
*
* Description: Sample polynomial with uniformly random coefficients
*              in [0,Q-1] by performing rejection sampling on the
*              output stream of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
#define POLY_UNIFORM_NBLOCKS ((768 + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
void poly_uniform_custom(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce)
{
  unsigned int i, ctr, off;
  unsigned int buflen = POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES + 2];

  stream128_init_custom(seed,nonce);
  stream128_squeezeblocks_custom(buf,POLY_UNIFORM_NBLOCKS,false);

  ctr=rej_uniform_custom(a->coeffs, N, buf, buflen);

  while(ctr<N){
    stream128_squeezeblocks_custom(buf,1,true);
    buflen=STREAM128_BLOCKBYTES;
    ctr+=rej_uniform_custom(a->coeffs + ctr, N - ctr, buf, buflen);
  }
}

/*************************************************
* Name:        rej_eta_custom
*
* Description: Sample uniformly random coefficients in [-ETA, ETA] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_eta_custom(int32_t *a,
                            unsigned int len,
                            const uint8_t *buf,
                            unsigned int buflen)
{
  unsigned int ctr,pos,vl,valid_num;
  vint32m1_t vreg_rej;
  vint32m1_t vtemp;
  vuint8m1_t vreg_m;

  csr_validnum_rw();
  ctr=pos=0;

  while((ctr<len)&&(pos<buflen)){
    if((pos + UNPACK_REJ_ETA_LOAD_BYTE) < buflen) {
      vl = vsetvl_e8m1(UNPACK_REJ_ETA_LOAD_BYTE);
    }
    else {
      vl = vsetvl_e8m1(buflen - pos);
    }

    vreg_m=vle8_v_u8m1(buf+pos,vl);
    pos+=vl;

    size_t unpacked_vl=vl<<1;

    vsetvl_e32m1_wrapper(unpacked_vl);

    vreg_rej=unpack_vx_i32m1(vreg_m,4);
    #if ETA==2
      vreg_rej=sample_rej_vx_i32m1(vreg_rej,15);
      valid_num=csr_validnum_rw();
      if(ctr+valid_num>len){
        valid_num=len-ctr;
      }
      vl=vsetvl_e32m1(valid_num);
      vtemp=vmul_vx_i32m1(vreg_rej,205,vl);
      vtemp=vsra_vx_i32m1(vtemp,10,vl);
      vtemp=vmul_vx_i32m1(vtemp,5,vl);
      vreg_rej=vsub_vv_i32m1(vreg_rej,vtemp,vl);
      vreg_rej=vrsub_vx_i32m1(vreg_rej,2,vl);
    #elif ETA==4
      vreg_rej=sample_rej_vx_i32m1(vreg_rej,9);
      valid_num=csr_validnum_rw();
      if(ctr+valid_num>len){
        valid_num=len-ctr;
      }
      vl=vsetvl_e32m1(valid_num);
      vreg_rej=vrsub_vx_i32m1(vreg_rej,4,vl);
    #endif
    vse32_v_i32m1(a+ctr,vreg_rej,vl);
    ctr+=valid_num;
  }
  return ctr;
}

/*************************************************
* Name:        poly_uniform_eta_custom
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-ETA,ETA] by performing rejection sampling on the
*              output stream from SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
#if ETA == 2
#define POLY_UNIFORM_ETA_NBLOCKS ((136 + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
#elif ETA == 4
#define POLY_UNIFORM_ETA_NBLOCKS ((227 + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
#endif
void poly_uniform_eta_custom(poly *a,
                      const uint8_t seed[SEEDBYTES],
                      uint16_t nonce)
{
  unsigned int ctr;
  unsigned int buflen = POLY_UNIFORM_ETA_NBLOCKS*STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_ETA_NBLOCKS*STREAM128_BLOCKBYTES];

  stream128_init_custom(seed, nonce);
  stream128_squeezeblocks_custom(buf, POLY_UNIFORM_ETA_NBLOCKS, false);

  ctr = rej_eta_custom(a->coeffs, N, buf, buflen);

  while(ctr<N){
    stream128_squeezeblocks_custom(buf, 1, true);
    ctr += rej_eta_custom(a->coeffs + ctr, N - ctr, buf, STREAM128_BLOCKBYTES);
  }
}

/*************************************************
* Name:        poly_uniform_gamma1m1_custom
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-(GAMMA1 - 1), GAMMA1] by unpacking output stream
*              of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length CRHBYTES
*              - uint16_t nonce: 16-bit nonce
**************************************************/
#if GAMMA1 == (1 << 17)
#define POLY_UNIFORM_GAMMA1_NBLOCKS ((576 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#elif GAMMA1 == (1 << 19)
#define POLY_UNIFORM_GAMMA1_NBLOCKS ((640 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#endif
void poly_uniform_gamma1_custom(poly *a,
                         const uint8_t seed[CRHBYTES],
                         uint16_t nonce)
{
  uint8_t buf[POLY_UNIFORM_GAMMA1_NBLOCKS*STREAM256_BLOCKBYTES];

  stream256_init_custom(seed, nonce);
  stream256_squeezeblocks_custom(buf, POLY_UNIFORM_GAMMA1_NBLOCKS, false);
  polyz_unpack_custom(a, buf);
}

/*************************************************
* Name:        poly_challenge_custom
*
* Description: Implementation of H. Samples polynomial with TAU nonzero
*              coefficients in {-1,1} using the output stream of
*              SHAKE256(seed).
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const uint8_t mu[]: byte array containing seed of length SEEDBYTES
**************************************************/
void poly_challenge_custom(poly *c, const uint8_t seed[SEEDBYTES]) {
  unsigned int i, b, pos;
  uint64_t signs;
  uint8_t buf[SHAKE256_RATE];
  keccak_state state;

  // shake256_init(&state);
  // shake256_absorb(&state, seed, SEEDBYTES);
  // shake256_finalize(&state);
  // shake256_squeezeblocks(buf, 1, &state);
  csr_keccakmode_rw(SHAKE_256_MODE);
  keccak_absorb_custom(seed,SEEDBYTES, SHAKE256_RATE);
  keccak_squeezeblocks_custom(buf, 1, SHAKE256_RATE, false);

  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (uint64_t)buf[i] << 8*i;
  pos = 8;

  for(i = 0; i < N; ++i)
    c->coeffs[i] = 0;
  for(i = N-TAU; i < N; ++i) {
    do {
      if(pos >= SHAKE256_RATE) {
        // shake256_squeezeblocks(buf, 1, &state);
        keccak_squeezeblocks_custom(buf, 1, SHAKE256_RATE, true);
        pos = 0;
      }

      b = buf[pos++];
    } while(b > i);

    c->coeffs[i] = c->coeffs[b];
    c->coeffs[b] = 1 - 2*(signs & 1);
    signs >>= 1;
  }
}

/*************************************************
* Name:        polyeta_pack_custom
*
* Description: Bit-pack polynomial with coefficients in [-ETA,ETA].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYETA_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyeta_pack_custom(uint8_t *r, const poly *a) {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  uint8_t* r_ptr=r;
  vint32m1_t vt;
  vuint8m1_t vx;
  while(avl>0){
    vl=vsetvl_e32m1(avl);
    vt=vle32_v_i32m1(a_ptr,vl);
    vt=vrsub_vx_i32m1(vt,ETA,vl);
  #if ETA == 2
    size_t vl_packed=(vl*3)>>3;
    vx=pack_vx_i32m1(vt,3);
  #elif ETA == 4
    size_t vl_packed=vl>>1;
    vx=pack_vx_i32m1(vt,4);
  #endif
    a_ptr+=vl;
    avl-=vl;
    vse8_v_u8m1(r_ptr,vx,vl_packed);
    r_ptr+=vl_packed;
  }
}

/*************************************************
* Name:        polyeta_unpack_custom
*
* Description: Unpack polynomial with coefficients in [-ETA,ETA].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyeta_unpack_custom(poly *r, const uint8_t *a) {
  size_t vl;
  size_t avl=N;
  int32_t* r_ptr=r->coeffs;
  uint8_t* a_ptr=a;
  vint32m1_t vt;
  vuint8m1_t vx;
  while(avl>0){
    vl=vsetvl_e32m1(avl);
  #if ETA == 2
    size_t vl_packed=(vl*3)>>3;
    vx=vle8_v_u8m1(a_ptr,vl_packed);
    vsetvl_e32m1_wrapper(avl);
    vt=unpack_vx_i32m1(vx,3);
  #elif ETA == 4
    size_t vl_packed=vl>>1;
    vx=vle8_v_u8m1(a_ptr,vl_packed);
    vsetvl_e32m1_wrapper(avl);
    vt=unpack_vx_i32m1(vx,4);
  #endif
    vt=vrsub_vx_i32m1(vt,ETA,vl);
    vse32_v_i32m1(r_ptr,vt,vl);
    a_ptr+=vl_packed;
    r_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        polyt1_pack_custom
*
* Description: Bit-pack polynomial t1 with coefficients fitting in 10 bits.
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYT1_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyt1_pack_custom(uint8_t *r, const poly *a) {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  uint8_t* r_ptr=r;
  vint32m1_t vt;
  vuint8m1_t vx;
  while(avl>0){
    vl=vsetvl_e32m1(avl);
    vt=vle32_v_i32m1(a_ptr,vl);
    size_t vl_packed=(vl*10)>>3;
    vx=pack_vx_i32m1(vt,10);
    a_ptr+=vl;
    avl-=vl;
    vse8_v_u8m1(r_ptr,vx,vl_packed);
    r_ptr+=vl_packed;
  }
}

/*************************************************
* Name:        polyt1_unpack_custom
*
* Description: Unpack polynomial t1 with 10-bit coefficients.
*              Output coefficients are standard representatives.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyt1_unpack_custom(poly *r, const uint8_t *a) {
  size_t vl;
  size_t avl=N;
  int32_t* r_ptr=r->coeffs;
  uint8_t* a_ptr=a;
  vint32m1_t vt;
  vuint8m1_t vx;
  while(avl>0){
    vl=vsetvl_e32m1(avl);
    size_t vl_packed=(vl*10)>>3;
    vx=vle8_v_u8m1(a_ptr,vl_packed);
    vsetvl_e32m1_wrapper(avl);
    vt=unpack_vx_i32m1(vx,10);
    vse32_v_i32m1(r_ptr,vt,vl);
    r_ptr+=vl;
    a_ptr+=vl_packed;
    avl-=vl;
  }
}

/*************************************************
* Name:        polyt0_pack_custom
*
* Description: Bit-pack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYT0_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyt0_pack_custom(uint8_t *r, const poly *a) {
  size_t vl;
  size_t avl=N;
  int32_t* a_ptr=a->coeffs;
  uint8_t* r_ptr=r;
  vint32m1_t vt;
  vuint8m1_t vx;
  int32_t cons=1 << (D-1);
  while(avl>0){
    vl=vsetvl_e32m1(avl);
    vt=vle32_v_i32m1(a_ptr,vl);
    vt=vrsub_vx_i32m1(vt,cons,vl);
    size_t vl_packed=(vl*13)>>3;
    vx=pack_vx_i32m1(vt,13);
    a_ptr+=vl;
    avl-=vl;
    vse8_v_u8m1(r_ptr,vx,vl_packed);
    r_ptr+=vl_packed;
  }
}

/*************************************************
* Name:        polyt0_unpack_custom
*
* Description: Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyt0_unpack_custom(poly *r, const uint8_t *a) {
  size_t vl;
  size_t avl=N;
  int32_t* r_ptr=r->coeffs;
  uint8_t* a_ptr=a;
  vint32m1_t vt;
  vuint8m1_t vx;
  int32_t cons=1 << (D-1);
  while(avl>0){
    vl=vsetvl_e32m1(avl);
    size_t vl_packed=(vl*13)>>3;
    vx=vle8_v_u8m1(a_ptr,vl_packed);
    vsetvl_e32m1_wrapper(avl);
    vt=unpack_vx_i32m1(vx,13);
    vt=vrsub_vx_i32m1(vt,cons,vl);
    vse32_v_i32m1(r_ptr,vt,vl);
    a_ptr+=vl_packed;
    r_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        polyz_pack_custom
*
* Description: Bit-pack polynomial with coefficients
*              in [-(GAMMA1 - 1), GAMMA1].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYZ_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyz_pack_custom(uint8_t *r, const poly *a) {
  size_t vl;
  size_t avl=N;
  uint8_t* r_ptr=r;
  int32_t* a_ptr=a->coeffs;
  vint32m1_t vt;
  vuint8m1_t vx;
  while(avl>0){
    vl=vsetvl_e32m1(avl);
    vt=vle32_v_i32m1(a_ptr,vl);
    vt=vrsub_vx_i32m1(vt,GAMMA1,vl);
    #if GAMMA1 == (1 << 17)
      size_t vl_packed=(vl*18)>>3;
      vx=pack_vx_i32m1(vt,18);
    #elif GAMMA1 == (1 << 19)
      size_t vl_packed=(vl*20)>>3;
      vx=pack_vx_i32m1(vt,20);
    #endif
    a_ptr+=vl;
    avl-=vl;
    vse8_v_u8m1(r_ptr,vx,vl_packed);
    r_ptr+=vl_packed;
  }
}

/*************************************************
* Name:        polyz_unpack_custom
*
* Description: Unpack polynomial z with coefficients
*              in [-(GAMMA1 - 1), GAMMA1].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyz_unpack_custom(poly *r, const uint8_t *a) {
  size_t vl;
  size_t avl=N;
  int32_t* r_ptr=r->coeffs;
  uint8_t* a_ptr=a;
  vint32m1_t vt;
  vuint8m1_t vx;
  while(avl>0){
    vl=vsetvl_e32m1(avl);
    #if GAMMA1 == (1 << 17)
      size_t vl_packed=(vl*18)>>3;
      vx=vle8_v_u8m1(a_ptr,vl_packed);
      vsetvl_e32m1_wrapper(avl);
      vt=unpack_vx_i32m1(vx,18);
    #elif GAMMA1 == (1 << 19)
      size_t vl_packed=(vl*20)>>3;
      vx=vle8_v_u8m1(a_ptr,vl_packed);
      vsetvl_e32m1_wrapper(avl);
      vt=unpack_vx_i32m1(vx,20);
    #endif
    vt=vrsub_vx_i32m1(vt,GAMMA1,vl);
    vse32_v_i32m1(r_ptr,vt,vl);
    a_ptr+=vl_packed;
    r_ptr+=vl;
    avl-=vl;
  }
}

/*************************************************
* Name:        polyw1_pack_custom
*
* Description: Bit-pack polynomial w1 with coefficients in [0,15] or [0,43].
*              Input coefficients are assumed to be standard representatives.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYW1_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyw1_pack_custom(uint8_t *r, const poly *a) {
  size_t vl;
  size_t avl=N;
  uint8_t* r_ptr=r;
  int32_t* a_ptr=a->coeffs;
  vint32m1_t vt;
  vuint8m1_t vx;
  while(avl>0){
    vl=vsetvl_e32m1(avl);
    vt=vle32_v_i32m1(a_ptr,vl);
    #if GAMMA2 == (Q-1)/88
      size_t vl_packed=(vl*6)>>3;
      vx=pack_vx_i32m1(vt,6);
    #elif GAMMA2 == (Q-1)/32
      size_t vl_packed=vl>>1;
      vx=pack_vx_i32m1(vt,4);
    #endif
    a_ptr+=vl;
    avl-=vl;
    vse8_v_u8m1(r_ptr,vx,vl_packed);
    r_ptr+=vl_packed;
  }
}

/******************************************************************
 *              ASM OPTIMIZED VERSIONS 
*****************************************************************/
/*************************************************
* Name:        poly_ntt_custom_asm
*
* Description: Outof-place NTT(Constant Geometry NTT). Input and output are both in normal order.
*              Coefficients of both input and output polynomial range in (0,Q).
*              Coefficients of both input and output polynomial are in normal domain (not montgomery domain).
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_ntt_custom_asm(poly *a) {
  ntt_cg_custom_reorder_asm(a);
}

/*************************************************
* Name:        poly_invntt_custom_asm
*
* Description: Outof-place inverse NTT. Input and output are both in normal order.
*              Coefficients of both input and output polynomial range in (0,Q).
*              Input and output polynomial are both in montgomery domain.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invntt_custom_asm(poly *a) {
  invntt_cg_custom_reorder_asm(a);
}