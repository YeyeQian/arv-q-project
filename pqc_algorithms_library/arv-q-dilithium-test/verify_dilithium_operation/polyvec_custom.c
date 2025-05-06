#include <stdint.h>
#include "param.h"
#include "polyvec.h"
#include "poly.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

/*************************************************
* Name:        expand_mat_custom
*
* Description: Implementation of ExpandA. Generates matrix A with uniformly
*              random coefficients a_{i,j} by performing rejection
*              sampling on the output stream of SHAKE128(rho|j|i)
*              or AES256CTR(rho,j|i).
*
* Arguments:   - polyvecl mat[K]: output matrix
*              - const uint8_t rho[]: byte array containing seed rho
**************************************************/
void polyvec_matrix_expand_custom(polyvecl mat[K], const uint8_t rho[SEEDBYTES]) {
  unsigned int i, j;

  for(i = 0; i < K; ++i)
    for(j = 0; j < L; ++j)
      poly_uniform_custom(&mat[i].vec[j], rho, (i << 8) + j);
}

void polyvec_matrix_pointwise_montgomery_custom(polyveck *t, const polyvecl mat[K], const polyvecl *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    polyvecl_pointwise_acc_montgomery_custom(&t->vec[i], &mat[i], v);
}

void polyvec_matrix_chorder(polyvecl mat[K]){
  unsigned int i;
  for(i = 0; i < K; ++i)
    polyvecl_chorder(&mat[i]);
}

void polyvec_matrix_print(const polyvecl mat[K]){
  unsigned int i;
  for(i = 0; i < K; ++i)
    polyvecl_print(&mat[i]);
}

/**************************************************************/
/************ Vectors of polynomials of length L **************/
/**************************************************************/

void polyvecl_print(const polyvecl *v){
  for(int i=0;i<L;i++){
    poly_print(&v->vec[i]);
  }
}

void polyvecl_chorder(polyvecl *v){
  for(int i=0;i<L;i++){
    poly_chorder(&v->vec[i]);
  }
}

void polyvecl_mon2nor_custom(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_mon2nor_custom(&v->vec[i]);
}

void polyvecl_uniform_eta_custom(polyvecl *v, const uint8_t seed[SEEDBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_uniform_eta_custom(&v->vec[i], seed, nonce++);
}

void polyvecl_uniform_gamma1_custom(polyvecl *v, const uint8_t seed[CRHBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_uniform_gamma1_custom(&v->vec[i], seed, L*nonce + i);
}

void polyvecl_reduce_rvv(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_reduce_rvv(&v->vec[i]);
}

void polyvecl_caddq_custom(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_caddq_custom(&v->vec[i]);
}

/*************************************************
* Name:        polyvecl_freeze_rvv
*
* Description: Reduce coefficients of polynomials in vector of length L
*              to standard representatives.
*
* Arguments:   - polyvecl *v: pointer to input/output vector
**************************************************/
void polyvecl_freeze_rvv(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_freeze_rvv(&v->vec[i]);
}

/*************************************************
* Name:        polyvecl_add_custom
*
* Description: Add vectors of polynomials of length L.
*              Modular reduction is performed.
*
* Arguments:   - polyvecl *w: pointer to output vector
*              - const polyvecl *u: pointer to first summand
*              - const polyvecl *v: pointer to second summand
**************************************************/
void polyvecl_add_custom(polyvecl *w, const polyvecl *u, const polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_add_custom(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyvecl_ntt_custom
*
* Description: Forward NTT of all polynomials in vector of length L. Output
*              coefficients range in (0,Q).
*
* Arguments:   - polyvecl *v: pointer to input/output vector
**************************************************/
void polyvecl_ntt_custom(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_ntt_custom(&v->vec[i]);
}

void polyvecl_ntt_custom_asm(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_ntt_custom_asm(&v->vec[i]);
}

void polyvecl_invntt_custom(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_invntt_custom(&v->vec[i]);
}

void polyvecl_invntt_custom_asm(polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_invntt_custom_asm(&v->vec[i]);
}

void polyvecl_pointwise_poly_montgomery_custom(polyvecl *r, const poly *a, const polyvecl *v) {
  unsigned int i;

  for(i = 0; i < L; ++i)
    poly_pointwise_montgomery_custom(&r->vec[i], a, &v->vec[i]);
}

/*************************************************
* Name:        polyvecl_pointwise_acc_montgomery_custom
*
* Description: Pointwise multiply vectors of polynomials of length L, multiply
*              resulting vector by 2^{-32} and add (accumulate) polynomials
*              in it. Input/output vectors are in NTT domain representation.
*
* Arguments:   - poly *w: output polynomial
*              - const polyvecl *u: pointer to first input vector
*              - const polyvecl *v: pointer to second input vector
**************************************************/
void polyvecl_pointwise_acc_montgomery_custom(poly *w,
                                       const polyvecl *u,
                                       const polyvecl *v)
{
  unsigned int i;
  poly t;

  poly_pointwise_montgomery_custom(w, &u->vec[0], &v->vec[0]);
  for(i = 1; i < L; ++i) {
    poly_pointwise_montgomery_custom(&t, &u->vec[i], &v->vec[i]);
    poly_add_custom(w, w, &t);
  }
}

/*************************************************
* Name:        polyvecl_chknorm_rvv
*
* Description: Check infinity norm of polynomials in vector of length L.
*              Assumes input polyvecl to be reduced by polyvecl_reduce().
*
* Arguments:   - const polyvecl *v: pointer to vector
*              - int32_t B: norm bound
*
* Returns 0 if norm of all polynomials is strictly smaller than B <= (Q-1)/8
* and 1 otherwise.
**************************************************/
int polyvecl_chknorm_rvv(const polyvecl *v, int32_t bound)  {
  unsigned int i;

  for(i = 0; i < L; ++i)
    if(poly_chknorm_rvv(&v->vec[i], bound))
      return 1;

  return 0;
}

/**************************************************************/
/************ Vectors of polynomials of length K **************/
/**************************************************************/

void polyveck_print(const polyveck *v){
  for(int i=0;i<K;i++){
    poly_print(&v->vec[i]);
  }
}

void polyveck_chorder(polyveck *v){
  for(int i=0;i<K;i++){
    poly_chorder(&v->vec[i]);
  }
}

void polyveck_copy(polyveck *dst_v, const polyveck *src_v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_copy(&dst_v->vec[i],&src_v->vec[i]);
}

void polyveck_mon2nor_custom(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_mon2nor_custom(&v->vec[i]);
}

void polyveck_uniform_eta_custom(polyveck *v, const uint8_t seed[SEEDBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_uniform_eta_custom(&v->vec[i], seed, nonce++);
}

/*************************************************
* Name:        polyveck_reduce_rvv
*
* Description: Reduce coefficients of polynomials in vector of length K
*              to representatives in [-6283009,6283007].
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_reduce_rvv(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_reduce_rvv(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_caddq_rvv
*
* Description: For all coefficients of polynomials in vector of length K
*              add Q if coefficient is negative.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_caddq_rvv(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_caddq_rvv(&v->vec[i]);
}

void polyveck_caddq_custom(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_caddq_custom(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_freeze_rvv
*
* Description: Reduce coefficients of polynomials in vector of length K
*              to standard representatives.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_freeze_rvv(polyveck *v)  {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_freeze_rvv(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_add_custom
*
* Description: Add vectors of polynomials of length K.
*              Modular reduction is performed.
*
* Arguments:   - polyveck *w: pointer to output vector
*              - const polyveck *u: pointer to first summand
*              - const polyveck *v: pointer to second summand
**************************************************/
void polyveck_add_custom(polyveck *w, const polyveck *u, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_add_custom(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyveck_add_rvv(polyveck *w, const polyveck *u, const polyveck *v) {//no modular reduction
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_add_rvv(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_sub_custom
*
* Description: Subtract vectors of polynomials of length K.
*              Modular reduction is performed.
*
* Arguments:   - polyveck *w: pointer to output vector
*              - const polyveck *u: pointer to first input vector
*              - const polyveck *v: pointer to second input vector to be
*                                   subtracted from first input vector
**************************************************/
void polyveck_sub_custom(polyveck *w, const polyveck *u, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_sub_custom(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_shiftl_rvv
*
* Description: Multiply vector of polynomials of Length K by 2^D without modular
*              reduction. Assumes input coefficients to be less than 2^{31-D}.
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_shiftl_rvv(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_shiftl_rvv(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_ntt_custom
*
* Description: Forward NTT of all polynomials in vector of length K. Output
*              coefficients range in (0,Q).
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_ntt_custom(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_ntt_custom(&v->vec[i]);
}

void polyveck_ntt_custom_asm(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_ntt_custom_asm(&v->vec[i]);
}

/*************************************************
* Name:        polyveck_invntt_tomont_custom
*
* Description: Inverse NTT and multiplication by 2^{32} of polynomials
*              in vector of length K. Input coefficients range in(0,Q).
*
* Arguments:   - polyveck *v: pointer to input/output vector
**************************************************/
void polyveck_invntt_custom(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_invntt_custom(&v->vec[i]);
}

void polyveck_invntt_custom_asm(polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_invntt_custom_asm(&v->vec[i]);
}


void polyveck_pointwise_poly_montgomery_custom(polyveck *r, const poly *a, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_pointwise_montgomery_custom(&r->vec[i], a, &v->vec[i]);
}

/*************************************************
* Name:        polyveck_chknorm_rvv
*
* Description: Check infinity norm of polynomials in vector of length K.
*              Assumes input polyveck to be reduced by polyveck_reduce().
*
* Arguments:   - const polyveck *v: pointer to vector
*              - int32_t B: norm bound
*
* Returns 0 if norm of all polynomials are strictly smaller than B <= (Q-1)/8
* and 1 otherwise.
**************************************************/
int polyveck_chknorm_rvv(const polyveck *v, int32_t bound) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    if(poly_chknorm_rvv(&v->vec[i], bound))
      return 1;

  return 0;
}

/*************************************************
* Name:        polyveck_power2round_rvv
*
* Description: For all coefficients a of polynomials in vector of length K,
*              compute a0, a1 such that a mod^+ Q = a1*2^D + a0
*              with -2^{D-1} < a0 <= 2^{D-1}. Assumes coefficients to be
*              standard representatives.
*
* Arguments:   - polyveck *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyveck *v0: pointer to output vector of polynomials with
*                              coefficients a0
*              - const polyveck *v: pointer to input vector
**************************************************/
void polyveck_power2round_rvv(polyveck *v1, polyveck *v0, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_power2round_rvv(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_decompose_rvv
*
* Description: For all coefficients a of polynomials in vector of length K,
*              compute high and low bits a0, a1 such a mod^+ Q = a1*ALPHA + a0
*              with -ALPHA/2 < a0 <= ALPHA/2 except a1 = (Q-1)/ALPHA where we
*              set a1 = 0 and -ALPHA/2 <= a0 = a mod Q - Q < 0.
*              Assumes coefficients to be standard representatives.
*
* Arguments:   - polyveck *v1: pointer to output vector of polynomials with
*                              coefficients a1
*              - polyveck *v0: pointer to output vector of polynomials with
*                              coefficients a0
*              - const polyveck *v: pointer to input vector
**************************************************/
void polyveck_decompose_rvv(polyveck *v1, polyveck *v0, const polyveck *v) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_decompose_rvv(&v1->vec[i], &v0->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyveck_make_hint_rvv
*
* Description: Compute hint vector.
*
* Arguments:   - polyveck *h: pointer to output vector
*              - const polyveck *v0: pointer to low part of input vector
*              - const polyveck *v1: pointer to high part of input vector
*
* Returns number of 1 bits.
**************************************************/
unsigned int polyveck_make_hint_rvv(polyveck *h,
                                const polyveck *v0,
                                const polyveck *v1)
{
  unsigned int i, s = 0;

  for(i = 0; i < K; ++i)
    s += poly_make_hint_rvv(&h->vec[i], &v0->vec[i], &v1->vec[i]);

  return s;
}

/*************************************************
* Name:        polyveck_use_hint_rvv
*
* Description: Use hint vector to correct the high bits of input vector.
*
* Arguments:   - polyveck *w: pointer to output vector of polynomials with
*                             corrected high bits
*              - const polyveck *u: pointer to input vector
*              - const polyveck *h: pointer to input hint vector
**************************************************/
void polyveck_use_hint_rvv(polyveck *w, const polyveck *u, const polyveck *h) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    poly_use_hint_rvv(&w->vec[i], &u->vec[i], &h->vec[i]);
}

void polyveck_pack_w1_custom(uint8_t r[K*POLYW1_PACKEDBYTES], const polyveck *w1) {
  unsigned int i;

  for(i = 0; i < K; ++i)
    polyw1_pack_custom(&r[i*POLYW1_PACKEDBYTES], &w1->vec[i]);
}