#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "param.h"
#include "sign.h"
#include "sign_custom_asm.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "poly_custom_rorg_asm.h"
#include "symmetric.h"
#include "fips202.h"

/*************************************************
* Name:        crypto_sign_keypair_custom_asm
*
* Description: Generates public and private key.
*
* Arguments:   - uint8_t *pk: pointer to output public key (allocated
*                             array of CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key (allocated
*                             array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_sign_keypair_custom_asm(uint8_t *pk, uint8_t *sk) {
  uint8_t seedbuf[3*SEEDBYTES];
  uint8_t tr[CRHBYTES];
  const uint8_t *rho, *rhoprime, *key;
  polyvecl mat[K];
  polyvecl s1;
  polyveck s2, t1, t0;
  uint16_t nonce=0;

  /* Get randomness for rho, rhoprime and key */
  //randombytes(seedbuf, SEEDBYTES);//temporarily replaced to avoid usage of AES
  for(int i = 0; i < SEEDBYTES; i++) {
    seedbuf[i] = i;
  }
  shake256_custom(seedbuf, 3*SEEDBYTES, seedbuf, SEEDBYTES);
  rho = seedbuf;
  rhoprime = seedbuf + SEEDBYTES;
  key = seedbuf + 2*SEEDBYTES;

  //s1: sample generation+pack+mod_add+ntt
  uint8_t* sk_s1=sk+SEEDBYTES+SEEDBYTES+CRHBYTES;//pointer to packed byte start position for s1 in sk
  for(int i=0;i<L;i++){
    poly_uniform_eta_custom(&s1.vec[i],rhoprime,nonce++);
    // polyeta_pack_custom(sk_s1+i*POLYETA_PACKEDBYTES,&s1.vec[i]);//pack s1 into sk
    // poly_caddq_custom(&s1.vec[i]);
    // poly_ntt_custom_asm(&s1.vec[i]);
    poly_eta_pack_caddq_ntt_asm(&s1.vec[i],sk_s1+i*POLYETA_PACKEDBYTES);
  }
  //Gradually generate matrix row and multiply it with s1
  for(int i=0;i<K;i++){
    //generate a row of the matrix
    for(int j=0;j<L;j++){
      poly_uniform_custom(&mat[i].vec[j],rho,(i << 8) + j);
    }
    //multiply this row with s1
    // polyvecl_pointwise_acc_montgomery_custom(&t1.vec[i],&mat[i],&s1);
    // poly_invntt_custom(&t1.vec[i]);
    // poly_mon2nor_custom(&t1.vec[i]);
    polyvecl_macc_invntt_mon2nor_asm(&t1.vec[i],&mat[i],&s1);
  }
  
  //s2:generation+pack+mod_add
  uint8_t* sk_s2=sk_s1+L*POLYETA_PACKEDBYTES;//pointer to packed byte start position for s2 in sk
  for(int i=0;i<K;i++){
    poly_uniform_eta_custom(&s2.vec[i],rhoprime,nonce++);
    // polyeta_pack_custom(sk_s2+i*POLYETA_PACKEDBYTES,&s2.vec[i]);//pack s2 into sk
    // poly_caddq_custom(&s2.vec[i]);
    poly_eta_pack_caddq_asm(&s2.vec[i],sk_s2+i*POLYETA_PACKEDBYTES);
  }

  polyveck_add_custom(&t1, &t1, &s2);

  polyveck_power2round_rvv(&t1, &t0, &t1);
  pack_pk_custom(pk, rho, &t1);

  crh_custom(tr, pk, CRYPTO_PUBLICKEYBYTES);

  //Finish the remaining sk pack
  for(int i=0;i<SEEDBYTES;i++) sk[i]=rho[i];
  sk+=SEEDBYTES;
  for(int i=0;i<SEEDBYTES;i++)sk[i]=key[i];
  sk+=SEEDBYTES;
  for(int i=0;i<CRHBYTES;i++)sk[i]=tr[i];
  sk+=CRHBYTES;
  uint8_t* sk_t0=sk_s2+K*POLYETA_PACKEDBYTES;//pointer to packed byte start position for t0 in sk
  for(int i=0;i<K;i++)polyt0_pack_custom(sk_t0+i*POLYT0_PACKEDBYTES,&t0.vec[i]);

  return 0;
}

/*************************************************
* Name:        crypto_sign_keypair_custom_asm_fortest
*
* Description: Generates public and private key.
*
* Arguments:   - uint8_t *pk: pointer to output public key (allocated
*                             array of CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key (allocated
*                             array of CRYPTO_SECRETKEYBYTES bytes)
*
* Returns 0 (success)
**************************************************/
int crypto_sign_keypair_custom_asm_fortest(uint8_t *pk, uint8_t *sk) {
  uint8_t seedbuf[3*SEEDBYTES];
  uint8_t tr[CRHBYTES];
  const uint8_t *rho, *rhoprime, *key;
  polyvecl mat[K];
  polyvecl s1;
  polyveck s2, t1, t0;
  uint16_t nonce=0;

  /* Get randomness for rho, rhoprime and key */
  //randombytes(seedbuf, SEEDBYTES);//temporarily replaced to avoid usage of AES
  for(int i = 0; i < SEEDBYTES; i++) {
    seedbuf[i] = i;
  }
  shake256_custom(seedbuf, 3*SEEDBYTES, seedbuf, SEEDBYTES);
  rho = seedbuf;
  rhoprime = seedbuf + SEEDBYTES;
  key = seedbuf + 2*SEEDBYTES;

  //s1: sample generation+pack+mod_add+ntt
  uint8_t* sk_s1=sk+SEEDBYTES+SEEDBYTES+CRHBYTES;//pointer to packed byte start position for s1 in sk
  for(int i=0;i<L;i++){
    poly_uniform_eta_custom(&s1.vec[i],rhoprime,nonce++);
    // polyeta_pack_custom(sk_s1+i*POLYETA_PACKEDBYTES,&s1.vec[i]);//pack s1 into sk
    // poly_caddq_custom(&s1.vec[i]);
    // poly_ntt_custom_asm(&s1.vec[i]);
    poly_eta_pack_caddq_ntt_asm(&s1.vec[i],sk_s1+i*POLYETA_PACKEDBYTES);
  }
  //Gradually generate matrix row and multiply it with s1
  for(int i=0;i<K;i++){
    //generate a row of the matrix
    for(int j=0;j<L;j++){
      poly_uniform_custom(&mat[i].vec[j],rho,(i << 8) + j);
      poly_chorder(&mat[i].vec[j]);//for consistency test only
    }
    //multiply this row with s1
    // polyvecl_pointwise_acc_montgomery_custom(&t1.vec[i],&mat[i],&s1);
    // poly_invntt_custom(&t1.vec[i]);
    // poly_mon2nor_custom(&t1.vec[i]);
    polyvecl_macc_invntt_mon2nor_asm(&t1.vec[i],&mat[i],&s1);
  }
  
  //s2:generation+pack+mod_add
  uint8_t* sk_s2=sk_s1+L*POLYETA_PACKEDBYTES;//pointer to packed byte start position for s2 in sk
  for(int i=0;i<K;i++){
    poly_uniform_eta_custom(&s2.vec[i],rhoprime,nonce++);
    // polyeta_pack_custom(sk_s2+i*POLYETA_PACKEDBYTES,&s2.vec[i]);//pack s2 into sk
    // poly_caddq_custom(&s2.vec[i]);
    poly_eta_pack_caddq_asm(&s2.vec[i],sk_s2+i*POLYETA_PACKEDBYTES);
  }

  polyveck_add_custom(&t1, &t1, &s2);

  polyveck_power2round_rvv(&t1, &t0, &t1);
  pack_pk_custom(pk, rho, &t1);

  crh_custom(tr, pk, CRYPTO_PUBLICKEYBYTES);

  //Finish the remaining sk pack
  for(int i=0;i<SEEDBYTES;i++) sk[i]=rho[i];
  sk+=SEEDBYTES;
  for(int i=0;i<SEEDBYTES;i++)sk[i]=key[i];
  sk+=SEEDBYTES;
  for(int i=0;i<CRHBYTES;i++)sk[i]=tr[i];
  sk+=CRHBYTES;
  uint8_t* sk_t0=sk_s2+K*POLYETA_PACKEDBYTES;//pointer to packed byte start position for t0 in sk
  for(int i=0;i<K;i++)polyt0_pack_custom(sk_t0+i*POLYT0_PACKEDBYTES,&t0.vec[i]);

  return 0;
}

/*************************************************
* Name:        crypto_sign_signature_custom_asm
*
* Description: Computes signature.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign_signature_custom_asm(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk)
{
  unsigned int n;
  uint8_t seedbuf[2*SEEDBYTES + 3*CRHBYTES];
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  polyvecl mat[K], s1, y, z;
  polyveck t0, s2, w1, w0, h;
  poly cp;

  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + CRHBYTES;
  mu = key + SEEDBYTES;
  rhoprime = mu + CRHBYTES;

  uint8_t* sk_tr=sk+(SEEDBYTES<<1);//pointer to packed byte start position for tr in sk

  memcpy(rho,sk,SEEDBYTES);//unpack for rho
  polyvec_matrix_expand_custom(mat, rho);//think of mat in normal order in the first place
  //s1:unpack+mod_add+ntt
  uint8_t* sk_s1=sk_tr+CRHBYTES;
  for(int i=0;i<L;i++){
    // polyeta_unpack_custom(&s1.vec[i],sk_s1+i*POLYETA_PACKEDBYTES);
    // poly_caddq_custom(&s1.vec[i]);
    // poly_ntt_custom_asm(&s1.vec[i]);
    poly_eta_unpack_caddq_ntt_asm(&s1.vec[i],sk_s1+i*POLYETA_PACKEDBYTES);
  }
  //s2:unpack+mod_add+ntt
  uint8_t* sk_s2=sk_s1+L*POLYETA_PACKEDBYTES;
  for(int i=0;i<K;i++){
    // polyeta_unpack_custom(&s2.vec[i],sk_s2+i*POLYETA_PACKEDBYTES);
    // poly_caddq_custom(&s2.vec[i]);
    // poly_ntt_custom_asm(&s2.vec[i]);
    poly_eta_unpack_caddq_ntt_asm(&s2.vec[i],sk_s2+i*POLYETA_PACKEDBYTES);
  }
  //t0:unpack+mod_add+ntt
  uint8_t* sk_t0=sk_s2+K*POLYETA_PACKEDBYTES;
  for(int i=0;i<K;i++){
    // polyt0_unpack_custom(&t0.vec[i],sk_t0+i*POLYT0_PACKEDBYTES);
    // poly_caddq_custom(&t0.vec[i]);
    // poly_ntt_custom_asm(&t0.vec[i]);
    poly_t0_unpack_caddq_ntt_asm(&t0.vec[i],sk_t0+i*POLYT0_PACKEDBYTES);
  }

  /* Compute CRH(tr, msg) */
  size_t tr_msg_len=mlen+CRHBYTES;
  uint8_t* tr_msg=(uint8_t*)malloc(tr_msg_len*sizeof(uint8_t));
  memcpy(tr_msg,sk_tr,CRHBYTES);//skip unpack for tr here
  memcpy(tr_msg+CRHBYTES,m,mlen);
  shake256_custom(mu,CRHBYTES,tr_msg,tr_msg_len);

  memcpy(key,sk+SEEDBYTES,SEEDBYTES);//unpack for key

#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else
  crh_custom(rhoprime, key, SEEDBYTES + CRHBYTES);//consumes both key and mu
#endif

rej:
  /* Sample intermediate vector y */
  //y:Generation+mod_add+store into y+ntt
  for(int i=0;i<L;i++){
    // poly_uniform_gamma1_custom(&y.vec[i],rhoprime,L*nonce+i);
    // poly_caddq_custom(&y.vec[i]);
    // z.vec[i]=y.vec[i];
    // poly_ntt_custom_asm(&z.vec[i]);
    poly_gamma1_sample_caddq_ntt_asm(&y.vec[i],&z.vec[i],rhoprime,L*nonce+i);
  }
  nonce++;

  /* Matrix-vector multiplication */
  //Gradually multiply each row of matrix with z
  //w1:macc+invntt+mon2nor
  for(int i=0;i<K;i++){
    // polyvecl_pointwise_acc_montgomery_custom(&w1.vec[i],&mat[i],&z);
    // poly_invntt_custom_asm(&w1.vec[i]);
    // poly_mon2nor_custom(&w1.vec[i]);
    polyvecl_macc_invntt_mon2nor_asm(&w1.vec[i],&mat[i],&z);
  }
  polyveck_decompose_rvv(&w1, &w0, &w1);
  polyveck_pack_w1_custom(sig, &w1);

  size_t mu_sig_len=CRHBYTES+K*POLYW1_PACKEDBYTES;
  uint8_t* mu_sig=(uint8_t*)malloc(mu_sig_len*sizeof(uint8_t));
  memcpy(mu_sig,mu,CRHBYTES);
  memcpy(mu_sig+CRHBYTES,sig,K*POLYW1_PACKEDBYTES);
  shake256_custom(sig,SEEDBYTES,mu_sig,mu_sig_len);
  poly_challenge_custom(&cp, sig);
  //cp:mod_add+ntt
  // poly_caddq_custom(&cp);//modular add Q on cp's coefficients
  // poly_ntt_custom(&cp);
  poly_caddq_ntt_asm(&cp);

  /* Compute z, reject if it reveals secret */
  //z:pointwise multiplication+invntt+mon2nor+mod add+reduce
  for(int i=0;i<L;i++){
    // poly_pointwise_montgomery_custom(&z.vec[i],&cp,&s1.vec[i]);
    // poly_invntt_custom_asm(&z.vec[i]);
    // poly_mon2nor_custom(&z.vec[i]);
    poly_pointwisemul_invntt_mon2nor_asm(&z.vec[i],&cp,&s1.vec[i]);
    // poly_add_custom(&z.vec[i],&z.vec[i],&y.vec[i]);
    // poly_reduce_rvv(&z.vec[i]);
    poly_add_reduce_asm(&z.vec[i],&z.vec[i],&y.vec[i]);
  }
  if(polyvecl_chknorm_rvv(&z, GAMMA1 - BETA)) {
    goto rej;
  }

  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  //h:pointwise multiplication+invntt+mon2nor ++ w0:mod_add+sub+reduce
  for(int i=0;i<K;i++){
    //h process
    // poly_pointwise_montgomery_custom(&h.vec[i],&cp,&s2.vec[i]);
    // poly_invntt_custom_asm(&h.vec[i]);
    // poly_mon2nor_custom(&h.vec[i]);
    poly_pointwisemul_invntt_mon2nor_asm(&h.vec[i],&cp,&s2.vec[i]);
    //w0 process
    // poly_caddq_custom(&w0.vec[i]);
    // poly_sub_custom(&w0.vec[i],&w0.vec[i],&h.vec[i]);
    // poly_reduce_rvv(&w0.vec[i]);
    poly_caddq_modsub_reduce_asm(&w0.vec[i],&h.vec[i]);
  }
  if(polyveck_chknorm_rvv(&w0, GAMMA2 - BETA)) {
    goto rej;
  }

  /* Compute hints for w1 */
  //h:pointwise multiplication+invntt+mon2nor+reduce
  for(int i=0;i<K;i++){
    // poly_pointwise_montgomery_custom(&h.vec[i],&cp,&t0.vec[i]);
    // poly_invntt_custom_asm(&h.vec[i]);
    // poly_mon2nor_custom(&h.vec[i]);
    poly_pointwisemul_invntt_mon2nor_asm(&h.vec[i],&cp,&t0.vec[i]);
    poly_reduce_rvv(&h.vec[i]);
  }
  if(polyveck_chknorm_rvv(&h, GAMMA2)) {
    goto rej;
  }

  //w0:add+mod add q
  for(int i=0;i<K;i++){
    // poly_add_rvv(&w0.vec[i],&w0.vec[i],&h.vec[i]);
    // poly_caddq_custom(&w0.vec[i]);
    poly_add_caddq_asm(&w0.vec[i],&w0.vec[i],&h.vec[i]);
  }
  n = polyveck_make_hint_rvv(&h, &w0, &w1);
  if(n > OMEGA) {
    goto rej;
  }

  /* Write signature */
  pack_sig_custom(sig, sig, &z, &h);
  *siglen = CRYPTO_BYTES;

  free(tr_msg);
  free(mu_sig);
  return 0;

  return 0;
}

/*************************************************
* Name:        crypto_sign_signature_custom_asm_fortest
*
* Description: Computes signature.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_BYTES)
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign_signature_custom_asm_fortest(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const uint8_t *sk)
{
  unsigned int n;
  uint8_t seedbuf[2*SEEDBYTES + 3*CRHBYTES];
  uint8_t *rho, *tr, *key, *mu, *rhoprime;
  uint16_t nonce = 0;
  polyvecl mat[K], s1, y, z;
  polyveck t0, s2, w1, w0, h;
  poly cp;

  rho = seedbuf;
  tr = rho + SEEDBYTES;
  key = tr + CRHBYTES;
  mu = key + SEEDBYTES;
  rhoprime = mu + CRHBYTES;

  uint8_t* sk_tr=sk+(SEEDBYTES<<1);//pointer to packed byte start position for tr in sk

  memcpy(rho,sk,SEEDBYTES);//unpack for rho
  polyvec_matrix_expand_custom(mat, rho);
  polyvec_matrix_chorder(mat);//for consistency comparison, change mat from bit-reverse order to normal order
  //s1:unpack+mod_add+ntt
  uint8_t* sk_s1=sk_tr+CRHBYTES;
  for(int i=0;i<L;i++){
    // polyeta_unpack_custom(&s1.vec[i],sk_s1+i*POLYETA_PACKEDBYTES);
    // poly_caddq_custom(&s1.vec[i]);
    // poly_ntt_custom_asm(&s1.vec[i]);
    poly_eta_unpack_caddq_ntt_asm(&s1.vec[i],sk_s1+i*POLYETA_PACKEDBYTES);
  }
  //s2:unpack+mod_add+ntt
  uint8_t* sk_s2=sk_s1+L*POLYETA_PACKEDBYTES;
  for(int i=0;i<K;i++){
    // polyeta_unpack_custom(&s2.vec[i],sk_s2+i*POLYETA_PACKEDBYTES);
    // poly_caddq_custom(&s2.vec[i]);
    // poly_ntt_custom_asm(&s2.vec[i]);
    poly_eta_unpack_caddq_ntt_asm(&s2.vec[i],sk_s2+i*POLYETA_PACKEDBYTES);
  }
  //t0:unpack+mod_add+ntt
  uint8_t* sk_t0=sk_s2+K*POLYETA_PACKEDBYTES;
  for(int i=0;i<K;i++){
    // polyt0_unpack_custom(&t0.vec[i],sk_t0+i*POLYT0_PACKEDBYTES);
    // poly_caddq_custom(&t0.vec[i]);
    // poly_ntt_custom_asm(&t0.vec[i]);
    poly_t0_unpack_caddq_ntt_asm(&t0.vec[i],sk_t0+i*POLYT0_PACKEDBYTES);
  }

  /* Compute CRH(tr, msg) */
  size_t tr_msg_len=mlen+CRHBYTES;
  uint8_t* tr_msg=(uint8_t*)malloc(tr_msg_len*sizeof(uint8_t));
  memcpy(tr_msg,sk_tr,CRHBYTES);//skip unpack for tr here
  memcpy(tr_msg+CRHBYTES,m,mlen);
  shake256_custom(mu,CRHBYTES,tr_msg,tr_msg_len);

  memcpy(key,sk+SEEDBYTES,SEEDBYTES);//unpack for key

#ifdef DILITHIUM_RANDOMIZED_SIGNING
  randombytes(rhoprime, CRHBYTES);
#else
  crh_custom(rhoprime, key, SEEDBYTES + CRHBYTES);//consumes both key and mu
#endif

rej:
  /* Sample intermediate vector y */
  //y:Generation+mod_add+store into y+ntt
  for(int i=0;i<L;i++){
    // poly_uniform_gamma1_custom(&y.vec[i],rhoprime,L*nonce+i);
    // poly_caddq_custom(&y.vec[i]);
    // z.vec[i]=y.vec[i];
    // poly_ntt_custom_asm(&z.vec[i]);
    poly_gamma1_sample_caddq_ntt_asm(&y.vec[i],&z.vec[i],rhoprime,L*nonce+i);
  }
  nonce++;

  /* Matrix-vector multiplication */
  //Gradually multiply each row of matrix with z
  //w1:macc+invntt+mon2nor
  for(int i=0;i<K;i++){
    // polyvecl_pointwise_acc_montgomery_custom(&w1.vec[i],&mat[i],&z);
    // poly_invntt_custom_asm(&w1.vec[i]);
    // poly_mon2nor_custom(&w1.vec[i]);
    polyvecl_macc_invntt_mon2nor_asm(&w1.vec[i],&mat[i],&z);
  }
  polyveck_decompose_rvv(&w1, &w0, &w1);
  polyveck_pack_w1_custom(sig, &w1);

  size_t mu_sig_len=CRHBYTES+K*POLYW1_PACKEDBYTES;
  uint8_t* mu_sig=(uint8_t*)malloc(mu_sig_len*sizeof(uint8_t));
  memcpy(mu_sig,mu,CRHBYTES);
  memcpy(mu_sig+CRHBYTES,sig,K*POLYW1_PACKEDBYTES);
  shake256_custom(sig,SEEDBYTES,mu_sig,mu_sig_len);
  poly_challenge_custom(&cp, sig);
  //cp:mod_add+ntt
  // poly_caddq_custom(&cp);//modular add Q on cp's coefficients
  // poly_ntt_custom(&cp);
  poly_caddq_ntt_asm(&cp);

  /* Compute z, reject if it reveals secret */
  //z:pointwise multiplication+invntt+mon2nor+mod add+reduce
  for(int i=0;i<L;i++){
    // poly_pointwise_montgomery_custom(&z.vec[i],&cp,&s1.vec[i]);
    // poly_invntt_custom_asm(&z.vec[i]);
    // poly_mon2nor_custom(&z.vec[i]);
    poly_pointwisemul_invntt_mon2nor_asm(&z.vec[i],&cp,&s1.vec[i]);
    // poly_add_custom(&z.vec[i],&z.vec[i],&y.vec[i]);
    // poly_reduce_rvv(&z.vec[i]);
    poly_add_reduce_asm(&z.vec[i],&z.vec[i],&y.vec[i]);
  }
  if(polyvecl_chknorm_rvv(&z, GAMMA1 - BETA)) {
    goto rej;
  }

  /* Check that subtracting cs2 does not change high bits of w and low bits
   * do not reveal secret information */
  //h:pointwise multiplication+invntt+mon2nor ++ w0:mod_add+sub+reduce
  for(int i=0;i<K;i++){
    //h process
    // poly_pointwise_montgomery_custom(&h.vec[i],&cp,&s2.vec[i]);
    // poly_invntt_custom_asm(&h.vec[i]);
    // poly_mon2nor_custom(&h.vec[i]);
    poly_pointwisemul_invntt_mon2nor_asm(&h.vec[i],&cp,&s2.vec[i]);
    // poly_caddq_custom(&w0.vec[i]);
    // poly_sub_custom(&w0.vec[i],&w0.vec[i],&h.vec[i]);
    // poly_reduce_rvv(&w0.vec[i]);
    poly_caddq_modsub_reduce_asm(&w0.vec[i],&h.vec[i]);
  }
  if(polyveck_chknorm_rvv(&w0, GAMMA2 - BETA)) {
    goto rej;
  }

  /* Compute hints for w1 */
  //h:pointwise multiplication+invntt+mon2nor+reduce
  for(int i=0;i<K;i++){
    // poly_pointwise_montgomery_custom(&h.vec[i],&cp,&t0.vec[i]);
    // poly_invntt_custom_asm(&h.vec[i]);
    // poly_mon2nor_custom(&h.vec[i]);
    poly_pointwisemul_invntt_mon2nor_asm(&h.vec[i],&cp,&t0.vec[i]);
    poly_reduce_rvv(&h.vec[i]);
  }
  if(polyveck_chknorm_rvv(&h, GAMMA2)) {
    goto rej;
  }

  //w0:add+mod add q
  for(int i=0;i<K;i++){
    // poly_add_rvv(&w0.vec[i],&w0.vec[i],&h.vec[i]);
    // poly_caddq_custom(&w0.vec[i]);
    poly_add_caddq_asm(&w0.vec[i],&w0.vec[i],&h.vec[i]);
  }
  n = polyveck_make_hint_rvv(&h, &w0, &w1);
  if(n > OMEGA) {
    goto rej;
  }

  /* Write signature */
  pack_sig_custom(sig, sig, &z, &h);
  *siglen = CRYPTO_BYTES;

  free(tr_msg);
  free(mu_sig);
  return 0;

  return 0;
}

/*************************************************
* Name:        crypto_sign_custom_asm
*
* Description: Compute signed message.
*
* Arguments:   - uint8_t *sm: pointer to output signed message (allocated
*                             array with CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smlen: pointer to output length of signed
*                               message
*              - const uint8_t *m: pointer to message to be signed
*              - size_t mlen: length of message
*              - const uint8_t *sk: pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign_custom_asm(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk)
{
  size_t i;
  for(i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  crypto_sign_signature_custom_asm(sm, smlen, sm + CRYPTO_BYTES, mlen, sk);
  *smlen += mlen;
  return 0;
}

/*************************************************
* Name:        crypto_sign_custom_asm_fortest
*
* Description: Compute signed message.
*
* Arguments:   - uint8_t *sm: pointer to output signed message (allocated
*                             array with CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smlen: pointer to output length of signed
*                               message
*              - const uint8_t *m: pointer to message to be signed
*              - size_t mlen: length of message
*              - const uint8_t *sk: pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign_custom_asm_fortest(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const uint8_t *sk)
{
  size_t i;
  for(i = 0; i < mlen; ++i)
    sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
  crypto_sign_signature_custom_asm_fortest(sm, smlen, sm + CRYPTO_BYTES, mlen, sk);
  *smlen += mlen;
  return 0;
}

/*************************************************
* Name:        crypto_sign_verify_custom_asm
*
* Description: Verifies signature.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_verify_custom_asm(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk)
{
  unsigned int i;
  uint8_t buf[K*POLYW1_PACKEDBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[SEEDBYTES];
  uint8_t c2[SEEDBYTES];
  poly cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h;

  if(siglen != CRYPTO_BYTES) {
    return -1;
  }

  if(unpack_sig_custom(c, &z, &h, sig)) {
    return -1;
  }

  if(polyvecl_chknorm_rvv(&z, GAMMA1 - BETA)) {
    return -1;
  }
  //z:mod_add+ntt
  for(int i=0;i<L;i++){
    // poly_caddq_custom(&z.vec[i]);
    // poly_ntt_custom_asm(&z.vec[i]);
    poly_caddq_ntt_asm(&z.vec[i]);
  }

  poly_challenge_custom(&cp, c);
  //cp:mod_add+ntt
  // poly_caddq_custom(&cp);
  // poly_ntt_custom(&cp);
  poly_caddq_ntt_asm(&cp);
  //t1:unpack+shiftl+ntt+montgomery mul
  uint8_t* pk_t1=pk+SEEDBYTES;
  for(int i=0;i<K;i++){
    // polyt1_unpack_custom(&t1.vec[i],pk_t1+i*POLYT1_PACKEDBYTES);
    // poly_shiftl_rvv(&t1.vec[i]);
    // poly_ntt_custom_asm(&t1.vec[i]);
    poly_unpackt1_shiftl_ntt_asm(&t1.vec[i],pk_t1+i*POLYT1_PACKEDBYTES);
    poly_pointwise_montgomery_custom(&t1.vec[i],&cp,&t1.vec[i]);
  }

  memcpy(rho,pk,SEEDBYTES);//unpack rho from pk
  //Gradually generate each row of matrix and multiply the row with z
  for(int i=0;i<K;i++){
    for(int j=0;j<L;j++){
      poly_uniform_custom(&mat[i].vec[j],rho,(i << 8) + j);
    }
    //w1:macc+sub+invntt+mon2nor
    // polyvecl_pointwise_acc_montgomery_custom(&w1.vec[i],&mat[i],&z);
    // poly_sub_custom(&w1.vec[i],&w1.vec[i],&t1.vec[i]);
    // poly_invntt_custom_asm(&w1.vec[i]);
    // poly_mon2nor_custom(&w1.vec[i]);
    poly_macc_sub_invntt_mon2nor_asm(&w1.vec[i],&mat[i],&z,&t1.vec[i]);
  }

  /* Reconstruct w1 */
  //polyveck_caddq(&w1);//no longer need caddq
  polyveck_use_hint_rvv(&w1, &w1, &h);
  polyveck_pack_w1_custom(buf, &w1);

  /* Compute CRH(CRH(rho, t1), msg) */
  crh_custom(mu, pk, CRYPTO_PUBLICKEYBYTES);
  size_t mu_m_len=CRHBYTES+mlen;
  uint8_t* mu_m=(uint8_t*)malloc(mu_m_len*sizeof(uint8_t));
  memcpy(mu_m,mu,CRHBYTES);
  memcpy(mu_m+CRHBYTES,m,mlen);
  shake256_custom(mu,CRHBYTES,mu_m,mu_m_len);

  /* Call random oracle and verify challenge */
  size_t mu_buf_len=CRHBYTES+K*POLYW1_PACKEDBYTES;
  uint8_t* mu_buf=(uint8_t*)malloc(mu_buf_len*sizeof(uint8_t));
  memcpy(mu_buf,mu,CRHBYTES);
  memcpy(mu_buf+CRHBYTES,buf, K*POLYW1_PACKEDBYTES);
  shake256_custom(c2,SEEDBYTES,mu_buf,mu_buf_len);

  for(i = 0; i < SEEDBYTES; ++i)
    if(c[i] != c2[i]) {
      return -1;
    }

  free(mu_m);
  free(mu_buf);

  return 0;
}

/*************************************************
* Name:        crypto_sign_verify_custom_asm_fortest
*
* Description: Verifies signature.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_verify_custom_asm_fortest(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const uint8_t *pk)
{
  unsigned int i;
  uint8_t buf[K*POLYW1_PACKEDBYTES];
  uint8_t rho[SEEDBYTES];
  uint8_t mu[CRHBYTES];
  uint8_t c[SEEDBYTES];
  uint8_t c2[SEEDBYTES];
  poly cp;
  polyvecl mat[K], z;
  polyveck t1, w1, h;

  if(siglen != CRYPTO_BYTES) {
    return -1;
  }

  if(unpack_sig_custom(c, &z, &h, sig)) {
    return -1;
  }

  if(polyvecl_chknorm_rvv(&z, GAMMA1 - BETA)) {
    return -1;
  }
  //z:mod_add+ntt
  for(int i=0;i<L;i++){
    // poly_caddq_custom(&z.vec[i]);
    // poly_ntt_custom_asm(&z.vec[i]);
    poly_caddq_ntt_asm(&z.vec[i]);
  }

  poly_challenge_custom(&cp, c);
  // poly_caddq_custom(&cp);
  // poly_ntt_custom(&cp);
  poly_caddq_ntt_asm(&cp);
  //t1:unpack+shiftl+ntt+montgomery mul
  uint8_t* pk_t1=pk+SEEDBYTES;
  for(int i=0;i<K;i++){
    // polyt1_unpack_custom(&t1.vec[i],pk_t1+i*POLYT1_PACKEDBYTES);
    // poly_shiftl_rvv(&t1.vec[i]);
    // poly_ntt_custom_asm(&t1.vec[i]);
    poly_unpackt1_shiftl_ntt_asm(&t1.vec[i],pk_t1+i*POLYT1_PACKEDBYTES);
    poly_pointwise_montgomery_custom(&t1.vec[i],&cp,&t1.vec[i]);
  }

  memcpy(rho,pk,SEEDBYTES);//unpack rho from pk
  //Gradually generate each row of matrix and multiply the row with z
  for(int i=0;i<K;i++){
    for(int j=0;j<L;j++){
      poly_uniform_custom(&mat[i].vec[j],rho,(i << 8) + j);
      poly_chorder(&mat[i].vec[j]);//for consistency comparison, change mat from bit-reverse order to normal order
    }
    //w1:macc+sub+invntt+mon2nor
    // polyvecl_pointwise_acc_montgomery_custom(&w1.vec[i],&mat[i],&z);
    // poly_sub_custom(&w1.vec[i],&w1.vec[i],&t1.vec[i]);
    // poly_invntt_custom_asm(&w1.vec[i]);
    // poly_mon2nor_custom(&w1.vec[i]);
    poly_macc_sub_invntt_mon2nor_asm(&w1.vec[i],&mat[i],&z,&t1.vec[i]);
  }

  /* Reconstruct w1 */
  //polyveck_caddq(&w1);//no longer need caddq
  polyveck_use_hint_rvv(&w1, &w1, &h);
  polyveck_pack_w1_custom(buf, &w1);

  /* Compute CRH(CRH(rho, t1), msg) */
  crh_custom(mu, pk, CRYPTO_PUBLICKEYBYTES);
  size_t mu_m_len=CRHBYTES+mlen;
  uint8_t* mu_m=(uint8_t*)malloc(mu_m_len*sizeof(uint8_t));
  memcpy(mu_m,mu,CRHBYTES);
  memcpy(mu_m+CRHBYTES,m,mlen);
  shake256_custom(mu,CRHBYTES,mu_m,mu_m_len);

  /* Call random oracle and verify challenge */
  size_t mu_buf_len=CRHBYTES+K*POLYW1_PACKEDBYTES;
  uint8_t* mu_buf=(uint8_t*)malloc(mu_buf_len*sizeof(uint8_t));
  memcpy(mu_buf,mu,CRHBYTES);
  memcpy(mu_buf+CRHBYTES,buf, K*POLYW1_PACKEDBYTES);
  shake256_custom(c2,SEEDBYTES,mu_buf,mu_buf_len);

  for(i = 0; i < SEEDBYTES; ++i)
    if(c[i] != c2[i]) {
      return -1;
    }

  free(mu_m);
  free(mu_buf);

  return 0;
}

/*************************************************
* Name:        crypto_sign_open_custom_asm
*
* Description: Verify signed message.
*
* Arguments:   - uint8_t *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mlen: pointer to output length of message
*              - const uint8_t *sm: pointer to signed message
*              - size_t smlen: length of signed message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_open_custom_asm(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *pk)
{
  size_t i;

  if(smlen < CRYPTO_BYTES) {
    goto badsig;
  }

  *mlen = smlen - CRYPTO_BYTES;
  if(crypto_sign_verify_custom_asm(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk)) {
    goto badsig;
  }
  else {
    /* All good, copy msg, return 0 */
    for(i = 0; i < *mlen; ++i)
      m[i] = sm[CRYPTO_BYTES + i];
    return 0;
  }

badsig:
  /* Signature verification failed */
  *mlen = -1;
  for(i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}

/*************************************************
* Name:        crypto_sign_open_custom_asm_fortest
*
* Description: Verify signed message.
*
* Arguments:   - uint8_t *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mlen: pointer to output length of message
*              - const uint8_t *sm: pointer to signed message
*              - size_t smlen: length of signed message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_open_custom_asm_fortest(uint8_t *m,
                     size_t *mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const uint8_t *pk)
{
  size_t i;

  if(smlen < CRYPTO_BYTES) {
    goto badsig;
  }

  *mlen = smlen - CRYPTO_BYTES;
  if(crypto_sign_verify_custom_asm_fortest(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, *mlen, pk)) {
    goto badsig;
  }
  else {
    /* All good, copy msg, return 0 */
    for(i = 0; i < *mlen; ++i)
      m[i] = sm[CRYPTO_BYTES + i];
    return 0;
  }

badsig:
  /* Signature verification failed */
  *mlen = -1;
  for(i = 0; i < smlen; ++i)
    m[i] = 0;

  return -1;
}