#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include "fips202.h"
#include "params.h"
#include "symmetric.h"
#include "indcpa.h"
#include "poly.h"
#include "indcpa.h"
#include "reduce.h"
#include "poly_reorg.h"

/*************************************************
* Name:        pack_pk_custom
*
* Description: Serialize the public key as concatenation of the
*              serialized vector of polynomials pk
*              and the public seed used to generate the matrix A.
*
* Arguments:   uint8_t *r:          pointer to the output serialized public key
*              polyvec *pk:         pointer to the input public-key polyvec
*              const uint8_t *seed: pointer to the input public seed
**************************************************/
static void pack_pk_custom(uint8_t r[KYBER_INDCPA_PUBLICKEYBYTES],
                    polyvec *pk,
                    const uint8_t seed[KYBER_SYMBYTES])
{
  size_t i;
  polyvec_tobytes_custom(r, pk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    r[i+KYBER_POLYVECBYTES] = seed[i];
}

/*************************************************
* Name:        unpack_pk_custom
*
* Description: De-serialize public key from a byte array;
*              approximate inverse of pack_pk
*
* Arguments:   - polyvec *pk:             pointer to output public-key
*                                         polynomial vector
*              - uint8_t *seed:           pointer to output seed to generate
*                                         matrix A
*              - const uint8_t *packedpk: pointer to input serialized public key
**************************************************/
static void unpack_pk_custom(polyvec *pk,
                      uint8_t seed[KYBER_SYMBYTES],
                      const uint8_t packedpk[KYBER_INDCPA_PUBLICKEYBYTES])
{
  size_t i;
  polyvec_frombytes_custom(pk, packedpk);
  for(i=0;i<KYBER_SYMBYTES;i++)
    seed[i] = packedpk[i+KYBER_POLYVECBYTES];
}

/*************************************************
* Name:        pack_sk_custom
*
* Description: Serialize the secret key
*
* Arguments:   - uint8_t *r:  pointer to output serialized secret key
*              - polyvec *sk: pointer to input vector of polynomials (secret key)
**************************************************/
static void      pack_sk_custom(uint8_t r[KYBER_INDCPA_SECRETKEYBYTES], polyvec *sk)
{
  polyvec_tobytes_custom(r, sk);
}

/*************************************************
* Name:        unpack_sk_custom
*
* Description: De-serialize the secret key;
*              inverse of pack_sk
*
* Arguments:   - polyvec *sk:             pointer to output vector of
*                                         polynomials (secret key)
*              - const uint8_t *packedsk: pointer to input serialized secret key
**************************************************/
static void unpack_sk_custom(polyvec *sk,
                      const uint8_t packedsk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec_frombytes_custom(sk, packedsk);
}

/*************************************************
* Name:        pack_ciphertext_custom
*
* Description: Serialize the ciphertext as concatenation of the
*              compressed and serialized vector of polynomials b
*              and the compressed and serialized polynomial v
*
* Arguments:   uint8_t *r: pointer to the output serialized ciphertext
*              poly *pk:   pointer to the input vector of polynomials b
*              poly *v:    pointer to the input polynomial v
**************************************************/
static void pack_ciphertext_custom(uint8_t r[KYBER_INDCPA_BYTES],
                            polyvec *b,
                            poly *v)
{
  polyvec_compress_custom(r, b);
  poly_compress_custom(r+KYBER_POLYVECCOMPRESSEDBYTES, v);
}

/*************************************************
* Name:        unpack_ciphertext_custom
*
* Description: De-serialize and decompress ciphertext from a byte array;
*              approximate inverse of pack_ciphertext
*
* Arguments:   - polyvec *b:       pointer to the output vector of polynomials b
*              - poly *v:          pointer to the output polynomial v
*              - const uint8_t *c: pointer to the input serialized ciphertext
**************************************************/
static void unpack_ciphertext_custom(polyvec *b,
                              poly *v,
                              const uint8_t c[KYBER_INDCPA_BYTES])
{
  polyvec_decompress_custom(b, c);
  poly_decompress_custom(v, c+KYBER_POLYVECCOMPRESSEDBYTES);
}

/*************************************************
* Name:        rej_uniform_custom
*
* Description: Run rejection sampling on uniform random bytes to generate
*              uniform random integers mod q
*
* Arguments:   - int16_t *r:          pointer to output buffer
*              - unsigned int len:    requested number of 16-bit integers
*                                     (uniform mod q)
*              - const uint8_t *buf:  pointer to input buffer
*                                     (assumed to be uniform random bytes)
*              - unsigned int buflen: length of input buffer in bytes
*
* Returns number of sampled 16-bit integers (at most len)
**************************************************/
unsigned int rej_uniform_custom(int16_t *r,
                                       unsigned int len,
                                       const uint8_t *buf,
                                       unsigned int buflen)
{
  size_t ctr, pos, vl, valid_num;
  vuint8m1_t vreg_m;
  vint16m1_t vreg_rej;

  csr_validnum_rw();
  ctr = pos = 0;

  while((ctr < len) && (pos < buflen)) {

    if((pos + UNPACK_REJ_LOAD_BYTE_KYBER) < buflen) {
      vl = vsetvl_e8m1(UNPACK_REJ_LOAD_BYTE_KYBER);
    }
    else {
      vl = vsetvl_e8m1(buflen - pos);
    }
    
    vreg_m = vle8_v_u8m1(buf + pos, vl);
    pos += vl;

    // in such case, vl is always divisible by 3
    // because vl is not used, user should use vsetvl_e16m1_wrapper instead of vsetvl_e16m1 
    // to avoid this vector configuration being omitted 
    vsetvl_e16m1_wrapper((vl / 3) << 1);
    vreg_rej = unpack_vx_i16m1(vreg_m, 12);
    vreg_rej = sample_rej_vx_i16m1(vreg_rej, KYBER_Q);
    valid_num = csr_validnum_rw();

    if(ctr + valid_num > len) {
      valid_num = len - ctr;
    }

    vl = vsetvl_e16m1(valid_num);
    vse16_v_i16m1(r + ctr, vreg_rej, valid_num);
    ctr += valid_num;
  }
  
  return ctr;
}

#define gen_a_custom(A,B)  gen_matrix_custom(A,B,0)
#define gen_at_custom(A,B) gen_matrix_custom(A,B,1)

/*************************************************
* Name:        gen_matrix_custom
*
* Description: Deterministically generate matrix A (or the transpose of A)
*              from a seed. Entries of the matrix are polynomials that look
*              uniformly random. Performs rejection sampling on output of
*              a XOF
*
* Arguments:   - polyvec *a:          pointer to ouptput matrix A
*              - const uint8_t *seed: pointer to input seed
*              - int transposed:      boolean deciding whether A or A^T
*                                     is generated
**************************************************/
#define GEN_MATRIX_NBLOCKS ((12*KYBER_N/8*(1 << 12)/KYBER_Q \
                             + XOF_BLOCKBYTES)/XOF_BLOCKBYTES)
// Not static for benchmarking
// gen_matrix_custom goes wrong when the optimize level set to O3, 
// because a vse16 inst executes before vset inst from a reason I dont know
#pragma GCC push_options
#pragma GCC optimize("O2")
void gen_matrix_custom(polyvec *a, const uint8_t seed[KYBER_SYMBYTES], int transposed)
{
  unsigned int ctr, i, j;
  unsigned int buflen;
  uint8_t buf[GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES+2];

  for(i=0;i<KYBER_K;i++) {
    for(j=0;j<KYBER_K;j++) {
      if(transposed) {
        xof_absorb_custom(seed, i, j);
      }
      else {
        xof_absorb_custom(seed, j, i);
      }

      xof_squeezeblocks_custom(buf, GEN_MATRIX_NBLOCKS, false);
      buflen = GEN_MATRIX_NBLOCKS*XOF_BLOCKBYTES;
      ctr = rej_uniform_custom(a[i].vec[j].coeffs, KYBER_N, buf, buflen);

      while(ctr < KYBER_N) {
        xof_squeezeblocks_custom(buf, 1, true);
        buflen = XOF_BLOCKBYTES;
        ctr += rej_uniform_custom(a[i].vec[j].coeffs + ctr, KYBER_N - ctr, buf, buflen);
      }
    }
  }
}

/*************************************************
* Name:        indcpa_keypair_custom_for_test
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*              Specially considered(change the coefficient permutation after ntt) 
*              so that the function result is the same as indcpa_keypair               
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
                              (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
**************************************************/
void indcpa_keypair_custom_for_test(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;
  polyvec a[KYBER_K], e, pkpv, skpv;

  csr_modulusq_rw(KYBER_Q);
  csr_qinv_rw(QINV_HW);

  // randombytes(buf, KYBER_SYMBYTES);
  for(int i = 0; i < KYBER_SYMBYTES; i++) {
    buf[i] = i;
  }
  hash_g_custom(buf, buf, KYBER_SYMBYTES);

  gen_a_custom(a, publicseed);

  for(i=0; i < KYBER_K; i++)
    poly_getnoise_eta1_custom(&skpv.vec[i], noiseseed, nonce++);
  for(i=0; i < KYBER_K; i++)
    poly_getnoise_eta1_custom(&e.vec[i], noiseseed, nonce++);

  /*
  * get non-negetive coefficients
  */
  polyvec_mod_add_q(&skpv);
  polyvec_mod_add_q(&e);

  polyvec_ntt_custom(&skpv);
  polyvec_ntt_custom(&e);

  /*
  * transfer the coefficient permutation in polynomils of a 
  * from bit-reversed order(the same as in indcpa_keypair) to normal order
  */ 
  for(i = 0; i < KYBER_K; i++) {
    // the coefficient permutation order of a is considered as bitreversed order by default here
    polyvec_bitreverse_to_standard_all(&a[i]);
  } 

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    polyvec_pointwise_acc_montgomery_custom(&pkpv.vec[i], &a[i], &skpv);
    poly_tomont_custom(&pkpv.vec[i]);
  }
  
  /*
  * transfer the coefficient permutation from normal order to bitreversed order (the same as in indcpa_keypair)
  */ 
  polyvec_standard_to_bitreverse_all(&skpv);
  polyvec_standard_to_bitreverse_all(&pkpv);
  polyvec_standard_to_bitreverse_all(&e);

  polyvec_add_custom(&pkpv, &pkpv, &e);

  pack_sk_custom(sk, &skpv);
  pack_pk_custom(pk, &pkpv, publicseed);
}

/*************************************************
* Name:        indcpa_enc_custom_for_test
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *c:           pointer to output ciphertext
*                                      (of length KYBER_INDCPA_BYTES bytes)
*              - const uint8_t *m:     pointer to input message
*                                      (of length KYBER_INDCPA_MSGBYTES bytes)
*              - const uint8_t *pk:    pointer to input public key
*                                      (of length KYBER_INDCPA_PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins
*                                      used as seed (of length KYBER_SYMBYTES)
*                                      to deterministically generate all
*                                      randomness
**************************************************/
void indcpa_enc_custom_for_test(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[KYBER_K], bp;
  poly v, k, epp;

  csr_modulusq_rw(KYBER_Q);
  csr_qinv_rw(QINV_HW);

  unpack_pk_custom(&pkpv, seed, pk);  // the coefficient permutation order of pkpv is bitreversed order here
  poly_frommsg_custom(&k, m);  // the coefficient permutation order of m is normal order here
  gen_at_custom(at, seed);  // the coefficient permutation order of at is considered as bitreversed order by default here

  /*
  * transfer the coefficient permutation in polynomils of pkpv and at
  * from bit-reversed order(the same as in indcpa_enc) to normal order
  */ 
  for(i = 0; i < KYBER_K; i++) {
    polyvec_bitreverse_to_standard_all(&at[i]);
  }
  polyvec_bitreverse_to_standard_all(&pkpv);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1_custom(sp.vec+i, coins, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta2_custom(ep.vec+i, coins, nonce++);
  poly_getnoise_eta2_custom(&epp, coins, nonce++);

  /*
  * get non-negetive coefficients
  */
  polyvec_mod_add_q(&sp);
  polyvec_mod_add_q(&ep);
  poly_mod_add_q(&epp);

  polyvec_ntt_custom(&sp);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_pointwise_acc_montgomery_custom(&bp.vec[i], &at[i], &sp);

  polyvec_pointwise_acc_montgomery_custom(&v, &pkpv, &sp);

  polyvec_invntt_tomont_custom(&bp);
  poly_invntt_tomont_custom(&v);

  polyvec_add_custom(&bp, &bp, &ep);
  poly_add_custom(&v, &v, &epp);
  poly_add_custom(&v, &v, &k);

  pack_ciphertext_custom(c, &bp, &v);// coefficients in v, bp are normal order here
}

/*************************************************
* Name:        indcpa_dec_custom_for_test
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *m:        pointer to output decrypted message
*                                   (of length KYBER_INDCPA_MSGBYTES)
*              - const uint8_t *c:  pointer to input ciphertext
*                                   (of length KYBER_INDCPA_BYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length KYBER_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec_custom_for_test(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec bp, skpv;
  poly v, mp;

  csr_modulusq_rw(KYBER_Q);
  csr_qinv_rw(QINV_HW);

  unpack_ciphertext_custom(&bp, &v, c);// coefficients in bp are normal order here
  unpack_sk_custom(&skpv, sk);  // the coefficient permutation order of skpv is bitreversed order here

  /*
  * transfer the coefficient permutation in polynomils of skpv
  * from bit-reversed order(the same as in indcpa_dec) to normal order
  */ 
  polyvec_bitreverse_to_standard_all(&skpv);

  polyvec_ntt_custom(&bp);
  polyvec_pointwise_acc_montgomery_custom(&mp, &skpv, &bp);
  poly_invntt_tomont_custom(&mp);

  poly_sub_custom(&mp, &v, &mp);

  poly_tomsg_custom(m, &mp);
}


/*************************************************
* Name:        indcpa_keypair_custom
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
                              (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
**************************************************/
void indcpa_keypair_custom(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                    uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;
  polyvec a[KYBER_K], e, pkpv, skpv;

  csr_modulusq_rw(KYBER_Q);
  csr_qinv_rw(QINV_HW);

  // randombytes(buf, KYBER_SYMBYTES);
  for(int i = 0; i < KYBER_SYMBYTES; i++) {
    buf[i] = i;
  }
  hash_g_custom(buf, buf, KYBER_SYMBYTES);

  gen_a_custom(a, publicseed);

  for(i=0; i < KYBER_K; i++)
    poly_getnoise_eta1_custom(&skpv.vec[i], noiseseed, nonce++);
  for(i=0; i < KYBER_K; i++)
    poly_getnoise_eta1_custom(&e.vec[i], noiseseed, nonce++);

  /*
  * get non-negetive coefficients
  */
  polyvec_mod_add_q(&skpv);
  polyvec_mod_add_q(&e);

  polyvec_ntt_custom(&skpv);
  polyvec_ntt_custom(&e);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    // the coefficient permutation order of a is considered as normal order by default here
    polyvec_pointwise_acc_montgomery_custom(&pkpv.vec[i], &a[i], &skpv);
    poly_tomont_custom(&pkpv.vec[i]);
  }

  polyvec_add_custom(&pkpv, &pkpv, &e);

  pack_sk_custom(sk, &skpv);
  pack_pk_custom(pk, &pkpv, publicseed);
}

/*************************************************
* Name:        indcpa_enc_custom
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *c:           pointer to output ciphertext
*                                      (of length KYBER_INDCPA_BYTES bytes)
*              - const uint8_t *m:     pointer to input message
*                                      (of length KYBER_INDCPA_MSGBYTES bytes)
*              - const uint8_t *pk:    pointer to input public key
*                                      (of length KYBER_INDCPA_PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins
*                                      used as seed (of length KYBER_SYMBYTES)
*                                      to deterministically generate all
*                                      randomness
**************************************************/
void indcpa_enc_custom(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[KYBER_K], bp;
  poly v, k, epp;

  csr_modulusq_rw(KYBER_Q);
  csr_qinv_rw(QINV_HW);

  unpack_pk_custom(&pkpv, seed, pk);  // the coefficient permutation order of pkpv is normal order here
  poly_frommsg_custom(&k, m);  // the coefficient permutation order of m is normal order here
  gen_at_custom(at, seed);  // the coefficient permutation order of at is considered as normal order by default here

  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta1_custom(sp.vec+i, coins, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise_eta2_custom(ep.vec+i, coins, nonce++);
  poly_getnoise_eta2_custom(&epp, coins, nonce++);

  /*
  * get non-negetive coefficients
  */
  polyvec_mod_add_q(&sp);
  polyvec_mod_add_q(&ep);
  poly_mod_add_q(&epp);

  polyvec_ntt_custom(&sp);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++)
    polyvec_pointwise_acc_montgomery_custom(&bp.vec[i], &at[i], &sp);

  polyvec_pointwise_acc_montgomery_custom(&v, &pkpv, &sp);

  polyvec_invntt_tomont_custom(&bp);
  poly_invntt_tomont_custom(&v);

  polyvec_add_custom(&bp, &bp, &ep);
  poly_add_custom(&v, &v, &epp);
  poly_add_custom(&v, &v, &k);

  pack_ciphertext_custom(c, &bp, &v);
}


/*************************************************
* Name:        indcpa_dec_custom
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *m:        pointer to output decrypted message
*                                   (of length KYBER_INDCPA_MSGBYTES)
*              - const uint8_t *c:  pointer to input ciphertext
*                                   (of length KYBER_INDCPA_BYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length KYBER_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec_custom(uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec bp, skpv;
  poly v, mp;

  csr_modulusq_rw(KYBER_Q);
  csr_qinv_rw(QINV_HW);

  unpack_ciphertext_custom(&bp, &v, c);
  unpack_sk_custom(&skpv, sk);
  
  polyvec_ntt_custom(&bp);
  polyvec_pointwise_acc_montgomery_custom(&mp, &skpv, &bp);
  poly_invntt_tomont_custom(&mp);
  
  poly_sub_custom(&mp, &v, &mp);

  poly_tomsg_custom(m, &mp);
}

/*************************************************
* Name:        indcpa_keypair_custom_asm
*
* Description: Generates public and private key for the CPA-secure
*              public-key encryption scheme underlying Kyber
*
* Arguments:   - uint8_t *pk: pointer to output public key
*                             (of length KYBER_INDCPA_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key
                              (of length KYBER_INDCPA_SECRETKEYBYTES bytes)
**************************************************/
void indcpa_keypair_custom_asm(uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                              uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  unsigned int i;
  uint8_t buf[2*KYBER_SYMBYTES];
  const uint8_t *publicseed = buf;
  const uint8_t *noiseseed = buf+KYBER_SYMBYTES;
  uint8_t nonce = 0;
  polyvec a[KYBER_K], e, pkpv, skpv;

  csr_modulusq_rw(KYBER_Q);
  csr_qinv_rw(QINV_HW);

  // randombytes(buf, KYBER_SYMBYTES);
  for(int i = 0; i < KYBER_SYMBYTES; i++) {
    buf[i] = i;
  }
  hash_g_custom(buf, buf, KYBER_SYMBYTES);

  /*
  ** skpv: polynomial generation + mod_add + ntt + pack
  */
  // skpv: polynomial generation + mod_add + ntt
  for(i=0; i < KYBER_K; i++)
    poly_eta1_add_ntt_asm(&skpv.vec[i], noiseseed, nonce++);
  // skpv: polynomial pack  
  pack_sk_custom(sk, &skpv);

  gen_a_custom(a, publicseed);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    // the coefficient permutation order of a is considered as normal order by default here
    polyvec_base_mul_acc_tomont_custom_asm(&pkpv.vec[i], &a[i], &skpv);
  }

  /*
  ** e: polynomial generation + mod_add + ntt
  */
  for(i=0; i < KYBER_K; i++)
    poly_eta1_add_ntt_asm(&e.vec[i], noiseseed, nonce++);

  polyvec_add_custom(&pkpv, &pkpv, &e);

  pack_pk_custom(pk, &pkpv, publicseed);
}

/*************************************************
* Name:        indcpa_enc_custom_asm
*
* Description: Encryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *c:           pointer to output ciphertext
*                                      (of length KYBER_INDCPA_BYTES bytes)
*              - const uint8_t *m:     pointer to input message
*                                      (of length KYBER_INDCPA_MSGBYTES bytes)
*              - const uint8_t *pk:    pointer to input public key
*                                      (of length KYBER_INDCPA_PUBLICKEYBYTES)
*              - const uint8_t *coins: pointer to input random coins
*                                      used as seed (of length KYBER_SYMBYTES)
*                                      to deterministically generate all
*                                      randomness
**************************************************/
void indcpa_enc_custom_asm(uint8_t c[KYBER_INDCPA_BYTES],
                const uint8_t m[KYBER_INDCPA_MSGBYTES],
                const uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES],
                const uint8_t coins[KYBER_SYMBYTES])
{
  unsigned int i;
  uint8_t seed[KYBER_SYMBYTES];
  uint8_t nonce = 0;
  polyvec sp, pkpv, ep, at[KYBER_K], bp;
  poly v, k, epp;

  csr_modulusq_rw(KYBER_Q);
  csr_qinv_rw(QINV_HW);

  unpack_pk_custom(&pkpv, seed, pk);  // the coefficient permutation order of pkpv is normal order here

  // for(i=0;i<KYBER_K;i++)
  //   poly_getnoise_eta1_custom(sp.vec+i, coins, nonce++);
  // polyvec_mod_add_q(&sp);
  // polyvec_ntt_custom(&sp);
  for(i=0;i<KYBER_K;i++)
    poly_eta1_add_ntt_asm(&sp.vec[i], coins, nonce++);
  
  gen_at_custom(at, seed);  // the coefficient permutation order of at is considered as normal order by default here
  
  // matrix-vector multiplication
  // for(i=0;i<KYBER_K;i++)
  //   polyvec_pointwise_acc_montgomery_custom(&bp.vec[i], &at[i], &sp);
  // polyvec_invntt_tomont_custom(&bp);
  for(i=0;i<KYBER_K;i++)
    polyvec_base_mul_acc_intt_tomont_custom_asm(&bp.vec[i], &at[i], &sp);

  // for(i=0;i<KYBER_K;i++)
  //   poly_getnoise_eta2_custom(ep.vec+i, coins, nonce++);
  // polyvec_mod_add_q(&ep);
  for(i=0;i<KYBER_K;i++) {
    // poly_eta2_add_asm(ep.vec+i, coins, nonce++);
    poly_eta2_addq_add_asm(bp.vec+i, coins, nonce++);
  }

  // poly_getnoise_eta2_custom(&epp, coins, nonce++);
  // poly_mod_add_q(&epp);
  poly_eta2_add_asm(&epp, coins, nonce++);

  // polyvec_add_custom(&bp, &bp, &ep);
  
  // polyvec_pointwise_acc_montgomery_custom(&v, &pkpv, &sp);
  // poly_invntt_tomont_custom(&v);
  polyvec_base_mul_acc_intt_tomont_custom_asm(&v, &pkpv, &sp);

  poly_add_custom(&v, &v, &epp);

  // poly_frommsg_custom(&k, m);  // the coefficient permutation order of m is normal order here
  // poly_add_custom(&v, &v, &k);
  poly_frommsg_add_custom(&v, m);

  pack_ciphertext_custom(c, &bp, &v);
}

/*************************************************
* Name:        indcpa_dec_custom_asm
*
* Description: Decryption function of the CPA-secure
*              public-key encryption scheme underlying Kyber.
*
* Arguments:   - uint8_t *m:        pointer to output decrypted message
*                                   (of length KYBER_INDCPA_MSGBYTES)
*              - const uint8_t *c:  pointer to input ciphertext
*                                   (of length KYBER_INDCPA_BYTES)
*              - const uint8_t *sk: pointer to input secret key
*                                   (of length KYBER_INDCPA_SECRETKEYBYTES)
**************************************************/
void indcpa_dec_custom_asm(uint8_t m[KYBER_INDCPA_MSGBYTES],
                          const uint8_t c[KYBER_INDCPA_BYTES],
                          const uint8_t sk[KYBER_INDCPA_SECRETKEYBYTES])
{
  polyvec bp, skpv;
  poly v, mp;

  csr_modulusq_rw(KYBER_Q);
  csr_qinv_rw(QINV_HW);

  unpack_ciphertext_custom(&bp, &v, c);
  unpack_sk_custom(&skpv, sk);

  // polyvec_ntt_custom(&bp);
  polyvec_ntt_custom_asm(&bp);
  // polyvec_pointwise_acc_montgomery_custom(&mp, &skpv, &bp);
  // poly_invntt_tomont_custom(&mp);
  polyvec_base_mul_acc_intt_tomont_custom_asm(&mp, &skpv, &bp);
  
  // poly_sub_custom(&mp, &v, &mp);
  // poly_tomsg_custom(m, &mp);
  poly_sub_tomsg_custom(m, &v, &mp);
}