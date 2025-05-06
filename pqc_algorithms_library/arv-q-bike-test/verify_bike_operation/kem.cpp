/******************************************************************************
 * BIKE -- Bit Flipping Key Encapsulation
 *
 * Copyright (c) 2021 Nir Drucker, Shay Gueron, Rafael Misoczki, Tobias Oder,
 * Tim Gueneysu, Jan Richter-Brockmann.
 * Contact: drucker.nir@gmail.com, shay.gueron@gmail.com,
 * rafaelmisoczki@google.com, tobias.oder@rub.de, tim.gueneysu@rub.de,
 * jan.richter-brockmann@rub.de.
 *
 * Permission to use this code for BIKE is granted.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * The names of the contributors may not be used to endorse or promote
 *   products derived from this software without specific prior written
 *   permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ""AS IS"" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS CORPORATION OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "hash_wrapper.h"
#include "decode.h"
#include "sampling.h"
#include "kem.h"
#include "conversions.h"
#include "shake_prng.h"
#include "poly_op_util.h"
#include "sampling.h"
#include "gf2x_port.h"
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"

 //Comparing value in a constant time manner.
_INLINE_ uint32_t safe_cmp(IN const uint8_t* a,
    IN const uint8_t* b,
    IN const uint32_t size)
{
    volatile uint8_t res = 0;

    for (uint32_t i = 0; i < size; ++i)
    {
        res |= (a[i] ^ b[i]);
    }

    return (res == 0);
}

// Function H. It uses the extract-then-expand paradigm based on SHA384 and
// AES256-CTR PRNG to produce e from m.
_INLINE_ status_t functionH(
        OUT uint8_t * e,
        IN const uint8_t * m)
{
    status_t res = SUCCESS;

    // format seed as a 32-bytes input:
    seed_t seed_for_hash;
    memcpy(seed_for_hash.raw, m, ELL_SIZE);

    // use the seed to generate sparse error vector e:
    DMSG("    Generating random error.\n");
    shake256_prng_state_t prng_state = {0};
    shake256_init(seed_for_hash.raw, ELL_SIZE, &prng_state);
    res = generate_sparse_rep_keccak(e, T1, N_BITS, &prng_state); CHECK_STATUS(res);

    EXIT:
    DMSG("  Exit functionH.\n");
    return res;
}

// Function L. Computes L(e0 || e1)
_INLINE_ status_t functionL(
        OUT uint8_t * output,
        IN const uint8_t * e)
{
    status_t res = SUCCESS;
    uint8_t hash_value[SHA384_HASH_SIZE] = {0};
    uint8_t e_split[2 * R_SIZE] = {0};

    poly_split_polynomial(e_split, &e_split[R_SIZE], e);

    // select hash function
    sha3_384(hash_value, e_split, 2*R_SIZE);

    memcpy(output, hash_value, ELL_SIZE);

    DMSG("  Exit functionL.\n");
    return res;
}

// Function K. Computes K(m || c0 || c1).
_INLINE_ status_t functionK(
        OUT uint8_t * output,
        IN const uint8_t * m,
        IN const uint8_t * c0,
        IN const uint8_t * c1)
{
    status_t res = SUCCESS;
    sha384_hash_t large_hash = {0};

    // preparing buffer with: [m || c0 || c1]
    uint8_t tmp1[ELL_SIZE + 2*R_SIZE] = {0};
    memcpy(tmp1, m, ELL_SIZE);
    memcpy(tmp1 + ELL_SIZE, c0, R_SIZE);
    memcpy(tmp1 + ELL_SIZE + R_SIZE, c1, ELL_SIZE);

    // shared secret =  K(m || c0 || c1)
    // select hash function
    sha3_384(large_hash.raw, tmp1, 2*ELL_SIZE + R_SIZE);
    memcpy(output, large_hash.raw, ELL_SIZE);
  
    DMSG("  Exit functionK.\n");
    return res;
}

_INLINE_ status_t compute_syndrome(OUT syndrome_t* syndrome,
        IN const ct_t* ct,
        IN const sk_t* sk)
{
    status_t res = SUCCESS;
    uint8_t s_tmp_bytes[R_BITS] = {0};
    uint8_t s0[R_SIZE] = {0};

    // syndrome: s = c0*h0
    // poly_mod_mul(s0, ct->val0, sk->val0);
    poly_mod_mul_in8(s0, ct->val0, sk->val0);

    // store the syndrome in a bit array
    convertByteToBinary(s_tmp_bytes, s0, R_BITS);
    transpose(syndrome->raw, s_tmp_bytes);

    DMSG("  Exit compute_syndrome.\n");

    return res;
}

////////////////////////////////////////////////////////////////
//The three APIs below (keypair, enc, dec) are defined by NIST:
//In addition there are two KAT versions of this API as defined.
////////////////////////////////////////////////////////////////
int crypto_kem_keypair(OUT unsigned char *pk, OUT unsigned char *sk)
{
    //Convert to these implementation types
    sk_t* l_sk = (sk_t*)sk;
    pk_t* l_pk = (pk_t*)pk;

    // return code
    status_t res = SUCCESS;

    //For NIST DRBG_CTR
    double_seed_t seeds = {0};
    shake256_prng_state_t h_prng_state = {0};

    //Get the entropy seeds
    get_seeds_fortest(&seeds, KEYGEN_SEEDS);

    ////print out seeds
    //printf("seeds.s1.raw:\n");
    //for (int i = 0; i < 32; i++) {
    //    printf("%d,", seeds.s1.raw[i]);
    //}
    //printf("\nseeds.s2.raw:\n");
    //for (int i = 0; i < 32; i++) {
    //    printf("%d,", seeds.s2.raw[i]);
    //}

    // sk = (h0, h1, sigma)
    uint8_t * h0 = l_sk->val0;
    uint8_t * h1 = l_sk->val1;
    uint8_t * sigma = l_sk->sigma;

    uint8_t inv_h0[R_SIZE] = {0};

    DMSG("  Enter crypto_kem_keypair.\n");
    DMSG("    Calculating the secret key.\n");

    shake256_init(seeds.s1.raw, ELL_SIZE, &h_prng_state);
    res = generate_sparse_rep_keccak(h0, DV, R_BITS, &h_prng_state); CHECK_STATUS(res);
    res = generate_sparse_rep_keccak(h1, DV, R_BITS, &h_prng_state); CHECK_STATUS(res);

    //print out some variables
    /*printf("h0:\n");
    for (int i = 0; i < R_SIZE; i++) {
        printf("%d,",h0[i]);
        if ((i + 1) % 64 == 0) {
            printf("\n");
        }
    }*/

    // use the second seed as sigma
    memcpy(sigma, seeds.s2.raw, ELL_SIZE);

    DMSG("    Calculating the public key.\n");

    // pk = (1, h1*h0^(-1)), the first pk component (1) is implicitly assumed
    // poly_mod_inv(inv_h0, h0);
    gf2x_mod_inv_wrapper(inv_h0,h0);
    // poly_mod_mul(l_pk->val, inv_h0, h1);
    poly_mod_mul_in8(l_pk->val, inv_h0, h1);

    EXIT:
    DMSG("  Exit crypto_kem_keypair.\n");
    return res;
}

//Encapsulate - pk is the public key,
//              ct is a key encapsulation message (ciphertext),
//              ss is the shared secret.
int crypto_kem_enc(OUT unsigned char *ct,
        OUT unsigned char *ss,
        IN  const unsigned char *pk)
{
    DMSG("  Enter crypto_kem_enc.\n");

    status_t res = SUCCESS;

    //Convert to these implementation types
    const pk_t* l_pk = (pk_t*)pk;
    ct_t* l_ct = (ct_t*)ct;
    ss_t* l_ss = (ss_t*)ss;

    //For NIST DRBG_CTR.
    double_seed_t seeds = {0};

    //Get the entropy seeds.
    get_seeds_fortest(&seeds, ENCAPS_SEEDS);

    // quantity m:
    uint8_t m[ELL_SIZE] = {0};

    // error vector:
    uint8_t e[N_SIZE] = {0};
    uint8_t e0[R_SIZE] = {0};
    uint8_t e1[R_SIZE] = {0};

    // temporary buffer:
    uint8_t tmp[ELL_SIZE] = {0};

    //random data generator; Using seed s1
    memcpy(m, seeds.s1.raw, ELL_SIZE);

    // (e0, e1) = H(m)
    functionH(e, m);
    poly_split_polynomial(e0, e1, e);

    // ct = (c0, c1) = (e0 + e1*h, L(e0, e1) \XOR m)
    // poly_mod_mul(l_ct->val0, l_pk->val,e1);
    poly_mod_mul_in8(l_ct->val0, l_pk->val,e1);
    poly_add(l_ct->val0, l_ct->val0, e0);
    functionL(tmp, e);
    for (uint32_t i = 0; i < ELL_SIZE; i++)
        l_ct->val1[i] = tmp[i] ^ m[i];

    // Function K:
    //shared secret =  K(m || c0 || c1)
    functionK(l_ss->raw, m, l_ct->val0, l_ct->val1);


    EXIT:

    DMSG("  Exit crypto_kem_enc.\n");
    return res;
}

//Decapsulate - ct is a key encapsulation message (ciphertext),
//              sk is the private key,
//              ss is the shared secret
int crypto_kem_dec(OUT unsigned char *ss,
        IN const unsigned char *ct,
        IN const unsigned char *sk)
{
    DMSG("  Enter crypto_kem_dec.\n");
    status_t res = SUCCESS;

    // convert to this implementation types
    const sk_t* l_sk = (sk_t*)sk;
    const ct_t* l_ct = (ct_t*)ct;
    ss_t* l_ss = (ss_t*)ss;

    int failed = 0;

    // for NIST DRBG_CTR
    double_seed_t seeds = {0};
  
    uint8_t e_recomputed[N_SIZE] = {0};

    uint8_t Le0e1[ELL_SIZE + 2*R_SIZE] = {0};
    uint8_t m_prime[ELL_SIZE] = {0};

    uint32_t h0_compact[DV] = {0};
    uint32_t h1_compact[DV] = {0};

    uint8_t e_prime[N_SIZE] = {0};
    uint8_t e_twoprime[R_BITS*2] = {0};
   
    uint8_t e_tmp1[R_BITS*2] = {0};
    uint8_t e_tmp2[N_SIZE] = {0};

    uint8_t e0rand[R_SIZE] = {0};
    uint8_t e1rand[R_SIZE] = {0};

    int rc;

    DMSG("  Converting to compact rep.\n");
    convert2compact(h0_compact, l_sk->val0);
    convert2compact(h1_compact, l_sk->val1);

    DMSG("  Computing s.\n");
    syndrome_t syndrome;

       // Step 1. computing syndrome:
    res = compute_syndrome(&syndrome, l_ct, l_sk); CHECK_STATUS(res);

    // Step 2. decoding:
    DMSG("  Decoding.\n");
    rc = BGF_decoder(e_tmp1, syndrome.raw, h0_compact, h1_compact);

    convertBinaryToByte(e_prime, e_tmp1, 2*R_BITS);

    // Step 3. compute L(e0 || e1)
    functionL(Le0e1, e_prime);

    // Step 4. retrieve m' = c1 \xor L(e0 || e1)
    for(uint32_t i = 0; i < ELL_SIZE; i++)
    {
        m_prime[i] = l_ct->val1[i] ^ Le0e1[i];
    }

    // Step 5. (e0, e1) = H(m)
    functionH(e_recomputed, m_prime);

    if (!safe_cmp(e_recomputed, e_prime, N_SIZE))
    {
        DMSG("recomputed error vector does not match decoded error vector\n");
        failed = 1;
    }

    // Step 6. compute shared secret k = K()
    if (failed) {
        // shared secret = K(sigma || c0 || c1)
        functionK(l_ss->raw, l_sk->sigma, l_ct->val0, l_ct->val1);
    }
    else
    {
       // shared secret = K(m' || c0 || c1)
        functionK(l_ss->raw, m_prime, l_ct->val0, l_ct->val1);
    }

    EXIT:

    DMSG("  Exit crypto_kem_dec.\n");
    return res;
}

//==================================================================
//                      Customized Version 
//==================================================================
 //Comparing value in a constant time manner.
uint32_t safe_cmp_custom(IN const uint8_t* a,
    IN const uint8_t* b,
    IN const uint32_t size)
{
    uint8_t res = 0;

    int32_t avl=size;
    size_t vl;
    const uint8_t* vsrc1_addr=a;
    const uint8_t* vsrc2_addr=b;
    vuint8m4_t vsrc1,vsrc2,vres1;
    vuint8m1_t vres2,vres3;
    vl=vsetvl_e8m1(avl);
    vres3=vmv_v_x_u8m1(0,vl);
    while(avl>0){
        vl=vsetvl_e8m4(avl);
        vsrc1=vle8_v_u8m4(vsrc1_addr,vl);
        vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
        vres1=vxor_vv_u8m4(vsrc1,vsrc2,vl);
        vres2=vredor_vs_u8m4_u8m1(vres2,vres1,vres2,vl);
        vres3=vor_vv_u8m1(vres3,vres2,1);
        
        vsrc1_addr+=vl,vsrc2_addr+=vl,avl-=vl;
    }
    res=vmv_x_s_u8m1_u8(vres3); 
    return (res == 0);
}

// Function H. It uses the extract-then-expand paradigm based on SHA384 and
// AES256-CTR PRNG to produce e from m.
_INLINE_ status_t functionH_custom(
        OUT uint8_t * e,
        IN const uint8_t * m)
{
    status_t res = SUCCESS;

    // use the seed to generate sparse error vector e:
    DMSG("    Generating random error.\n");
    // shake256_prng_state_t prng_state = {0};
    // shake256_init(seed_for_hash.raw, ELL_SIZE, &prng_state);
    // res = generate_sparse_rep_keccak(e, T1, N_BITS, &prng_state); CHECK_STATUS(res);
    setZero_custom(e,N_SIZE);

    uint8_t buffer[T1<<2]={0};
    shake256_prng_custom(buffer,T1<<2,m,ELL_SIZE);
    uint32_t* buffer_in32=(uint32_t*)buffer;

    uint32_t idx[T1]={0};
    for(int i=0;i<T1;i++)idx[i]=T1-i-1;
    uint32_t* idx_addr=idx;

    int32_t avl=T1;
    size_t vl=0;
    while(avl>0){
        vl=vsetvl_e32m4(avl);
        vuint32m4_t vidx=vle32_v_u32m4(idx_addr,vl);
        vuint32m4_t vtemp=vrsub_vx_u32m4(vidx,N_BITS,vl);
        vuint32m4_t vbuffer=vle32_v_u32m4(buffer_in32,vl);
        vuint64m8_t vmul=vwmulu_vv_u64m8(vbuffer,vtemp,vl);
        vbuffer=vnsrl_wx_u32m4(vmul,32,vl);
        vbuffer=vadd_vv_u32m4(vbuffer,vidx,vl);
        vse32_v_u32m4(buffer_in32,vbuffer,vl);

        //Address update
        buffer_in32+=vl,idx_addr+=vl,avl-=vl;
    }
    buffer_in32=(uint32_t*)buffer;
    for(int i=T1-1;i>=0;i--){
        uint32_t rand_pos=buffer_in32[T1-i-1];
        if (CHECK_BIT(e, rand_pos))
        {
            // If collision, then select index i instead of rand_pos
            rand_pos = i;
        }
        SET_BIT(e, rand_pos);
    }

    EXIT:
    DMSG("  Exit functionH.\n");
    return res;
}

// Function L. Computes L(e0 || e1)
_INLINE_ status_t functionL_custom(
        OUT uint8_t * output,
        IN const uint8_t * e)
{
    status_t res = SUCCESS;
    uint8_t hash_value[SHA384_HASH_SIZE] = {0};
    uint8_t e_split[2 * R_SIZE] = {0};

    poly_split_polynomial_custom(e_split, &e_split[R_SIZE], e);

    // select hash function
    sha3_384_custom(hash_value, e_split, 2*R_SIZE);

    memcpy(output, hash_value, ELL_SIZE);

    DMSG("  Exit functionL.\n");
    return res;
}

// Function K. Computes K(m || c0 || c1).
_INLINE_ status_t functionK_custom(
        OUT uint8_t * output,
        IN const uint8_t * m,
        IN const uint8_t * c0,
        IN const uint8_t * c1)
{
    status_t res = SUCCESS;
    sha384_hash_t large_hash = {0};

    // preparing buffer with: [m || c0 || c1]
    uint8_t tmp1[ELL_SIZE + 2*R_SIZE] = {0};
    memcpy(tmp1, m, ELL_SIZE);
    memcpy(tmp1 + ELL_SIZE, c0, R_SIZE);
    memcpy(tmp1 + ELL_SIZE + R_SIZE, c1, ELL_SIZE);

    // shared secret =  K(m || c0 || c1)
    // select hash function
    sha3_384_custom(large_hash.raw, tmp1, 2*ELL_SIZE + R_SIZE);
    memcpy(output, large_hash.raw, ELL_SIZE);
  
    DMSG("  Exit functionK.\n");
    return res;
}

_INLINE_ status_t compute_syndrome_custom(OUT syndrome_t* syndrome,
        IN const ct_t* ct,
        IN const sk_t* sk)
{
    status_t res = SUCCESS;
    uint8_t s_tmp_bytes[R_BITS] = {0};
    uint8_t s0[R_SIZE] = {0};

    // syndrome: s = c0*h0
    poly_mul_binary_custom(s0, ct->val0, sk->val0);

    // store the syndrome in a bit array
    convertByteToBinary_custom(s_tmp_bytes, s0, R_BITS);
    transpose_custom(syndrome->raw, s_tmp_bytes);

    DMSG("  Exit compute_syndrome.\n");

    return res;
}

int crypto_kem_keypair_custom(OUT unsigned char *pk, OUT unsigned char *sk)
{
    //Convert to these implementation types
    sk_t* l_sk = (sk_t*)sk;
    pk_t* l_pk = (pk_t*)pk;

    // return code
    status_t res = SUCCESS;

    //For NIST DRBG_CTR
    double_seed_t seeds = {0};

    //Get the entropy seeds
    get_seeds_fortest(&seeds, KEYGEN_SEEDS);

    // sk = (h0, h1, sigma)
    uint8_t * h0 = l_sk->val0;
    uint8_t * h1 = l_sk->val1;
    uint8_t * sigma = l_sk->sigma;

    uint8_t inv_h0[R_SIZE] = {0};

    DMSG("  Enter crypto_kem_keypair.\n");
    DMSG("    Calculating the secret key.\n");

    //Generate h0 and h1
    setZero_custom(h0,R_SIZE);
    setZero_custom(h1,R_SIZE);

    uint8_t buffer[DV<<3]={0};
    shake256_prng_custom(buffer,DV<<3,seeds.s1.raw,ELL_SIZE);
    uint32_t* buffer_in32=(uint32_t*)buffer;

    uint32_t idx[DV]={0};
    for(int i=0;i<DV;i++)idx[i]=DV-i-1;

    for(int i=0;i<2;i++){
        uint32_t* idx_addr=idx;
        int32_t avl=DV;
        size_t vl=0;
        while(avl>0){
            vl=vsetvl_e32m4(avl);
            vuint32m4_t vidx=vle32_v_u32m4(idx_addr,vl);
            vuint32m4_t vtemp=vrsub_vx_u32m4(vidx,R_BITS,vl);
            vuint32m4_t vbuffer=vle32_v_u32m4(buffer_in32,vl);
            vuint64m8_t vmul=vwmulu_vv_u64m8(vbuffer,vtemp,vl);
            vbuffer=vnsrl_wx_u32m4(vmul,32,vl);
            vbuffer=vadd_vv_u32m4(vbuffer,vidx,vl);
            vse32_v_u32m4(buffer_in32,vbuffer,vl);

            //Address update
            buffer_in32+=vl,idx_addr+=vl,avl-=vl;
        }
    }
    buffer_in32=(uint32_t*)buffer;
    for(int i=DV-1;i>=0;i--){
        uint32_t rand_pos=buffer_in32[DV-i-1];
        if (CHECK_BIT(h0, rand_pos))
        {
            // If collision, then select index i instead of rand_pos
            rand_pos = i;
        }
        SET_BIT(h0, rand_pos);
    }
    for(int i=DV-1;i>=0;i--){
        uint32_t rand_pos=buffer_in32[DV+DV-i-1];
        if (CHECK_BIT(h1, rand_pos))
        {
            // If collision, then select index i instead of rand_pos
            rand_pos = i;
        }
        SET_BIT(h1, rand_pos);
    }

    // use the second seed as sigma
    memcpy(sigma, seeds.s2.raw, ELL_SIZE);

    //printf("    Calculating the public key.\n");//for test

    // pk = (1, h1*h0^(-1)), the first pk component (1) is implicitly assumed
    //poly_inv_binary_custom64_wrapper(inv_h0, h0);
    gf2x_mod_inv_wrapper_custom(inv_h0, h0);
    poly_mul_binary_custom(l_pk->val, inv_h0, h1);

    EXIT:
    DMSG("  Exit crypto_kem_keypair.\n");
    return res;
}

//Encapsulate - pk is the public key,
//              ct is a key encapsulation message (ciphertext),
//              ss is the shared secret.
int crypto_kem_enc_custom(OUT unsigned char *ct,
        OUT unsigned char *ss,
        IN  const unsigned char *pk)
{
    DMSG("  Enter crypto_kem_enc.\n");

    status_t res = SUCCESS;

    //Convert to these implementation types
    const pk_t* l_pk = (pk_t*)pk;
    ct_t* l_ct = (ct_t*)ct;
    ss_t* l_ss = (ss_t*)ss;

    //For NIST DRBG_CTR.
    double_seed_t seeds = {0};

    //Get the entropy seeds.
    get_seeds_fortest(&seeds, ENCAPS_SEEDS);

    // quantity m:
    uint8_t m[ELL_SIZE] = {0};

    // error vector:
    uint8_t e[N_SIZE] = {0};
    uint8_t e0[R_SIZE] = {0};
    uint8_t e1[R_SIZE] = {0};

    // temporary buffer:
    uint8_t tmp[ELL_SIZE] = {0};

    //random data generator; Using seed s1
    memcpy(m, seeds.s1.raw, ELL_SIZE);

    // (e0, e1) = H(m)
    functionH_custom(e, m);
    poly_split_polynomial_custom(e0, e1, e);

    // ct = (c0, c1) = (e0 + e1*h, L(e0, e1) \XOR m)
    poly_mul_binary_custom(l_ct->val0, l_pk->val,e1);
    poly_add_custom(l_ct->val0, l_ct->val0, e0,R_SIZE);
    functionL_custom(tmp, e);
    poly_add_custom(l_ct->val1,tmp,m,ELL_SIZE);

    // Function K:
    //shared secret =  K(m || c0 || c1)
    functionK_custom(l_ss->raw, m, l_ct->val0, l_ct->val1);


    EXIT:

    DMSG("  Exit crypto_kem_enc.\n");
    return res;
}

//Decapsulate - ct is a key encapsulation message (ciphertext),
//              sk is the private key,
//              ss is the shared secret
int crypto_kem_dec_custom(OUT unsigned char *ss,
        IN const unsigned char *ct,
        IN const unsigned char *sk)
{
    DMSG("  Enter crypto_kem_dec.\n");
    status_t res = SUCCESS;

    // convert to this implementation types
    const sk_t* l_sk = (sk_t*)sk;
    const ct_t* l_ct = (ct_t*)ct;
    ss_t* l_ss = (ss_t*)ss;

    int failed = 0;

    // for NIST DRBG_CTR
    double_seed_t seeds = {0};
  
    uint8_t e_recomputed[N_SIZE] = {0};

    uint8_t Le0e1[ELL_SIZE + 2*R_SIZE] = {0};
    uint8_t m_prime[ELL_SIZE] = {0};

    uint32_t h0_compact[DV] = {0};
    uint32_t h1_compact[DV] = {0};

    uint8_t e_prime[N_SIZE] = {0};
    uint8_t e_twoprime[R_BITS*2] = {0};
   
    uint8_t e_tmp1[R_BITS*2] = {0};
    uint8_t e_tmp2[N_SIZE] = {0};

    uint8_t e0rand[R_SIZE] = {0};
    uint8_t e1rand[R_SIZE] = {0};

    int rc;

    DMSG("  Converting to compact rep.\n");
    convert2compact(h0_compact, l_sk->val0);
    convert2compact(h1_compact, l_sk->val1);

    DMSG("  Computing s.\n");
    syndrome_t syndrome;

       // Step 1. computing syndrome:
    res = compute_syndrome_custom(&syndrome, l_ct, l_sk); CHECK_STATUS(res);
    
    // Step 2. decoding:
    DMSG("  Decoding.\n");
    rc = BGF_decoder_custom(e_tmp1, syndrome.raw, h0_compact, h1_compact);

    convertBinaryToByte_custom(e_prime, e_tmp1, 2*R_BITS);
    
    // Step 3. compute L(e0 || e1)
    functionL_custom(Le0e1, e_prime);
    
    // Step 4. retrieve m' = c1 \xor L(e0 || e1)
    poly_add_custom(m_prime,l_ct->val1,Le0e1,ELL_SIZE);
    
    // Step 5. (e0, e1) = H(m)
    functionH_custom(e_recomputed, m_prime);
    
    if (!safe_cmp(e_recomputed, e_prime, N_SIZE))
    {
        DMSG("recomputed error vector does not match decoded error vector\n");
        failed = 1;
    }

    // Step 6. compute shared secret k = K()
    if (failed) {
        // shared secret = K(sigma || c0 || c1)
        functionK_custom(l_ss->raw, l_sk->sigma, l_ct->val0, l_ct->val1);
    }
    else
    {
       // shared secret = K(m' || c0 || c1)
        functionK_custom(l_ss->raw, m_prime, l_ct->val0, l_ct->val1);
    }

    printf("\n");
    EXIT:

    DMSG("  Exit crypto_kem_dec.\n");
    return res;
}

//==================================================================
//                     Test For Customization 
//==================================================================
#include <time.h>
#include <stdbool.h>
void test_safe_cmp_custom(){
    size_t vl=vl=vsetvl_e8m1(ELEMENT_SEW8_PER_VECREG);
    uint8_t a[10]={0};
    uint8_t b[10]={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<10;i++){
        a[i]=b[i]=rand()&255;
    }
    b[0]=123;
    uint32_t ref,res;

    // ref=safe_cmp(a,b,N_SIZE);
    res=safe_cmp_custom(a,b,10);

    printf("res=%d\n",res);
    return;
}

void test_functionH_custom(){
    uint8_t m[ELL_SIZE]={0};
    uint8_t e_ref[N_SIZE];
    uint8_t e_res[N_SIZE];
    srand((unsigned)time(NULL));
    for(int i=0;i<ELL_SIZE;i++)m[i]=rand()&255;
    functionH(e_ref,m);
    functionH_custom(e_res,m);
    bool flag=true;
    for(int i=0;i<N_SIZE;i++){
        if(e_ref[i]!=e_res[i]){
            // printf("e_ref[%d]=%d,e_res[%d]=%d\n",i,e_ref[i],i,e_res[i]);
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_functionH_custom");
    else printf("fail in test_functionH_custom");
    return;
}

void test_functionL_custom(){
    uint8_t tmp_ref[ELL_SIZE]={0};
    uint8_t tmp_res[ELL_SIZE]={0};
    uint8_t e[N_SIZE]={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<N_SIZE;i++)e[i]=rand()&255;

    functionL(tmp_ref,e);
    functionL_custom(tmp_res,e);

    bool flag=true;
    for(int i=0;i<ELL_SIZE;i++){
        if(tmp_ref[i]!=tmp_res[i]){
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_functionL_custom");
    else printf("fail in test_functionL_custom");
    return;
}

void test_functionK_custom(){
    uint8_t c0[R_SIZE]={0};
    uint8_t c1[ELL_SIZE]={0};
    uint8_t m[ELL_SIZE]={0};

    uint8_t ref[ELL_SIZE]={0};
    uint8_t res[ELL_SIZE]={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<ELL_SIZE;i++){
        c1[i]=rand()&255;
        m[i]=rand()&255;
    }
    for(int i=0;i<R_SIZE;i++)c0[i]=rand()&255;

    functionK(ref,m,c0,c1);
    functionK_custom(res,m,c0,c1);

    bool flag=true;
    for(int i=0;i<ELL_SIZE;i++){
        if(ref[i]!=res[i]){
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_functionK_custom");
    else printf("fail in test_functionK_custom");
    return;
}

void test_compute_syndrome_custom(){
    sk_t l_sk;
    ct_t l_ct;

    syndrome_t syndrome_ref,syndrome_res;

    srand((unsigned)time(NULL));
    for(int i=0;i<R_SIZE;i++){
        l_ct.val0[i]=rand()&255;
        l_sk.val0[i]=0;
    }
    l_ct.val0[R_SIZE-1]&=((1 << (R_BITS & 7)) - 1);

    for(int i=0;i<DV;i++){
        uint16_t idx=(((uint32_t)rand()<<16)|(uint32_t)rand())%R_BITS;
        SET_BIT(l_sk.val0,idx);
    }

    compute_syndrome(&syndrome_ref,&l_ct,&l_sk);
    compute_syndrome_custom(&syndrome_res,&l_ct,&l_sk);

    bool flag=true;
    for(int i=0;i<R_BITS;i++){
        if(syndrome_ref.raw[i]!=syndrome_res.raw[i]){
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_compute_syndrome_custom");
    else printf("fail in test_compute_syndrome_custom");
    return;
}

void test_crypto_kem_keypair_custom(){
    uint8_t pk_ref[CRYPTO_PUBLICKEYBYTES];
    uint8_t pk_res[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk_ref[CRYPTO_SECRETKEYBYTES];
    uint8_t sk_res[CRYPTO_SECRETKEYBYTES];

    crypto_kem_keypair(pk_ref,sk_ref);
    crypto_kem_keypair_custom(pk_res,sk_res);

    bool flag=true;
    for(int i=0;i<CRYPTO_PUBLICKEYBYTES;i++){
        if(pk_ref[i]!=pk_res[i]){
            flag=false;
            break;
        }
    }
    for(int i=0;i<CRYPTO_SECRETKEYBYTES;i++){
        if(sk_ref[i]!=sk_res[i]){
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_crypto_kem_keypair_custom");
    else printf("fail in test_crypto_kem_keypair_custom");
    return;
}

void test_crypto_kem_keypair_custom_level1(){
    assert(R_BITS==12323);
    uint8_t pk_ref[CRYPTO_PUBLICKEYBYTES]={
        0,41,36,24,59,191,168,145,16,2,123,118,228,16,124,231,18,196,179,98,181,99,19,158,221,7,33,171,153,208,176,140,242,190,208,3,150,76,170,212,173,248,142,111,46,70,15,123,42,58,29,183,180,141,85,145,136,139,220,174,135,226,237,224,
        174,239,214,232,85,13,33,14,103,242,135,179,146,4,174,0,134,245,178,243,245,236,167,160,13,215,98,201,89,226,26,195,66,255,53,77,186,91,156,245,109,65,32,55,8,158,74,252,249,168,102,229,155,22,78,158,130,76,53,106,50,191,236,193,
        6,89,215,38,203,108,46,141,134,85,191,226,60,93,250,138,229,31,216,110,223,91,143,0,56,250,203,52,5,82,76,114,194,210,158,187,147,17,122,118,250,165,25,59,99,96,39,160,225,38,48,139,94,34,157,216,89,212,75,123,220,97,78,77,
        217,143,119,93,117,214,171,184,92,125,192,87,162,18,163,160,235,77,171,86,41,110,82,177,84,68,227,79,233,57,65,19,13,164,160,211,91,35,26,17,246,193,127,127,171,200,16,234,230,33,147,218,227,209,106,35,198,160,112,136,3,238,85,209,
        77,167,18,35,84,45,192,156,191,72,58,100,45,207,154,187,192,63,184,17,98,229,5,118,244,237,243,81,180,18,3,128,84,196,38,227,44,57,114,96,249,125,232,179,23,6,22,190,28,236,108,68,53,206,82,5,7,180,176,219,221,49,149,35,
        135,17,139,93,120,29,255,113,160,205,46,241,67,91,106,254,166,73,203,55,66,253,124,120,198,238,84,96,113,11,239,97,167,148,131,32,76,93,199,140,253,163,156,28,80,148,84,163,205,131,149,101,40,45,240,156,129,189,210,237,222,231,117,54,
        215,120,76,32,152,58,89,147,133,184,187,212,121,65,84,62,116,83,117,155,89,193,248,190,99,82,123,57,190,48,244,201,43,60,114,160,37,109,199,119,70,102,94,71,251,5,99,41,88,29,231,50,251,69,76,136,215,251,96,164,153,149,58,239,
        34,212,3,103,182,254,11,90,7,176,234,23,147,247,154,233,69,239,81,67,41,39,168,146,211,111,160,96,217,191,171,136,255,235,2,223,154,100,57,200,84,254,133,75,208,251,59,141,225,117,152,33,85,77,72,89,46,206,220,168,32,102,70,114,
        94,18,159,119,116,65,38,219,108,29,255,230,169,54,143,50,81,116,241,188,203,30,71,221,174,56,55,185,117,83,87,138,135,60,134,249,145,215,164,134,131,172,56,247,86,207,3,211,7,59,24,81,214,243,239,155,196,174,239,194,33,127,91,225,
        30,232,202,229,184,164,39,33,12,155,105,161,117,109,5,203,65,127,130,112,144,153,220,141,220,24,178,254,120,230,161,161,165,207,210,217,129,5,62,30,26,199,251,160,196,0,75,165,198,252,122,214,54,117,82,101,48,239,234,220,144,230,121,216,
        110,205,149,237,22,49,102,215,28,210,86,248,59,159,6,169,161,50,84,99,251,2,231,187,201,40,91,12,108,119,213,213,193,73,128,106,100,73,255,169,242,47,40,186,207,195,251,161,5,221,85,87,183,159,48,108,163,70,59,201,37,50,185,236,
        86,20,158,135,228,191,169,2,52,86,169,3,238,223,143,23,252,78,230,193,45,179,152,249,106,14,158,210,221,22,199,52,130,240,10,2,196,241,63,25,204,163,201,170,2,218,176,160,35,209,87,68,58,215,20,57,92,248,242,19,236,109,227,204,
        252,235,213,208,161,27,222,2,33,152,56,41,14,108,239,243,33,37,239,175,110,124,41,240,44,65,99,135,45,184,235,13,5,224,232,48,250,237,115,156,164,121,146,156,252,197,233,1,9,10,61,143,185,9,218,59,149,0,16,103,146,209,173,225,
        189,23,126,129,48,144,69,157,215,187,159,185,18,82,233,11,24,196,206,59,104,194,158,3,27,200,124,123,69,137,97,19,65,143,112,196,100,5,67,56,63,146,68,213,183,130,61,138,2,253,164,18,188,176,101,136,175,243,179,79,139,39,132,145,
        183,146,129,97,35,2,248,239,72,128,31,248,210,102,53,191,161,21,44,244,129,132,34,79,127,19,210,110,214,129,105,140,92,255,74,119,141,40,22,231,56,133,21,218,82,102,71,136,118,209,242,124,224,77,165,148,34,40,55,97,33,210,229,153,
        30,240,241,67,44,1,69,148,125,196,206,227,158,210,160,0,109,248,247,150,246,152,29,7,94,100,125,8,238,242,101,185,140,196,162,205,198,213,94,226,120,60,78,207,230,163,41,166,129,95,4,253,15,35,103,5,15,142,160,253,209,29,92,230,
        246,189,196,193,7,216,188,228,243,23,145,43,216,50,182,209,191,39,186,215,195,118,105,20,126,219,254,238,55,133,146,255,181,58,58,254,165,46,132,163,205,225,5,149,14,102,48,68,55,46,152,228,170,175,139,96,211,190,42,84,116,99,186,0,
        151,27,218,4,226,60,234,158,131,66,242,184,173,159,140,148,235,246,57,145,32,133,107,40,198,169,95,183,60,81,236,183,190,220,98,137,202,116,144,115,178,58,44,14,217,72,112,68,182,179,46,139,194,180,229,160,68,251,59,70,96,151,89,28,
        34,221,49,221,192,26,228,113,9,158,229,105,127,114,152,167,58,57,2,95,34,186,70,50,57,59,10,155,29,77,13,39,97,73,37,136,155,214,76,46,145,177,111,210,70,5,131,64,118,133,58,19,101,194,247,199,167,11,0,134,236,185,169,32,
        208,122,40,22,12,13,94,45,195,56,209,27,26,36,143,42,183,131,121,207,54,112,248,209,6,123,57,64,252,193,114,112,137,84,190,254,138,173,172,166,215,183,63,43,159,10,0,124,47,114,223,160,240,172,35,44,86,19,13,30,163,91,3,95,
        64,146,162,61,188,143,30,40,245,141,173,195,79,232,150,56,69,48,204,234,122,125,69,207,217,128,31,186,170,67,220,48,231,251,201,98,106,122,104,59,77,151,148,22,2,178,158,222,165,2,158,174,138,230,233,37,185,233,125,105,76,1,214,35,
        77,4,32,144,81,185,149,133,128,96,130,205,232,15,46,101,144,0,89,217,48,175,66,165,59,35,21,88,40,26,191,66,236,226,129,82,171,132,247,83,77,169,124,35,9,237,2,250,63,41,43,235,181,145,106,104,210,245,10,75,220,96,97,150,
        180,211,184,76,237,125,81,215,126,51,176,2,228,98,132,254,43,194,188,206,85,141,188,169,48,243,55,8,243,226,68,163,8,192,84,77,167,249,244,38,212,53,84,70,117,102,22,52,229,251,144,216,217,150,10,177,17,20,215,62,95,22,240,252,
        178,108,241,107,101,8,186,41,222,95,188,201,214,139,66,75,220,188,242,112,127,179,57,246,118,51,242,56,98,123,143,162,79,175,157,220,196,98,128,78,224,122,132,147,21,76,113,160,91,191,104,152,56,89,167,191,206,141,84,144,100,130,14,125,
        24,70,168,129,6
    };
    uint8_t pk_res[CRYPTO_PUBLICKEYBYTES];
    uint8_t sk_ref[CRYPTO_SECRETKEYBYTES]={
        0,0,0,0,0,128,1,0,0,0,64,0,2,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,32,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,8,0,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,
        128,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,128,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,16,16,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,1,0,0,0,64,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,32,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,4,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,16,128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        4,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,64,0,2,0,0,
        0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,16,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,64,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,4,0,0,0,4,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,64,32,0,0,0,0,0,0,8,0,0,0,0,32,0,2,0,0,2,0,0,0,0,0,0,0,0,0,4,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,64,0,0,0,
        0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,
        0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,32,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,64,0,128,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,2,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,
        0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,32,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,4,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,16,0,0,4,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31
    };
    uint8_t sk_res[CRYPTO_SECRETKEYBYTES];

    crypto_kem_keypair_custom(pk_res,sk_res);

    bool flag=true;
    for(int i=0;i<CRYPTO_PUBLICKEYBYTES;i++){
        if(pk_ref[i]!=pk_res[i]){
            flag=false;
            break;
        }
    }
    for(int i=0;i<CRYPTO_SECRETKEYBYTES;i++){
        if(sk_ref[i]!=sk_res[i]){
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_crypto_kem_keypair_custom");
    else printf("fail in test_crypto_kem_keypair_custom");
    return;
}

void test_crypto_kem_enc_custom_level1(){
    assert(R_BITS==12323);
    uint8_t pk_ref[CRYPTO_PUBLICKEYBYTES]={
        0,41,36,24,59,191,168,145,16,2,123,118,228,16,124,231,18,196,179,98,181,99,19,158,221,7,33,171,153,208,176,140,242,190,208,3,150,76,170,212,173,248,142,111,46,70,15,123,42,58,29,183,180,141,85,145,136,139,220,174,135,226,237,224,
        174,239,214,232,85,13,33,14,103,242,135,179,146,4,174,0,134,245,178,243,245,236,167,160,13,215,98,201,89,226,26,195,66,255,53,77,186,91,156,245,109,65,32,55,8,158,74,252,249,168,102,229,155,22,78,158,130,76,53,106,50,191,236,193,
        6,89,215,38,203,108,46,141,134,85,191,226,60,93,250,138,229,31,216,110,223,91,143,0,56,250,203,52,5,82,76,114,194,210,158,187,147,17,122,118,250,165,25,59,99,96,39,160,225,38,48,139,94,34,157,216,89,212,75,123,220,97,78,77,
        217,143,119,93,117,214,171,184,92,125,192,87,162,18,163,160,235,77,171,86,41,110,82,177,84,68,227,79,233,57,65,19,13,164,160,211,91,35,26,17,246,193,127,127,171,200,16,234,230,33,147,218,227,209,106,35,198,160,112,136,3,238,85,209,
        77,167,18,35,84,45,192,156,191,72,58,100,45,207,154,187,192,63,184,17,98,229,5,118,244,237,243,81,180,18,3,128,84,196,38,227,44,57,114,96,249,125,232,179,23,6,22,190,28,236,108,68,53,206,82,5,7,180,176,219,221,49,149,35,
        135,17,139,93,120,29,255,113,160,205,46,241,67,91,106,254,166,73,203,55,66,253,124,120,198,238,84,96,113,11,239,97,167,148,131,32,76,93,199,140,253,163,156,28,80,148,84,163,205,131,149,101,40,45,240,156,129,189,210,237,222,231,117,54,
        215,120,76,32,152,58,89,147,133,184,187,212,121,65,84,62,116,83,117,155,89,193,248,190,99,82,123,57,190,48,244,201,43,60,114,160,37,109,199,119,70,102,94,71,251,5,99,41,88,29,231,50,251,69,76,136,215,251,96,164,153,149,58,239,
        34,212,3,103,182,254,11,90,7,176,234,23,147,247,154,233,69,239,81,67,41,39,168,146,211,111,160,96,217,191,171,136,255,235,2,223,154,100,57,200,84,254,133,75,208,251,59,141,225,117,152,33,85,77,72,89,46,206,220,168,32,102,70,114,
        94,18,159,119,116,65,38,219,108,29,255,230,169,54,143,50,81,116,241,188,203,30,71,221,174,56,55,185,117,83,87,138,135,60,134,249,145,215,164,134,131,172,56,247,86,207,3,211,7,59,24,81,214,243,239,155,196,174,239,194,33,127,91,225,
        30,232,202,229,184,164,39,33,12,155,105,161,117,109,5,203,65,127,130,112,144,153,220,141,220,24,178,254,120,230,161,161,165,207,210,217,129,5,62,30,26,199,251,160,196,0,75,165,198,252,122,214,54,117,82,101,48,239,234,220,144,230,121,216,
        110,205,149,237,22,49,102,215,28,210,86,248,59,159,6,169,161,50,84,99,251,2,231,187,201,40,91,12,108,119,213,213,193,73,128,106,100,73,255,169,242,47,40,186,207,195,251,161,5,221,85,87,183,159,48,108,163,70,59,201,37,50,185,236,
        86,20,158,135,228,191,169,2,52,86,169,3,238,223,143,23,252,78,230,193,45,179,152,249,106,14,158,210,221,22,199,52,130,240,10,2,196,241,63,25,204,163,201,170,2,218,176,160,35,209,87,68,58,215,20,57,92,248,242,19,236,109,227,204,
        252,235,213,208,161,27,222,2,33,152,56,41,14,108,239,243,33,37,239,175,110,124,41,240,44,65,99,135,45,184,235,13,5,224,232,48,250,237,115,156,164,121,146,156,252,197,233,1,9,10,61,143,185,9,218,59,149,0,16,103,146,209,173,225,
        189,23,126,129,48,144,69,157,215,187,159,185,18,82,233,11,24,196,206,59,104,194,158,3,27,200,124,123,69,137,97,19,65,143,112,196,100,5,67,56,63,146,68,213,183,130,61,138,2,253,164,18,188,176,101,136,175,243,179,79,139,39,132,145,
        183,146,129,97,35,2,248,239,72,128,31,248,210,102,53,191,161,21,44,244,129,132,34,79,127,19,210,110,214,129,105,140,92,255,74,119,141,40,22,231,56,133,21,218,82,102,71,136,118,209,242,124,224,77,165,148,34,40,55,97,33,210,229,153,
        30,240,241,67,44,1,69,148,125,196,206,227,158,210,160,0,109,248,247,150,246,152,29,7,94,100,125,8,238,242,101,185,140,196,162,205,198,213,94,226,120,60,78,207,230,163,41,166,129,95,4,253,15,35,103,5,15,142,160,253,209,29,92,230,
        246,189,196,193,7,216,188,228,243,23,145,43,216,50,182,209,191,39,186,215,195,118,105,20,126,219,254,238,55,133,146,255,181,58,58,254,165,46,132,163,205,225,5,149,14,102,48,68,55,46,152,228,170,175,139,96,211,190,42,84,116,99,186,0,
        151,27,218,4,226,60,234,158,131,66,242,184,173,159,140,148,235,246,57,145,32,133,107,40,198,169,95,183,60,81,236,183,190,220,98,137,202,116,144,115,178,58,44,14,217,72,112,68,182,179,46,139,194,180,229,160,68,251,59,70,96,151,89,28,
        34,221,49,221,192,26,228,113,9,158,229,105,127,114,152,167,58,57,2,95,34,186,70,50,57,59,10,155,29,77,13,39,97,73,37,136,155,214,76,46,145,177,111,210,70,5,131,64,118,133,58,19,101,194,247,199,167,11,0,134,236,185,169,32,
        208,122,40,22,12,13,94,45,195,56,209,27,26,36,143,42,183,131,121,207,54,112,248,209,6,123,57,64,252,193,114,112,137,84,190,254,138,173,172,166,215,183,63,43,159,10,0,124,47,114,223,160,240,172,35,44,86,19,13,30,163,91,3,95,
        64,146,162,61,188,143,30,40,245,141,173,195,79,232,150,56,69,48,204,234,122,125,69,207,217,128,31,186,170,67,220,48,231,251,201,98,106,122,104,59,77,151,148,22,2,178,158,222,165,2,158,174,138,230,233,37,185,233,125,105,76,1,214,35,
        77,4,32,144,81,185,149,133,128,96,130,205,232,15,46,101,144,0,89,217,48,175,66,165,59,35,21,88,40,26,191,66,236,226,129,82,171,132,247,83,77,169,124,35,9,237,2,250,63,41,43,235,181,145,106,104,210,245,10,75,220,96,97,150,
        180,211,184,76,237,125,81,215,126,51,176,2,228,98,132,254,43,194,188,206,85,141,188,169,48,243,55,8,243,226,68,163,8,192,84,77,167,249,244,38,212,53,84,70,117,102,22,52,229,251,144,216,217,150,10,177,17,20,215,62,95,22,240,252,
        178,108,241,107,101,8,186,41,222,95,188,201,214,139,66,75,220,188,242,112,127,179,57,246,118,51,242,56,98,123,143,162,79,175,157,220,196,98,128,78,224,122,132,147,21,76,113,160,91,191,104,152,56,89,167,191,206,141,84,144,100,130,14,125,
        24,70,168,129,6
    };
    uint8_t ct_ref[CRYPTO_CIPHERTEXTBYTES]={
        62,128,124,76,119,113,248,122,178,232,155,28,102,59,108,211,47,231,150,72,36,197,0,242,157,48,95,74,86,116,60,95,166,221,74,176,21,63,245,255,141,184,64,52,137,211,198,80,80,19,228,132,101,84,136,83,201,2,59,98,175,151,192,166,
        171,185,69,102,6,186,68,229,235,221,152,167,159,222,24,245,130,46,64,9,107,14,205,79,187,178,119,208,125,98,0,71,18,207,215,211,188,70,70,49,9,166,219,235,232,168,104,214,13,198,186,159,129,128,13,37,63,50,41,223,136,61,114,47,
        194,126,27,179,29,128,221,57,173,58,113,152,198,119,199,44,33,76,148,16,25,85,116,156,146,12,177,157,158,120,60,99,12,72,232,201,133,9,146,212,87,81,112,222,211,244,74,137,251,22,101,94,172,8,53,27,130,32,213,105,109,42,27,241,
        64,93,45,60,101,251,178,230,49,32,124,199,1,66,137,191,73,221,103,121,163,73,84,16,125,12,34,103,241,95,215,29,137,232,5,81,123,151,168,167,97,119,207,20,151,214,46,40,144,36,203,151,228,90,210,118,53,73,153,7,175,54,79,95,
        28,114,49,34,252,77,240,235,51,198,51,74,16,207,203,241,135,42,89,211,19,157,186,68,2,17,110,65,184,51,47,31,72,211,189,50,87,228,128,45,206,138,128,18,184,249,150,120,21,95,109,35,99,210,139,35,99,87,144,40,163,206,162,55,
        25,186,245,213,70,37,145,136,219,198,101,254,148,152,174,150,42,156,151,141,37,193,17,46,84,162,210,181,7,241,114,4,46,151,232,199,242,36,188,80,254,159,121,216,36,13,109,146,146,80,36,244,191,2,29,101,19,136,228,73,203,80,204,192,
        67,134,239,124,51,89,223,185,70,52,154,42,130,213,101,219,255,255,108,150,48,4,87,88,217,227,51,130,167,101,78,125,179,56,190,16,244,214,34,166,39,233,143,128,132,225,126,211,83,236,137,133,42,5,77,90,192,107,253,5,183,89,17,17,
        9,149,172,255,200,154,57,123,58,71,148,103,112,58,88,123,154,157,86,138,122,82,129,250,99,209,144,169,230,202,71,205,94,48,196,35,214,23,88,152,181,152,190,145,147,182,199,146,182,213,204,141,164,123,45,99,51,187,1,19,230,19,63,143,
        137,198,241,161,63,8,137,244,220,151,220,125,14,87,154,121,154,177,72,74,100,93,57,188,124,76,183,13,63,58,8,194,191,170,65,113,222,248,225,77,89,127,145,62,15,145,31,31,180,248,209,139,99,23,129,153,158,186,141,148,79,165,194,85,
        132,73,23,186,23,27,178,152,118,9,215,35,149,167,21,67,182,200,92,204,179,188,60,130,103,6,12,169,76,132,144,195,29,127,38,165,215,47,51,236,237,110,40,33,154,202,151,172,20,142,85,209,47,66,101,8,173,180,232,127,75,255,217,63,
        181,111,214,149,76,87,32,137,15,36,224,208,157,124,87,132,119,147,29,56,233,3,66,99,244,48,129,115,94,254,237,55,204,181,114,251,243,60,16,22,20,239,231,183,37,68,139,247,238,250,113,241,28,221,177,249,140,45,28,182,15,245,90,44,
        81,211,186,108,17,46,9,147,96,95,175,219,180,207,176,171,228,146,21,191,130,57,230,149,252,62,53,4,147,43,186,62,33,37,25,116,12,20,8,100,194,115,249,179,186,3,219,7,86,15,40,236,45,191,65,170,207,1,194,162,7,7,29,4,
        28,48,77,119,68,42,237,100,157,15,18,95,26,216,170,242,112,202,213,192,180,21,43,164,173,194,159,237,84,188,246,129,235,176,91,224,69,226,134,116,140,118,187,137,56,142,16,35,11,219,22,40,203,217,110,65,134,167,3,243,235,205,106,249,
        68,244,20,193,12,163,183,242,85,216,103,68,252,218,65,32,15,95,120,26,204,253,4,142,46,194,40,136,197,78,227,151,203,13,49,79,211,120,175,2,148,173,153,201,64,29,62,60,130,114,218,201,241,205,251,159,71,88,42,74,116,177,109,233,
        194,232,243,71,218,8,139,130,7,202,137,85,114,105,137,232,126,173,178,158,56,219,67,34,95,41,210,31,226,240,38,190,138,72,120,244,146,80,163,38,225,212,196,46,16,193,211,203,49,202,123,248,228,106,238,150,128,104,58,121,94,178,32,89,
        199,214,232,133,109,242,2,205,175,224,62,146,223,73,234,205,191,17,119,209,176,231,222,193,236,177,3,92,79,181,125,212,204,98,63,58,65,241,138,65,35,202,50,160,162,7,122,38,66,88,27,17,50,88,220,53,190,131,210,110,236,58,199,0,
        6,223,92,83,5,245,27,223,248,106,156,140,93,59,117,242,169,132,249,95,125,226,240,100,147,52,155,95,206,186,185,0,74,102,213,15,89,167,97,3,87,144,217,208,80,207,61,20,137,72,230,83,114,191,245,75,111,26,43,38,152,61,5,142,
        34,226,96,100,125,180,164,125,164,178,232,73,192,207,174,25,54,2,0,6,73,22,70,18,194,52,100,161,156,154,156,17,11,140,155,129,37,74,207,36,57,179,79,116,193,77,187,244,3,84,185,26,144,80,132,31,10,36,36,141,86,171,37,198,
        204,25,215,151,118,111,220,29,253,90,96,80,244,17,25,211,147,147,103,32,72,66,10,44,14,250,125,210,150,179,84,18,68,198,249,9,102,186,179,242,246,47,23,148,48,199,189,203,253,49,78,171,177,200,205,156,236,84,102,236,90,119,155,81,
        36,181,11,167,77,152,21,122,151,156,110,168,189,209,188,138,202,85,107,248,220,125,79,69,176,220,109,193,237,109,122,20,107,121,125,3,100,104,130,24,248,84,3,36,7,104,231,130,108,11,141,117,34,164,197,252,130,213,221,162,206,200,73,219,
        3,199,144,146,217,104,31,121,29,158,236,239,199,74,16,148,142,151,0,92,130,252,57,56,210,46,3,174,33,205,178,236,208,16,29,123,228,21,25,218,54,244,110,81,150,231,127,100,50,17,36,242,15,119,80,45,223,166,101,75,184,107,141,175,
        178,82,139,151,164,220,131,95,56,118,10,140,231,209,26,160,38,140,30,139,232,231,13,181,50,103,40,155,51,31,131,69,53,68,127,255,136,77,178,143,39,73,79,80,2,244,43,59,249,183,132,57,80,98,97,233,20,219,181,226,113,145,221,100,
        170,43,9,51,83,15,177,253,248,140,217,174,50,176,189,57,204,36,10,123,99,153,75,9,248,194,20,132,8,247,211,173,202,165,193,6,149,235,204,80,165,13,90,249,65,202,100,90,16,208,51,129,27,43,214,95,200,0,165,200,56,137,99,218,
        216,216,53,19,229,213,98,77,129,180,198,191,244,244,33,89,213,186,121,36,190,58,10,144,241,106,136,28,157,62,96,243,10,210,244,212,224,250,145,226,63,168,19,159,34,162,249,175,57,22,192,17,120,186,204,20,145,125,129,28,172,170,116,84,
        180,103,32,175,0,6,146,117,18,215,223,126,195,114,24,65,57,233,226,219,140,36,238,107,41,126,131,123,154,24,157,72,64,7,10,64,16
    };
    uint8_t ss_ref[CRYPTO_BYTES]={
        197,114,75,90,23,55,251,197,171,63,74,206,215,24,6,58,228,23,221,37,48,160,107,176,240,135,139,129,192,102,61,26
    };
    uint8_t ct_res[CRYPTO_CIPHERTEXTBYTES]={0};
    uint8_t ss_res[CRYPTO_BYTES]={0};

    crypto_kem_enc_custom(ct_res,ss_res,pk_ref);

    bool flag=true;
    for(int i=0;i<CRYPTO_CIPHERTEXTBYTES;i++){
        if(ct_ref[i]!=ct_res[i]){
            flag=false;
            break;
        }
    }
    for(int i=0;i<CRYPTO_BYTES;i++){
        if(ss_ref[i]!=ss_res[i]){
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_crypto_kem_enc_custom_level1");
    else printf("fail in test_crypto_kem_enc_custom_level1");
    return;
}

void test_crypto_kem_dec_custom_level1(){
    assert(R_BITS==12323);
    uint8_t sk_ref[CRYPTO_SECRETKEYBYTES]={
        0,0,0,0,0,128,1,0,0,0,64,0,2,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,32,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,8,0,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,
        128,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,128,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,16,16,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,1,0,0,0,64,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,32,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,4,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,16,128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        4,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,64,0,2,0,0,
        0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,16,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,64,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,4,0,0,0,4,0,0,0,0,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,64,32,0,0,0,0,0,0,8,0,0,0,0,32,0,2,0,0,2,0,0,0,0,0,0,0,0,0,4,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,64,0,0,0,
        0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,8,0,0,0,0,0,0,0,0,
        0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,0,0,0,0,0,0,32,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,32,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,64,0,128,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,2,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,0,0,0,0,0,8,0,0,0,0,
        0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,32,64,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,16,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,32,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,4,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,128,16,0,0,4,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,128,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,16,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31
    };
    uint8_t ct_ref[CRYPTO_CIPHERTEXTBYTES]={
        62,128,124,76,119,113,248,122,178,232,155,28,102,59,108,211,47,231,150,72,36,197,0,242,157,48,95,74,86,116,60,95,166,221,74,176,21,63,245,255,141,184,64,52,137,211,198,80,80,19,228,132,101,84,136,83,201,2,59,98,175,151,192,166,
        171,185,69,102,6,186,68,229,235,221,152,167,159,222,24,245,130,46,64,9,107,14,205,79,187,178,119,208,125,98,0,71,18,207,215,211,188,70,70,49,9,166,219,235,232,168,104,214,13,198,186,159,129,128,13,37,63,50,41,223,136,61,114,47,
        194,126,27,179,29,128,221,57,173,58,113,152,198,119,199,44,33,76,148,16,25,85,116,156,146,12,177,157,158,120,60,99,12,72,232,201,133,9,146,212,87,81,112,222,211,244,74,137,251,22,101,94,172,8,53,27,130,32,213,105,109,42,27,241,
        64,93,45,60,101,251,178,230,49,32,124,199,1,66,137,191,73,221,103,121,163,73,84,16,125,12,34,103,241,95,215,29,137,232,5,81,123,151,168,167,97,119,207,20,151,214,46,40,144,36,203,151,228,90,210,118,53,73,153,7,175,54,79,95,
        28,114,49,34,252,77,240,235,51,198,51,74,16,207,203,241,135,42,89,211,19,157,186,68,2,17,110,65,184,51,47,31,72,211,189,50,87,228,128,45,206,138,128,18,184,249,150,120,21,95,109,35,99,210,139,35,99,87,144,40,163,206,162,55,
        25,186,245,213,70,37,145,136,219,198,101,254,148,152,174,150,42,156,151,141,37,193,17,46,84,162,210,181,7,241,114,4,46,151,232,199,242,36,188,80,254,159,121,216,36,13,109,146,146,80,36,244,191,2,29,101,19,136,228,73,203,80,204,192,
        67,134,239,124,51,89,223,185,70,52,154,42,130,213,101,219,255,255,108,150,48,4,87,88,217,227,51,130,167,101,78,125,179,56,190,16,244,214,34,166,39,233,143,128,132,225,126,211,83,236,137,133,42,5,77,90,192,107,253,5,183,89,17,17,
        9,149,172,255,200,154,57,123,58,71,148,103,112,58,88,123,154,157,86,138,122,82,129,250,99,209,144,169,230,202,71,205,94,48,196,35,214,23,88,152,181,152,190,145,147,182,199,146,182,213,204,141,164,123,45,99,51,187,1,19,230,19,63,143,
        137,198,241,161,63,8,137,244,220,151,220,125,14,87,154,121,154,177,72,74,100,93,57,188,124,76,183,13,63,58,8,194,191,170,65,113,222,248,225,77,89,127,145,62,15,145,31,31,180,248,209,139,99,23,129,153,158,186,141,148,79,165,194,85,
        132,73,23,186,23,27,178,152,118,9,215,35,149,167,21,67,182,200,92,204,179,188,60,130,103,6,12,169,76,132,144,195,29,127,38,165,215,47,51,236,237,110,40,33,154,202,151,172,20,142,85,209,47,66,101,8,173,180,232,127,75,255,217,63,
        181,111,214,149,76,87,32,137,15,36,224,208,157,124,87,132,119,147,29,56,233,3,66,99,244,48,129,115,94,254,237,55,204,181,114,251,243,60,16,22,20,239,231,183,37,68,139,247,238,250,113,241,28,221,177,249,140,45,28,182,15,245,90,44,
        81,211,186,108,17,46,9,147,96,95,175,219,180,207,176,171,228,146,21,191,130,57,230,149,252,62,53,4,147,43,186,62,33,37,25,116,12,20,8,100,194,115,249,179,186,3,219,7,86,15,40,236,45,191,65,170,207,1,194,162,7,7,29,4,
        28,48,77,119,68,42,237,100,157,15,18,95,26,216,170,242,112,202,213,192,180,21,43,164,173,194,159,237,84,188,246,129,235,176,91,224,69,226,134,116,140,118,187,137,56,142,16,35,11,219,22,40,203,217,110,65,134,167,3,243,235,205,106,249,
        68,244,20,193,12,163,183,242,85,216,103,68,252,218,65,32,15,95,120,26,204,253,4,142,46,194,40,136,197,78,227,151,203,13,49,79,211,120,175,2,148,173,153,201,64,29,62,60,130,114,218,201,241,205,251,159,71,88,42,74,116,177,109,233,
        194,232,243,71,218,8,139,130,7,202,137,85,114,105,137,232,126,173,178,158,56,219,67,34,95,41,210,31,226,240,38,190,138,72,120,244,146,80,163,38,225,212,196,46,16,193,211,203,49,202,123,248,228,106,238,150,128,104,58,121,94,178,32,89,
        199,214,232,133,109,242,2,205,175,224,62,146,223,73,234,205,191,17,119,209,176,231,222,193,236,177,3,92,79,181,125,212,204,98,63,58,65,241,138,65,35,202,50,160,162,7,122,38,66,88,27,17,50,88,220,53,190,131,210,110,236,58,199,0,
        6,223,92,83,5,245,27,223,248,106,156,140,93,59,117,242,169,132,249,95,125,226,240,100,147,52,155,95,206,186,185,0,74,102,213,15,89,167,97,3,87,144,217,208,80,207,61,20,137,72,230,83,114,191,245,75,111,26,43,38,152,61,5,142,
        34,226,96,100,125,180,164,125,164,178,232,73,192,207,174,25,54,2,0,6,73,22,70,18,194,52,100,161,156,154,156,17,11,140,155,129,37,74,207,36,57,179,79,116,193,77,187,244,3,84,185,26,144,80,132,31,10,36,36,141,86,171,37,198,
        204,25,215,151,118,111,220,29,253,90,96,80,244,17,25,211,147,147,103,32,72,66,10,44,14,250,125,210,150,179,84,18,68,198,249,9,102,186,179,242,246,47,23,148,48,199,189,203,253,49,78,171,177,200,205,156,236,84,102,236,90,119,155,81,
        36,181,11,167,77,152,21,122,151,156,110,168,189,209,188,138,202,85,107,248,220,125,79,69,176,220,109,193,237,109,122,20,107,121,125,3,100,104,130,24,248,84,3,36,7,104,231,130,108,11,141,117,34,164,197,252,130,213,221,162,206,200,73,219,
        3,199,144,146,217,104,31,121,29,158,236,239,199,74,16,148,142,151,0,92,130,252,57,56,210,46,3,174,33,205,178,236,208,16,29,123,228,21,25,218,54,244,110,81,150,231,127,100,50,17,36,242,15,119,80,45,223,166,101,75,184,107,141,175,
        178,82,139,151,164,220,131,95,56,118,10,140,231,209,26,160,38,140,30,139,232,231,13,181,50,103,40,155,51,31,131,69,53,68,127,255,136,77,178,143,39,73,79,80,2,244,43,59,249,183,132,57,80,98,97,233,20,219,181,226,113,145,221,100,
        170,43,9,51,83,15,177,253,248,140,217,174,50,176,189,57,204,36,10,123,99,153,75,9,248,194,20,132,8,247,211,173,202,165,193,6,149,235,204,80,165,13,90,249,65,202,100,90,16,208,51,129,27,43,214,95,200,0,165,200,56,137,99,218,
        216,216,53,19,229,213,98,77,129,180,198,191,244,244,33,89,213,186,121,36,190,58,10,144,241,106,136,28,157,62,96,243,10,210,244,212,224,250,145,226,63,168,19,159,34,162,249,175,57,22,192,17,120,186,204,20,145,125,129,28,172,170,116,84,
        180,103,32,175,0,6,146,117,18,215,223,126,195,114,24,65,57,233,226,219,140,36,238,107,41,126,131,123,154,24,157,72,64,7,10,64,16
    };
    uint8_t ss_ref[CRYPTO_BYTES]={
        197,114,75,90,23,55,251,197,171,63,74,206,215,24,6,58,228,23,221,37,48,160,107,176,240,135,139,129,192,102,61,26
    };
    uint8_t ss_res[CRYPTO_BYTES]={0};
    crypto_kem_dec_custom(ss_res,ct_ref,sk_ref);

    bool flag=true;
    for(int i=0;i<CRYPTO_BYTES;i++){
        if(ss_ref[i]!=ss_res[i]){
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_crypto_kem_dec_custom_level1");
    else printf("fail in test_crypto_kem_dec_custom_level1");
    return;
}

// int main(){
//     // test_safe_cmp_custom();
//     // test_functionH_custom();
//     test_functionL_custom();
//     // test_functionK_custom();
//     // test_compute_syndrome_custom();
//     // test_crypto_kem_keypair_custom();
//     // test_crypto_kem_keypair_custom_level1();
//     // test_crypto_kem_enc_custom_level1();
//     // test_crypto_kem_dec_custom_level1();
//     return 0;
// }