#include <string.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include "SABER_indcpa.h"
#include "poly.h"
#include "pack_unpack.h"
#include "poly_mul.h"
#include "fips202.h"
#include "SABER_params.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

#define h1 (1 << (SABER_EQ - SABER_EP - 1))
#define h2 ((1 << (SABER_EP - 2)) - (1 << (SABER_EP - SABER_ET - 1)) + (1 << (SABER_EQ - SABER_EP - 1)))

/***************************************************************************
 * 			        Using NTT Customized Versions
***************************************************************************/

void indcpa_kem_keypair_bf_custom(uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], uint8_t sk[SABER_INDCPA_SECRETKEYBYTES])
{
	csr_modulusq_rw(SABER_Q_EXT);
  	csr_qinv_rw(SABER_Q_EXT_INV);
	uint16_t A[SABER_L][SABER_L][SABER_N];
	uint16_t s[SABER_L][SABER_N];
	uint16_t b[SABER_L][SABER_N] = {0};

	uint8_t seed_A[SABER_SEEDBYTES];
	uint8_t seed_s[SABER_NOISE_SEEDBYTES];
	int i, j;

	srand(time(NULL));
	for (int i = 0; i < SABER_SEEDBYTES; i++){
		#if TEST==1
			seed_A[i] = i;
		#else
			seed_A[i] = rand() & 255;
		#endif
	}
	shake128_custom(seed_A, SABER_SEEDBYTES, seed_A, SABER_SEEDBYTES); // for not revealing system RNG state
	for (int i = 0; i < SABER_NOISE_SEEDBYTES; i++){
		#if TEST==1
			seed_s[i] = i;
		#else
			seed_s[i] = rand() & 255;
		#endif
	}

	GenMatrix_custom(A, seed_A);
	GenSecret_custom(s, seed_s);
	MatrixVectorMul_bf_custom(A, s, b, 1);
	for(i=0;i<SABER_L;i++){
        size_t vl;
        size_t avl=SABER_N;
        uint16_t* addr=b[i];
        vuint16m8_t va;
        while(avl>0){
            vl=vsetvl_e16m8(avl);
            va=vle16_v_u16m8(addr,vl);
            va=vadd_vx_u16m8(va,h1,vl);
            va=vsrl_vx_u16m8(va,SABER_EQ - SABER_EP,vl);
            vse16_v_u16m8(addr,va,vl);
            avl-=vl,addr+=vl;
        }
    }
	POLVECq2BS_custom(sk, s);
	POLVECp2BS_custom(pk, b);
	memcpy(pk + SABER_POLYVECCOMPRESSEDBYTES, seed_A, sizeof(seed_A));
}

void indcpa_kem_enc_bf_custom(const uint8_t m[SABER_KEYBYTES], const uint8_t seed_sp[SABER_NOISE_SEEDBYTES], const uint8_t pk[SABER_INDCPA_PUBLICKEYBYTES], uint8_t ciphertext[SABER_BYTES_CCA_DEC])
{
	csr_modulusq_rw(SABER_Q_EXT);
  	csr_qinv_rw(SABER_Q_EXT_INV);
	uint16_t A[SABER_L][SABER_L][SABER_N];
	uint16_t sp[SABER_L][SABER_N];
	uint16_t bp[SABER_L][SABER_N] = {0};
	uint16_t vp[SABER_N] = {0};
	uint16_t mp[SABER_N];
	uint32_t b[SABER_L][SABER_N];
	int i, j;
	const uint8_t *seed_A = pk + SABER_POLYVECCOMPRESSEDBYTES;

	GenMatrix_custom(A, seed_A);
	GenSecret_custom(sp, seed_sp);
	BS2POLVECp_custom(pk, b);
	MatrixVectorMul_bf_InnerPrd_custom(A, sp, bp, 0, b, vp);

	for (i = 0; i < SABER_L; i++)
	{
		size_t vl;
        size_t avl=SABER_N;
        uint16_t* addr=bp[i];
        vuint16m8_t va;
        while(avl>0){
            vl=vsetvl_e16m8(avl);
            va=vle16_v_u16m8(addr,vl);
            va=vadd_vx_u16m8(va,h1,vl);
            va=vsrl_vx_u16m8(va,SABER_EQ - SABER_EP,vl);
            vse16_v_u16m8(addr,va,vl);
            avl-=vl,addr+=vl;
        }
	}

	POLVECp2BS_custom(ciphertext, bp);

	BS2POLmsg_custom(m, mp);

	size_t vl;
    size_t avl=SABER_N;
    uint16_t* addr_vp=vp;
    uint16_t* addr_mp=mp;
    vuint16m8_t va,vb;
    while(avl>0){
        vl=vsetvl_e16m8(avl);
        va=vle16_v_u16m8(addr_mp,vl);
        vb=vle16_v_u16m8(addr_vp,vl);
        va=vsll_vx_u16m8(va,SABER_EP - 1,vl);
        vb=vadd_vx_u16m8(vb,h1,vl);
        vb=vsub_vv_u16m8(vb,va,vl);
        vb=vsrl_vx_u16m8(vb,SABER_EP - SABER_ET,vl);
        vse16_v_u16m8(addr_vp,vb,vl);
        addr_mp+=vl,addr_vp+=vl,avl-=vl;
    }

	POLT2BS_custom(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, vp);
}

void indcpa_kem_dec_bf_custom(const uint8_t sk[SABER_INDCPA_SECRETKEYBYTES], const uint8_t ciphertext[SABER_BYTES_CCA_DEC], uint8_t m[SABER_KEYBYTES])
{
	csr_modulusq_rw(SABER_Q_EXT);
  	csr_qinv_rw(SABER_Q_EXT_INV);
	uint16_t s[SABER_L][SABER_N];
    uint32_t s_tmp[SABER_L][SABER_N];
	uint32_t b[SABER_L][SABER_N];
	uint16_t v[SABER_N] = {0};
	uint16_t cm[SABER_N];
	int i;
    size_t vl;
    size_t avl;

	BS2POLVECq_custom(sk, s);
	for (int i = 0; i < SABER_L; i++) {
		vuint16m4_t va;
        vint16m4_t va_prime;
        vint32m8_t vb;
        vuint32m8_t vb_prime;
        vbool4_t vflag;
        //process s
        avl=SABER_N;
        uint16_t* src_addr=s[i];
        uint32_t* dst_addr=s_tmp[i];
        while(avl>0){
            vl=vsetvl_e16m4(avl);
            va=vle16_v_u16m4(src_addr,vl);
            va_prime=vreinterpret_v_u16m4_i16m4(va);
            vflag=vmsgt_vx_i16m4_b4(va_prime,8192 >> 1,vl);
            va_prime=vsub_vx_i16m4_m(vflag,va_prime,va_prime,8192,vl);
            vb=vsext_vf2_i32m8(va_prime,vl);
            vb=vadd_vx_i32m8_m(vflag,vb,vb,SABER_Q_EXT,vl);
            vb_prime=vreinterpret_v_i32m8_u32m8(vb);
            vse32_v_u32m8(dst_addr,vb_prime,vl);
            src_addr+=vl,dst_addr+=vl,avl-=vl;
        }
	}

	BS2POLVECp_custom(ciphertext, b);
	InnerProd_bf_custom(b, s_tmp, v);
	BS2POLT_custom(ciphertext + SABER_POLYVECCOMPRESSEDBYTES, cm);

	for (i = 0; i < SABER_N; i++)
	{
		v[i] = (v[i] + h2 - (cm[i] << (SABER_EP - SABER_ET))) >> (SABER_EP - 1);
	}

	POLmsg2BS_custom(m, v);
}