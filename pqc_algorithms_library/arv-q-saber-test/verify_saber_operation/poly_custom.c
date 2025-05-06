#include <stdio.h>
#include "api.h"
#include "poly.h"
#include "poly_mul.h"
#include "pack_unpack.h"
#include "cbd.h"
#include "fips202.h"
#include "ntt.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

#define USE_NTT_ASM

void MatrixVectorMul_bf_custom(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose)
{
	uint32_t temp_res[SABER_L][SABER_N]={0};
	uint32_t temp_A[SABER_L][SABER_L][SABER_N];
	uint32_t temp_s[SABER_L][SABER_N];
    size_t vl;
    size_t avl;
	for (int i = 0; i < SABER_L; i++) {
		for (int j = 0; j < SABER_L; j++) {
            vuint16m4_t va;
            vuint32m8_t vb;
			avl=SABER_N;
            uint16_t* src_addr=A[i][j];
            uint32_t* dst_addr=temp_A[i][j];
            while(avl>0){
                vl=vsetvl_e16m4(avl);
                va=vle16_v_u16m4(src_addr,vl);
                vb=vwaddu_vx_u32m8(va,0,vl);
                vse32_v_u32m8(dst_addr,vb,vl);
                avl-=vl,src_addr+=vl,dst_addr+=vl;
            }
			#ifdef USE_NTT_ASM
            ntt_cg_custom_reorder_asm(temp_A[i][j]);
            #else
            ntt_cg_custom_reorder(temp_A[i][j]);
            #endif
		}
	}
	for (int i = 0; i < SABER_L; i++) {
		vuint16m4_t va;
        vint16m4_t va_prime;
        vint32m8_t vb;
        vuint32m8_t vb_prime;
        vbool4_t vflag;
        avl=SABER_N;
        uint16_t* src_addr=s[i];
        uint32_t* dst_addr=temp_s[i];
        while(avl>0){
            vl=vsetvl_e16m4(avl);
            va=vle16_v_u16m4(src_addr,vl);
            va_prime=vreinterpret_v_u16m4_i16m4(va);
            vb=vsext_vf2_i32m8(va_prime,vl);
            vflag=vmslt_vx_i32m8_b4(vb,0,vl);
            vb=vadd_vx_i32m8_m(vflag,vb,vb,SABER_Q_EXT,vl);
            vb_prime=vreinterpret_v_i32m8_u32m8(vb);
            vse32_v_u32m8(dst_addr,vb_prime,vl);
            avl-=vl,src_addr+=vl,dst_addr+=vl;
        }
		#ifdef USE_NTT_ASM
        ntt_cg_custom_reorder_asm(temp_s[i]);
        #else
        ntt_cg_custom_reorder(temp_s[i]);
        #endif
	}
	int i, j;
	for (i = 0; i < SABER_L; i++)
	{
		for (j = 0; j < SABER_L; j++)
		{
			if (transpose == 1)
			{
				poly_mul_acc_bf_custom(temp_A[j][i], temp_s[j], temp_res[i]);
			}
			else
			{
				poly_mul_acc_bf_custom(temp_A[i][j], temp_s[j], temp_res[i]);
			}
		}
	}

	for (int i = 0; i < SABER_L; i++) {
		#ifdef USE_NTT_ASM
        invntt_cg_custom_reorder_asm(temp_res[i]);
        #else
        invntt_cg_custom_reorder(temp_res[i]);
        #endif
        //mon2nor and mod-down
        vuint32m8_t va;
        vint32m8_t va_prime;
        vuint16m4_t vb;
        vint16m4_t vb_prime;
        vbool4_t vflag;
        avl=SABER_N;
        uint32_t* src_addr=temp_res[i];
        uint16_t* dst_addr=res[i];
        while(avl>0){
            vl=vsetvl_e32m8(avl);
            va=vle32_v_u32m8(src_addr,vl);
            va=vmod_mul_vx_u32m8(va,22495266);
            va_prime=vreinterpret_v_u32m8_i32m8(va);
            vflag=vmsgt_vx_i32m8_b4(va_prime,25171457>>1,vl);
            va_prime=vsub_vx_i32m8_m(vflag,va_prime,va_prime,25171457,vl);
            vb_prime=vnsra_wx_i16m4(va_prime,0,vl);
            vb_prime=vand_vx_i16m4(vb_prime,8191,vl);
            vb=vreinterpret_v_i16m4_u16m4(vb_prime);
            vse16_v_u16m4(dst_addr,vb,vl);
            avl-=vl,src_addr+=vl,dst_addr+=vl;
        }
	}
}

void MatrixVectorMul_bf_InnerPrd_custom(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose, 
    uint32_t b[SABER_L][SABER_N], uint16_t v[SABER_N])
{
    //begin MatrixVectorMul_bf part
	uint32_t temp_res[SABER_L][SABER_N]={0};
	uint32_t temp_A[SABER_L][SABER_L][SABER_N];
	uint32_t temp_s[SABER_L][SABER_N];
    size_t vl;
    size_t avl;
	for (int i = 0; i < SABER_L; i++) {
		for (int j = 0; j < SABER_L; j++) {
            vuint16m4_t va;
            vuint32m8_t vb;
			avl=SABER_N;
            uint16_t* src_addr=A[i][j];
            uint32_t* dst_addr=temp_A[i][j];
            while(avl>0){
                vl=vsetvl_e16m4(avl);
                va=vle16_v_u16m4(src_addr,vl);
                vb=vwaddu_vx_u32m8(va,0,vl);
                vse32_v_u32m8(dst_addr,vb,vl);
                avl-=vl,src_addr+=vl,dst_addr+=vl;
            }
			#ifdef USE_NTT_ASM
            ntt_cg_custom_reorder_asm(temp_A[i][j]);
            #else
            ntt_cg_custom_reorder(temp_A[i][j]);
            #endif
		}
	}
	for (int i = 0; i < SABER_L; i++) {
		vuint16m4_t va;
        vint16m4_t va_prime;
        vint32m8_t vb;
        vuint32m8_t vb_prime;
        vbool4_t vflag;
        avl=SABER_N;
        uint16_t* src_addr=s[i];
        uint32_t* dst_addr=temp_s[i];
        while(avl>0){
            vl=vsetvl_e16m4(avl);
            va=vle16_v_u16m4(src_addr,vl);
            va_prime=vreinterpret_v_u16m4_i16m4(va);
            vb=vsext_vf2_i32m8(va_prime,vl);
            vflag=vmslt_vx_i32m8_b4(vb,0,vl);
            vb=vadd_vx_i32m8_m(vflag,vb,vb,SABER_Q_EXT,vl);
            vb_prime=vreinterpret_v_i32m8_u32m8(vb);
            vse32_v_u32m8(dst_addr,vb_prime,vl);
            avl-=vl,src_addr+=vl,dst_addr+=vl;
        }
		#ifdef USE_NTT_ASM
        ntt_cg_custom_reorder_asm(temp_s[i]);
        #else
        ntt_cg_custom_reorder(temp_s[i]);
        #endif
	}
	int i, j;
	for (i = 0; i < SABER_L; i++)
	{
		for (j = 0; j < SABER_L; j++)
		{
			if (transpose == 1)
			{
				poly_mul_acc_bf_custom(temp_A[j][i], temp_s[j], temp_res[i]);
			}
			else
			{
				poly_mul_acc_bf_custom(temp_A[i][j], temp_s[j], temp_res[i]);
			}
		}
	}

	for (int i = 0; i < SABER_L; i++) {
		#ifdef USE_NTT_ASM
        invntt_cg_custom_reorder_asm(temp_res[i]);
        #else
        invntt_cg_custom_reorder(temp_res[i]);
        #endif
        //mon2nor and mod-down
        vuint32m8_t va;
        vint32m8_t va_prime;
        vuint16m4_t vb;
        vint16m4_t vb_prime;
        vbool4_t vflag;
        avl=SABER_N;
        uint32_t* src_addr=temp_res[i];
        uint16_t* dst_addr=res[i];
        while(avl>0){
            vl=vsetvl_e32m8(avl);
            va=vle32_v_u32m8(src_addr,vl);
            va=vmod_mul_vx_u32m8(va,22495266);
            va_prime=vreinterpret_v_u32m8_i32m8(va);
            vflag=vmsgt_vx_i32m8_b4(va_prime,25171457>>1,vl);
            va_prime=vsub_vx_i32m8_m(vflag,va_prime,va_prime,25171457,vl);
            vb_prime=vnsra_wx_i16m4(va_prime,0,vl);
            vb_prime=vand_vx_i16m4(vb_prime,8191,vl);
            vb=vreinterpret_v_i16m4_u16m4(vb_prime);
            vse16_v_u16m4(dst_addr,vb,vl);
            avl-=vl,src_addr+=vl,dst_addr+=vl;
        }
	}

    //begin InnerProd Part
    uint32_t temp_v[SABER_N]={0};
    for(int i = 0; i < SABER_L; i++){
        #ifdef USE_NTT_ASM
        ntt_cg_custom_reorder_asm(b[i]);
        #else
        ntt_cg_custom_reorder(b[i]);
        #endif
        poly_mul_acc_bf_custom(b[i],temp_s[i],temp_v);
    }
    #ifdef USE_NTT_ASM
    invntt_cg_custom_reorder_asm(temp_v);
    #else
    invntt_cg_custom_reorder(temp_v);
    #endif
    vuint32m8_t va;
    vint32m8_t va_prime;
    vuint16m4_t vb;
    vint16m4_t vb_prime;
    vbool4_t vflag;
    avl=SABER_N;
    uint32_t* src_addr=temp_v;
    uint16_t* dst_addr=v;
    while(avl>0){
        vl=vsetvl_e32m8(avl);
        va=vle32_v_u32m8(src_addr,vl);
        va=vmod_mul_vx_u32m8(va,22495266);
        va_prime=vreinterpret_v_u32m8_i32m8(va);
        vflag=vmsgt_vx_i32m8_b4(va_prime,25171457>>1,vl);
        va_prime=vsub_vx_i32m8_m(vflag,va_prime,va_prime,25171457,vl);
        vb_prime=vnsra_wx_i16m4(va_prime,0,vl);
        vb_prime=vand_vx_i16m4(vb_prime,8191,vl);
        vb=vreinterpret_v_i16m4_u16m4(vb_prime);
        vse16_v_u16m4(dst_addr,vb,vl);
        avl-=vl,src_addr+=vl,dst_addr+=vl;
    }
}

void InnerProd_bf_custom(uint32_t b[SABER_L][SABER_N], uint32_t s[SABER_L][SABER_N], uint16_t res[SABER_N]) {
	uint32_t temp_s[SABER_L][SABER_N];
    uint32_t temp_res[SABER_N]={0};
    size_t vl;
    size_t avl;

    for(int i=0;i<SABER_L;i++){
        #ifdef USE_NTT_ASM
        ntt_cg_custom_reorder_asm(s[i]);
        //process b
        ntt_cg_custom_reorder_asm(b[i]);
        #else
        ntt_cg_custom_reorder(s[i]);
        //process b
        ntt_cg_custom_reorder(b[i]);
        #endif
        //mul and acc
        poly_mul_acc_bf_custom(b[i],s[i],temp_res);
    }
    #ifdef USE_NTT_ASM
    invntt_cg_custom_reorder_asm(temp_res);
    #else
    invntt_cg_custom_reorder(temp_res);
    #endif
    vuint32m8_t va;
    vint32m8_t va_prime;
    vuint16m4_t vb;
    vint16m4_t vb_prime;
    vbool4_t vflag;
    avl=SABER_N;
    uint32_t* src_addr=temp_res;
    uint16_t* dst_addr=res;
    while(avl>0){
        vl=vsetvl_e32m8(avl);
        va=vle32_v_u32m8(src_addr,vl);
        va=vmod_mul_vx_u32m8(va,22495266);
        va_prime=vreinterpret_v_u32m8_i32m8(va);
        vflag=vmsgt_vx_i32m8_b4(va_prime,25171457>>1,vl);
        va_prime=vsub_vx_i32m8_m(vflag,va_prime,va_prime,25171457,vl);
        vb_prime=vnsra_wx_i16m4(va_prime,0,vl);
        vb_prime=vand_vx_i16m4(vb_prime,8191,vl);
        vb=vreinterpret_v_i16m4_u16m4(vb_prime);
        vse16_v_u16m4(dst_addr,vb,vl);
        avl-=vl,src_addr+=vl,dst_addr+=vl;
    }
}

void GenMatrix_custom(uint16_t A[SABER_L][SABER_L][SABER_N], const uint8_t seed[SABER_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYVECBYTES];
	int i;

	shake128_custom(buf, sizeof(buf), seed, SABER_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		BS2POLVECq_custom(buf + i * SABER_POLYVECBYTES, A[i]);
	}
}

void GenSecret_custom(uint16_t s[SABER_L][SABER_N], const uint8_t seed[SABER_NOISE_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYCOINBYTES];
	size_t i;

	shake128_custom(buf, sizeof(buf), seed, SABER_NOISE_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		cbd_custom(s[i], buf + i * SABER_POLYCOINBYTES);
	}
}