#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include "SABER_params.h"
#include "reduce.h"
#include "poly.h"
#include "ntt.h"

#if (VLEN==256)
/*************************************************
* Name:        ntt_cg_custom_reorder
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
*
* Requirements: -VLEN=256
**************************************************/
void ntt_cg_custom_reorder(int32_t r[SABER_N]){
    int32_t r_temp0[256];
    unsigned int i, j;
    int32_t* r_ptr0 = NULL;
    int32_t* r_ptr1 = NULL;
    int32_t* r_ptr2 = NULL;

    vint32m1_t v_zetas;
    //for stage 1 to stage 4
    vint32m2_t v_in0,v_in1;
    vint32m4_t v_out0;
    //for stage 5 to stage 8
    vint32m1_t v_in2,v_in3;
    vint32m2_t v_out1;

    //stage1, repeat bound=1
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
    v_zetas=vle32_v_i32m1(&zetas[1],1);
    for(i=0;i<8;i++){
        v_in0=vle32_v_i32m2(r_ptr0,16);
        v_in1=vle32_v_i32m2(r_ptr1,16);
        r_ptr0+=16;
        r_ptr1+=16;
        v_out0=vbutterfly_ct_vvm_i32m4(v_in0,v_in1,v_zetas);
        vse32_v_i32m4(r_ptr2,v_out0,32);
        r_ptr2+=32;
    }

    //stage2, repeat bound=2
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));
    v_zetas=vle32_v_i32m1(&zetas[2],2);
    for(i=0;i<8;i++){
        v_in0=vle32_v_i32m2(r_ptr0,16);
        v_in1=vle32_v_i32m2(r_ptr1,16);
        r_ptr0+=16;
        r_ptr1+=16;
        v_out0=vbutterfly_ct_vvm_i32m4(v_in0,v_in1,v_zetas);
        vse32_v_i32m4(r_ptr2,v_out0,32);
        r_ptr2+=32;
    }

    //stage3, repeat bound=4
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
    v_zetas=vle32_v_i32m1(&zetas[4],4);
    for(i=0;i<8;i++){
        v_in0=vle32_v_i32m2(r_ptr0,16);
        v_in1=vle32_v_i32m2(r_ptr1,16);
        r_ptr0+=16;
        r_ptr1+=16;
        v_out0=vbutterfly_ct_vvm_i32m4(v_in0,v_in1,v_zetas);
        vse32_v_i32m4(r_ptr2,v_out0,32);
        r_ptr2+=32;
    }

    //stage4, repeat bound=8
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
    v_zetas=vle32_v_i32m1(&zetas[8],8);
    for(i=0;i<8;i++){
        v_in0=vle32_v_i32m2(r_ptr0,16);
        v_in1=vle32_v_i32m2(r_ptr1,16);
        r_ptr0+=16;
        r_ptr1+=16;
        v_out0=vbutterfly_ct_vvm_i32m4(v_in0,v_in1,v_zetas);
        vse32_v_i32m4(r_ptr2,v_out0,32);
        r_ptr2+=32;
    }

    //stage5, actual repeat bound=16, exceed v0's capacity, so set repeat bound to 8
    for(i=0;i<2;i++){
        int offset=8*i;
        r_ptr0=r+offset;
        r_ptr1=&r[128+offset];
        r_ptr2=r_temp0+(offset<<1);
        v_zetas=vle32_v_i32m1(&zetas[16+offset],8);
        for(j=0;j<8;j++){
            v_in2=vle32_v_i32m1(r_ptr0,8);
            v_in3=vle32_v_i32m1(r_ptr1,8);
            r_ptr0+=16;
            r_ptr1+=16;
            v_out1=vbutterfly_ct_vvm_i32m2(v_in2,v_in3,v_zetas);
            vse32_v_i32m2(r_ptr2,v_out1,16);
            r_ptr2+=32;
        }
    }

    //stage6, actual repeat bound=32, exceed v0's capacity, so set repeat bound to 8
    for(i=0;i<4;i++){
        int offset=8*i;
        r_ptr0=r_temp0+offset;
        r_ptr1=&r_temp0[128+offset];
        r_ptr2=r+(offset<<1);
        v_zetas=vle32_v_i32m1(&zetas[32+offset],8);
        for(j=0;j<4;j++){
            v_in2=vle32_v_i32m1(r_ptr0,8);
            v_in3=vle32_v_i32m1(r_ptr1,8);
            r_ptr0+=32;
            r_ptr1+=32;
            v_out1=vbutterfly_ct_vvm_i32m2(v_in2,v_in3,v_zetas);
            vse32_v_i32m2(r_ptr2,v_out1,16);
            r_ptr2+=64;
        }
    }

    //stage7, actual repeat bound=64, exceed v0's capacity, so set repeat bound to 8
    for(i=0;i<8;i++){
        int offset=8*i;
        r_ptr0=r+offset;
        r_ptr1=&r[128+offset];
        r_ptr2=r_temp0+(offset<<1);
        v_zetas=vle32_v_i32m1(&zetas[64+offset],8);
        for(j=0;j<2;j++){
            v_in2=vle32_v_i32m1(r_ptr0,8);
            v_in3=vle32_v_i32m1(r_ptr1,8);
            r_ptr0+=64;
            r_ptr1+=64;
            v_out1=vbutterfly_ct_vvm_i32m2(v_in2,v_in3,v_zetas);
            vse32_v_i32m2(r_ptr2,v_out1,16);
            r_ptr2+=128;
        }
    }

    vuint32m2_t v_offset;
    uint32_t* tree_ptr=NULL;
    //stage8, actual repeat bound=128, exceed v0's capacity, so set repeat bound to 8
    for(i=0;i<16;i++){
        int offset=8*i;
        r_ptr0=r_temp0+offset;
        r_ptr1=&r_temp0[128+offset];
        tree_ptr=tree_byteoffset+(offset<<1);
        v_offset=vle32_v_u32m2(tree_ptr,16);
        v_zetas=vle32_v_i32m1(&zetas[128+offset],8);
        v_in2=vle32_v_i32m1(r_ptr0,8);
        v_in3=vle32_v_i32m1(r_ptr1,8);
        v_out1=vbutterfly_ct_vvm_i32m2(v_in2,v_in3,v_zetas);
        vsuxei32_v_i32m2(r,v_offset,v_out1,16);
    }
}

/*************************************************
* Name:        invntt_cg_custom_reorder
*
* Description: Inverse NTT.
*              Out-of-place, constant geometry.
*              Input coefficients are in standard order and normal domain.
*              Output coefficient are smaller than Q in absolute value;
*              and output coefficients are in standard order and normal domain.
*              zeta in compact form, that zetas in v0 are not the same
*
* Arguments:   - int32_t p[N]: input/output coefficient array
*
* Requirements: -VLEN=256
**************************************************/
void invntt_cg_custom_reorder(int32_t r[SABER_N]){
    int32_t r_temp0[256];
    unsigned int i;
    int32_t* r_ptr0 = NULL;
    int32_t* r_ptr1 = NULL;
    int32_t* r_ptr2 = NULL;
    int32_t* zeta_ptr = inv_zetas_inorder;

    vint32m1_t v_zetas;
    //for stage1 and stage2
    vint32m1_t v_in0,v_in1;
    vint32m2_t v_out0;
    //for stage3 to stage8
    vint32m4_t v_in2,v_in3;
    vint32m8_t v_out1;

    //stage1,same num=1
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 6));//repeat time here should be set large enough
    for(i=0;i<16;i++){//totally use 128zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,8);//8bf use 8zetas
        zeta_ptr+=8;
        v_in0=vle32_v_i32m1(r_ptr0,8);
        v_in1=vle32_v_i32m1(r_ptr1,8);
        r_ptr0+=8;
        r_ptr1+=8;
        v_out0=vbutterfly_gs_vvm_i32m2(v_in0,v_in1,v_zetas);
        vse32_v_i32m2(r_ptr2,v_out0,16);
        r_ptr2+=16;
    }

    //stage2, same num=2
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(1, 6));
    for(i=0;i<16;i++){//totally use 64zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,4);//8bf use 4zetas
        zeta_ptr+=4;
        v_in0=vle32_v_i32m1(r_ptr0,8);
        v_in1=vle32_v_i32m1(r_ptr1,8);
        r_ptr0+=8;
        r_ptr1+=8;
        v_out0=vbutterfly_gs_vvm_i32m2(v_in0,v_in1,v_zetas);
        vse32_v_i32m2(r_ptr2,v_out0,16);
        r_ptr2+=16;
    }

    //stage3, same num=4
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(2, 6));
    for(i=0;i<4;i++){//totally use 32zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,8);//32bf use 8zetas
        zeta_ptr+=8;
        v_in2=vle32_v_i32m4(r_ptr0,32);
        v_in3=vle32_v_i32m4(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        vsetvl_e32m4_wrapper(32);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,64);
        r_ptr2+=64;
    }

    //stage4, same num=8
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(3, 6));
    for(i=0;i<4;i++){//totally use 16zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,4);//32bf use 4zetas
        zeta_ptr+=4;
        v_in2=vle32_v_i32m4(r_ptr0,32);
        v_in3=vle32_v_i32m4(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        vsetvl_e32m4_wrapper(32);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,64);
        r_ptr2+=64;
    }

    //stage5, same num=16
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(4, 6));
    for(i=0;i<4;i++){//totally use 8zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,2);//32bf use 2zetas
        zeta_ptr+=2;
        v_in2=vle32_v_i32m4(r_ptr0,32);
        v_in3=vle32_v_i32m4(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        vsetvl_e32m4_wrapper(32);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,64);
        r_ptr2+=64;
    }

    //stage6, same num=32
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(5, 6));
    for(i=0;i<4;i++){//totally use 4zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,1);//32bf use 1zetas
        zeta_ptr+=1;
        v_in2=vle32_v_i32m4(r_ptr0,32);
        v_in3=vle32_v_i32m4(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        vsetvl_e32m4_wrapper(32);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,64);
        r_ptr2+=64;
    }

    //stage7, same num=64
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(6, 6));
    for(i=0;i<2;i++){//totally use 2zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,1);
        zeta_ptr+=1;

        v_in2=vle32_v_i32m4(r_ptr0,32);
        v_in3=vle32_v_i32m4(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        vsetvl_e32m4_wrapper(32);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,64);
        r_ptr2+=64;

        v_in2=vle32_v_i32m4(r_ptr0,32);
        v_in3=vle32_v_i32m4(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        vsetvl_e32m4_wrapper(32);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,64);
        r_ptr2+=64;
    }

    vuint32m8_t v_offset;
    uint32_t* tree_ptr=tree_byteoffset;

    //stage8, same num=128, all use the same zeta
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    csr_zeta_selMode_rw(GET_ZETA_MODE(7, 6));
    v_zetas=vle32_v_i32m1(zeta_ptr,1);
    for(i=0;i<4;i++){
        v_in2=vle32_v_i32m4(r_ptr0,32);
        v_in3=vle32_v_i32m4(r_ptr1,32);
        v_offset=vle32_v_u32m8(tree_ptr,64);
        tree_ptr+=64;
        r_ptr0+=32;
        r_ptr1+=32;
        vsetvl_e32m4_wrapper(32);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vsetvl_e32m8_wrapper(64);
        vsuxei32_v_i32m8(r,v_offset,v_out1,64);
    }
}

#elif (VLEN==512)
/*************************************************
* Name:        ntt_cg_custom_reorder
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
*
* Requirements: -VLEN=512
**************************************************/
void ntt_cg_custom_reorder(int32_t r[SABER_N]){
    int32_t r_temp0[256];
    unsigned int i, j;
    int32_t* r_ptr0 = NULL;
    int32_t* r_ptr1 = NULL;
    int32_t* r_ptr2 = NULL;

    vint32m1_t v_zetas;
    //for stage 1 to stage 5
    vint32m2_t v_in0,v_in1;
    vint32m4_t v_out0;
    //for stage 6 to stage 8
    vint32m1_t v_in2,v_in3;
    vint32m2_t v_out1;

    //stage1, repeat bound=1
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
    v_zetas=vle32_v_i32m1(&zetas[1],1);
    for(i=0;i<4;i++){
        v_in0=vle32_v_i32m2(r_ptr0,32);
        v_in1=vle32_v_i32m2(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        v_out0=vbutterfly_ct_vvm_i32m4(v_in0,v_in1,v_zetas);
        vse32_v_i32m4(r_ptr2,v_out0,64);
        r_ptr2+=64;
    }

    //stage2, repeat bound=2
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));
    v_zetas=vle32_v_i32m1(&zetas[2],2);
    for(i=0;i<4;i++){
        v_in0=vle32_v_i32m2(r_ptr0,32);
        v_in1=vle32_v_i32m2(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        v_out0=vbutterfly_ct_vvm_i32m4(v_in0,v_in1,v_zetas);
        vse32_v_i32m4(r_ptr2,v_out0,64);
        r_ptr2+=64;
    }

    //stage3, repeat bound=4
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
    v_zetas=vle32_v_i32m1(&zetas[4],4);
    for(i=0;i<4;i++){
        v_in0=vle32_v_i32m2(r_ptr0,32);
        v_in1=vle32_v_i32m2(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        v_out0=vbutterfly_ct_vvm_i32m4(v_in0,v_in1,v_zetas);
        vse32_v_i32m4(r_ptr2,v_out0,64);
        r_ptr2+=64;
    }

    //stage4, repeat bound=8
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
    v_zetas=vle32_v_i32m1(&zetas[8],8);
    for(i=0;i<4;i++){
        v_in0=vle32_v_i32m2(r_ptr0,32);
        v_in1=vle32_v_i32m2(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        v_out0=vbutterfly_ct_vvm_i32m4(v_in0,v_in1,v_zetas);
        vse32_v_i32m4(r_ptr2,v_out0,64);
        r_ptr2+=64;
    }

    //stage5, repeat bound=16
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
    v_zetas=vle32_v_i32m1(&zetas[16],16);
    for(i=0;i<4;i++){
        v_in0=vle32_v_i32m2(r_ptr0,32);
        v_in1=vle32_v_i32m2(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        v_out0=vbutterfly_ct_vvm_i32m4(v_in0,v_in1,v_zetas);
        vse32_v_i32m4(r_ptr2,v_out0,64);
        r_ptr2+=64;
    }

    //stage6, actual repeat bound=32, exceed v0's capacity, so set repeat bound to 16
    for(i=0;i<2;i++){
        int offset=16*i;
        r_ptr0=r_temp0+offset;
        r_ptr1=&r_temp0[128+offset];
        r_ptr2=r+(offset<<1);
        v_zetas=vle32_v_i32m1(&zetas[32+offset],16);
        for(j=0;j<4;j++){
            v_in2=vle32_v_i32m1(r_ptr0,16);
            v_in3=vle32_v_i32m1(r_ptr1,16);
            r_ptr0+=32;
            r_ptr1+=32;
            v_out1=vbutterfly_ct_vvm_i32m2(v_in2,v_in3,v_zetas);
            vse32_v_i32m2(r_ptr2,v_out1,32);
            r_ptr2+=64;
        }
    }

    //stage7, actual repeat bound=64, exceed v0's capacity, so set repeat bound to 16
    for(i=0;i<4;i++){
        int offset=16*i;
        r_ptr0=r+offset;
        r_ptr1=&r[128+offset];
        r_ptr2=r_temp0+(offset<<1);
        v_zetas=vle32_v_i32m1(&zetas[64+offset],16);
        for(j=0;j<2;j++){
            v_in2=vle32_v_i32m1(r_ptr0,16);
            v_in3=vle32_v_i32m1(r_ptr1,16);
            r_ptr0+=64;
            r_ptr1+=64;
            v_out1=vbutterfly_ct_vvm_i32m2(v_in2,v_in3,v_zetas);
            vse32_v_i32m2(r_ptr2,v_out1,32);
            r_ptr2+=128;
        }
    }

    vuint32m2_t v_offset;
    uint32_t* tree_ptr=NULL;
    //stage8, actual repeat bound=128, exceed v0's capacity, so set repeat bound to 16
    for(i=0;i<8;i++){
        int offset=16*i;
        r_ptr0=r_temp0+offset;
        r_ptr1=&r_temp0[128+offset];
        tree_ptr=tree_byteoffset+(offset<<1);
        v_offset=vle32_v_u32m2(tree_ptr,32);
        v_zetas=vle32_v_i32m1(&zetas[128+offset],16);
        v_in2=vle32_v_i32m1(r_ptr0,16);
        v_in3=vle32_v_i32m1(r_ptr1,16);
        v_out1=vbutterfly_ct_vvm_i32m2(v_in2,v_in3,v_zetas);
        vsuxei32_v_i32m2(r,v_offset,v_out1,32);
    }
}

/*************************************************
* Name:        ntt_cg_custom_reorder_asm
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
*
* Requirements: -VLEN=512
**************************************************/
void ntt_cg_custom_reorder_asm(int32_t r[SABER_N]){
    size_t load_zeta_vl;
    size_t avl=32;//m2

    //stage1,repeat bound=1
    //v4 and v6 with m2 stores input coefficients
    //v16, v20, v24, v28 with m4 store output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
    load_zeta_vl=1;
    __asm__ __volatile__ ( 
    "vsetvli zero, %[zeta_vl], e32, m1, ta, mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli zero, %[avl], e32, m2, ta, mu\n"
    "vle32.v v4, (%[top_base_0])\n"
    "vle32.v v6, (%[btm_base_0])\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    "vle32.v v4, (%[top_base_1])\n"
    "vle32.v v6, (%[btm_base_1])\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    "vle32.v v4, (%[top_base_2])\n"
    "vle32.v v6, (%[btm_base_2])\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    "vle32.v v4, (%[top_base_3])\n"
    "vle32.v v6, (%[btm_base_3])\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base]"r"(&zetas[1]), [avl]"r"(avl),
     [top_base_0]"r"(&r[0]), [btm_base_0]"r"(&r[128]), [top_base_1]"r"(&r[32]), [btm_base_1]"r"(&r[160]),
     [top_base_2]"r"(&r[64]), [btm_base_2]"r"(&r[192]), [top_base_3]"r"(&r[96]), [btm_base_3]"r"(&r[224])
  );

    //stage2, repeat bound=2
    // (v16, v24), (v18, v26), (v20, v28), (v22, v30) with m2 are input coefficients
    // v4, v8, v12, v16 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));
    load_zeta_vl=2;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.ct.vvm v4, v16, v24, v0\n"
    "vbutterfly.ct.vvm v8, v18, v26, v0\n"
    "vbutterfly.ct.vvm v12, v20, v28, v0\n"
    "vbutterfly.ct.vvm v16, v22, v30, v0\n"
    :
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas[2]), [avl]"r"(avl)
  );

    //stage3, repeat bound=4
    // (v4, v12), (v6, v14), (v8, v16), (v10, v18) with m2 are input coefficients
    // v20, v24, v28, v4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
    load_zeta_vl=4;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.ct.vvm v20, v4, v12, v0\n"
    "vbutterfly.ct.vvm v24, v6, v14, v0\n"
    "vbutterfly.ct.vvm v28, v8, v16, v0\n"
    "vbutterfly.ct.vvm v4, v10, v18, v0\n"
    :
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas[4]), [avl]"r"(avl)
  );

    //stage4, repeat bound=8
    // (v20, v28), (v22, v30), (v24, v4), (v26, v6) with m2 are input coefficients
    // v8, v12, v16, v20 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
    load_zeta_vl=8;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.ct.vvm v8, v20, v28, v0\n"
    "vbutterfly.ct.vvm v12, v22, v30, v0\n"
    "vbutterfly.ct.vvm v16, v24, v4, v0\n"
    "vbutterfly.ct.vvm v20, v26, v6, v0\n"
    :
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas[8]), [avl]"r"(avl)
  );

    //stage5, repeat bound=16
    // (v8, v16), (v10, v18), (v12, v20), (v14, v22) with m2 are input coefficients
    // v24, v28, v4, v8 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
    load_zeta_vl=16;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.ct.vvm v24, v8, v16, v0\n"
    "vbutterfly.ct.vvm v28, v10, v18, v0\n"
    "vbutterfly.ct.vvm v4, v12, v20, v0\n"
    "vbutterfly.ct.vvm v8, v14, v22, v0\n"
    :
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas[16]), [avl]"r"(avl)
  );

    //stage6, repeat bound=32
    // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m1 are input coefficients, with same zetas
    // (v25, v5), (v27, v7), (v29, v9), (v31, v11) with m1 are input coefficients, with same zetas
    // v12, v14, v16, v18, v20, v22, v2, v4 with m2 are output coefficients
    avl=16;//m1
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v12, v24, v4, v0\n"
    "vbutterfly.ct.vvm v16, v26, v6, v0\n"
    "vbutterfly.ct.vvm v20, v28, v8, v0\n"
    "vbutterfly.ct.vvm v2, v30, v10, v0\n"
    "vle32.v v0, (%[zeta_base1])\n"
    "vbutterfly.ct.vvm v14, v25, v5, v0\n"
    "vbutterfly.ct.vvm v18, v27, v7, v0\n"
    "vbutterfly.ct.vvm v22, v29, v9, v0\n"
    "vbutterfly.ct.vvm v4, v31, v11, v0\n"
    :
    : [zeta_vl] "r"(avl), [zeta_base0] "r"(&zetas[32]), [zeta_base1] "r"(&zetas[48])
  );

    //stage7, repeat bound=64
    // (v12, v20), (v13, v21), (v14, v22), (v15, v23) with m1 are input coefficients
    // (v16, v2), (v17, v3), (v18, v4), (v19, v5) with m1 are input coefficients
    // ((v12, v20), (v16, v2)) same zetas, ((v13, v21), (v17, v3)) same zetas
    // ((v14, v22), (v18, v4)) same zetas, ((v15, v23), (v19, v5)) same zetas
    // v24, v26, v28, v30, v6, v8, v10, v12 with m2 is output coefficients
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v24, v12, v20, v0\n"
    "vbutterfly.ct.vvm v6, v16, v2, v0\n"
    "vle32.v v0, (%[zeta_base1])\n" 
    "vbutterfly.ct.vvm v26, v13, v21, v0\n"
    "vbutterfly.ct.vvm v8, v17, v3, v0\n"
    "vle32.v v0, (%[zeta_base2])\n"
    "vbutterfly.ct.vvm v28, v14, v22, v0\n"
    "vbutterfly.ct.vvm v10, v18, v4, v0\n"
    "vle32.v v0, (%[zeta_base3])\n"
    "vbutterfly.ct.vvm v30, v15, v23, v0\n"
    "vbutterfly.ct.vvm v12, v19, v5, v0\n"
    :
    : [zeta_vl] "r"(avl), [zeta_base0] "r"(&zetas[64]), [zeta_base1] "r"(&zetas[80]), [zeta_base2] "r"(&zetas[96]), [zeta_base3] "r"(&zetas[112])
  );

    //stage8, repeat bound=128
    //(v24, v6), (v25, v7), (v26, v8), (v27, v9), (v28, v10), (v29, v11), (v30, v12), (v31, v13) with m1 are input coefficients
    //all have their own zetas
    //v16, v18, v20, v22, v24, v26, v28, v30 with m2 is output coefficients
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v16, v24, v6, v0\n"
    "vle32.v v0, (%[zeta_base1])\n"
    "vbutterfly.ct.vvm v18, v25, v7, v0\n"
    "vle32.v v0, (%[zeta_base2])\n" 
    "vbutterfly.ct.vvm v20, v26, v8, v0\n"
    "vle32.v v0, (%[zeta_base3])\n"
    "vbutterfly.ct.vvm v22, v27, v9, v0\n"
    "vle32.v v0, (%[zeta_base4])\n"
    "vbutterfly.ct.vvm v24, v28, v10, v0\n"
    "vle32.v v0, (%[zeta_base5])\n"
    "vbutterfly.ct.vvm v26, v29, v11, v0\n"
    "vle32.v v0, (%[zeta_base6])\n"
    "vbutterfly.ct.vvm v28, v30, v12, v0\n"
    "vle32.v v0, (%[zeta_base7])\n"
    "vbutterfly.ct.vvm v30, v31, v13, v0\n"
    :
    : [zeta_vl] "r"(avl), [zeta_base0] "r"(&zetas[128]), [zeta_base1] "r"(&zetas[144]), [zeta_base2] "r"(&zetas[160]), [zeta_base3] "r"(&zetas[176]),
      [zeta_base4] "r"(&zetas[192]), [zeta_base5] "r"(&zetas[208]), [zeta_base6] "r"(&zetas[224]), [zeta_base7] "r"(&zetas[240])
  );

    //store back, from bit-reverse order to normal-order
    //v8-v15 used to store index, v16-v31 contains coefficients
    uint32_t* tree_ptr=NULL;
    avl=128;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m8,ta,mu\n"
    "vle32.v v8, (%[index_base0])\n"
    "vsuxei32.v v16, (%[coeff_base]), v8\n"     // indexed store v16-v23
    "vle32.v v8, (%[index_base1])\n"
    "vsuxei32.v v24, (%[coeff_base]), v8\n"     // indexed store v24-v31
    :
    : [avl]"r"(avl), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&r[0]), [index_base1]"r"(&tree_byteoffset[128])
  );
}

/*************************************************
* Name:        invntt_cg_custom_reorder
*
* Description: Inverse NTT.
*              Out-of-place, constant geometry.
*              Input coefficients are in standard order and normal domain.
*              Output coefficient are smaller than Q in absolute value;
*              and output coefficients are in standard order and normal domain.
*              zeta in compact form, that zetas in v0 are not the same
*
* Arguments:   - int32_t p[SABER_N]: input/output coefficient array
*
* Requirements: -VLEN=512
**************************************************/
void invntt_cg_custom_reorder(int32_t r[SABER_N]){
    int32_t r_temp0[256];
    unsigned int i;
    int32_t* r_ptr0 = NULL;
    int32_t* r_ptr1 = NULL;
    int32_t* r_ptr2 = NULL;
    int32_t* zeta_ptr = inv_zetas_inorder;

    vint32m1_t v_zetas;
    //for stage1 and stage2
    vint32m1_t v_in0,v_in1;
    vint32m2_t v_out0;
    //for stage3 to stage8
    vint32m4_t v_in2,v_in3;
    vint32m8_t v_out1;

    //stage1,same num=1
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 7));//repeat time here should be set large enough
    for(i=0;i<8;i++){//totally use 128zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,16);//16bf use 16zetas
        zeta_ptr+=16;
        v_in0=vle32_v_i32m1(r_ptr0,16);
        v_in1=vle32_v_i32m1(r_ptr1,16);
        r_ptr0+=16;
        r_ptr1+=16;
        v_out0=vbutterfly_gs_vvm_i32m2(v_in0,v_in1,v_zetas);
        vse32_v_i32m2(r_ptr2,v_out0,32);
        r_ptr2+=32;
    }

    //stage2, same num=2
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(1, 7));
    for(i=0;i<8;i++){//totally use 64zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,8);//16bf use 8zetas
        zeta_ptr+=8;
        v_in0=vle32_v_i32m1(r_ptr0,16);
        v_in1=vle32_v_i32m1(r_ptr1,16);
        r_ptr0+=16;
        r_ptr1+=16;
        v_out0=vbutterfly_gs_vvm_i32m2(v_in0,v_in1,v_zetas);
        vse32_v_i32m2(r_ptr2,v_out0,32);
        r_ptr2+=32;
    }

    //stage3, same num=4
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(2, 7));
    for(i=0;i<2;i++){//totally use 32zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,16);//64bf use 16zetas
        zeta_ptr+=16;
        v_in2=vle32_v_i32m4(r_ptr0,64);
        v_in3=vle32_v_i32m4(r_ptr1,64);
        r_ptr0+=64;
        r_ptr1+=64;
        vsetvl_e32m4_wrapper(64);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,128);
        r_ptr2+=128;
    }

    //stage4, same num=8
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(3, 7));
    for(i=0;i<2;i++){//totally use 16zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,8);//64bf use 8zetas
        zeta_ptr+=8;
        v_in2=vle32_v_i32m4(r_ptr0,64);
        v_in3=vle32_v_i32m4(r_ptr1,64);
        r_ptr0+=64;
        r_ptr1+=64;
        vsetvl_e32m4_wrapper(64);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,128);
        r_ptr2+=128;
    }

    //stage5, same num=16
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(4, 7));
    for(i=0;i<2;i++){//totally use 8zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,4);//64bf use 4zetas
        zeta_ptr+=4;
        v_in2=vle32_v_i32m4(r_ptr0,64);
        v_in3=vle32_v_i32m4(r_ptr1,64);
        r_ptr0+=64;
        r_ptr1+=64;
        vsetvl_e32m4_wrapper(64);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,128);
        r_ptr2+=128;
    }

    //stage6, same num=32
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(5, 7));
    for(i=0;i<2;i++){//totally use 4zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,2);//64bf use 2zetas
        zeta_ptr+=2;
        v_in2=vle32_v_i32m4(r_ptr0,64);
        v_in3=vle32_v_i32m4(r_ptr1,64);
        r_ptr0+=64;
        r_ptr1+=64;
        vsetvl_e32m4_wrapper(64);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,128);
        r_ptr2+=128;
    }

    //stage7, same num=64
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(6, 7));
    for(i=0;i<2;i++){//totally use 2zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,1);//64bf use 1zetas
        zeta_ptr+=1;
        v_in2=vle32_v_i32m4(r_ptr0,64);
        v_in3=vle32_v_i32m4(r_ptr1,64);
        r_ptr0+=64;
        r_ptr1+=64;
        vsetvl_e32m4_wrapper(64);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vse32_v_i32m8(r_ptr2,v_out1,128);
        r_ptr2+=128;
    }

    vuint32m8_t v_offset;
    uint32_t* tree_ptr=tree_byteoffset;

    //stage8, same num=128, all use the same zeta
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    csr_zeta_selMode_rw(GET_ZETA_MODE(7, 7));
    v_zetas=vle32_v_i32m1(zeta_ptr,1);
    for(i=0;i<2;i++){
        v_in2=vle32_v_i32m4(r_ptr0,64);
        v_in3=vle32_v_i32m4(r_ptr1,64);
        v_offset=vle32_v_u32m8(tree_ptr,128);
        tree_ptr+=128;
        r_ptr0+=64;
        r_ptr1+=64;
        vsetvl_e32m4_wrapper(64);
        v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
        vsetvl_e32m8_wrapper(128);
        vsuxei32_v_i32m8(r,v_offset,v_out1,128);
    }
}

/*************************************************
* Name:        invntt_cg_custom_reorder_asm
*
* Description: Inverse NTT.
*              Out-of-place, constant geometry.
*              Input coefficients are in standard order and normal domain.
*              Output coefficient are smaller than Q in absolute value;
*              and output coefficients are in standard order and normal domain.
*              zeta in compact form, that zetas in v0 are not the same
*
* Arguments:   - int32_t p[SABER_N]: input/output coefficient array
*
* Requirements: -VLEN=512
**************************************************/
void invntt_cg_custom_reorder_asm(int32_t r[SABER_N]){
    size_t load_zeta_vl;
    size_t avl;//m1

    //stage1, same num=1
    //v4 and v6 with m1 stores input coefficients
    //v16, v18, v20, v22, v24, v26, v28, v30 with m2 stores output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 7));//repeat time here should be set large enough
    load_zeta_vl=16;
    __asm__ __volatile__ ( 
    "vsetvli zero, %[zeta_vl], e32, m1, ta, mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vle32.v v4, (%[top_base_0])\n"
    "vle32.v v6, (%[btm_base_0])\n"
    "vbutterfly.gs.vvm v16, v4, v6, v0\n"
    "vle32.v v0, (%[zeta_base_1])\n"
    "vle32.v v4, (%[top_base_1])\n"
    "vle32.v v6, (%[btm_base_1])\n"
    "vbutterfly.gs.vvm v18, v4, v6, v0\n"
    "vle32.v v0, (%[zeta_base_2])\n"
    "vle32.v v4, (%[top_base_2])\n"
    "vle32.v v6, (%[btm_base_2])\n"
    "vbutterfly.gs.vvm v20, v4, v6, v0\n"
    "vle32.v v0, (%[zeta_base_3])\n"
    "vle32.v v4, (%[top_base_3])\n"
    "vle32.v v6, (%[btm_base_3])\n"
    "vbutterfly.gs.vvm v22, v4, v6, v0\n"
    "vle32.v v0, (%[zeta_base_4])\n"
    "vle32.v v4, (%[top_base_4])\n"
    "vle32.v v6, (%[btm_base_4])\n"
    "vbutterfly.gs.vvm v24, v4, v6, v0\n"
    "vle32.v v0, (%[zeta_base_5])\n"
    "vle32.v v4, (%[top_base_5])\n"
    "vle32.v v6, (%[btm_base_5])\n"
    "vbutterfly.gs.vvm v26, v4, v6, v0\n"
    "vle32.v v0, (%[zeta_base_6])\n"
    "vle32.v v4, (%[top_base_6])\n"
    "vle32.v v6, (%[btm_base_6])\n"
    "vbutterfly.gs.vvm v28, v4, v6, v0\n"
    "vle32.v v0, (%[zeta_base_7])\n"
    "vle32.v v4, (%[top_base_7])\n"
    "vle32.v v6, (%[btm_base_7])\n"
    "vbutterfly.gs.vvm v30, v4, v6, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&inv_zetas_inorder[0]), [top_base_0]"r"(&r[0]), [btm_base_0]"r"(&r[128]),
      [zeta_base_1]"r"(&inv_zetas_inorder[16]), [top_base_1]"r"(&r[16]), [btm_base_1]"r"(&r[144]),
      [zeta_base_2]"r"(&inv_zetas_inorder[32]), [top_base_2]"r"(&r[32]), [btm_base_2]"r"(&r[160]),
      [zeta_base_3]"r"(&inv_zetas_inorder[48]), [top_base_3]"r"(&r[48]), [btm_base_3]"r"(&r[176]),
      [zeta_base_4]"r"(&inv_zetas_inorder[64]), [top_base_4]"r"(&r[64]), [btm_base_4]"r"(&r[192]),
      [zeta_base_5]"r"(&inv_zetas_inorder[80]), [top_base_5]"r"(&r[80]), [btm_base_5]"r"(&r[208]),
      [zeta_base_6]"r"(&inv_zetas_inorder[96]), [top_base_6]"r"(&r[96]), [btm_base_6]"r"(&r[224]),
      [zeta_base_7]"r"(&inv_zetas_inorder[112]), [top_base_7]"r"(&r[112]), [btm_base_7]"r"(&r[240])
  );
  
    //stage2, same num=2
    // (v16, v24), (v18, v26), (v20, v28), (v22, v30) with m2 are input coefficients
    // v4, v8, v12, v16 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(1, 7));
    avl=32;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v4, v16, v24, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v8, v18, v26, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v12, v20, v28, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v16, v22, v30, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&inv_zetas_inorder[128]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&inv_zetas_inorder[144]), [zeta_base_2]"r"(&inv_zetas_inorder[160]), [zeta_base_3]"r"(&inv_zetas_inorder[176])
  );

    //stage3, same num=4
    // (v4, v12), (v6, v14), (v8, v16), (v10, v18) with m2 are input coefficients
    // v20, v24, v28, v4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(2, 7));
    load_zeta_vl=8;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v20, v4, v12, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v24, v6, v14, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v28, v8, v16, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v4, v10, v18, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&inv_zetas_inorder[192]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&inv_zetas_inorder[200]), [zeta_base_2]"r"(&inv_zetas_inorder[208]), [zeta_base_3]"r"(&inv_zetas_inorder[216])
  );

    //stage4, same num=8
    // (v20, v28), (v22, v30), (v24, v4), (v26, v6) with m2 are input coefficients
    // v8, v12, v16, v20 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(3, 7));
    load_zeta_vl=4;
     __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v8, v20, v28, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v12, v22, v30, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v16, v24, v4, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v20, v26, v6, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&inv_zetas_inorder[224]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&inv_zetas_inorder[228]), [zeta_base_2]"r"(&inv_zetas_inorder[232]), [zeta_base_3]"r"(&inv_zetas_inorder[236])
  );

    //stage5, same num=16
    // (v8, v16), (v10, v18), (v12, v20), (v14, v22) with m2 are input coefficients
    // v24, v28, v4, v8 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(4, 7));
    load_zeta_vl=2;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v24, v8, v16, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v28, v10, v18, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v4, v12, v20, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v8, v14, v22, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&inv_zetas_inorder[240]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&inv_zetas_inorder[242]), [zeta_base_2]"r"(&inv_zetas_inorder[244]), [zeta_base_3]"r"(&inv_zetas_inorder[246])
  );

    //stage6, same num=32
    // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m2 are input coefficients
    //  v12, v16, v20, v24 with m4 are output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(5, 7));
    load_zeta_vl=1;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v12, v24, v4, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v16, v26, v6, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v20, v28, v8, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v24, v30, v10, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&inv_zetas_inorder[248]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&inv_zetas_inorder[249]), [zeta_base_2]"r"(&inv_zetas_inorder[250]), [zeta_base_3]"r"(&inv_zetas_inorder[251])
  );

    //stage7, same num=64
    // (v12, v20), (v14, v22), (v16, v24), (v18, v26) with m2 are input coefficients
    // v28, v4, v8, v12 with m4 are output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(6, 7));
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v28, v12, v20, v0\n"
    "vbutterfly.gs.vvm v4, v14, v22, v0\n"
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v8, v16, v24, v0\n"
    "vbutterfly.gs.vvm v12, v18, v26, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&inv_zetas_inorder[252]), [avl]"r"(avl), [zeta_base_1]"r"(&inv_zetas_inorder[253])
  );

    //stage8, same num=128
    // (v28, v8), (v30, v10), (v4, v12), (v6, v14) with m2 are input coefficients
    // v16, v20, v24, v28 with m4 are output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(7, 7));
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
    "vbutterfly.gs.vvm v16, v28, v8, v0\n"
    "vbutterfly.gs.vvm v20, v30, v10, v0\n"
    "vbutterfly.gs.vvm v24, v4, v12, v0\n"
    "vbutterfly.gs.vvm v28, v6, v14, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&inv_zetas_inorder[254]), [avl]"r"(avl)
  );

    //store back, back, from bit-reverse order to normal-order
    //v8-v15 used to store index, v16-v31 contains coefficients
    uint32_t* tree_ptr=NULL;
    avl=128;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m8,ta,mu\n"
    "vle32.v v8, (%[index_base0])\n"
    "vsuxei32.v v16, (%[coeff_base]), v8\n"     // indexed store v16-v23
    "vle32.v v8, (%[index_base1])\n"
    "vsuxei32.v v24, (%[coeff_base]), v8\n"     // indexed store v24-v31
    :
    : [avl]"r"(avl), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&r[0]), [index_base1]"r"(&tree_byteoffset[128])
  );
}
#elif (VLEN==1024)
/*************************************************
* Name:        ntt_cg_custom_reorder
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
*
* Requirements: -VLEN=512
**************************************************/
void ntt_cg_custom_reorder(int32_t r[SABER_N]){
    int32_t r_temp0[256];
    unsigned int i, j;
    int32_t* r_ptr0 = NULL;
    int32_t* r_ptr1 = NULL;
    int32_t* r_ptr2 = NULL;

    vint32m1_t v_zetas;
    //for stage 1 to stage 6
    vint32m4_t v_in0,v_in1;
    vint32m8_t v_out0;
    //for stage 7 to stage 8
    vint32m1_t v_in2,v_in3;
    vint32m2_t v_out1;

    //stage1, repeat bound=1
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
    v_zetas=vle32_v_i32m1(&zetas[1],1);
    v_in0=vle32_v_i32m4(r_ptr0,128);
    v_in1=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out0=vbutterfly_ct_vvm_i32m8(v_in0,v_in1,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out0,256);

    //stage2, repeat bound=2
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));
    v_zetas=vle32_v_i32m1(&zetas[2],2);
    v_in0=vle32_v_i32m4(r_ptr0,128);
    v_in1=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out0=vbutterfly_ct_vvm_i32m8(v_in0,v_in1,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out0,256);

    //stage3, repeat bound=4
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
    v_zetas=vle32_v_i32m1(&zetas[4],4);
    v_in0=vle32_v_i32m4(r_ptr0,128);
    v_in1=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out0=vbutterfly_ct_vvm_i32m8(v_in0,v_in1,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out0,256);

    //stage4, repeat bound=8
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
    v_zetas=vle32_v_i32m1(&zetas[8],8);
    v_in0=vle32_v_i32m4(r_ptr0,128);
    v_in1=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out0=vbutterfly_ct_vvm_i32m8(v_in0,v_in1,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out0,256);

    //stage5, repeat bound=16
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
    v_zetas=vle32_v_i32m1(&zetas[16],16);
    v_in0=vle32_v_i32m4(r_ptr0,128);
    v_in1=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out0=vbutterfly_ct_vvm_i32m8(v_in0,v_in1,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out0,256);

    //stage6, actual repeat bound=32
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));
    v_zetas=vle32_v_i32m1(&zetas[32],32);
    v_in0=vle32_v_i32m4(r_ptr0,128);
    v_in1=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out0=vbutterfly_ct_vvm_i32m8(v_in0,v_in1,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out0,256);

    //stage7, actual repeat bound=64, exceed v0's capacity, so set repeat bound to 32
    for(i=0;i<2;i++){
        int offset=32*i;
        r_ptr0=r+offset;
        r_ptr1=&r[128+offset];
        r_ptr2=r_temp0+(offset<<1);
        v_zetas=vle32_v_i32m1(&zetas[64+offset],32);
        for(j=0;j<2;j++){
            v_in2=vle32_v_i32m1(r_ptr0,32);
            v_in3=vle32_v_i32m1(r_ptr1,32);
            r_ptr0+=64;
            r_ptr1+=64;
            v_out1=vbutterfly_ct_vvm_i32m2(v_in2,v_in3,v_zetas);
            vse32_v_i32m2(r_ptr2,v_out1,64);
            r_ptr2+=128;
        }
    }

    vuint32m2_t v_offset;
    uint32_t* tree_ptr=NULL;
    //stage8, actual repeat bound=128, exceed v0's capacity, so set repeat bound to 32
    for(i=0;i<4;i++){
        int offset=32*i;
        r_ptr0=r_temp0+offset;
        r_ptr1=&r_temp0[128+offset];
        tree_ptr=tree_byteoffset+(offset<<1);
        v_offset=vle32_v_u32m2(tree_ptr,64);
        v_zetas=vle32_v_i32m1(&zetas[128+offset],32);
        v_in2=vle32_v_i32m1(r_ptr0,32);
        v_in3=vle32_v_i32m1(r_ptr1,32);
        v_out1=vbutterfly_ct_vvm_i32m2(v_in2,v_in3,v_zetas);
        vsuxei32_v_i32m2(r,v_offset,v_out1,64);
    }
}

/*************************************************
* Name:        invntt_cg_custom_reorder
*
* Description: Inverse NTT.
*              Out-of-place, constant geometry.
*              Input coefficients are in standard order and normal domain.
*              Output coefficient are smaller than Q in absolute value;
*              and output coefficients are in standard order and normal domain.
*              zeta in compact form, that zetas in v0 are not the same
*
* Arguments:   - int32_t p[N]: input/output coefficient array
*
* Requirements: -VLEN=512
**************************************************/
void invntt_cg_custom_reorder(int32_t r[SABER_N]){
    int32_t r_temp0[256];
    unsigned int i;
    int32_t* r_ptr0 = NULL;
    int32_t* r_ptr1 = NULL;
    int32_t* r_ptr2 = NULL;
    int32_t* zeta_ptr = inv_zetas_inorder;

    vint32m1_t v_zetas;
    //for stage1 to stage2
    vint32m1_t v_in0,v_in1;
    vint32m2_t v_out0;
    //for stage3 to stage8
    vint32m4_t v_in2,v_in3;
    vint32m8_t v_out1;

    //stage1,same num=1
    r_ptr0=r;
    r_ptr1=&r[128];
    r_ptr2=r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 7));//repeat time here should be set large enough
    for(i=0;i<4;i++){//totally use 128zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,32);//32bf use 32zetas
        zeta_ptr+=32;
        v_in0=vle32_v_i32m1(r_ptr0,32);
        v_in1=vle32_v_i32m1(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        v_out0=vbutterfly_gs_vvm_i32m2(v_in0,v_in1,v_zetas);
        vse32_v_i32m2(r_ptr2,v_out0,64);
        r_ptr2+=64;
    }

    //stage2, same num=2
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(1, 7));
    for(i=0;i<4;i++){//totally use 64zetas
        v_zetas=vle32_v_i32m1(zeta_ptr,16);//32bf use 16zetas
        zeta_ptr+=16;
        v_in0=vle32_v_i32m1(r_ptr0,32);
        v_in1=vle32_v_i32m1(r_ptr1,32);
        r_ptr0+=32;
        r_ptr1+=32;
        v_out0=vbutterfly_gs_vvm_i32m2(v_in0,v_in1,v_zetas);
        vse32_v_i32m2(r_ptr2,v_out0,64);
        r_ptr2+=64;
    }

    //stage3, same num=4
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(2, 7));
    v_zetas=vle32_v_i32m1(zeta_ptr,32);//128bf use 32zetas
    zeta_ptr+=32;
    v_in2=vle32_v_i32m4(r_ptr0,128);
    v_in3=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out1,256);

    //stage4, same num=8
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(3, 7));
    v_zetas=vle32_v_i32m1(zeta_ptr,16);//128bf use 16zetas
    zeta_ptr+=16;
    v_in2=vle32_v_i32m4(r_ptr0,128);
    v_in3=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out1,256);

    //stage5, same num=16
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(4, 7));
    v_zetas=vle32_v_i32m1(zeta_ptr,8);//128bf use 8zetas
    zeta_ptr+=8;
    v_in2=vle32_v_i32m4(r_ptr0,128);
    v_in3=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out1,256);

    //stage6, same num=32
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    r_ptr2=r;
    csr_zeta_selMode_rw(GET_ZETA_MODE(5, 7));
    v_zetas=vle32_v_i32m1(zeta_ptr,4);//128bf use 4zetas
    zeta_ptr+=4;
    v_in2=vle32_v_i32m4(r_ptr0,128);
    v_in3=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out1,256);

    //stage7, same num=64
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(6, 7));
    v_zetas=vle32_v_i32m1(zeta_ptr,2);//128bf use 2zetas
    zeta_ptr+=2;
    v_in2=vle32_v_i32m4(r_ptr0,128);
    v_in3=vle32_v_i32m4(r_ptr1,128);
    vsetvl_e32m4_wrapper(128);
    v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
    vse32_v_i32m8(r_ptr2,v_out1,256);

    vuint32m8_t v_offset;
    uint32_t* tree_ptr=tree_byteoffset;

    //stage8, same num=128, all use the same zeta
    r_ptr0=r_temp0;
    r_ptr1=&r_temp0[128];
    csr_zeta_selMode_rw(GET_ZETA_MODE(7, 7));
    v_zetas=vle32_v_i32m1(zeta_ptr,1);
    v_in2=vle32_v_i32m4(r_ptr0,128);
    v_in3=vle32_v_i32m4(r_ptr1,128);
    v_offset=vle32_v_u32m8(tree_ptr,256);
    vsetvl_e32m4_wrapper(128);
    v_out1=vbutterfly_gs_vvm_i32m8(v_in2,v_in3,v_zetas);
    vsuxei32_v_i32m8(r,v_offset,v_out1,256);
}

#else
#error "VLEN must be 256/512/1024"
#endif