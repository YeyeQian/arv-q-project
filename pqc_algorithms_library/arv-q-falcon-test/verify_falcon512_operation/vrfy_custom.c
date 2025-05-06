#include "inner.h"
#include <riscv_vector.h>
#include <assert.h>
#include "../apis/custom_inst_api.h"
#include "params_custom.h"

/* ===================================================================== */
/*
 * Constants for NTT.
 *
 *   n = 2^logn  (2 <= n <= 1024)
 *   phi = X^n + 1
 *   q = 12289
 *   q0i = -1/q mod 2^16
 *   R = 2^16 mod q
 *   R2 = 2^32 mod q
 */

// #define Q     12289
// #define Q0I   12287
// #define R      4091
// #define R2    10952

/*************************************************
* Name:        mq_NTT_custom
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in bitreversed order
*              zeta in compact form, that zetas in v0 are not the same.
*              
*
* Arguments:   - uint16_t a[512]: pointer to input/output vector of elements
*                                of Zq
*
* Requirements: -VLEN=512
**************************************************/
void
mq_NTT_custom(uint16_t *a, unsigned logn)
{
    assert(logn == 9);

    const unsigned int elem_per_vreg = VLEN / 16;               // VLEN=512, elem_per_vreg=32
    const unsigned int vbf_times = (512 / 2) / elem_per_vreg;   // VLEN=512, vbf_times=

    uint16_t r_temp0[512];
    unsigned int i, j;
    uint16_t* r_ptr0 = NULL;
    uint16_t* r_ptr1 = NULL;
    uint16_t* r_ptr2 = NULL;

    vuint16m1_t v_zetas;
    vuint16m1_t v_in0, v_in1;
    vuint16m2_t v_out0;

    //stage1, repeat bound=1
    r_ptr0 = a;
    r_ptr1 = &a[256];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
    v_zetas = vle16_v_u16m1(&zetas_cg_512[1],1);
    for(i = 0; i < vbf_times; i++) {
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_ct_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage2, repeat bound=2
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[256];
    r_ptr2 = a;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));
    v_zetas = vle16_v_u16m1(&zetas_cg_512[2],2);
    for(i = 0; i < vbf_times; i++) {
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_ct_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage3, repeat bound=4
    r_ptr0 = a;
    r_ptr1 = &a[256];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
    v_zetas = vle16_v_u16m1(&zetas_cg_512[4],4);
    for(i = 0; i < vbf_times; i++) {
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_ct_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage4, repeat bound=8
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[256];
    r_ptr2 = a;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
    v_zetas = vle16_v_u16m1(&zetas_cg_512[8],8);
    for(i = 0; i < vbf_times; i++) {
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_ct_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage5, repeat bound=16
    r_ptr0 = a;
    r_ptr1 = &a[256];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
    v_zetas = vle16_v_u16m1(&zetas_cg_512[16],16);
    for(i = 0; i < vbf_times; i++) {
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_ct_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage6, repeat bound=32
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[256];
    r_ptr2 = a;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));
    v_zetas = vle16_v_u16m1(&zetas_cg_512[32],32);
    for(i = 0; i < vbf_times; i++) {
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_ct_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage7, actual repeat bound=64, exceed v0's capacity, so set repeat bound to elem_per_vreg = 32
    // Two-layer loop, total loop number is also vbf_times
    for(i = 0; i < (64 / elem_per_vreg); i++) {                 // 2 iteration
        uint32_t offset = 32 * i;
        r_ptr0 = a + offset;
        r_ptr1 = &a[256 + offset];
        r_ptr2 = r_temp0 + (offset<<1);
        v_zetas = vle16_v_u16m1(&zetas_cg_512[64 + offset], elem_per_vreg);
        for(j = 0; j < vbf_times / (64 / elem_per_vreg); j++) {      // 4 iteration     
            v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
            v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
            r_ptr0 += (2*elem_per_vreg);
            r_ptr1 += (2*elem_per_vreg);
            v_out0 = vbutterfly_ct_vvm_u16m2(v_in0, v_in1, v_zetas);
            vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
            r_ptr2 += (4*elem_per_vreg);
        }   
    }

    //stage8, actual repeat bound=128, exceed v0's capacity, so set repeat bound to elem_per_vreg = 32
    // Two-layer loop, total loop number is also vbf_times
    for(i = 0; i < (128 / elem_per_vreg); i++) {                 // 4 iteration
        uint32_t offset = 32 * i;
        r_ptr0 = r_temp0 + offset;
        r_ptr1 = &r_temp0[256 + offset];
        r_ptr2 = a + (offset<<1);
        v_zetas = vle16_v_u16m1(&zetas_cg_512[128 + offset], elem_per_vreg);
        for(j = 0; j < vbf_times / (128 / elem_per_vreg); j++) {      // 2 iteration
            v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
            v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
            r_ptr0 += (4*elem_per_vreg);
            r_ptr1 += (4*elem_per_vreg);
            v_out0 = vbutterfly_ct_vvm_u16m2(v_in0, v_in1, v_zetas);
            vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
            r_ptr2 += (8*elem_per_vreg);          
        }   
    }

    //stage9, actual repeat bound=256, exceed v0's capacity, so set repeat bound to elem_per_vreg = 32
    for(i = 0; i < vbf_times; i++) {
        uint32_t offset = 32 * i;
        r_ptr0 = a + offset;
        r_ptr1 = &a[256 + offset];
        r_ptr2 = r_temp0 + (offset<<1);
        v_zetas = vle16_v_u16m1(&zetas_cg_512[256 + offset], elem_per_vreg);
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        v_out0 = vbutterfly_ct_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
    }

    copy_array(a, r_temp0);
}

/*************************************************
* Name:        mq_iNTT_custom
*
* Description: Inverse NTT.
*              Out-of-place, constant geometry.
*              Input coefficients are in bit-reversed order and normal domain.
*              Output coefficient are smaller than Q in absolute value;
*              and output coefficients are in standard order and normal domain.
*
* Arguments:   - uint16_t p[FALCON_N]: input/output coefficient array
*
* Requirements: -VLEN=512
**************************************************/
void
mq_iNTT_custom(uint16_t *a, unsigned logn)
{
    assert(logn == 9);

    bitreverse_standard_pos_transfer(a);

    const unsigned int elem_per_vreg = VLEN / 16;
    const unsigned int vbf_times = (512 / 2) / elem_per_vreg;

    uint16_t r_temp0[512];
    unsigned int i, j;
    uint16_t* r_ptr0 = NULL;
    uint16_t* r_ptr1 = NULL;
    uint16_t* r_ptr2 = NULL;
    const uint16_t* zeta_ptr = zetas_inv_cg_inorder_512;

    vuint16m1_t v_zetas;
    vuint16m1_t v_in0, v_in1;
    vuint16m2_t v_out0;

    //stage1, same num=1
    r_ptr0 = a;
    r_ptr1 = &a[256];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 7));       // repeat time here should be set large enough
    for(i = 0; i < vbf_times; i++) {                // totally use 256 zetas
        v_zetas = vle16_v_u16m1(zeta_ptr, elem_per_vreg);      // 32 bf use 32 zetas
        zeta_ptr += elem_per_vreg;
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_gs_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage2, same num=2
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[256];
    r_ptr2 = a;
    csr_zeta_selMode_rw(GET_ZETA_MODE(1, 7));       // repeat time here should be set large enough
    for(i = 0; i < vbf_times; i++) {                // totally use 128 zetas
        v_zetas = vle16_v_u16m1(zeta_ptr, elem_per_vreg>>1);      // 32 bf use 16 zetas
        zeta_ptr += (elem_per_vreg>>1);
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_gs_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage3, same num=4
    r_ptr0 = a;
    r_ptr1 = &a[256];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(2, 7));       // repeat time here should be set large enough
    for(i = 0; i < vbf_times; i++) {                // totally use 64 zetas
        v_zetas = vle16_v_u16m1(zeta_ptr, elem_per_vreg>>2);      // 32 bf use 8 zetas
        zeta_ptr += (elem_per_vreg>>2);
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_gs_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage4, same num=8
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[256];
    r_ptr2 = a;
    csr_zeta_selMode_rw(GET_ZETA_MODE(3, 7));       // repeat time here should be set large enough
    for(i = 0; i < vbf_times; i++) {                // totally use 32 zetas
        v_zetas = vle16_v_u16m1(zeta_ptr, elem_per_vreg>>3);      // 32 bf use 4 zetas
        zeta_ptr += (elem_per_vreg>>3);
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_gs_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage5, same num=16
    r_ptr0 = a;
    r_ptr1 = &a[256];
    r_ptr2 = r_temp0;
    csr_zeta_selMode_rw(GET_ZETA_MODE(4, 7));       // repeat time here should be set large enough
    for(i = 0; i < vbf_times; i++) {                // totally use 16 zetas
        v_zetas = vle16_v_u16m1(zeta_ptr, elem_per_vreg>>4);      // 32 bf use 2 zetas
        zeta_ptr += (elem_per_vreg>>4);
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_gs_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage6, same num=32
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[256];
    r_ptr2 = a;
    csr_zeta_selMode_rw(GET_ZETA_MODE(5, 7));       // repeat time here should be set large enough
    for(i = 0; i < vbf_times; i++) {                // totally use 8 zetas
        v_zetas = vle16_v_u16m1(zeta_ptr, elem_per_vreg>>5);      // 32 bf use 1 zetas
        zeta_ptr += (elem_per_vreg>>5);
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_gs_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);      
    }

    //stage7, same num=64
    r_ptr0 = a;
    r_ptr1 = &a[256];
    r_ptr2 = r_temp0;
    for(i = 0; i < vbf_times / 2; i++) {            // totally use 4 zetas
        v_zetas = vle16_v_u16m1(zeta_ptr, 1);       // 64 bf use 1 zetas
        zeta_ptr += 1;
        for(j = 0; j < 2; j++) {
            v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
            v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
            r_ptr0 += elem_per_vreg;
            r_ptr1 += elem_per_vreg;
            v_out0 = vbutterfly_gs_vvm_u16m2(v_in0, v_in1, v_zetas);
            vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
            r_ptr2 += (2*elem_per_vreg);   
        }   
    }

    //stage8, same num=128
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[256];
    r_ptr2 = a;
    for(i = 0; i < vbf_times / 4; i++) {            // totally use 2 zetas
        v_zetas = vle16_v_u16m1(zeta_ptr, 1);       // 128 bf use 1 zetas
        zeta_ptr += 1;
        for(j = 0; j < 4; j++) {
            v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
            v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
            r_ptr0 += elem_per_vreg;
            r_ptr1 += elem_per_vreg;
            v_out0 = vbutterfly_gs_vvm_u16m2(v_in0, v_in1, v_zetas);
            vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
            r_ptr2 += (2*elem_per_vreg);   
        }   
    }

    //stage9, same num=256
    r_ptr0 = a;
    r_ptr1 = &a[256];
    r_ptr2 = r_temp0;
    v_zetas = vle16_v_u16m1(zeta_ptr, 1);           // 256 bf use 1 zetas
    for(i = 0; i < vbf_times; i++) {                // totally use 1 zetas
        v_in0 = vle16_v_u16m1(r_ptr0, elem_per_vreg);
        v_in1 = vle16_v_u16m1(r_ptr1, elem_per_vreg);
        r_ptr0 += elem_per_vreg;
        r_ptr1 += elem_per_vreg;
        v_out0 = vbutterfly_gs_vvm_u16m2(v_in0, v_in1, v_zetas);
        vse16_v_u16m2(r_ptr2, v_out0, 2*elem_per_vreg);
        r_ptr2 += (2*elem_per_vreg);
    }

    copy_array(a, r_temp0);
    bitreverse_standard_pos_transfer(a);           
}

/*************************************************
* Name:        mq_NTT_custom_asm
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in bitreversed order
*              zeta in compact form, that zetas in v0 are not the same.
*              
*
* Arguments:   - uint16_t a[512]: pointer to input/output vector of elements
*                                of Zq
*
* Requirements: -VLEN=512
**************************************************/
void
mq_NTT_custom_asm(uint16_t *a, unsigned logn)
{
    assert(logn == 9);

    size_t load_zeta_vl;
    size_t avl = 64;    //m2

    //stage1, repeat bound=1
    //v4 and v6 with m2 stores input coefficients
    //v16, v20, v24, v28 with m4 store output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
    load_zeta_vl=1;
    __asm__ __volatile__ ( 
    "vsetvli zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli zero, %[avl], e16, m2, ta, mu\n"
    "vle16.v v4, (%[top_base_0])\n"
    "vle16.v v6, (%[btm_base_0])\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    "vle16.v v4, (%[top_base_1])\n"
    "vle16.v v6, (%[btm_base_1])\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    "vle16.v v4, (%[top_base_2])\n"
    "vle16.v v6, (%[btm_base_2])\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    "vle16.v v4, (%[top_base_3])\n"
    "vle16.v v6, (%[btm_base_3])\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base]"r"(&zetas_cg_512[1]), [avl]"r"(avl),
     [top_base_0]"r"(&a[0]), [btm_base_0]"r"(&a[256]), [top_base_1]"r"(&a[64]), [btm_base_1]"r"(&a[320]),
     [top_base_2]"r"(&a[128]), [btm_base_2]"r"(&a[384]), [top_base_3]"r"(&a[192]), [btm_base_3]"r"(&a[448])
    );

    //stage2, repeat bound=2
    // (v16, v24), (v18, v26), (v20, v28), (v22, v30) with m2 are input coefficients
    // v4, v8, v12, v16 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));
    load_zeta_vl=2;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v4, v16, v24, v0\n"
    "vbutterfly.ct.vvm v8, v18, v26, v0\n"
    "vbutterfly.ct.vvm v12, v20, v28, v0\n"
    "vbutterfly.ct.vvm v16, v22, v30, v0\n"
    :
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_cg_512[2]), [avl]"r"(avl)
  );

    //stage3, repeat bound=4
    // (v4, v12), (v6, v14), (v8, v16), (v10, v18) with m2 are input coefficients
    // v20, v24, v28, v4 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
    load_zeta_vl=4;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v20, v4, v12, v0\n"
    "vbutterfly.ct.vvm v24, v6, v14, v0\n"
    "vbutterfly.ct.vvm v28, v8, v16, v0\n"
    "vbutterfly.ct.vvm v4, v10, v18, v0\n"
    :
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_cg_512[4]), [avl]"r"(avl)
  );

    //stage4, repeat bound=8
    // (v20, v28), (v22, v30), (v24, v4), (v26, v6) with m2 are input coefficients
    // v8, v12, v16, v20 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
    load_zeta_vl=8;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v8, v20, v28, v0\n"
    "vbutterfly.ct.vvm v12, v22, v30, v0\n"
    "vbutterfly.ct.vvm v16, v24, v4, v0\n"
    "vbutterfly.ct.vvm v20, v26, v6, v0\n"
    :
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_cg_512[8]), [avl]"r"(avl)
  );

    //stage5, repeat bound=16
    // (v8, v16), (v10, v18), (v12, v20), (v14, v22) with m2 are input coefficients
    // v24, v28, v4, v8 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
    load_zeta_vl=16;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v24, v8, v16, v0\n"
    "vbutterfly.ct.vvm v28, v10, v18, v0\n"
    "vbutterfly.ct.vvm v4, v12, v20, v0\n"
    "vbutterfly.ct.vvm v8, v14, v22, v0\n"
    :
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_cg_512[16]), [avl]"r"(avl)
  );

    //stage6, repeat bound=32
    // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m1 are input coefficients, with same zetas
    // (v25, v5), (v27, v7), (v29, v9), (v31, v11) with m1 are input coefficients, with same zetas
    // v12, v14, v16, v18, v20, v22, v2, v4 with m2 are output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));
    avl = 32;               // m1
    load_zeta_vl = 32;      // load_zeta_vl is equal to avl
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vbutterfly.ct.vvm v12, v24, v4, v0\n"
    "vbutterfly.ct.vvm v16, v26, v6, v0\n"
    "vbutterfly.ct.vvm v20, v28, v8, v0\n"
    "vbutterfly.ct.vvm v2, v30, v10, v0\n"
    "vbutterfly.ct.vvm v14, v25, v5, v0\n"
    "vbutterfly.ct.vvm v18, v27, v7, v0\n"
    "vbutterfly.ct.vvm v22, v29, v9, v0\n"
    "vbutterfly.ct.vvm v4, v31, v11, v0\n"
    :
    : [zeta_vl] "r"(avl), [zeta_base] "r"(&zetas_cg_512[32])
  );

    //stage7, repeat bound=64
    // (v12, v20), (v13, v21), (v14, v22), (v15, v23),
    // (v16, v2), (v17, v3), (v18, v4), (v19, v5) with m1 are input coefficients
    // ((v12, v20), (v14, v22), (v16, v2), (v18, v4)) same zetas
    // ((v13, v21), (v15, v23), (v17, v3), (v19, v5)) same zetas
    // v24, v26, v28, v30, v6, v8, v10, v12 with m2 is output coefficients
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v24, v12, v20, v0\n"
    "vbutterfly.ct.vvm v28, v14, v22, v0\n"
    "vbutterfly.ct.vvm v6, v16, v2, v0\n"
    "vbutterfly.ct.vvm v10, v18, v4, v0\n"
    "vle16.v v0, (%[zeta_base1])\n" 
    "vbutterfly.ct.vvm v26, v13, v21, v0\n"
    "vbutterfly.ct.vvm v30, v15, v23, v0\n"
    "vbutterfly.ct.vvm v8, v17, v3, v0\n"
    "vbutterfly.ct.vvm v12, v19, v5, v0\n"
    :
    : [zeta_vl] "r"(avl), [zeta_base0] "r"(&zetas_cg_512[64]), [zeta_base1] "r"(&zetas_cg_512[96])
  );

    //stage8, repeat bound=128
    // (v24, v6), (v25, v7), (v26, v8), (v27, v9),
    // (v28, v10), (v29, v11), (v30, v12), (v31, v13) with m1 are input coefficients
    // (v24, v6)->v2 and (v28, v10)->v18 same zeta
    // (v25, v7)->v4 and (v29, v11)->v20 same zeta
    // (v26, v8)->14 and (v30, v12)->v22 same zeta  
    // (v27, v9)->v16 and (v31, v13)->v24 same zeta
    // v2, v4, v14, v16, v18, v20, v22, v24 with m2 is output coefficients
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v2, v24, v6, v0\n"
    "vbutterfly.ct.vvm v18, v28, v10, v0\n"
    "vle16.v v0, (%[zeta_base1])\n"
    "vbutterfly.ct.vvm v4, v25, v7, v0\n"    
    "vbutterfly.ct.vvm v20, v29, v11, v0\n"
    "vle16.v v0, (%[zeta_base2])\n" 
    "vbutterfly.ct.vvm v14, v26, v8, v0\n"      
    "vbutterfly.ct.vvm v22, v30, v12, v0\n"
    "vle16.v v0, (%[zeta_base3])\n"
    "vbutterfly.ct.vvm v16, v27, v9, v0\n"
    "vbutterfly.ct.vvm v24, v31, v13, v0\n"
    :
    : [zeta_vl] "r"(avl), 
      [zeta_base0] "r"(&zetas_cg_512[128]), [zeta_base1] "r"(&zetas_cg_512[160]), 
      [zeta_base2] "r"(&zetas_cg_512[192]), [zeta_base3] "r"(&zetas_cg_512[224])
  );

    //stage9, repeat bound=256
    // (v2, v18)->v6, (v3, v19)->v8, (v4, v20)->v10, (v5, v21)->v12, 
    // (v14, v22)->v26, (v15, v23)->v28, (v16, v24)->v30, (v17, v25)->v2 with m1 are input coefficients
    // all have their own zetas
    // v6, v8, v10, v12, v26, v28, v30, v2 with m2 is output coefficients
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v6, v2, v18, v0\n"
    "vle16.v v0, (%[zeta_base1])\n"
    "vbutterfly.ct.vvm v8, v3, v19, v0\n"
    "vle16.v v0, (%[zeta_base2])\n" 
    "vbutterfly.ct.vvm v10, v4, v20, v0\n"
    "vle16.v v0, (%[zeta_base3])\n"
    "vbutterfly.ct.vvm v12, v5, v21, v0\n"
    "vle16.v v0, (%[zeta_base4])\n"
    "vbutterfly.ct.vvm v26, v14, v22, v0\n"
    "vle16.v v0, (%[zeta_base5])\n"
    "vbutterfly.ct.vvm v28, v15, v23, v0\n"
    "vle16.v v0, (%[zeta_base6])\n"
    "vbutterfly.ct.vvm v30, v16, v24, v0\n"
    "vle16.v v0, (%[zeta_base7])\n"
    "vbutterfly.ct.vvm v2, v17, v25, v0\n"
    :
    : [zeta_vl] "r"(avl), 
      [zeta_base0] "r"(&zetas_cg_512[256]), [zeta_base1] "r"(&zetas_cg_512[288]), 
      [zeta_base2] "r"(&zetas_cg_512[320]), [zeta_base3] "r"(&zetas_cg_512[352]),
      [zeta_base4] "r"(&zetas_cg_512[384]), [zeta_base5] "r"(&zetas_cg_512[416]), 
      [zeta_base6] "r"(&zetas_cg_512[448]), [zeta_base7] "r"(&zetas_cg_512[480])
  );

    // store coefficients from vector register v2~v17 with m2 to memory
    avl = 64;   // m2
      __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vse16.v v6, (%[coeff_base0])\n"
    "vse16.v v8, (%[coeff_base1])\n"
    "vse16.v v10, (%[coeff_base2])\n"
    "vse16.v v12, (%[coeff_base3])\n"
    "vse16.v v26, (%[coeff_base4])\n"
    "vse16.v v28, (%[coeff_base5])\n"
    "vse16.v v30, (%[coeff_base6])\n"
    "vse16.v v2, (%[coeff_base7])\n"
    :
    : [avl] "r"(avl), 
      [coeff_base0] "r"(&a[0]), [coeff_base1] "r"(&a[64]),
      [coeff_base2] "r"(&a[128]), [coeff_base3] "r"(&a[192]),
      [coeff_base4] "r"(&a[256]), [coeff_base5] "r"(&a[320]),
      [coeff_base6] "r"(&a[384]), [coeff_base7] "r"(&a[448])
  );
}

/*************************************************
* Name:        mq_iNTT_custom_asm
*
* Description: Inverse NTT.
*              Out-of-place, constant geometry.
*              Input coefficients are in bit-reversed order and normal domain.
*              Output coefficient are smaller than Q in absolute value;
*              and output coefficients are in standard order and normal domain.
*
* Arguments:   - uint16_t p[FALCON_N]: input/output coefficient array
*
* Requirements: -VLEN=512
**************************************************/
void
mq_iNTT_custom_asm(uint16_t *a, unsigned logn)
{
    assert(logn == 9);
    size_t avl, load_zeta_vl;

    //load coefficients, from bit-reverse order to normal-order
    //v8-v15 contains offsets, v16-v31 contains coefficients
    avl = 256;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e16,m8,ta,mu\n"
    "vle16.v v8, (%[index_base0])\n"
    "vluxei16.v v16, (%[coeff_base]), v8\n"     // indexed load to v16-v23
    "vle16.v v8, (%[index_base1])\n"
    "vluxei16.v v24, (%[coeff_base]), v8\n"     // indexed load to v24-v31
    :
    : [avl]"r"(avl), 
      [index_base0]"r"(&tree_byteoffset[0]), [index_base1]"r"(&tree_byteoffset[256]),
      [coeff_base]"r"(&a[0])
  );

    //stage1, same num=1
    // (v16, v24)->v4, (v17, v25)->v6, (v18, v26)->v8, (v19, v27)->v10
    // (v20, v28)->v12, (v21, v29)->v14, (v22, v30)->v16, (v23, v31)->18 
    // with m1 are input coefficients, all have their own zetas
    // v4, v6, v8, v10, v12, v14, v16, v18 with m2 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 7));       // repeat time here should be set large enough
    avl = 32;           // m1
    load_zeta_vl = 32;  // load_zeta_vl is equal to avl
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base0])\n"
    "vbutterfly.gs.vvm v4, v16, v24, v0\n"
    "vle16.v v0, (%[zeta_base1])\n"
    "vbutterfly.gs.vvm v6, v17, v25, v0\n"
    "vle16.v v0, (%[zeta_base2])\n" 
    "vbutterfly.gs.vvm v8, v18, v26, v0\n"
    "vle16.v v0, (%[zeta_base3])\n"
    "vbutterfly.gs.vvm v10, v19, v27, v0\n"
    "vle16.v v0, (%[zeta_base4])\n"
    "vbutterfly.gs.vvm v12, v20, v28, v0\n"
    "vle16.v v0, (%[zeta_base5])\n"
    "vbutterfly.gs.vvm v14, v21, v29, v0\n"
    "vle16.v v0, (%[zeta_base6])\n"
    "vbutterfly.gs.vvm v16, v22, v30, v0\n"
    "vle16.v v0, (%[zeta_base7])\n"
    "vbutterfly.gs.vvm v18, v23, v31, v0\n"
    :
    : [zeta_vl] "r"(avl), 
      [zeta_base0] "r"(&zetas_inv_cg_inorder_512[0]), [zeta_base1] "r"(&zetas_inv_cg_inorder_512[32]), 
      [zeta_base2] "r"(&zetas_inv_cg_inorder_512[64]), [zeta_base3] "r"(&zetas_inv_cg_inorder_512[96]),
      [zeta_base4] "r"(&zetas_inv_cg_inorder_512[128]), [zeta_base5] "r"(&zetas_inv_cg_inorder_512[160]), 
      [zeta_base6] "r"(&zetas_inv_cg_inorder_512[192]), [zeta_base7] "r"(&zetas_inv_cg_inorder_512[224])
  );

    //stage2, same num=2
    // (v4, v12)->v20, (v6, v14)->v24,
    // (v8, v16)->v28, (v10, v18)->v4,
    // with m2 are input coefficients
    // v20, v24, v28, v4 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(1, 7));       // repeat time here should be set large enough
    avl = 64;           // m2
    load_zeta_vl = 32;  // avl / 2 = 32
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base0])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v20, v4, v12, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base1])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v24, v6, v14, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base2])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v28, v8, v16, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base3])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v4, v10, v18, v0\n"
    :
    : [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl),
      [zeta_base0] "r"(&zetas_inv_cg_inorder_512[256]), [zeta_base1] "r"(&zetas_inv_cg_inorder_512[288]), 
      [zeta_base2] "r"(&zetas_inv_cg_inorder_512[320]), [zeta_base3] "r"(&zetas_inv_cg_inorder_512[352])
  );    

    //stage3, same num=4
    // (v20, v28), (v22, v30), (v24, v4), (v26, v6) with m2 are input coefficients
    // v8, v12, v16, v20 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(2, 7));
    load_zeta_vl = 16;  // avl / 4 = 16
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v8, v20, v28, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v12, v22, v30, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v16, v24, v4, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v20, v26, v6, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [avl]"r"(avl),
      [zeta_base_0]"r"(&zetas_inv_cg_inorder_512[384]), [zeta_base_1]"r"(&zetas_inv_cg_inorder_512[400]),   
      [zeta_base_2]"r"(&zetas_inv_cg_inorder_512[416]), [zeta_base_3]"r"(&zetas_inv_cg_inorder_512[432])
  );

    //stage4, same num=8
    // (v8, v16), (v10, v18), (v12, v20), (v14, v22) with m2 are input coefficients
    // v24, v28, v4, v8 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(3, 7));
    load_zeta_vl = 8;  // avl / 8 = 8
     __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v24, v8, v16, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v28, v10, v18, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v4, v12, v20, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v8, v14, v22, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [avl]"r"(avl),
      [zeta_base_0]"r"(&zetas_inv_cg_inorder_512[448]), [zeta_base_1]"r"(&zetas_inv_cg_inorder_512[456]),  
      [zeta_base_2]"r"(&zetas_inv_cg_inorder_512[464]), [zeta_base_3]"r"(&zetas_inv_cg_inorder_512[472])
  );

    //stage5, same num=16
    // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m2 are input coefficients
    // v12, v16, v20, v4 with m4 are output coefficients    
    csr_zeta_selMode_rw(GET_ZETA_MODE(4, 7));
    load_zeta_vl = 4;   // avl / 16 = 4
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v12, v24, v4, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v16, v26, v6, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v20, v28, v8, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v4, v30, v10, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [avl]"r"(avl), 
      [zeta_base_0]"r"(&zetas_inv_cg_inorder_512[480]), [zeta_base_1]"r"(&zetas_inv_cg_inorder_512[484]),
      [zeta_base_2]"r"(&zetas_inv_cg_inorder_512[488]), [zeta_base_3]"r"(&zetas_inv_cg_inorder_512[492])
  );

    //stage6, same num=32
    // (v12, v20), (v14, v22), (v16, v4), (v18, v6) with m2 are input coefficients
    //  v8, v24, v28, v20 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(5, 7));
    load_zeta_vl = 2;   // avl / 32 = 2     
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v8, v12, v20, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v24, v14, v22, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v28, v16, v4, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v20, v18, v6, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [avl]"r"(avl), 
      [zeta_base_0]"r"(&zetas_inv_cg_inorder_512[496]), [zeta_base_1]"r"(&zetas_inv_cg_inorder_512[498]),
      [zeta_base_2]"r"(&zetas_inv_cg_inorder_512[500]), [zeta_base_3]"r"(&zetas_inv_cg_inorder_512[502])
  );

    // stage7, same num=64
    // (v8, v28), (v10, v30), (v24, v20), (v26, v22) with m2 are input coefficients
    // v4, v12, v16, v28 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(6, 7));
    load_zeta_vl = 1;   // avl / 64 = 1    
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v4, v8, v28, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v12, v10, v30, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_2])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v16, v24, v20, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_3])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v28, v26, v22, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [avl]"r"(avl), 
      [zeta_base_0]"r"(&zetas_inv_cg_inorder_512[504]), [zeta_base_1]"r"(&zetas_inv_cg_inorder_512[505]),
      [zeta_base_2]"r"(&zetas_inv_cg_inorder_512[506]), [zeta_base_3]"r"(&zetas_inv_cg_inorder_512[507])
  );  

    // stage8, same num=128
    // (v4, v16) and (v6, v18) with m2 are input coefficients with same zetas
    // (v12, v28) and (v14, v30) with m2 are input coefficients with same zetas
    // v8, v20, v24, v16 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(7, 7));
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_0])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v8, v4, v16, v0\n"
    "vbutterfly.gs.vvm v20, v6, v18, v0\n"
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base_1])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v24, v12, v28, v0\n"
    "vbutterfly.gs.vvm v16, v14, v30, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [avl]"r"(avl), 
      [zeta_base_0]"r"(&zetas_inv_cg_inorder_512[508]), [zeta_base_1]"r"(&zetas_inv_cg_inorder_512[509])
  );    

    // stage9, same num=256
    // (v8, v24), (v10, v26), (v20, v16), (v22, v18) with m2 are input coefficients with same zetas
    // v4, v12, v24, v28 with m4 is output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(8, 7));
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[zeta_vl],e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero,%[avl],e16,m2,ta,mu\n"
    "vbutterfly.gs.vvm v4, v8, v24, v0\n"
    "vbutterfly.gs.vvm v12, v10, v26, v0\n"
    "vbutterfly.gs.vvm v24, v20, v16, v0\n"
    "vbutterfly.gs.vvm v28, v22, v18, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [avl]"r"(avl), 
      [zeta_base]"r"(&zetas_inv_cg_inorder_512[510])
  );


    //store back, from bit-reverse order to normal-order
    //v16 with m4 used to store index, v4, v12, v24, v28 with m4 contains coefficients
    avl = 128;  // m4
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e16,m4,ta,mu\n"
    "vle16.v v16, (%[index_base0])\n"
    "vsuxei16.v v4, (%[coeff_base]), v16\n"
    "vle16.v v16, (%[index_base1])\n"
    "vsuxei16.v v12, (%[coeff_base]), v16\n"
    "vle16.v v16, (%[index_base2])\n"
    "vsuxei16.v v24, (%[coeff_base]), v16\n"
    "vle16.v v16, (%[index_base3])\n"
    "vsuxei16.v v28, (%[coeff_base]), v16\n"         
    :
    : [avl]"r"(avl), 
      [index_base0]"r"(&tree_byteoffset[0]),
      [index_base1]"r"(&tree_byteoffset[128]),
      [index_base2]"r"(&tree_byteoffset[256]),
      [index_base3]"r"(&tree_byteoffset[384]),
      [coeff_base]"r"(a)
  );

}
/*
 * Convert a polynomial (mod q) to Montgomery representation.
 */
static void
mq_poly_tomonty_custom(uint16_t *f, unsigned logn)
{
	size_t n = (size_t)1 << logn;
    for (size_t vl; n > 0; n -= vl, f += vl) {
        vl = vsetvl_e16m1(n);
        vuint16m1_t vec_f = vle16_v_u16m1(f, vl);
        vec_f = vmod_mul_vx_u16m1(vec_f, R2);
        vse16_v_u16m1(f, vec_f, vl);
    }
}


/*
 * Multiply two polynomials together (NTT representation, and using
 * a Montgomery multiplication). Result f*g is written over f.
 */
static void
mq_poly_montymul_ntt_custom(uint16_t *f, const uint16_t *g, unsigned logn)
{
	size_t n = (size_t)1 << logn;
    for (size_t vl; n > 0; n -= vl, f += vl, g += vl) {
        vl = vsetvl_e16m1(n);
        vuint16m1_t vec_f = vle16_v_u16m1(f, vl);
        vuint16m1_t vec_g = vle16_v_u16m1(g, vl);
        vec_f = vmod_mul_vv_u16m1(vec_f, vec_g);
        vse16_v_u16m1(f, vec_f, vl);
    }
}


/*
 * Subtract polynomial g from polynomial f.
 */
static void
mq_poly_sub_custom(uint16_t *f, const uint16_t *g, unsigned logn)
{
	size_t n = (size_t)1 << logn;
    for (size_t vl; n > 0; n -= vl, f += vl, g += vl) {
        vl = vsetvl_e16m1(n);
        vuint16m1_t vec_f = vle16_v_u16m1(f, vl);
        vuint16m1_t vec_g = vle16_v_u16m1(g, vl);
        // f = f - g
        vec_f = vmod_sub_vv_u16m1(vec_f, vec_g);
        vse16_v_u16m1(f, vec_f, vl);
    }
}


/* see inner.h */
void
Zf(to_ntt_monty_custom)(uint16_t *h, unsigned logn)
{
	mq_NTT_custom_asm(h, logn);
	mq_poly_tomonty_custom(h, logn);
}


/* see inner.h */
int
Zf(verify_raw_custom)(const uint16_t *c0, const int16_t *s2,
	const uint16_t *h, unsigned logn, uint8_t *tmp)
{
	size_t n;
	uint16_t *tt;

	n = (size_t)1 << logn;
	tt = (uint16_t *)tmp;

	/*
	 * Reduce s2 elements modulo q ([0..q-1] range).
	 */
    const int16_t *pi = s2;
    uint16_t *p = tt;
    for (size_t vl; n > 0; n -= vl, p += vl, pi += vl) {
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_value = vle16_v_i16m1(pi, vl);
        vuint16m1_t vec_out;
        vec_out = vmod_add_vx_u16m1_i16m1(vec_value, Q);
        vse16_v_u16m1(p, vec_out, vl);
    }

	/*
	 * Compute -s1 = s2*h - c0 mod phi mod q (in tt[]).
	 */
	mq_NTT_custom_asm(tt, logn);
	mq_poly_montymul_ntt_custom(tt, h, logn);
	mq_iNTT_custom_asm(tt, logn);
	mq_poly_sub_custom(tt, c0, logn);

	/*
	 * Normalize -s1 elements into the [-q/2..q/2] range.
	 */
	n = (size_t)1 << logn;
    p = tt;
    int16_t *ti = (int16_t *)tt;
    for (size_t vl; n > 0; n -= vl, p += vl, ti += vl) {
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_value = vle16_v_i16m1((int16_t*)p, vl);
        vint16m1_t vec_out;
        vbool16_t vmask;
        // check if the element > (Q>>1) or not
        vmask = vmsgt_vx_i16m1_b16(vec_value, Q>>1, vl);
        // if element > (Q>>1), then element -= Q
        vec_out = vsub_vx_i16m1_m(vmask, vec_value, vec_value, Q, vl);
        vse16_v_i16m1(ti, vec_out, vl);
    }

	/*
	 * Signature is valid if and only if the aggregate (-s1,s2) vector
	 * is short enough.
	 */
	return Zf(is_short_custom)((int16_t *)tt, s2, logn);
}


/* see inner.h */
int
Zf(compute_public_custom)(uint16_t *h,
	const int8_t *f, const int8_t *g, unsigned logn, uint8_t *tmp)
{
	uint16_t *tt;

	size_t n = (size_t)1 << logn;
	tt = (uint16_t *)tmp;
    uint16_t *p = tt;
    uint16_t *ph = h;

    for (size_t vl; n > 0; n -= vl, p += vl, ph += vl, f += vl, g += vl) {
        vl = vsetvl_e8mf2(n);
        vint8mf2_t vec_f = vle8_v_i8mf2(f, vl);
        vint8mf2_t vec_g = vle8_v_i8mf2(g, vl);
		vint16m1_t vec_value_f = vsext_vf2_i16m1(vec_f, vl);
		vint16m1_t vec_value_g = vsext_vf2_i16m1(vec_g, vl);
        vuint16m1_t vec_out;
        vec_out = vmod_add_vx_u16m1_i16m1(vec_value_f, Q);
        vse16_v_u16m1(p, vec_out, vl);
        vec_out = vmod_add_vx_u16m1_i16m1(vec_value_g, Q);
        vse16_v_u16m1(ph, vec_out, vl);
    }

	mq_NTT_custom_asm(h, logn);
	mq_NTT_custom_asm(tt, logn);

	n = (size_t)1 << logn;
    p = tt;
    ph = h;
    for (size_t vl; n > 0; n -= vl, p += vl, ph += vl) {
        vl = vsetvl_e16m1(n);
		vuint16m1_t vec_y0, vec_y9, vec_y10, vec_tu;
        vuint16m1_t vec_tt = vle16_v_u16m1(p, vl);
		vbool16_t mask = vmseq_vx_u16m1_b16(vec_tt, 0, vl);
		size_t num = vcpop_m_b16(mask, vl);
		if (num > 0) {
			return 0;
		}
        vec_tt = vmod_mul_vx_u16m1(vec_tt, R2);         // y0
		vec_y0 = vmv_v_v_u16m1(vec_tt, vl);
        vec_tu = vmod_mul_vv_u16m1(vec_tt, vec_tt);     // y1
        vec_tt = vmod_mul_vv_u16m1(vec_tu, vec_tt);     // y2
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tt);     // y3
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y4
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y5
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y6
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y7
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y8
        vec_y9 = vmod_mul_vv_u16m1(vec_tu, vec_tt);     // y9
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y9);     // y10
		vec_y10 = vmv_v_v_u16m1(vec_tu, vl);
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y11
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y12
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y9);     // y13    
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y14
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y15
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y10);    // y16
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y17
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y0);     // y18
        vec_tt = vle16_v_u16m1(ph, vl);
        vec_tt = vmod_mul_vv_u16m1(vec_tu, vec_tt);
        vse16_v_u16m1(ph, vec_tt, vl);
    }

	mq_iNTT_custom_asm(h, logn);

	return 1;
}


/* see inner.h */
int
Zf(complete_private_custom)(int8_t *G,
	const int8_t *f, const int8_t *g, const int8_t *F,
	unsigned logn, uint8_t *tmp)
{
	size_t n;
	uint16_t *t1, *t2;

	n = (size_t)1 << logn;
	t1 = (uint16_t *)tmp;
	t2 = t1 + n;
    uint16_t *pt1 = t1;
    uint16_t *pt2 = t2;

    for (size_t vl; n > 0; n -= vl, pt1 += vl, pt2 += vl, g += vl, F += vl) {
        vl = vsetvl_e8mf2(n);
        vint8mf2_t vec_g = vle8_v_i8mf2(g, vl);
        vint8mf2_t vec_F = vle8_v_i8mf2(F, vl);
		vint16m1_t vec_value_g = vsext_vf2_i16m1(vec_g, vl);
		vint16m1_t vec_value_F = vsext_vf2_i16m1(vec_F, vl);
        vint16m1_t vec_out;
        vec_out = vmod_add_vx_i16m1(vec_value_g, Q);
        vse16_v_i16m1((int16_t*)pt1, vec_out, vl);
        vec_out = vmod_add_vx_i16m1(vec_value_F, Q);
        vse16_v_i16m1((int16_t*)pt2, vec_out, vl);
    }

	mq_NTT_custom_asm(t1, logn);
	mq_NTT_custom_asm(t2, logn);
	mq_poly_tomonty_custom(t1, logn);
	mq_poly_montymul_ntt_custom(t1, t2, logn);

	n = (size_t)1 << logn;
    pt2 = t2;
    for (size_t vl; n > 0; n -= vl, pt2 += vl, f += vl) {
        vl = vsetvl_e8mf2(n);
        vint8mf2_t vec_f = vle8_v_i8mf2(f, vl);
		vint16m1_t vec_value_f = vsext_vf2_i16m1(vec_f, vl);
        vint16m1_t vec_out;
        vec_out = vmod_add_vx_i16m1(vec_value_f, Q);
        vse16_v_i16m1((int16_t*)pt2, vec_out, vl);
    }

	mq_NTT_custom_asm(t2, logn);

	n = (size_t)1 << logn;
    pt1 = t1;
    pt2 = t2;
    for (size_t vl; n > 0; n -= vl, pt1 += vl, pt2 += vl) {
        vl = vsetvl_e16m1(n);
		vuint16m1_t vec_y0, vec_y9, vec_y10, vec_tu;
        vuint16m1_t vec_tt = vle16_v_u16m1(pt2, vl);
		vbool16_t mask = vmseq_vx_u16m1_b16(vec_tt, 0, vl);
		size_t num = vcpop_m_b16(mask, vl);
		if (num > 0) {
			return 0;
		}
        vec_tt = vmod_mul_vx_u16m1(vec_tt, R2);         // y0
		vec_y0 = vmv_v_v_u16m1(vec_tt, vl);
        vec_tu = vmod_mul_vv_u16m1(vec_tt, vec_tt);     // y1
        vec_tt = vmod_mul_vv_u16m1(vec_tu, vec_tt);     // y2
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tt);     // y3
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y4
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y5
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y6
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y7
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y8
        vec_y9 = vmod_mul_vv_u16m1(vec_tu, vec_tt);     // y9
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y9);     // y10
		vec_y10 = vmv_v_v_u16m1(vec_tu, vl);
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y11
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y12
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y9);     // y13
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y14
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y15
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y10);    // y16
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y17
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y0);     // y18
        vec_tt = vle16_v_u16m1(pt1, vl);
        vec_tt = vmod_mul_vv_u16m1(vec_tu, vec_tt);
        vse16_v_u16m1(pt1, vec_tt, vl);
    }
	mq_iNTT_custom_asm(t1, logn);
	n = (size_t)1 << logn;
    pt1 = t1;
    for (size_t vl; n > 0; n -= vl, pt1 += vl, G += vl) {
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_value = vle16_v_i16m1((int16_t*)pt1, vl);
        vint16m1_t vec_ivalue;
        vbool16_t mask;

        // check if the element > (Q>>1) or not
        mask = vmsgt_vx_i16m1_b16(vec_value, Q>>1, vl);
        // if element > (Q>>1), then element -= Q
        vec_ivalue = vsub_vx_i16m1_m(mask, vec_value, vec_value, Q, vl);

		mask = vmsgt_vx_i16m1_b16(vec_ivalue, 127, vl);
		size_t num = vcpop_m_b16(mask, vl);
		mask = vmslt_vx_i16m1_b16(vec_ivalue, -127, vl);
		num += vcpop_m_b16(mask, vl);
		if (num > 0) {
			return 0;
		}
		vint8mf2_t vec_out = vnsra_wx_i8mf2(vec_ivalue, 0, vl);
        vse8_v_i8mf2(G, vec_out, vl);
    }
             
	return 1;
}


/* see inner.h */
int
Zf(is_invertible_custom)(
	const int16_t *s2, unsigned logn, uint8_t *tmp)
{
	size_t n;
	uint16_t *tt;
	uint32_t r;

	n = (size_t)1 << logn;
	tt = (uint16_t *)tmp;
    uint16_t *pt = tt;
    for (size_t vl; n > 0; n -= vl, pt += vl, s2 += vl) {
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_value = vle16_v_i16m1(s2, vl);
        vuint16m1_t vec_out;
        vec_out = vmod_add_vx_u16m1_i16m1(vec_value, Q);
        vse16_v_u16m1(pt, vec_out, vl);
    }
	mq_NTT_custom_asm(tt, logn);
	r = 0;
	n = (size_t)1 << logn;
    pt = tt;
    for (size_t vl; n > 0; n -= vl, pt += vl) {
        vl = vsetvl_e16m1(n);
        vuint16m1_t vec_tt = vle16_v_u16m1(pt, vl);
		vbool16_t mask = vmseq_vx_u16m1_b16(vec_tt, 0, vl);
		r += vcpop_m_b16(mask, vl);
    }

	return (int)(r == 0);
}


/* see inner.h */
int
Zf(verify_recover_custom)(uint16_t *h,
	const uint16_t *c0, const int16_t *s1, const int16_t *s2,
	unsigned logn, uint8_t *tmp)
{
	size_t n;
	uint16_t *tt;
	uint32_t r;

	n = (size_t)1 << logn;

	/*
	 * Reduce elements of s1 and s2 modulo q; then write s2 into tt[]
	 * and c0 - s1 into h[].
	 */
	tt = (uint16_t *)tmp;
    uint16_t *p = tt;
    uint16_t *ph = h;
    const int16_t *ps1 = s1;
    const int16_t *ps2 = s2;
    for (size_t vl; n > 0; n -= vl, p += vl, ph += vl, ps1 += vl, ps2 += vl, c0 += vl) {
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_s = vle16_v_i16m1(ps2, vl);
        vuint16m1_t vec_out;
        vec_out = vmod_add_vx_u16m1_i16m1(vec_s, Q);
        vse16_v_u16m1(p, vec_out, vl);
        vec_s = vle16_v_i16m1(ps1, vl);
        vec_out = vmod_add_vx_u16m1_i16m1(vec_s, Q);
        vuint16m1_t vec_c0 = vle16_v_u16m1(c0, vl);
        vec_out = vmod_sub_vv_u16m1(vec_c0, vec_out);
        vse16_v_u16m1(ph, vec_out, vl);
    }

	/*
	 * Compute h = (c0 - s1) / s2. If one of the coefficients of s2
	 * is zero (in NTT representation) then the operation fails. We
	 * keep that information into a flag so that we do not deviate
	 * from strict constant-time processing; if all coefficients of
	 * s2 are non-zero, then the high bit of r will be zero.
	 */
	mq_NTT_custom_asm(tt, logn);
	mq_NTT_custom_asm(h, logn);
	r = 0;
	n = (size_t)1 << logn;
    p = tt;
	ph = h;
    for (size_t vl; n > 0; n -= vl, p += vl, ph += vl) {
        vl = vsetvl_e16m1(n);
        vuint16m1_t vec_tt = vle16_v_u16m1(p, vl);
		vbool16_t mask = vmseq_vx_u16m1_b16(vec_tt, 0, vl);
		r += vcpop_m_b16(mask, vl);
		vuint16m1_t vec_y0, vec_y9, vec_y10, vec_tu;
        vec_tt = vmod_mul_vx_u16m1(vec_tt, R2);         // y0
		vec_y0 = vmv_v_v_u16m1(vec_tt, vl);
        vec_tu = vmod_mul_vv_u16m1(vec_tt, vec_tt);     // y1
        vec_tt = vmod_mul_vv_u16m1(vec_tu, vec_tt);     // y2
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tt);     // y3
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y4
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y5
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y6
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y7
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y8
        vec_y9 = vmod_mul_vv_u16m1(vec_tu, vec_tt);     // y9
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y9);     // y10
		vec_y10 = vmv_v_v_u16m1(vec_tu, vl);
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y11
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y12
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y9);     // y13
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y14
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y15
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y10);    // y16
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_tu);     // y17
        vec_tu = vmod_mul_vv_u16m1(vec_tu, vec_y0);     // y18
        vec_tt = vle16_v_u16m1(ph, vl);
        vec_tt = vmod_mul_vv_u16m1(vec_tu, vec_tt);
        vse16_v_u16m1(ph, vec_tt, vl);
    }

	mq_iNTT_custom_asm(h, logn);

	/*
	 * Signature is acceptable if and only if it is short enough,
	 * and s2 was invertible mod phi mod q. The caller must still
	 * check that the rebuilt public key matches the expected
	 * value (e.g. through a hash).
	 */
	return (int)((r == 0) && Zf(is_short_custom)(s1, s2, logn));
}


/* see inner.h */
int
Zf(count_nttzero_custom)(const int16_t *sig, unsigned logn, uint8_t *tmp)
{
	uint16_t *s2;
	size_t n;
	uint32_t r;

	s2 = (uint16_t *)tmp;
	n = (size_t)1 << logn;
    uint16_t *pt = s2;
    for (size_t vl; n > 0; n -= vl, pt += vl, sig += vl) {
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_value = vle16_v_i16m1(sig, vl);
        vuint16m1_t vec_out;
        vec_out = vmod_add_vx_u16m1_i16m1(vec_value, Q);
        vse16_v_u16m1(pt, vec_out, vl);
    }
	mq_NTT_custom_asm(s2, logn);
	r = 0;
	pt = s2;
    for (size_t vl; n > 0; n -= vl, pt += vl) {
        vl = vsetvl_e16m1(n);
        vuint16m1_t vec_tt = vle16_v_u16m1(pt, vl);
		vbool16_t mask = vmseq_vx_u16m1_b16(vec_tt, 0, vl);
		r += vcpop_m_b16(mask, vl);
    }

	return (int)r;
}