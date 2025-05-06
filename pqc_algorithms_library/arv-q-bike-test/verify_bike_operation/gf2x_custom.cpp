#include "gf2x_port.h"
#include "poly_op_util.h"
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"
/*****************************************
 * 
 *    Karatsuba Multiplication Related
 *    Ref: bike_kem: gf2x_mul_base_vpclmul.c
 * 
*****************************************/
//At the end of recursion we perform Karatsuba on 1024x1024bit
//Using RVV, VLEN at least 256
#if (VLEN>=512)

vuint64m2_t mul2_512_custom(IN vuint64m1_t va, IN vuint64m1_t vb){
    uint64_t buffer[8]={0};//for test
    uint64_t shuffle_idx[8]={1,0,3,2,5,4,7,6};
    vuint64m1_t vshuffle_idx=vle64_v_u64m1(shuffle_idx,8);
    vuint64m1_t va_tmp=vrgather_vv_u64m1(va,vshuffle_idx,8);
    vuint64m1_t vb_tmp=vrgather_vv_u64m1(vb,vshuffle_idx,8);
    vuint64m1_t vs1=vxor_vv_u64m1(va,va_tmp,8);
    vuint64m1_t vs2=vxor_vv_u64m1(vb,vb_tmp,8);

    uint64_t shuffle_idx1[8]={0,2,4,6,1,3,5,7};
    vuint64m1_t vshuffle_idx1=vle64_v_u64m1(shuffle_idx1,8);
    va=vrgather_vv_u64m1(va,vshuffle_idx1,8);
    vb=vrgather_vv_u64m1(vb,vshuffle_idx1,8);
    vs1=vrgather_vv_u64m1(vs1,vshuffle_idx1,8);
    vs2=vrgather_vv_u64m1(vs2,vshuffle_idx1,8);
    
    __asm__ __volatile__ ("clmul.vv v2, %[vt], %[vs]"::[vt]"vr"(va),[vs]"vr"(vb));//Occupies v2 and v3, v2 stores lq results, v3 stores hq results
    __asm__ __volatile__ ("clmul.vv v4, %[vt], %[vs]"::[vt]"vr"(vs1),[vs]"vr"(vs2));//Only v4 is actually used

    vsetvl_e64m1(8);
    vuint64m1_t vabq;
    __asm__ __volatile__ ("vxor.vv %[vd], v2, v3":[vd]"=&vr"(vabq):);// lq xor hq, do not overwrite
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v4":[vd]"=vr"(vabq):[vt]"vr"(vabq));//lq xor hq xor CLMUL(s1, s2, 0x00)
    
    vabq=vrgather_vv_u64m1(vabq,vshuffle_idx,8);//abq         = PERMXVAR_I64(mask_abq, abq);
    
    uint8_t mask=0xaa;
    vbool64_t vmask=vlm_v_b64(&mask,8);//One byte is enough
    __asm__ __volatile__ ("vsetivli	zero,8,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v2, v2, %[vt], %[vmask].t"::[vt]"vr"(vabq),[vmask]"vm"(vmask));//*l          = MXOR_I64(lq, 0xaa, lq, abq);
 
    mask=0x55;
    vmask=vlm_v_b64(&mask,8);//One byte is enough
    __asm__ __volatile__ ("vsetivli	zero,8,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v3, v3, %[vt], %[vmask].t"::[vt]"vr"(vabq),[vmask]"vm"(vmask));//*h          = MXOR_I64(hq, 0x55, hq, abq);

    __asm__ __volatile__ ("vsetivli	zero,16,e64,m2,ta,mu"::);
    vuint64m2_t vres;
    __asm__ __volatile__ ("vmv.v.v %[vd], v2":[vd]"=vr"(vres):);

    return vres;
}

// 8x8 Karatsuba multiplication
vuint64m2_t gf2x_mul8_512_int_custom(IN vuint64m1_t va, IN vuint64m1_t vb){
    uint64_t mask_s1[8]={2,3,0,1,4,5,6,7};
    uint64_t mask_s2[8]={0,1,4,5,6,7,2,3};
    uint64_t buffer[8]={0};//for test
    // Calculate:
    // AX1^AX3|| AX2^AX3 || AX0^AX2 || AX0^AX1
    // BX1^BX3|| BX2^BX3 || BX0^BX2 || BX0^BX1
    // Where (AX1^AX3 || AX0^AX2) stands for (AX1 || AX0)^(AX3 || AX2) = AY0^AY1
    vuint64m1_t vtmp1,vtmp2,vmask1,vmask2;
    vmask1=vle64_v_u64m1(mask_s1,8);
    vmask2=vle64_v_u64m1(mask_s2,8);
    vtmp1=vrgather_vv_u64m1(va,vmask1,8);
    vtmp2=vrgather_vv_u64m1(va,vmask2,8);
    vuint64m1_t vt0=vxor_vv_u64m1(vtmp1,vtmp2,8);
    vtmp1=vrgather_vv_u64m1(vb,vmask1,8);
    vtmp2=vrgather_vv_u64m1(vb,vmask2,8);
    vuint64m1_t vt1=vxor_vv_u64m1(vtmp1,vtmp2,8);

    // Calculate:
    // Don't care || AX1^AX3^AX0^AX2
    // Don't care || BX1^BX3^BX0^BX2
    uint64_t mask_tmp[8]={4,5,6,7,0,1,2,3};
    vmask1=vle64_v_u64m1(mask_tmp,8);
    vtmp1=vrgather_vv_u64m1(vt0,vmask1,8);
    vuint64m1_t vt2=vxor_vv_u64m1(vtmp1,vt0,8);
    vtmp1=vrgather_vv_u64m1(vt1,vmask1,8);
    vuint64m1_t vt3=vxor_vv_u64m1(vtmp1,vt1,8);

    vuint64m2_t vxhl=mul2_512_custom(va,vb);
    vuint64m2_t vxabhl=mul2_512_custom(vt0,vt1);
    vuint64m2_t vyabhl=mul2_512_custom(vt2,vt3);


    uint64_t mask01[16]={0,1,8,9,4,5,12,13,2,3,10,11,6,7,14,15};
    vuint64m2_t vmask_m2=vle64_v_u64m2(mask01,16);
    __asm__ __volatile__ ("vrgather.vv v18, %[vt], %[vs]"::[vt]"vr"(vxabhl),[vs]"vr"(vmask_m2));//Now v18 stores PERMX2VAR_I64(xabl, mask0, xabh), v19 stores oxh = PERMX2VAR_I64(xabl, mask1, xabh)
    __asm__ __volatile__ ("vrgather.vv v22, %[vt], %[vs]"::[vt]"vr"(vyabhl),[vs]"vr"(vmask_m2));//Now v22 stores PERMX2VAR_I64(yabl, mask0, yabh)

    __asm__ __volatile__ ("vmv.v.v v16, %[vt]"::[vt]"vr"(vxhl));//Now v16 stores xl and v17 stores xh
    __asm__ __volatile__ ("vsetivli	zero,8,e64,m1,ta,mu"::);
    vuint64m1_t vxab;
    __asm__ __volatile__ ("vxor.vv %[vd], v16, v17":[vd]"=&vr"(vxab):);//xl ^ xh
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v18":[vd]"=vr"(vxab):[vt]"vr"(vxab));//xab  = xl ^ xh ^ PERMX2VAR_I64(xabl, mask0, xabh);

    uint64_t mask34[16]={0,1,2,3,8,9,10,11,4,5,6,7,12,13,14,15};
    vmask_m2=vle64_v_u64m2(mask34,16);
    __asm__ __volatile__ ("vrgather.vv v20, v16, %[vs]"::[vs]"vr"(vmask_m2));//Now v20 stores yl, v21 stores yh

    uint64_t mask_tmp1[8]={6,7,0,1,2,3,4,5};
    vmask1=vle64_v_u64m1(mask_tmp1,8);
    vuint64m1_t vxab1=vrgather_vv_u64m1(vxab,vmask1,8);//xab1 = VALIGN(xab, xab, 6);
    uint64_t mask_tmp2[8]={2,3,4,5,6,7,0,1};
    vmask1=vle64_v_u64m1(mask_tmp2,8);
    vuint64m1_t vxab2=vrgather_vv_u64m1(vxab,vmask1,8);//xab2 = VALIGN(xab, xab, 2);

    uint8_t mask=0x3c;
    vbool64_t vmask=vlm_v_b64(&mask,8);//One byte is enough
    __asm__ __volatile__ ("vsetivli	zero,8,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v20, v20, %[vt], %[vmask].t"::[vt]"vr"(vxab1),[vmask]"vm"(vmask));//yl   = MXOR_I64(yl, 0x3c, yl, xab1);
    __asm__ __volatile__ ("vxor.vv v21, v21, %[vt], %[vmask].t"::[vt]"vr"(vxab2),[vmask]"vm"(vmask));//yh   = MXOR_I64(yh, 0x3c, yh, xab2);

    uint64_t mask_tmp3[8]={4,5,6,7,0,1,2,3};
    vmask1=vle64_v_u64m1(mask_tmp3,8);
    __asm__ __volatile__ ("vrgather.vv v18, v19, %[vs]"::[vs]"vr"(vmask1));//Now v18 stores oxl, oxl = VALIGN(oxh, oxh, 4);

    __asm__ __volatile__ ("vsetivli	zero,8,e64,m1,ta,mu"::);
    vuint64m1_t vyab;
    __asm__ __volatile__ ("vxor.vv %[vd], v18, v19":[vd]"=&vr"(vyab):);//oxl ^ oxh
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v22":[vd]"=vr"(vyab):[vt]"vr"(vyab));//yab         = oxl ^ oxh ^ PERMX2VAR_I64(yabl, mask0, yabh);
    
    vmask1=vle64_v_u64m1(mask_tmp2,8);
    vyab=vrgather_vv_u64m1(vyab,vmask1,8);
    __asm__ __volatile__ ("vxor.vv v19, v19, %[vt], %[vmask].t"::[vt]"vr"(vyab),[vmask]"vm"(vmask));//yab         = MXOR_I64(oxh, 0x3c, oxh, VALIGN(yab, yab, 2)); Now v19 stores yab
    __asm__ __volatile__ ("vxor.vv %[vd], v19, v20":[vd]"=vr"(vyab):);//yab ^ yl
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v21":[vd]"=vr"(vyab):[vt]"vr"(vyab));//yab ^ yl ^ yh
    
    vmask1=vle64_v_u64m1(mask_tmp,8);
    vyab=vrgather_vv_u64m1(vyab,vmask1,8);//yab = PERMXVAR_I64(mask2, yab);
    mask=0xf0;
    vmask=vlm_v_b64(&mask,8);
    __asm__ __volatile__ ("vsetivli	zero,8,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v20, v20, %[vs], %[vmask].t"::[vs]"vr"(vyab),[vmask]"vm"(vmask));
    mask=0x0f;
    vmask=vlm_v_b64(&mask,8);
    __asm__ __volatile__ ("vsetivli	zero,8,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v21, v21, %[vs], %[vmask].t"::[vs]"vr"(vyab),[vmask]"vm"(vmask));

    __asm__ __volatile__ ("vsetivli	zero,16,e64,m2,ta,mu"::);
    vuint64m2_t vres;
    __asm__ __volatile__ ("vmv.v.v %[vd], v20":[vd]"=vr"(vres):);

    return vres;
}

// 1024x1024 bit multiplication performed by Karatsuba algorithm.
// Here, a and b are considered as having 16 digits of size 64 bits.
void gf2x_mul_base_vpclmul_custom(OUT uint64_t *c,IN const uint64_t *a,IN const uint64_t *b){
    vuint64m1_t va0=vle64_v_u64m1(a,8);
    vuint64m1_t va1=vle64_v_u64m1(&a[8],8);
    vuint64m1_t vb0=vle64_v_u64m1(b,8);
    vuint64m1_t vb1=vle64_v_u64m1(&b[8],8);

    vuint64m2_t vlo=gf2x_mul8_512_int_custom(va0,vb0);
    vuint64m2_t vhi=gf2x_mul8_512_int_custom(va1,vb1);
    vuint64m1_t vtemp1=vxor_vv_u64m1(va0,va1,8);
    vuint64m1_t vtemp2=vxor_vv_u64m1(vb0,vb1,8);
    vuint64m2_t vmi=gf2x_mul8_512_int_custom(vtemp1,vtemp2);

    __asm__ __volatile__ ("vsetivli	zero,16,e64,m2,ta,mu"::);
    __asm__ __volatile__ ("vmv.v.v v16, %[vt]"::[vt]"vr"(vlo));//Now v16 stores vlo_0, v17 stores vlo_1
    __asm__ __volatile__ ("vmv.v.v v18, %[vt]"::[vt]"vr"(vhi));//Now v18 stores vhi_0, v19 stores vhi_1
    __asm__ __volatile__ ("vmv.v.v v20, %[vt]"::[vt]"vr"(vmi));//Now v20 stroes vmi_0, v21 stores vmi_1

    __asm__ __volatile__ ("vsetivli	zero,8,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v22, v17, v18"::);//__m512i m = lo[1] ^ hi[0];Now v22 stores m

    __asm__ __volatile__ ("vse64.v v16, (%[addr])"::[addr]"r"(&c[0]));//STORE(&c[0 * QWORDS_IN_ZMM], lo[0]);


    __asm__ __volatile__ ("vxor.vv %[vd], v20, v16":[vd]"=&vr"(vtemp1):);//mi[0] ^ lo[0]
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v22":[vd]"=vr"(vtemp1):[vt]"vr"(vtemp1));//mi[0] ^ lo[0] ^ m
    vse64_v_u64m1(&c[8],vtemp1,8);//STORE(&c[1 * QWORDS_IN_ZMM], mi[0] ^ lo[0] ^ m);

    __asm__ __volatile__ ("vxor.vv %[vd], v19, v21":[vd]"=&vr"(vtemp1):);//mi[1] ^ hi[1]
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v22":[vd]"=vr"(vtemp1):[vt]"vr"(vtemp1));//mi[1] ^ hi[1] ^ m
    vse64_v_u64m1(&c[16],vtemp1,8);//STORE(&c[2 * QWORDS_IN_ZMM], mi[1] ^ hi[1] ^ m);

    __asm__ __volatile__ ("vse64.v v19, (%[addr])"::[addr]"r"(&c[24]));//STORE(&c[3 * QWORDS_IN_ZMM], hi[1]);
}

#elif (VLEN==256)
vuint64m2_t gf2x_mul_256_int_custom(IN vuint64m1_t va, IN vuint64m1_t vb){
    //va=(a3,a2,a1,a0) vb=(b3,b2,b1,b0)
    uint64_t shuffle_idx[4]={1,0,3,2};
    vuint64m1_t vshuffle_idx=vle64_v_u64m1(shuffle_idx,4);
    vuint64m1_t va_tmp=vrgather_vv_u64m1(va,vshuffle_idx,4);
    vuint64m1_t vb_tmp=vrgather_vv_u64m1(vb,vshuffle_idx,4);
    vuint64m1_t vs1=vxor_vv_u64m1(va,va_tmp,4);
    vuint64m1_t vs2=vxor_vv_u64m1(vb,vb_tmp,4);

    uint64_t shuffle_idx2[4]={2,3,0,1};
    vuint64m1_t vshuffle_idx2=vle64_v_u64m1(shuffle_idx2,4);
    va_tmp=vrgather_vv_u64m1(va,vshuffle_idx2,4);
    va_tmp=vxor_vv_u64m1(va,va_tmp,4);
    vb_tmp=vrgather_vv_u64m1(vb,vshuffle_idx2,4);
    vb_tmp=vxor_vv_u64m1(vb,vb_tmp,4);

    uint64_t shuffle_idx1[4]={0,2,1,3};
    vuint64m1_t vshuffle_idx1=vle64_v_u64m1(shuffle_idx1,4);
    va=vrgather_vv_u64m1(va,vshuffle_idx1,4);
    vb=vrgather_vv_u64m1(vb,vshuffle_idx1,4);
    vs1=vrgather_vv_u64m1(vs1,vshuffle_idx1,4);
    vs2=vrgather_vv_u64m1(vs2,vshuffle_idx1,4);
    
    __asm__ __volatile__ ("clmul.vv v2, %[vt], %[vs]"::[vt]"vr"(va),[vs]"vr"(vb));//Occupies v2 and v3, v2 stores lq results, v3 stores hq results
    __asm__ __volatile__ ("clmul.vv v4, %[vt], %[vs]"::[vt]"vr"(vs1),[vs]"vr"(vs2));//Only v4 is actually used

    vsetvl_e64m1(4);
    vuint64m1_t vabq;
    __asm__ __volatile__ ("vxor.vv %[vd], v2, v3":[vd]"=&vr"(vabq):);// lq xor hq, do not overwrite
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v4":[vd]"=vr"(vabq):[vt]"vr"(vabq));//lq xor hq xor CLMUL(s1, s2, 0x00)
    
    vabq=vrgather_vv_u64m1(vabq,vshuffle_idx,4);//abq         = PERMXVAR_I64(mask_abq, abq);
    
    uint8_t mask=0xa;
    vbool64_t vmask=vlm_v_b64(&mask,4);//One byte is enough
    __asm__ __volatile__ ("vsetivli	zero,4,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v2, v2, %[vt], %[vmask].t"::[vt]"vr"(vabq),[vmask]"vm"(vmask));//*l          = MXOR_I64(lq, 0xaa, lq, abq);
 
    mask=0x5;
    vmask=vlm_v_b64(&mask,4);//One byte is enough
    __asm__ __volatile__ ("vsetivli	zero,4,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v3, v3, %[vt], %[vmask].t"::[vt]"vr"(vabq),[vmask]"vm"(vmask));//*h          = MXOR_I64(hq, 0x55, hq, abq);

    //Now v2 stores (a1,a0)*(b1,b0), v3 stores (a3,a2)*(b3,b2)

    /********************************************************
     * //treat va_tmp vb_tmp as va and vb again
    *********************************************************/
    va=vrgather_vv_u64m1(va_tmp,vshuffle_idx,4);
    vb=vrgather_vv_u64m1(vb_tmp,vshuffle_idx,4);
    vs1=vxor_vv_u64m1(va,va_tmp,4);
    vs2=vxor_vv_u64m1(vb,vb_tmp,4);

    va=vrgather_vv_u64m1(va_tmp,vshuffle_idx1,4);
    vb=vrgather_vv_u64m1(vb_tmp,vshuffle_idx1,4);
    vs1=vrgather_vv_u64m1(vs1,vshuffle_idx1,4);
    vs2=vrgather_vv_u64m1(vs2,vshuffle_idx1,4);

    __asm__ __volatile__ ("clmul.vv v6, %[vt], %[vs]"::[vt]"vr"(va),[vs]"vr"(vb));//Occupies v6 and v7, v6 stores lq results, v7 stores hq results
    __asm__ __volatile__ ("clmul.vv v8, %[vt], %[vs]"::[vt]"vr"(vs1),[vs]"vr"(vs2));//Only v8 is actually used

    vsetvl_e64m1(4);
    __asm__ __volatile__ ("vxor.vv %[vd], v6, v7":[vd]"=&vr"(vabq):);// lq xor hq, do not overwrite
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v8":[vd]"=vr"(vabq):[vt]"vr"(vabq));//lq xor hq xor CLMUL(s1, s2, 0x00)

    vabq=vrgather_vv_u64m1(vabq,vshuffle_idx,4);//abq         = PERMXVAR_I64(mask_abq, abq);
    
    mask=0xa;
    vmask=vlm_v_b64(&mask,4);//One byte is enough
    __asm__ __volatile__ ("vsetivli	zero,4,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v6, v6, %[vt], %[vmask].t"::[vt]"vr"(vabq),[vmask]"vm"(vmask));//*l          = MXOR_I64(lq, 0xaa, lq, abq);
 
    mask=0x5;
    vmask=vlm_v_b64(&mask,4);//One byte is enough
    __asm__ __volatile__ ("vsetivli	zero,4,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v7, v7, %[vt], %[vmask].t"::[vt]"vr"(vabq),[vmask]"vm"(vmask));//*h          = MXOR_I64(hq, 0x55, hq, abq);

    mask=0xc;
    vmask=vlm_v_b64(&mask,4);
    __asm__ __volatile__ ("vsetivli	zero,4,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v2, v6, %[vt], %[vmask].t"::[vt]"vr"(vabq),[vmask]"vm"(vmask));

    mask=0x3;
    vmask=vlm_v_b64(&mask,4);
    __asm__ __volatile__ ("vsetivli	zero,4,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v3, v6, %[vt], %[vmask].t"::[vt]"vr"(vabq),[vmask]"vm"(vmask));

    __asm__ __volatile__ ("vsetivli	zero,16,e64,m2,ta,mu"::);
    vuint64m2_t vres;
    __asm__ __volatile__ ("vmv.v.v %[vd], v2":[vd]"=vr"(vres):);

    return vres;
}

vuint64m4_t gf2x_mul8_512_int_custom(IN vuint64m1_t va_0, IN vuint64m1_t va_1, IN vuint64m1_t vb_0, IN vuint64m1_t vb_1){
    vuint64m1_t vtemp1=vxor_vv_u64m1(va_0,va_1,4);
    vuint64m1_t vtemp2=vxor_vv_u64m1(vb_0,vb_1,4);
    vuint64m2_t va0b0=gf2x_mul_256_int_custom(va_0,vb_0);
    vuint64m2_t va1b1=gf2x_mul_256_int_custom(va_1,vb_1);
    vuint64m2_t vmi=gf2x_mul_256_int_custom(vtemp1,vtemp2);
    vmi=vxor_vv_u64m2(vmi,va0b0,8);
    vmi=vxor_vv_u64m2(vmi,va1b1,8);

    __asm__ __volatile__ ("vsetivli	zero,8,e64,m2,ta,mu"::);
    __asm__ __volatile__ ("vmv.v.v v4, %[vt]"::[vt]"vr"(va0b0));//Now v4 stores va0b0_l, v5 stores va0b0_h
    __asm__ __volatile__ ("vmv.v.v v6, %[vt]"::[vt]"vr"(va1b1));//Now v6 stores va1b1_l, v7 stores va1b1_l
    __asm__ __volatile__ ("vmv.v.v v8, %[vt]"::[vt]"vr"(vmi));//Now v8 stores vmi_l, v9 stores vmi_h

    __asm__ __volatile__ ("vsetivli	zero,4,e64,m1,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v5, v5, v8"::);
    __asm__ __volatile__ ("vxor.vv v6, v6, v9"::);
    vuint64m4_t vres;
    __asm__ __volatile__ ("vmv.v.v %[vres], v4":[vres]"=vr"(vres):);
    return vres;
}


// 1024x1024 bit multiplication performed by Karatsuba algorithm.
// Here, a and b are considered as having 16 digits of size 64 bits.
void gf2x_mul_base_vpclmul_custom(OUT uint64_t *c,IN const uint64_t *a,IN const uint64_t *b){
    vuint64m1_t va0_0=vle64_v_u64m1(&a[0],4); //256bit
    vuint64m1_t va0_1=vle64_v_u64m1(&a[4],4);
    vuint64m1_t va1_0=vle64_v_u64m1(&a[8],4);
    vuint64m1_t va1_1=vle64_v_u64m1(&a[12],4);

    vuint64m1_t vb0_0=vle64_v_u64m1(&b[0],4);
    vuint64m1_t vb0_1=vle64_v_u64m1(&b[4],4);
    vuint64m1_t vb1_0=vle64_v_u64m1(&b[8],4);
    vuint64m1_t vb1_1=vle64_v_u64m1(&b[12],4);

    vuint64m4_t vlo=gf2x_mul8_512_int_custom(va0_0,va0_1,vb0_0,vb0_1); //1024bit
    vuint64m4_t vhi=gf2x_mul8_512_int_custom(va1_0,va1_1,vb1_0,vb1_1);
    vuint64m1_t vtemp1_0=vxor_vv_u64m1(va0_0,va1_0,4);
    vuint64m1_t vtemp1_1=vxor_vv_u64m1(va0_1,va1_1,4);
    vuint64m1_t vtemp2_0=vxor_vv_u64m1(vb0_0,vb1_0,4);
    vuint64m1_t vtemp2_1=vxor_vv_u64m1(vb0_1,vb1_1,4);
    vuint64m4_t vmi=gf2x_mul8_512_int_custom(vtemp1_0,vtemp1_1,vtemp2_0,vtemp2_1);
    vuint64m2_t temp3;

    __asm__ __volatile__ ("vsetivli	zero,16,e64,m4,ta,mu"::);
    __asm__ __volatile__ ("vmv.v.v v16, %[vt]"::[vt]"vr"(vlo));//Now v16 v17 stores vlo_0, v18 v19 stores vlo_1
    __asm__ __volatile__ ("vmv.v.v v20, %[vt]"::[vt]"vr"(vhi));//Now v20 v21 stores vhi_0, v22 v23 stores vhi_1
    __asm__ __volatile__ ("vmv.v.v v24, %[vt]"::[vt]"vr"(vmi));//Now v24 v25 stroes vmi_0, v26 v27 stores vmi_1

    __asm__ __volatile__ ("vsetivli	zero,8,e64,m2,ta,mu"::);
    __asm__ __volatile__ ("vxor.vv v28, v18, v20"::);//__m512i m = lo[1] ^ hi[0];Now v28 stores m

    __asm__ __volatile__ ("vse64.v v16, (%[addr])"::[addr]"r"(&c[0]));//STORE(&c[0 * QWORDS_IN_ZMM], lo[0]);


    __asm__ __volatile__ ("vxor.vv %[vd], v24, v16":[vd]"=&vr"(temp3):);//mi[0] ^ lo[0]
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v28":[vd]"=vr"(temp3):[vt]"vr"(temp3));//mi[0] ^ lo[0] ^ m
    vse64_v_u64m2(&c[8],temp3,8);//STORE(&c[1 * QWORDS_IN_ZMM], mi[0] ^ lo[0] ^ m);

    __asm__ __volatile__ ("vxor.vv %[vd], v26, v22":[vd]"=&vr"(temp3):);//mi[1] ^ hi[1]
    __asm__ __volatile__ ("vxor.vv %[vd], %[vt], v28":[vd]"=vr"(temp3):[vt]"vr"(temp3));//mi[1] ^ hi[1] ^ m
    vse64_v_u64m2(&c[16],temp3,8);//STORE(&c[2 * QWORDS_IN_ZMM], mi[1] ^ hi[1] ^ m);

    __asm__ __volatile__ ("vse64.v v22, (%[addr])"::[addr]"r"(&c[24]));//STORE(&c[3 * QWORDS_IN_ZMM], hi[1]);
}
#endif

// The secure buffer size required for Karatsuba is computed by:
//    size(n) = 3*n/2 + size(n/2) = 3*sum_{i}{n/2^i} < 3n
#define SECURE_BUFFER_QWORDS (3 * R_PADDED_QWORDS)

void karatzuba_add1_custom(OUT uint64_t *alah,
                         OUT uint64_t *blbh,
                         IN const uint64_t *a,
                         IN const uint64_t *b,
                         IN const size_t    qwords_len)
{
    size_t vl=0;
    int32_t avl=qwords_len;
    const uint64_t* ain_addr1=a;
    const uint64_t* ain_addr2=a+qwords_len;
    const uint64_t* bin_addr1=b;
    const uint64_t* bin_addr2=b+qwords_len;
    uint64_t* aout_addr=alah;
    uint64_t* bout_addr=blbh;

    vuint64m8_t v0,v1;

    while(avl>0){
        vl=vsetvl_e64m8(avl);
        v0=vle64_v_u64m8(ain_addr1,vl);
        v1=vle64_v_u64m8(ain_addr2,vl);
        v0=vxor_vv_u64m8(v0,v1,vl);
        vse64_v_u64m8(aout_addr,v0,vl);
        ain_addr1+=vl;
        ain_addr2+=vl;
        aout_addr+=vl;

        v0=vle64_v_u64m8(bin_addr1,vl);
        v1=vle64_v_u64m8(bin_addr2,vl);
        v0=vxor_vv_u64m8(v0,v1,vl);
        vse64_v_u64m8(bout_addr,v0,vl);
        bin_addr1+=vl;
        bin_addr2+=vl;
        bout_addr+=vl;

        avl-=vl;
    }
}

void karatzuba_add2_custom(OUT uint64_t *z,
                         IN const uint64_t *x,
                         IN const uint64_t *y,
                         IN const size_t    qwords_len)
{
    size_t vl=0;
    int32_t avl=qwords_len;
    const uint64_t* in_addr1=x;
    const uint64_t* in_addr2=y;
    uint64_t* out_addr=z;

    vuint64m8_t v0,v1;

    while(avl>0){
        vl=vsetvl_e64m8(avl);
        v0=vle64_v_u64m8(in_addr1,vl);
        v1=vle64_v_u64m8(in_addr2,vl);
        v0=vxor_vv_u64m8(v0,v1,vl);
        vse64_v_u64m8(out_addr,v0,vl);
        in_addr1+=vl;
        in_addr2+=vl;
        out_addr+=vl;
        avl-=vl;
    }
}

//Should not use u64m8 here, will cause register file overflow and result in page table fault in GEM5
void karatzuba_add3_custom(OUT uint64_t *c,
                         IN const uint64_t *mid,
                         IN const size_t    qwords_len)
{
    uint64_t* c0_addr=c;
    uint64_t* c1_addr=c+qwords_len;
    uint64_t* c2_addr=c+2*qwords_len;
    uint64_t* c3_addr=c+3*qwords_len;
    const uint64_t* mid_addr=mid;
    size_t vl=0;
    int32_t avl=qwords_len;
    vuint64m2_t vr0,vr1,vr2,vr3,vt;

    while(avl>0){
        vl=vsetvl_e64m2(avl);
        vr0=vle64_v_u64m2(c0_addr,vl);
        vr1=vle64_v_u64m2(c1_addr,vl);
        vr2=vle64_v_u64m2(c2_addr,vl);
        vr3=vle64_v_u64m2(c3_addr,vl);
        vt=vle64_v_u64m2(mid_addr,vl);
        vr0=vxor_vv_u64m2(vr0,vr1,vl);
        vr0=vxor_vv_u64m2(vr0,vt,vl);
        vr2=vxor_vv_u64m2(vr2,vr3,vl);
        vr2=vxor_vv_u64m2(vr2,vt,vl);
        vse64_v_u64m2(c1_addr,vr0,vl);
        vse64_v_u64m2(c2_addr,vr2,vl);

        c0_addr+=vl;
        c1_addr+=vl;
        c2_addr+=vl;
        c3_addr+=vl;
        mid_addr+=vl;
        avl-=vl;
    }
}

_INLINE_ void karatzuba_custom(OUT uint64_t* c,
    IN const uint64_t* a,
    IN const uint64_t* b,
    IN const size_t    qwords_len,
    IN const size_t    qwords_len_pad,
    uint64_t* sec_buf)
{
    // if (qwords_len <= 1) {
    //     gf2x_mul_base_port(c, a, b);//Portable Version
    //     return;
    // }

    if (qwords_len <= 16) {
        gf2x_mul_base_vpclmul_custom(c, a, b);//RVV Version
        return;
    }

    const size_t half_qw_len = qwords_len_pad >> 1;

    // Split a and b into low and high parts of size n_padded/2
    const uint64_t* a_lo = a;
    const uint64_t* b_lo = b;
    const uint64_t* a_hi = &a[half_qw_len];
    const uint64_t* b_hi = &b[half_qw_len];

    // Split c into 4 parts of size n_padded/2 (the last ptr is not needed)
    uint64_t* c0 = c;
    uint64_t* c1 = &c[half_qw_len];
    uint64_t* c2 = &c[half_qw_len * 2];

    // Allocate 3 ptrs of size n_padded/2  on sec_buf
    uint64_t* alah = sec_buf;
    uint64_t* blbh = &sec_buf[half_qw_len];
    uint64_t* tmp = &sec_buf[half_qw_len * 2];

    // Move sec_buf ptr to the first free location for the next recursion call
    sec_buf = &sec_buf[half_qw_len * 3];

    // Compute a_lo*b_lo and store the result in (c1|c0)
    karatzuba_custom(c0, a_lo, b_lo, half_qw_len, half_qw_len, sec_buf);

    // If the real number of digits n is less or equal to n_padded/2 then:
    //     a_hi = 0 and b_hi = 0
    // and
    //     (a_hi|a_lo)*(b_hi|b_lo) = a_lo*b_lo
    // so we can skip the remaining two multiplications
    if (qwords_len > half_qw_len) {
        // Compute a_hi*b_hi and store the result in (c3|c2)
        karatzuba_custom(c2, a_hi, b_hi, qwords_len - half_qw_len, half_qw_len, sec_buf);

        // Compute alah = (a_lo + a_hi) and blbh = (b_lo + b_hi)
        karatzuba_add1_custom(alah, blbh, a, b, half_qw_len);

        // Compute (c1 + c2) and store the result in tmp
        karatzuba_add2_custom(tmp, c1, c2, half_qw_len);

        // Compute alah*blbh and store the result in (c2|c1)
        karatzuba_custom(c1, alah, blbh, half_qw_len, half_qw_len, sec_buf);

        // Add (tmp|tmp) and (c3|c0) to (c2|c1)
        karatzuba_add3_custom(c0, tmp, half_qw_len);
    }
}

/*****************************************
 * 
 *    GF2X ITI INVERSION Related
 * 
*****************************************/

// c = a^2
void gf2x_sqr_custom(OUT dbl_pad_r_t* c, IN const pad_r_t* a){
    const uint64_t* a64 = (const uint64_t*)a;
    uint64_t* c64 = (uint64_t*)c;

    int32_t avl=R_QWORDS;
    size_t vl;
    while(avl>0){
        vl=vsetvl_e64m4(avl);
        vuint64m4_t va=vle64_v_u64m4(a64,vl);
        vuint64m8_t vc;
        __asm__ __volatile__ ("clmul.vv %[vd],%[vt],%[vs]":[vd]"=&vr"(vc):[vt]"vr"(va),[vs]"vr"(va));
        size_t tmp_vl=vl<<1;
        vse64_v_u64m8(c64,vc,tmp_vl);
        //address update
        a64+=vl,c64+=tmp_vl,avl-=vl;
    }
}

// c = a mod (x^r - 1)
void gf2x_red_custom(OUT pad_r_t* c, IN const dbl_pad_r_t* a){
    uint8_t temp[R_SIZE]={0};
    poly_shiftleft_custom(&(a->raw[R_SIZE]),(8-(R_BITS&7)),R_SIZE,a->raw[R_SIZE-1],temp);
    const uint8_t* a_addr=&a->raw[0];
    uint8_t* temp_addr=temp;
    uint8_t* c_addr=&c->val.raw[0];
    int32_t avl=R_SIZE;
    size_t vl;
    while(avl>0){
        vl=vsetvl_e8m8(avl);
        vuint8m8_t va=vle8_v_u8m8(a_addr,vl);
        vuint8m8_t vtemp=vle8_v_u8m8(temp_addr,vl);
        vuint8m8_t vc=vxor_vv_u8m8(va,vtemp,vl);
        vse8_v_u8m8(c_addr,vc,vl);
        //address update
        a_addr+=vl,temp_addr+=vl,c_addr+=vl,avl-=vl;
    }
    c->val.raw[R_SIZE-1]&=LAST_R_BYTE_MASK;
}

// a = a^2 mod (x^r - 1)
_INLINE_ void gf2x_mod_sqr_in_place_custom(IN OUT pad_r_t* a,
    OUT dbl_pad_r_t* secure_buffer)
{
    gf2x_sqr_custom(secure_buffer, a);
    gf2x_red_custom(a, secure_buffer);
}

// c = a^2^2^num_sqrs
_INLINE_ void repeated_squaring_custom(OUT pad_r_t* c,
    IN pad_r_t* a,
    IN const size_t num_sqrs,
    OUT dbl_pad_r_t* sec_buf)
{
    c->val = a->val;

    for (size_t i = 0; i < num_sqrs; i++) {
        gf2x_mod_sqr_in_place_custom(c, sec_buf);
    }
}

void gf2x_mod_mul_with_ctx_custom(OUT pad_r_t* c,
    IN const pad_r_t* a,
    IN const pad_r_t* b)
{
    static_assert((R_PADDED_BYTES % 2 == 0), "Karatsuba is odd");

    dbl_pad_r_t t = { 0 };
    uint64_t secure_buffer[SECURE_BUFFER_QWORDS];

    karatzuba_custom((uint64_t*)&t, (const uint64_t*)a, (const uint64_t*)b, R_QWORDS,
        R_PADDED_QWORDS, secure_buffer);

    gf2x_red_custom(c, &t);

    memset((uint8_t*)secure_buffer,0, sizeof(secure_buffer));
}

_INLINE_ void generate_map_custom(OUT uint16_t *map, IN const uint16_t l_param){
    // The permutation map is generated in the following way:
    //   1. for i = 0 to map size:
    //   2.  map[i] = (i * l_param) % r
    // However, to avoid the expensive multiplication and modulo operations
    // we modify the algorithm to:
    //   1. map[0] = l_param
    //   2. for i = 1 to map size:
    //   3.   map[i] = map[i - 1] + l_param
    //   4.   if map[i] >= r:
    //   5.     map[i] = map[i] - r
    // This algorithm is parallelized with vector instructions by processing
    // certain number of values (NUM_OF_VALS) in parallel. Therefore,
    // in the beginning we need to initialize the first NUM_OF_VALS elements.

    size_t vl=vsetvl_e16m4(256);//Obtained in VLEN=1024bit
    for(size_t i = 0; i < vl; i++) {
        map[i] = (i * l_param) % R_BITS;
    }

    vuint16m4_t vinc=vmv_v_x_u16m4((l_param * vl) % R_BITS,vl);
    vuint16m4_t vmap=vle16_v_u16m4(map,vl);

    int32_t avl=R_BITS-vl;
    uint16_t* addr=map+vl;
    while(avl>0){
        vl=vsetvl_e16m4(avl);
        vmap=vadd_vv_u16m4(vmap,vinc,vl);
        vbool4_t vmask=vmsgeu_vx_u16m4_b4(vmap,R_BITS,vl);
        vmap=vsub_vx_u16m4_m(vmask,vmap,vmap,R_BITS,vl);
        vse16_v_u16m4(addr,vmap,vl);

        addr+=vl,avl-=vl;
    }
}

// The k-squaring function computes c = a^(2^k) % (x^r - 1),
// By [1](Observation 1), if
//     a = sum_{j in supp(a)} x^j,
// then
//     a^(2^k) % (x^r - 1) = sum_{j in supp(a)} x^((j * 2^k) % r).
// Therefore, k-squaring can be computed as permutation of the bits of "a":
//     pi0 : j --> (j * 2^k) % r.
// For improved performance, we compute the result by inverted permutation pi1:
//     pi1 : (j * 2^-k) % r --> j.
// Input argument l_param is defined as the value (2^-k) % r.
void k_sqr_custom(OUT pad_r_t *c, IN const pad_r_t *a, IN const size_t l_param){
    uint16_t map[R_BITS]={0};
    uint8_t a_bytes[R_BITS]={0};
    uint8_t c_bytes[R_BITS]={0};

    generate_map_custom(map,l_param);

    convertByteToBinary_custom(a_bytes,&(a->val.raw[0]),R_BITS);

    // Permute "a" using the generated permutation map.
    for(size_t i = 0; i < R_BITS; i++) {
        c_bytes[i] = a_bytes[map[i]];
    }

    convertBinaryToByte_custom(&(c->val.raw[0]),c_bytes,R_BITS);
}

// Inversion in F_2[x]/(x^R - 1), [1](Algorithm 2).
// c = a^{-1} mod x^r-1
void gf2x_mod_inv_custom(OUT pad_r_t *c, IN const pad_r_t *a)
{
  // Note that exp0/1_k/l are predefined constants that depend only on the value
  // of R. This value is public. Therefore, branches in this function, which
  // depends on R, are also "public". Code that releases these branches
  // (taken/not-taken) does not leak secret information.
  const size_t exp0_k[MAX_I] = {EXP0_K_VALS};
  const size_t exp0_l[MAX_I] = {EXP0_L_VALS};
  const size_t exp1_k[MAX_I] = {EXP1_K_VALS};
  const size_t exp1_l[MAX_I] = {EXP1_L_VALS};

  pad_r_t f = {0};
  pad_r_t g = {0};
  pad_r_t t = {0};
  dbl_pad_r_t sec_buf = {0};

  // Steps 2 and 3 in [1](Algorithm 2)
  f.val = a->val;
  t.val = a->val;

  for(size_t i = 1; i < MAX_I; i++) {
    // Step 5 in [1](Algorithm 2), exponentiation 0: g = f^2^2^(i-1)
    if(exp0_k[i - 1] <= K_SQR_THR) {
      repeated_squaring_custom(&g, &f, exp0_k[i - 1], &sec_buf);
    } else {
      k_sqr_custom(&g, &f, exp0_l[i - 1]);
    }

    // Step 6, [1](Algorithm 2): f = f*g
    gf2x_mod_mul_with_ctx_custom(&f, &g, &f);

    if(exp1_k[i] != 0) {
      // Step 8, [1](Algorithm 2), exponentiation 1: g = f^2^((r-2) % 2^i)
      if(exp1_k[i] <= K_SQR_THR) {
        repeated_squaring_custom(&g, &f, exp1_k[i], &sec_buf);
      } else {
        k_sqr_custom(&g, &f, exp1_l[i]);
      }

      // Step 9, [1](Algorithm 2): t = t*g;
      gf2x_mod_mul_with_ctx_custom(&t, &g, &t);
    }
  }

  // Step 10, [1](Algorithm 2): c = t^2
  gf2x_mod_sqr_in_place_custom(&t, &sec_buf);
  c->val = t.val;
}

void gf2x_mod_inv_wrapper_custom(OUT uint8_t res_bin[R_SIZE], IN const uint8_t a_bin[R_SIZE]) {
    pad_r_t a_bin_padded = { 0 };
    memcpy(a_bin_padded.val.raw,a_bin,R_SIZE);
    pad_r_t res_bin_padded = { 0 };
    gf2x_mod_inv_custom(&res_bin_padded,&a_bin_padded);
    memcpy(res_bin,res_bin_padded.val.raw,R_SIZE);
}