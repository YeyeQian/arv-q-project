#include <stdio.h>
#include <stdint.h>
#include "ntt.h"
#include "reduce.h"
#include "poly_custom_rorg_asm.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

#if (VLEN == 256)

#elif (VLEN == 512)

/*************************************************
* Name:        poly_eta_pack_caddq_ntt_asm
*
* Description: 1. pack a polynomial P with coefficients in [-ETA,ETA] into sk
*              2. Modular add q to each coefficients in polynomial P
*              in this way turn each coefficients to non-negetive number
*              3. Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in normal order  
*
* Arguments:   - poly *r:             pointer to both input and output polynomial
*              - uint8_t *packed_addr: pointer to position to store the packed data
**************************************************/
void poly_eta_pack_caddq_ntt_asm(poly* r, uint8_t* packed_addr){
    size_t load_zeta_vl,avl_m1,avl_m2,avl_m4,packed_vl,packed_vlx2,actual_width;
    uint8_t* packed_addr0;
    uint8_t* packed_addr1;
    uint8_t* packed_addr2;
    uint8_t* packed_addr3;
    int32_t* coeff_ptr=r->coeffs;
    packed_addr0=packed_addr;
    avl_m1=16,avl_m2=32,avl_m4=64;
#if ETA==2
    actual_width=3;
    packed_addr1=packed_addr0+6;
    packed_addr2=packed_addr0+48;
    packed_addr3=packed_addr0+54;
    packed_vl=6;
    packed_vlx2=12;
#elif ETA==4
    actual_width=4;
    packed_addr1=packed_addr0+8;
    packed_addr2=packed_addr0+64;
    packed_addr3=packed_addr0+72;
    packed_vl=8;
    packed_vlx2=16;
#endif
    csr_modulusq_rw(Q);
    csr_qinv_rw(QINV_HW);
    /* 
    ** Firstly, pack process+mod add q+ first stage of NTT
    ** 1. Load part of polynomial into v4 and v6 with m2
    ** 2. Pack v4,v5 v6,v7 with m1 into appropriate position with respect to packed_addr, v3 used to store temporary data
    ** 3. Mod add q on v4 and v6 with m2
    ** 4. Compute the first stage of NTT
    ** 5. Finally, v16, v20, v24, v28 with m4 stores the result of NTT after first stage
    */
    //elements 0-31 in v4 and v5, elements 128-159 in v6 and v7
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
    load_zeta_vl=1;
        //Load polynomial and perform pack
    __asm__ __volatile__ ( 
    "vsetvli zero, %[avl_m2], e32, m2, ta, mu\n"
    "vle32.v v4, (%[top_base])\n"
    "vle32.v v6, (%[btm_base])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v4, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_0])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v5, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_1])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v6, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_2])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v7, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_3])\n"
    :
    : [avl_m2]"r"(avl_m2), [top_base]"r"(&coeff_ptr[0]), [btm_base]"r"(&coeff_ptr[128]), [avl_m1]"r"(avl_m1), [eta]"r"(ETA), 
      [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3)
  );
        //Mod add q on v4 with m4
    __asm__ __volatile__ ( 
    "vsetvli zero, %[avl_m4], e32, m4, ta, mu\n"
    "vmod.add.vx v4, v4, %[q]\n"
    :
    : [avl_m4]"r"(avl_m4), [q]"r"(Q)
  );
        //Stage1 NTT
    __asm__ __volatile__ ( 
    "vsetvli zero, %[zeta_vl], e32, m1, ta, mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[1]), [avl_m2]"r"(avl_m2)
  );
    packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

    //elements 32-63 in v4 and v5, elements 160-191 in v6 and v7
        //Load polynomial and perform pack
    __asm__ __volatile__ ( 
    "vle32.v v4, (%[top_base])\n"
    "vle32.v v6, (%[btm_base])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v4, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_0])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v5, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_1])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v6, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_2])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v7, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_3])\n"
    :
    : [avl_m2]"r"(avl_m2), [top_base]"r"(&coeff_ptr[32]), [btm_base]"r"(&coeff_ptr[160]), [avl_m1]"r"(avl_m1), [eta]"r"(ETA), 
      [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3)
  );
        //Mod add q on v4 with m4
    __asm__ __volatile__ ( 
    "vsetvli zero, %[avl_m4], e32, m4, ta, mu\n"
    "vmod.add.vx v4, v4, %[q]\n"
    :
    : [avl_m4]"r"(avl_m4), [q]"r"(Q)
  );
        //Stage1 NTT
    __asm__ __volatile__ ( 
    "vsetvli zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    :
    : [avl_m2]"r"(avl_m2)
  );
    packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

    //elements 64-95 in v4 and v5, elements 192-223 in v6 and v7
        //Load polynomial and perform pack
    __asm__ __volatile__ ( 
    "vle32.v v4, (%[top_base])\n"
    "vle32.v v6, (%[btm_base])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v4, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_0])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v5, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_1])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v6, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_2])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v7, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_3])\n"
    :
    : [avl_m2]"r"(avl_m2), [top_base]"r"(&coeff_ptr[64]), [btm_base]"r"(&coeff_ptr[192]), [avl_m1]"r"(avl_m1), [eta]"r"(ETA), 
      [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3)
  );
        //Mod add q on v4 with m4
    __asm__ __volatile__ ( 
    "vsetvli zero, %[avl_m4], e32, m4, ta, mu\n"
    "vmod.add.vx v4, v4, %[q]\n"
    :
    : [avl_m4]"r"(avl_m4), [q]"r"(Q)
  );
        //Stage1 NTT
    __asm__ __volatile__ ( 
    "vsetvli zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    :
    : [avl_m2]"r"(avl_m2)
  );
    packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

    //elements 96-127 in v4 and v5, elements 224-255 in v6 and v7
        //Load polynomial and perform pack
    __asm__ __volatile__ ( 
    "vle32.v v4, (%[top_base])\n"
    "vle32.v v6, (%[btm_base])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v4, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_0])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v5, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_1])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v6, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_2])\n"
    "vsetvli zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v7, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr_3])\n"
    :
    : [avl_m2]"r"(avl_m2), [top_base]"r"(&coeff_ptr[96]), [btm_base]"r"(&coeff_ptr[224]), [avl_m1]"r"(avl_m1), [eta]"r"(ETA), 
      [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3)
  );
        //Mod add q on v4 with m4
    __asm__ __volatile__ ( 
    "vsetvli zero, %[avl_m4], e32, m4, ta, mu\n"
    "vmod.add.vx v4, v4, %[q]\n"
    :
    : [avl_m4]"r"(avl_m4), [q]"r"(Q)
  );
        //Stage1 NTT
    __asm__ __volatile__ ( 
    "vsetvli zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [avl_m2]"r"(avl_m2)
  );
    
    /* 
    ** Secondly, performing the remaining stages of NTT
    */
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[2]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[4]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[8]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[16]), [avl]"r"(avl_m2)
  );

    //stage6, repeat bound=32
    // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m1 are input coefficients, with same zetas
    // (v25, v5), (v27, v7), (v29, v9), (v31, v11) with m1 are input coefficients, with same zetas
    // v12, v14, v16, v18, v20, v22, v2, v4 with m2 are output coefficients
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[32]), [zeta_base1] "r"(&zetas_pos[48])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[64]), [zeta_base1] "r"(&zetas_pos[80]), [zeta_base2] "r"(&zetas_pos[96]), [zeta_base3] "r"(&zetas_pos[112])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[128]), [zeta_base1] "r"(&zetas_pos[144]), [zeta_base2] "r"(&zetas_pos[160]), [zeta_base3] "r"(&zetas_pos[176]),
      [zeta_base4] "r"(&zetas_pos[192]), [zeta_base5] "r"(&zetas_pos[208]), [zeta_base6] "r"(&zetas_pos[224]), [zeta_base7] "r"(&zetas_pos[240])
  );

    //store back, from bit-reverse order to normal-order
    //v8-v15 used to store index, v16-v31 contains coefficients
    uint32_t* tree_ptr=NULL;
    size_t avl=128;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m8,ta,mu\n"
    "vle32.v v8, (%[index_base0])\n"
    "vsuxei32.v v16, (%[coeff_base]), v8\n"     // indexed store v16-v23
    "vle32.v v8, (%[index_base1])\n"
    "vsuxei32.v v24, (%[coeff_base]), v8\n"     // indexed store v24-v31
    :
    : [avl]"r"(avl), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&coeff_ptr[0]), [index_base1]"r"(&tree_byteoffset[128])
  );
}

/*************************************************
* Name:        polyvecl_macc_invntt_mon2nor_asm
*
* Description: 1. Pointwise multiply two polyvecl and accumulate into a polynomial
*              2. Invntt on the polynomial
*              3. Mon2nor on the polynomial  
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const polyvecl *a:   pointer to input polyvecl 
*              - const polyvecl *b:   pointer to input polyvecl
**************************************************/
void polyvecl_macc_invntt_mon2nor_asm(poly* r, const polyvecl *a, const polyvecl *b){
  size_t load_zeta_vl,avl_m1,avl_m2,avl_m4,avl_m8;
  avl_m1=16,avl_m2=32,avl_m4=64,avl_m8=128;
  const int32_t* coeff_ptr_a=a->vec[0].coeffs;
  const int32_t* coeff_ptr_b=b->vec[0].coeffs;

  csr_modulusq_rw(Q);
  csr_qinv_rw(QINV_HW);

  /*
  ** Firstly, compute polynomial multiplication of a->vec[0] and b->vec[0]
  ** Result stored in v16-v31, which form the accumulation base
  */
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl_m8],e32,m8,ta,mu\n"
    "vle32.v v8, (%[a_ptr0])\n"
    "vle32.v v16, (%[b_ptr0])\n"
    "vmod.mul.vv v16, v16, v8\n"
    "vle32.v v8, (%[a_ptr1])\n"
    "vle32.v v24, (%[b_ptr1])\n"
    "vmod.mul.vv v24, v24, v8\n"
    :
    : [avl_m8]"r"(avl_m8), [a_ptr0]"r"(&coeff_ptr_a[0]), [b_ptr0]"r"(&coeff_ptr_b[0]), 
      [a_ptr1]"r"(&coeff_ptr_a[128]), [b_ptr1]"r"(&coeff_ptr_b[128])
  );

  /*
  ** Secondly, compute polynomial multiplication of a->vec[i] and b->vec[i], i = 1, 2, ..., L-1
  ** 1. Load partial of a->vec[i] into v4 and partial of b->vec[i] into v8 with m4 set
  ** 2. Mod mul between v4 and v8 and results stored in v8 with m4
  ** 3. Mod add v8 with corresponding part in v16-v31
  ** 4. Finally, v16-v31 stores the results of macc
  */
  for(int i=1;i<L;i++){
    coeff_ptr_a=a->vec[i].coeffs;
    coeff_ptr_b=b->vec[i].coeffs;
    __asm__ __volatile__ (
    "vsetvli	zero,%[avl_m4],e32,m4,ta,mu\n"
    "vle32.v v4, (%[a_ptr0])\n"
    "vle32.v v8, (%[b_ptr0])\n"
    "vmod.mul.vv v8, v8, v4\n"
    "vmod.add.vv v16,v16,v8\n"
    "vle32.v v4, (%[a_ptr1])\n"
    "vle32.v v8, (%[b_ptr1])\n"
    "vmod.mul.vv v8, v8, v4\n"
    "vmod.add.vv v20,v20,v8\n"
    "vle32.v v4, (%[a_ptr2])\n"
    "vle32.v v8, (%[b_ptr2])\n"
    "vmod.mul.vv v8, v8, v4\n"
    "vmod.add.vv v24,v24,v8\n"
    "vle32.v v4, (%[a_ptr3])\n"
    "vle32.v v8, (%[b_ptr3])\n"
    "vmod.mul.vv v8, v8, v4\n"
    "vmod.add.vv v28,v28,v8\n"
    :
    : [avl_m4]"r"(avl_m4), [a_ptr0]"r"(&coeff_ptr_a[0]), [b_ptr0]"r"(&coeff_ptr_b[0]), 
      [a_ptr1]"r"(&coeff_ptr_a[64]), [b_ptr1]"r"(&coeff_ptr_b[64]), 
      [a_ptr2]"r"(&coeff_ptr_a[128]), [b_ptr2]"r"(&coeff_ptr_b[128]), 
      [a_ptr3]"r"(&coeff_ptr_a[192]), [b_ptr3]"r"(&coeff_ptr_b[192])
  );
  }

  /*
  ** Thirdly, compute invntt on coefficients in v16-v31
  */

  //stage1, same num=1
  //v16-v31 with m1 stores input coefficients
  //v4, v6, v8, v10, v12, v14, v16, v18 with m2 stores output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 7));//repeat time here should be set large enough
  load_zeta_vl=16;
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl_m1],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vbutterfly.gs.vvm v4, v16, v24, v0\n"
    "vle32.v v0, (%[zeta_base_1])\n"
    "vbutterfly.gs.vvm v6, v17, v25, v0\n"
    "vle32.v v0, (%[zeta_base_2])\n"
    "vbutterfly.gs.vvm v8, v18, v26, v0\n"
    "vle32.v v0, (%[zeta_base_3])\n"
    "vbutterfly.gs.vvm v10, v19, v27, v0\n"
    "vle32.v v0, (%[zeta_base_4])\n"
    "vbutterfly.gs.vvm v12, v20, v28, v0\n"
    "vle32.v v0, (%[zeta_base_5])\n"
    "vbutterfly.gs.vvm v14, v21, v29, v0\n"
    "vle32.v v0, (%[zeta_base_6])\n"
    "vbutterfly.gs.vvm v16, v22, v30, v0\n"
    "vle32.v v0, (%[zeta_base_7])\n"
    "vbutterfly.gs.vvm v18, v23, v31, v0\n"
    :
    : [avl_m1]"r"(avl_m1), [zeta_base_0]"r"(&zetas_inv_in_order[0]), [zeta_base_1]"r"(&zetas_inv_in_order[16]),
      [zeta_base_2]"r"(&zetas_inv_in_order[32]), [zeta_base_3]"r"(&zetas_inv_in_order[48]), 
      [zeta_base_4]"r"(&zetas_inv_in_order[64]), [zeta_base_5]"r"(&zetas_inv_in_order[80]),
      [zeta_base_6]"r"(&zetas_inv_in_order[96]), [zeta_base_7]"r"(&zetas_inv_in_order[112])
  );

  //stage2, same num=2
  // (v4, v12), (v6, v14), (v8, v16), (v10, v18) with m2 are input coefficients
  // v20, v24, v28, v4 with m4 is output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(1, 7));
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
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[128]), [avl]"r"(avl_m2), 
      [zeta_base_1]"r"(&zetas_inv_in_order[144]), [zeta_base_2]"r"(&zetas_inv_in_order[160]), [zeta_base_3]"r"(&zetas_inv_in_order[176])
  );

  //stage3, same num=4
  // (v20, v28), (v22, v30), (v24, v4), (v26, v6) with m2 are input coefficients
  // v8, v12, v16, v20 is output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(2, 7));
  load_zeta_vl=8;
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
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[192]), [avl]"r"(avl_m2), 
    [zeta_base_1]"r"(&zetas_inv_in_order[200]), [zeta_base_2]"r"(&zetas_inv_in_order[208]), [zeta_base_3]"r"(&zetas_inv_in_order[216])
);

  //stage4, same num=8
  // (v8, v16), (v10, v18), (v12, v20), (v14, v22) with m2 are input coefficients
  // v24, v28, v4, v8 is output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(3, 7));
  load_zeta_vl=4;
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
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[224]), [avl]"r"(avl_m2), 
    [zeta_base_1]"r"(&zetas_inv_in_order[228]), [zeta_base_2]"r"(&zetas_inv_in_order[232]), [zeta_base_3]"r"(&zetas_inv_in_order[236])
);

  //stage5, same num=16
  // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m2 are input coefficients
  // v12, v16, v20, v24 is output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(4, 7));
  load_zeta_vl=2;
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
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[240]), [avl]"r"(avl_m2), 
    [zeta_base_1]"r"(&zetas_inv_in_order[242]), [zeta_base_2]"r"(&zetas_inv_in_order[244]), [zeta_base_3]"r"(&zetas_inv_in_order[246])
);

  //stage6, same num=32
  // (v12, v20), (v14, v22), (v16, v24), (v18, v26) with m2 are input coefficients
  //  v28, v4, v8, v12 with m4 are output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(5, 7));
  load_zeta_vl=1;
  __asm__ __volatile__ ( 
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_0])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v28, v12, v20, v0\n"
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_1])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v4, v14, v22, v0\n"
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_2])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v8, v16, v24, v0\n"
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_3])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v12, v18, v26, v0\n"
  :
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[248]), [avl]"r"(avl_m2), 
    [zeta_base_1]"r"(&zetas_inv_in_order[249]), [zeta_base_2]"r"(&zetas_inv_in_order[250]), [zeta_base_3]"r"(&zetas_inv_in_order[251])
);

  //stage7, same num=64
  // (v28, v8), (v30, v10), (v4, v12), (v6, v14) with m2 are input coefficients
  // v16, v20, v24, v28 with m4 are output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(6, 7));
  __asm__ __volatile__ ( 
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_0])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v16, v28, v8, v0\n"
  "vbutterfly.gs.vvm v20, v30, v10, v0\n"
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_1])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v24, v4, v12, v0\n"
  "vbutterfly.gs.vvm v28, v6, v14, v0\n"
  :
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[252]), [avl]"r"(avl_m2), [zeta_base_1]"r"(&zetas_inv_in_order[253])
);

  //stage8, same num=128
  // (v16, v24), (v18, v26), (v20, v28), (v22, v30) with m2 are input coefficients
  // v4, v8, v12, v16 with m4 are output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(7, 7));
  __asm__ __volatile__ ( 
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_0])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v4, v16, v24, v0\n"
  "vbutterfly.gs.vvm v8, v18, v26, v0\n"
  "vbutterfly.gs.vvm v12, v20, v28, v0\n"
  "vbutterfly.gs.vvm v16, v22, v30, v0\n"
  :
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[254]), [avl]"r"(avl_m2)
);

  /*
  ** Fourthly, compute mon2nor on coefficients in v4-v19 with m4
  */
  int32_t cons=2365951;
  __asm__ __volatile__ ( 
  "vsetvli	zero,%[avl],e32,m4,ta,mu\n"
  "vmod.mul.vx v4,v4,%[cons]\n"
  "vmod.mul.vx v8,v8,%[cons]\n"
  "vmod.mul.vx v12,v12,%[cons]\n"
  "vmod.mul.vx v16,v16,%[cons]\n"
  :
  : [avl]"r"(avl_m4), [cons]"r"(cons)
);

  //store back, back, from bit-reverse order to normal-order
  __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m4,ta,mu\n"
    "vle32.v v28, (%[index_base0])\n"
    "vsuxei32.v v4, (%[coeff_base]), v28\n"     
    "vle32.v v28, (%[index_base1])\n"
    "vsuxei32.v v8, (%[coeff_base]), v28\n" 
    "vle32.v v28, (%[index_base2])\n"
    "vsuxei32.v v12, (%[coeff_base]), v28\n"
    "vle32.v v28, (%[index_base3])\n"
    "vsuxei32.v v16, (%[coeff_base]), v28\n"    
    :
    : [avl]"r"(avl_m4), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&(r->coeffs[0])), [index_base1]"r"(&tree_byteoffset[64]),
      [index_base2]"r"(&tree_byteoffset[128]), [index_base3]"r"(&tree_byteoffset[192])
  );
}

/*************************************************
* Name:        poly_eta_pack_caddq_asm
*
* Description: 1. pack a polynomial P with coefficients in [-ETA,ETA] into sk
*              2. Modular add q to each coefficients in polynomial P
*              in this way turn each coefficients to non-negetive number
*
* Arguments:   - poly *r:             pointer to both input and output polynomial
*              - uint8_t *packed_addr: pointer to position to store the packed data
**************************************************/
void poly_eta_pack_caddq_asm(poly* r, uint8_t* packed_addr){
  size_t packed_vl,avl_m1,avl_m8,actual_width;
  avl_m1=16,avl_m8=128;
  int32_t* coeff_ptr=r->coeffs;
#if ETA==2
  actual_width=3;
  packed_vl=6;
#elif ETA==4
  actual_width=4;
  packed_vl=8;
#endif
  csr_modulusq_rw(Q);
  csr_qinv_rw(QINV_HW);
  /*
  ** 1.Load coefficients of r into v16-v31 with m8
  ** 2.Pack v16-v31 with m1 into v3 and store v3 back with respect packed_addr
  ** 3.Mod add q on v16-v31 with m8
  */
  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m8], e32, m8, ta, mu\n"
    "vle32.v v16, (%[r_base0])\n"
    "vle32.v v24, (%[r_base1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v16, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m8]"r"(avl_m8), [r_base0]"r"(&coeff_ptr[0]), [r_base1]"r"(&coeff_ptr[128]), 
      [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v17, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v18, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v19, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v20, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v21, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v22, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v23, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v24, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v25, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v26, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v27, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v28, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v29, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v30, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );
  packed_addr+=packed_vl;

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "vrsub.vx v3, v31, %[eta]\n"
    "pack.vx v3, v3, %[actual_width]\n"
    "vsetvli zero, %[packed_vl], e8, m1, ta, mu\n"
    "vse8.v v3, (%[packed_addr])\n"
    :
    : [avl_m1]"r"(avl_m1), [eta]"r"(ETA), [actual_width]"r"(actual_width), [packed_vl]"r"(packed_vl), [packed_addr]"r"(packed_addr)
  );

  __asm__ __volatile__ (
    "vsetvli	zero, %[avl_m8], e32, m8, ta, mu\n"
    "vmod.add.vx v16, v16, %[q]\n"
    "vmod.add.vx v24, v24, %[q]\n"
    "vse32.v v16, (%[r_base0])\n"
    "vse32.v v24, (%[r_base1])\n"
    :
    : [avl_m8]"r"(avl_m8), [q]"r"(Q), [r_base0]"r"(&coeff_ptr[0]), [r_base1]"r"(&coeff_ptr[128])
  );
}

/*************************************************
* Name:        poly_eta_unpack_caddq_ntt_asm
*
* Description: 1. Unpack a polynomial P with coefficients in [-ETA,ETA] from sk
*              2. Modular add q to each coefficients in polynomial P
*              in this way turn each coefficients to non-negetive number
*              3. Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in normal order  
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const uint8_t *packed_addr: pointer to position of bytes to start the unpack
**************************************************/
void poly_eta_unpack_caddq_ntt_asm(poly* r, const uint8_t* packed_addr){
  size_t load_zeta_vl,avl_m1,avl_m2,avl_m4,packed_vl,packed_vlx2,actual_width;
    uint8_t* packed_addr0;
    uint8_t* packed_addr1;
    uint8_t* packed_addr2;
    uint8_t* packed_addr3;
    int32_t* coeff_ptr=r->coeffs;
    packed_addr0=packed_addr;
    avl_m1=16,avl_m2=32,avl_m4=64;
#if ETA==2
    actual_width=3;
    packed_addr1=packed_addr0+6;
    packed_addr2=packed_addr0+48;
    packed_addr3=packed_addr0+54;
    packed_vl=6;
    packed_vlx2=12;
#elif ETA==4
    actual_width=4;
    packed_addr1=packed_addr0+8;
    packed_addr2=packed_addr0+64;
    packed_addr3=packed_addr0+72;
    packed_vl=8;
    packed_vlx2=16;
#endif
    csr_modulusq_rw(Q);
    csr_qinv_rw(QINV_HW);
  /* 
    ** Firstly, Unpack process+mod add q+ first stage of NTT
    ** 1. Load part of packed data into v3 and gradually unpack into v4,v5,v6,v7 with m1
    ** 3. Mod add q on v4 with m4
    ** 4. Compute the first stage of NTT
    ** 5. Finally, v16, v20, v24, v28 with m4 stores the result of NTT after first stage
    */
  //elements 0-31 in v4 and v5, elements 128-159 in v6 and v7
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  load_zeta_vl=1;
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[eta]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli zero, %[zeta_vl], e32, m1, ta, mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [eta]"r"(ETA), [q]"r"(Q), [zeta_vl]"r"(load_zeta_vl), [avl_m2]"r"(avl_m2), [zeta_base]"r"(&zetas_pos[1])
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;
  
  //elements 32-63 in v4 and v5, elements 160-191 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[eta]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [eta]"r"(ETA), [q]"r"(Q), [avl_m2]"r"(avl_m2)
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 64-95 in v4 and v5, elements 192-223 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[eta]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [eta]"r"(ETA), [q]"r"(Q), [avl_m2]"r"(avl_m2)
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 96-127 in v4 and v5, elements 224-255 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[eta]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [eta]"r"(ETA), [q]"r"(Q), [avl_m2]"r"(avl_m2)
  );

  /* 
  ** Secondly, performing the remaining stages of NTT
  */
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[2]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[4]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[8]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[16]), [avl]"r"(avl_m2)
  );

    //stage6, repeat bound=32
    // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m1 are input coefficients, with same zetas
    // (v25, v5), (v27, v7), (v29, v9), (v31, v11) with m1 are input coefficients, with same zetas
    // v12, v14, v16, v18, v20, v22, v2, v4 with m2 are output coefficients
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[32]), [zeta_base1] "r"(&zetas_pos[48])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[64]), [zeta_base1] "r"(&zetas_pos[80]), [zeta_base2] "r"(&zetas_pos[96]), [zeta_base3] "r"(&zetas_pos[112])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[128]), [zeta_base1] "r"(&zetas_pos[144]), [zeta_base2] "r"(&zetas_pos[160]), [zeta_base3] "r"(&zetas_pos[176]),
      [zeta_base4] "r"(&zetas_pos[192]), [zeta_base5] "r"(&zetas_pos[208]), [zeta_base6] "r"(&zetas_pos[224]), [zeta_base7] "r"(&zetas_pos[240])
  );

    //store back, from bit-reverse order to normal-order
    //v8-v15 used to store index, v16-v31 contains coefficients
    uint32_t* tree_ptr=NULL;
    size_t avl=128;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m8,ta,mu\n"
    "vle32.v v8, (%[index_base0])\n"
    "vsuxei32.v v16, (%[coeff_base]), v8\n"     // indexed store v16-v23
    "vle32.v v8, (%[index_base1])\n"
    "vsuxei32.v v24, (%[coeff_base]), v8\n"     // indexed store v24-v31
    :
    : [avl]"r"(avl), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&coeff_ptr[0]), [index_base1]"r"(&tree_byteoffset[128])
  );
}

/*************************************************
* Name:        poly_t0_unpack_caddq_ntt_asm
*
* Description: 1. Unpack a polynomial P with coefficients in [-2^{D-1}, 2^{D-1}] from sk
*              2. Modular add q to each coefficients in polynomial P
*              in this way turn each coefficients to non-negetive number
*              3. Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in normal order  
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const uint8_t *packed_addr: pointer to position of bytes to start the unpack
**************************************************/
void poly_t0_unpack_caddq_ntt_asm(poly* r, const uint8_t* packed_addr){
  size_t load_zeta_vl,avl_m1,avl_m2,avl_m4,packed_vl,packed_vlx2,actual_width;
  uint8_t* packed_addr0;
  uint8_t* packed_addr1;
  uint8_t* packed_addr2;
  uint8_t* packed_addr3;
  int32_t* coeff_ptr=r->coeffs;
  packed_addr0=packed_addr;
  avl_m1=16,avl_m2=32,avl_m4=64;
  int32_t cons=1 << (D-1);

  actual_width=13;
  packed_addr1=packed_addr0+26;
  packed_addr2=packed_addr0+208;
  packed_addr3=packed_addr0+234;
  packed_vl=26;
  packed_vlx2=52;

  csr_modulusq_rw(Q);
  csr_qinv_rw(QINV_HW);
  /* 
    ** Firstly, Unpack process+mod add q+ first stage of NTT
    ** 1. Load part of packed data into v3 and gradually unpack into v4,v5,v6,v7 with m1
    ** 3. Mod add q on v4 with m4
    ** 4. Compute the first stage of NTT
    ** 5. Finally, v16, v20, v24, v28 with m4 stores the result of NTT after first stage
    */
  //elements 0-31 in v4 and v5, elements 128-159 in v6 and v7
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  load_zeta_vl=1;
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[cons]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli zero, %[zeta_vl], e32, m1, ta, mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [cons]"r"(cons), [q]"r"(Q), [zeta_vl]"r"(load_zeta_vl), [avl_m2]"r"(avl_m2), [zeta_base]"r"(&zetas_pos[1])
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;
  
  //elements 32-63 in v4 and v5, elements 160-191 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[cons]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [cons]"r"(cons), [q]"r"(Q), [avl_m2]"r"(avl_m2)
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 64-95 in v4 and v5, elements 192-223 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[cons]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [cons]"r"(cons), [q]"r"(Q), [avl_m2]"r"(avl_m2)
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 96-127 in v4 and v5, elements 224-255 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[cons]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [cons]"r"(cons), [q]"r"(Q), [avl_m2]"r"(avl_m2)
  );

  /* 
  ** Secondly, performing the remaining stages of NTT
  */
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[2]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[4]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[8]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[16]), [avl]"r"(avl_m2)
  );

    //stage6, repeat bound=32
    // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m1 are input coefficients, with same zetas
    // (v25, v5), (v27, v7), (v29, v9), (v31, v11) with m1 are input coefficients, with same zetas
    // v12, v14, v16, v18, v20, v22, v2, v4 with m2 are output coefficients
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[32]), [zeta_base1] "r"(&zetas_pos[48])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[64]), [zeta_base1] "r"(&zetas_pos[80]), [zeta_base2] "r"(&zetas_pos[96]), [zeta_base3] "r"(&zetas_pos[112])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[128]), [zeta_base1] "r"(&zetas_pos[144]), [zeta_base2] "r"(&zetas_pos[160]), [zeta_base3] "r"(&zetas_pos[176]),
      [zeta_base4] "r"(&zetas_pos[192]), [zeta_base5] "r"(&zetas_pos[208]), [zeta_base6] "r"(&zetas_pos[224]), [zeta_base7] "r"(&zetas_pos[240])
  );

    //store back, from bit-reverse order to normal-order
    //v8-v15 used to store index, v16-v31 contains coefficients
    uint32_t* tree_ptr=NULL;
    size_t avl=128;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m8,ta,mu\n"
    "vle32.v v8, (%[index_base0])\n"
    "vsuxei32.v v16, (%[coeff_base]), v8\n"     // indexed store v16-v23
    "vle32.v v8, (%[index_base1])\n"
    "vsuxei32.v v24, (%[coeff_base]), v8\n"     // indexed store v24-v31
    :
    : [avl]"r"(avl), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&coeff_ptr[0]), [index_base1]"r"(&tree_byteoffset[128])
  );
}

/*************************************************
* Name:        poly_gamma1_sample_caddq_ntt_asm
*
* Description: 1. Use SHA3 to fill buf
*              2. Perform polyz_unpack to get polynomial with coefficients in [-(GAMMA1 - 1), GAMMA1]
*              3. Mod Add q on all coefficients
*              4. Store the coefficients into y, consider the coefficients in vreg to come from z
*              5. Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in normal order  
*
* Arguments:   - poly *y:             pointer to output polynomial y
*              - poly *z:             pointer to output polynomial z
*              - const uint8_t seed[CRHBYTES]: byte array with seed of length CRHBYTES
*              - uint16_t nonce:      16-bit nonce
**************************************************/
#if GAMMA1 == (1 << 17)
#define POLY_UNIFORM_GAMMA1_NBLOCKS ((576 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#elif GAMMA1 == (1 << 19)
#define POLY_UNIFORM_GAMMA1_NBLOCKS ((640 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
#endif
void poly_gamma1_sample_caddq_ntt_asm(poly* y, poly* z, const uint8_t seed[CRHBYTES], uint16_t nonce){
  /*
  ** Firstly, SHA3 fill the buf
  */
  uint8_t buf[POLY_UNIFORM_GAMMA1_NBLOCKS*STREAM256_BLOCKBYTES];

  stream256_init_custom(seed, nonce);
  stream256_squeezeblocks_custom(buf, POLY_UNIFORM_GAMMA1_NBLOCKS, false);

  size_t load_zeta_vl,avl_m1,avl_m2,avl_m4,packed_vl,packed_vlx2,actual_width;
  uint8_t* packed_addr0;
  uint8_t* packed_addr1;
  uint8_t* packed_addr2;
  uint8_t* packed_addr3;
  int32_t* coeff_ptr_y=y->coeffs;
  int32_t* coeff_ptr_z=z->coeffs;
  packed_addr0=buf;
  avl_m1=16,avl_m2=32,avl_m4=64;

#if GAMMA1 == (1 << 17)
    actual_width=18;
    packed_addr1=packed_addr0+36;
    packed_addr2=packed_addr0+288;
    packed_addr3=packed_addr0+324;
    packed_vl=36;
    packed_vlx2=72;
#elif GAMMA1 == (1 << 19)
    actual_width=20;
    packed_addr1=packed_addr0+40;
    packed_addr2=packed_addr0+320;
    packed_addr3=packed_addr0+360;
    packed_vl=40;
    packed_vlx2=80;
#endif

  csr_modulusq_rw(Q);
  csr_qinv_rw(QINV_HW);

  /* 
  ** Secondly, Unpack gamma1 process+mod add q+ first stage of NTT
  ** 1. Load part of packed data in buf into v3 and gradually unpack into v4,v5,v6,v7 with m1
  ** 2. Mod add q on v4 with m4
  ** 3. Store v4 and v6 with m2 into polynomial y
  ** 4. Compute the first stage of NTT
  ** 5. Finally, v16, v20, v24, v28 with m4 stores the result of NTT after first stage
  */
  //elements 0-31 in v4 and v5, elements 128-159 in v6 and v7
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  load_zeta_vl=1;
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[cons]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vse32.v v4, (%[y_tp_base])\n"
    "vse32.v v6, (%[y_btm_base])\n"
    "vsetvli zero, %[zeta_vl], e32, m1, ta, mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [cons]"r"(GAMMA1), [q]"r"(Q), [zeta_vl]"r"(load_zeta_vl), [avl_m2]"r"(avl_m2), [zeta_base]"r"(&zetas_pos[1]),
      [y_tp_base]"r"(&coeff_ptr_y[0]), [y_btm_base]"r"(&coeff_ptr_y[128])
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 32-63 in v4 and v5, elements 160-191 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[cons]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vse32.v v4, (%[y_tp_base])\n"
    "vse32.v v6, (%[y_btm_base])\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [cons]"r"(GAMMA1), [q]"r"(Q), [avl_m2]"r"(avl_m2), [y_tp_base]"r"(&coeff_ptr_y[32]), [y_btm_base]"r"(&coeff_ptr_y[160])
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 64-95 in v4 and v5, elements 192-223 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[cons]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vse32.v v4, (%[y_tp_base])\n"
    "vse32.v v6, (%[y_btm_base])\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [cons]"r"(GAMMA1), [q]"r"(Q), [avl_m2]"r"(avl_m2), [y_tp_base]"r"(&coeff_ptr_y[64]), [y_btm_base]"r"(&coeff_ptr_y[192])
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 96-127 in v4 and v5, elements 224-255 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vrsub.vx v4, v4, %[cons]\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vse32.v v4, (%[y_tp_base])\n"
    "vse32.v v6, (%[y_btm_base])\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [cons]"r"(GAMMA1), [q]"r"(Q), [avl_m2]"r"(avl_m2), [y_tp_base]"r"(&coeff_ptr_y[96]), [y_btm_base]"r"(&coeff_ptr_y[224])
  );


  /* 
  ** Secondly, performing the remaining stages of NTT
  */
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[2]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[4]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[8]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[16]), [avl]"r"(avl_m2)
  );

    //stage6, repeat bound=32
    // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m1 are input coefficients, with same zetas
    // (v25, v5), (v27, v7), (v29, v9), (v31, v11) with m1 are input coefficients, with same zetas
    // v12, v14, v16, v18, v20, v22, v2, v4 with m2 are output coefficients
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[32]), [zeta_base1] "r"(&zetas_pos[48])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[64]), [zeta_base1] "r"(&zetas_pos[80]), [zeta_base2] "r"(&zetas_pos[96]), [zeta_base3] "r"(&zetas_pos[112])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[128]), [zeta_base1] "r"(&zetas_pos[144]), [zeta_base2] "r"(&zetas_pos[160]), [zeta_base3] "r"(&zetas_pos[176]),
      [zeta_base4] "r"(&zetas_pos[192]), [zeta_base5] "r"(&zetas_pos[208]), [zeta_base6] "r"(&zetas_pos[224]), [zeta_base7] "r"(&zetas_pos[240])
  );

    //store back, from bit-reverse order to normal-order
    //v8-v15 used to store index, v16-v31 contains coefficients
    uint32_t* tree_ptr=NULL;
    size_t avl=128;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m8,ta,mu\n"
    "vle32.v v8, (%[index_base0])\n"
    "vsuxei32.v v16, (%[coeff_base]), v8\n"     // indexed store v16-v23
    "vle32.v v8, (%[index_base1])\n"
    "vsuxei32.v v24, (%[coeff_base]), v8\n"     // indexed store v24-v31
    :
    : [avl]"r"(avl), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&coeff_ptr_z[0]), [index_base1]"r"(&tree_byteoffset[128])
  );
}

/*************************************************
* Name:        poly_caddq_ntt_asm
*
* Description: 1. Mod Add q on all coefficients
*              2. Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in normal order  
*
* Arguments:   - poly *r:             pointer to both input and output polynomial
**************************************************/
void poly_caddq_ntt_asm(poly* r){
    size_t load_zeta_vl;
    size_t avl=32;//m2
    int32_t* coeff_ptr=r->coeffs;

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
    "vmod.add.vx v4, v4, %[q]\n"
    "vmod.add.vx v6, v6, %[q]\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    "vle32.v v4, (%[top_base_1])\n"
    "vle32.v v6, (%[btm_base_1])\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vmod.add.vx v6, v6, %[q]\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    "vle32.v v4, (%[top_base_2])\n"
    "vle32.v v6, (%[btm_base_2])\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vmod.add.vx v6, v6, %[q]\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    "vle32.v v4, (%[top_base_3])\n"
    "vle32.v v6, (%[btm_base_3])\n"
    "vmod.add.vx v4, v4, %[q]\n"
    "vmod.add.vx v6, v6, %[q]\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[1]), [avl]"r"(avl), [q]"r"(Q),
      [top_base_0]"r"(&coeff_ptr[0]), [btm_base_0]"r"(&coeff_ptr[128]), [top_base_1]"r"(&coeff_ptr[32]), [btm_base_1]"r"(&coeff_ptr[160]),
      [top_base_2]"r"(&coeff_ptr[64]), [btm_base_2]"r"(&coeff_ptr[192]), [top_base_3]"r"(&coeff_ptr[96]), [btm_base_3]"r"(&coeff_ptr[224])
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[2]), [avl]"r"(avl)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[4]), [avl]"r"(avl)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[8]), [avl]"r"(avl)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[16]), [avl]"r"(avl)
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
    : [zeta_vl] "r"(avl), [zeta_base0] "r"(&zetas_pos[32]), [zeta_base1] "r"(&zetas_pos[48])
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
    : [zeta_vl] "r"(avl), [zeta_base0] "r"(&zetas_pos[64]), [zeta_base1] "r"(&zetas_pos[80]), [zeta_base2] "r"(&zetas_pos[96]), [zeta_base3] "r"(&zetas_pos[112])
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
    : [zeta_vl] "r"(avl), [zeta_base0] "r"(&zetas_pos[128]), [zeta_base1] "r"(&zetas_pos[144]), [zeta_base2] "r"(&zetas_pos[160]), [zeta_base3] "r"(&zetas_pos[176]),
      [zeta_base4] "r"(&zetas_pos[192]), [zeta_base5] "r"(&zetas_pos[208]), [zeta_base6] "r"(&zetas_pos[224]), [zeta_base7] "r"(&zetas_pos[240])
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
    : [avl]"r"(avl), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&coeff_ptr[0]), [index_base1]"r"(&tree_byteoffset[128])
  );
}

/*************************************************
* Name:        poly_pointwisemul_invntt_mon2nor_asm
*
* Description: 1. Montgomery multiply two polynomials (already in NTT field)
*              2. Compute INTT
*              3. Perform Mon2nor 
*
* Arguments:   - poly *z:            pointer to output polynomial
*              - const poly *a:      pointer to input polynomial a
*              - const poly *b:      pointer to input polynomial b
* Function:    z=Mon2nor(Invntt(a*b))
**************************************************/
void poly_pointwisemul_invntt_mon2nor_asm(poly *z, const poly *a, const poly *b){
    size_t load_zeta_vl;
    size_t avl;//m1
    int32_t* coeff_ptr_a=a->coeffs;
    int32_t* coeff_ptr_b=b->coeffs;
    int32_t* coeff_ptr_z=z->coeffs;

    //stage1, same num=1
    //v4 and v6 with m1 stores input coefficients
    //v16, v18, v20, v22, v24, v26, v28, v30 with m2 stores output coefficients
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 7));//repeat time here should be set large enough
    load_zeta_vl=16;
    __asm__ __volatile__ ( 
    "vsetvli zero, %[zeta_vl], e32, m1, ta, mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vle32.v v4, (%[a_top_base_0])\n"
    "vle32.v v6, (%[a_btm_base_0])\n"
    "vle32.v v5, (%[b_top_base_0])\n"
    "vle32.v v7, (%[b_btm_base_0])\n"
    "vmod.mul.vv v4, v4, v5\n"
    "vmod.mul.vv v6, v6, v7\n"
    "vbutterfly.gs.vvm v16, v4, v6, v0\n"
    :
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[0]), [a_top_base_0]"r"(&coeff_ptr_a[0]), [a_btm_base_0]"r"(&coeff_ptr_a[128]),
      [b_top_base_0]"r"(&coeff_ptr_b[0]), [b_btm_base_0]"r"(&coeff_ptr_b[128])
  );

    __asm__ __volatile__ (
    "vle32.v v0, (%[zeta_base_1])\n"
    "vle32.v v4, (%[a_top_base_1])\n"
    "vle32.v v6, (%[a_btm_base_1])\n"
    "vle32.v v5, (%[b_top_base_1])\n"
    "vle32.v v7, (%[b_btm_base_1])\n"
    "vmod.mul.vv v4, v4, v5\n"
    "vmod.mul.vv v6, v6, v7\n"
    "vbutterfly.gs.vvm v18, v4, v6, v0\n"
    :
    : [zeta_base_1]"r"(&zetas_inv_in_order[16]), [a_top_base_1]"r"(&coeff_ptr_a[16]), [a_btm_base_1]"r"(&coeff_ptr_a[144]),
      [b_top_base_1]"r"(&coeff_ptr_b[16]), [b_btm_base_1]"r"(&coeff_ptr_b[144])
  );

    __asm__ __volatile__ (
    "vle32.v v0, (%[zeta_base_2])\n"
    "vle32.v v4, (%[a_top_base_2])\n"
    "vle32.v v6, (%[a_btm_base_2])\n"
    "vle32.v v5, (%[b_top_base_2])\n"
    "vle32.v v7, (%[b_btm_base_2])\n"
    "vmod.mul.vv v4, v4, v5\n"
    "vmod.mul.vv v6, v6, v7\n"
    "vbutterfly.gs.vvm v20, v4, v6, v0\n"
    :
    : [zeta_base_2]"r"(&zetas_inv_in_order[32]), [a_top_base_2]"r"(&coeff_ptr_a[32]), [a_btm_base_2]"r"(&coeff_ptr_a[160]),
      [b_top_base_2]"r"(&coeff_ptr_b[32]), [b_btm_base_2]"r"(&coeff_ptr_b[160])
  );

    __asm__ __volatile__ (
    "vle32.v v0, (%[zeta_base_3])\n"
    "vle32.v v4, (%[a_top_base_3])\n"
    "vle32.v v6, (%[a_btm_base_3])\n"
    "vle32.v v5, (%[b_top_base_3])\n"
    "vle32.v v7, (%[b_btm_base_3])\n"
    "vmod.mul.vv v4, v4, v5\n"
    "vmod.mul.vv v6, v6, v7\n"
    "vbutterfly.gs.vvm v22, v4, v6, v0\n"
    :
    : [zeta_base_3]"r"(&zetas_inv_in_order[48]), [a_top_base_3]"r"(&coeff_ptr_a[48]), [a_btm_base_3]"r"(&coeff_ptr_a[176]),
      [b_top_base_3]"r"(&coeff_ptr_b[48]), [b_btm_base_3]"r"(&coeff_ptr_b[176])
  );

    __asm__ __volatile__ (
    "vle32.v v0, (%[zeta_base_4])\n"
    "vle32.v v4, (%[a_top_base_4])\n"
    "vle32.v v6, (%[a_btm_base_4])\n"
    "vle32.v v5, (%[b_top_base_4])\n"
    "vle32.v v7, (%[b_btm_base_4])\n"
    "vmod.mul.vv v4, v4, v5\n"
    "vmod.mul.vv v6, v6, v7\n"
    "vbutterfly.gs.vvm v24, v4, v6, v0\n"
    :
    : [zeta_base_4]"r"(&zetas_inv_in_order[64]), [a_top_base_4]"r"(&coeff_ptr_a[64]), [a_btm_base_4]"r"(&coeff_ptr_a[192]),
      [b_top_base_4]"r"(&coeff_ptr_b[64]), [b_btm_base_4]"r"(&coeff_ptr_b[192])
  );

    __asm__ __volatile__ (
    "vle32.v v0, (%[zeta_base_5])\n"
    "vle32.v v4, (%[a_top_base_5])\n"
    "vle32.v v6, (%[a_btm_base_5])\n"
    "vle32.v v5, (%[b_top_base_5])\n"
    "vle32.v v7, (%[b_btm_base_5])\n"
    "vmod.mul.vv v4, v4, v5\n"
    "vmod.mul.vv v6, v6, v7\n"
    "vbutterfly.gs.vvm v26, v4, v6, v0\n"
    :
    : [zeta_base_5]"r"(&zetas_inv_in_order[80]), [a_top_base_5]"r"(&coeff_ptr_a[80]), [a_btm_base_5]"r"(&coeff_ptr_a[208]),
      [b_top_base_5]"r"(&coeff_ptr_b[80]), [b_btm_base_5]"r"(&coeff_ptr_b[208])
  );

    __asm__ __volatile__ (
    "vle32.v v0, (%[zeta_base_6])\n"
    "vle32.v v4, (%[a_top_base_6])\n"
    "vle32.v v6, (%[a_btm_base_6])\n"
    "vle32.v v5, (%[b_top_base_6])\n"
    "vle32.v v7, (%[b_btm_base_6])\n"
    "vmod.mul.vv v4, v4, v5\n"
    "vmod.mul.vv v6, v6, v7\n"
    "vbutterfly.gs.vvm v28, v4, v6, v0\n"
    :
    : [zeta_base_6]"r"(&zetas_inv_in_order[96]), [a_top_base_6]"r"(&coeff_ptr_a[96]), [a_btm_base_6]"r"(&coeff_ptr_a[224]),
      [b_top_base_6]"r"(&coeff_ptr_b[96]), [b_btm_base_6]"r"(&coeff_ptr_b[224])
  );

    __asm__ __volatile__ (
    "vle32.v v0, (%[zeta_base_7])\n"
    "vle32.v v4, (%[a_top_base_7])\n"
    "vle32.v v6, (%[a_btm_base_7])\n"
    "vle32.v v5, (%[b_top_base_7])\n"
    "vle32.v v7, (%[b_btm_base_7])\n"
    "vmod.mul.vv v4, v4, v5\n"
    "vmod.mul.vv v6, v6, v7\n"
    "vbutterfly.gs.vvm v30, v4, v6, v0\n"
    :
    : [zeta_base_7]"r"(&zetas_inv_in_order[112]), [a_top_base_7]"r"(&coeff_ptr_a[112]), [a_btm_base_7]"r"(&coeff_ptr_a[240]),
      [b_top_base_7]"r"(&coeff_ptr_b[112]), [b_btm_base_7]"r"(&coeff_ptr_b[240])
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
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[128]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&zetas_inv_in_order[144]), [zeta_base_2]"r"(&zetas_inv_in_order[160]), [zeta_base_3]"r"(&zetas_inv_in_order[176])
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
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[192]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&zetas_inv_in_order[200]), [zeta_base_2]"r"(&zetas_inv_in_order[208]), [zeta_base_3]"r"(&zetas_inv_in_order[216])
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
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[224]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&zetas_inv_in_order[228]), [zeta_base_2]"r"(&zetas_inv_in_order[232]), [zeta_base_3]"r"(&zetas_inv_in_order[236])
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
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[240]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&zetas_inv_in_order[242]), [zeta_base_2]"r"(&zetas_inv_in_order[244]), [zeta_base_3]"r"(&zetas_inv_in_order[246])
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
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[248]), [avl]"r"(avl), 
      [zeta_base_1]"r"(&zetas_inv_in_order[249]), [zeta_base_2]"r"(&zetas_inv_in_order[250]), [zeta_base_3]"r"(&zetas_inv_in_order[251])
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
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[252]), [avl]"r"(avl), [zeta_base_1]"r"(&zetas_inv_in_order[253])
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
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[254]), [avl]"r"(avl)
  );

    //Mon2nor
    int32_t cons=2365951;
    size_t avl_m8=128;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m8,ta,mu\n"
    "vmod.mul.vx v16,v16,%[cons]\n"
    "vmod.mul.vx v24,v24,%[cons]\n"
    :
    : [avl]"r"(avl_m8), [cons]"r"(cons)
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
    : [avl]"r"(avl), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&coeff_ptr_z[0]), [index_base1]"r"(&tree_byteoffset[128])
  );
}

/*************************************************
* Name:        poly_add_reduce_asm
*
* Description: 1. Add two polynomials
*              2. Perform reduce
*
* Arguments:   - poly *z:            pointer to output polynomial
*              - const poly *a:      pointer to input polynomial a
*              - const poly *b:      pointer to input polynomial b
* Function:    z=Reduce(a+b)
**************************************************/
void poly_add_reduce_asm(poly *z, const poly *a, const poly *b){
  int32_t* coeff_ptr_z=z->coeffs;
  int32_t* coeff_ptr_a=a->coeffs;
  int32_t* coeff_ptr_b=b->coeffs;
  size_t avl_m8=128;
  __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl_m8],e32,m8,ta,mu\n"
    "vle32.v v8, (%[a_addr_0])\n"
    "vle32.v v16, (%[b_addr_0])\n"
    "vmod.add.vv v16, v16, v8\n"
    "vle32.v v8, (%[a_addr_1])\n"
    "vle32.v v24, (%[b_addr_1])\n"
    "vmod.add.vv v24, v24, v8\n"
    :
    : [avl_m8]"r"(avl_m8), [a_addr_0]"r"(&coeff_ptr_a[0]), [b_addr_0]"r"(&coeff_ptr_b[0]), [a_addr_1]"r"(&coeff_ptr_a[128]), [b_addr_1]"r"(&coeff_ptr_b[128])
  );

  __asm__ __volatile__ ( 
    "vadd.vx v8, v16, %[cons]\n"
    "vsra.vx v8, v8, %[offset]\n"
    "vmul.vx v8, v8, %[q]\n"
    "vsub.vv v16, v16, v8\n"
    "vadd.vx v8, v24, %[cons]\n"
    "vsra.vx v8, v8, %[offset]\n"
    "vmul.vx v8, v8, %[q]\n"
    "vsub.vv v24, v24, v8\n"
    "vse32.v v16, (%[z_addr_0])\n"
    "vse32.v v24, (%[z_addr_1])\n"
    :
    : [cons]"r"(1<<22), [offset]"r"(23), [q]"r"(Q), [z_addr_0]"r"(&coeff_ptr_z[0]), [z_addr_1]"r"(&coeff_ptr_z[128])
  );
}

/*************************************************
* Name:        poly_caddq_modsub_reduce_asm
*
* Description: 1. Mod add q on coefficients of polynomial a
*              2. Sub polynomial a with polynomial b
*              3. Reduce polynomial a
*
* Arguments:   - poly *a:            pointer to both input and output polynomial
*              - const poly *b:      pointer to input polynomial b
* Function:    a=Reduce(caddq(a)-b)
**************************************************/
void poly_caddq_modsub_reduce_asm(poly *a, const poly *b){
  int32_t* coeff_ptr_a=a->coeffs;
  int32_t* coeff_ptr_b=b->coeffs;
  size_t avl_m8=128;
  //Handling elements 0-127
  __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl_m8],e32,m8,ta,mu\n"
    "vle32.v v8, (%[a_addr])\n"
    "vle32.v v16, (%[b_addr])\n"
    "vmod.add.vx v8, v8, %[q]\n"
    "vmod.sub.vv v16, v8, v16\n"
    "vadd.vx v8, v16, %[cons]\n"
    "vsra.vx v8, v8, %[offset]\n"
    "vmul.vx v8, v8, %[q]\n"
    "vsub.vv v16, v16, v8\n"
    :
    : [avl_m8]"r"(avl_m8), [a_addr]"r"(&coeff_ptr_a[0]), [b_addr]"r"(&coeff_ptr_b[0]), [q]"r"(Q), [cons]"r"(1<<22), [offset]"r"(23)
  );

  //Handling elements 128-255
  __asm__ __volatile__ ( 
    "vle32.v v8, (%[a_addr])\n"
    "vle32.v v24, (%[b_addr])\n"
    "vmod.add.vx v8, v8, %[q]\n"
    "vmod.sub.vv v24, v8, v24\n"
    "vadd.vx v8, v24, %[cons]\n"
    "vsra.vx v8, v8, %[offset]\n"
    "vmul.vx v8, v8, %[q]\n"
    "vsub.vv v24, v24, v8\n"
    :
    : [a_addr]"r"(&coeff_ptr_a[128]), [b_addr]"r"(&coeff_ptr_b[128]), [q]"r"(Q), [cons]"r"(1<<22), [offset]"r"(23)
  );

  //store back
  __asm__ __volatile__ ( 
    "vse32.v v16, (%[a_addr0])\n"
    "vse32.v v24, (%[a_addr1])\n"
    :
    : [a_addr0]"r"(&coeff_ptr_a[0]), [a_addr1]"r"(&coeff_ptr_a[128]) 
  );
}

/*************************************************
* Name:        poly_add_caddq_asm
*
* Description: 1. Add polynomial a and polynomial b (No mudular reduction)
*              2. Mod add q on the results and store back to polynomial c
*
* Arguments:   - poly *c:            pointer to output polynomial
*              - const poly *a:      pointer to input polynomial a
*              - const poly *b:      pointer to input polynomial b
* Function:    c=caddq(a+b)
**************************************************/
void poly_add_caddq_asm(poly *c, const poly *a, const poly *b){
  int32_t* coeff_ptr_c=c->coeffs;
  int32_t* coeff_ptr_a=a->coeffs;
  int32_t* coeff_ptr_b=b->coeffs;
  size_t avl_m8=128;
  __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl_m8],e32,m8,ta,mu\n"
    "vle32.v v8, (%[a_addr_0])\n"
    "vle32.v v16, (%[b_addr_0])\n"
    "vadd.vv v16, v16, v8\n"
    "vle32.v v8, (%[a_addr_1])\n"
    "vle32.v v24, (%[b_addr_1])\n"
    "vadd.vv v24, v24, v8\n"
    "vmod.add.vx v16, v16, %[q]\n"
    "vmod.add.vx v24, v24, %[q]\n"
    "vse32.v v16, (%[c_addr_0])\n"
    "vse32.v v24, (%[c_addr_1])\n"
    :
    : [avl_m8]"r"(avl_m8), [a_addr_0]"r"(&coeff_ptr_a[0]), [b_addr_0]"r"(&coeff_ptr_b[0]), 
      [a_addr_1]"r"(&coeff_ptr_a[128]), [b_addr_1]"r"(&coeff_ptr_b[128]), [q]"r"(Q), 
      [c_addr_0]"r"(&coeff_ptr_c[0]), [c_addr_1]"r"(&coeff_ptr_c[128])
  );
}

/*************************************************
* Name:        poly_unpackt1_shiftl_ntt_asm
*
* Description: 1. Unpack polynomial t1 with 10-bit coefficients from pk
*              2. Perform shiftl
*              3. Perform ntt and store back
*
* Arguments:   - poly *r:            pointer to output polynomial
*              - uint8_t *packed_addr: pointer to position to store the packed data      
**************************************************/
void poly_unpackt1_shiftl_ntt_asm(poly *r, const uint8_t* packed_addr){
  size_t load_zeta_vl,avl_m1,avl_m2,avl_m4,packed_vl,packed_vlx2,actual_width;
  uint8_t* packed_addr0;
  uint8_t* packed_addr1;
  uint8_t* packed_addr2;
  uint8_t* packed_addr3;
  int32_t* coeff_ptr=r->coeffs;
  packed_addr0=packed_addr;
  avl_m1=16,avl_m2=32,avl_m4=64;

  actual_width=10;
  packed_addr1=packed_addr0+20;
  packed_addr2=packed_addr0+160;
  packed_addr3=packed_addr0+180;
  packed_vl=20;
  packed_vlx2=40;

  /* 
    ** Firstly, Unpack process+shiftl+ first stage of NTT
    ** 1. Load part of packed data into v3 and gradually unpack into v4,v5,v6,v7 with m1
    ** 3. Shiftl on v4 with m4
    ** 4. Compute the first stage of NTT
    ** 5. Finally, v16, v20, v24, v28 with m4 stores the result of NTT after first stage
    */
  //elements 0-31 in v4 and v5, elements 128-159 in v6 and v7
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  load_zeta_vl=1;
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vsll.vx v4, v4, %[cons]\n"
    "vsetvli zero, %[zeta_vl], e32, m1, ta, mu\n"
    "vle32.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [cons]"r"(D), [zeta_vl]"r"(load_zeta_vl), [avl_m2]"r"(avl_m2), [zeta_base]"r"(&zetas_pos[1])
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 32-63 in v4 and v5, elements 160-191 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vsll.vx v4, v4, %[cons]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [avl_m2]"r"(avl_m2), [cons]"r"(D)
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 64-95 in v4 and v5, elements 192-223 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vsll.vx v4, v4, %[cons]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [avl_m2]"r"(avl_m2), [cons]"r"(D)
  );
  packed_addr0+=packed_vlx2,packed_addr1+=packed_vlx2,packed_addr2+=packed_vlx2,packed_addr3+=packed_vlx2;

  //elements 96-127 in v4 and v5, elements 224-255 in v6 and v7
  __asm__ __volatile__ (
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_0])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v4, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_1])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v5, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_2])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v6, v3, %[actual_width]\n"
    "vsetvli	zero, %[packed_vl], e8, m1, ta, mu\n"
    "vle8.v v3, (%[packed_addr_3])\n"
    "vsetvli	zero, %[avl_m1], e32, m1, ta, mu\n"
    "unpack.vx v7, v3, %[actual_width]\n"
    "vsetvli	zero, %[avl_m4], e32, m4, ta, mu\n"
    "vsll.vx v4, v4, %[cons]\n"
    "vsetvli	zero, %[avl_m2], e32, m2, ta, mu\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [packed_vl]"r"(packed_vl), [packed_addr_0]"r"(packed_addr0), [packed_addr_1]"r"(packed_addr1),
      [packed_addr_2]"r"(packed_addr2), [packed_addr_3]"r"(packed_addr3), [avl_m1]"r"(avl_m1), [actual_width]"r"(actual_width), 
      [avl_m4]"r"(avl_m4), [avl_m2]"r"(avl_m2), [cons]"r"(D)
  );
  
  /* 
  ** Secondly, performing the remaining stages of NTT
  */
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[2]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[4]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[8]), [avl]"r"(avl_m2)
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
    : [zeta_vl] "r"(load_zeta_vl), [zeta_base]"r"(&zetas_pos[16]), [avl]"r"(avl_m2)
  );

    //stage6, repeat bound=32
    // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m1 are input coefficients, with same zetas
    // (v25, v5), (v27, v7), (v29, v9), (v31, v11) with m1 are input coefficients, with same zetas
    // v12, v14, v16, v18, v20, v22, v2, v4 with m2 are output coefficients
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[32]), [zeta_base1] "r"(&zetas_pos[48])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[64]), [zeta_base1] "r"(&zetas_pos[80]), [zeta_base2] "r"(&zetas_pos[96]), [zeta_base3] "r"(&zetas_pos[112])
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
    : [zeta_vl] "r"(avl_m1), [zeta_base0] "r"(&zetas_pos[128]), [zeta_base1] "r"(&zetas_pos[144]), [zeta_base2] "r"(&zetas_pos[160]), [zeta_base3] "r"(&zetas_pos[176]),
      [zeta_base4] "r"(&zetas_pos[192]), [zeta_base5] "r"(&zetas_pos[208]), [zeta_base6] "r"(&zetas_pos[224]), [zeta_base7] "r"(&zetas_pos[240])
  );

    //store back, from bit-reverse order to normal-order
    //v8-v15 used to store index, v16-v31 contains coefficients
    uint32_t* tree_ptr=NULL;
    size_t avl=128;
    __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m8,ta,mu\n"
    "vle32.v v8, (%[index_base0])\n"
    "vsuxei32.v v16, (%[coeff_base]), v8\n"     // indexed store v16-v23
    "vle32.v v8, (%[index_base1])\n"
    "vsuxei32.v v24, (%[coeff_base]), v8\n"     // indexed store v24-v31
    :
    : [avl]"r"(avl), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&coeff_ptr[0]), [index_base1]"r"(&tree_byteoffset[128])
  );
}

/*************************************************
* Name:        poly_macc_sub_invntt_mon2nor_asm
*
* Description: 1. Pointwise multiply two polyvecl and accumulate into a polynomial
*              2. Sub the polynomial
*              3. Perform INVNTT
*              4. Perform Mon2nor
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const polyvecl *a:   pointer to input polyvecl 
*              - const polyvecl *b:   pointer to input polyvecl
*              - const poly *c:       pointer to input polynomial
* Function:    r=mon2nor(invntt((a*b)-c))
**************************************************/
void poly_macc_sub_invntt_mon2nor_asm(poly *r, const polyvecl *a, const polyvecl *b, const poly *c){
  size_t load_zeta_vl,avl_m1,avl_m2,avl_m4,avl_m8;
  avl_m1=16,avl_m2=32,avl_m4=64,avl_m8=128;
  const int32_t* coeff_ptr_a=a->vec[0].coeffs;
  const int32_t* coeff_ptr_b=b->vec[0].coeffs;
  const int32_t* coeff_ptr_c=c->coeffs;

  csr_modulusq_rw(Q);
  csr_qinv_rw(QINV_HW);

  /*
  ** Firstly, compute polynomial multiplication of a->vec[0] and b->vec[0]
  ** Result stored in v16-v31, which form the accumulation base
  */
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl_m8],e32,m8,ta,mu\n"
    "vle32.v v8, (%[a_ptr0])\n"
    "vle32.v v16, (%[b_ptr0])\n"
    "vmod.mul.vv v16, v16, v8\n"
    "vle32.v v8, (%[a_ptr1])\n"
    "vle32.v v24, (%[b_ptr1])\n"
    "vmod.mul.vv v24, v24, v8\n"
    :
    : [avl_m8]"r"(avl_m8), [a_ptr0]"r"(&coeff_ptr_a[0]), [b_ptr0]"r"(&coeff_ptr_b[0]), 
      [a_ptr1]"r"(&coeff_ptr_a[128]), [b_ptr1]"r"(&coeff_ptr_b[128])
  );

  /*
  ** Secondly, compute polynomial multiplication of a->vec[i] and b->vec[i], i = 1, 2, ..., L-1
  ** 1. Load partial of a->vec[i] into v4 and partial of b->vec[i] into v8 with m4 set
  ** 2. Mod mul between v4 and v8 and results stored in v8 with m4
  ** 3. Mod add v8 with corresponding part in v16-v31
  ** 4. Finally, v16-v31 stores the results of macc
  */
  for(int i=1;i<L;i++){
    coeff_ptr_a=a->vec[i].coeffs;
    coeff_ptr_b=b->vec[i].coeffs;
    __asm__ __volatile__ (
    "vsetvli	zero,%[avl_m4],e32,m4,ta,mu\n"
    "vle32.v v4, (%[a_ptr0])\n"
    "vle32.v v8, (%[b_ptr0])\n"
    "vmod.mul.vv v8, v8, v4\n"
    "vmod.add.vv v16,v16,v8\n"
    "vle32.v v4, (%[a_ptr1])\n"
    "vle32.v v8, (%[b_ptr1])\n"
    "vmod.mul.vv v8, v8, v4\n"
    "vmod.add.vv v20,v20,v8\n"
    "vle32.v v4, (%[a_ptr2])\n"
    "vle32.v v8, (%[b_ptr2])\n"
    "vmod.mul.vv v8, v8, v4\n"
    "vmod.add.vv v24,v24,v8\n"
    "vle32.v v4, (%[a_ptr3])\n"
    "vle32.v v8, (%[b_ptr3])\n"
    "vmod.mul.vv v8, v8, v4\n"
    "vmod.add.vv v28,v28,v8\n"
    :
    : [avl_m4]"r"(avl_m4), [a_ptr0]"r"(&coeff_ptr_a[0]), [b_ptr0]"r"(&coeff_ptr_b[0]), 
      [a_ptr1]"r"(&coeff_ptr_a[64]), [b_ptr1]"r"(&coeff_ptr_b[64]), 
      [a_ptr2]"r"(&coeff_ptr_a[128]), [b_ptr2]"r"(&coeff_ptr_b[128]), 
      [a_ptr3]"r"(&coeff_ptr_a[192]), [b_ptr3]"r"(&coeff_ptr_b[192])
  );
  }

  /*
  ** Thirdly, subtract coefficients in v16-v31 with polynomial t1
  */
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl_m8],e32,m8,ta,mu\n"
    "vle32.v v8, (%[c_addr0])\n"
    "vmod.sub.vv v16, v16, v8\n"
    "vle32.v v8, (%[c_addr1])\n"
    "vmod.sub.vv v24, v24, v8\n"
    :
    : [avl_m8]"r"(avl_m8), [c_addr0]"r"(&coeff_ptr_c[0]), [c_addr1]"r"(&coeff_ptr_c[128])
  );

  /*
  ** Fourthly, compute invntt on coefficients in v16-v31
  */

  //stage1, same num=1
  //v16-v31 with m1 stores input coefficients
  //v4, v6, v8, v10, v12, v14, v16, v18 with m2 stores output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 7));//repeat time here should be set large enough
  load_zeta_vl=16;
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl_m1],e32,m1,ta,mu\n"
    "vle32.v v0, (%[zeta_base_0])\n"
    "vbutterfly.gs.vvm v4, v16, v24, v0\n"
    "vle32.v v0, (%[zeta_base_1])\n"
    "vbutterfly.gs.vvm v6, v17, v25, v0\n"
    "vle32.v v0, (%[zeta_base_2])\n"
    "vbutterfly.gs.vvm v8, v18, v26, v0\n"
    "vle32.v v0, (%[zeta_base_3])\n"
    "vbutterfly.gs.vvm v10, v19, v27, v0\n"
    "vle32.v v0, (%[zeta_base_4])\n"
    "vbutterfly.gs.vvm v12, v20, v28, v0\n"
    "vle32.v v0, (%[zeta_base_5])\n"
    "vbutterfly.gs.vvm v14, v21, v29, v0\n"
    "vle32.v v0, (%[zeta_base_6])\n"
    "vbutterfly.gs.vvm v16, v22, v30, v0\n"
    "vle32.v v0, (%[zeta_base_7])\n"
    "vbutterfly.gs.vvm v18, v23, v31, v0\n"
    :
    : [avl_m1]"r"(avl_m1), [zeta_base_0]"r"(&zetas_inv_in_order[0]), [zeta_base_1]"r"(&zetas_inv_in_order[16]),
      [zeta_base_2]"r"(&zetas_inv_in_order[32]), [zeta_base_3]"r"(&zetas_inv_in_order[48]), 
      [zeta_base_4]"r"(&zetas_inv_in_order[64]), [zeta_base_5]"r"(&zetas_inv_in_order[80]),
      [zeta_base_6]"r"(&zetas_inv_in_order[96]), [zeta_base_7]"r"(&zetas_inv_in_order[112])
  );

  //stage2, same num=2
  // (v4, v12), (v6, v14), (v8, v16), (v10, v18) with m2 are input coefficients
  // v20, v24, v28, v4 with m4 is output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(1, 7));
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
    : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[128]), [avl]"r"(avl_m2), 
      [zeta_base_1]"r"(&zetas_inv_in_order[144]), [zeta_base_2]"r"(&zetas_inv_in_order[160]), [zeta_base_3]"r"(&zetas_inv_in_order[176])
  );

  //stage3, same num=4
  // (v20, v28), (v22, v30), (v24, v4), (v26, v6) with m2 are input coefficients
  // v8, v12, v16, v20 is output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(2, 7));
  load_zeta_vl=8;
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
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[192]), [avl]"r"(avl_m2), 
    [zeta_base_1]"r"(&zetas_inv_in_order[200]), [zeta_base_2]"r"(&zetas_inv_in_order[208]), [zeta_base_3]"r"(&zetas_inv_in_order[216])
);

  //stage4, same num=8
  // (v8, v16), (v10, v18), (v12, v20), (v14, v22) with m2 are input coefficients
  // v24, v28, v4, v8 is output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(3, 7));
  load_zeta_vl=4;
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
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[224]), [avl]"r"(avl_m2), 
    [zeta_base_1]"r"(&zetas_inv_in_order[228]), [zeta_base_2]"r"(&zetas_inv_in_order[232]), [zeta_base_3]"r"(&zetas_inv_in_order[236])
);

  //stage5, same num=16
  // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m2 are input coefficients
  // v12, v16, v20, v24 is output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(4, 7));
  load_zeta_vl=2;
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
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[240]), [avl]"r"(avl_m2), 
    [zeta_base_1]"r"(&zetas_inv_in_order[242]), [zeta_base_2]"r"(&zetas_inv_in_order[244]), [zeta_base_3]"r"(&zetas_inv_in_order[246])
);

  //stage6, same num=32
  // (v12, v20), (v14, v22), (v16, v24), (v18, v26) with m2 are input coefficients
  //  v28, v4, v8, v12 with m4 are output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(5, 7));
  load_zeta_vl=1;
  __asm__ __volatile__ ( 
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_0])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v28, v12, v20, v0\n"
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_1])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v4, v14, v22, v0\n"
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_2])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v8, v16, v24, v0\n"
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_3])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v12, v18, v26, v0\n"
  :
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[248]), [avl]"r"(avl_m2), 
    [zeta_base_1]"r"(&zetas_inv_in_order[249]), [zeta_base_2]"r"(&zetas_inv_in_order[250]), [zeta_base_3]"r"(&zetas_inv_in_order[251])
);

  //stage7, same num=64
  // (v28, v8), (v30, v10), (v4, v12), (v6, v14) with m2 are input coefficients
  // v16, v20, v24, v28 with m4 are output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(6, 7));
  __asm__ __volatile__ ( 
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_0])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v16, v28, v8, v0\n"
  "vbutterfly.gs.vvm v20, v30, v10, v0\n"
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_1])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v24, v4, v12, v0\n"
  "vbutterfly.gs.vvm v28, v6, v14, v0\n"
  :
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[252]), [avl]"r"(avl_m2), [zeta_base_1]"r"(&zetas_inv_in_order[253])
);

  //stage8, same num=128
  // (v16, v24), (v18, v26), (v20, v28), (v22, v30) with m2 are input coefficients
  // v4, v8, v12, v16 with m4 are output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(7, 7));
  __asm__ __volatile__ ( 
  "vsetvli	zero,%[zeta_vl],e32,m1,ta,mu\n"
  "vle32.v v0, (%[zeta_base_0])\n"
  "vsetvli	zero,%[avl],e32,m2,ta,mu\n"
  "vbutterfly.gs.vvm v4, v16, v24, v0\n"
  "vbutterfly.gs.vvm v8, v18, v26, v0\n"
  "vbutterfly.gs.vvm v12, v20, v28, v0\n"
  "vbutterfly.gs.vvm v16, v22, v30, v0\n"
  :
  : [zeta_vl]"r"(load_zeta_vl), [zeta_base_0]"r"(&zetas_inv_in_order[254]), [avl]"r"(avl_m2)
);

  /*
  ** Fifthly, compute mon2nor on coefficients in v4-v19 with m4
  */
  int32_t cons=2365951;
  __asm__ __volatile__ ( 
  "vsetvli	zero,%[avl],e32,m4,ta,mu\n"
  "vmod.mul.vx v4,v4,%[cons]\n"
  "vmod.mul.vx v8,v8,%[cons]\n"
  "vmod.mul.vx v12,v12,%[cons]\n"
  "vmod.mul.vx v16,v16,%[cons]\n"
  :
  : [avl]"r"(avl_m4), [cons]"r"(cons)
);

  //store back, back, from bit-reverse order to normal-order
  __asm__ __volatile__ ( 
    "vsetvli	zero,%[avl],e32,m4,ta,mu\n"
    "vle32.v v28, (%[index_base0])\n"
    "vsuxei32.v v4, (%[coeff_base]), v28\n"     
    "vle32.v v28, (%[index_base1])\n"
    "vsuxei32.v v8, (%[coeff_base]), v28\n" 
    "vle32.v v28, (%[index_base2])\n"
    "vsuxei32.v v12, (%[coeff_base]), v28\n"
    "vle32.v v28, (%[index_base3])\n"
    "vsuxei32.v v16, (%[coeff_base]), v28\n"    
    :
    : [avl]"r"(avl_m4), [index_base0]"r"(&tree_byteoffset[0]), [coeff_base]"r"(&(r->coeffs[0])), [index_base1]"r"(&tree_byteoffset[64]),
      [index_base2]"r"(&tree_byteoffset[128]), [index_base3]"r"(&tree_byteoffset[192])
  );
}
#elif (VLEN == 1024)

# else
#error "VLEN must be 256/512/1024"
#endif