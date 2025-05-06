#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"

#if (VLEN == 256)
const int16_t zetas_cg[11][16] = 
{
  {2571, 2571, 2571, 2571, 2571, 2571, 2571, 2571, 2571, 2571, 2571, 2571, 2571, 2571, 2571, 2571},
  {2970, 1812, 2970, 1812, 2970, 1812, 2970, 1812, 2970, 1812, 2970, 1812, 2970, 1812, 2970, 1812},
  {1493, 1422, 287, 202, 1493, 1422, 287, 202, 1493, 1422, 287, 202, 1493, 1422, 287, 202},
  {3158, 622, 1577, 182, 962, 2127, 1855, 1468, 3158, 622, 1577, 182, 962, 2127, 1855, 1468},
  {573, 2004, 264, 383, 2500, 1458, 1727, 3199, 2648, 1017, 732, 608, 1787, 411, 3124, 1758},
  {1223, 652, 2777, 1015, 2036, 1491, 3047, 1785, 516, 3321, 3009, 2663, 1711, 2167, 126, 1469},
  {2476, 3239, 3058, 830, 107, 1908, 3082, 2378, 2931, 961, 1821, 2604, 448, 2264, 677, 2054},
  {2226, 430, 555, 843, 2078, 871, 1550, 105, 422, 587, 177, 3094, 3038, 2869, 1574, 1653},
  {3083, 778, 1159, 3182, 2552, 1483, 2727, 1119, 1739, 644, 2457, 349, 418, 329, 3173, 3254},
  {817, 1097, 603, 610, 1322, 2044, 1864, 384, 2114, 3193, 1218, 1994, 2455, 220, 2142, 1670},
  {2144, 1799, 2051, 794, 1819, 2475, 2459, 478, 3221, 3021, 996, 991, 958, 1869, 1522, 1628}
};
#elif (VLEN == 512)

#elif (VLEN == 1024)

# else
#error "VLEN must be 256/512/1024"
#endif

#if (VLEN == 256)
/*************************************************
* Name:        ntt_cg_custom_asm
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in bitreversed order
*              assembly optimized
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_asm(int16_t r[256]) {
  // stage 1
  // a2 is avl, a3 is coefficient start address, a4 is coefficient middle address
  // v4 and v6 with m2 is input coefficients
  // v16, v20, v24, v28 with m4 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "mv a3, %[coeff_base0]\n"
    "mv a4, %[coeff_base1]\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vle16.v v4, (a3)\n"
    "vle16.v v6, (a4)\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    "slli a2, a2, 1\n"
    "add a3, a3, a2\n"
    "add a4, a4, a2\n"
    "vle16.v v4, (a3)\n"
    "vle16.v v6, (a4)\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    "add a3, a3, a2\n"
    "add a4, a4, a2\n"
    "vle16.v v4, (a3)\n"
    "vle16.v v6, (a4)\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    "add a3, a3, a2\n"
    "add a4, a4, a2\n"
    "vle16.v v4, (a3)\n"
    "vle16.v v6, (a4)\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [coeff_base0] "r"(r), [coeff_base1] "r"(&r[128]), [zeta_base] "r"(zetas_cg[0])
  );

  // stage 2
  // a2 is avl
  // (v16, v24), (v18, v26), (v20, v28), (v22, v30) with m2 are input coefficients
  // v4, v8, v12, v16 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v4, v16, v24, v0\n"
    "vbutterfly.ct.vvm v8, v18, v26, v0\n"
    "vbutterfly.ct.vvm v12, v20, v28, v0\n"
    "vbutterfly.ct.vvm v16, v22, v30, v0\n"
    :
    : [zeta_base] "r"(zetas_cg[1])
  );

  // stage 3
  // a2 is avl
  // (v4, v12), (v6, v14), (v8, v16), (v10, v18) with m2 are input coefficients
  // v20, v24, v28, v4 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v20, v4, v12, v0\n"
    "vbutterfly.ct.vvm v24, v6, v14, v0\n"
    "vbutterfly.ct.vvm v28, v8, v16, v0\n"
    "vbutterfly.ct.vvm v4, v10, v18, v0\n"
    :
    : [zeta_base] "r"(zetas_cg[2])
  );

  // stage 4
  // a2 is avl
  // (v20, v28), (v22, v30), (v24, v4), (v26, v6) with m2 are input coefficients
  // v8, v12, v16, v20 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v8, v20, v28, v0\n"
    "vbutterfly.ct.vvm v12, v22, v30, v0\n"
    "vbutterfly.ct.vvm v16, v24, v4, v0\n"
    "vbutterfly.ct.vvm v20, v26, v6, v0\n"
    :
    : [zeta_base] "r"(zetas_cg[3])
  );  
  
  // stage 5
  // a2 is avl
  // (v8, v16), (v10, v18), (v12, v20), (v14, v22) with m2 are input coefficients
  // v24, v28, v4, v8 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v24, v8, v16, v0\n"
    "vbutterfly.ct.vvm v28, v10, v18, v0\n"
    "vbutterfly.ct.vvm v4, v12, v20, v0\n"
    "vbutterfly.ct.vvm v8, v14, v22, v0\n"
    :
    : [zeta_base] "r"(zetas_cg[4])
  );

  // stage 6
  // a2 is avl = 16, a3 is avl = 32
  // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m1 are input coefficients, with same zetas
  // (v25, v5), (v27, v7), (v29, v9), (v31, v11) with m1 are input coefficients, with same zetas
  // v12, v14, v16, v18, v20, v22, v24 v26 with m2 are output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "li a3, 32\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v12, v24, v4, v0\n"
    "vbutterfly.ct.vvm v16, v26, v6, v0\n"
    "vbutterfly.ct.vvm v20, v28, v8, v0\n"
    "vbutterfly.ct.vvm v24, v30, v10, v0\n"
    "vle16.v v0, (%[zeta_base1])\n"
    "vbutterfly.ct.vvm v14, v25, v5, v0\n"
    "vbutterfly.ct.vvm v18, v27, v7, v0\n"
    "vbutterfly.ct.vvm v22, v29, v9, v0\n"
    "vbutterfly.ct.vvm v26, v31, v11, v0\n"
    :
    : [zeta_base0] "r"(zetas_cg[5]), [zeta_base1] "r"(zetas_cg[6])
  );

  // stage 7
  // a2 is avl = 16, a3 is avl = 32
  // (v12, v20), (v13, v21), (v14, v22), (v15, v23) with m2 are input coefficients
  // (v16, v24), (v17, v25), (v18, v26), (v19, v27) with m2 are input coefficients
  // ((v12, v20), (v16, v24)) same zetas, ((v13, v21), (v17, v25)) same zetas
  // ((v14, v22), (v18, v26)) same zetas, ((v15, v23), (v19, v27)) same zetas
  // v4, v6, v8, v10, v12, v14, v28, v30 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "li a3, 32\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v4, v12, v20, v0\n"
    "vbutterfly.ct.vvm v12, v16, v24, v0\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base1])\n"
    "vsetvli	zero,a3,e16,m2,ta,mu\n"    
    "vbutterfly.ct.vvm v6, v13, v21, v0\n"
    "vbutterfly.ct.vvm v14, v17, v25, v0\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base2])\n"
    "vbutterfly.ct.vvm v18, v14, v22, v0\n"
    "vbutterfly.ct.vvm v28, v18, v26, v0\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base3])\n"
    "vbutterfly.ct.vvm v10, v15, v23, v0\n"
    "vbutterfly.ct.vvm v30, v19, v27, v0\n"
    :
    : [zeta_base0] "r"(zetas_cg[7]), [zeta_base1] "r"(zetas_cg[8]), [zeta_base2] "r"(zetas_cg[9]), [zeta_base3] "r"(zetas_cg[10])
  );

  // store back
  // a2 is avl = 64 with m4
  // a3 is coefficient store address
  // v4, v8, v12, v28 to be stored back
  __asm__ __volatile__ ( 
    "li a2, 64\n"
    "mv a3, %[coeff_base]\n"
    "vsetvli	zero,a2,e16,m4,ta,mu\n"
    "vse16.v v4, (a3)\n"
    "slli a2, a2, 1\n"
    "add a3, a3, a2\n"
    "vse16.v v8, (a3)\n"
    "add a3, a3, a2\n"
    "vse16.v v12, (a3)\n"
    "add a3, a3, a2\n"
    "vse16.v v28, (a3)\n"             
    :
    : [coeff_base] "r"(r)
  );  
}

/*************************************************
* Name:        ntt_cg_custom_reorder_asm
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*              assembly optimized
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_reorder_asm(int16_t r[256]) {
  // stage 1
  // a2 is avl, a3 is coefficient start address, a4 is coefficient middle address
  // v4 and v6 with m2 is input coefficients
  // v16, v20, v24, v28 with m4 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "mv a3, %[coeff_base0]\n"
    "mv a4, %[coeff_base1]\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vle16.v v4, (a3)\n"
    "vle16.v v6, (a4)\n"
    "vbutterfly.ct.vvm v16, v4, v6, v0\n"
    "slli a2, a2, 1\n"
    "add a3, a3, a2\n"
    "add a4, a4, a2\n"
    "vle16.v v4, (a3)\n"
    "vle16.v v6, (a4)\n"
    "vbutterfly.ct.vvm v20, v4, v6, v0\n"
    "add a3, a3, a2\n"
    "add a4, a4, a2\n"
    "vle16.v v4, (a3)\n"
    "vle16.v v6, (a4)\n"
    "vbutterfly.ct.vvm v24, v4, v6, v0\n"
    "add a3, a3, a2\n"
    "add a4, a4, a2\n"
    "vle16.v v4, (a3)\n"
    "vle16.v v6, (a4)\n"
    "vbutterfly.ct.vvm v28, v4, v6, v0\n"
    :
    : [coeff_base0] "r"(r), [coeff_base1] "r"(&r[128]), [zeta_base] "r"(zetas_cg[0])
  );

  // stage 2
  // a2 is avl
  // (v16, v24), (v18, v26), (v20, v28), (v22, v30) with m2 are input coefficients
  // v4, v8, v12, v16 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v4, v16, v24, v0\n"
    "vbutterfly.ct.vvm v8, v18, v26, v0\n"
    "vbutterfly.ct.vvm v12, v20, v28, v0\n"
    "vbutterfly.ct.vvm v16, v22, v30, v0\n"
    :
    : [zeta_base] "r"(zetas_cg[1])
  );

  // stage 3
  // a2 is avl
  // (v4, v12), (v6, v14), (v8, v16), (v10, v18) with m2 are input coefficients
  // v20, v24, v28, v4 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v20, v4, v12, v0\n"
    "vbutterfly.ct.vvm v24, v6, v14, v0\n"
    "vbutterfly.ct.vvm v28, v8, v16, v0\n"
    "vbutterfly.ct.vvm v4, v10, v18, v0\n"
    :
    : [zeta_base] "r"(zetas_cg[2])
  );

  // stage 4
  // a2 is avl
  // (v20, v28), (v22, v30), (v24, v4), (v26, v6) with m2 are input coefficients
  // v8, v12, v16, v20 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v8, v20, v28, v0\n"
    "vbutterfly.ct.vvm v12, v22, v30, v0\n"
    "vbutterfly.ct.vvm v16, v24, v4, v0\n"
    "vbutterfly.ct.vvm v20, v26, v6, v0\n"
    :
    : [zeta_base] "r"(zetas_cg[3])
  );  
  
  // stage 5
  // a2 is avl
  // (v8, v16), (v10, v18), (v12, v20), (v14, v22) with m2 are input coefficients
  // v24, v28, v4, v8 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "addi a2, a2, 16\n"
    "vsetvli	zero,a2,e16,m2,ta,mu\n"
    "vbutterfly.ct.vvm v24, v8, v16, v0\n"
    "vbutterfly.ct.vvm v28, v10, v18, v0\n"
    "vbutterfly.ct.vvm v4, v12, v20, v0\n"
    "vbutterfly.ct.vvm v8, v14, v22, v0\n"
    :
    : [zeta_base] "r"(zetas_cg[4])
  );

  // stage 6
  // a2 is avl = 16, a3 is avl = 32
  // (v24, v4), (v26, v6), (v28, v8), (v30, v10) with m1 are input coefficients, with same zetas
  // (v25, v5), (v27, v7), (v29, v9), (v31, v11) with m1 are input coefficients, with same zetas
  // v12, v14, v16, v18, v20, v22, v24 v26 with m2 are output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "li a3, 32\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v12, v24, v4, v0\n"
    "vbutterfly.ct.vvm v16, v26, v6, v0\n"
    "vbutterfly.ct.vvm v20, v28, v8, v0\n"
    "vbutterfly.ct.vvm v24, v30, v10, v0\n"
    "vle16.v v0, (%[zeta_base1])\n"
    "vbutterfly.ct.vvm v14, v25, v5, v0\n"
    "vbutterfly.ct.vvm v18, v27, v7, v0\n"
    "vbutterfly.ct.vvm v22, v29, v9, v0\n"
    "vbutterfly.ct.vvm v26, v31, v11, v0\n"
    :
    : [zeta_base0] "r"(zetas_cg[5]), [zeta_base1] "r"(zetas_cg[6])
  );

  // stage 7
  // a2 is avl = 16, a3 is avl = 32
  // (v12, v20), (v13, v21), (v14, v22), (v15, v23) with m2 are input coefficients
  // (v16, v24), (v17, v25), (v18, v26), (v19, v27) with m2 are input coefficients
  // ((v12, v20), (v16, v24)) same zetas, ((v13, v21), (v17, v25)) same zetas
  // ((v14, v22), (v18, v26)) same zetas, ((v15, v23), (v19, v27)) same zetas
  // v4, v6, v8, v10, v12, v14, v28, v30 is output coefficients
  __asm__ __volatile__ ( 
    "li a2, 16\n"
    "li a3, 32\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base0])\n"
    "vbutterfly.ct.vvm v4, v12, v20, v0\n"
    "vbutterfly.ct.vvm v12, v16, v24, v0\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base1])\n"
    "vsetvli	zero,a3,e16,m2,ta,mu\n"    
    "vbutterfly.ct.vvm v6, v13, v21, v0\n"
    "vbutterfly.ct.vvm v14, v17, v25, v0\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base2])\n"
    "vbutterfly.ct.vvm v18, v14, v22, v0\n"
    "vbutterfly.ct.vvm v28, v18, v26, v0\n"
    "vsetvli	zero,a2,e16,m1,ta,mu\n"
    "vle16.v v0, (%[zeta_base3])\n"
    "vbutterfly.ct.vvm v10, v15, v23, v0\n"
    "vbutterfly.ct.vvm v30, v19, v27, v0\n"
    :
    : [zeta_base0] "r"(zetas_cg[7]), [zeta_base1] "r"(zetas_cg[8]), [zeta_base2] "r"(zetas_cg[9]), [zeta_base3] "r"(zetas_cg[10])
  );

  // store back
  // v4, v8, v12, v28 to be stored back
  // indexed store coefficients from vector register v16~v23 to memory
  // const uint16_t* tree_prt = NULL;

  // tree_prt = byteoffset_even;
  // __asm__ __volatile__ ( 
  //   "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
  //   "vle16.v v16, (%[offset_base])\n"
  //   "vsuxei16.v v4, (%[coeff_base]), v16\n"     // indexed store v4-7
  //   "vadd.vi v16, v16, %[offset_plus]\n"
  //   "vsuxei16.v v12, (%[coeff_base]), v16\n"     // indexed store v12-15
  //   :
  //   : [coeff_base] "r"(r), [avl] "r"(64), [offset_plus] "i"(2), [offset_base] "r"(tree_prt)
  // );

  // tree_prt = byteoffset_even + 64;
  // __asm__ __volatile__ ( 
  //   "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
  //   "vle16.v v16, (%[offset_base])\n"
  //   "vsuxei16.v v8, (%[coeff_base]), v16\n"     // indexed store v8-11
  //   "vadd.vi v16, v16, %[offset_plus]\n"
  //   "vsuxei16.v v28, (%[coeff_base]), v16\n"     // indexed store v28-31
  //   :
  //   : [coeff_base] "r"(r), [avl] "r"(64), [offset_plus] "i"(2), [offset_base] "r"(tree_prt)
  // );
  __asm__ __volatile__ ( 
    "li a2, 64\n"
    "mv a3, %[coeff_base]\n"
    "vsetvli	zero,a2,e16,m4,ta,mu\n"
    "vse16.v v4, (a3)\n"
    "slli a2, a2, 1\n"
    "add a3, a3, a2\n"
    "vse16.v v8, (a3)\n"
    "add a3, a3, a2\n"
    "vse16.v v12, (a3)\n"
    "add a3, a3, a2\n"
    "vse16.v v28, (a3)\n"             
    :
    : [coeff_base] "r"(r)
  );  
}

/*************************************************
* Name:        ntt_cg_custom
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in bitreversed order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom(int16_t r[256]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i, j;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;

  // each bfu uses a single zeta
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));

  vint16m1_t v_zetas;
  // for stage 1 to stage 5
  vint16m2_t v_in0, v_in1;
  vint16m4_t v_out0;  
  // for stage 6 to stage 7
  vint16m1_t v_in2, v_in3;
  vint16m2_t v_out1;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  v_zetas = vle16_v_i16m1(zetas_cg[0], 16);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  v_zetas = vle16_v_i16m1(zetas_cg[1], 16);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }  

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  v_zetas = vle16_v_i16m1(zetas_cg[2], 16);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  v_zetas = vle16_v_i16m1(zetas_cg[3], 16);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  v_zetas = vle16_v_i16m1(zetas_cg[4], 16);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }


  // Vector LMUL Truncation
  v_in2 = vget_v_i16m2_i16m1(v_in0,0);
  v_in3 = vget_v_i16m2_i16m1(v_in1,0);
  v_out1 = vget_v_i16m4_i16m2(v_out0,0);

  // stage 6
  for(i = 0; i < 2; i++) {
    r_ptr0 = r_temp0 + 16*i;
    r_ptr1 = &r_temp0[128 + 16*i];
    r_ptr2 = r_temp1 + 32*i;

    v_zetas = vle16_v_i16m1(zetas_cg[5+i], 16);
    for(j = 0; j < 4; j++) {
      v_in2 = vle16_v_i16m1(r_ptr0, 16);
      v_in3 = vle16_v_i16m1(r_ptr1, 16);
      r_ptr0 += 32;
      r_ptr1 += 32;
      v_out1 = vbutterfly_ct_vvm_i16m2(v_in2, v_in3, v_zetas);
      vse16_v_i16m2(r_ptr2, v_out1, 32);
      r_ptr2 += 64;
    }
  }

  // stage 7
  for(i = 0; i < 4; i++) {
    r_ptr0 = r_temp1 + 16*i;
    r_ptr1 = &r_temp1[128 + 16*i];
    r_ptr2 = r + 32*i;

    v_zetas = vle16_v_i16m1(zetas_cg[7+i], 16);
    for(j = 0; j < 2; j++) {
      v_in2 = vle16_v_i16m1(r_ptr0, 16);
      v_in3 = vle16_v_i16m1(r_ptr1, 16);
      r_ptr0 += 64;
      r_ptr1 += 64;
      v_out1 = vbutterfly_ct_vvm_i16m2(v_in2, v_in3, v_zetas);
      vse16_v_i16m2(r_ptr2, v_out1, 32);
      r_ptr2 += 128;
    }
  }


}


/*************************************************
* Name:        ntt_cg_custom_zeta_compact
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in bitreversed order
*              zeta in compact form, that zetas in v0 are not the same
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_zeta_compact(int16_t r[256]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i, j;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;

  vint16m1_t v_zetas;
  // for stage 1 to stage 5
  vint16m2_t v_in0, v_in1;
  vint16m4_t v_out0;
  // for stage 6 to stage 7
  vint16m1_t v_in2, v_in3;
  vint16m2_t v_out1;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 1
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  v_zetas = vle16_v_i16m1(&zetas[1], 1);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));  
  v_zetas = vle16_v_i16m1(&zetas[2], 2);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }  

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
  v_zetas = vle16_v_i16m1(&zetas[4], 4);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
  v_zetas = vle16_v_i16m1(&zetas[8], 8);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
  v_zetas = vle16_v_i16m1(&zetas[16], 16);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }


  // Vector LMUL Truncation
  v_in2 = vget_v_i16m2_i16m1(v_in0,0);
  v_in3 = vget_v_i16m2_i16m1(v_in1,0);
  v_out1 = vget_v_i16m4_i16m2(v_out0,0);

  // stage 6
  for(i = 0; i < 2; i++) {
    r_ptr0 = r_temp0 + 16*i;
    r_ptr1 = &r_temp0[128 + 16*i];
    r_ptr2 = r_temp1 + 32*i;

    // repeat bound = 16
    v_zetas = vle16_v_i16m1(&zetas[32+16*i], 16);
    for(j = 0; j < 4; j++) {
      v_in2 = vle16_v_i16m1(r_ptr0, 16);
      v_in3 = vle16_v_i16m1(r_ptr1, 16);
      r_ptr0 += 32;
      r_ptr1 += 32;
      v_out1 = vbutterfly_ct_vvm_i16m2(v_in2, v_in3, v_zetas);
      vse16_v_i16m2(r_ptr2, v_out1, 32);
      r_ptr2 += 64;
    }
  }

  // stage 7
  for(i = 0; i < 4; i++) {
    r_ptr0 = r_temp1 + 16*i;
    r_ptr1 = &r_temp1[128 + 16*i];
    r_ptr2 = r + 32*i;

    // repeat bound = 16
    v_zetas = vle16_v_i16m1(&zetas[64+16*i], 16);
    for(j = 0; j < 2; j++) {
      v_in2 = vle16_v_i16m1(r_ptr0, 16);
      v_in3 = vle16_v_i16m1(r_ptr1, 16);
      r_ptr0 += 64;
      r_ptr1 += 64;
      v_out1 = vbutterfly_ct_vvm_i16m2(v_in2, v_in3, v_zetas);
      vse16_v_i16m2(r_ptr2, v_out1, 32);
      r_ptr2 += 128;
    }
  }


}


/*************************************************
* Name:        ntt_cg_custom_reorder
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_reorder(int16_t r[256]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i, j;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;

  vint16m1_t v_zetas;
  // for stage 1 to stage 5
  vint16m2_t v_in0, v_in1;
  vint16m4_t v_out0;
  // for stage 6 to stage 7
  vint16m1_t v_in2, v_in3;
  vint16m2_t v_out1;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 1
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));  
  v_zetas = vle16_v_i16m1(&zetas[1], 1);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));  
  v_zetas = vle16_v_i16m1(&zetas[2], 2);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }  

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
  v_zetas = vle16_v_i16m1(&zetas[4], 4);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
  v_zetas = vle16_v_i16m1(&zetas[8], 8);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
  v_zetas = vle16_v_i16m1(&zetas[16], 16);
  for(i = 0; i < 4; i++) {  // 4 * vectorized butterfly (each 32 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 32);
    v_in1 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }


  // Vector LMUL Truncation
  v_in2 = vget_v_i16m2_i16m1(v_in0,0);
  v_in3 = vget_v_i16m2_i16m1(v_in1,0);
  v_out1 = vget_v_i16m4_i16m2(v_out0,0);

  // stage 6
  for(i = 0; i < 2; i++) {
    r_ptr0 = r_temp0 + 16*i;
    r_ptr1 = &r_temp0[128 + 16*i];
    r_ptr2 = r_temp1 + 32*i;

    // repeat bound = 16
    v_zetas = vle16_v_i16m1(&zetas[32+16*i], 16);
    for(j = 0; j < 4; j++) {
      v_in2 = vle16_v_i16m1(r_ptr0, 16);
      v_in3 = vle16_v_i16m1(r_ptr1, 16);
      r_ptr0 += 32;
      r_ptr1 += 32;
      v_out1 = vbutterfly_ct_vvm_i16m2(v_in2, v_in3, v_zetas);
      vse16_v_i16m2(r_ptr2, v_out1, 32);
      r_ptr2 += 64;
    }
  }

  // indexed store: used as offsets
  vuint16m2_t v_offset;
  // tree pointer
  const uint16_t* tree_prt = NULL;

  // stage 7
  for(i = 0; i < 4; i++) {
    r_ptr0 = r_temp1 + 16*i;
    r_ptr1 = &r_temp1[128 + 16*i];
    tree_prt = byteoffset_even + 32*i;

    // repeat bound = 16
    v_zetas = vle16_v_i16m1(&zetas[64+16*i], 16);
    v_offset = vle16_v_u16m2(tree_prt, 32);
    for(j = 0; j < 2; j++) {
      v_in2 = vle16_v_i16m1(r_ptr0, 16);
      v_in3 = vle16_v_i16m1(r_ptr1, 16);
      r_ptr0 += 64;
      r_ptr1 += 64;
      v_out1 = vbutterfly_ct_vvm_i16m2(v_in2, v_in3, v_zetas);
      // reorder, from special bit-reversed order to standard order
      vsuxei16_v_i16m2(r, v_offset, v_out1, 32);
      v_offset = vadd_vx_u16m2(v_offset, 2, 32);
    }
  }

}

/*************************************************
* Name:        invntt_cg_zeta_compact
*
* Description: Inverse NTT.
*              Out-of-place, constant geometry.
*              Input coefficients are in standard order and normal domain.
*              Output coefficient are smaller than Q in absolute value;
*              and output coefficients are in bit-reversed order and normal domain.
*              zeta in compact form, that zetas in v0 are not the same
*
* Arguments:   - int16_t p[N]: input/output coefficient array
**************************************************/
void invntt_cg_zeta_compact(int16_t r[KYBER_N]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;
  const int16_t* zeta_ptr = zetas_inv_in_order;

  vint16m1_t v_zetas;
  vint16m1_t v_in0, v_in1;
  vint16m2_t v_out0;
  vint16m2_t v_in2, v_in3;
  vint16m4_t v_out1;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(1, 4));
  for(i = 0; i < 8; i++) {                // each 16* butterfly, 16/2=8 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 8);
    zeta_ptr += 8;
    v_in0 = vle16_v_i16m1(r_ptr0, 16);
    v_in1 = vle16_v_i16m1(r_ptr1, 16);
    r_ptr0 += 16;
    r_ptr1 += 16;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 32);
    r_ptr2 += 32;
  }

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // same_num = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(2, 4));  
  for(i = 0; i < 8; i++) {                // each 16* butterfly, 16/4=4 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 4);
    zeta_ptr += 4;
    v_in0 = vle16_v_i16m1(r_ptr0, 16);
    v_in1 = vle16_v_i16m1(r_ptr1, 16);
    r_ptr0 += 16;
    r_ptr1 += 16;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 32);
    r_ptr2 += 32;
  }

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(3, 4));
  for(i = 0; i < 8; i++) {                // each 16* butterfly, 16/8=2 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 2);
    zeta_ptr += 2;
    v_in0 = vle16_v_i16m1(r_ptr0, 16);
    v_in1 = vle16_v_i16m1(r_ptr1, 16);
    r_ptr0 += 16;
    r_ptr1 += 16;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 32);
    r_ptr2 += 32;
  }

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // same_num = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(4, 4));
  for(i = 0; i < 8; i++) {                // each 16* butterfly, 16/16=1 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    zeta_ptr += 1;
    v_in0 = vle16_v_i16m1(r_ptr0, 16);
    v_in1 = vle16_v_i16m1(r_ptr1, 16);
    r_ptr0 += 16;
    r_ptr1 += 16;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 32);
    r_ptr2 += 32;
  }

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(5, 4));
  for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/32=1 zeta
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    zeta_ptr += 1;
    v_in2 = vle16_v_i16m2(r_ptr0, 32);
    v_in3 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out1, 64);
    r_ptr2 += 64;
  }


  // stage 6
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r_temp1;
  // same_num = 64
  csr_zeta_selMode_rw(GET_ZETA_MODE(6, 4));
  for(i = 0; i < 2; i++) {                // each 64* butterfly, 64/64=1 zeta
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    zeta_ptr += 1;
    
    v_in2 = vle16_v_i16m2(r_ptr0, 32);
    v_in3 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out1, 64);
    r_ptr2 += 64;

    v_in2 = vle16_v_i16m2(r_ptr0, 32);
    v_in3 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out1, 64);
    r_ptr2 += 64;    
  }

  // stage 7
  r_ptr0 = r_temp1;
  r_ptr1 = &r_temp1[128];
  r_ptr2 = r;
  // same_num = 128
  csr_zeta_selMode_rw(GET_ZETA_MODE(7, 4));
  v_zetas = vle16_v_i16m1(zeta_ptr, 1);
  for(i = 0; i < 4; i++) {                // each 32* butterfly, all same zeta
    v_in2 = vle16_v_i16m2(r_ptr0, 32);
    v_in3 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out1, 64);
    r_ptr2 += 64;
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
* Arguments:   - int16_t p[N]: input/output coefficient array
**************************************************/
void invntt_cg_custom_reorder(int16_t r[KYBER_N]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i, j;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;
  const int16_t* zeta_ptr = zetas_inv_in_order;

  vint16m1_t v_zetas;
  vint16m1_t v_in0, v_in1;
  vint16m2_t v_out0;
  vint16m2_t v_in2, v_in3;
  vint16m4_t v_out1;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(1, 4));
  for(i = 0; i < 8; i++) {                // each 16* butterfly, 16/2=8 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 8);
    zeta_ptr += 8;
    v_in0 = vle16_v_i16m1(r_ptr0, 16);
    v_in1 = vle16_v_i16m1(r_ptr1, 16);
    r_ptr0 += 16;
    r_ptr1 += 16;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 32);
    r_ptr2 += 32;
  }

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // same_num = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(2, 4));
  for(i = 0; i < 8; i++) {                // each 16* butterfly, 16/4=4 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 4);
    zeta_ptr += 4;
    v_in0 = vle16_v_i16m1(r_ptr0, 16);
    v_in1 = vle16_v_i16m1(r_ptr1, 16);
    r_ptr0 += 16;
    r_ptr1 += 16;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 32);
    r_ptr2 += 32;
  }

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(3, 4));  
  for(i = 0; i < 8; i++) {                // each 16* butterfly, 16/8=2 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 2);
    zeta_ptr += 2;
    v_in0 = vle16_v_i16m1(r_ptr0, 16);
    v_in1 = vle16_v_i16m1(r_ptr1, 16);
    r_ptr0 += 16;
    r_ptr1 += 16;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 32);
    r_ptr2 += 32;
  }

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // same_num = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(4, 4));  
  for(i = 0; i < 8; i++) {                // each 16* butterfly, 16/16=1 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    zeta_ptr += 1;
    v_in0 = vle16_v_i16m1(r_ptr0, 16);
    v_in1 = vle16_v_i16m1(r_ptr1, 16);
    r_ptr0 += 16;
    r_ptr1 += 16;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 32);
    r_ptr2 += 32;
  }

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(5, 4));
  for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/32=1 zeta
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    zeta_ptr += 1;
    v_in2 = vle16_v_i16m2(r_ptr0, 32);
    v_in3 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out1, 64);
    r_ptr2 += 64;
  }


  // stage 6
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r_temp1;
  // same_num = 64
  csr_zeta_selMode_rw(GET_ZETA_MODE(6, 4));
  for(i = 0; i < 2; i++) {                // each 64* butterfly, 64/64=1 zeta
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    zeta_ptr += 1;
    
    v_in2 = vle16_v_i16m2(r_ptr0, 32);
    v_in3 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out1, 64);
    r_ptr2 += 64;

    v_in2 = vle16_v_i16m2(r_ptr0, 32);
    v_in3 = vle16_v_i16m2(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out1, 64);
    r_ptr2 += 64;    
  }

  // indexed store: used as offsets
  vuint16m4_t v_offset;
  // tree pointer
  const uint16_t* tree_prt = NULL;

  // stage 7
  // same_num = 128
  csr_zeta_selMode_rw(GET_ZETA_MODE(7, 4));
  v_zetas = vle16_v_i16m1(zeta_ptr, 1);
  for(i = 0; i < 2; i++) {              // each 32* butterfly, all same zeta
    r_ptr0 = r_temp1 + 32*i;
    r_ptr1 = &r_temp1[128 + 32*i];
    tree_prt = byteoffset_even + 64*i;

    v_offset = vle16_v_u16m4(tree_prt, 64);
    for(j = 0; j < 2; j++) {
      v_in2 = vle16_v_i16m2(r_ptr0, 32);
      v_in3 = vle16_v_i16m2(r_ptr1, 32);
      r_ptr0 += 64;
      r_ptr1 += 64;
      v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
      // reorder, from special bit-reversed order to standard order
      vsuxei16_v_i16m4(r, v_offset, v_out1, 64);
      v_offset = vadd_vx_u16m4(v_offset, 2, 64);      
    }    
  }
  
}

#elif (VLEN == 512)
/*************************************************
* Name:        ntt_cg_custom_zeta_compact
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in bitreversed order
*              zeta in compact form, that zetas in v0 are not the same
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_zeta_compact(int16_t r[256]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i, j;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;

  vint16m1_t v_zetas;
  // for stage 1 to stage 6
  vint16m2_t v_in0, v_in1;
  vint16m4_t v_out0;
  // for stage 7
  vint16m1_t v_in2, v_in3;
  vint16m2_t v_out1;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 1
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  v_zetas = vle16_v_i16m1(&zetas[1], 1);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));  
  v_zetas = vle16_v_i16m1(&zetas[2], 2);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
  v_zetas = vle16_v_i16m1(&zetas[4], 4);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
  v_zetas = vle16_v_i16m1(&zetas[8], 8);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
  v_zetas = vle16_v_i16m1(&zetas[16], 16);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }


  // stage 6
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r_temp1;
  // repeat bound = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));
  v_zetas = vle16_v_i16m1(&zetas[32], 32);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }

  // Vector LMUL Truncation
  v_in2 = vget_v_i16m2_i16m1(v_in0, 0);
  v_in3 = vget_v_i16m2_i16m1(v_in1, 0);
  v_out1 = vget_v_i16m4_i16m2(v_out0, 0);

  // stage 7
  for(i = 0; i < 2; i++) {
    r_ptr0 = r_temp1 + 32*i;
    r_ptr1 = &r_temp1[128 + 32*i];
    r_ptr2 = r + 64*i;

    // repeat bound = 32
    v_zetas = vle16_v_i16m1(&zetas[64+32*i], 32);
    for(j = 0; j < 2; j++) {
      v_in2 = vle16_v_i16m1(r_ptr0, 32);
      v_in3 = vle16_v_i16m1(r_ptr1, 32);
      r_ptr0 += 64;
      r_ptr1 += 64;
      v_out1 = vbutterfly_ct_vvm_i16m2(v_in2, v_in3, v_zetas);
      vse16_v_i16m2(r_ptr2, v_out1, 64);
      r_ptr2 += 128;
    }
  }

}

/*************************************************
* Name:        ntt_cg_custom_reorder
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_reorder(int16_t r[256]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i, j;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;

  vint16m1_t v_zetas;
  // for stage 1 to stage 6
  vint16m2_t v_in0, v_in1;
  vint16m4_t v_out0;
  // for stage 7
  vint16m1_t v_in2, v_in3;
  vint16m2_t v_out1;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 1
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  v_zetas = vle16_v_i16m1(&zetas[1], 1);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));  
  v_zetas = vle16_v_i16m1(&zetas[2], 2);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
  v_zetas = vle16_v_i16m1(&zetas[4], 4);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
  v_zetas = vle16_v_i16m1(&zetas[8], 8);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
  v_zetas = vle16_v_i16m1(&zetas[16], 16);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }


  // stage 6
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r_temp1;
  // repeat bound = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));
  v_zetas = vle16_v_i16m1(&zetas[32], 32);
  for(i = 0; i < 2; i++) {  // 2 * vectorized butterfly (each 64 set of bf) for 128 sets of bf totally
    v_in0 = vle16_v_i16m2(r_ptr0, 64);
    v_in1 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out0, 128);
    r_ptr2 += 128;
  }


  // Vector LMUL Truncation
  v_in2 = vget_v_i16m2_i16m1(v_in0, 0);
  v_in3 = vget_v_i16m2_i16m1(v_in1, 0);
  v_out1 = vget_v_i16m4_i16m2(v_out0, 0);

  // indexed store: used as offsets
  vuint16m2_t v_offset;
  // tree pointer
  const uint16_t* tree_prt = NULL;

  // stage 7
  for(i = 0; i < 2; i++) {
    r_ptr0 = r_temp1 + 32*i;
    r_ptr1 = &r_temp1[128 + 32*i];
    tree_prt = byteoffset_even + 64*i;

    // repeat bound = 32
    v_zetas = vle16_v_i16m1(&zetas[64+32*i], 32);
    v_offset = vle16_v_u16m2(tree_prt, 64);
    for(j = 0; j < 2; j++) {
      v_in2 = vle16_v_i16m1(r_ptr0, 32);
      v_in3 = vle16_v_i16m1(r_ptr1, 32);
      r_ptr0 += 64;
      r_ptr1 += 64;
      v_out1 = vbutterfly_ct_vvm_i16m2(v_in2, v_in3, v_zetas);
      // reorder, from special bit-reversed order to standard order
      vsuxei16_v_i16m2(r, v_offset, v_out1, 64);
      v_offset = vadd_vx_u16m2(v_offset, 2, 64);
    }
  }

}

/*************************************************
* Name:        ntt_cg_custom_asm
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in bitreversed order
*              assembly optimized
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_asm(int16_t r[256]) {
  size_t avl = KYBER_N;
  size_t load_zeta_vl;

  // load coefficients from memory to vector register v8~v15
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
    "vle16.v v8, (%[coeff_base])\n"
    :
    : [coeff_base] "r"(r), [avl] "r"(avl)
  );  

  // stage 1
  // v8 to v15 with m4 is input coefficients, v8~v11 is top half, v12~v15 is down half
  // v16 to v23 is output coefficients
  // repeat bound = 1
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  load_zeta_vl = 1;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v16, v8, v12, v0\n"
    :
    : [zeta_base] "r"(&zetas[1]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 2
  // v16 to v23 with m4 is input coefficients, v16~v19 is top half, v20~v23 is down half
  // v8 to v15 is output coefficients
  // repeat bound = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));  
  load_zeta_vl = 2;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v8, v16, v20, v0\n"
    :
    : [zeta_base] "r"(&zetas[2]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 3
  // v8 to v15 with m4 is input coefficients, v8~v11 is top half, v12~v15 is down half
  // v16 to v23 is output coefficients
  // repeat bound = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
  load_zeta_vl = 4;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v16, v8, v12, v0\n"
    :
    : [zeta_base] "r"(&zetas[4]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 4
  // v16 to v23 with m4 is input coefficients, v16~v19 is top half, v20~v23 is down half
  // v8 to v15 is output coefficients
  // repeat bound = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));  
  load_zeta_vl = 8;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v8, v16, v20, v0\n"
    :
    : [zeta_base] "r"(&zetas[8]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );
  
  // stage 5
  // v8 to v15 with m4 is input coefficients, v8~v11 is top half, v12~v15 is down half
  // v16 to v23 is output coefficients
  // repeat bound = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));  
  load_zeta_vl = 16;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v16, v8, v12, v0\n"
    :
    : [zeta_base] "r"(&zetas[16]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 6
  // v16 to v23 with m4 is input coefficients, v16~v19 is top half, v20~v23 is down half
  // v8 to v15 is output coefficients
  // repeat bound = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));  
  load_zeta_vl = 32;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v8, v16, v20, v0\n"
    :
    : [zeta_base] "r"(&zetas[32]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 7
  // v8 to v15 with m1 is input coefficients, v8~v11 is top half, v12~v15 is down half
  // (v8, v12)->v16 and (v10, v14)->v20 shares the same zetas
  // (v9, v13)->v18 and (v11, v15)->v22 shares the same zetas
  // v16 to v23 is output coefficients
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vbutterfly.ct.vvm v16, v8, v12, v0\n"
    "vbutterfly.ct.vvm v20, v10, v14, v0\n"
    :
    : [zeta_base] "r"(&zetas[64]), [zeta_vl] "r"(load_zeta_vl)
  );  
  __asm__ __volatile__ (
    "vle16.v v0, (%[zeta_base])\n"
    "vbutterfly.ct.vvm v18, v9, v13, v0\n"
    "vbutterfly.ct.vvm v22, v11, v15, v0\n"
    :
    : [zeta_base] "r"(&zetas[96]), [zeta_vl] "r"(load_zeta_vl)
  );

  // store coefficients from vector register v16~v23 to memory
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
    "vse16.v v16, (%[coeff_base])\n"
    :
    : [coeff_base] "r"(r), [avl] "r"(avl)
  );
}

/*************************************************
* Name:        ntt_cg_custom_reorder_asm
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*              assembly optimized
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_reorder_asm(int16_t r[256]) {
  size_t avl = KYBER_N;
  size_t load_zeta_vl;

  // load coefficients from memory to vector register v8~v15
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
    "vle16.v v8, (%[coeff_base])\n"
    :
    : [coeff_base] "r"(r), [avl] "r"(avl)
  );

  // stage 1
  // v8 to v15 with m4 is input coefficients, v8~v11 is top half, v12~v15 is down half
  // v16 to v23 is output coefficients
  // repeat bound = 1
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  load_zeta_vl = 1;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v16, v8, v12, v0\n"
    :
    : [zeta_base] "r"(&zetas[1]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 2
  // v16 to v23 with m4 is input coefficients, v16~v19 is top half, v20~v23 is down half
  // v8 to v15 is output coefficients
  // repeat bound = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));  
  load_zeta_vl = 2;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v8, v16, v20, v0\n"
    :
    : [zeta_base] "r"(&zetas[2]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 3
  // v8 to v15 with m4 is input coefficients, v8~v11 is top half, v12~v15 is down half
  // v16 to v23 is output coefficients
  // repeat bound = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
  load_zeta_vl = 4;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v16, v8, v12, v0\n"
    :
    : [zeta_base] "r"(&zetas[4]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 4
  // v16 to v23 with m4 is input coefficients, v16~v19 is top half, v20~v23 is down half
  // v8 to v15 is output coefficients
  // repeat bound = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));  
  load_zeta_vl = 8;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v8, v16, v20, v0\n"
    :
    : [zeta_base] "r"(&zetas[8]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );
  
  // stage 5
  // v8 to v15 with m4 is input coefficients, v8~v11 is top half, v12~v15 is down half
  // v16 to v23 is output coefficients
  // repeat bound = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));  
  load_zeta_vl = 16;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v16, v8, v12, v0\n"
    :
    : [zeta_base] "r"(&zetas[16]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 6
  // v16 to v23 with m4 is input coefficients, v16~v19 is top half, v20~v23 is down half
  // v8 to v15 is output coefficients
  // repeat bound = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));  
  load_zeta_vl = 32;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vbutterfly.ct.vvm v8, v16, v20, v0\n"
    :
    : [zeta_base] "r"(&zetas[32]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 7
  // v8 to v15 with m1 is input coefficients, v8~v11 is top half, v12~v15 is down half
  // (v8, v12)->v16 and (v10, v14)->v20 shares the same zetas
  // (v9, v13)->v18 and (v11, v15)->v22 shares the same zetas
  // v16 to v23 is output coefficients
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vbutterfly.ct.vvm v16, v8, v12, v0\n"
    "vbutterfly.ct.vvm v20, v10, v14, v0\n"
    :
    : [zeta_base] "r"(&zetas[64]), [zeta_vl] "r"(load_zeta_vl)
  );  
  __asm__ __volatile__ (
    "vle16.v v0, (%[zeta_base])\n"
    "vbutterfly.ct.vvm v18, v9, v13, v0\n"
    "vbutterfly.ct.vvm v22, v11, v15, v0\n"
    :
    : [zeta_base] "r"(&zetas[96]), [zeta_vl] "r"(load_zeta_vl)
  );

  // indexed store coefficients from vector register v16~v23 to memory
  const uint16_t* tree_prt = NULL;

  tree_prt = byteoffset_even;
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vle16.v v2, (%[offset_base])\n"
    "vsuxei16.v v16, (%[coeff_base]), v2\n"     // indexed store v16 and v17
    "vadd.vi v2, v2, %[offset_plus]\n"
    "vsuxei16.v v20, (%[coeff_base]), v2\n"     // indexed store v20 and v21
    :
    : [coeff_base] "r"(r), [avl] "r"(64), [offset_plus] "i"(2), [offset_base] "r"(tree_prt)
  );

  tree_prt = byteoffset_even + 64;
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vle16.v v2, (%[offset_base])\n"
    "vsuxei16.v v18, (%[coeff_base]), v2\n"     // indexed store v18 and v19
    "vadd.vi v2, v2, %[offset_plus]\n"
    "vsuxei16.v v22, (%[coeff_base]), v2\n"     // indexed store v22 and v23
    :
    : [coeff_base] "r"(r), [avl] "r"(64), [offset_plus] "i"(2), [offset_base] "r"(tree_prt)
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
* Arguments:   - int16_t p[N]: input/output coefficient array
**************************************************/
void invntt_cg_custom_reorder(int16_t r[KYBER_N]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i, j;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;
  const int16_t* zeta_ptr = zetas_inv_in_order;

  vint16m1_t v_zetas;
  vint16m1_t v_in0, v_in1;
  vint16m2_t v_out0;
  vint16m2_t v_in2, v_in3;
  vint16m4_t v_out1;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(1, 5));
  for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/2=16 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 16);
    zeta_ptr += 16;
    v_in0 = vle16_v_i16m1(r_ptr0, 32);
    v_in1 = vle16_v_i16m1(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // same_num = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(2, 5));
  for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/4=8 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 8);
    zeta_ptr += 8;
    v_in0 = vle16_v_i16m1(r_ptr0, 32);
    v_in1 = vle16_v_i16m1(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(3, 5));  
  for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/8=4 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 4);
    zeta_ptr += 4;
    v_in0 = vle16_v_i16m1(r_ptr0, 32);
    v_in1 = vle16_v_i16m1(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // same_num = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(4, 5));  
  for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/16=2 zetas
    v_zetas = vle16_v_i16m1(zeta_ptr, 2);
    zeta_ptr += 2;
    v_in0 = vle16_v_i16m1(r_ptr0, 32);
    v_in1 = vle16_v_i16m1(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(5, 5));
  for(i = 0; i < 4; i++) {                // each 32* butterfly, all same zeta
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    zeta_ptr += 1;
    v_in0 = vle16_v_i16m1(r_ptr0, 32);
    v_in1 = vle16_v_i16m1(r_ptr1, 32);
    r_ptr0 += 32;
    r_ptr1 += 32;
    v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
    vse16_v_i16m2(r_ptr2, v_out0, 64);
    r_ptr2 += 64;
  }

  // stage 6
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r_temp1;
  // same_num = 64
  csr_zeta_selMode_rw(GET_ZETA_MODE(6, 5));
  for(i = 0; i < 2; i++) {              // each 64* butterfly, all same zeta
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    zeta_ptr += 1;
    
    v_in2 = vle16_v_i16m2(r_ptr0, 64);
    v_in3 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
    vse16_v_i16m4(r_ptr2, v_out1, 128);
    r_ptr2 += 128;
  }

  // indexed store: used as offsets
  vuint16m4_t v_offset;
  // tree pointer
  const uint16_t* tree_prt = NULL;

  // stage 7
  // same_num = 128
  csr_zeta_selMode_rw(GET_ZETA_MODE(7, 5));
  v_zetas = vle16_v_i16m1(zeta_ptr, 1);
  r_ptr0 = r_temp1;
  r_ptr1 = &r_temp1[128];
  tree_prt = byteoffset_even;
  v_offset = vle16_v_u16m4(tree_prt, 128);
  for(j = 0; j < 2; j++) {                // each 64* butterfly, all same zeta
    v_in2 = vle16_v_i16m2(r_ptr0, 64);
    v_in3 = vle16_v_i16m2(r_ptr1, 64);
    r_ptr0 += 64;
    r_ptr1 += 64;
    v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
    // reorder, from special bit-reversed order to standard order
    vsuxei16_v_i16m4(r, v_offset, v_out1, 128);
    v_offset = vadd_vx_u16m4(v_offset, 2, 128);
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
* Arguments:   - int16_t p[N]: input/output coefficient array
**************************************************/
void invntt_cg_custom_reorder_asm(int16_t r[256]) {
  size_t avl = KYBER_N;
  size_t load_zeta_vl, bf_vl;
  const int16_t* zeta_ptr = zetas_inv_in_order;

  // load coefficients from memory to vector register v24~v31
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
    "vle16.v v24, (%[coeff_base])\n"
    :
    : [coeff_base] "r"(r), [avl] "r"(avl)
  );

  // from stage 1 to stage 5, bf_vl=32 
  bf_vl = 32;

  // stage 1
  // v24 to v31 with m1 is input coefficients, v24~v27 is top half, v28~v31 is down half
  // v16 to v23 is output coefficients
  // same_num = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(1, 5));
  load_zeta_vl = 16;
  // (v24, v28) -> v16~v17
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v16, v24, v28, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v25, v29) -> v18~v19
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v18, v25, v29, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v26, v30) -> v20~v21
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v20, v26, v30, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v27, v31) -> v22~v23
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v22, v27, v31, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;  


  // stage 2
  // v16 to v23 with m1 is input coefficients, v16~v19 is top half, v20~v23 is down half
  // v24 to v31 is output coefficients
  // same_num = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(2, 5));
  load_zeta_vl = 8;
  // (v16, v20) -> v24~v25
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v24, v16, v20, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v17, v21) -> v26~v27
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v26, v17, v21, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;  
  // (v18, v22) -> v28~v29
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v28, v18, v22, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;  
  // (v19, v23) -> v30~v31
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v30, v19, v23, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;


  // stage 3
  // v24 to v31 with m1 is input coefficients, v24~v27 is top half, v28~v31 is down half
  // v16 to v23 is output coefficients
  // same_num = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(3, 5));  
  load_zeta_vl = 4;
  // (v24, v28) -> v16~v17
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v16, v24, v28, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v25, v29) -> v18~v19
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v18, v25, v29, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v26, v30) -> v20~v21
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v20, v26, v30, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v27, v31) -> v22~v23
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v22, v27, v31, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;

  // stage 4
  // v16 to v23 with m1 is input coefficients, v16~v19 is top half, v20~v23 is down half
  // v24 to v31 is output coefficients
  // same_num = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(4, 5));  
  load_zeta_vl = 2;
  // (v16, v20) -> v24~v25
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v24, v16, v20, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v17, v21) -> v26~v27
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v26, v17, v21, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;  
  // (v18, v22) -> v28~v29
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v28, v18, v22, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;  
  // (v19, v23) -> v30~v31
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v30, v19, v23, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  
  // stage 5
  // v24 to v31 with m1 is input coefficients, v24~v27 is top half, v28~v31 is down half
  // v16 to v23 is output coefficients
  // same_num = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(5, 5)); 
  load_zeta_vl = 1;
  // (v24, v28) -> v16~v17
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v16, v24, v28, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v25, v29) -> v18~v19
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v18, v25, v29, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v26, v30) -> v20~v21
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v20, v26, v30, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v27, v31) -> v22~v23
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m1, ta, mu\n"
    "vbutterfly.gs.vvm v22, v27, v31, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;

  // from stage 6 to stage 7, bf_vl=64
  bf_vl = 64;

  // stage 6
  // v16 to v23 with m2 is input coefficients, v16~v19 is top half, v20~v23 is down half
  // v24 to v31 is output coefficients
  // same_num = 64
  csr_zeta_selMode_rw(GET_ZETA_MODE(6, 5)); 
  load_zeta_vl = 1;
  // (v16~v17, v20~v21) -> v24~v27
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m2, ta, mu\n"
    "vbutterfly.gs.vvm v24, v16, v20, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;
  // (v18~v19, v22~v23) -> v28~v31
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m2, ta, mu\n"
    "vbutterfly.gs.vvm v28, v18, v22, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );
  zeta_ptr += load_zeta_vl;

  // stage 7
  // v24 to v31 with m2 is input coefficients, v24~v27 is top half, v28~v31 is down half
  // v16 to v23 is output coefficients
  // same_num = 128
  csr_zeta_selMode_rw(GET_ZETA_MODE(7, 5));
  // (v24~v25, v28~v29) -> v16~v19
  // (v26~v27, v30~v31) -> v20~v23
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[bf_vl], e16, m2, ta, mu\n"
    "vbutterfly.gs.vvm v16, v24, v28, v0\n"
    "vbutterfly.gs.vvm v20, v26, v30, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(load_zeta_vl), [bf_vl] "r"(bf_vl)
  );    


  // Finally,
  // indexed store coefficients from vector register v16~v23 to memory
  const uint16_t* tree_prt = NULL;

  tree_prt = byteoffset_even;
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vle16.v v2, (%[offset_base])\n"
    "vsuxei16.v v16, (%[coeff_base]), v2\n"     // indexed store v16 and v17
    "vadd.vi v2, v2, %[offset_plus]\n"
    "vsuxei16.v v20, (%[coeff_base]), v2\n"     // indexed store v20 and v21
    :
    : [coeff_base] "r"(r), [avl] "r"(64), [offset_plus] "i"(2), [offset_base] "r"(tree_prt)
  );

  tree_prt = byteoffset_even + 64;
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vle16.v v2, (%[offset_base])\n"
    "vsuxei16.v v18, (%[coeff_base]), v2\n"     // indexed store v18 and v19
    "vadd.vi v2, v2, %[offset_plus]\n"
    "vsuxei16.v v22, (%[coeff_base]), v2\n"     // indexed store v22 and v23
    :
    : [coeff_base] "r"(r), [avl] "r"(64), [offset_plus] "i"(2), [offset_base] "r"(tree_prt)
  );
}

#elif (VLEN == 1024)
/*************************************************
* Name:        ntt_cg_custom_reorder
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_reorder(int16_t r[256]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i, j;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;

  vint16m1_t v_zetas;
  // for stage 1 to stage 7
  vint16m2_t v_in0, v_in1;
  vint16m4_t v_out0;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 1
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  v_zetas = vle16_v_i16m1(&zetas[1], 1);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));  
  v_zetas = vle16_v_i16m1(&zetas[2], 2);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
  v_zetas = vle16_v_i16m1(&zetas[4], 4);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // repeat bound = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));
  v_zetas = vle16_v_i16m1(&zetas[8], 8);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // repeat bound = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));
  v_zetas = vle16_v_i16m1(&zetas[16], 16);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);


  // stage 6
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r_temp1;
  // repeat bound = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));
  v_zetas = vle16_v_i16m1(&zetas[32], 32);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);

  // stage 7
  // indexed store: used as offsets
  vuint16m4_t v_offset;
  // tree pointer
  const uint16_t* tree_prt = NULL;
  tree_prt = byteoffset_total;
  r_ptr0 = r_temp1;
  r_ptr1 = &r_temp1[128];
  //repeat bound=64
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 6));
  v_zetas = vle16_v_i16m1(&zetas[64], 64);
  v_offset = vle16_v_u16m4(tree_prt, 256);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_ct_vvm_i16m4(v_in0, v_in1, v_zetas);
  vsuxei16_v_i16m4(r, v_offset, v_out0, 256);

}

/*************************************************
* Name:        ntt_cg_custom_reorder_asm
*
* Description: Out-of-place number-theoretic transform (NTT) in Rq
*              input is in standard order, output is in standard order
*              assembly optimized
*
* Arguments:   - int16_t r[256]: pointer to input/output vector of elements
*                                of Zq
**************************************************/
void ntt_cg_custom_reorder_asm(int16_t r[256]) {
  size_t avl = KYBER_N;
  size_t load_zeta_vl;

  // load coefficients from memory to vector register v4~v7
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vle16.v v4, (%[coeff_base])\n"
    :
    : [coeff_base] "r"(r), [avl] "r"(avl)
  );

  // stage 1
  // v4 to v7 with m2 is input coefficients, v4 v5 is top half, v6 v7 is down half
  // v8 to v11 is output coefficients
  // repeat bound = 1
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 0));
  load_zeta_vl = 1;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vbutterfly.ct.vvm v8, v4, v6, v0\n"
    :
    : [zeta_base] "r"(&zetas[1]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 2
  // v8 to v11 with m2 is input coefficients, v8~v9 is top half, v10~v11 is down half
  // v4 to v7 is output coefficients
  // repeat bound = 2
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 1));  
  load_zeta_vl = 2;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vbutterfly.ct.vvm v4, v8, v10, v0\n"
    :
    : [zeta_base] "r"(&zetas[2]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 3
  // v4 to v7 with m2 is input coefficients, v4~v5 is top half, v6~v7 is down half
  // v8 to v11 is output coefficients
  // repeat bound = 4
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 2));
  load_zeta_vl = 4;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vbutterfly.ct.vvm v8, v4, v6, v0\n"
    :
    : [zeta_base] "r"(&zetas[4]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 4
  // v8 to v11 with m2 is input coefficients, v8~v9 is top half, v10~v11 is down half
  // v4 to v7 is output coefficients
  // repeat bound = 8
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 3));  
  load_zeta_vl = 8;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vbutterfly.ct.vvm v4, v8, v10, v0\n"
    :
    : [zeta_base] "r"(&zetas[8]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );
  
  // stage 5
  // v4 to v7 with m2 is input coefficients, v4~v5 is top half, v6~v7 is down half
  // v8 to v11 is output coefficients
  // repeat bound = 16
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 4));  
  load_zeta_vl = 16;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vbutterfly.ct.vvm v8, v4, v6, v0\n"
    :
    : [zeta_base] "r"(&zetas[16]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 6
  // v8 to v11 with m2 is input coefficients, v8~v9 is top half, v10~v11 is down half
  // v4 to v7 is output coefficients
  // repeat bound = 32
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));  
  load_zeta_vl = 32;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vbutterfly.ct.vvm v4, v8, v10, v0\n"
    :
    : [zeta_base] "r"(&zetas[32]), [zeta_vl] "r"(load_zeta_vl), [avl] "r"(avl)
  );

  // stage 7
  // v4 to v7 with m1 is input coefficients, v4~v5 is top half, v6~v7 is down half
  // v8 to v11 is output coefficients
  csr_zeta_selMode_rw(GET_ZETA_MODE(0, 6)); 
  load_zeta_vl = 64;
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vbutterfly.ct.vvm v8, v4, v6, v0\n"
    :
    : [zeta_base] "r"(&zetas[64]), [avl] "r"(avl), [zeta_vl] "r"(load_zeta_vl)
  );  

  // indexed store coefficients from vector register v8-v11 to memory
  const uint16_t* tree_prt = NULL;

  tree_prt = byteoffset_total;
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m4, ta, mu\n"
    "vle16.v v4, (%[offset_base])\n"
    "vsuxei16.v v8, (%[coeff_base]), v4\n"
    :
    : [coeff_base] "r"(r), [avl] "r"(avl), [offset_base] "r"(tree_prt)
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
* Arguments:   - int16_t p[N]: input/output coefficient array
**************************************************/
void invntt_cg_custom_reorder(int16_t r[KYBER_N]) {
  int16_t r_temp0[256];
  int16_t r_temp1[256];
  unsigned int i, j;
  int16_t* r_ptr0 = NULL;
  int16_t* r_ptr1 = NULL;
  int16_t* r_ptr2 = NULL;
  const int16_t* zeta_ptr = zetas_inv_in_order;

  vint16m1_t v_zetas;
  vint16m2_t v_in0, v_in1;
  vint16m4_t v_out0;

  // stage 1
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 2, 128 butterfly, 128/2=64 zetas
  csr_zeta_selMode_rw(GET_ZETA_MODE(1, 8));
  v_zetas = vle16_v_i16m1(zeta_ptr, 64);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_gs_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);
  zeta_ptr+=64;

  // stage 2
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // same_num = 4, 128 butterfly, 128/4=32 zetas
  csr_zeta_selMode_rw(GET_ZETA_MODE(2, 8));
  v_zetas = vle16_v_i16m1(zeta_ptr, 32);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_gs_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);
  zeta_ptr+=32;

  // stage 3
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 8, 128 butterfly, 128/8=16 zetas
  csr_zeta_selMode_rw(GET_ZETA_MODE(3, 8)); 
  v_zetas = vle16_v_i16m1(zeta_ptr, 16);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_gs_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);
  zeta_ptr+=16;

  // stage 4
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r;
  // same_num = 16, 128 butterfly, 128/16=8 zetas
  csr_zeta_selMode_rw(GET_ZETA_MODE(4, 8));  
  v_zetas = vle16_v_i16m1(zeta_ptr, 8);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_gs_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);
  zeta_ptr+=8;

  // stage 5
  r_ptr0 = r;
  r_ptr1 = &r[128];
  r_ptr2 = r_temp0;
  // same_num = 32, 128 butterfly, 128/32=4 zetas
  csr_zeta_selMode_rw(GET_ZETA_MODE(5, 8));
  v_zetas = vle16_v_i16m1(zeta_ptr, 4);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_gs_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);
  zeta_ptr+=4;

  // stage 6
  r_ptr0 = r_temp0;
  r_ptr1 = &r_temp0[128];
  r_ptr2 = r_temp1;
  // same_num = 64
  csr_zeta_selMode_rw(GET_ZETA_MODE(6, 8));
  v_zetas = vle16_v_i16m1(zeta_ptr, 2);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_gs_vvm_i16m4(v_in0, v_in1, v_zetas);
  vse16_v_i16m4(r_ptr2, v_out0, 256);
  zeta_ptr+=2;

  // indexed store: used as offsets
  vuint16m4_t v_offset;
  // tree pointer
  const uint16_t* tree_prt = NULL;

  // stage 7
  // same_num = 128
  csr_zeta_selMode_rw(GET_ZETA_MODE(7, 8));
  v_zetas = vle16_v_i16m1(zeta_ptr, 1);
  r_ptr0 = r_temp1;
  r_ptr1 = &r_temp1[128];
  tree_prt = byteoffset_total;
  v_offset = vle16_v_u16m4(tree_prt, 256);
  v_in0 = vle16_v_i16m2(r_ptr0, 128);
  v_in1 = vle16_v_i16m2(r_ptr1, 128);
  v_out0 = vbutterfly_gs_vvm_i16m4(v_in0, v_in1, v_zetas);
  vsuxei16_v_i16m4(r, v_offset, v_out0, 256);  
}
# else
#error "VLEN must be 256/512/1024"
#endif