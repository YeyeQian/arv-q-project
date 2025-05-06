#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include "fips202.h"
#include "params.h"
#include "symmetric.h"
#include "poly.h"
#include "polyvec.h"
#include "ntt.h"
#include "reduce.h"
#include "cbd.h"
#include "poly_reorg.h"

/*************************************************
* Name:        poly_eta1_add_ntt
*
* Description: 1. Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA1
*              2. Modular add q to each coefficients in a polynomial
*              in this way turn each coefficients to non-negetive number
*              3. Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in normal order  
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length KYBER_SYMBYTES bytes)
*              - uint8_t nonce:       one-byte input nonce
**************************************************/
void poly_eta1_add_ntt(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce)
{
	uint8_t buf[KYBER_ETA1*KYBER_N/4];
  kyber_shake256_prf_custom(buf, sizeof(buf), seed, nonce);
  cbd_eta1_custom(r, buf);
	poly_mod_add_q(r);
	poly_ntt_custom(r);
}

/*************************************************
* Name:        poly_eta1_add_ntt_asm
*
* Description: 1. Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA1
*              2. Modular add q to each coefficients in a polynomial
*              in this way turn each coefficients to non-negetive number
*              3. Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place;
*              inputs assumed to be in normal order, output in normal order  
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length KYBER_SYMBYTES bytes)
*              - uint8_t nonce:       one-byte input nonce
**************************************************/
void poly_eta1_add_ntt_asm(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce)
{
  uint8_t buf[KYBER_ETA1*KYBER_N/4];
  kyber_shake256_prf_custom(buf, sizeof(buf), seed, nonce);

  /*
	** Firstly, binomial sampling process
  ** Use vector register v8 ~ v15 to contain a polynomial in binomial sampling
	** Use vector register v2 and v4 to contian cbd_even_index and cbd_odd_index
	** Use v16 to contain buf, use v17 to contian unpacked data and vcpopv
	** Use v18 to contain even vrgather, v19 to contain odd vrgather
  */
 	uint8_t* buf_ptr = buf;
#if KYBER_ETA1 == 2
	size_t vl_packed = ELEMENT_SEW8_PER_VECREG >> 2;
	// load v_index0 and v_index1 to v2 and v4
  __asm__ __volatile__ (
    "vsetvli	zero, %[avl], e8, m1, ta, mu\n"
    "vle16.v v2, (%[even_addr])\n"
    "vle16.v v4, (%[odd_addr])\n"
    :
    : [even_addr] "r"(cbd_even_index), [odd_addr] "r"(cbd_odd_index), [avl] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );

	// get v8
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v8, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v9
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v9, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v10
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v10, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v11
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v11, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v12
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v12, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v13
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v13, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v14
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v14, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v15
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v15, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );

#elif KYBER_ETA1 == 3
	size_t vl_packed = (ELEMENT_SEW8_PER_VECREG * 3) >> 3;
	// load v_index0 and v_index1 to v2 and v4
  __asm__ __volatile__ (
    "vsetvli	zero, %[avl], e8, m1, ta, mu\n"
    "vle16.v v2, (%[even_addr])\n"
    "vle16.v v4, (%[odd_addr])\n"
    :
    : [even_addr] "r"(cbd_even_index), [odd_addr] "r"(cbd_odd_index), [avl] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );

	// get v8
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v8, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(3), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v9
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v9, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(3), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v10
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v10, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(3), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v11
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v11, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(3), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v12
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v12, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(3), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v13
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v13, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(3), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v14
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v14, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(3), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v15
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v15, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(3), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
#else
#error "This implementation requires eta1 in {2,3}"
#endif

  /*
	** Secondly, modular add q process
	** Input coefficients are in v8 ~ v15
	** Result are also in v8 ~ v15
  */
  __asm__ __volatile__ (
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
		"vmod.add.vx v8, v8, %[q]\n"
    :
    : [avl] "r"(KYBER_N), [q] "r"(KYBER_Q)
  ); 	

 	/*
	** Third, NTT with reordered store
	** Input coefficients are in v8 ~ v15
	** Result are also in v16 ~ v23
  */
  size_t avl = KYBER_N;
  size_t load_zeta_vl;

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
	int16_t* coeff_ptr = r->coeffs;

  tree_prt = byteoffset_even;
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vle16.v v2, (%[offset_base])\n"
    "vsuxei16.v v16, (%[coeff_base]), v2\n"     // indexed store v16 and v17
    "vadd.vi v2, v2, %[offset_plus]\n"
    "vsuxei16.v v20, (%[coeff_base]), v2\n"     // indexed store v20 and v21
    :
    : [coeff_base] "r"(coeff_ptr), [avl] "r"(64), [offset_plus] "i"(2), [offset_base] "r"(tree_prt)
  );

  tree_prt = byteoffset_even + 64;
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m2, ta, mu\n"
    "vle16.v v2, (%[offset_base])\n"
    "vsuxei16.v v18, (%[coeff_base]), v2\n"     // indexed store v18 and v19
    "vadd.vi v2, v2, %[offset_plus]\n"
    "vsuxei16.v v22, (%[coeff_base]), v2\n"     // indexed store v22 and v23
    :
    : [coeff_base] "r"(coeff_ptr), [avl] "r"(64), [offset_plus] "i"(2), [offset_base] "r"(tree_prt)
  );
}

/*************************************************
* Name:        polyvec_base_mul_acc_tomont_custom
*
* Description: 1. Pointwise multiply elements of a and b, accumulate into r,
*              and multiply by 2^-16.
*              2. Inplace conversion of all coefficients of a polynomial
*              from normal domain to Montgomery domain
*
* Arguments: - poly *r:          pointer to output polynomial
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_base_mul_acc_tomont_custom(poly *r,
                                      const polyvec *a,
                                      const polyvec *b)
{
  unsigned int i;
  poly t;

  poly_basemul_montgomery_custom(r, &a->vec[0], &b->vec[0]);
  for(i=1;i<KYBER_K;i++) {
    poly_basemul_montgomery_custom(&t, &a->vec[i], &b->vec[i]);
    poly_add_custom(r, r, &t);
  }

  poly_tomont_custom(r);
}

/*************************************************
* Name:        polyvec_base_mul_acc_tomont_custom_asm
*
* Description: 1. Pointwise multiply elements of a and b, accumulate into r,
*              and multiply by 2^-16.
*              2. Inplace conversion of all coefficients of a polynomial
*              from normal domain to Montgomery domain
*
* Arguments: - poly *r:          pointer to output polynomial
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_base_mul_acc_tomont_custom_asm(poly *r,
                                            const polyvec *a,
                                            const polyvec *b)
{
  size_t avl = KYBER_N;
  size_t vl = ELEMENT_SEW16_PER_VECREG;
  size_t zeta_vl = ELEMENT_SEW16_PER_VECREG/2;

  /*
  ** Firstly, compute polynomial basemul of a->vec[0] and b->vec[0]
  ** a->vec[0] is in v8~v15, b->vec[0] is in v16~v23
  ** result polynomial are in v24~v31
  */
  const int16_t* a_0_ptr = a->vec[0].coeffs;
  const int16_t* b_0_ptr = b->vec[0].coeffs;
  const int16_t* zeta_ptr = zetas_basemul_cg;
  // v8 basemul v16 -> v24
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v8, (%[a_ptr])\n"
    "vle16.v v16, (%[b_ptr])\n"
    "vmod.basemul.vvm v24, v8, v16, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v9 basemul v17 -> v25
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v9, (%[a_ptr])\n"
    "vle16.v v17, (%[b_ptr])\n"
    "vmod.basemul.vvm v25, v9, v17, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v10 basemul v18 -> v26
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v10, (%[a_ptr])\n"
    "vle16.v v18, (%[b_ptr])\n"
    "vmod.basemul.vvm v26, v10, v18, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;    

  // v11 basemul v19 -> v27
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v11, (%[a_ptr])\n"
    "vle16.v v19, (%[b_ptr])\n"
    "vmod.basemul.vvm v27, v11, v19, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;      

  // v12 basemul v20 -> v28
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v12, (%[a_ptr])\n"
    "vle16.v v20, (%[b_ptr])\n"
    "vmod.basemul.vvm v28, v12, v20, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v13 basemul v21 -> v29
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v13, (%[a_ptr])\n"
    "vle16.v v21, (%[b_ptr])\n"
    "vmod.basemul.vvm v29, v13, v21, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v14 basemul v22 -> v30
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v14, (%[a_ptr])\n"
    "vle16.v v22, (%[b_ptr])\n"
    "vmod.basemul.vvm v30, v14, v22, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v15 basemul v23 -> v31
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v15, (%[a_ptr])\n"
    "vle16.v v23, (%[b_ptr])\n"
    "vmod.basemul.vvm v31, v15, v23, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );

  /*
  ** Secondly, compute polynomial basemul of a->vec[i] and b->vec[i], i = 1, 2, ..., KYBER_K-1
  ** and accumulate the result to r
  ** a->vec[i] is in v8~v15, b->vec[i] is in v16~v23
  ** result polynomial are in v8~v15, and then accumulated to v24~v31
  */
  uint32_t i;
  const int16_t* a_i_ptr = NULL;
  const int16_t* b_i_ptr = NULL;
  for(i = 1; i < KYBER_K; i++) {
    a_i_ptr = a->vec[i].coeffs;
    b_i_ptr = b->vec[i].coeffs;
    zeta_ptr = zetas_basemul_cg;

    /*2.1. Basemul Multiplication*/
    // v8 basemul v16 -> v8
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v8, (%[a_ptr])\n"
      "vle16.v v16, (%[b_ptr])\n"
      "vmod.basemul.vvm v8, v8, v16, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v9 basemul v17 -> v9
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v9, (%[a_ptr])\n"
      "vle16.v v17, (%[b_ptr])\n"
      "vmod.basemul.vvm v9, v9, v17, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v10 basemul v18 -> v10
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v10, (%[a_ptr])\n"
      "vle16.v v18, (%[b_ptr])\n"
      "vmod.basemul.vvm v10, v10, v18, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;    

    // v11 basemul v19 -> v27
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v11, (%[a_ptr])\n"
      "vle16.v v19, (%[b_ptr])\n"
      "vmod.basemul.vvm v11, v11, v19, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;      

    // v12 basemul v20 -> v12
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v12, (%[a_ptr])\n"
      "vle16.v v20, (%[b_ptr])\n"
      "vmod.basemul.vvm v12, v12, v20, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v13 basemul v21 -> v13
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v13, (%[a_ptr])\n"
      "vle16.v v21, (%[b_ptr])\n"
      "vmod.basemul.vvm v13, v13, v21, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v14 basemul v22 -> v14
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v14, (%[a_ptr])\n"
      "vle16.v v22, (%[b_ptr])\n"
      "vmod.basemul.vvm v14, v14, v22, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v15 basemul v23 -> v15
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v15, (%[a_ptr])\n"
      "vle16.v v23, (%[b_ptr])\n"
      "vmod.basemul.vvm v15, v15, v23, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );

    /*2.2. Accumulate Process*/
    // add v8~v15 and v24~v31 t0 v24~v31
    __asm__ __volatile__ (
      "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
      "vmod.add.vv v24, v8, v24\n"
      :
      : [avl] "r"(avl)
    );    
  }

  /*
  ** Third, compute poly_tomont
  ** input polynomial is in v24~v31
  ** result polynomial r is in v24~v31
  */
  const int16_t f = (1ULL << 32) % KYBER_Q;
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
    "vmod.mul.vx v24, v24, %[src1]"
    :
    : [avl] "r"(avl), [src1] "r"(f)
  );  

  /*
  ** Finally, store result from vector register v24~v31 to memory
  */ 
  __asm__ __volatile__ ( 
    "vse16.v v24, (%[coeff_base])\n"
    :
    : [coeff_base] "r"(r->coeffs)
  );  
}

/*************************************************
* Name:        polyvec_base_mul_acc_intt_tomont_custom
*
* Description: 1. Pointwise multiply elements of a and b, accumulate into r,
*              and multiply by 2^-16.
*              2. Inverse NTT.
*              3. Inplace conversion of all coefficients of a polynomial
*              from normal domain to Montgomery domain.
*
* Arguments: - poly *r:          pointer to output polynomial
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_base_mul_acc_intt_tomont_custom(poly *r,
                                                const polyvec *a,
                                                const polyvec *b)
{
  polyvec_pointwise_acc_montgomery_custom(r, a, b);
  poly_invntt_tomont_custom(r);
}                                                

/*************************************************
* Name:        polyvec_base_mul_acc_intt_tomont_custom_asm
*
* Description: 1. Pointwise multiply elements of a and b, accumulate into r,
*              and multiply by 2^-16.
*              2. Inverse NTT.
*              3. Inplace conversion of all coefficients of a polynomial
*              from normal domain to Montgomery domain.
*
* Arguments: - poly *r:          pointer to output polynomial
*            - const polyvec *a: pointer to first input vector of polynomials
*            - const polyvec *b: pointer to second input vector of polynomials
**************************************************/
void polyvec_base_mul_acc_intt_tomont_custom_asm(poly *r,
                                                const polyvec *a,
                                                const polyvec *b)
{
  size_t avl = KYBER_N;
  size_t vl = ELEMENT_SEW16_PER_VECREG;
  size_t zeta_vl = ELEMENT_SEW16_PER_VECREG/2;

  /*
  ** Firstly, compute polynomial basemul of a->vec[0] and b->vec[0]
  ** a->vec[0] is in v8~v15, b->vec[0] is in v16~v23
  ** result polynomial are in v24~v31
  */
  const int16_t* a_0_ptr = a->vec[0].coeffs;
  const int16_t* b_0_ptr = b->vec[0].coeffs;
  const int16_t* zeta_ptr = zetas_basemul_cg;
  // v8 basemul v16 -> v24
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v8, (%[a_ptr])\n"
    "vle16.v v16, (%[b_ptr])\n"
    "vmod.basemul.vvm v24, v8, v16, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v9 basemul v17 -> v25
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v9, (%[a_ptr])\n"
    "vle16.v v17, (%[b_ptr])\n"
    "vmod.basemul.vvm v25, v9, v17, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v10 basemul v18 -> v26
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v10, (%[a_ptr])\n"
    "vle16.v v18, (%[b_ptr])\n"
    "vmod.basemul.vvm v26, v10, v18, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;    

  // v11 basemul v19 -> v27
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v11, (%[a_ptr])\n"
    "vle16.v v19, (%[b_ptr])\n"
    "vmod.basemul.vvm v27, v11, v19, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;      

  // v12 basemul v20 -> v28
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v12, (%[a_ptr])\n"
    "vle16.v v20, (%[b_ptr])\n"
    "vmod.basemul.vvm v28, v12, v20, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v13 basemul v21 -> v29
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v13, (%[a_ptr])\n"
    "vle16.v v21, (%[b_ptr])\n"
    "vmod.basemul.vvm v29, v13, v21, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v14 basemul v22 -> v30
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v14, (%[a_ptr])\n"
    "vle16.v v22, (%[b_ptr])\n"
    "vmod.basemul.vvm v30, v14, v22, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );
  zeta_ptr += zeta_vl;
  a_0_ptr += vl;
  b_0_ptr += vl;

  // v15 basemul v23 -> v31
  __asm__ __volatile__ (
    "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
    "vle16.v v0, (%[zeta_base])\n"
    "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
    "vle16.v v15, (%[a_ptr])\n"
    "vle16.v v23, (%[b_ptr])\n"
    "vmod.basemul.vvm v31, v15, v23, v0\n"
    :
    : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_0_ptr), [b_ptr] "r"(b_0_ptr)
  );



  /*
  ** Secondly, compute polynomial basemul of a->vec[i] and b->vec[i], i = 1, 2, ..., KYBER_K-1
  ** and accumulate the result to r
  ** a->vec[i] is in v8~v15, b->vec[i] is in v16~v23
  ** result polynomial are in v8~v15, and then accumulated to v24~v31
  */
  uint32_t i;
  const int16_t* a_i_ptr = NULL;
  const int16_t* b_i_ptr = NULL;
  for(i = 1; i < KYBER_K; i++) {
    a_i_ptr = a->vec[i].coeffs;
    b_i_ptr = b->vec[i].coeffs;
    zeta_ptr = zetas_basemul_cg;

    /*2.1. Basemul Multiplication*/
    // v8 basemul v16 -> v8
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v8, (%[a_ptr])\n"
      "vle16.v v16, (%[b_ptr])\n"
      "vmod.basemul.vvm v8, v8, v16, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v9 basemul v17 -> v9
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v9, (%[a_ptr])\n"
      "vle16.v v17, (%[b_ptr])\n"
      "vmod.basemul.vvm v9, v9, v17, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v10 basemul v18 -> v10
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v10, (%[a_ptr])\n"
      "vle16.v v18, (%[b_ptr])\n"
      "vmod.basemul.vvm v10, v10, v18, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;    

    // v11 basemul v19 -> v27
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v11, (%[a_ptr])\n"
      "vle16.v v19, (%[b_ptr])\n"
      "vmod.basemul.vvm v11, v11, v19, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;      

    // v12 basemul v20 -> v12
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v12, (%[a_ptr])\n"
      "vle16.v v20, (%[b_ptr])\n"
      "vmod.basemul.vvm v12, v12, v20, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v13 basemul v21 -> v13
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v13, (%[a_ptr])\n"
      "vle16.v v21, (%[b_ptr])\n"
      "vmod.basemul.vvm v13, v13, v21, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v14 basemul v22 -> v14
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v14, (%[a_ptr])\n"
      "vle16.v v22, (%[b_ptr])\n"
      "vmod.basemul.vvm v14, v14, v22, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );
    zeta_ptr += zeta_vl;
    a_i_ptr += vl;
    b_i_ptr += vl;

    // v15 basemul v23 -> v15
    __asm__ __volatile__ (
      "vsetvli	zero, %[zeta_vl], e16, m1, ta, mu\n"
      "vle16.v v0, (%[zeta_base])\n"
      "vsetvli	zero, %[vl], e16, m1, ta, mu\n"
      "vle16.v v15, (%[a_ptr])\n"
      "vle16.v v23, (%[b_ptr])\n"
      "vmod.basemul.vvm v15, v15, v23, v0\n"
      :
      : [zeta_base] "r"(zeta_ptr), [zeta_vl] "r"(zeta_vl), [vl] "r"(vl), [a_ptr] "r"(a_i_ptr), [b_ptr] "r"(b_i_ptr)
    );

    /*2.2. Accumulate Process*/
    // add v8~v15 and v24~v31 t0 v24~v31
    __asm__ __volatile__ (
      "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
      "vmod.add.vv v24, v8, v24\n"
      :
      : [avl] "r"(avl)
    );    
  }



  /*
  ** Third, compute INTT
  ** input polynomial is in v24~v31
  ** result polynomial r is in v24~v31
  */  
  avl = KYBER_N;
  size_t load_zeta_vl, bf_vl;
  zeta_ptr = zetas_inv_in_order;

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


  /*
  ** Fourth, compute poly_tomont
  ** input polynomial is in v16~v23
  ** result polynomial r is in v16~v23
  */
  const int16_t f = (1ULL << 32) % KYBER_Q;
  __asm__ __volatile__ ( 
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
    "vmod.mul.vx v16, v16, %[src1]"
    :
    : [avl] "r"(avl), [src1] "r"(f)
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


/*************************************************
* Name:        poly_eta2_add_ntt
*
* Description: 1. Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA2
*              2. Modular add q to each coefficients in a polynomial
*              in this way turn each coefficients to non-negetive number
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length KYBER_SYMBYTES bytes)
*              - uint8_t nonce:       one-byte input nonce
**************************************************/
void poly_eta2_add(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce)
{
  uint8_t buf[KYBER_ETA2*KYBER_N/4];
  kyber_shake256_prf_custom(buf, sizeof(buf), seed, nonce);
  cbd_eta2_custom(r, buf);
	poly_mod_add_q(r);
}

/*************************************************
* Name:        poly_eta2_add_asm
*
* Description: 1. Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA2
*              2. Modular add q to each coefficients in a polynomial
*              in this way turn each coefficients to non-negetive number 
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length KYBER_SYMBYTES bytes)
*              - uint8_t nonce:       one-byte input nonce
**************************************************/
void poly_eta2_add_asm(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce)
{
  uint8_t buf[KYBER_ETA2*KYBER_N/4];
  kyber_shake256_prf_custom(buf, sizeof(buf), seed, nonce);

  /*
	** Firstly, binomial sampling process
  ** Use vector register v8 ~ v15 to contain a polynomial in binomial sampling
	** Use vector register v2 and v4 to contian cbd_even_index and cbd_odd_index
	** Use v16 to contain buf, use v17 to contian unpacked data and vcpopv
	** Use v18 to contain even vrgather, v19 to contain odd vrgather
  */
 	uint8_t* buf_ptr = buf;

#if KYBER_ETA2 != 2
#error "This implementation requires eta2 = 2"
#else
	size_t vl_packed = ELEMENT_SEW8_PER_VECREG >> 2;
	// load v_index0 and v_index1 to v2 and v4
  __asm__ __volatile__ (
    "vsetvli	zero, %[avl], e8, m1, ta, mu\n"
    "vle16.v v2, (%[even_addr])\n"
    "vle16.v v4, (%[odd_addr])\n"
    :
    : [even_addr] "r"(cbd_even_index), [odd_addr] "r"(cbd_odd_index), [avl] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );

	// get v8
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v8, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v9
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v9, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v10
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v10, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v11
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v11, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v12
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v12, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v13
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v13, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v14
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v14, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v15
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v15, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
#endif


  /*
	** Secondly, modular add q process
	** Input coefficients are in v8 ~ v15
	** Result are also in v8 ~ v15
  */
  __asm__ __volatile__ (
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
		"vmod.add.vx v8, v8, %[q]\n"
    :
    : [avl] "r"(KYBER_N), [q] "r"(KYBER_Q)
  ); 	

  /*
  ** Finally, store result from vector register v8~v15 to memory
  */ 
  __asm__ __volatile__ (
    "vse16.v v8, (%[coeff_base])\n"
    :
    : [coeff_base] "r"(r->coeffs)
  );   
}

/*************************************************
* Name:        poly_eta2_addq_add_asm
*
* Description: 1. Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA2
*              2. Modular add q to each coefficients in a polynomial
*              in this way turn each coefficients to non-negetive number
*              3. Modular add with the polynomial r, and store the result to r 
*
* Arguments:   - poly *r:             pointer to output polynomial
*              - const uint8_t *seed: pointer to input seed
*                                     (of length KYBER_SYMBYTES bytes)
*              - uint8_t nonce:       one-byte input nonce
**************************************************/
void poly_eta2_addq_add_asm(poly *r, const uint8_t seed[KYBER_SYMBYTES], uint8_t nonce)
{
  uint8_t buf[KYBER_ETA2*KYBER_N/4];
  kyber_shake256_prf_custom(buf, sizeof(buf), seed, nonce);

  /*
	** Firstly, binomial sampling process
  ** Use vector register v8 ~ v15 to contain a polynomial in binomial sampling
	** Use vector register v2 and v4 to contian cbd_even_index and cbd_odd_index
	** Use v16 to contain buf, use v17 to contian unpacked data and vcpopv
	** Use v18 to contain even vrgather, v19 to contain odd vrgather
  */
 	uint8_t* buf_ptr = buf;

#if KYBER_ETA2 != 2
#error "This implementation requires eta2 = 2"
#else
	size_t vl_packed = ELEMENT_SEW8_PER_VECREG >> 2;
	// load v_index0 and v_index1 to v2 and v4
  __asm__ __volatile__ (
    "vsetvli	zero, %[avl], e8, m1, ta, mu\n"
    "vle16.v v2, (%[even_addr])\n"
    "vle16.v v4, (%[odd_addr])\n"
    :
    : [even_addr] "r"(cbd_even_index), [odd_addr] "r"(cbd_odd_index), [avl] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );

	// get v8
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v8, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v9
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v9, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v10
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v10, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v11
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v11, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v12
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v12, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v13
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v13, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v14
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v14, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
	buf_ptr += vl_packed;

	// get v15
  __asm__ __volatile__ (
    "vsetvli	zero, %[vl_packed], e8, m1, ta, mu\n"
    "vle16.v v16, (%[buf_ptr])\n"
    "vsetvli	zero, %[avl0], e8, m1, ta, mu\n"
    "unpack.vx v17, v16, %[valid_num]\n"
		"vcpop.v v17, v17\n"
		"vsetvli	zero, %[avl1], e8, m1, ta, mu\n"
		"vrgather.vv v18, v17, v2\n"
		"vrgather.vv v19, v17, v4\n"
		"vsetvli	zero, %[avl1], e8, mf2, ta, mu\n"
		"vsub.vv v18, v18, v19\n"
		"vwadd.vx v15, v18, x0\n"
    :
    : [buf_ptr] "r"(buf_ptr), [vl_packed] "r"(vl_packed), [avl0] "r"(ELEMENT_SEW8_PER_VECREG), [valid_num] "r"(2), [avl1] "r"(ELEMENT_SEW8_PER_VECREG/2)
  );
#endif


  /*
	** Secondly, modular add q process
	** Input coefficients are in v8 ~ v15
	** Result are also in v8 ~ v15
  */
  __asm__ __volatile__ (
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
		"vmod.add.vx v8, v8, %[q]\n"
    :
    : [avl] "r"(KYBER_N), [q] "r"(KYBER_Q)
  ); 	

  /*
	** Third, modular add with poly r
	** Input coefficients are in v8 ~ v15
  ** Input coefficients of poly r are loaded to v16 ~ v23
	** Result are also in v16 ~ v23
  */
  __asm__ __volatile__ (
    "vsetvli	zero, %[avl], e16, m8, ta, mu\n"
    "vle16.v v16, (%[coeff_base])\n"
		"vmod.add.vv v16, v8, v16\n"
    :
    : [avl] "r"(KYBER_N), [coeff_base] "r"(r->coeffs)
  );

  /*
  ** Finally, store result from vector register v16~v23 to memory
  */ 
  __asm__ __volatile__ (
    "vse16.v v16, (%[coeff_base])\n"
    :
    : [coeff_base] "r"(r->coeffs)
  );   
}

/*************************************************
* Name:        polyvec_ntt_custom_asm
*
* Description: Apply forward NTT to all elements of a vector of polynomials
*
* Arguments:   - polyvec *r: pointer to in/output vector of polynomials
**************************************************/
void polyvec_ntt_custom_asm(polyvec *r)
{
  unsigned int i;
  for(i=0;i<KYBER_K;i++)
    ntt_cg_custom_reorder_asm(r->vec[i].coeffs);
}