#include <stdint.h>
#include "params.h"
#include "cbd.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

#if (VLEN == 256)
const uint8_t cbd_even_index[ELEMENT_SEW8_PER_VECREG/2] = {
  0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30
};
const uint8_t cbd_odd_index[ELEMENT_SEW8_PER_VECREG/2] = {
  1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31
};
#elif (VLEN == 512)
const uint8_t cbd_even_index[ELEMENT_SEW8_PER_VECREG/2] = {
  0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
  32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62
};
const uint8_t cbd_odd_index[ELEMENT_SEW8_PER_VECREG/2] = {
  1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31,
  33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63
};
#elif (VLEN == 1024)
const uint8_t cbd_even_index[ELEMENT_SEW8_PER_VECREG/2] = {
  0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30
};
const uint8_t cbd_odd_index[ELEMENT_SEW8_PER_VECREG/2] = {
  1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31
};
# else
#error "VLEN must be 256/512/1024"
#endif 

/*************************************************
* Name:        cbd2_custom
*
* Description: Given an array of uniformly random bytes, compute
*              polynomial with coefficients distributed according to
*              a centered binomial distribution with parameter eta=2
*
* Arguments:   - poly *r:            pointer to output polynomial
*              - const uint8_t *buf: pointer to input byte array
**************************************************/
static void cbd2_custom(poly *r, const uint8_t buf[2*KYBER_N/4])
{
  int16_t *coeffs_ptr = r->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  vuint8m1_t vt;
  vint8m1_t vx;
  vint8m1_t v_op0, v_op1;
  vint16m2_t vy;
  vint16m1_t vz;

  vuint8m1_t v_index0, v_index1;
  v_index0 = vle8_v_u8m1(cbd_even_index, ELEMENT_SEW8_PER_VECREG/2);
  v_index1 = vle8_v_u8m1(cbd_odd_index, ELEMENT_SEW8_PER_VECREG/2);

  while(avl > 0) {
    // in this case, ELEMENT_SEW8_PER_VECREG >> 2 bytes packed data precisely 
    // unpacked to a whole vector register file with ELEMENT_SEW8_PER_VECREG elements
    // each element has 2 valid bits
    size_t vl_packed = ELEMENT_SEW8_PER_VECREG >> 2;
    vt = vle8_v_u8m1(buf, vl_packed);
    buf += vl_packed;
    //should use asm volatile otherwise will be optimized out and triger panic if
    vl = vsetvl_e8m1_wrapper(ELEMENT_SEW8_PER_VECREG);
    vx = unpack_vx_i8m1(vt, 2);
    vx = vcpop_v_i8m1(vx);

    v_op0 = vrgather_vv_i8m1(vx, v_index0, ELEMENT_SEW8_PER_VECREG/2);
    v_op1 = vrgather_vv_i8m1(vx, v_index1, ELEMENT_SEW8_PER_VECREG/2);
    v_op0 = vsub_vv_i8m1(v_op0, v_op1, ELEMENT_SEW8_PER_VECREG/2);
    vy = vsext_vf2_i16m2(v_op0, ELEMENT_SEW8_PER_VECREG/2);
    vz = vget_v_i16m2_i16m1(vy, 0);
    vse16_v_i16m1(coeffs_ptr, vz, ELEMENT_SEW8_PER_VECREG/2);
    coeffs_ptr += ELEMENT_SEW8_PER_VECREG/2;
    avl -= ELEMENT_SEW8_PER_VECREG/2;
  }
}

/*************************************************
* Name:        cbd3_custom
*
* Description: Given an array of uniformly random bytes, compute
*              polynomial with coefficients distributed according to
*              a centered binomial distribution with parameter eta=3
*              This function is only needed for Kyber-512
*
* Arguments:   - poly *r:            pointer to output polynomial
*              - const uint8_t *buf: pointer to input byte array
**************************************************/
static void cbd3_custom(poly *r, const uint8_t buf[3*KYBER_N/4])
{
  int16_t *coeffs_ptr = r->coeffs;
  size_t vl;
  size_t avl = KYBER_N;
  vuint8m1_t vt;
  vint8m1_t vx;
  vint8m1_t v_op0, v_op1;
  vint16m2_t vy;
  vint16m1_t vz;

  vuint8m1_t v_index0, v_index1;
  v_index0 = vle8_v_u8m1(cbd_even_index, ELEMENT_SEW8_PER_VECREG/2);
  v_index1 = vle8_v_u8m1(cbd_odd_index, ELEMENT_SEW8_PER_VECREG/2);
  while( avl > 0) {
    // in this case, (ELEMENT_SEW8_PER_VECREG * 3) >> 3 bytes packed data precisely 
    // unpacked to a whole vector register file with ELEMENT_SEW8_PER_VECREG elements
    // each element has 3 valid bits
    size_t vl_packed = (ELEMENT_SEW8_PER_VECREG * 3) >> 3;
    vt = vle8_v_u8m1(buf, vl_packed);
    buf += vl_packed;
    //should use asm volatile otherwise will be optimized out and tirger panic if
    vl = vsetvl_e8m1_wrapper(ELEMENT_SEW8_PER_VECREG);
    vx = unpack_vx_i8m1(vt,3);
    vx = vcpop_v_i8m1(vx);

    v_op0 = vrgather_vv_i8m1(vx, v_index0, ELEMENT_SEW8_PER_VECREG/2);
    v_op1 = vrgather_vv_i8m1(vx, v_index1, ELEMENT_SEW8_PER_VECREG/2);
    v_op0 = vsub_vv_i8m1(v_op0, v_op1, ELEMENT_SEW8_PER_VECREG/2);
    vy = vsext_vf2_i16m2(v_op0, ELEMENT_SEW8_PER_VECREG/2);
    vz = vget_v_i16m2_i16m1(vy, 0);
    vse16_v_i16m1(coeffs_ptr, vz, ELEMENT_SEW8_PER_VECREG/2);
    coeffs_ptr += ELEMENT_SEW8_PER_VECREG/2;
    avl -= ELEMENT_SEW8_PER_VECREG/2;
  }
}

void cbd_eta1_custom(poly *r, const uint8_t buf[KYBER_ETA1*KYBER_N/4])
{
#if KYBER_ETA1 == 2
  cbd2_custom(r, buf);
#elif KYBER_ETA1 == 3
  cbd3_custom(r, buf);
#else
#error "This implementation requires eta1 in {2,3}"
#endif
}

void cbd_eta2_custom(poly *r, const uint8_t buf[KYBER_ETA1*KYBER_N/4])
{
#if KYBER_ETA2 != 2
#error "This implementation requires eta2 = 2"
#else
  cbd2_custom(r, buf);
#endif
}