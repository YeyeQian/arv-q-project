#include "inner.h"
#include <riscv_vector.h>
#include <assert.h>
#include "../apis/custom_inst_api.h"
#include "params_custom.h"
#include "../apis/vsetvl_wrapper.h"

/*************************************************
* Name:        rej_uniform_custom
*
* Description: Sample uniformly random coefficients in [0, Q-1] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - uint16_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
#define UNPACK_REJ_LOAD_BYTE	56	// 14*32/8 = 56 bytes
#define UNPACKED_E16_VL 32				// vsetvlmax_e16m1
static unsigned int rej_uniform_custom(uint16_t *a,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  size_t ctr, pos, vl, valid_num;
	vuint8m1_t vreg_m;
  vuint16m1_t vreg_rej;

  csr_validnum_rw();
  ctr = pos = 0;

  while((ctr < len) && (pos < buflen)) {
    if((pos + UNPACK_REJ_LOAD_BYTE) < buflen) {
      vl = vsetvl_e8m1(UNPACK_REJ_LOAD_BYTE);
    }
    else {
      vl = vsetvl_e8m1(buflen - pos);
    }

    vreg_m = vle8_v_u8m1(buf + pos, vl);
    pos += vl;
		
		vsetvl_e16m1_wrapper(UNPACKED_E16_VL);

    vreg_rej = unpack_vx_u16m1(vreg_m, 14);

    vreg_rej = sample_rej_vx_u16m1(vreg_rej, Q);
    valid_num = csr_validnum_rw();
		
    if(ctr + valid_num > len) {
      valid_num = len - ctr;
    }

    vl = vsetvl_e16m1(valid_num);
    vse16_v_u16m1(a + ctr, vreg_rej, valid_num);
    ctr += valid_num;
  }

  return ctr;
}

/* see inner.h */
#if FALCON_N == 512
#define POLY_UNIFORM_NBLOCKS_SHAKE256 ((1195 + SHAKE256_RATE - 1) / SHAKE256_RATE)	// equal to 9
#elif FALCON_N == 1024
#define POLY_UNIFORM_NBLOCKS_SHAKE256 ((2390 + SHAKE256_RATE - 1) / SHAKE256_RATE)
#endif
void
Zf(hash_to_point_vartime_custom)(uint16_t *x, unsigned logn)
{
  unsigned int ctr;
  unsigned int buflen = POLY_UNIFORM_NBLOCKS_SHAKE256 * SHAKE256_RATE;
  uint8_t buf[POLY_UNIFORM_NBLOCKS_SHAKE256 * SHAKE256_RATE];
	uint32_t n = 1 << logn;

	keccak_squeezeblocks_custom(buf, POLY_UNIFORM_NBLOCKS_SHAKE256, SHAKE256_RATE, false);

	ctr = rej_uniform_custom(x, n, buf, buflen);

  while(ctr < n){
    keccak_squeezeblocks_custom(buf, 1, SHAKE256_RATE, true);
    buflen = SHAKE256_RATE;
    ctr += rej_uniform_custom(x + ctr, n - ctr, buf, buflen);
  }	
}


/* see inner.h */
void
Zf(hash_to_point_ct_custom)(
	size_t *sc, uint16_t *x, unsigned logn, uint8_t *tmp)
{
	/*
	 * Each 16-bit sample is a value in 0..65535. The value is
	 * kept if it falls in 0..61444 (because 61445 = 5*12289)
	 * and rejected otherwise; thus, each sample has probability
	 * about 0.93758 of being selected.
	 *
	 * We want to oversample enough to be sure that we will
	 * have enough values with probability at least 1 - 2^(-256).
	 * Depending on degree N, this leads to the following
	 * required oversampling:
	 *
	 *   logn     n  oversampling
	 *     1      2     65
	 *     2      4     67
	 *     3      8     71
	 *     4     16     77
	 *     5     32     86
	 *     6     64    100
	 *     7    128    122
	 *     8    256    154
	 *     9    512    205
	 *    10   1024    287
	 *
	 * If logn >= 7, then the provided temporary buffer is large
	 * enough. Otherwise, we use a stack buffer of 63 entries
	 * (i.e. 126 bytes) for the values that do not fit in tmp[].
	 */

	static const uint16_t overtab[] = {
		0, /* unused */
		65,
		67,
		71,
		77,
		86,
		100,
		122,
		154,
		205,
		287
	};

	size_t vlmax = vsetvlmax_e16m1();
	uint8_t buf[vlmax << 1];
	uint16_t *tt = (uint16_t *)buf;
	uint16_t n = (size_t)1 << logn;
	uint16_t m = n + overtab[logn];
	uint16_t *t = (uint16_t *)tmp;
	uint16_t fiveq = 61445;
    for (uint16_t vlrej; m > 0; t += vlrej) {
		// inner_shake256_extract_custom(sc, (void *)buf, sizeof(buf));
        size_t vl = vsetvl_e16m1(m);
        vuint16m1_t vec_value = vle16_v_u16m1(tt, vl);
        vuint16m1_t vec_out;
        vec_out = sample_rej_vx_u16m1(vec_value, fiveq);
        vlrej = csr_validnum_rw();
        vse16_v_u16m1(t, vec_out, vlrej);
		m -= vl;
    }
	t = (uint16_t *)tmp;
    for (size_t vl; n > 0; n -= vl, t += vl, x += vl) {
        vl = vsetvl_e16m1(n);
        vuint16m1_t vec_value = vle16_v_u16m1(t, vl);
        vec_value = vmod_mul_vx_u16m1(vec_value, R);
        vse16_v_u16m1(x, vec_value, vl);
	}
}


/*
 * Acceptance bound for the (squared) l2-norm of the signature depends
 * on the degree. This array is indexed by logn (1 to 10). These bounds
 * are _inclusive_ (they are equal to floor(beta^2)).
 */
static const uint32_t l2bound[] = {
	0,    /* unused */
	101498,
	208714,
	428865,
	892039,
	1852696,
	3842630,
	7959734,
	16468416,
	34034726,
	70265242
};


/* see inner.h */
int
Zf(is_short_custom)(
	const int16_t *s1, const int16_t *s2, unsigned logn)
{
	/*
	 * We use the l2-norm. Code below uses only 32-bit operations to
	 * compute the square of the norm with saturation to 2^32-1 if
	 * the value exceeds 2^31-1.
	 */

	size_t n = (size_t)1 << logn;
	vint64m1_t scalar = vmv_v_x_i64m1(0, 1);
    for (size_t vl; n > 0; n -= vl, s1 += vl, s2 += vl) {
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_s1 = vle16_v_i16m1(s1, vl);
        vint16m1_t vec_s2 = vle16_v_i16m1(s2, vl);
		vint32m2_t vec_s1mul = vwmul_vv_i32m2(vec_s1, vec_s1, vl);
		vint32m2_t vec_s2mul = vwmul_vv_i32m2(vec_s2, vec_s2, vl);
		scalar = vwredsum_vs_i32m2_i64m1(scalar, vec_s1mul, scalar, vl);
		scalar = vwredsum_vs_i32m2_i64m1(scalar, vec_s2mul, scalar, vl);
    }
	uint64_t s = vmv_x_s_i64m1_i64(scalar);

	return s <= l2bound[logn];
}


/* see inner.h */
int
Zf(is_short_half_custom)(
	uint64_t sqn, const int16_t *s2, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	vint64m1_t scalar = vmv_v_x_i64m1(sqn, 1);
    for (size_t vl; n > 0; n -= vl, s2 += vl) {
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_s2 = vle16_v_i16m1(s2, vl);
		vint32m2_t vec_s2mul = vwmul_vv_i32m2(vec_s2, vec_s2, vl);
		scalar = vwredsum_vs_i32m2_i64m1(scalar, vec_s2mul, scalar, vl);
    }
	uint64_t s = vmv_x_s_i64m1_i64(scalar);

	return s <= l2bound[logn];
}


// void
// Zf(hash_to_point_vartime_custom)(uint16_t *x, unsigned logn)
// {
// 	/*
// 	 * This is the straightforward per-the-spec implementation. It
// 	 * is not constant-time, thus it might reveal information on the
// 	 * plaintext (at least, enough to check the plaintext against a
// 	 * list of potential plaintexts) in a scenario where the
// 	 * attacker does not have access to the signature value or to
// 	 * the public key, but knows the nonce (without knowledge of the
// 	 * nonce, the hashed output cannot be matched against potential
// 	 * plaintexts).
// 	 */
// 	size_t n = (size_t)1 << logn;
// 	size_t vlmax = vsetvlmax_e16m1();
// 	uint8_t buf[vlmax << 1];
// 	uint16_t *tt = (uint16_t *)buf;
// 	uint16_t fiveq = 61445;
//     for (uint16_t vlx; n > 0; n -= vlx, x += vlx) {
// 		inner_shake256_extract_custom(sc, (void *)buf, sizeof(buf));
//         size_t vl = vsetvl_e16m1(vlmax);
//         vuint16m1_t vec_value = vle16_v_u16m1(tt, vl);
//         vuint16m1_t vec_out;
// 		uint16_t vlrej;
//         vec_out = sample_rej_vx_u16m1(vec_value, fiveq);
//         vlrej = csr_validnum_rw();
//         vl = vsetvl_e16m1(vlrej);
//         vec_out = vmod_mul_vx_u16m1(vec_out, R);
// 		vlx = vlrej > n ? n : vlrej;
//         vse16_v_u16m1(x, vec_out, vlx);
// 	}
// }