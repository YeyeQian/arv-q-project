#include "inner.h"
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"
#include "params_custom.h"

#if (VLEN == 256)
#define VLMAXFFT 4
#elif (VLEN == 512)
#define VLMAXFFT 8
#elif (VLEN == 1024)
#define VLMAXFFT 16
# else
#error "VLEN must be 256/512/1024"
#endif

/*
 * Rules for complex number macros:
 * --------------------------------
 *
 * Operand order is: destination, source1, source2...
 *
 * Each operand is a real and an imaginary part.
 *
 * All overlaps are allowed.
 */

/*
 * Addition of two complex numbers (d = a + b).
 */
#define FPC_ADD(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_re, fpct_im; \
		fpct_re = fpr_add(a_re, b_re); \
		fpct_im = fpr_add(a_im, b_im); \
		(d_re) = fpct_re; \
		(d_im) = fpct_im; \
	} while (0)

/*
 * Subtraction of two complex numbers (d = a - b).
 */
#define FPC_SUB(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_re, fpct_im; \
		fpct_re = fpr_sub(a_re, b_re); \
		fpct_im = fpr_sub(a_im, b_im); \
		(d_re) = fpct_re; \
		(d_im) = fpct_im; \
	} while (0)

/*
 * Multplication of two complex numbers (d = a * b).
 */
#define FPC_MUL(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_a_re, fpct_a_im; \
		fpr fpct_b_re, fpct_b_im; \
		fpr fpct_d_re, fpct_d_im; \
		fpct_a_re = (a_re); \
		fpct_a_im = (a_im); \
		fpct_b_re = (b_re); \
		fpct_b_im = (b_im); \
		fpct_d_re = fpr_sub( \
			fpr_mul(fpct_a_re, fpct_b_re), \
			fpr_mul(fpct_a_im, fpct_b_im)); \
		fpct_d_im = fpr_add( \
			fpr_mul(fpct_a_re, fpct_b_im), \
			fpr_mul(fpct_a_im, fpct_b_re)); \
		(d_re) = fpct_d_re; \
		(d_im) = fpct_d_im; \
	} while (0)

/*
 * Division of complex numbers (d = a / b).
 */
#define FPC_DIV(d_re, d_im, a_re, a_im, b_re, b_im)   do { \
		fpr fpct_a_re, fpct_a_im; \
		fpr fpct_b_re, fpct_b_im; \
		fpr fpct_d_re, fpct_d_im; \
		fpr fpct_m; \
		fpct_a_re = (a_re); \
		fpct_a_im = (a_im); \
		fpct_b_re = (b_re); \
		fpct_b_im = (b_im); \
		fpct_m = fpr_add(fpr_sqr(fpct_b_re), fpr_sqr(fpct_b_im)); \
		fpct_m = fpr_inv(fpct_m); \
		fpct_b_re = fpr_mul(fpct_b_re, fpct_m); \
		fpct_b_im = fpr_mul(fpr_neg(fpct_b_im), fpct_m); \
		fpct_d_re = fpr_sub( \
			fpr_mul(fpct_a_re, fpct_b_re), \
			fpr_mul(fpct_a_im, fpct_b_im)); \
		fpct_d_im = fpr_add( \
			fpr_mul(fpct_a_re, fpct_b_im), \
			fpr_mul(fpct_a_im, fpct_b_re)); \
		(d_re) = fpct_d_re; \
		(d_im) = fpct_d_im; \
	} while (0)


/* see inner.h */
void
Zf(FFT_custom)(fpr *f, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t m;
	unsigned u;
	size_t t = hn;

	uint16_t vlmax = vsetvlmax_e64m1();
	uint16_t log2vlmax;
	for (log2vlmax = 1; log2vlmax <= 6; log2vlmax ++)
		if ((vlmax >> log2vlmax) & 1)
			break;
	
	for (u = 1, m = 2; u < logn; u ++, m <<= 1) {
		size_t ht, hm, i1, j1;
		ht = t >> 1;
		hm = m >> 1;
		for (i1 = 0, j1 = 0; i1 < hm; i1 ++, j1 += t) {
			size_t j, j2;
			j2 = j1 + ht;
			fpr s_re, s_im;
			s_re = fpr_gm_tab_reduce[m + (i1 << 1) + 0];
			s_im = fpr_gm_tab_reduce[m + (i1 << 1) + 1];
			if(ht >= VLMAXFFT) {
				double *pc = (double *)(f + j1);
				for(j = 0; j < (ht >> log2vlmax); j++) {
					vfloat64m1_t vec_hi_re = vle64_v_f64m1(pc + ht, vlmax);
					vfloat64m1_t vec_hi_im = vle64_v_f64m1(pc + ht + hn, vlmax);
					vfloat64m1_t vec_lo_re = vle64_v_f64m1(pc, vlmax);
					vfloat64m1_t vec_lo_im = vle64_v_f64m1(pc + hn, vlmax);
					vfloat64m1_t vec_mul_re = vfmul_vf_f64m1(vec_hi_re, s_re.v, vlmax);
					vfloat64m1_t vec_mul_im = vfmul_vf_f64m1(vec_hi_re, s_im.v, vlmax);
					vec_mul_re = vfnmsac_vf_f64m1 (vec_mul_re, s_im.v, vec_hi_im, vlmax);
					vec_mul_im = vfmacc_vf_f64m1 (vec_mul_im, s_re.v, vec_hi_im, vlmax);
					vfloat64m1_t vec_add_re = vfadd_vv_f64m1(vec_lo_re, vec_mul_re, vlmax);
					vfloat64m1_t vec_sub_re = vfsub_vv_f64m1(vec_lo_re, vec_mul_re, vlmax);
					vfloat64m1_t vec_add_im = vfadd_vv_f64m1(vec_lo_im, vec_mul_im, vlmax);
					vfloat64m1_t vec_sub_im = vfsub_vv_f64m1(vec_lo_im, vec_mul_im, vlmax);
					vse64_v_f64m1(pc, vec_add_re, vlmax);
					vse64_v_f64m1(pc + ht, vec_sub_re, vlmax);
					vse64_v_f64m1(pc + hn, vec_add_im, vlmax);
					vse64_v_f64m1(pc + ht + hn, vec_sub_im, vlmax);
					pc += vlmax;
				}
			}
			else {
				for (j = j1; j < j2; j ++) {
					fpr x_re, x_im, y_re, y_im;
					x_re = f[j];
					x_im = f[j + hn];
					y_re = f[j + ht];
					y_im = f[j + ht + hn];
					FPC_MUL(y_re, y_im, y_re, y_im, s_re, s_im);
					FPC_ADD(f[j], f[j + hn],
						x_re, x_im, y_re, y_im);
					FPC_SUB(f[j + ht], f[j + ht + hn],
						x_re, x_im, y_re, y_im);
				}
			}
		}
		t = ht;
	}
}

/* see inner.h */
void
Zf(iFFT_custom)(fpr *f, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t t = 1;
	size_t m = n;

	uint16_t vlmax = vsetvlmax_e64m1();
	uint16_t log2vlmax;
	for (log2vlmax = 1; log2vlmax <= 6; log2vlmax ++)
		if ((vlmax >> log2vlmax) & 1)
			break;
	
	for (size_t u = logn; u > 1; u --) {
		size_t hm, dt, i1, j1;
		hm = m >> 1;
		dt = t << 1;
		for (i1 = 0, j1 = 0; j1 < hn; i1 ++, j1 += dt) {
			size_t j, j2;
			j2 = j1 + t;
			fpr s_re, s_im;
			s_re = fpr_gm_tab_reduce[hm + (i1 << 1) + 0];
			s_im = fpr_gm_tab_reduce[hm + (i1 << 1) + 1];
			if (t >= VLMAXFFT) {
				double *pc = (double *)(f + j1);
				for(uint16_t j = 0; j < (t >> log2vlmax); j++) {
					vfloat64m1_t vec_lo_re = vle64_v_f64m1(pc, vlmax);
					vfloat64m1_t vec_hi_re = vle64_v_f64m1(pc + t, vlmax);
					vfloat64m1_t vec_lo_im = vle64_v_f64m1(pc + hn, vlmax);
					vfloat64m1_t vec_hi_im = vle64_v_f64m1(pc + t + hn, vlmax);
					vfloat64m1_t vec_sub_re = vfsub_vv_f64m1(vec_lo_re, vec_hi_re, vlmax);
					vfloat64m1_t vec_add_re = vfadd_vv_f64m1(vec_lo_re, vec_hi_re, vlmax);
					vfloat64m1_t vec_sub_im = vfsub_vv_f64m1(vec_lo_im, vec_hi_im, vlmax);
					vfloat64m1_t vec_add_im = vfadd_vv_f64m1(vec_lo_im, vec_hi_im, vlmax);
					vfloat64m1_t vec_mul_re = vfmul_vf_f64m1(vec_sub_re, s_re.v, vlmax);
					vfloat64m1_t vec_mul_im = vfmul_vf_f64m1(vec_sub_re, s_im.v, vlmax);
					vec_mul_re = vfmacc_vf_f64m1 (vec_mul_re, s_im.v, vec_sub_im, vlmax);
					vec_mul_im = vfmsac_vf_f64m1 (vec_mul_im, s_re.v, vec_sub_im, vlmax);
					vse64_v_f64m1(pc, vec_add_re, vlmax);
					vse64_v_f64m1(pc + hn, vec_add_im, vlmax);
					vse64_v_f64m1(pc + t, vec_mul_re, vlmax);
					vse64_v_f64m1(pc + t + hn, vec_mul_im, vlmax);
					pc += vlmax;
				}
			}
			else {
				s_im = fpr_neg(s_im);
				for (j = j1; j < j2; j ++) {
					fpr x_re, x_im, y_re, y_im;
					x_re = f[j];
					x_im = f[j + hn];
					y_re = f[j + t];
					y_im = f[j + t + hn];
					FPC_ADD(f[j], f[j + hn],
						x_re, x_im, y_re, y_im);
					FPC_SUB(x_re, x_im, x_re, x_im, y_re, y_im);
					FPC_MUL(f[j + t], f[j + t + hn],
						x_re, x_im, s_re, s_im);
				}
			}
		}
		t = dt;
		m = hm;
	}
	/*
	* Last iteration is a no-op, provided that we divide by N/2
	* instead of N. We need to make a special case for logn = 0.
	*/
	if (logn > 0) {
		fpr ni;
		ni = fpr_p2_tab[logn];
		if(n >= vlmax) {
			vfloat64m1_t vec_ni = vfmv_v_f_f64m1(ni.v, vlmax);
			double *a = (double *)f;
			for (size_t vl; n > 0; n -= vl, a += vl) {
				vl = vsetvl_e64m1(n);
				vfloat64m1_t vec_a = vle64_v_f64m1(a, vl);
				vfloat64m1_t vec_mul = vfmul_vv_f64m1(vec_a, vec_ni, vl);
				vse64_v_f64m1(a, vec_mul, vl);
			}
		}
		else {
			for (size_t u = 0; u < n; u ++) {
				f[u] = fpr_mul(f[u], ni);
			}
		}
	}
}


/* see inner.h */
void
Zf(poly_add_custom)(
	fpr *restrict a, const fpr *restrict b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	if (n >= VLMAXFFT) {
		for (size_t vl; n > 0; n -= vl, a += vl, b += vl) {
			vl = vsetvl_e64m1(n);
			vfloat64m1_t vec_a = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_b = vle64_v_f64m1((double *)b, vl);
			vfloat64m1_t vec_add = vfadd_vv_f64m1(vec_a, vec_b, vl);
			vse64_v_f64m1((double *)a, vec_add, vl);
		}
		// n = (size_t)1 << logn;
		// a -= n;
	}
	else {
		for (size_t u = 0; u < n; u ++) {
			a[u] = fpr_add(a[u], b[u]);
		}
	}
}


/* see inner.h */
void
Zf(poly_sub_custom)(
	fpr *restrict a, const fpr *restrict b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	if (n >= VLMAXFFT) {
		for (size_t vl; n > 0; n -= vl, a += vl, b += vl) {
			vl = vsetvl_e64m1(n);
			vfloat64m1_t vec_a = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_b = vle64_v_f64m1((double *)b, vl);
			vfloat64m1_t vec_sub = vfsub_vv_f64m1(vec_a, vec_b, vl);
			vse64_v_f64m1((double *)a, vec_sub, vl);
		}
	}
	else {
		for (size_t u = 0; u < n; u ++) {
			a[u] = fpr_sub(a[u], b[u]);
		}
	}
}


/* see inner.h */
void
Zf(poly_neg_custom)(fpr *a, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	if (n >= VLMAXFFT) {
		for (size_t vl; n > 0; n -= vl, a += vl) {
			vl = vsetvl_e64m1(n);
			vfloat64m1_t vec_a = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_sub = vfrsub_vf_f64m1(vec_a, 0, vl);
			vse64_v_f64m1((double *)a, vec_sub, vl);
		}
	}
	else {
		for (size_t u = 0; u < n; u ++) {
			a[u] = fpr_neg(a[u]);
		}
	}
}


/* see inner.h */
void
Zf(poly_adj_fft_custom)(fpr *a, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		double *pa = (double *)a + hn;
		for (size_t vl; hn > 0; hn -= vl, pa += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_a = vle64_v_f64m1(pa, vl);
			vfloat64m1_t vec_sub = vfrsub_vf_f64m1(vec_a, 0, vl);
			vse64_v_f64m1(pa, vec_sub, vl);
		}
	}
	else {
		for (size_t u = hn; u < n; u ++) {
			a[u] = fpr_neg(a[u]);
		}
	}
}


/* see inner.h */
void
Zf(poly_mul_fft_custom)(
	fpr *restrict a, const fpr *restrict b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		for (size_t vl; hn > 0; hn -= vl, a += vl, b += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_b_re = vle64_v_f64m1((double *)b, vl);
			vfloat64m1_t vec_a_re = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_b_im = vle64_v_f64m1((double *)b + hN, vl);
			vfloat64m1_t vec_a_im = vle64_v_f64m1((double *)a + hN, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_a_re, vec_b_re, vl);
			vfloat64m1_t vec_mul_im = vfmul_vv_f64m1(vec_a_re, vec_b_im, vl);
			vec_mul_re = vfnmsac_vv_f64m1(vec_mul_re, vec_a_im, vec_b_im, vl);
			vec_mul_im = vfmacc_vv_f64m1(vec_mul_im, vec_a_im, vec_b_re, vl);
			vse64_v_f64m1((double *)a, vec_mul_re, vl);
			vse64_v_f64m1((double *)a + hN, vec_mul_im, vl);
		}
		// a -= hN;
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			fpr a_re, a_im, b_re, b_im;

			a_re = a[u];
			a_im = a[u + hn];
			b_re = b[u];
			b_im = b[u + hn];
			FPC_MUL(a[u], a[u + hn], a_re, a_im, b_re, b_im);
		}
	}
}


/* see inner.h */
void
Zf(poly_muladj_fft_custom)(
	fpr *restrict a, const fpr *restrict b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		for (size_t vl; hn > 0; hn -= vl, a += vl, b += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_b_re = vle64_v_f64m1((double *)b, vl);
			vfloat64m1_t vec_a_re = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_b_im = vle64_v_f64m1((double *)b + hN, vl);
			vfloat64m1_t vec_a_im = vle64_v_f64m1((double *)a + hN, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_a_re, vec_b_re, vl);
			vfloat64m1_t vec_mul_im = vfmul_vv_f64m1(vec_a_re, vec_b_im, vl);
			vec_mul_re = vfmacc_vv_f64m1(vec_mul_re, vec_a_im, vec_b_im, vl);
			vec_mul_im = vfmsac_vv_f64m1(vec_mul_im, vec_a_im, vec_b_re, vl);
			vse64_v_f64m1((double *)a, vec_mul_re, vl);
			vse64_v_f64m1((double *)a + hN, vec_mul_im, vl);
		}
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			fpr a_re, a_im, b_re, b_im;

			a_re = a[u];
			a_im = a[u + hn];
			b_re = b[u];
			b_im = fpr_neg(b[u + hn]);
			FPC_MUL(a[u], a[u + hn], a_re, a_im, b_re, b_im);
		}
	}
}


/* see inner.h */
void
Zf(poly_mulselfadj_fft_custom)(fpr *a, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		vfloat64m1_t vec_zero = vfmv_v_f_f64m1(0, VLMAXFFT);
		for (size_t vl; hn > 0; hn -= vl, a += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_a_re = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_a_im = vle64_v_f64m1((double *)a + hN, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_a_re, vec_a_re, vl);
			vec_mul_re = vfmacc_vv_f64m1(vec_mul_re, vec_a_im, vec_a_im, vl);
			vse64_v_f64m1((double *)a, vec_mul_re, vl);
			vse64_v_f64m1((double *)a + hN, vec_zero, vl);
		}
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			fpr a_re, a_im;

			a_re = a[u];
			a_im = a[u + hn];
			a[u] = fpr_add(fpr_sqr(a_re), fpr_sqr(a_im));
			a[u + hn] = fpr_zero;
		}
	}
}


/* see inner.h */
void
Zf(poly_mulconst_custom)(fpr *a, fpr x, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	if (n >= VLMAXFFT) {
		for (size_t vl; n > 0; n -= vl, a += vl) {
			vl = vsetvl_e64m1(n);
			vfloat64m1_t vec_a = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_mul = vfmul_vf_f64m1(vec_a, x.v, vl);
			vse64_v_f64m1((double *)a, vec_mul, vl);
		}
		// n = (size_t)1 << logn;
		// a -= n;
	}
	else {
		for (size_t u = 0; u < n; u ++) {
			a[u] = fpr_mul(a[u], x);
		}
	}
}


/* see inner.h */
void
Zf(poly_div_fft_custom)(
	fpr *restrict a, const fpr *restrict b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		for (size_t vl; hn > 0; hn -= vl, a += vl, b += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_b_re = vle64_v_f64m1((double *)b, vl);
			vfloat64m1_t vec_a_re = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_b_im = vle64_v_f64m1((double *)b + hN, vl);
			vfloat64m1_t vec_a_im = vle64_v_f64m1((double *)a + hN, vl);
			vfloat64m1_t vec_sqr = vfmul_vv_f64m1(vec_b_re, vec_b_re, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_a_re, vec_b_re, vl);
			vfloat64m1_t vec_mul_im = vfmul_vv_f64m1(vec_a_re, vec_b_im, vl);
			vec_sqr = vfmacc_vv_f64m1(vec_sqr, vec_b_im, vec_b_im, vl);
			vec_mul_re = vfmacc_vv_f64m1(vec_mul_re, vec_a_im, vec_b_im, vl);
			vec_mul_im = vfmsac_vv_f64m1(vec_mul_im, vec_a_im, vec_b_re, vl);
			vfloat64m1_t vec_sqr_inv = vfrdiv_vf_f64m1(vec_sqr, 1.0f, vl);
			vec_mul_re = vfmul_vv_f64m1(vec_mul_re, vec_sqr_inv, vl);
			vec_mul_im = vfmul_vv_f64m1(vec_mul_im, vec_sqr_inv, vl);
			vse64_v_f64m1((double *)a, vec_mul_re, vl);
			vse64_v_f64m1((double *)a + hN, vec_mul_im, vl);
		}
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			fpr a_re, a_im, b_re, b_im;
			a_re = a[u];
			a_im = a[u + hn];
			b_re = b[u];
			b_im = b[u + hn];
			FPC_DIV(a[u], a[u + hn], a_re, a_im, b_re, b_im);
		}
	}
}


/* see inner.h */
void
Zf(poly_invnorm2_fft_custom)(fpr *restrict d,
	const fpr *restrict a, const fpr *restrict b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		for (size_t vl; hn > 0; hn -= vl, a += vl, b += vl, d += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_a_re = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_b_re = vle64_v_f64m1((double *)b, vl);
			vfloat64m1_t vec_a_im = vle64_v_f64m1((double *)a + hN, vl);
			vfloat64m1_t vec_b_im = vle64_v_f64m1((double *)b + hN, vl);
			vfloat64m1_t vec_a_sqr = vfmul_vv_f64m1(vec_a_re, vec_a_re, vl);
			vfloat64m1_t vec_b_sqr = vfmul_vv_f64m1(vec_b_re, vec_b_re, vl);
			vec_a_sqr = vfmacc_vv_f64m1(vec_a_sqr, vec_a_im, vec_a_im, vl);
			vec_b_sqr = vfmacc_vv_f64m1(vec_b_sqr, vec_b_im, vec_b_im, vl);
			vfloat64m1_t vec_add = vfadd_vv_f64m1(vec_a_sqr, vec_b_sqr, vl);
			vec_add = vfrdiv_vf_f64m1(vec_add, 1.0f, vl);
			vse64_v_f64m1((double *)d, vec_add, vl);
		}
		// d -= hN;
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			fpr a_re, a_im;
			fpr b_re, b_im;
			a_re = a[u];
			a_im = a[u + hn];
			b_re = b[u];
			b_im = b[u + hn];
			d[u] = fpr_inv(fpr_add(
				fpr_add(fpr_sqr(a_re), fpr_sqr(a_im)),
				fpr_add(fpr_sqr(b_re), fpr_sqr(b_im))));
		}
	}
}


/* see inner.h */
void
Zf(poly_add_muladj_fft_custom)(fpr *restrict d,
	const fpr *restrict F, const fpr *restrict G,
	const fpr *restrict f, const fpr *restrict g, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		for (size_t vl; hn > 0; hn -= vl, F += vl, G += vl, f += vl, g += vl, d += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_f_re = vle64_v_f64m1((double *)f, vl);
			vfloat64m1_t vec_F_re = vle64_v_f64m1((double *)F, vl);
			vfloat64m1_t vec_f_im = vle64_v_f64m1((double *)f + hN, vl);
			vfloat64m1_t vec_F_im = vle64_v_f64m1((double *)F + hN, vl);
			vfloat64m1_t vec_a_re = vfmul_vv_f64m1(vec_F_re, vec_f_re, vl);
			vfloat64m1_t vec_a_im = vfmul_vv_f64m1(vec_F_re, vec_f_im, vl);
			vec_a_re = vfmacc_vv_f64m1(vec_a_re, vec_F_im, vec_f_im, vl);
			vec_a_im = vfmsac_vv_f64m1(vec_a_im, vec_F_im, vec_f_re, vl);
			vfloat64m1_t vec_g_re = vle64_v_f64m1((double *)g, vl);
			vfloat64m1_t vec_G_re = vle64_v_f64m1((double *)G, vl);
			vfloat64m1_t vec_g_im = vle64_v_f64m1((double *)g + hN, vl);
			vfloat64m1_t vec_G_im = vle64_v_f64m1((double *)G + hN, vl);
			vfloat64m1_t vec_b_re = vfmul_vv_f64m1(vec_G_re, vec_g_re, vl);
			vfloat64m1_t vec_b_im = vfmul_vv_f64m1(vec_G_re, vec_g_im, vl);
			vec_b_re = vfmacc_vv_f64m1(vec_b_re, vec_G_im, vec_g_im, vl);
			vec_b_im = vfmsac_vv_f64m1(vec_b_im, vec_G_im, vec_g_re, vl);
			vfloat64m1_t vec_add_re = vfadd_vv_f64m1(vec_a_re, vec_b_re, vl);
			vfloat64m1_t vec_add_im = vfadd_vv_f64m1(vec_a_im, vec_b_im, vl);
			vse64_v_f64m1((double *)d, vec_add_re, vl);
			vse64_v_f64m1((double *)d + hN, vec_add_im, vl);
		}
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			fpr F_re, F_im, G_re, G_im;
			fpr f_re, f_im, g_re, g_im;
			fpr a_re, a_im, b_re, b_im;
			F_re = F[u];
			F_im = F[u + hn];
			G_re = G[u];
			G_im = G[u + hn];
			f_re = f[u];
			f_im = f[u + hn];
			g_re = g[u];
			g_im = g[u + hn];
			FPC_MUL(a_re, a_im, F_re, F_im, f_re, fpr_neg(f_im));
			FPC_MUL(b_re, b_im, G_re, G_im, g_re, fpr_neg(g_im));
			d[u] = fpr_add(a_re, b_re);
			d[u + hn] = fpr_add(a_im, b_im);
		}
	}
}


/* see inner.h */
void
Zf(poly_mul_autoadj_fft_custom)(
	fpr *restrict a, const fpr *restrict b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		for (size_t vl; hn > 0; hn -= vl, a += vl, b += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_b_re = vle64_v_f64m1((double *)b, vl);
			vfloat64m1_t vec_a_re = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_a_im = vle64_v_f64m1((double *)a + hN, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_a_re, vec_b_re, vl);
			vfloat64m1_t vec_mul_im = vfmul_vv_f64m1(vec_a_im, vec_b_re, vl);
			vse64_v_f64m1((double *)a, vec_mul_re, vl);
			vse64_v_f64m1((double *)a + hN, vec_mul_im, vl);
		}
		// a -= hN;
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			a[u] = fpr_mul(a[u], b[u]);
			a[u + hn] = fpr_mul(a[u + hn], b[u]);
		}
	}
}


/* see inner.h */
void
Zf(poly_div_autoadj_fft_custom)(
	fpr *restrict a, const fpr *restrict b, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		for (size_t vl; hn > 0; hn -= vl, a += vl, b += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_b_re = vle64_v_f64m1((double *)b, vl);
			vfloat64m1_t vec_a_re = vle64_v_f64m1((double *)a, vl);
			vfloat64m1_t vec_a_im = vle64_v_f64m1((double *)a + hN, vl);
			vfloat64m1_t vec_b_inv = vfrdiv_vf_f64m1(vec_b_re, 1.0f, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_a_re, vec_b_inv, vl);
			vfloat64m1_t vec_mul_im = vfmul_vv_f64m1(vec_a_im, vec_b_inv, vl);
			vse64_v_f64m1((double *)a, vec_mul_re, vl);
			vse64_v_f64m1((double *)a + hN, vec_mul_im, vl);
		}
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			fpr ib;
			ib = fpr_inv(b[u]);
			a[u] = fpr_mul(a[u], ib);
			a[u + hn] = fpr_mul(a[u + hn], ib);
		}
	}
}


/* see inner.h */
void
Zf(poly_LDL_fft_custom)(
	const fpr *restrict g00,
	fpr *restrict g01, fpr *restrict g11, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		for (size_t vl; hn > 0; hn -= vl, g00 += vl, g01 += vl, g11 += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_g00_re = vle64_v_f64m1((double *)g00, vl);
			vfloat64m1_t vec_g01_re = vle64_v_f64m1((double *)g01, vl);
			vfloat64m1_t vec_g00_im = vle64_v_f64m1((double *)g00 + hN, vl);
			vfloat64m1_t vec_g01_im = vle64_v_f64m1((double *)g01 + hN, vl);
			vfloat64m1_t vec_sqr = vfmul_vv_f64m1(vec_g00_re, vec_g00_re, vl);
			vfloat64m1_t vec_mu_re = vfmul_vv_f64m1(vec_g01_re, vec_g00_re, vl);
			vfloat64m1_t vec_mu_im = vfmul_vv_f64m1(vec_g01_re, vec_g00_im, vl);
			vec_sqr = vfmacc_vv_f64m1(vec_sqr, vec_g00_im, vec_g00_im, vl);
			vec_mu_re = vfmacc_vv_f64m1(vec_mu_re, vec_g01_im, vec_g00_im, vl);
			vec_mu_im = vfnmsac_vv_f64m1(vec_mu_im, vec_g01_im, vec_g00_re, vl);
			vfloat64m1_t vec_sqr_inv = vfrdiv_vf_f64m1(vec_sqr, 1.0f, vl);
			vec_mu_re = vfmul_vv_f64m1(vec_mu_re, vec_sqr_inv, vl);
			vec_mu_im = vfmul_vv_f64m1(vec_mu_im, vec_sqr_inv, vl);
			vfloat64m1_t vec_g11_re = vle64_v_f64m1((double *)g11, vl);
			vfloat64m1_t vec_g11_im = vle64_v_f64m1((double *)g11 + hN, vl);
			vse64_v_f64m1((double *)g01, vec_mu_re, vl);
			vse64_v_f64m1((double *)g01 + hN, vec_mu_im, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_mu_re, vec_g01_re, vl);
			vfloat64m1_t vec_mul_im = vfmul_vv_f64m1(vec_mu_re, vec_g01_im, vl);
			vec_mul_re = vfnmsac_vv_f64m1(vec_mul_re, vec_mu_im, vec_g01_im, vl);
			vec_mul_im = vfnmacc_vv_f64m1(vec_mul_im, vec_mu_im, vec_g01_re, vl);
			vfloat64m1_t vec_sub_re = vfsub_vv_f64m1(vec_g11_re, vec_mul_re, vl);
			vfloat64m1_t vec_sub_im = vfsub_vv_f64m1(vec_g11_im, vec_mul_im, vl);
			vse64_v_f64m1((double *)g11, vec_sub_re, vl);
			vse64_v_f64m1((double *)g11 + hN, vec_sub_im, vl);
		}
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			fpr g00_re, g00_im, g01_re, g01_im, g11_re, g11_im;
			fpr mu_re, mu_im;
			g00_re = g00[u];
			g00_im = g00[u + hn];
			g01_re = g01[u];
			g01_im = g01[u + hn];
			g11_re = g11[u];
			g11_im = g11[u + hn];
			FPC_DIV(mu_re, mu_im, g01_re, g01_im, g00_re, g00_im);
			FPC_MUL(g01_re, g01_im, mu_re, mu_im, g01_re, fpr_neg(g01_im));
			FPC_SUB(g11[u], g11[u + hn], g11_re, g11_im, g01_re, g01_im);
			g01[u] = mu_re;
			g01[u + hn] = fpr_neg(mu_im);
		}
	}
}


/* see inner.h */
void
Zf(poly_LDLmv_fft_custom)(
	fpr *restrict d11, fpr *restrict l10,
	const fpr *restrict g00, const fpr *restrict g01,
	const fpr *restrict g11, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	if (hn >= VLMAXFFT) {
		size_t hN = hn;
		for (size_t vl; hn > 0; hn -= vl, g00 += vl, g01 += vl, g11 += vl, d11 += vl, l10 += vl) {
			vl = vsetvl_e64m1(hn);
			vfloat64m1_t vec_g00_re = vle64_v_f64m1((double *)g00, vl);
			vfloat64m1_t vec_g01_re = vle64_v_f64m1((double *)g01, vl);
			vfloat64m1_t vec_g00_im = vle64_v_f64m1((double *)g00 + hN, vl);
			vfloat64m1_t vec_g01_im = vle64_v_f64m1((double *)g01 + hN, vl);
			vfloat64m1_t vec_sqr = vfmul_vv_f64m1(vec_g00_re, vec_g00_re, vl);
			vfloat64m1_t vec_mu_re = vfmul_vv_f64m1(vec_g01_re, vec_g00_re, vl);
			vfloat64m1_t vec_mu_im = vfmul_vv_f64m1(vec_g01_re, vec_g00_im, vl);
			vec_sqr = vfmacc_vv_f64m1(vec_sqr, vec_g00_im, vec_g00_im, vl);
			vec_mu_re = vfmacc_vv_f64m1(vec_mu_re, vec_g01_im, vec_g00_im, vl);
			vec_mu_im = vfnmsac_vv_f64m1(vec_mu_im, vec_g01_im, vec_g00_re, vl);
			vfloat64m1_t vec_sqr_inv = vfrdiv_vf_f64m1(vec_sqr, 1.0f, vl);
			vec_mu_re = vfmul_vv_f64m1(vec_mu_re, vec_sqr_inv, vl);
			vec_mu_im = vfmul_vv_f64m1(vec_mu_im, vec_sqr_inv, vl);
			vfloat64m1_t vec_g11_re = vle64_v_f64m1((double *)g11, vl);
			vfloat64m1_t vec_g11_im = vle64_v_f64m1((double *)g11 + hN, vl);
			vse64_v_f64m1((double *)l10, vec_mu_re, vl);
			vse64_v_f64m1((double *)l10 + hN, vec_mu_im, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_mu_re, vec_g01_re, vl);
			vfloat64m1_t vec_mul_im = vfmul_vv_f64m1(vec_mu_re, vec_g01_im, vl);
			vec_mul_re = vfnmsac_vv_f64m1(vec_mul_re, vec_mu_im, vec_g01_im, vl);
			vec_mul_im = vfnmacc_vv_f64m1(vec_mul_im, vec_mu_im, vec_g01_re, vl);
			vfloat64m1_t vec_sub_re = vfsub_vv_f64m1(vec_g11_re, vec_mul_re, vl);
			vfloat64m1_t vec_sub_im = vfsub_vv_f64m1(vec_g11_im, vec_mul_im, vl);
			vse64_v_f64m1((double *)d11, vec_sub_re, vl);
			vse64_v_f64m1((double *)d11 + hN, vec_sub_im, vl);
		}
	}
	else {
		for (size_t u = 0; u < hn; u ++) {
			fpr g00_re, g00_im, g01_re, g01_im, g11_re, g11_im;
			fpr mu_re, mu_im;
			g00_re = g00[u];
			g00_im = g00[u + hn];
			g01_re = g01[u];
			g01_im = g01[u + hn];
			g11_re = g11[u];
			g11_im = g11[u + hn];
			FPC_DIV(mu_re, mu_im, g01_re, g01_im, g00_re, g00_im);
			FPC_MUL(g01_re, g01_im, mu_re, mu_im, g01_re, fpr_neg(g01_im));
			FPC_SUB(d11[u], d11[u + hn], g11_re, g11_im, g01_re, g01_im);
			l10[u] = mu_re;
			l10[u + hn] = fpr_neg(mu_im);
		}
	}
}


/* see inner.h */
void
Zf(poly_split_fft_custom)(
	fpr *restrict f0, fpr *restrict f1,
	const fpr *restrict f, unsigned logn)
{
	/*
	 * The FFT representation we use is in bit-reversed order
	 * (element i contains f(w^(rev(i))), where rev() is the
	 * bit-reversal function over the ring degree. This changes
	 * indexes with regards to the Falcon specification.
	 */
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t qn = hn >> 1;
	if (qn >= VLMAXFFT) {
		size_t qN = qn;
    	double *gm = (double *)&fpr_gm_tab_reduce[hn];
		for (size_t vl; qn > 0; qn -= vl, f += (vl << 1), gm += (vl << 1), f0 += vl, f1 += vl) {
			vl = vsetvl_e64m1(qn);
			vfloat64m1_t vec_a_re = vlse64_v_f64m1((double *)f, 16, vl);
			vfloat64m1_t vec_b_re = vlse64_v_f64m1((double *)f + 1, 16, vl);
			vfloat64m1_t vec_a_im = vlse64_v_f64m1((double *)f + hn, 16, vl);
			vfloat64m1_t vec_b_im = vlse64_v_f64m1((double *)f + 1 + hn, 16, vl);
			vfloat64m1_t vec_add_re = vfadd_vv_f64m1(vec_a_re, vec_b_re, vl);
			vfloat64m1_t vec_add_im = vfadd_vv_f64m1(vec_a_im, vec_b_im, vl);
			vfloat64m1_t vec_sub_re = vfsub_vv_f64m1(vec_a_re, vec_b_re, vl);
			vfloat64m1_t vec_sub_im = vfsub_vv_f64m1(vec_a_im, vec_b_im, vl);
			vfloat64m1_t vec_add_half_re = vfmul_vf_f64m1(vec_add_re, 0.5, vl);
			vfloat64m1_t vec_add_half_im = vfmul_vf_f64m1(vec_add_im, 0.5, vl);
			vfloat64m1_t vec_gm_re = vlse64_v_f64m1(gm, 16, vl);
			vfloat64m1_t vec_gm_im = vlse64_v_f64m1(gm + 1, 16, vl);
			vse64_v_f64m1((double *)f0, vec_add_half_re, vl);
			vse64_v_f64m1((double *)f0 + qN, vec_add_half_im, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_sub_re, vec_gm_re, vl);
			vfloat64m1_t vec_mul_im = vfmul_vv_f64m1(vec_sub_re, vec_gm_im, vl);
			vec_mul_re = vfmacc_vv_f64m1(vec_mul_re, vec_sub_im, vec_gm_im, vl);
			vec_mul_im = vfmsac_vv_f64m1(vec_mul_im, vec_sub_im, vec_gm_re, vl);
			vfloat64m1_t vec_mul_half_re = vfmul_vf_f64m1(vec_mul_re, 0.5, vl);
			vfloat64m1_t vec_mul_half_im = vfmul_vf_f64m1(vec_mul_im, 0.5, vl);
			vse64_v_f64m1((double *)f1, vec_mul_half_re, vl);
			vse64_v_f64m1((double *)f1 + qN, vec_mul_half_im, vl);
		}
	}
	else {
		/*
		* We process complex values by pairs. For logn = 1, there is only
		* one complex value (the other one is the implicit conjugate),
		* so we add the two lines below because the loop will be
		* skipped.
		*/
		f0[0] = f[0];
		f1[0] = f[hn];
		for (size_t u = 0; u < qn; u ++) {
			fpr a_re, a_im, b_re, b_im;
			fpr t_re, t_im;
			a_re = f[(u << 1) + 0];
			a_im = f[(u << 1) + 0 + hn];
			b_re = f[(u << 1) + 1];
			b_im = f[(u << 1) + 1 + hn];
			FPC_ADD(t_re, t_im, a_re, a_im, b_re, b_im);
			f0[u] = fpr_half(t_re);
			f0[u + qn] = fpr_half(t_im);
			FPC_SUB(t_re, t_im, a_re, a_im, b_re, b_im);
			FPC_MUL(t_re, t_im, t_re, t_im,
				fpr_gm_tab_reduce[hn + (u << 1) + 0],
				fpr_neg(fpr_gm_tab_reduce[hn + (u << 1) + 1]));
			f1[u] = fpr_half(t_re);
			f1[u + qn] = fpr_half(t_im);
		}
	}
}


/* see inner.h */
void
Zf(poly_merge_fft_custom)(
	fpr *restrict f,
	const fpr *restrict f0, const fpr *restrict f1, unsigned logn)
{
	size_t n = (size_t)1 << logn;
	size_t hn = n >> 1;
	size_t qn = hn >> 1;
	if (qn >= VLMAXFFT) {
		size_t qN = qn;
    	double *gm = (double *)&fpr_gm_tab_reduce[hn];
		for (size_t vl; qn > 0; qn -= vl, f += (vl << 1), gm += (vl << 1), f0 += vl, f1 += vl) {
			vl = vsetvl_e64m1(qn);
			vfloat64m1_t vec_gm_re = vlse64_v_f64m1(gm, 16, vl);
			vfloat64m1_t vec_gm_im = vlse64_v_f64m1(gm + 1, 16, vl);
			vfloat64m1_t vec_b_re = vle64_v_f64m1((double *)f1, vl);
			vfloat64m1_t vec_b_im = vle64_v_f64m1((double *)f1 + qN, vl);
			vfloat64m1_t vec_a_re = vle64_v_f64m1((double *)f0, vl);
			vfloat64m1_t vec_a_im = vle64_v_f64m1((double *)f0 + qN, vl);
			vfloat64m1_t vec_mul_re = vfmul_vv_f64m1(vec_b_re, vec_gm_re, vl);
			vfloat64m1_t vec_mul_im = vfmul_vv_f64m1(vec_b_re, vec_gm_im, vl);
			vec_mul_re = vfnmsac_vv_f64m1(vec_mul_re, vec_b_im, vec_gm_im, vl);
			vec_mul_im = vfmacc_vv_f64m1(vec_mul_im, vec_b_im, vec_gm_re, vl);
			vfloat64m1_t vec_add_re = vfadd_vv_f64m1(vec_a_re, vec_mul_re, vl);
			vfloat64m1_t vec_sub_re = vfsub_vv_f64m1(vec_a_re, vec_mul_re, vl);
			vfloat64m1_t vec_add_im = vfadd_vv_f64m1(vec_a_im, vec_mul_im, vl);
			vfloat64m1_t vec_sub_im = vfsub_vv_f64m1(vec_a_im, vec_mul_im, vl);
			vsse64_v_f64m1((double *)f, 16, vec_add_re, vl);
			vsse64_v_f64m1((double *)f + 1, 16, vec_sub_re, vl);
			vsse64_v_f64m1((double *)f + hn, 16, vec_add_im, vl);
			vsse64_v_f64m1((double *)f + 1 + hn, 16, vec_sub_im, vl);
		}
	}
	else {
		/*
		* An extra copy to handle the special case logn = 1.
		*/
		f[0] = f0[0];
		f[hn] = f1[0];
		for (size_t u = 0; u < qn; u ++) {
			fpr a_re, a_im, b_re, b_im;
			fpr t_re, t_im;

			a_re = f0[u];
			a_im = f0[u + qn];
			FPC_MUL(b_re, b_im, f1[u], f1[u + qn],
				fpr_gm_tab_reduce[hn + (u << 1) + 0],
				fpr_gm_tab_reduce[hn + (u << 1) + 1]);
			FPC_ADD(t_re, t_im, a_re, a_im, b_re, b_im);
			f[(u << 1) + 0] = t_re;
			f[(u << 1) + 0 + hn] = t_im;
			FPC_SUB(t_re, t_im, a_re, a_im, b_re, b_im);
			f[(u << 1) + 1] = t_re;
			f[(u << 1) + 1 + hn] = t_im;
		}
	}
}