#include "inner.h"
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"
#include "params_custom.h"
#include <stdio.h>
/* =================================================================== */

/*
 * Compute degree N from logarithm 'logn'.
 */
#define MKN(logn)   ((size_t)1 << (logn))

/* =================================================================== */
/*
 * Binary case:
 *   N = 2^logn
 *   phi = X^N+1
 */

/*
 * Get the size of the LDL tree for an input with polynomials of size
 * 2^logn. The size is expressed in the number of elements.
 */
static inline unsigned
ffLDL_treesize(unsigned logn)
{
	/*
	 * For logn = 0 (polynomials are constant), the "tree" is a
	 * single element. Otherwise, the tree node has size 2^logn, and
	 * has two child trees for size logn-1 each. Thus, treesize s()
	 * must fulfill these two relations:
	 *
	 *   s(0) = 1
	 *   s(logn) = (2^logn) + 2*s(logn-1)
	 */
	return (logn + 1) << logn;
}

/*
 * Inner function for ffLDL_fft(). It expects the matrix to be both
 * auto-adjoint and quasicyclic; also, it uses the source operands
 * as modifiable temporaries.
 *
 * tmp[] must have room for at least one polynomial.
 */
static void
ffLDL_fft_inner_custom_with_rvvf(fpr *restrict tree,
	fpr *restrict g0, fpr *restrict g1, unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;

	n = MKN(logn);
	if (n == 1) {
		tree[0] = g0[0];
		return;
	}
	hn = n >> 1;

	/*
	 * The LDL decomposition yields L (which is written in the tree)
	 * and the diagonal of D. Since d00 = g0, we just write d11
	 * into tmp.
	 */
	Zf(poly_LDLmv_fft_custom)(tmp, tree, g0, g1, g0, logn);

	/*
	 * Split d00 (currently in g0) and d11 (currently in tmp). We
	 * reuse g0 and g1 as temporary storage spaces:
	 *   d00 splits into g1, g1+hn
	 *   d11 splits into g0, g0+hn
	 */
	Zf(poly_split_fft_custom)(g1, g1 + hn, g0, logn);
	Zf(poly_split_fft_custom)(g0, g0 + hn, tmp, logn);

	/*
	 * Each split result is the first row of a new auto-adjoint
	 * quasicyclic matrix for the next recursive step.
	 */
	ffLDL_fft_inner_custom_with_rvvf(tree + n,
		g1, g1 + hn, logn - 1, tmp);
	ffLDL_fft_inner_custom_with_rvvf(tree + n + ffLDL_treesize(logn - 1),
		g0, g0 + hn, logn - 1, tmp);
}

/*
 * Compute the ffLDL tree of an auto-adjoint matrix G. The matrix
 * is provided as three polynomials (FFT representation).
 *
 * The "tree" array is filled with the computed tree, of size
 * (logn+1)*(2^logn) elements (see ffLDL_treesize()).
 *
 * Input arrays MUST NOT overlap, except possibly the three unmodified
 * arrays g00, g01 and g11. tmp[] should have room for at least three
 * polynomials of 2^logn elements each.
 */
static void
ffLDL_fft_custom_with_rvvf(fpr *restrict tree, const fpr *restrict g00,
	const fpr *restrict g01, const fpr *restrict g11,
	unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;
	fpr *d00, *d11;

	n = MKN(logn);
	if (n == 1) {
		tree[0] = g00[0];
		return;
	}
	hn = n >> 1;
	d00 = tmp;
	d11 = tmp + n;
	tmp += n << 1;

	memcpy(d00, g00, n * sizeof *g00);
	Zf(poly_LDLmv_fft_custom)(d11, tree, g00, g01, g11, logn);

	Zf(poly_split_fft_custom)(tmp, tmp + hn, d00, logn);
	Zf(poly_split_fft_custom)(d00, d00 + hn, d11, logn);
	memcpy(d11, tmp, n * sizeof *tmp);
	ffLDL_fft_inner_custom_with_rvvf(tree + n,
		d11, d11 + hn, logn - 1, tmp);
	ffLDL_fft_inner_custom_with_rvvf(tree + n + ffLDL_treesize(logn - 1),
		d00, d00 + hn, logn - 1, tmp);
}

/*
 * Normalize an ffLDL tree: each leaf of value x is replaced with
 * sigma / sqrt(x).
 */
static void
ffLDL_binary_normalize(fpr *tree, unsigned orig_logn, unsigned logn)
{
	/*
	 * TODO: make an iterative version.
	 */
	size_t n;

	n = MKN(logn);
	if (n == 1) {
		/*
		 * We actually store in the tree leaf the inverse of
		 * the value mandated by the specification: this
		 * saves a division both here and in the sampler.
		 */
		tree[0] = fpr_mul(fpr_sqrt(tree[0]), fpr_inv_sigma[orig_logn]);
	} else {
		ffLDL_binary_normalize(tree + n, orig_logn, logn - 1);
		ffLDL_binary_normalize(tree + n + ffLDL_treesize(logn - 1),
			orig_logn, logn - 1);
	}
}

/* =================================================================== */

/*
 * Convert an integer polynomial (with small values) into the
 * representation with complex numbers.
 */
static void
smallints_to_fpr(fpr *r, const int8_t *t, unsigned logn)
{
	size_t n, u;

	n = MKN(logn);
	for (u = 0; u < n; u ++) {
		r[u] = fpr_of(t[u]);
	}
}

/*
 * The expanded private key contains:
 *  - The B0 matrix (four elements)
 *  - The ffLDL tree
 */

static inline size_t
skoff_b00(unsigned logn)
{
	(void)logn;
	return 0;
}

static inline size_t
skoff_b01(unsigned logn)
{
	return MKN(logn);
}

static inline size_t
skoff_b10(unsigned logn)
{
	return 2 * MKN(logn);
}

static inline size_t
skoff_b11(unsigned logn)
{
	return 3 * MKN(logn);
}

static inline size_t
skoff_tree(unsigned logn)
{
	return 4 * MKN(logn);
}

/* see inner.h */
void
Zf(expand_privkey_custom_with_rvvf)(fpr *restrict expanded_key,
	const int8_t *f, const int8_t *g,
	const int8_t *F, const int8_t *G,
	unsigned logn, uint8_t *restrict tmp)
{
	size_t n;
	fpr *rf, *rg, *rF, *rG;
	fpr *b00, *b01, *b10, *b11;
	fpr *g00, *g01, *g11, *gxx;
	fpr *tree;

	n = MKN(logn);
	b00 = expanded_key + skoff_b00(logn);
	b01 = expanded_key + skoff_b01(logn);
	b10 = expanded_key + skoff_b10(logn);
	b11 = expanded_key + skoff_b11(logn);
	tree = expanded_key + skoff_tree(logn);

	/*
	 * We load the private key elements directly into the B0 matrix,
	 * since B0 = [[g, -f], [G, -F]].
	 */
	rf = b01;
	rg = b00;
	rF = b11;
	rG = b10;

	smallints_to_fpr(rf, f, logn);
	smallints_to_fpr(rg, g, logn);
	smallints_to_fpr(rF, F, logn);
	smallints_to_fpr(rG, G, logn);

	/*
	 * Compute the FFT for the key elements, and negate f and F.
	 */
	Zf(FFT_custom)(rf, logn);
	Zf(FFT_custom)(rg, logn);
	Zf(FFT_custom)(rF, logn);
	Zf(FFT_custom)(rG, logn);
	Zf(poly_neg_custom)(rf, logn);
	Zf(poly_neg_custom)(rF, logn);

	/*
	 * The Gram matrix is G = B·B*. Formulas are:
	 *   g00 = b00*adj(b00) + b01*adj(b01)
	 *   g01 = b00*adj(b10) + b01*adj(b11)
	 *   g10 = b10*adj(b00) + b11*adj(b01)
	 *   g11 = b10*adj(b10) + b11*adj(b11)
	 *
	 * For historical reasons, this implementation uses
	 * g00, g01 and g11 (upper triangle).
	 */
	g00 = (fpr *)tmp;
	g01 = g00 + n;
	g11 = g01 + n;
	gxx = g11 + n;

	memcpy(g00, b00, n * sizeof *b00);
	Zf(poly_mulselfadj_fft_custom)(g00, logn);
	memcpy(gxx, b01, n * sizeof *b01);
	Zf(poly_mulselfadj_fft_custom)(gxx, logn);
	Zf(poly_add_custom)(g00, gxx, logn);

	memcpy(g01, b00, n * sizeof *b00);
	Zf(poly_muladj_fft_custom)(g01, b10, logn);
	memcpy(gxx, b01, n * sizeof *b01);
	Zf(poly_muladj_fft_custom)(gxx, b11, logn);
	Zf(poly_add_custom)(g01, gxx, logn);

	memcpy(g11, b10, n * sizeof *b10);
	Zf(poly_mulselfadj_fft_custom)(g11, logn);
	memcpy(gxx, b11, n * sizeof *b11);
	Zf(poly_mulselfadj_fft_custom)(gxx, logn);
	Zf(poly_add_custom)(g11, gxx, logn);

	/*
	 * Compute the Falcon tree.
	 */
	ffLDL_fft_custom_with_rvvf(tree, g00, g01, g11, logn, gxx);

	/*
	 * Normalize tree.
	 */
	ffLDL_binary_normalize(tree, logn, logn);
}

typedef int (*samplerZ)(void *ctx, fpr mu, fpr sigma);
typedef fpr (*samplerZ_custom)(fpr mu, fpr sigma);

/*
 * Perform Fast Fourier Sampling for target vector t. The Gram matrix
 * is provided (G = [[g00, g01], [adj(g01), g11]]). The sampled vector
 * is written over (t0,t1). The Gram matrix is modified as well. The
 * tmp[] buffer must have room for four polynomials.
 */
static void
ffSampling_fft_dyntree_custom_with_rvvf(samplerZ_custom samp, 
	fpr *restrict t0, fpr *restrict t1,
	fpr *restrict g00, fpr *restrict g01, fpr *restrict g11,
	unsigned orig_logn, unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;
	fpr *z0, *z1;

	/*
	 * Deepest level: the LDL tree leaf value is just g00 (the
	 * array has length only 1 at this point); we normalize it
	 * with regards to sigma, then use it for sampling.
	 */
	if (logn == 0) {
		fpr leaf;

		leaf = g00[0];
		leaf = fpr_mul(fpr_sqrt(leaf), fpr_inv_sigma[orig_logn]);
		t0[0] = samp(t0[0], leaf);
		t1[0] = samp(t1[0], leaf);
		return;
	}

	n = (size_t)1 << logn;
	hn = n >> 1;

	/*
	 * Decompose G into LDL. We only need d00 (identical to g00),
	 * d11, and l10; we do that in place.
	 */
	Zf(poly_LDL_fft_custom)(g00, g01, g11, logn);

	/*
	 * Split d00 and d11 and expand them into half-size quasi-cyclic
	 * Gram matrices. We also save l10 in tmp[].
	 */
	Zf(poly_split_fft_custom)(tmp, tmp + hn, g00, logn);
	memcpy(g00, tmp, n * sizeof *tmp);
	Zf(poly_split_fft_custom)(tmp, tmp + hn, g11, logn);
	memcpy(g11, tmp, n * sizeof *tmp);
	memcpy(tmp, g01, n * sizeof *g01);
	memcpy(g01, g00, hn * sizeof *g00);
	memcpy(g01 + hn, g11, hn * sizeof *g00);

	/*
	 * The half-size Gram matrices for the recursive LDL tree
	 * building are now:
	 *   - left sub-tree: g00, g00+hn, g01
	 *   - right sub-tree: g11, g11+hn, g01+hn
	 * l10 is in tmp[].
	 */

	/*
	 * We split t1 and use the first recursive call on the two
	 * halves, using the right sub-tree. The result is merged
	 * back into tmp + 2*n.
	 */
	z1 = tmp + n;
	Zf(poly_split_fft_custom)(z1, z1 + hn, t1, logn);
	ffSampling_fft_dyntree_custom_with_rvvf(samp, z1, z1 + hn,
		g11, g11 + hn, g01 + hn, orig_logn, logn - 1, z1 + n);
	Zf(poly_merge_fft_custom)(tmp + (n << 1), z1, z1 + hn, logn);

	/*
	 * Compute tb0 = t0 + (t1 - z1) * l10.
	 * At that point, l10 is in tmp, t1 is unmodified, and z1 is
	 * in tmp + (n << 1). The buffer in z1 is free.
	 *
	 * In the end, z1 is written over t1, and tb0 is in t0.
	 */
	memcpy(z1, t1, n * sizeof *t1);
	Zf(poly_sub_custom)(z1, tmp + (n << 1), logn);
	memcpy(t1, tmp + (n << 1), n * sizeof *tmp);
	Zf(poly_mul_fft_custom)(tmp, z1, logn);
	Zf(poly_add_custom)(t0, tmp, logn);

	/*
	 * Second recursive invocation, on the split tb0 (currently in t0)
	 * and the left sub-tree.
	 */
	z0 = tmp;
	Zf(poly_split_fft_custom)(z0, z0 + hn, t0, logn);
	ffSampling_fft_dyntree_custom_with_rvvf(samp, z0, z0 + hn,
		g00, g00 + hn, g01, orig_logn, logn - 1, z0 + n);
	Zf(poly_merge_fft_custom)(t0, z0, z0 + hn, logn);
}


/*
 * Perform Fast Fourier Sampling for target vector t and LDL tree T.
 * tmp[] must have size for at least two polynomials of size 2^logn.
 */
static void
ffSampling_fft_custom_with_rvvf(samplerZ_custom samp,
	fpr *restrict z0, fpr *restrict z1,
	const fpr *restrict tree,
	const fpr *restrict t0, const fpr *restrict t1, unsigned logn,
	fpr *restrict tmp)
{
	size_t n, hn;
	const fpr *tree0, *tree1;

	/*
	 * When logn == 2, we inline the last two recursion levels.
	 */
	if (logn == 2) {
		fpr x0, x1, y0, y1, w0, w1, w2, w3, sigma;
		fpr a_re, a_im, b_re, b_im, c_re, c_im;

		tree0 = tree + 4;
		tree1 = tree + 8;

		/*
		 * We split t1 into w*, then do the recursive invocation,
		 * with output in w*. We finally merge back into z1.
		 */
		a_re = t1[0];
		a_im = t1[2];
		b_re = t1[1];
		b_im = t1[3];
		c_re = fpr_add(a_re, b_re);
		c_im = fpr_add(a_im, b_im);
		w0 = fpr_half(c_re);
		w1 = fpr_half(c_im);
		c_re = fpr_sub(a_re, b_re);
		c_im = fpr_sub(a_im, b_im);
		w2 = fpr_mul(fpr_add(c_re, c_im), fpr_invsqrt8);
		w3 = fpr_mul(fpr_sub(c_im, c_re), fpr_invsqrt8);

		x0 = w2;
		x1 = w3;
		sigma = tree1[3];
		w2 = samp(x0, sigma);
		w3 = samp(x1, sigma);
		a_re = fpr_sub(x0, w2);
		a_im = fpr_sub(x1, w3);
		b_re = tree1[0];
		b_im = tree1[1];
		c_re = fpr_sub(fpr_mul(a_re, b_re), fpr_mul(a_im, b_im));
		c_im = fpr_add(fpr_mul(a_re, b_im), fpr_mul(a_im, b_re));
		x0 = fpr_add(c_re, w0);
		x1 = fpr_add(c_im, w1);
		sigma = tree1[2];
		w0 = samp(x0, sigma);
		w1 = samp(x1, sigma);

		a_re = w0;
		a_im = w1;
		b_re = w2;
		b_im = w3;
		c_re = fpr_mul(fpr_sub(b_re, b_im), fpr_invsqrt2);
		c_im = fpr_mul(fpr_add(b_re, b_im), fpr_invsqrt2);
		z1[0] = w0 = fpr_add(a_re, c_re);
		z1[2] = w2 = fpr_add(a_im, c_im);
		z1[1] = w1 = fpr_sub(a_re, c_re);
		z1[3] = w3 = fpr_sub(a_im, c_im);

		/*
		 * Compute tb0 = t0 + (t1 - z1) * L. Value tb0 ends up in w*.
		 */
		w0 = fpr_sub(t1[0], w0);
		w1 = fpr_sub(t1[1], w1);
		w2 = fpr_sub(t1[2], w2);
		w3 = fpr_sub(t1[3], w3);

		a_re = w0;
		a_im = w2;
		b_re = tree[0];
		b_im = tree[2];
		w0 = fpr_sub(fpr_mul(a_re, b_re), fpr_mul(a_im, b_im));
		w2 = fpr_add(fpr_mul(a_re, b_im), fpr_mul(a_im, b_re));
		a_re = w1;
		a_im = w3;
		b_re = tree[1];
		b_im = tree[3];
		w1 = fpr_sub(fpr_mul(a_re, b_re), fpr_mul(a_im, b_im));
		w3 = fpr_add(fpr_mul(a_re, b_im), fpr_mul(a_im, b_re));

		w0 = fpr_add(w0, t0[0]);
		w1 = fpr_add(w1, t0[1]);
		w2 = fpr_add(w2, t0[2]);
		w3 = fpr_add(w3, t0[3]);

		/*
		 * Second recursive invocation.
		 */
		a_re = w0;
		a_im = w2;
		b_re = w1;
		b_im = w3;
		c_re = fpr_add(a_re, b_re);
		c_im = fpr_add(a_im, b_im);
		w0 = fpr_half(c_re);
		w1 = fpr_half(c_im);
		c_re = fpr_sub(a_re, b_re);
		c_im = fpr_sub(a_im, b_im);
		w2 = fpr_mul(fpr_add(c_re, c_im), fpr_invsqrt8);
		w3 = fpr_mul(fpr_sub(c_im, c_re), fpr_invsqrt8);

		x0 = w2;
		x1 = w3;
		sigma = tree0[3];
		w2 = y0 = samp(x0, sigma);
		w3 = y1 = samp(x1, sigma);
		a_re = fpr_sub(x0, y0);
		a_im = fpr_sub(x1, y1);
		b_re = tree0[0];
		b_im = tree0[1];
		c_re = fpr_sub(fpr_mul(a_re, b_re), fpr_mul(a_im, b_im));
		c_im = fpr_add(fpr_mul(a_re, b_im), fpr_mul(a_im, b_re));
		x0 = fpr_add(c_re, w0);
		x1 = fpr_add(c_im, w1);
		sigma = tree0[2];
		w0 = samp(x0, sigma);
		w1 = samp(x1, sigma);

		a_re = w0;
		a_im = w1;
		b_re = w2;
		b_im = w3;
		c_re = fpr_mul(fpr_sub(b_re, b_im), fpr_invsqrt2);
		c_im = fpr_mul(fpr_add(b_re, b_im), fpr_invsqrt2);
		z0[0] = fpr_add(a_re, c_re);
		z0[2] = fpr_add(a_im, c_im);
		z0[1] = fpr_sub(a_re, c_re);
		z0[3] = fpr_sub(a_im, c_im);

		return;
	}

	/*
	 * Case logn == 1 is reachable only when using Falcon-2 (the
	 * smallest size for which Falcon is mathematically defined, but
	 * of course way too insecure to be of any use).
	 */
	if (logn == 1) {
		fpr x0, x1, y0, y1, sigma;
		fpr a_re, a_im, b_re, b_im, c_re, c_im;

		x0 = t1[0];
		x1 = t1[1];
		sigma = tree[3];
		z1[0] = y0 = samp(x0, sigma);
		z1[1] = y1 = samp(x1, sigma);
		a_re = fpr_sub(x0, y0);
		a_im = fpr_sub(x1, y1);
		b_re = tree[0];
		b_im = tree[1];
		c_re = fpr_sub(fpr_mul(a_re, b_re), fpr_mul(a_im, b_im));
		c_im = fpr_add(fpr_mul(a_re, b_im), fpr_mul(a_im, b_re));
		x0 = fpr_add(c_re, t0[0]);
		x1 = fpr_add(c_im, t0[1]);
		sigma = tree[2];
		z0[0] = samp(x0, sigma);
		z0[1] = samp(x1, sigma);

		return;
	}

	/*
	 * Normal end of recursion is for logn == 0. Since the last
	 * steps of the recursions were inlined in the blocks above
	 * (when logn == 1 or 2), this case is not reachable, and is
	 * retained here only for documentation purposes.

	if (logn == 0) {
		fpr x0, x1, sigma;

		x0 = t0[0];
		x1 = t1[0];
		sigma = tree[0];
		z0[0] = samp(x0, sigma);
		z1[0] = samp(x1, sigma);
		return;
	}

	 */

	/*
	 * General recursive case (logn >= 3).
	 */

	n = (size_t)1 << logn;
	hn = n >> 1;
	tree0 = tree + n;
	tree1 = tree + n + ffLDL_treesize(logn - 1);

	/*
	 * We split t1 into z1 (reused as temporary storage), then do
	 * the recursive invocation, with output in tmp. We finally
	 * merge back into z1.
	 */
	Zf(poly_split_fft_custom)(z1, z1 + hn, t1, logn);
	ffSampling_fft_custom_with_rvvf(samp, tmp, tmp + hn,
		tree1, z1, z1 + hn, logn - 1, tmp + n);
	Zf(poly_merge_fft_custom)(z1, tmp, tmp + hn, logn);

	/*
	 * Compute tb0 = t0 + (t1 - z1) * L. Value tb0 ends up in tmp[].
	 */
	memcpy(tmp, t1, n * sizeof *t1);
	Zf(poly_sub_custom)(tmp, z1, logn);
	Zf(poly_mul_fft_custom)(tmp, tree, logn);
	Zf(poly_add_custom)(tmp, t0, logn);

	/*
	 * Second recursive invocation.
	 */
	Zf(poly_split_fft_custom)(z0, z0 + hn, tmp, logn);
	ffSampling_fft_custom_with_rvvf(samp, tmp, tmp + hn,
		tree0, z0, z0 + hn, logn - 1, tmp + n);
	Zf(poly_merge_fft_custom)(z0, tmp, tmp + hn, logn);
}


/*
 * Compute a signature: the signature contains two vectors, s1 and s2.
 * The s1 vector is not returned. The squared norm of (s1,s2) is
 * computed, and if it is short enough, then s2 is returned into the
 * s2[] buffer, and 1 is returned; otherwise, s2[] is untouched and 0 is
 * returned; the caller should then try again. This function uses an
 * expanded key.
 *
 * tmp[] must have room for at least six polynomials.
 */
static int
do_sign_tree_custom_with_rvvf(samplerZ_custom samp, int16_t *s2,
	const fpr *restrict expanded_key,
	const uint16_t *hm,
	unsigned logn, fpr *restrict tmp)
{
	size_t n, u;
	fpr *t0, *t1, *tx, *ty;
	const fpr *b00, *b01, *b10, *b11, *tree;
	fpr ni;
	int16_t *s1tmp, *s2tmp;

	n = MKN(logn);
	t0 = tmp;
	t1 = t0 + n;
	b00 = expanded_key + skoff_b00(logn);
	b01 = expanded_key + skoff_b01(logn);
	b10 = expanded_key + skoff_b10(logn);
	b11 = expanded_key + skoff_b11(logn);
	tree = expanded_key + skoff_tree(logn);

	/*
	 * Set the target vector to [hm, 0] (hm is the hashed message).
	 */
	for (u = 0; u < n; u ++) {
		t0[u] = fpr_of(hm[u]);
		/* This is implicit.
		t1[u] = fpr_zero;
		*/
	}

	/*
	 * Apply the lattice basis to obtain the real target
	 * vector (after normalization with regards to modulus).
	 */
	Zf(FFT_custom)(t0, logn);
	ni = fpr_inverse_of_q;
	memcpy(t1, t0, n * sizeof *t0);
	Zf(poly_mul_fft_custom)(t1, b01, logn);
	Zf(poly_mulconst_custom)(t1, fpr_neg(ni), logn);
	Zf(poly_mul_fft_custom)(t0, b11, logn);
	Zf(poly_mulconst_custom)(t0, ni, logn);

	tx = t1 + n;
	ty = tx + n;

	/*
	 * Apply sampling. Output is written back in [tx, ty].
	 */
	ffSampling_fft_custom_with_rvvf(samp, tx, ty, tree, t0, t1, logn, ty + n);

	/*
	 * Get the lattice point corresponding to that tiny vector.
	 */
	memcpy(t0, tx, n * sizeof *tx);
	memcpy(t1, ty, n * sizeof *ty);
	Zf(poly_mul_fft_custom)(tx, b00, logn);
	Zf(poly_mul_fft_custom)(ty, b10, logn);
	Zf(poly_add_custom)(tx, ty, logn);
	memcpy(ty, t0, n * sizeof *t0);
	Zf(poly_mul_fft_custom)(ty, b01, logn);

	memcpy(t0, tx, n * sizeof *tx);
	Zf(poly_mul_fft_custom)(t1, b11, logn);
	Zf(poly_add_custom)(t1, ty, logn);

	Zf(iFFT_custom)(t0, logn);
	Zf(iFFT_custom)(t1, logn);

	/*
	 * Compute the signature.
	 */
	int16_t *pt0 = (int16_t *)tx;
	int16_t *pt1 = (int16_t *)tmp;
	vint64m1_t scalar = vmv_v_x_i64m1(0, 1);
	size_t N = n;
    for (size_t vl; N > 0; N -= vl, t0 += vl, t1 += vl, pt0 += vl, pt1 += vl, hm += vl) {
        vl = vsetvl_e16m1(N);
        vuint16m1_t vec_hm = vle16_v_u16m1(hm, vl);
		vuint32m2_t vec_hm_uint32 = vzext_vf2_u32m2(vec_hm, vl);
		vint32m2_t vec_hm_int32 = vreinterpret_v_u32m2_i32m2(vec_hm_uint32);
        vfloat64m4_t vec_t0_f = vle64_v_f64m4((double *)t0, vl);
        vfloat64m4_t vec_t1_f = vle64_v_f64m4((double *)t1, vl);
        vint64m4_t vec_t0_int64 = vfcvt_x_f_v_i64m4(vec_t0_f, vl);
        vint64m4_t vec_t1_int64 = vfcvt_x_f_v_i64m4(vec_t1_f, vl);
		vint32m2_t vec_t0_int32 = vnsra_wx_i32m2(vec_t0_int64, 0, vl);
		vint32m2_t vec_t1_int32 = vnsra_wx_i32m2(vec_t1_int64, 0, vl);
		vint32m2_t vec_sub = vsub_vv_i32m2(vec_hm_int32, vec_t0_int32, vl);
		vint32m2_t vec_rsub = vrsub_vx_i32m2(vec_t1_int32, 0, vl);
		vint16m1_t vec_z = vnsra_wx_i16m1(vec_sub, 0, vl);
		vint16m1_t vec_y = vnsra_wx_i16m1(vec_rsub, 0, vl);
		vse16_v_i16m1(pt0, vec_z, vl);
		vse16_v_i16m1(pt1, vec_y, vl);
		vint32m2_t vec_mul = vwmul_vv_i32m2(vec_z, vec_z, vl);
		scalar = vwredsum_vs_i32m2_i64m1(scalar, vec_mul, scalar, vl);
    }
	uint64_t sqn = vmv_x_s_i64m1_i64(scalar);

	/*
	 * With "normal" degrees (e.g. 512 or 1024), it is very
	 * improbable that the computed vector is not short enough;
	 * however, it may happen in practice for the very reduced
	 * versions (e.g. degree 16 or below). In that case, the caller
	 * will loop, and we must not write anything into s2[] because
	 * s2[] may overlap with the hashed message hm[] and we need
	 * hm[] for the next iteration.
	 */
	s1tmp = (int16_t *)tx;
	s2tmp = (int16_t *)tmp;
	if (Zf(is_short_half_custom)(sqn, s2tmp, logn)) {
		memcpy(s2, s2tmp, n * sizeof *s2);
		memcpy(tmp, s1tmp, n * sizeof *s1tmp);
		return 1;
	}
	return 0;
}


/*
 * Compute a signature: the signature contains two vectors, s1 and s2.
 * The s1 vector is not returned. The squared norm of (s1,s2) is
 * computed, and if it is short enough, then s2 is returned into the
 * s2[] buffer, and 1 is returned; otherwise, s2[] is untouched and 0 is
 * returned; the caller should then try again.
 *
 * tmp[] must have room for at least nine polynomials.
 */
static int
do_sign_dyn_custom_with_rvvf(samplerZ_custom samp, int16_t *s2,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint16_t *hm, unsigned logn, fpr *restrict tmp)
{
	size_t n, u;
	fpr *t0, *t1, *tx, *ty;
	fpr *b00, *b01, *b10, *b11, *g00, *g01, *g11;
	fpr ni;
	int16_t *s1tmp, *s2tmp;

	n = MKN(logn);

	/*
	 * Lattice basis is B = [[g, -f], [G, -F]]. We convert it to FFT.
	 */
	b00 = tmp;
	b01 = b00 + n;
	b10 = b01 + n;
	b11 = b10 + n;
	smallints_to_fpr(b01, f, logn);
	smallints_to_fpr(b00, g, logn);
	smallints_to_fpr(b11, F, logn);
	smallints_to_fpr(b10, G, logn);
	Zf(FFT_custom)(b01, logn);
	Zf(FFT_custom)(b00, logn);
	Zf(FFT_custom)(b11, logn);
	Zf(FFT_custom)(b10, logn);
	Zf(poly_neg_custom)(b01, logn);
	Zf(poly_neg_custom)(b11, logn);

	/*
	 * Compute the Gram matrix G = B·B*. Formulas are:
	 *   g00 = b00*adj(b00) + b01*adj(b01)
	 *   g01 = b00*adj(b10) + b01*adj(b11)
	 *   g10 = b10*adj(b00) + b11*adj(b01)
	 *   g11 = b10*adj(b10) + b11*adj(b11)
	 *
	 * For historical reasons, this implementation uses
	 * g00, g01 and g11 (upper triangle). g10 is not kept
	 * since it is equal to adj(g01).
	 *
	 * We _replace_ the matrix B with the Gram matrix, but we
	 * must keep b01 and b11 for computing the target vector.
	 */
	t0 = b11 + n;
	t1 = t0 + n;

	memcpy(t0, b01, n * sizeof *b01);
	Zf(poly_mulselfadj_fft_custom)(t0, logn);    // t0 <- b01*adj(b01)

	memcpy(t1, b00, n * sizeof *b00);
	Zf(poly_muladj_fft_custom)(t1, b10, logn);   // t1 <- b00*adj(b10)
	Zf(poly_mulselfadj_fft_custom)(b00, logn);   // b00 <- b00*adj(b00)
	Zf(poly_add_custom)(b00, t0, logn);      // b00 <- g00
	memcpy(t0, b01, n * sizeof *b01);
	Zf(poly_muladj_fft_custom)(b01, b11, logn);  // b01 <- b01*adj(b11)
	Zf(poly_add_custom)(b01, t1, logn);      // b01 <- g01

	Zf(poly_mulselfadj_fft_custom)(b10, logn);   // b10 <- b10*adj(b10)
	memcpy(t1, b11, n * sizeof *b11);
	Zf(poly_mulselfadj_fft_custom)(t1, logn);    // t1 <- b11*adj(b11)
	Zf(poly_add_custom)(b10, t1, logn);      // b10 <- g11

	/*
	 * We rename variables to make things clearer. The three elements
	 * of the Gram matrix uses the first 3*n slots of tmp[], followed
	 * by b11 and b01 (in that order).
	 */
	g00 = b00;
	g01 = b01;
	g11 = b10;
	b01 = t0;
	t0 = b01 + n;
	t1 = t0 + n;

	/*
	 * Memory layout at that point:
	 *   g00 g01 g11 b11 b01 t0 t1
	 */

	/*
	 * Set the target vector to [hm, 0] (hm is the hashed message).
	 */
	for (u = 0; u < n; u ++) {
		t0[u] = fpr_of(hm[u]);
		/* This is implicit.
		t1[u] = fpr_zero;
		*/
	}

	/*
	 * Apply the lattice basis to obtain the real target
	 * vector (after normalization with regards to modulus).
	 */
	Zf(FFT_custom)(t0, logn);
	ni = fpr_inverse_of_q;
	memcpy(t1, t0, n * sizeof *t0);
	Zf(poly_mul_fft_custom)(t1, b01, logn);
	Zf(poly_mulconst_custom)(t1, fpr_neg(ni), logn);
	Zf(poly_mul_fft_custom)(t0, b11, logn);
	Zf(poly_mulconst_custom)(t0, ni, logn);

	/*
	 * b01 and b11 can be discarded, so we move back (t0,t1).
	 * Memory layout is now:
	 *      g00 g01 g11 t0 t1
	 */
	memcpy(b11, t0, n * 2 * sizeof *t0);
	t0 = g11 + n;
	t1 = t0 + n;

	/*
	 * Apply sampling; result is written over (t0,t1).
	 */
	ffSampling_fft_dyntree_custom_with_rvvf(samp,
		t0, t1, g00, g01, g11, logn, logn, t1 + n);

	/*
	 * We arrange the layout back to:
	 *     b00 b01 b10 b11 t0 t1
	 *
	 * We did not conserve the matrix basis, so we must recompute
	 * it now.
	 */
	b00 = tmp;
	b01 = b00 + n;
	b10 = b01 + n;
	b11 = b10 + n;
	memmove(b11 + n, t0, n * 2 * sizeof *t0);
	t0 = b11 + n;
	t1 = t0 + n;
	smallints_to_fpr(b01, f, logn);
	smallints_to_fpr(b00, g, logn);
	smallints_to_fpr(b11, F, logn);
	smallints_to_fpr(b10, G, logn);
	Zf(FFT_custom)(b01, logn);
	Zf(FFT_custom)(b00, logn);
	Zf(FFT_custom)(b11, logn);
	Zf(FFT_custom)(b10, logn);
	Zf(poly_neg_custom)(b01, logn);
	Zf(poly_neg_custom)(b11, logn);
	tx = t1 + n;
	ty = tx + n;

	/*
	 * Get the lattice point corresponding to that tiny vector.
	 */
	memcpy(tx, t0, n * sizeof *t0);
	memcpy(ty, t1, n * sizeof *t1);
	Zf(poly_mul_fft_custom)(tx, b00, logn);
	Zf(poly_mul_fft_custom)(ty, b10, logn);
	Zf(poly_add_custom)(tx, ty, logn);
	memcpy(ty, t0, n * sizeof *t0);
	Zf(poly_mul_fft_custom)(ty, b01, logn);

	memcpy(t0, tx, n * sizeof *tx);
	Zf(poly_mul_fft_custom)(t1, b11, logn);
	Zf(poly_add_custom)(t1, ty, logn);
	Zf(iFFT_custom)(t0, logn);
	Zf(iFFT_custom)(t1, logn);
	
	int16_t *pt0 = (int16_t *)tx;
	int16_t *pt1 = (int16_t *)tmp;
	vint64m1_t scalar = vmv_v_x_i64m1(0, 1);
	size_t N = n;
    for (size_t vl; N > 0; N -= vl, t0 += vl, t1 += vl, pt0 += vl, pt1 += vl, hm += vl) {
        vl = vsetvl_e16m1(N);
        vuint16m1_t vec_hm = vle16_v_u16m1(hm, vl);
		vuint32m2_t vec_hm_uint32 = vzext_vf2_u32m2(vec_hm, vl);
		vint32m2_t vec_hm_int32 = vreinterpret_v_u32m2_i32m2(vec_hm_uint32);
        vfloat64m4_t vec_t0_f = vle64_v_f64m4((double *)t0, vl);
        vfloat64m4_t vec_t1_f = vle64_v_f64m4((double *)t1, vl);
        vint64m4_t vec_t0_int64 = vfcvt_x_f_v_i64m4(vec_t0_f, vl);
        vint64m4_t vec_t1_int64 = vfcvt_x_f_v_i64m4(vec_t1_f, vl);
		vint32m2_t vec_t0_int32 = vnsra_wx_i32m2(vec_t0_int64, 0, vl);
		vint32m2_t vec_t1_int32 = vnsra_wx_i32m2(vec_t1_int64, 0, vl);
		vint32m2_t vec_sub = vsub_vv_i32m2(vec_hm_int32, vec_t0_int32, vl);
		vint32m2_t vec_rsub = vrsub_vx_i32m2(vec_t1_int32, 0, vl);
		vint16m1_t vec_z = vnsra_wx_i16m1(vec_sub, 0, vl);
		vint16m1_t vec_y = vnsra_wx_i16m1(vec_rsub, 0, vl);
		vse16_v_i16m1(pt0, vec_z, vl);
		vse16_v_i16m1(pt1, vec_y, vl);
		vint32m2_t vec_mul = vwmul_vv_i32m2(vec_z, vec_z, vl);
		scalar = vwredsum_vs_i32m2_i64m1(scalar, vec_mul, scalar, vl);
    }
	uint64_t sqn = vmv_x_s_i64m1_i64(scalar);

	/*
	 * With "normal" degrees (e.g. 512 or 1024), it is very
	 * improbable that the computed vector is not short enough;
	 * however, it may happen in practice for the very reduced
	 * versions (e.g. degree 16 or below). In that case, the caller
	 * will loop, and we must not write anything into s2[] because
	 * s2[] may overlap with the hashed message hm[] and we need
	 * hm[] for the next iteration.
	 */
	s1tmp = (int16_t *)tx;
	s2tmp = (int16_t *)tmp;
	if (Zf(is_short_half_custom)(sqn, s2tmp, logn)) {
		memcpy(s2, s2tmp, n * sizeof *s2);
		memcpy(tmp, s1tmp, n * sizeof *s1tmp);
		return 1;
	}
	return 0;
}

/* see inner.h */
void
Zf(sign_tree_custom_with_rvvf)(int16_t *sig,
	const fpr *restrict expanded_key,
	const uint16_t *hm, unsigned logn, uint8_t *tmp)
{
	fpr *ftmp;

	ftmp = (fpr *)tmp;
	for (;;) {
		/*
		 * Signature produces short vectors s1 and s2. The
		 * signature is acceptable only if the aggregate vector
		 * s1,s2 is short; we must use the same bound as the
		 * verifier.
		 *
		 * If the signature is acceptable, then we return only s2
		 * (the verifier recomputes s1 from s2, the hashed message,
		 * and the public key).
		 */
		samplerZ_custom samp;

		/*
		 * Normal sampling. We use a fast PRNG seeded from our
		 * SHAKE context ('rng').
		 */
		Zf(prng_init_custom)();
		csr_selGaussMode_rw(9);
		samp = Zf(sampler_custom);

		/*
		 * Do the actual signature.
		 */
		if (do_sign_tree_custom_with_rvvf(samp, sig,
			expanded_key, hm, logn, ftmp))
		{
			break;
		}
	}
}

/* see inner.h */
void
Zf(sign_dyn_custom_with_rvvf)(int16_t *sig,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint16_t *hm, unsigned logn, uint8_t *tmp)
{
	fpr *ftmp;

	ftmp = (fpr *)tmp;
	for (;;) {
		/*
		 * Signature produces short vectors s1 and s2. The
		 * signature is acceptable only if the aggregate vector
		 * s1,s2 is short; we must use the same bound as the
		 * verifier.
		 *
		 * If the signature is acceptable, then we return only s2
		 * (the verifier recomputes s1 from s2, the hashed message,
		 * and the public key).
		 */
		samplerZ_custom samp;

		/*
		 * Normal sampling. We use a fast PRNG seeded from our
		 * SHAKE context ('rng').
		 */
		Zf(prng_init_custom)();
		csr_selGaussMode_rw(9);
		samp = Zf(sampler_custom);

		/*
		 * Do the actual signature.
		 */
		if (do_sign_dyn_custom_with_rvvf(samp, sig,
			f, g, F, G, hm, logn, ftmp))
		{
			break;
		}
	}
}

/*
** modified from ffSampling_fft_dyntree_custom
** but without custom samplerz
*/
static void
ffSampling_fft_dyntree_custom_no_cus_spz(samplerZ samp, void *samp_ctx,
	fpr *restrict t0, fpr *restrict t1,
	fpr *restrict g00, fpr *restrict g01, fpr *restrict g11,
	unsigned orig_logn, unsigned logn, fpr *restrict tmp)
{
	size_t n, hn;
	fpr *z0, *z1;

	/*
	 * Deepest level: the LDL tree leaf value is just g00 (the
	 * array has length only 1 at this point); we normalize it
	 * with regards to sigma, then use it for sampling.
	 */
	if (logn == 0) {
		fpr leaf;

		leaf = g00[0];
		leaf = fpr_mul(fpr_sqrt(leaf), fpr_inv_sigma[orig_logn]);
		t0[0] = fpr_of(samp(samp_ctx, t0[0], leaf));
		t1[0] = fpr_of(samp(samp_ctx, t1[0], leaf));
		return;
	}

	n = (size_t)1 << logn;
	hn = n >> 1;

	/*
	 * Decompose G into LDL. We only need d00 (identical to g00),
	 * d11, and l10; we do that in place.
	 */
	Zf(poly_LDL_fft_custom)(g00, g01, g11, logn);

	/*
	 * Split d00 and d11 and expand them into half-size quasi-cyclic
	 * Gram matrices. We also save l10 in tmp[].
	 */
	Zf(poly_split_fft_custom)(tmp, tmp + hn, g00, logn);
	memcpy(g00, tmp, n * sizeof *tmp);
	Zf(poly_split_fft_custom)(tmp, tmp + hn, g11, logn);
	memcpy(g11, tmp, n * sizeof *tmp);
	memcpy(tmp, g01, n * sizeof *g01);
	memcpy(g01, g00, hn * sizeof *g00);
	memcpy(g01 + hn, g11, hn * sizeof *g00);

	/*
	 * The half-size Gram matrices for the recursive LDL tree
	 * building are now:
	 *   - left sub-tree: g00, g00+hn, g01
	 *   - right sub-tree: g11, g11+hn, g01+hn
	 * l10 is in tmp[].
	 */

	/*
	 * We split t1 and use the first recursive call on the two
	 * halves, using the right sub-tree. The result is merged
	 * back into tmp + 2*n.
	 */
	z1 = tmp + n;
	Zf(poly_split_fft_custom)(z1, z1 + hn, t1, logn);
	ffSampling_fft_dyntree_custom_no_cus_spz(samp, samp_ctx, z1, z1 + hn,
		g11, g11 + hn, g01 + hn, orig_logn, logn - 1, z1 + n);
	Zf(poly_merge_fft_custom)(tmp + (n << 1), z1, z1 + hn, logn);

	/*
	 * Compute tb0 = t0 + (t1 - z1) * l10.
	 * At that point, l10 is in tmp, t1 is unmodified, and z1 is
	 * in tmp + (n << 1). The buffer in z1 is free.
	 *
	 * In the end, z1 is written over t1, and tb0 is in t0.
	 */
	memcpy(z1, t1, n * sizeof *t1);
	Zf(poly_sub_custom)(z1, tmp + (n << 1), logn);
	memcpy(t1, tmp + (n << 1), n * sizeof *tmp);
	Zf(poly_mul_fft_custom)(tmp, z1, logn);
	Zf(poly_add_custom)(t0, tmp, logn);

	/*
	 * Second recursive invocation, on the split tb0 (currently in t0)
	 * and the left sub-tree.
	 */
	z0 = tmp;
	Zf(poly_split_fft_custom)(z0, z0 + hn, t0, logn);
	ffSampling_fft_dyntree_custom_no_cus_spz(samp, samp_ctx, z0, z0 + hn,
		g00, g00 + hn, g01, orig_logn, logn - 1, z0 + n);
	Zf(poly_merge_fft_custom)(t0, z0, z0 + hn, logn);
}

/*
** modified from do_sign_dyn_custom
** but without custom samplerz
*/
static int
do_sign_dyn_custom_no_cus_spz(samplerZ samp, void *samp_ctx, int16_t *s2,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint16_t *hm, unsigned logn, fpr *restrict tmp)
{
	size_t n, u;
	fpr *t0, *t1, *tx, *ty;
	fpr *b00, *b01, *b10, *b11, *g00, *g01, *g11;
	fpr ni;
	int16_t *s1tmp, *s2tmp;

	n = MKN(logn);

	/*
	 * Lattice basis is B = [[g, -f], [G, -F]]. We convert it to FFT.
	 */
	b00 = tmp;
	b01 = b00 + n;
	b10 = b01 + n;
	b11 = b10 + n;
	smallints_to_fpr(b01, f, logn);
	smallints_to_fpr(b00, g, logn);
	smallints_to_fpr(b11, F, logn);
	smallints_to_fpr(b10, G, logn);
	Zf(FFT_custom)(b01, logn);
	Zf(FFT_custom)(b00, logn);
	Zf(FFT_custom)(b11, logn);
	Zf(FFT_custom)(b10, logn);
	Zf(poly_neg_custom)(b01, logn);
	Zf(poly_neg_custom)(b11, logn);

	/*
	 * Compute the Gram matrix G = B·B*. Formulas are:
	 *   g00 = b00*adj(b00) + b01*adj(b01)
	 *   g01 = b00*adj(b10) + b01*adj(b11)
	 *   g10 = b10*adj(b00) + b11*adj(b01)
	 *   g11 = b10*adj(b10) + b11*adj(b11)
	 *
	 * For historical reasons, this implementation uses
	 * g00, g01 and g11 (upper triangle). g10 is not kept
	 * since it is equal to adj(g01).
	 *
	 * We _replace_ the matrix B with the Gram matrix, but we
	 * must keep b01 and b11 for computing the target vector.
	 */
	t0 = b11 + n;
	t1 = t0 + n;

	memcpy(t0, b01, n * sizeof *b01);
	Zf(poly_mulselfadj_fft_custom)(t0, logn);    // t0 <- b01*adj(b01)

	memcpy(t1, b00, n * sizeof *b00);
	Zf(poly_muladj_fft_custom)(t1, b10, logn);   // t1 <- b00*adj(b10)
	Zf(poly_mulselfadj_fft_custom)(b00, logn);   // b00 <- b00*adj(b00)
	Zf(poly_add_custom)(b00, t0, logn);      // b00 <- g00
	memcpy(t0, b01, n * sizeof *b01);
	Zf(poly_muladj_fft_custom)(b01, b11, logn);  // b01 <- b01*adj(b11)
	Zf(poly_add_custom)(b01, t1, logn);      // b01 <- g01

	Zf(poly_mulselfadj_fft_custom)(b10, logn);   // b10 <- b10*adj(b10)
	memcpy(t1, b11, n * sizeof *b11);
	Zf(poly_mulselfadj_fft_custom)(t1, logn);    // t1 <- b11*adj(b11)
	Zf(poly_add_custom)(b10, t1, logn);      // b10 <- g11

	/*
	 * We rename variables to make things clearer. The three elements
	 * of the Gram matrix uses the first 3*n slots of tmp[], followed
	 * by b11 and b01 (in that order).
	 */
	g00 = b00;
	g01 = b01;
	g11 = b10;
	b01 = t0;
	t0 = b01 + n;
	t1 = t0 + n;

	/*
	 * Memory layout at that point:
	 *   g00 g01 g11 b11 b01 t0 t1
	 */

	/*
	 * Set the target vector to [hm, 0] (hm is the hashed message).
	 */
	for (u = 0; u < n; u ++) {
		t0[u] = fpr_of(hm[u]);
		/* This is implicit.
		t1[u] = fpr_zero;
		*/
	}

	/*
	 * Apply the lattice basis to obtain the real target
	 * vector (after normalization with regards to modulus).
	 */
	Zf(FFT_custom)(t0, logn);
	ni = fpr_inverse_of_q;
	memcpy(t1, t0, n * sizeof *t0);
	Zf(poly_mul_fft_custom)(t1, b01, logn);
	Zf(poly_mulconst_custom)(t1, fpr_neg(ni), logn);
	Zf(poly_mul_fft_custom)(t0, b11, logn);
	Zf(poly_mulconst_custom)(t0, ni, logn);

	/*
	 * b01 and b11 can be discarded, so we move back (t0,t1).
	 * Memory layout is now:
	 *      g00 g01 g11 t0 t1
	 */
	memcpy(b11, t0, n * 2 * sizeof *t0);
	t0 = g11 + n;
	t1 = t0 + n;

	/*
	 * Apply sampling; result is written over (t0,t1).
	 */
	ffSampling_fft_dyntree_custom_no_cus_spz(samp, samp_ctx,
		t0, t1, g00, g01, g11, logn, logn, t1 + n);

	/*
	 * We arrange the layout back to:
	 *     b00 b01 b10 b11 t0 t1
	 *
	 * We did not conserve the matrix basis, so we must recompute
	 * it now.
	 */
	b00 = tmp;
	b01 = b00 + n;
	b10 = b01 + n;
	b11 = b10 + n;
	memmove(b11 + n, t0, n * 2 * sizeof *t0);
	t0 = b11 + n;
	t1 = t0 + n;
	smallints_to_fpr(b01, f, logn);
	smallints_to_fpr(b00, g, logn);
	smallints_to_fpr(b11, F, logn);
	smallints_to_fpr(b10, G, logn);
	Zf(FFT_custom)(b01, logn);
	Zf(FFT_custom)(b00, logn);
	Zf(FFT_custom)(b11, logn);
	Zf(FFT_custom)(b10, logn);
	Zf(poly_neg_custom)(b01, logn);
	Zf(poly_neg_custom)(b11, logn);
	tx = t1 + n;
	ty = tx + n;

	/*
	 * Get the lattice point corresponding to that tiny vector.
	 */
	memcpy(tx, t0, n * sizeof *t0);
	memcpy(ty, t1, n * sizeof *t1);
	Zf(poly_mul_fft_custom)(tx, b00, logn);
	Zf(poly_mul_fft_custom)(ty, b10, logn);
	Zf(poly_add_custom)(tx, ty, logn);
	memcpy(ty, t0, n * sizeof *t0);
	Zf(poly_mul_fft_custom)(ty, b01, logn);

	memcpy(t0, tx, n * sizeof *tx);
	Zf(poly_mul_fft_custom)(t1, b11, logn);
	Zf(poly_add_custom)(t1, ty, logn);
	Zf(iFFT_custom)(t0, logn);
	Zf(iFFT_custom)(t1, logn);
	
	int16_t *pt0 = (int16_t *)tx;
	int16_t *pt1 = (int16_t *)tmp;
	vint64m1_t scalar = vmv_v_x_i64m1(0, 1);
	size_t N = n;
    for (size_t vl; N > 0; N -= vl, t0 += vl, t1 += vl, pt0 += vl, pt1 += vl, hm += vl) {
        vl = vsetvl_e16m1(N);
        vuint16m1_t vec_hm = vle16_v_u16m1(hm, vl);
		vuint32m2_t vec_hm_uint32 = vzext_vf2_u32m2(vec_hm, vl);
		vint32m2_t vec_hm_int32 = vreinterpret_v_u32m2_i32m2(vec_hm_uint32);
        vfloat64m4_t vec_t0_f = vle64_v_f64m4((double *)t0, vl);
        vfloat64m4_t vec_t1_f = vle64_v_f64m4((double *)t1, vl);
        vint64m4_t vec_t0_int64 = vfcvt_x_f_v_i64m4(vec_t0_f, vl);
        vint64m4_t vec_t1_int64 = vfcvt_x_f_v_i64m4(vec_t1_f, vl);
		vint32m2_t vec_t0_int32 = vnsra_wx_i32m2(vec_t0_int64, 0, vl);
		vint32m2_t vec_t1_int32 = vnsra_wx_i32m2(vec_t1_int64, 0, vl);
		vint32m2_t vec_sub = vsub_vv_i32m2(vec_hm_int32, vec_t0_int32, vl);
		vint32m2_t vec_rsub = vrsub_vx_i32m2(vec_t1_int32, 0, vl);
		vint16m1_t vec_z = vnsra_wx_i16m1(vec_sub, 0, vl);
		vint16m1_t vec_y = vnsra_wx_i16m1(vec_rsub, 0, vl);
		vse16_v_i16m1(pt0, vec_z, vl);
		vse16_v_i16m1(pt1, vec_y, vl);
		vint32m2_t vec_mul = vwmul_vv_i32m2(vec_z, vec_z, vl);
		scalar = vwredsum_vs_i32m2_i64m1(scalar, vec_mul, scalar, vl);
    }
	uint64_t sqn = vmv_x_s_i64m1_i64(scalar);

	/*
	 * With "normal" degrees (e.g. 512 or 1024), it is very
	 * improbable that the computed vector is not short enough;
	 * however, it may happen in practice for the very reduced
	 * versions (e.g. degree 16 or below). In that case, the caller
	 * will loop, and we must not write anything into s2[] because
	 * s2[] may overlap with the hashed message hm[] and we need
	 * hm[] for the next iteration.
	 */
	s1tmp = (int16_t *)tx;
	s2tmp = (int16_t *)tmp;
	if (Zf(is_short_half_custom)(sqn, s2tmp, logn)) {
		memcpy(s2, s2tmp, n * sizeof *s2);
		memcpy(tmp, s1tmp, n * sizeof *s1tmp);
		return 1;
	}
	return 0;
}


static void prng_init_no_cus_spz(prng *p)
{
	/*
	 * To ensure reproducibility for a given seed, we
	 * must enforce little-endian interpretation of
	 * the state words.
	 */
	uint8_t tmp[56];
	uint64_t th, tl;
	int i;

	uint8_t buf[SHAKE256_RATE];
	keccak_squeezeblocks_custom(buf, 1, SHAKE256_RATE, false);

	memcpy(tmp, buf, 56);
	for (i = 0; i < 14; i ++) {
		uint32_t w;

		w = (uint32_t)tmp[(i << 2) + 0]
			| ((uint32_t)tmp[(i << 2) + 1] << 8)
			| ((uint32_t)tmp[(i << 2) + 2] << 16)
			| ((uint32_t)tmp[(i << 2) + 3] << 24);
		*(uint32_t *)(p->state.d + (i << 2)) = w;
	}
	tl = *(uint32_t *)(p->state.d + 48);
	th = *(uint32_t *)(p->state.d + 52);
	*(uint64_t *)(p->state.d + 48) = tl + (th << 32);
	Zf(prng_refill)(p);
}

/* see inner.h */
void
Zf(sign_dyn_custom_no_cus_spz)(int16_t *sig,
	const int8_t *restrict f, const int8_t *restrict g,
	const int8_t *restrict F, const int8_t *restrict G,
	const uint16_t *hm, unsigned logn, uint8_t *tmp)
{
	fpr *ftmp;

	ftmp = (fpr *)tmp;
	for (;;) {
		/*
		 * Signature produces short vectors s1 and s2. The
		 * signature is acceptable only if the aggregate vector
		 * s1,s2 is short; we must use the same bound as the
		 * verifier.
		 *
		 * If the signature is acceptable, then we return only s2
		 * (the verifier recomputes s1 from s2, the hashed message,
		 * and the public key).
		 */
		sampler_context spc;
		samplerZ samp;
		void *samp_ctx;

		/*
		 * Normal sampling. We use a fast PRNG seeded from our
		 * SHAKE context ('rng').
		 */
		spc.sigma_min = fpr_sigma_min[logn];
		prng_init_no_cus_spz(&spc.p);
		samp = Zf(sampler);
		samp_ctx = &spc;

		/*
		 * Do the actual signature.
		 */
		if (do_sign_dyn_custom_no_cus_spz(samp, samp_ctx, sig,
			f, g, F, G, hm, logn, ftmp))
		{
			break;
		}
	}
}