#ifndef __ARCH_RISCV_GAUSSIAN_SAMPLER_HH_
#define __ARCH_RISCV_GAUSSIAN_SAMPLER_HH_

#include <cstddef>
#include <cstdint>
#include <string>

namespace gem5
{

namespace RiscvISA
{

typedef struct {
	union {
		uint8_t d[512]; /* MUST be 512, exactly */
		uint32_t d32[128];
	} buf;
	size_t ptr;
	union {
		uint8_t d[256];
		uint32_t d32[64];
	} state;
	int type;
} prng;

static prng chacha_prng;
static prng *chacha = &chacha_prng;

static uint8_t chacha20_buffer[56];
static uint8_t chacha20_count = 0;

static inline void set_chacha20_count (uint8_t count) {
	chacha20_count = count;
}

static inline uint8_t get_chacha20_count () {
	return chacha20_count;
}

static inline void set_chacha20_buffer (uint8_t count, uint8_t num) {
	chacha20_buffer[count] = num;
}

static inline uint8_t get_chacha20_buf (uint32_t index) {
	return chacha->buf.d[index];
}

/*
 * PRNG based on ChaCha20.
 *
 * State consists in key (32 bytes) then IV (16 bytes) and block counter
 * (8 bytes). Normally, we should not care about local endianness (this
 * is for a PRNG), but for the NIST competition we need reproducible KAT
 * vectors that work across architectures, so we enforce little-endian
 * interpretation where applicable. Moreover, output words are "spread
 * out" over the output buffer with the interleaving pattern that is
 * naturally obtained from the AVX2 implementation that runs eight
 * ChaCha20 instances in parallel.
 *
 * The block counter is XORed into the first 8 bytes of the IV.
 */
static void prng_refill ()
{

	static const uint32_t CW[] = {
		0x61707865, 0x3320646e, 0x79622d32, 0x6b206574
	};

	uint64_t cc;
	size_t u;

	/*
	 * State uses local endianness. Only the output bytes must be
	 * converted to little endian (if used on a big-endian machine).
	 */
	cc = *(uint64_t *)(chacha->state.d + 48);
	// uint64_t cclow = chacha->state.d32[12];
	// uint64_t cchigh = chacha->state.d32[13];
	// cc = (cchigh << 32) | cclow;
	for (u = 0; u < 8; u ++) {
		uint32_t state[16];
		size_t v;
		int i;

		memcpy(&state[0], CW, sizeof CW);
		memcpy(&state[4], chacha->state.d, 48);
		state[14] ^= (uint32_t)cc;
		state[15] ^= (uint32_t)(cc >> 32);
		for (i = 0; i < 10; i ++) {

#define QROUND(a, b, c, d)   do { \
		state[a] += state[b]; \
		state[d] ^= state[a]; \
		state[d] = (state[d] << 16) | (state[d] >> 16); \
		state[c] += state[d]; \
		state[b] ^= state[c]; \
		state[b] = (state[b] << 12) | (state[b] >> 20); \
		state[a] += state[b]; \
		state[d] ^= state[a]; \
		state[d] = (state[d] <<  8) | (state[d] >> 24); \
		state[c] += state[d]; \
		state[b] ^= state[c]; \
		state[b] = (state[b] <<  7) | (state[b] >> 25); \
	} while (0)

			QROUND( 0,  4,  8, 12);
			QROUND( 1,  5,  9, 13);
			QROUND( 2,  6, 10, 14);
			QROUND( 3,  7, 11, 15);
			QROUND( 0,  5, 10, 15);
			QROUND( 1,  6, 11, 12);
			QROUND( 2,  7,  8, 13);
			QROUND( 3,  4,  9, 14);

#undef QROUND

		}

		for (v = 0; v < 4; v ++) {
			state[v] += CW[v];
		}
		for (v = 4; v < 14; v ++) {
			state[v] += ((uint32_t *)chacha->state.d)[v - 4];
			// state[v] += chacha->state.d32[v - 4];
		}
		state[14] += ((uint32_t *)chacha->state.d)[10]
			^ (uint32_t)cc;
		// state[14] += chacha->state.d32[10] ^ (uint32_t)cc;
		state[15] += ((uint32_t *)chacha->state.d)[11]
			^ (uint32_t)(cc >> 32);
		// state[15] += chacha->state.d32[11] ^ (uint32_t)(cc >> 32);
		cc ++;

		/*
		 * We mimic the interleaving that is used in the AVX2
		 * implementation.
		 */
		// this version is the same as the original falcon c code that when left random number size < 8 byte, these left random numbers are not discarded
		// for (v = 0; v < 16; v ++) {
		// 	chacha->buf.d[(u << 2) + (v << 5) + 0] =
		// 		(uint8_t)state[v];
		// 	chacha->buf.d[(u << 2) + (v << 5) + 1] =
		// 		(uint8_t)(state[v] >> 8);
		// 	chacha->buf.d[(u << 2) + (v << 5) + 2] =
		// 		(uint8_t)(state[v] >> 16);
		// 	chacha->buf.d[(u << 2) + (v << 5) + 3] =
		// 		(uint8_t)(state[v] >> 24);
		// }

		// this version is the same as the optimized falcon c code that the left random numbers are not discarded, when left size < 8
		// this version is the same as the actual hardware design
		for (v = 0; v < 16; v ++) {
			chacha->buf.d[(u << 6) + (v << 2) + 0] =
				(uint8_t)state[v];
			chacha->buf.d[(u << 6) + (v << 2) + 1] =
				(uint8_t)(state[v] >> 8);
			chacha->buf.d[(u << 6) + (v << 2) + 2] =
				(uint8_t)(state[v] >> 16);
			chacha->buf.d[(u << 6) + (v << 2) + 3] =
				(uint8_t)(state[v] >> 24);
		}		
	}
	*(uint64_t *)(chacha->state.d + 48) = cc;
	// chacha->state.d32[12] = (uint32_t)cc;
	// chacha->state.d32[13] = (uint32_t)(cc >> 32);

	chacha->ptr = 0;
}

__attribute__((unused)) static void prng_init ()
{
	uint64_t th, tl;
	uint16_t i;

	for (i = 0; i < 14; i ++) {
		uint32_t w;
		// chacha->state.d[(i << 2) + 0] = chacha20_buffer[(i << 2) + 0];
		// chacha->state.d[(i << 2) + 1] = chacha20_buffer[(i << 2) + 1];
		// chacha->state.d[(i << 2) + 2] = chacha20_buffer[(i << 2) + 2];
		// chacha->state.d[(i << 2) + 3] = chacha20_buffer[(i << 2) + 3];
		w = (uint32_t)chacha20_buffer[(i << 2) + 0]
			| ((uint32_t)chacha20_buffer[(i << 2) + 1] << 8)
			| ((uint32_t)chacha20_buffer[(i << 2) + 2] << 16)
			| ((uint32_t)chacha20_buffer[(i << 2) + 3] << 24);
		*(uint32_t *)(chacha->state.d + (i << 2)) = w;
	}
	tl = *(uint32_t *)(chacha->state.d + 48);
	th = *(uint32_t *)(chacha->state.d + 52);
	*(uint64_t *)(chacha->state.d + 48) = tl + (th << 32);
	prng_refill();
}

/*
 * Get a 64-bit random value from a PRNG.
 */
static inline uint64_t prng_get_u64() {
	size_t u;

	// this version is the same as the original falcon c code that when left random number size < 8 byte, these left random numbers are not discarded
	// /*
	//  * If there are less than 9 bytes in the buffer, we refill it.
	//  * This means that we may drop the last few bytes, but this allows
	//  * for faster extraction code. Also, it means that we never leave
	//  * an empty buffer.
	//  */
	// u = chacha->ptr;
	// if (u >= (sizeof chacha->buf.d) - 9) {
	// 	prng_refill();
	// 	u = 0;
	// }
	// chacha->ptr = u + 8;

	// /*
	//  * On systems that use little-endian encoding and allow
	//  * unaligned accesses, we can simply read the data where it is.
	//  */
	// return (uint64_t)chacha->buf.d[u + 0]
	// 	| ((uint64_t)chacha->buf.d[u + 1] << 8)
	// 	| ((uint64_t)chacha->buf.d[u + 2] << 16)
	// 	| ((uint64_t)chacha->buf.d[u + 3] << 24)
	// 	| ((uint64_t)chacha->buf.d[u + 4] << 32)
	// 	| ((uint64_t)chacha->buf.d[u + 5] << 40)
	// 	| ((uint64_t)chacha->buf.d[u + 6] << 48)
	// 	| ((uint64_t)chacha->buf.d[u + 7] << 56);

	// this version is the same as the optimized falcon c code that the left random numbers are not discarded, when left size < 8
	// this version is the same as the actual hardware design
	u = chacha->ptr;
	uint64_t result = 0;
	uint32_t left = (sizeof chacha->buf.d) - u;
	if (left < 8) {
		for (uint8_t i = 0; i < left; i ++)
			result |= ((uint64_t)chacha->buf.d[u + i] << (i << 3));
		prng_refill();
		u = 8 - left;
		for (uint8_t i = 0; i < u; i ++)
			result |= ((uint64_t)chacha->buf.d[i] << ((i + left) << 3));
		chacha->ptr = u;
	} else {
		result = (uint64_t)chacha->buf.d[u + 0]
		| ((uint64_t)chacha->buf.d[u + 1] << 8)
		| ((uint64_t)chacha->buf.d[u + 2] << 16)
		| ((uint64_t)chacha->buf.d[u + 3] << 24)
		| ((uint64_t)chacha->buf.d[u + 4] << 32)
		| ((uint64_t)chacha->buf.d[u + 5] << 40)
		| ((uint64_t)chacha->buf.d[u + 6] << 48)
		| ((uint64_t)chacha->buf.d[u + 7] << 56);
		chacha->ptr = u + 8;
		if (chacha->ptr == sizeof chacha->buf.d) {
			prng_refill();
		}
	}

	/*
	 * On systems that use little-endian encoding and allow
	 * unaligned accesses, we can simply read the data where it is.
	 */
	return (uint64_t)result;
}

/*
 * Get an 8-bit random value from a PRNG.
 */
static inline unsigned prng_get_u8() {
	unsigned v;

	v = chacha->buf.d[chacha->ptr ++];
	if (chacha->ptr == sizeof chacha->buf.d) {
		prng_refill();
	}
	return v;
}

static inline uint64_t fpr_expm_p63(double x, double ccs) {
	/*
	 * Polynomial approximation of exp(-x) is taken from FACCT:
	 *   https://eprint.iacr.org/2018/1234
	 * Specifically, values are extracted from the implementation
	 * referenced from the FACCT article, and available at:
	 *   https://github.com/raykzhao/gaussian
	 * Tests over more than 24 billions of random inputs in the
	 * 0..log(2) range have never shown a deviation larger than
	 * 2^(-50) from the true mathematical value.
	 */


	/*
	 * Normal implementation uses Horner's method, which minimizes
	 * the number of operations.
	 */

	double d, y;

	d = x;
	y = 0.000000002073772366009083061987;
	y = 0.000000025299506379442070029551 - y * d;
	y = 0.000000275607356160477811864927 - y * d;
	y = 0.000002755586350219122514855659 - y * d;
	y = 0.000024801566833585381209939524 - y * d;
	y = 0.000198412739277311890541063977 - y * d;
	y = 0.001388888894063186997887560103 - y * d;
	y = 0.008333333327800835146903501993 - y * d;
	y = 0.041666666666110491190622155955 - y * d;
	y = 0.166666666666984014666397229121 - y * d;
	y = 0.500000000000019206858326015208 - y * d;
	y = 0.999999999999994892974086724280 - y * d;
	y = 1.000000000000000000000000000000 - y * d;
	y *= ccs;
	const double fpr_ptwo63 = 9223372036854775808.0;
	return (uint64_t)(y * fpr_ptwo63);

}

int gaussian0_sampler () {
	static const uint32_t dist[] = {
		10745844u,  3068844u,  3741698u,
		 5559083u,  1580863u,  8248194u,
		 2260429u, 13669192u,  2736639u,
		  708981u,  4421575u, 10046180u,
		  169348u,  7122675u,  4136815u,
		   30538u, 13063405u,  7650655u,
		    4132u, 14505003u,  7826148u,
		     417u, 16768101u, 11363290u,
		      31u,  8444042u,  8086568u,
		       1u, 12844466u,   265321u,
		       0u,  1232676u, 13644283u,
		       0u,    38047u,  9111839u,
		       0u,      870u,  6138264u,
		       0u,       14u, 12545723u,
		       0u,        0u,  3104126u,
		       0u,        0u,    28824u,
		       0u,        0u,      198u,
		       0u,        0u,        1u
	};

	uint32_t v0, v1, v2, hi;
	uint64_t lo;
	size_t u;
	int z;

	/*
	 * Get a random 72-bit value, into three 24-bit limbs v0..v2.
	 */
	lo = prng_get_u64();
	hi = prng_get_u8();
	v0 = (uint32_t)lo & 0xFFFFFF;
	v1 = (uint32_t)(lo >> 24) & 0xFFFFFF;
	v2 = (uint32_t)(lo >> 48) | (hi << 16);

	/*
	 * Sampled value is z, such that v0..v2 is lower than the first
	 * z elements of the table.
	 */
	z = 0;
	for (u = 0; u < (sizeof dist) / sizeof(dist[0]); u += 3) {
		uint32_t w0, w1, w2, cc;

		w0 = dist[u + 2];
		w1 = dist[u + 1];
		w2 = dist[u + 0];
		cc = (v0 - w0) >> 31;
		cc = (v1 - w1 - cc) >> 31;
		cc = (v2 - w2 - cc) >> 31;
		z += (int)cc;
	}
	return z;
}

/*
 * Sample a bit with probability exp(-x) for some x >= 0.
 */
static int BerExp (double x, double ccs)
{
	int s, i;
	double r;
	uint32_t sw, w;
	uint64_t z;
	
	const double fpr_log2 = 0.69314718055994530941723212146;
	const double fpr_inv_log2 = 1.4426950408889634073599246810;

	/*
	 * Reduce x modulo log(2): x = s*log(2) + r, with s an integer,
	 * and 0 <= r < log(2). Since x >= 0, we can use fpr_trunc().
	 */
	s = (int)(x * fpr_inv_log2);
	r = x - s * fpr_log2;

	/*
	 * It may happen (quite rarely) that s >= 64; if sigma = 1.2
	 * (the minimum value for sigma), r = 0 and b = 1, then we get
	 * s >= 64 if the half-Gaussian produced a z >= 13, which happens
	 * with probability about 0.000000000230383991, which is
	 * approximatively equal to 2^(-32). In any case, if s >= 64,
	 * then BerExp will be non-zero with probability less than
	 * 2^(-64), so we can simply saturate s at 63.
	 */
	sw = (uint32_t)s;
	sw ^= (sw ^ 63) & -((63 - sw) >> 31);
	s = (int)sw;

	/*
	 * Compute exp(-r); we know that 0 <= r < log(2) at this point, so
	 * we can use fpr_expm_p63(), which yields a result scaled to 2^63.
	 * We scale it up to 2^64, then right-shift it by s bits because
	 * we really want exp(-x) = 2^(-s)*exp(-r).
	 *
	 * The "-1" operation makes sure that the value fits on 64 bits
	 * (i.e. if r = 0, we may get 2^64, and we prefer 2^64-1 in that
	 * case). The bias is negligible since fpr_expm_p63() only computes
	 * with 51 bits of precision or so.
	 */
	z = ((fpr_expm_p63(r, ccs) << 1) - 1) >> s;

	/*
	 * Sample a bit with probability exp(-x). Since x = s*log(2) + r,
	 * exp(-x) = 2^-s * exp(-r), we compare lazily exp(-x) with the
	 * PRNG output to limit its consumption, the sign of the difference
	 * yields the expected result.
	 */
	i = 64;
	do {
		i -= 8;
		w = prng_get_u8() - ((uint32_t)(z >> i) & 0xFF);
	} while (!w && i > 0);
	return (int)(w >> 31);
}

/*
 * The sampler produces a random integer that follows a discrete Gaussian
 * distribution, centered on mu, and with standard deviation sigma. The
 * provided parameter isigma is equal to 1/sigma.
 *
 * The value of sigma MUST lie between 1 and 2 (i.e. isigma lies between
 * 0.5 and 1); in Falcon, sigma should always be between 1.2 and 1.9.
 */
double sampler (double mu, double isigma, uint16_t logn)
{
	int s;
	double r, dss, ccs;

	const double fpr_inv_2sqrsigma0 = 0.150865048875372721532312163019;

	/*
	 * Center is mu. We compute mu = s + r where s is an integer
	 * and 0 <= r < 1.
	 */
	int64_t floor = (int64_t)mu;
	s = (int)(floor - (mu < (double)floor));
	r = mu - (double)s;

	/*
	 * dss = 1/(2*sigma^2) = 0.5*(isigma^2).
	 */
	dss = isigma * isigma * 0.5;

	// #define FALCONLOG2N 	9
	const double fpr_sigma_min[] = {
		1.2778336969128335860256340575729042,
		1.2982803343442918539708792538826807
	};

	/*
	 * ccs = sigma_min / sigma = sigma_min * isigma.
	 */
	ccs = isigma * fpr_sigma_min[logn - 9];

	/*
	 * We now need to sample on center r.
	 */
	for (;;) {
		int z0, z, b;
		double x;

		/*
		 * Sample z for a Gaussian distribution. Then get a
		 * random bit b to turn the sampling into a bimodal
		 * distribution: if b = 1, we use z+1, otherwise we
		 * use -z. We thus have two situations:
		 *
		 *  - b = 1: z >= 1 and sampled against a Gaussian
		 *    centered on 1.
		 *  - b = 0: z <= 0 and sampled against a Gaussian
		 *    centered on 0.
		 */
		z0 = gaussian0_sampler();
		b = (int)prng_get_u8() & 1;
		z = b + ((b << 1) - 1) * z0;

		/*
		 * Rejection sampling. We want a Gaussian centered on r;
		 * but we sampled against a Gaussian centered on b (0 or
		 * 1). But we know that z is always in the range where
		 * our sampling distribution is greater than the Gaussian
		 * distribution, so rejection works.
		 *
		 * We got z with distribution:
		 *    G(z) = exp(-((z-b)^2)/(2*sigma0^2))
		 * We target distribution:
		 *    S(z) = exp(-((z-r)^2)/(2*sigma^2))
		 * Rejection sampling works by keeping the value z with
		 * probability S(z)/G(z), and starting again otherwise.
		 * This requires S(z) <= G(z), which is the case here.
		 * Thus, we simply need to keep our z with probability:
		 *    P = exp(-x)
		 * where:
		 *    x = ((z-r)^2)/(2*sigma^2) - ((z-b)^2)/(2*sigma0^2)
		 *
		 * Here, we scale up the Bernouilli distribution, which
		 * makes rejection more probable, but makes rejection
		 * rate sufficiently decorrelated from the Gaussian
		 * center and standard deviation that the whole sampler
		 * can be said to be constant-time.
		 */
		double r1 = ((double)z - r);
		x = r1 * r1 * dss;
		// x = fpr_mul(fpr_sqr(fpr_sub(fpr_of(z), r)), dss);
		x = x - ((double)(z0 * z0) * fpr_inv_2sqrsigma0);
		// x = fpr_sub(x, fpr_mul(fpr_of(z0 * z0), fpr_inv_2sqrsigma0));
		if (BerExp(x, ccs)) {
			/*
			 * Rejection sampling was centered on r, but the
			 * actual center is mu = s + r.
			 */
			return double(s + z);
		}
	}
}

}

}

#endif  // __ARCH_RISCV_GAUSSIAN_SAMPLER_HH__