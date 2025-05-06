#include <riscv_vector.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include "inner.h"
#include "params_custom.h"
#include "../apis/custom_inst_api.h"

#define Np 20
#define N_SAMPLED 1

typedef int (*samplerZ)(void *ctx, fpr mu, fpr sigma);

void gen_rand_uint8_array(uint8_t* a, uint8_t* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
    {
        // srand(0xdeadbee8);
        uint8_t num = rand();
        a[i] = num;
        b[i] = num;
    }
}

bool comparehu_1d(uint8_t* a, uint8_t* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

void gen_rand_fpr_mu_array(fpr* a, fpr* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
    {
        // srand(0xdeadbee8);
        double random = ((double) rand()) /(double) RAND_MAX;
        double range = 150.9428f - 10.5231f;
        double num = (random * range) + 10.5231f;
        a[i].v = num;
        b[i].v = num;
    }
    for (size_t i = 0; i < len; i += 2) {
        a[i].v = -a[i].v;
        b[i].v = -b[i].v;
    }
}

void gen_rand_fpr_isigma_array(fpr* a, fpr* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
    {
        // srand(0xdeadbee8);
        double random = ((double) rand()) /(double) RAND_MAX;
        double range = 1.0f - 0.5f;
        double num = (random * range) + 0.5f;
        a[i].v = num;
        b[i].v = num;
    }
    // for (size_t i = 0; i < len; i += 2) {
    //     a[i].v = -a[i].v;
    //     b[i].v = -b[i].v;
    // }
}

void gen_fixed_fpr_mu_array(fpr* a, fpr* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
    {
        double range = 150.9428f - 10.5231f;
        double num = 85.0f;
        a[i].v = num;
        b[i].v = num;
    }
}

void gen_fixed_fpr_isigma_array(fpr* a, fpr* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
    {
        double range = 1.0f - 0.5f;
        double num = 0.55f;
        a[i].v = num;
        b[i].v = num;
    }
}

bool comparefpr_1d(fpr* a, fpr* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
        if (a[i].v != b[i].v)
            return false;
    return true;
}

bool check_same_fpr(fpr* r_ref, fpr* r, size_t length)
{
	size_t i;
	bool same = true;
	for (i = 0; i < length; i++) {
		if (r_ref[i].v != r[i].v) {
			printf("pos=%ld, r_ref=%d, while r=%d\n", i, (int)r_ref[i].v, (int)r[i].v);
			same = false;
		}
	}
	if (same == true) {
		printf("check pass, r_ref == r!\n");
	}
	else {
		printf("check failed, r_ref != r!\n");
	}

	return same;
}

/* see inner.h */
void
Zf(prng_init_for_test)(prng *p, uint8_t* buf)
{
	/*
	 * To ensure reproducibility for a given seed, we
	 * must enforce little-endian interpretation of
	 * the state words.
	 */
	uint8_t tmp[56];
	uint64_t th, tl;
	int i;

    memcpy(tmp, buf, sizeof(tmp));

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
Zf(prng_init_custom_for_test)(uint8_t* buf)
{
	/*
	 * To ensure reproducibility for a given seed, we
	 * must enforce little-endian interpretation of
	 * the state words.
	 */
	size_t vl = vsetvl_e8m2(56);
	vuint8m2_t vec_value = vle8_v_u8m2(buf, vl);
    chacha20_init_v_u8m2(vec_value);
}

bool test_samplerZ_custom()
{
    // prepare mu and isigma
	fpr mu1[FALCON_N], mu2[FALCON_N];
	fpr isigma1[FALCON_N],isigma2[FALCON_N];

	gen_rand_fpr_mu_array(mu1, mu2, FALCON_N);
	gen_rand_fpr_isigma_array(isigma1, isigma2, FALCON_N);

    // guass is sdandard output, testg is custom output
	fpr testg[FALCON_N], guass[FALCON_N];

    // prepare random seed
    unsigned char seed[48];
	generate_random_uint8(seed, sizeof seed, NOT_NEED_BOUND);

    /*
    ** standard process of samplerz
    */
    // in nist.c: crypto_sign
    inner_shake256_context sc;
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, seed, sizeof seed);
	inner_shake256_flip(&sc);

    // in sign.c: sign_dyn
    inner_shake256_context *rng = &sc;
	sampler_context spc;
	samplerZ samp;
	void *samp_ctx;

    spc.sigma_min = fpr_sigma_min[LOGN];

	Zf(prng_init)(&spc.p, rng);
	samp = Zf(sampler);
	samp_ctx = &spc;

	for(uint32_t i = 0; i < FALCON_N; i++)
		guass[i] = fpr_of(samp(samp_ctx, mu1[i], isigma1[i]));

    /*
    ** custom process of samplerz
    */
    csr_keccakmode_rw(SHAKE_256_MODE);
    keccak_absorb_custom(seed, sizeof seed, SHAKE256_RATE);
    csr_selGaussMode_rw(9);
    Zf(prng_init_custom)();
	for(uint32_t i = 0; i < FALCON_N; i++)
		testg[i] = Zf(sampler_custom)(mu2[i], isigma2[i]);

    for(uint32_t i = 0; i < Np; i++) {
        printf("guass[%d] = %27.26f \n", i, guass[i].v);
        printf("testg[%d] = %27.26f \n", i, testg[i].v);
    }

    puts(comparefpr_1d(guass, testg, Np) ? "pass " : "fail ");

    return true;    
}


bool samplerz_debug_test()
{
    // prepare mu and isigma
	fpr mu1[N_SAMPLED], mu2[N_SAMPLED];
	fpr isigma1[N_SAMPLED],isigma2[N_SAMPLED];

	gen_fixed_fpr_mu_array(mu1, mu2, N_SAMPLED);
	gen_fixed_fpr_isigma_array(isigma1, isigma2, N_SAMPLED);

    // guass is sdandard output, testg is custom output
	fpr testg[N_SAMPLED], guass[N_SAMPLED];

    // prepare random seed
    unsigned char seed[56];
	generate_inorder_uint8(seed, sizeof(seed), NOT_NEED_BOUND);

    /*
    ** standard process of samplerz
    */
    // in sign.c: sign_dyn
	sampler_context spc;
	samplerZ samp;
	void *samp_ctx;

    spc.sigma_min = fpr_sigma_min[LOGN];

	Zf(prng_init_for_test)(&spc.p, seed);
	samp = Zf(sampler);
	samp_ctx = &spc;

	for(uint32_t i = 0; i < N_SAMPLED; i++)
		guass[i] = fpr_of(samp(samp_ctx, mu1[i], isigma1[i]));

    /*
    ** custom process of samplerz
    */
    Zf(prng_init_custom_for_test)(seed);
    csr_selGaussMode_rw(9);
	for(uint32_t i = 0; i < N_SAMPLED; i++)
		testg[i] = Zf(sampler_custom)(mu2[i], isigma2[i]);

    for(uint32_t i = 0; i < Np; i++) {
        printf("guass[%d] = %.3f \n", i, guass[i].v);
        printf("testg[%d] = %.3f \n", i, testg[i].v);
    }

    bool result = check_same_fpr(guass, testg, N_SAMPLED);

    if(result) {
        printf("samplerz_debug_test pass!\n");
    }
    else {
        printf("samplerz_debug_test failed..\n");
    }

    return true;    
}

int main()
{
    samplerz_debug_test();

    return 0;
}