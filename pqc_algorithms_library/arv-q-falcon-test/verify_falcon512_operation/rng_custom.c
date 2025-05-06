#include "inner.h"
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"
#include "params_custom.h"

/* see inner.h */
void
Zf(prng_init_custom)(void)
{
	/*
	 * To ensure reproducibility for a given seed, we
	 * must enforce little-endian interpretation of
	 * the state words.
	 */
	// uint8_t tmp[56];
	// inner_shake256_extract_custom(src, tmp, 56);
	uint8_t tmp[SHAKE256_RATE];
	keccak_squeezeblocks_custom(tmp, 1, SHAKE256_RATE, false);
	size_t vl = vsetvl_e8m2(56);
	vuint8m2_t vec_value = vle8_v_u8m2(tmp, vl);
    chacha20_init_v_u8m2(vec_value);
}