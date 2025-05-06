#include <riscv_vector.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include "inner.h"
#include "params_custom.h"
#include "../apis/custom_inst_api.h"

#define NONCELEN   40
#define MLEN 32

void gen_rand_int8_array(int8_t* a, int8_t* b, size_t len) {
    srand((unsigned int)time(NULL));
    for (size_t i = 0; i < len; ++i)
    {
        uint8_t num = (rand() & 0xff);
        while(num == 0 || num == 128) {
        	num = (rand() & 0x0ff);
			break;
        }
		a[i] = num - 128;
		b[i] = num - 128;
    }
}

void gen_rand_uint16_array(uint16_t* a, uint16_t* b, size_t len) {
    srand((unsigned int)time(NULL));
    for (size_t i = 0; i < len; ++i)
    {
        uint16_t num = rand() % 12289;
        a[i] = num;
        b[i] = num;
    }
}

void gen_rand_int16_array(int16_t* a, int16_t* b, size_t len) {
    srand((unsigned int)time(NULL));
    for (size_t i = 0; i < len; ++i)
    {
        uint16_t num = (rand() & 0xfff);
        while(num == 0 || num == 2048) {
        	num = (rand() & 0xfff);
			break;
        }
		a[i] = num - 2048;
		b[i] = num - 2048;
    }
}

void gen_rand_small_int16_array(int16_t* a, int16_t* b, size_t len) {
    srand((unsigned int)time(NULL));
    for (size_t i = 0; i < len; ++i)
    {
        uint16_t num = (rand() & 0x0ff);
        while(num == 0 || num == 128) {
        	num = (rand() & 0x0ff);
			break;
        }
		a[i] = num - 128;
		b[i] = num - 128;
    }
}

bool test_hash_to_point_vartime_custom()
{
    uint16_t r_hm_ref[512];
    uint16_t r_hm[512];
    unsigned char nonce[NONCELEN], m[MLEN];

    generate_random_uint8(nonce, sizeof(nonce), NOT_NEED_BOUND);
    generate_random_uint8(m, sizeof(m), NOT_NEED_BOUND);

    inner_shake256_context sc;
	inner_shake256_init(&sc);
	inner_shake256_inject(&sc, nonce, sizeof nonce);
	inner_shake256_inject(&sc, m, MLEN);
	inner_shake256_flip(&sc);
	Zf(hash_to_point_vartime)(&sc, r_hm_ref, 9);

    uint8_t nonce_m[NONCELEN+MLEN];
    memcpy(nonce_m, nonce, NONCELEN);
    memcpy(nonce_m+NONCELEN, m, MLEN);

    csr_keccakmode_rw(SHAKE_256_MODE);
    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    keccak_absorb_custom(nonce_m, NONCELEN + MLEN, SHAKE256_RATE);
    Zf(hash_to_point_vartime_custom)(r_hm, 9);

    bool result = check_same_uint16(r_hm_ref, r_hm, 512);

    if(result) {
        printf("test_hash_to_point_vartime_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_hash_to_point_vartime_custom fail with result=%d..\n", result);
    }

    return result;    
}

bool test_is_short_custom()
{
    int16_t s1_cust[FALCON_N], s1_orig[FALCON_N];
    gen_rand_small_int16_array(s1_cust, s1_orig, FALCON_N);
    int16_t s2_cust[FALCON_N], s2_orig[FALCON_N];
    gen_rand_small_int16_array(s2_cust, s2_orig, FALCON_N);
    int sqn_cust, sqn_orig;

    sqn_orig = Zf(is_short)(s1_orig, s2_orig, LOGN);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    csr_keccakmode_rw(SHAKE_256_MODE);
    sqn_cust = Zf(is_short_custom)(s1_cust, s2_cust, LOGN);

	printf("sqn_cust = %d \n", sqn_cust);
	printf("sqn_orig = %d \n", sqn_orig);

    bool result = (sqn_cust == sqn_orig);

    if(result) {
        printf("test_is_short_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_is_short_custom fail with result=%d..\n", result);
    }

    return result;
}

bool test_is_short_half_custom()
{
    int16_t s2_cust[FALCON_N], s2_orig[FALCON_N];
    gen_rand_small_int16_array(s2_cust, s2_orig, FALCON_N);
    int sqn_cust, sqn_orig;

    sqn_orig = Zf(is_short_half)(9, s2_orig, LOGN);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    csr_keccakmode_rw(SHAKE_256_MODE);
    sqn_cust = Zf(is_short_half_custom)(9, s2_cust, LOGN);

	printf("sqn_cust = %d \n", sqn_cust);
	printf("sqn_orig = %d \n", sqn_orig);

    bool result = (sqn_cust == sqn_orig);

    if(result) {
        printf("test_is_short_half_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_is_short_half_custom fail with result=%d..\n", result);
    }

    return result;
}

int main()
{
    // test_hash_to_point_vartime_custom();

    // test_is_short_custom();
    test_is_short_half_custom();
    return 0;
}