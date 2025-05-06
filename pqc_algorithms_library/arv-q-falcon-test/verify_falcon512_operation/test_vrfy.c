#include <riscv_vector.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include "inner.h"
#include "params_custom.h"
#include "../apis/custom_inst_api.h"

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

bool test_mq_NTT_custom()
{
    uint16_t r_ref[512];
    uint16_t r[512];

    generate_inorder_uint16(r_ref, 512, Q);
    generate_inorder_uint16(r, 512, Q);

    mq_NTT(r_ref, LOGN);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    mq_NTT_custom(r, LOGN);

    bool result = check_same_uint16(r_ref, r, 512);

    if(result) {
        printf("test_mq_NTT_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_mq_NTT_custom fail with result=%d!\n", result);
    }

    return result;
}

bool test_mq_iNTT_custom()
{
    uint16_t r_ref[512];
    uint16_t r[512];

    generate_inorder_uint16(r_ref, 512, Q);
    generate_inorder_uint16(r, 512, Q);

    mq_NTT(r_ref, LOGN);
    mq_iNTT(r_ref, LOGN);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    mq_NTT_custom(r, LOGN);
    mq_iNTT_custom(r, LOGN);

    bool result = check_same_uint16(r_ref, r, 512);

    if(result) {
        printf("test_mq_iNTT_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_mq_iNTT_custom fail with result=%d!\n", result);
    }

    return result;
}

bool test_mq_NTT_custom_asm()
{
    uint16_t r_ref[512];
    uint16_t r[512];

    generate_inorder_uint16(r_ref, 512, Q);
    generate_inorder_uint16(r, 512, Q);

    ntt_CG_512(r_ref);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    mq_NTT_custom_asm(r, LOGN);

    bool result = check_same_uint16(r_ref, r, 512);

    if(result) {
        printf("test_mq_NTT_custom_asm pass with result=%d!\n", result);
    }
    else {
        printf("test_mq_NTT_custom_asm fail with result=%d!\n", result);
    }

    return result;
}

bool test_mq_iNTT_custom_asm()
{
    uint16_t r_ref[512];
    uint16_t r[512];

    generate_inorder_uint16(r_ref, 512, Q);
    generate_inorder_uint16(r, 512, Q);

    ntt_CG_512(r_ref);
    bitreverse_standard_pos_transfer(r_ref);
    invntt_CG_512_in_order_zetas_inv(r_ref);
    bitreverse_standard_pos_transfer(r_ref);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    mq_NTT_custom(r, LOGN);
    mq_iNTT_custom_asm(r, LOGN);
    // print_array_uint16(r, 512);
    bool result = check_same_uint16(r_ref, r, 512);

    if(result) {
        printf("test_mq_iNTT_custom_asm pass with result=%d!\n", result);
    }
    else {
        printf("test_mq_iNTT_custom_asm fail with result=%d!\n", result);
    }

    return result;
}

void get_tree_byteoffset()
{
    int i;
    uint16_t tree_byteoffset[FALCON_N];
    for(i = 0; i < FALCON_N; i++) {
        tree_byteoffset[i] = reverse_bit(i, LOGN) * sizeof(uint16_t);
    }
    print_array_uint16(tree_byteoffset, FALCON_N);
}

bool test_to_ntt_monty_custom()
{
    uint16_t r_ref[512];
    uint16_t r[512];

    generate_inorder_uint16(r_ref, 512, Q);
    generate_inorder_uint16(r, 512, Q);

    Zf(to_ntt_monty)(r_ref, LOGN);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    Zf(to_ntt_monty_custom)(r, LOGN);

    bool result = check_same_uint16(r_ref, r, 512);

    if(result) {
        printf("test_to_ntt_monty_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_to_ntt_monty_custom fail with result=%d!\n", result);
    }

    return result;
}

bool test_verify_raw_custom()
{
    int16_t s1_cust[FALCON_N], s1_orig[FALCON_N];
    
    int16_t s2_cust[FALCON_N], s2_orig[FALCON_N];
    gen_rand_small_int16_array(s2_cust, s2_orig, FALCON_N);

    int sqn_cust, sqn_orig;

    uint16_t c[FALCON_N], d[FALCON_N];
    gen_rand_uint16_array(c, d, FALCON_N);

    uint16_t a[FALCON_N], b[FALCON_N];
    gen_rand_uint16_array(a, b, FALCON_N);

    sqn_orig = Zf(verify_raw)(d, s2_orig, b, LOGN, (uint8_t *)s1_orig);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    
    sqn_cust = Zf(verify_raw_custom)(c, s2_cust, a, LOGN, (uint8_t *)s1_cust);

	printf("sqn_cust = %d \n", sqn_cust);
	printf("sqn_orig = %d \n", sqn_orig);

    bool result = check_same_int16(s1_orig, s1_cust, 512);

    if(result) {
        printf("test_verify_raw_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_verify_raw_custom fail with result=%d!\n", result);
    }

    return result;
}

bool test_compute_public_custom()
{
    int16_t s1_cust[FALCON_N], s1_orig[FALCON_N];
    gen_rand_small_int16_array(s1_cust, s1_orig, FALCON_N);
    
    int8_t f_cust[FALCON_N], f_orig[FALCON_N];
    gen_rand_int8_array(f_cust, f_orig, FALCON_N);

    int8_t g_cust[FALCON_N], g_orig[FALCON_N];
    gen_rand_int8_array(g_cust, g_orig, FALCON_N);

    uint16_t c[FALCON_N], d[FALCON_N];

    int sqn_cust, sqn_orig;

    sqn_orig = Zf(compute_public)(d, f_orig, g_orig, LOGN, (uint8_t *)s1_orig);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    
    sqn_cust = Zf(compute_public_custom)(c, f_cust, g_cust, LOGN, (uint8_t *)s1_cust);

	printf("sqn_cust = %d \n", sqn_cust);
	printf("sqn_orig = %d \n", sqn_orig);

    bool result = check_same_uint16(d, c, 512);

    if(result) {
        printf("test_compute_public_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_compute_public_custom fail with result=%d!\n", result);
    }

    return result;
}


bool test_complete_private_custom()
{
    int8_t f_cust[FALCON_N], f_orig[FALCON_N];
    gen_rand_int8_array(f_cust, f_orig, FALCON_N);
    int8_t g_cust[FALCON_N], g_orig[FALCON_N];
    gen_rand_int8_array(g_cust, g_orig, FALCON_N);

    int8_t G_cust[FALCON_N], G_orig[FALCON_N];
    int8_t F_cust[FALCON_N], F_orig[FALCON_N];
    gen_rand_int8_array(F_cust, F_orig, FALCON_N);
	uint8_t tmp_cust[FALCON_N<<2], tmp_orig[FALCON_N<<2];

    int sqn_cust, sqn_orig;

    sqn_orig = Zf(complete_private)(G_orig, f_orig, g_orig, F_orig, LOGN, tmp_orig);

    csr_modulusq_rw(Q);
    csr_qinv_rw(Q0I);
    
    sqn_cust = Zf(complete_private_custom)(G_cust, f_cust, g_cust, F_cust, LOGN, tmp_cust);

	printf("sqn_cust = %d \n", sqn_cust);
	printf("sqn_orig = %d \n", sqn_orig);

    bool result = check_same_int8(G_orig, G_cust, 512);

    if(result) {
        printf("test_complete_private_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_complete_private_custom fail with result=%d!\n", result);
    }

    return result;
}

int main()
{
    // test_mq_NTT_custom();
    // test_mq_iNTT_custom();

    // test_mq_NTT_custom_asm();

    // test_to_ntt_monty_custom();
    // test_verify_raw_custom();

    // test_compute_public_custom();
    // test_complete_private_custom();

    test_mq_iNTT_custom_asm();
    return 0;
}