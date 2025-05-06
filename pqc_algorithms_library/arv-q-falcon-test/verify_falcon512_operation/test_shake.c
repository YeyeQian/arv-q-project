#include <riscv_vector.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include "inner.h"

#define Np 20

#define SHAKE_256_IN0_LEN 134
#define SHAKE_256_IN1_LEN 200
#define SHAKE_256_OUT0_LEN 134
#define SHAKE_256_OUT1_LEN 200

void gen_rand_uint8_array(uint8_t* a, uint8_t* b, size_t len) {
    srand((unsigned int)time(NULL));
    for (size_t i = 0; i < len; ++i)
    {
        uint8_t num = rand() & 0xff;
        a[i] = num;
        b[i] = num;
    }
}

bool compareub_1d(uint8_t* a, uint8_t* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

bool test_shake_custom()
{
    uint8_t a0[SHAKE_256_IN0_LEN], b0[SHAKE_256_IN0_LEN];
    uint8_t a1[SHAKE_256_IN1_LEN], b1[SHAKE_256_IN1_LEN];
    uint8_t a01[SHAKE_256_IN0_LEN+SHAKE_256_IN1_LEN];
    uint8_t a01_out[SHAKE_256_OUT0_LEN+SHAKE_256_OUT1_LEN];

    uint8_t a0_out[SHAKE_256_OUT0_LEN], b0_out[SHAKE_256_OUT0_LEN];
    uint8_t a1_out[SHAKE_256_OUT1_LEN], b1_out[SHAKE_256_OUT1_LEN];

    gen_rand_uint8_array(a0, b0, SHAKE_256_IN0_LEN);
    gen_rand_uint8_array(a1, b1, SHAKE_256_IN1_LEN);

    memcpy(a01, a0, SHAKE_256_IN0_LEN);
    memcpy(a01+SHAKE_256_IN0_LEN, a1, SHAKE_256_IN1_LEN);

    shake256_custom(a01_out, SHAKE_256_OUT0_LEN+SHAKE_256_OUT1_LEN, a01, SHAKE_256_IN0_LEN+SHAKE_256_IN1_LEN);
    memcpy(a0_out, a01_out, SHAKE_256_OUT0_LEN);
    memcpy(a1_out, a01_out+SHAKE_256_OUT0_LEN, SHAKE_256_OUT1_LEN);

    inner_shake256_context scc;
	inner_shake256_init(&scc);
	inner_shake256_inject(&scc, a0, SHAKE_256_IN0_LEN);
	inner_shake256_inject(&scc, a1, SHAKE_256_IN1_LEN);
	inner_shake256_flip(&scc);
	inner_shake256_extract(&scc, b0_out, SHAKE_256_OUT0_LEN);
	inner_shake256_extract(&scc, b1_out, SHAKE_256_OUT1_LEN);

    bool result0 = check_same_uint8(a0_out, b0_out, SHAKE_256_OUT0_LEN);
    printf("\n");
    bool result1 = check_same_uint8(a1_out, b1_out, SHAKE_256_OUT1_LEN);

    if(result0 && result1) {
        printf("test_shake_custom pass!\n");
    }
    else {
        printf("test_shake_custom fail..\n");
    }

    return (result0 && result1);
}

int main()
{
    test_shake_custom();

    return 0;
}