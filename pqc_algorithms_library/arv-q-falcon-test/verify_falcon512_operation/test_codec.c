#include <riscv_vector.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include "inner.h"

#define N       512
#define LOGN    9
#define Np      20

void gen_rand_uint16_array(uint16_t* a, uint16_t* b, size_t len) {
    srand((unsigned int)time(NULL));
    for (size_t i = 0; i < len; ++i)
    {
        uint16_t num = rand() % 12289;
        a[i] = num;
        b[i] = num;
    }
}

bool compareu_1d(uint16_t* a, uint16_t* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

void gen_rand_int16_array(int16_t* a, int16_t* b, size_t len) {
    srand((unsigned int)time(NULL));
    for (size_t i = 0; i < len; ++i)
    {
        uint16_t num = (rand() & 0x07ff);
        if(num != 0) {
            a[i] = num - 1024;
            b[i] = num - 1024;
        }
    }
}

bool comparei_1d(int16_t* a, int16_t* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

void gen_rand_int8_array(int8_t* a, int8_t* b, size_t len) {
    srand((unsigned int)time(NULL));
    for (size_t i = 0; i < len; ++i)
    {
        uint8_t num = (rand() & 0x1f);
        if(num != 0) {
            a[i] = num - 16;
            b[i] = num - 16;
        }
    }
}

bool compareib_1d(int8_t* a, int8_t* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

bool compareub_1d(uint8_t* a, uint8_t* b, size_t len) {
    for (size_t i = 0; i < len; ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

bool test_modq_encode_decode_custom()
{
    size_t cust_len;
    uint8_t custom_packed[N * 3];

    uint16_t au[N];
    uint16_t au_final[N];

    generate_random_uint16(au, N, 12289);

    size_t max_out_len = (N * 14) >> 3;

    cust_len = Zf(modq_encode_custom)(custom_packed, max_out_len, au, LOGN);

    cust_len = Zf(modq_decode_custom)(au_final, LOGN, custom_packed, cust_len);

    bool result_custom = check_same_uint16(au_final, au, N);

    if(result_custom) {
        printf("test_modq_encode_decode_custom pass!\n");
    }
    else {
        printf("test_modq_encode_decode_custom fail..\n");
    }

    return result_custom;
        
}

bool test_trim_i16_encode_decode_custom()
{
    size_t cust_len;
    uint8_t custom_packed[N * 3];

    int16_t ai[N];
    int16_t bi[N];
    int16_t ai_final[N];

    gen_rand_int16_array(ai, bi, N);

    size_t max_out_len = (N * 11) >> 3;

    cust_len = Zf(trim_i16_encode_custom)(custom_packed, max_out_len, ai, LOGN, 11);

    cust_len = Zf(trim_i16_decode_custom)(ai_final, LOGN, 11, custom_packed, cust_len);

    bool result_custom = check_same_int16(ai_final, ai, N);

    if(result_custom) {
        printf("test_trim_i16_encode_decode_custom pass!\n");
    }
    else {
        printf("test_trim_i16_encode_decode_custom fail..\n");
    }

    return result_custom;
        
}

bool test_trim_i8_encode_decode_custom()
{
    size_t cust_len;
    uint8_t custom_packed[N * 3];

    int8_t ai[N];
    int8_t bi[N];
    int8_t ai_final[N];

    gen_rand_int8_array(ai, bi, N);

    size_t max_out_len = (N * 5) >> 3;

    cust_len = Zf(trim_i8_encode_custom)(custom_packed, max_out_len, ai, LOGN, 5);

    cust_len = Zf(trim_i8_decode_custom)(ai_final, LOGN, 5, custom_packed, cust_len);

    bool result_custom = check_same_int8(ai_final, ai, N);

    if(result_custom) {
        printf("test_trim_i8_encode_decode_custom pass!\n");
    }
    else {
        printf("test_trim_i8_encode_decode_custom fail..\n");
    }

    return result_custom;
        
}

int main()
{
    // test_modq_encode_decode_custom();
    // test_trim_i16_encode_decode_custom();
    test_trim_i8_encode_decode_custom();

    return 0;
}