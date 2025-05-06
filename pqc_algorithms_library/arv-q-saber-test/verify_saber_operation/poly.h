#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "SABER_params.h"
void print_poly32(uint32_t poly[SABER_N]);
void print_poly16(uint16_t poly[SABER_N]);
void MatrixVectorMul(const uint16_t a[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose);
void InnerProd(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N]);
void GenMatrix(uint16_t a[SABER_L][SABER_L][SABER_N], const uint8_t seed[SABER_SEEDBYTES]);
void GenSecret(uint16_t s[SABER_L][SABER_N], const uint8_t seed[SABER_NOISE_SEEDBYTES]);
void MatrixVectorMul_plain(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose);
void MatrixVectorMul_bf(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose);
void InnerProd_bf(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N]);

//Customized version
void MatrixVectorMul_bf_custom(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose);
void MatrixVectorMul_bf_InnerPrd_custom(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose, 
    uint32_t b[SABER_L][SABER_N], uint16_t v[SABER_N]);
void InnerProd_bf_custom(uint32_t b[SABER_L][SABER_N], uint32_t s[SABER_L][SABER_N], uint16_t res[SABER_N]);
void GenMatrix_custom(uint16_t A[SABER_L][SABER_L][SABER_N], const uint8_t seed[SABER_SEEDBYTES]);
void GenSecret_custom(uint16_t s[SABER_L][SABER_N], const uint8_t seed[SABER_NOISE_SEEDBYTES]);
#endif
