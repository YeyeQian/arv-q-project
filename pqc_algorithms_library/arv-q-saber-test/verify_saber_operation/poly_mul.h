#ifndef POLY_MUL_H
#define POLY_MUL_H

#include "SABER_params.h"
#include <stdint.h>

void poly_mul_acc(const uint16_t a[SABER_N], const uint16_t b[SABER_N], uint16_t res[SABER_N]);
void poly_mul_acc_plain(const uint16_t a[SABER_N], const uint16_t b[SABER_N], uint16_t res[SABER_N]);
void poly_mul_acc_bf(uint32_t a[SABER_N], uint32_t b[SABER_N], uint32_t res[SABER_N]);

//Customized version
void poly_mul_acc_bf_custom(uint32_t a[SABER_N], uint32_t b[SABER_N], uint32_t res[SABER_N]);
#endif