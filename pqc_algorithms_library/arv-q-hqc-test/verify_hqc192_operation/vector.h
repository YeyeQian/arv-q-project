#ifndef VECTOR_H
#define VECTOR_H

/**
 * @file vector.h
 * @brief Header file for vector.cpp
 */

#include "shake_prng.h"
#include <stdint.h>

void vect_set_random_fixed_weight(seedexpander_state* ctx, uint64_t* v, uint16_t weight);
void vect_set_random(seedexpander_state* ctx, uint64_t* v);
void vect_set_random_from_prng(uint64_t* v, uint32_t size_v);

void vect_add(uint64_t* o, const uint64_t* v1, const uint64_t* v2, uint32_t size);
uint8_t vect_compare(const uint8_t* v1, const uint8_t* v2, uint32_t size);
void vect_resize(uint64_t* o, uint32_t size_o, const uint64_t* v, uint32_t size_v);

void vect_print(const uint64_t* v, const uint32_t size);
void vect_print_sparse(const uint32_t* v, const uint16_t weight);


//Re-Organized version
void vect_set_random_fixed_weight_xy(uint64_t* x, uint64_t* y, uint32_t* support, const unsigned char* sk_seed);
void vect_set_random_h(unsigned char* pk_seed, uint64_t* v);
void vect_set_random_fixed_weight_r1r2e(uint64_t* r1, uint64_t* r2, uint32_t* support, uint64_t* e, const unsigned char* theta);

void vect_set_random_fixed_weight_xy_custom(uint64_t* x, uint64_t* y, uint32_t* support, const unsigned char* sk_seed);
void vect_set_random_fixed_weight_r1r2e_custom(uint64_t* r1, uint64_t* r2, uint32_t* support,uint64_t* e, const unsigned char* theta);
void vect_set_random_h_custom(unsigned char* pk_seed, uint64_t* v);
void vect_add_custom(uint64_t *o, const uint64_t *v1, const uint64_t *v2, uint32_t size);
uint8_t vect_compare_custom(const uint8_t *v1, const uint8_t *v2, uint32_t size);
#endif
