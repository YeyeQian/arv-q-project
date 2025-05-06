#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include "poly_op_util.h"
#include "conversions.h"

void poly_mod_mul(OUT uint8_t res_bin[R_SIZE],
    IN const uint8_t a_bin[R_SIZE],
    IN const uint8_t b_bin[R_SIZE])
{
    uint32_t b_compact[DV] = { 0 };
    convert2compact(b_compact, b_bin);

    T_TYPE a_inSEW[R_SIZE_PRIME + 1] = { 0 };
    General_sparsepoly_pack<T_TYPE>(a_bin, a_inSEW + 1);
    General_sparsepoly_compensate<T_TYPE, U_TYPE>(a_inSEW);

    T_TYPE myres_mul_inSEW[R_SIZE_PRIME] = { 0 };
    General_sparse_poly_mul<T_TYPE, U_TYPE>(a_inSEW, b_compact, myres_mul_inSEW);

    General_sparsepoly_unpack<T_TYPE>(myres_mul_inSEW, res_bin);
}

void cyclic_shift_left_in8(uint8_t res_out[R_SIZE], const uint8_t dense_padded[R_SIZE+1], const uint32_t offset) {
    uint8_t chunk_mask = 255;
    uint32_t chunk_solved = 0; //track the number of chunk already dealt with
    uint32_t chunk_shift = offset >> 3;
    uint32_t bit_shift1 = offset & (8 - 1);
    uint32_t bit_shift2 = bit_shift1 + 8 - (R_BITS & (8 - 1));
    bool flag0 = bit_shift2 > 8;
    //loop1: do stripmining on a_in, each round of loop deal with 16(VLMAX) chunk and there are 386(R_SIZE_) chunks
    while (chunk_solved < R_SIZE) {
        uint32_t len = (chunk_solved + 32) <= R_SIZE ? 32 : (R_SIZE - chunk_solved);
        uint32_t temp = chunk_solved + chunk_shift;
        bool flag1 = temp >= R_SIZE;
        uint32_t bunch_offset = flag1 ? (temp - R_SIZE) : temp;
        //loop2: simulate vector processing on a bunch of chunk
        for (int j = 0; j < len; j++) {
            uint32_t temp = bunch_offset + j;
            bool flag2 = temp >= R_SIZE;
            uint32_t chunk_offset = flag2 ? temp - R_SIZE : temp;
            uint32_t bit_shift = flag1 || flag2 ? (flag0 ? bit_shift2 - 8 : bit_shift2) : bit_shift1;
            int32_t sp_offset = (flag1 || flag2) && flag0 ? -1 : 0;
            uint8_t chunk_shifted = ((((uint16_t)dense_padded[chunk_solved + j + 1 + sp_offset] << 8 | (uint16_t)dense_padded[chunk_solved + j + sp_offset]) << bit_shift) >> 8) & chunk_mask;
            res_out[chunk_offset] ^= chunk_shifted;
        }
        chunk_solved += len;
    }
    uint8_t final_chunk_mask = ((1ULL << (R_BITS & (8 - 1))) - 1);
    res_out[R_SIZE - 1] &= final_chunk_mask;
}

void poly_mod_mul_in8(OUT uint8_t res_bin[R_SIZE],
    IN const uint8_t a_bin[R_SIZE],
    IN const uint8_t b_bin[R_SIZE])
{
    uint32_t b_compact[DV] = { 0 };
    convert2compact(b_compact, b_bin);

    uint8_t a_padded[R_SIZE + 1] = { 0 };
    memcpy(a_padded+1,a_bin,R_SIZE);
    uint8_t additional_chunk = 0;
    uint8_t last_chunk_bits = R_BITS & (8 - 1);
    additional_chunk = (((uint16_t)a_padded[R_SIZE] << 8 | (uint16_t)a_padded[R_SIZE - 1]) >> last_chunk_bits) & 255;
    a_padded[0] = additional_chunk;

    for (int i = 0; i < DV; i++) {
		cyclic_shift_left_in8(res_bin,a_padded,b_compact[i]);
	}
}

void poly_mod_inv(OUT uint8_t res_bin[R_SIZE],
    IN const uint8_t a_bin[R_SIZE])
{
    T_TYPE a_inSEW[R_SIZE_PRIME] = { 0 };
    General_sparsepoly_pack<T_TYPE>(a_bin, a_inSEW);

    T_TYPE myres_inv_inSEW[R_SIZE_PRIME] = { 0 };
    General_poly_inverse<T_TYPE, U_TYPE>(S, a_inSEW, myres_inv_inSEW);

    General_sparsepoly_unpack<T_TYPE>(myres_inv_inSEW, res_bin);
}

void poly_split_polynomial(OUT uint8_t e0[R_SIZE],
    OUT uint8_t e1[R_SIZE],
    IN const uint8_t e[2 * R_SIZE])
{
    uint8_t bits_remained = R_BITS & 7;
    uint8_t mask = (1 << bits_remained) - 1;

    //fill e0
    for (int i = 0; i < R_SIZE; i++)
        e0[i] = e[i];
    e0[R_SIZE - 1] &= mask;

    //fill e1
    for (int i = 0; i < R_SIZE-1; i++) {
        e1[i] = (((uint16_t)e[i + R_SIZE] << 8) | ((uint16_t)e[i + R_SIZE - 1])) >> bits_remained;
    }
    if (bits_remained <= 4) {
        e1[R_SIZE - 1] = e[2 * R_SIZE - 2] >> bits_remained;
    }
    else {
        e1[R_SIZE - 1] = (((uint16_t)e[2 * R_SIZE - 1] << 8) | ((uint16_t)e[2 * R_SIZE - 2])) >> bits_remained;
    }
}

void poly_add(OUT uint8_t res_bin[R_SIZE],
    IN const uint8_t a_bin[R_SIZE],
    IN const uint8_t b_bin[R_SIZE])
{
    for (int i = 0; i < R_SIZE; i++)
        res_bin[i] = a_bin[i] ^ b_bin[i];
}