#include <stdio.h>
#include "api.h"
#include "poly.h"
#include "poly_mul.h"
#include "pack_unpack.h"
#include "cbd.h"
#include "fips202.h"
#include "ntt.h"

void print_poly32(uint32_t poly[SABER_N]) {
	for (int i = 0; i < SABER_N; i++) {
		printf("%d,", poly[i]);
		if ((i + 1) % 16 == 0)printf("\n");
	}
}

void print_poly16(uint16_t poly[SABER_N]) {
	for (int i = 0; i < SABER_N; i++) {
		printf("%d,", poly[i]);
		if ((i + 1) % 16 == 0)printf("\n");
	}
}

void MatrixVectorMul(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose)
{
	int i, j;
	for (i = 0; i < SABER_L; i++)
	{
		for (j = 0; j < SABER_L; j++)
		{
			if (transpose == 1)
			{
				poly_mul_acc(A[j][i], s[j], res[i]);
			}
			else
			{
				poly_mul_acc(A[i][j], s[j], res[i]);
			}	
		}
	}
}

void MatrixVectorMul_plain(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose)
{
	int i, j;
	for (i = 0; i < SABER_L; i++)
	{
		for (j = 0; j < SABER_L; j++)
		{
			if (transpose == 1)
			{
				poly_mul_acc_plain(A[j][i], s[j], res[i]);
			}
			else
			{
				poly_mul_acc_plain(A[i][j], s[j], res[i]);
			}
		}
	}
	for (i = 0; i < SABER_L; i++) {
		for (j = 0; j < SABER_N; j++) {
			res[i][j] &= (8191);
		}
	}
}

void MatrixVectorMul_bf(const uint16_t A[SABER_L][SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_L][SABER_N], int16_t transpose)
{
	uint32_t temp_res[SABER_L][SABER_N]={0};
	uint32_t temp_A[SABER_L][SABER_L][SABER_N];
	uint32_t temp_s[SABER_L][SABER_N];
	for (int i = 0; i < SABER_L; i++) {
		for (int j = 0; j < SABER_L; j++) {
			for (int k = 0; k < SABER_N; k++) {
				temp_A[i][j][k] = (uint32_t)(A[i][j][k]);
			}
			ntt(temp_A[i][j]);
		}
	}
	for (int i = 0; i < SABER_L; i++) {
		for (int j = 0; j < SABER_N; j++) {
			if (((int16_t)s[i][j]) < 0) {
				temp_s[i][j] = (uint32_t)(((int16_t)s[i][j]) + SABER_Q_EXT);
			}
			else {
				temp_s[i][j] = (uint32_t)s[i][j];
			}
		}
		ntt(temp_s[i]);
	}
	int i, j;
	for (i = 0; i < SABER_L; i++)
	{
		for (j = 0; j < SABER_L; j++)
		{
			if (transpose == 1)
			{
				poly_mul_acc_bf(temp_A[j][i], temp_s[j], temp_res[i]);
			}
			else
			{
				poly_mul_acc_bf(temp_A[i][j], temp_s[j], temp_res[i]);
			}
		}
	}

	for (int i = 0; i < SABER_L; i++) {
		invntt_tomont(temp_res[i]);
	}
	for (i = 0; i < SABER_L; i++) {
		for (j = 0; j < SABER_N; j++) {
			if (temp_res[i][j] > 25171457 / 2) {
				res[i][j] = (uint16_t)((temp_res[i][j] - 25171457) & 8191);
			}
			else {
				res[i][j] =(uint16_t)(temp_res[i][j] & 8191);
			}
		}
	}
}

void InnerProd(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N])
{
	int j;
	for (j = 0; j < SABER_L; j++)
	{
		poly_mul_acc(b[j], s[j], res);
	}
}

void InnerProd_bf(const uint16_t b[SABER_L][SABER_N], const uint16_t s[SABER_L][SABER_N], uint16_t res[SABER_N]) {
	int j;
	uint32_t res_temp[SABER_N] = { 0 };
	for (j = 0; j < SABER_L; j++)
	{
		uint32_t a_tmp[SABER_N];
		uint32_t b_tmp[SABER_N];
		for (int k = 0; k < SABER_N; k++) {
			a_tmp[k] = (uint32_t)b[j][k];
			if (((int16_t)s[j][k]) < 0) {
				b_tmp[k] = (uint32_t)(((int16_t)s[j][k]) + SABER_Q_EXT);
			}
			else {
				b_tmp[k] = (uint32_t)s[j][k];
			}
		}
		ntt(a_tmp);
		ntt(b_tmp);
		poly_mul_acc_bf(a_tmp,b_tmp,res_temp);
	}
	invntt_tomont(res_temp);
	for (j = 0; j < SABER_N; j++) {
		if (res_temp[j] > 25171457 / 2) {
			res[j] = (uint16_t)((res_temp[j] - 25171457) & 8191);
		}
		else {
			res[j] = (uint16_t)(res_temp[j] & 8191);
		}
	}
}

void GenMatrix(uint16_t A[SABER_L][SABER_L][SABER_N], const uint8_t seed[SABER_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYVECBYTES];
	int i;

	shake128(buf, sizeof(buf), seed, SABER_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		BS2POLVECq(buf + i * SABER_POLYVECBYTES, A[i]);
	}
}

void GenSecret(uint16_t s[SABER_L][SABER_N], const uint8_t seed[SABER_NOISE_SEEDBYTES])
{
	uint8_t buf[SABER_L * SABER_POLYCOINBYTES];
	size_t i;

	shake128(buf, sizeof(buf), seed, SABER_NOISE_SEEDBYTES);

	for (i = 0; i < SABER_L; i++)
	{
		cbd(s[i], buf + i * SABER_POLYCOINBYTES);
	}
}
