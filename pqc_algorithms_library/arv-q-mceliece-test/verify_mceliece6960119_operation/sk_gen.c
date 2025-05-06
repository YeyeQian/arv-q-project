/*
  This file is for secret-key generation
*/

#include "sk_gen.h"

#include "controlbits.h"
#include "params.h"
#include "util.h"
#include "gf.h"
#include "crypto_declassify.h"
#include "crypto_uint16.h"
#include <string.h>
#include <stdbool.h>

static inline crypto_uint16 gf_is_zero_declassify(gf t)
{
  crypto_uint16 mask = crypto_uint16_zero_mask(t);
  crypto_declassify(&mask,sizeof mask);
  return mask;
}

/* input: f, element in GF((2^m)^t) */
/* output: out, minimal polynomial of f */
/* return: 0 for success and -1 for failure */
int genpoly_gen(gf *out, gf *f)
{
	int i, j, k, c;

	gf mat[ SYS_T+1 ][ SYS_T ];
	gf mask, inv, t;

	// fill matrix

	mat[0][0] = 1;

	for (i = 1; i < SYS_T; i++)
		mat[0][i] = 0;

	for (i = 0; i < SYS_T; i++)
		mat[1][i] = f[i];

	for (j = 2; j <= SYS_T; j++)
		GF_mul(mat[j], mat[j-1], f);

	// gaussian //Gauss-Jordan Elimination (Column Major)

	for (j = 0; j < SYS_T; j++)
	{
		for (k = j + 1; k < SYS_T; k++)
		{
			mask = gf_iszero(mat[ j ][ j ]);

			for (c = j; c < SYS_T + 1; c++)
				mat[ c ][ j ] ^= mat[ c ][ k ] & mask;

		}

		if ( gf_is_zero_declassify(mat[ j ][ j ]) ) // return if not systematic
		{
			return -1;
		}

		inv = gf_inv(mat[j][j]);

		for (c = j; c < SYS_T + 1; c++)
			mat[ c ][ j ] = gf_mul(mat[ c ][ j ], inv) ;

		for (c = j; c < SYS_T + 1; c++) {
			for (k = 0; k < SYS_T; k++)
			{
				if (k != j)
				{
					t = mat[j][k];
					mat[c][k] ^= gf_mul(mat[c][j], t);
				}
			}
		}
	}

	for (i = 0; i < SYS_T; i++)
		out[i] = mat[ SYS_T ][ i ];

	return 0;
}

//Vectorized version from RePQC
//#define V 32
//#define TOTAL_PHASE CEIL_DIVIDE(SYS_T+1,V)
//int genpoly_gen_vectorized(gf* out, gf* f)			//reduce in raw major
//{
//	int i, j, k, c;
//
//	gf mat[SYS_T + 1][SYS_T];
//	gf mask, inv, t;
//
//	// fill matrix
//
//	mat[0][0] = 1;
//
//	for (i = 1; i < SYS_T; i++)
//		mat[0][i] = 0;
//
//	for (i = 0; i < SYS_T; i++)
//		mat[1][i] = f[i];
//
//	for (j = 2; j <= SYS_T; j++)
//		GF_mul(mat[j], mat[j - 1], f);
//
//	//Gauss-Jordan Elimination
//	uint8_t order_list[V] = { 0 };
//	gf scale_list[SYS_T + 1][V];
//	for (int phase = 0; phase <= TOTAL_PHASE; phase++) {
//		uint8_t phase_addr = phase * V;
//		memset(order_list,0xFF,sizeof(order_list));
//		memset(scale_list,0,sizeof(gf)*(SYS_T+1)*V);
//		for (int step = phase_addr; step < SYS_T; step += V) {
//			uint8_t columnVec_len = MIN(V,SYS_T-step);
//			uint8_t columnVec_idx[V] = { 0 };
//			for (int i = 0; i < columnVec_len; i++)columnVec_idx[i] = step + i;
//			for (int row = phase_addr; row <= (phase + 1) * V - 1; row++) {
//				uint8_t pivotloc = row - phase * V;
//				gf Rrow[V];
//				memset(Rrow, 0, sizeof(Rrow));
//				for (int elirow = row; elirow <= SYS_T + row; elirow++) {
//					uint8_t eli = elirow % (SYS_T + 1);
//					gf A_in[V];
//					for (int i = 0; i < columnVec_len; i++)A_in[i] = mat[eli][columnVec_idx[i]];
//					gf pivot = A_in[pivotloc];
//					if (step == phase_addr) {//the first step of each pahse
//						if (order_list[row - phase_addr]==0xFF && pivot) {
//							order_list[row - phase_addr] = eli;
//							gf inv_pivot = gf_inv(pivot);
//							scale_list[eli][pivotloc] = inv_pivot;
//							for (int i = 0; i < columnVec_len;i++) {
//								Rrow[i] = gf_mul(A_in[i], inv_pivot);
//								mat[eli][columnVec_idx[i]] = Rrow[i];
//							}
//						}
//						else {
//							gf scale = A_in[pivotloc];
//							scale_list[eli][pivotloc] = scale;
//							for (int i = 0; i < columnVec_len; i++) {
//								mat[eli][columnVec_idx[i]] = gf_add(A_in[i],gf_mul(Rrow[i],scale));
//							}
//						}
//					}
//					else {
//						if (order_list[row - phase_addr] == eli) {
//							gf scale = scale_list[eli][pivotloc];
//							for (int i = 0; i < columnVec_len; i++) {
//								Rrow[i] = gf_mul(A_in[i], scale);
//								mat[eli][columnVec_idx[i]] = Rrow[i];
//							}
//						}
//						else {
//							gf scale = scale_list[eli][pivotloc];
//							for (int i = 0; i < columnVec_len; i++) {
//								mat[eli][columnVec_idx[i]] = gf_add(A_in[i], gf_mul(Rrow[i], scale));
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//	
//
//	for (i = 0; i < SYS_T; i++)
//		out[i] = mat[SYS_T][i];
//
//	return 0;
//}

#define V 32
#define TOTAL_PHASE CEIL_DIVIDE(SYS_T,V)
int genpoly_gen_vectorized(gf* out, gf* f)			//reduce in column major
{
	int i, j, k, c;

	gf mat[SYS_T + 1][SYS_T];
	gf mask, inv, t;

	// fill matrix

	mat[0][0] = 1;

	for (i = 1; i < SYS_T; i++)
		mat[0][i] = 0;

	for (i = 0; i < SYS_T; i++)
		mat[1][i] = f[i];

	for (j = 2; j <= SYS_T; j++)
		GF_mul(mat[j], mat[j - 1], f);

	//Gauss-Jordan Elimination
	uint8_t order_list[V] = { 0 };
	gf scale_list[V][SYS_T];
	for (int phase = 0; phase < TOTAL_PHASE; phase++) {
		uint8_t phase_addr = phase * V;
		memset(order_list,0xFF,sizeof(order_list));
		memset(scale_list,0,sizeof(gf)*SYS_T*V);
		for (int step = phase_addr; step < SYS_T+1; step += V) {
			uint8_t rowVec_len = MIN(V,SYS_T+1-step);						//rowVec means at same column but multiple rows
			uint8_t rowVec_idx[V] = { 0 };
			for (int i = 0; i < rowVec_len; i++)rowVec_idx[i] = step + i;
			for (int col = phase_addr; col <= (phase + 1) * V - 1; col++) {
				uint8_t pivotloc = col - phase * V;
				gf Ccol[V];
				memset(Ccol, 0, sizeof(Ccol));
				for (int elicol = col; elicol < SYS_T + col; elicol++) {
					uint8_t eli = elicol % SYS_T;
					gf A_in[V];
					for (int i = 0; i < rowVec_len; i++)A_in[i] = mat[rowVec_idx[i]][eli];
					gf pivot = A_in[pivotloc];
					if (step == phase_addr) {//the first step of each pahse
						if (order_list[col - phase_addr]==0xFF && pivot) {
							order_list[col - phase_addr] = eli;
							gf inv_pivot = gf_inv(pivot);
							scale_list[pivotloc][eli] = inv_pivot;
							for (int i = 0; i < rowVec_len;i++) {
								Ccol[i] = gf_mul(A_in[i], inv_pivot);
								mat[rowVec_idx[i]][eli] = Ccol[i];
							}
						}
						else {
							gf scale = A_in[pivotloc];
							scale_list[pivotloc][eli] = scale;
							for (int i = 0; i < rowVec_len; i++) {
								mat[rowVec_idx[i]][eli] = gf_add(A_in[i],gf_mul(Ccol[i],scale));
							}
						}
					}
					else {
						if (order_list[col - phase_addr] == eli) {
							gf scale = scale_list[pivotloc][eli];
							for (int i = 0; i < rowVec_len; i++) {
								Ccol[i] = gf_mul(A_in[i], scale);
								mat[rowVec_idx[i]][eli] = Ccol[i];
							}
						}
						else {
							gf scale = scale_list[pivotloc][eli];
							for (int i = 0; i < rowVec_len; i++) {
								mat[rowVec_idx[i]][eli] = gf_add(A_in[i], gf_mul(Ccol[i], scale));
							}
						}
					}
				}
			}
		}
	}
	

	for (i = 0; i < SYS_T; i++)
		out[i] = mat[SYS_T][i];

	return 0;
}