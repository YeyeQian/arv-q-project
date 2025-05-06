#include <stdlib.h>
#include "extGCD_polyinv.h"

uint64_t clmul_sew(uint32_t a, uint32_t b) {
	uint64_t temp[32] = { 0 };
	for (int i = 0; i < 32; i++) {
		if ((b & 1) == 1) {
			temp[i] = (uint64_t)a << i;
		}
		b >>= 1;
	}
	uint64_t res = 0;
	for (int i = 0; i < 32; i++) res ^= temp[i];
	return res;
}

//d_in of length len and dout of length len+1
void clmul_array(uint32_t* d_in, int len, uint32_t h, uint32_t* d_out) {
	uint64_t temp0 = clmul_sew(d_in[0], h);
	d_out[0] = temp0 & CHUNK_MASK;
	uint32_t temp1 = (temp0 >> SEW) & CHUNK_MASK;
	for (int i = 1; i < len; i++) {
		temp0 = clmul_sew(d_in[i], h);
		d_out[i] = (temp0 & CHUNK_MASK) ^ temp1;
		temp1 = (temp0 >> SEW) & CHUNK_MASK;
	}
	d_out[len] = temp1;
}

void xor_array(uint32_t* din1, uint32_t* din2, uint32_t* dout, int len) {
	for (int i = 0; i < len; i++)dout[i] = din1[i] ^ din2[i];
}

//should know that for if hmatrix's elements are 32bit width, the step_size should be 31(the largest degree number)
///hx[0]-hx[4] stands for h00,h01,h10,h11 in 2x2 matrix
void gen_hmatrix(int32_t* delta_ptr, uint32_t fm, uint32_t gm,uint32_t step_size,uint32_t h[4]) {
	for (int j = 1; j <= step_size; j++) {
		uint32_t g_msb = gm >> (SEW - 1);
		if (g_msb == 0) {
			(*delta_ptr) += 1;
			gm <<= 1;
			h[2] <<= 1, h[3] <<= 1;
		}
		else {
			if ((*delta_ptr) > 0) {
				(*delta_ptr) = -(*delta_ptr) + 1;
				uint32_t temp = fm;
				fm = gm, gm = (temp ^ gm) << 1;
				uint32_t h0_temp = h[0];
				uint32_t h1_temp = h[1];
				h[0] = h[2], h[1] = h[3], h[2] = (h[2] ^ h0_temp) << 1, h[3] = (h[3] ^ h1_temp) << 1;
			}
			else {
				(*delta_ptr) += 1;
				gm = (gm ^ fm) << 1;
				h[2] = (h[2] ^ h[0]) << 1, h[3] = (h[3] ^ h[1]) << 1;
			}
		}
	}
}

void trans(uint32_t *in, uint64_t* out,int len) {
	for (int i = 0; i < len; i += 2) out[i >> 1] = ((uint64_t)in[i + 1] << 32) | in[i];
}

void extGCD_polyinv(uint32_t gin[R_SIZE_PRIME], uint32_t gout[R_SIZE_PRIME]) {
	uint32_t f[R_SIZE_PRIME] = { 0 };
	f[0] = 1 << (SEW - 1 - (R_BITS & (SEW - 1))), f[R_SIZE_PRIME - 1] = 1 << (SEW - 1);
	uint32_t g[R_SIZE_PRIME] = { 0 };
	for (int i = R_SIZE_PRIME - 1; i > 0; i--) {
		g[i] = (((uint64_t)gin[i] << SEW) | (uint64_t)gin[i - 1]) >> (1 + (R_BITS & (SEW - 1)));
	}
	g[0] = gin[0]<<(SEW - 1 - (R_BITS & (SEW - 1)));

	uint32_t w[R_SIZE_PRIME<<1] = { 0 };
	w[0] = 1;
	uint32_t v[R_SIZE_PRIME<<1] = { 0 };

	int32_t delta = 0;
	uint32_t Tau = 2 * R_BITS - 1;

	uint32_t step_size = SEW - 1;

	//some stack spaces
	uint32_t* h0f = (uint32_t*)calloc(R_SIZE_PRIME + 1, sizeof(uint32_t));
	uint32_t* h1g = (uint32_t*)calloc(R_SIZE_PRIME + 1, sizeof(uint32_t));
	uint32_t* h2f = (uint32_t*)calloc(R_SIZE_PRIME + 1, sizeof(uint32_t));
	uint32_t* h3g = (uint32_t*)calloc(R_SIZE_PRIME + 1, sizeof(uint32_t));
	uint32_t* temp_f = (uint32_t*)calloc(R_SIZE_PRIME + 1, sizeof(uint32_t));
	uint32_t* temp_g = (uint32_t*)calloc(R_SIZE_PRIME + 1, sizeof(uint32_t));
	uint32_t* h0v = (uint32_t*)calloc(R_SIZE_PRIME * 2 + 1, sizeof(uint32_t));
	uint32_t* h1w = (uint32_t*)calloc(R_SIZE_PRIME * 2 + 1, sizeof(uint32_t));
	uint32_t* h2v = (uint32_t*)calloc(R_SIZE_PRIME * 2 + 1, sizeof(uint32_t));
	uint32_t* h3w = (uint32_t*)calloc(R_SIZE_PRIME * 2 + 1, sizeof(uint32_t));
	uint32_t* temp_v = (uint32_t*)calloc(R_SIZE_PRIME * 2 + 1, sizeof(uint32_t));
	uint32_t* temp_w = (uint32_t*)calloc(R_SIZE_PRIME * 2 + 1, sizeof(uint32_t));

	while (Tau >= step_size) {
		uint32_t h[4] = {1,0,0,1};
		uint32_t gm = g[R_SIZE_PRIME - 1];
		uint32_t fm = f[R_SIZE_PRIME - 1];
		gen_hmatrix(&delta,fm,gm,step_size,h);
		//do carry-less multiplication
		clmul_array(f, R_SIZE_PRIME, h[0], h0f);
		clmul_array(g, R_SIZE_PRIME, h[1], h1g);
		clmul_array(f, R_SIZE_PRIME, h[2], h2f);
		clmul_array(g, R_SIZE_PRIME, h[3], h3g);

		xor_array(h0f, h1g, temp_f, R_SIZE_PRIME + 1);
		xor_array(h2f, h3g, temp_g, R_SIZE_PRIME + 1);

		for (int i = 0; i < R_SIZE_PRIME; i++) { f[i] = temp_f[i]; g[i] = temp_g[i]; }

		clmul_array(v, R_SIZE_PRIME * 2, h[0], h0v);
		clmul_array(w, R_SIZE_PRIME * 2, h[1], h1w);
		clmul_array(v, R_SIZE_PRIME * 2, h[2], h2v);
		clmul_array(w, R_SIZE_PRIME * 2, h[3], h3w);

		xor_array(h0v, h1w, temp_v, R_SIZE_PRIME * 2 + 1);
		xor_array(h2v, h3w, temp_w, R_SIZE_PRIME * 2 + 1);

		for (int i = 0; i < (2* R_SIZE_PRIME); i++) { v[i] = temp_v[i]; w[i] = temp_w[i]; }

		Tau -= step_size;
	}
	if (Tau > 0) {
		uint32_t h[4] = { 1,0,0,1 };
		uint32_t gm = g[R_SIZE_PRIME - 1];
		uint32_t fm = f[R_SIZE_PRIME - 1];
		gen_hmatrix(&delta, fm, gm, Tau, h);
		//do carry-less multiplication
		clmul_array(v, R_SIZE_PRIME * 2, h[0], h0v);
		clmul_array(w, R_SIZE_PRIME * 2, h[1], h1w);
		clmul_array(v, R_SIZE_PRIME * 2, h[2], h2v);
		clmul_array(w, R_SIZE_PRIME * 2, h[3], h3w);

		xor_array(h0v, h1w, temp_v, R_SIZE_PRIME * 2 + 1);
		xor_array(h2v, h3w, temp_w, R_SIZE_PRIME * 2 + 1);

		for (int i = 0; i < (2 * R_SIZE_PRIME); i++) { v[i] = temp_v[i]; w[i] = temp_w[i]; }
	}

	//check point
	uint64_t temp_array[R_SIZE_PRIME];
	trans(v,temp_array,2*R_SIZE_PRIME);
	//shift v R_BITS right
	for (int i = 0; i < R_SIZE_PRIME - 1; i++) {
		gout[i] = (((uint64_t)v[i + R_SIZE_PRIME] << SEW) | ((uint64_t)v[i + R_SIZE_PRIME - 1])) >> (R_BITS & (SEW - 1));
	}
	gout[R_SIZE_PRIME - 1] = v[2 * R_SIZE_PRIME - 2] >> (R_BITS & (SEW - 1));
	//free the stack spaces
	free(h0f), free(h1g), free(h2f), free(h3g), free(h0v), free(h1w), free(h2v), free(h3w);
	free(temp_f), free(temp_g), free(temp_v), free(temp_w);
}

