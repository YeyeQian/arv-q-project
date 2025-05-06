#include "poly_op_helper.h"
#include "defs.h"
#include <stdio.h>

template <typename T, typename U>
void General_sparse_poly_mul(IN const T a_in[R_SIZE_PRIME + 1], IN const uint32_t b_in[DV], OUT T res_bin[R_SIZE_PRIME]);

template <typename T, typename U>
void General_poly_inverse(IN const uint32_t step_size, IN const T g_in[R_SIZE_PRIME], OUT T g_out[R_SIZE_PRIME]);

/****************************************************************
* Generalization of the polynomial multiplication and inversion:
* 1. SEW support 8,16,32(defined below)
* 2. VLEN support 128,256,512,1024(defined below)
* 3. r support 12323,24659,40973(defined in defs.h)
****************************************************************/

//U is twice the width as T
template <typename T, typename U>
void General_cyclic_shift_left(IN const T a_in[R_SIZE_PRIME + 1], IN const uint32_t current_shift, OUT T res_bin[R_SIZE_PRIME]) {
	T chunk_mask = (T)CHUNK_MASK;
	uint32_t chunk_solved = 0; //track the number of chunk already dealt with
	uint32_t chunk_shift = current_shift >> SEW_LOG;
	uint32_t bit_shift1 = current_shift & (SEW - 1);
	uint32_t bit_shift2 = bit_shift1 + SEW - (R_BITS & (SEW - 1));
	bool flag0 = bit_shift2 > SEW;
	//loop1: do stripmining on a_in, each round of loop deal with 16(VLMAX) chunk and there are 386(R_SIZE_) chunks
	while (chunk_solved < R_SIZE_PRIME) {
		uint32_t len = (chunk_solved + VLMAX_PRIME) <= R_SIZE_PRIME ? VLMAX_PRIME : (R_SIZE_PRIME - chunk_solved);
		uint32_t temp = chunk_solved + chunk_shift;
		bool flag1 = temp >= R_SIZE_PRIME;
		uint32_t bunch_offset = flag1 ? (temp - R_SIZE_PRIME) : temp;
		//loop2: simulate vector processing on a bunch of chunk
		for (int j = 0; j < len; j++) {
			uint32_t temp = bunch_offset + j;
			bool flag2 = temp >= R_SIZE_PRIME;
			uint32_t chunk_offset = flag2 ? temp - R_SIZE_PRIME : temp;
			uint32_t bit_shift = flag1 || flag2 ? (flag0 ? bit_shift2 - SEW : bit_shift2) : bit_shift1;
			int32_t sp_offset = (flag1 || flag2) && flag0 ? -1 : 0;
			T chunk_shifted = ((((U)a_in[chunk_solved + j + 1 + sp_offset] << SEW | (U)a_in[chunk_solved + j + sp_offset]) << bit_shift) >> SEW) & chunk_mask;
			res_bin[chunk_offset] ^= chunk_shifted;
		}
		chunk_solved += len;
	}
	T final_chunk_mask = (T)FINAL_CHUNK_MASK;
	res_bin[R_SIZE_PRIME - 1] &= final_chunk_mask;
}

template <typename T, typename U>
void General_sparse_poly_mul(IN const T a_in[R_SIZE_PRIME + 1], IN const uint32_t b_in[DV], OUT T res_bin[R_SIZE_PRIME]) {
	//outer loop: go through all the coefficient of sparse polynomial
	for (int i = 0; i < DV; i++) {
		General_cyclic_shift_left<T, U>(a_in, b_in[i], res_bin);
	}
}

template <typename T, typename U>
void General_bitreverse_coeff(IN const T a_in[R_SIZE_PRIME], OUT T res[R_SIZE_PRIME]) {
	//right shift a_in to the rightmost
	uint32_t shift_mount = SEW - (R_BITS & (SEW-1));
	for (int i = R_SIZE_PRIME - 1; i >= 0; i--) {
		T f1 = a_in[i];
		T f0 = i == 0 ? 0 : a_in[i - 1];
		res[i] = (((((U)f1) << SEW | (U)f0) << shift_mount) >> SEW) & CHUNK_MASK;
	}
	//do bitreverse on every element
	for (int i = 0; i < R_SIZE_PRIME; i++)
		reverse(&res[i]);
	//reverse element order
	int left = 0;
	int right = R_SIZE_PRIME - 1;
	T temp;
	while (left < right) {
		temp = res[right];
		res[right] = res[left];
		res[left] = temp;
		left++;
		right--;
	}
}

template <typename T>
void General_compute_control_bits(int32_t* delta, IN const uint32_t step_size, IN const T f[R_SIZE_PRIME], IN const T g[R_SIZE_PRIME], OUT bool c[2 * S]) {
	T f_ = f[0];
	T g_ = g[0];
	for (int i = 0; i < step_size; i++) {
		bool alpha = g_ & 1;
		bool swap = (*delta > 0 ? true : false) && alpha;
		c[2 * i] = swap;
		c[2 * i + 1] = alpha;
		*delta = swap ? -(*delta) + 1 : (*delta) + 1;
		T f_temp = f_;
		f_ = swap ? g_ : f_;
		g_ = alpha ? (g_ ^ f_temp) >> 1 : g_ >> 1;
	}
}

template <typename T,typename U>
void General_update_fg_or_vw(IN const bool c[2 * S], IN const uint32_t step_size, IN const T f1, IN const T f0, IN const T g1, IN const T g0, IN const bool is_updating_fg, OUT T* r0, OUT T* r1) {
	U f = ((U)f1) << SEW | (U)f0;
	U g = ((U)g1) << SEW | (U)g0;
	for (int i = 0; i < step_size; i++) {
		U f_temp = f;
		f = c[2 * i] ? g : f;
		g = c[2 * i + 1] ? g ^ f_temp : g;
		if (is_updating_fg)
			g >>= 1;
		else
			f <<= 1;
	}
	if (is_updating_fg) {
		*r0 = f & CHUNK_MASK;
		*r1 = g & CHUNK_MASK;
	}
	else {
		*r0 = (f >> SEW) & CHUNK_MASK;
		*r1 = (g >> SEW) & CHUNK_MASK;
	}
}

template <typename T, typename U>
void General_poly_inverse(IN const uint32_t step_size, IN const T g_in[R_SIZE_PRIME], OUT T g_out[R_SIZE_PRIME]) {
	T f[R_SIZE_PRIME] = { 0 };
	T g[R_SIZE_PRIME] = { 0 };
	T* v = g_out;//g_out should be initialized to 0 before passed to this function
	T w[R_SIZE_PRIME] = { 0 };
	w[0] = 1, f[0] = 1, f[R_SIZE_PRIME - 1] = 1<<(R_BITS&(SEW-1));
	//reverse g_in to get g
	General_bitreverse_coeff<T,U>(g_in, g);
	int32_t delta = 1;
	uint32_t Tau = 2 * R_BITS - 1;
	bool c[2 * S] = { 0 };

	//define a temp array for test
	uint32_t temp[R_SIZE_PRIME] = { 0 };

	while (Tau >= step_size) {
		General_compute_control_bits<T>(&delta, step_size, f, g, c);
		//do stripmining on the chunks
		//update f g
		uint32_t chunk_solved = 0;
		while (chunk_solved < R_SIZE_PRIME) {
			uint32_t len = (chunk_solved + VLMAX_PRIME) <= R_SIZE_PRIME ? VLMAX_PRIME : (R_SIZE_PRIME - chunk_solved);
			for (int j = 0; j < len; j++) {
				T f0 = f[chunk_solved + j];
				T f1 = chunk_solved + j + 1 >= R_SIZE_PRIME ? 0 : f[chunk_solved + j + 1];
				T g0 = g[chunk_solved + j];
				T g1 = chunk_solved + j + 1 >= R_SIZE_PRIME ? 0 : g[chunk_solved + j + 1];
				T r0, r1;
				General_update_fg_or_vw<T,U>(c, step_size, f1, f0, g1, g0, true, &r0, &r1);
				f[chunk_solved + j] = r0;
				g[chunk_solved + j] = r1;
			}
			chunk_solved += len;
		}

		//temp code for inspecting f
		General_bitreverse_coeff<uint32_t,uint64_t>(g,temp);

		//update v w
		chunk_solved = 0;
		while (chunk_solved < R_SIZE_PRIME) {
			uint32_t len = (chunk_solved + VLMAX_PRIME) <= R_SIZE_PRIME ? VLMAX_PRIME : (R_SIZE_PRIME - chunk_solved);
			for (int j = 0; j < len; j++) {
				T v0 = R_SIZE_PRIME - 1 - (chunk_solved + j) == 0 ? 0 : v[R_SIZE_PRIME - (chunk_solved + j) - 2];
				T v1 = v[R_SIZE_PRIME - 1 - (chunk_solved + j)];
				T w0 = R_SIZE_PRIME - 1 - (chunk_solved + j) == 0 ? 0 : w[R_SIZE_PRIME - (chunk_solved + j) - 2];
				T w1 = w[R_SIZE_PRIME - 1 - (chunk_solved + j)];
				T r0, r1;
				General_update_fg_or_vw<T,U>(c, step_size, v1, v0, w1, w0, false, &r0, &r1);
				v[R_SIZE_PRIME - 1 - (chunk_solved + j)] = r0;
				w[R_SIZE_PRIME - 1 - (chunk_solved + j)] = r1;
			}
			chunk_solved += len;
		}
		Tau -= step_size;
	}
	if (Tau > 0) {
		General_compute_control_bits<T>(&delta, Tau, f, g, c);
		//update v w
		uint32_t chunk_solved = 0;
		while (chunk_solved < R_SIZE_PRIME) {
			uint32_t len = (chunk_solved + VLMAX_PRIME) <= R_SIZE_PRIME ? VLMAX_PRIME : (R_SIZE_PRIME - chunk_solved);
			for (int j = 0; j < len; j++) {
				T v0 = R_SIZE_PRIME - 1 - (chunk_solved + j) == 0 ? 0 : v[R_SIZE_PRIME - (chunk_solved + j) - 2];
				T v1 = v[R_SIZE_PRIME - 1 - (chunk_solved + j)];
				T w0 = R_SIZE_PRIME - 1 - (chunk_solved + j) == 0 ? 0 : w[R_SIZE_PRIME - (chunk_solved + j) - 2];
				T w1 = w[R_SIZE_PRIME - 1 - (chunk_solved + j)];
				T r0, r1;
				General_update_fg_or_vw<T,U>(c, Tau, v1, v0, w1, w0, false, &r0, &r1);
				v[R_SIZE_PRIME - 1 - (chunk_solved + j)] = r0;
				w[R_SIZE_PRIME - 1 - (chunk_solved + j)] = r1;
			}
			chunk_solved += len;
		}
	}
	//shift v one bit right
	T temp_v[R_SIZE_PRIME] = { 0 };
	for (int i = 0; i < R_SIZE_PRIME; i++) {
		T v0 = v[i];
		T v1 = i + 1 >= R_SIZE_PRIME ? 0 : v[i + 1];
		temp_v[i] = ((((U)v1) << SEW | (U)v0) >> 1) & CHUNK_MASK;
	}
	//reverse v
	General_bitreverse_coeff<T,U>(temp_v, v);
}


/**************************************************************************
* Define some other functions for checking the correctness of replacement
**************************************************************************/

//General_sparsepoly_pack, transform 8-bit chunks into SEW-bit chunks
template <typename T>
void General_sparsepoly_pack(IN const uint8_t g_in[R_SIZE], OUT T g_out[R_SIZE_PRIME]) {
	uint8_t pack_size = sizeof(T) / sizeof(uint8_t);
	uint32_t chunk_solved = 0;
	uint32_t des_idx = 0;
	while (chunk_solved + pack_size <= R_SIZE) {
		T temp = 0;
		for (int i = 0; i < pack_size; i++) {
			temp |= (T)g_in[chunk_solved + i] << (i * 8);
		}
		g_out[des_idx] = temp;
		chunk_solved += pack_size;
		des_idx += 1;
	}
	T temp = 0;
	for (int i =0; i+chunk_solved < R_SIZE; i++) {
		temp |= (T)g_in[i + chunk_solved] << (i * 8);
	}
	g_out[des_idx] = temp;
}

//General_sparsepoly_compensate, add a compensation chunk in the first element
template <typename T, typename U>
void General_sparsepoly_compensate(T g[R_SIZE_PRIME + 1]) {
	T additional_chunk = 0;
	uint8_t last_chunk_bits = R_BITS & (SEW - 1);
	additional_chunk = (((U)g[R_SIZE_PRIME] << SEW | (U)g[R_SIZE_PRIME - 1]) >> last_chunk_bits) & CHUNK_MASK;
	g[0] = additional_chunk;
}

//General_sparsepoly_unpack, transform SEW-bit chunks into 8-bit chunks
template <typename T>
void General_sparsepoly_unpack(IN const T g_in[R_SIZE_PRIME], OUT uint8_t g_out[R_SIZE]) {
	uint8_t pack_size = sizeof(T) / sizeof(uint8_t);
	uint32_t chunk_solved = 0;
	uint32_t src_idx = 0;
	uint8_t len = 0;
	while (chunk_solved < R_SIZE) {
		len = (chunk_solved + pack_size <= R_SIZE) ? pack_size : R_SIZE - chunk_solved;
		for (int i = 0; i < len; i++) {
			g_out[chunk_solved + i] = (g_in[src_idx] >> (i * 8)) & 0xff;
		}
		chunk_solved += len;
		src_idx += 1;
	}
}

void poly_mod_mul(OUT uint8_t res_bin[R_SIZE],
	IN const uint8_t a_bin[R_SIZE],
	IN const uint8_t b_bin[R_SIZE]);

void poly_mod_mul_in8(OUT uint8_t res_bin[R_SIZE],
    IN const uint8_t a_bin[R_SIZE],
    IN const uint8_t b_bin[R_SIZE]);

void poly_mod_inv(OUT uint8_t res_bin[R_SIZE],
	IN const uint8_t a_bin[R_SIZE]);

void poly_split_polynomial(OUT uint8_t e0[R_SIZE],
	OUT uint8_t e1[R_SIZE],
	IN const uint8_t e[2 * R_SIZE]);

void poly_add(OUT uint8_t res_bin[R_SIZE],
	IN const uint8_t a_bin[R_SIZE],
	IN const uint8_t b_bin[R_SIZE]);

//=====================================================
//					Customized Version
//=====================================================
void poly_add_custom(OUT uint8_t res_bin[],
    IN const uint8_t a_bin[],
    IN const uint8_t b_bin[], IN const uint32_t len);

void poly_split_polynomial_custom(OUT uint8_t e0[R_SIZE],
    OUT uint8_t e1[R_SIZE],
    IN const uint8_t e[2 * R_SIZE]);

void poly_mul_binary_custom(uint8_t poly_res[R_SIZE],const uint8_t poly_dense[R_SIZE], const uint8_t poly_sparse[R_SIZE]);

void poly_mul_binary_shuffling_custom(uint8_t poly_res[R_SIZE],const uint8_t poly_dense[R_SIZE], const uint8_t poly_sparse[R_SIZE], const uint8_t shuffle_idx[DV]);

void poly_inv_binary_custom8(uint8_t gout[R_SIZE],uint8_t gin[R_SIZE]);

void poly_inv_binary_custom64_wrapper(uint8_t gout[R_SIZE],uint8_t gin[R_SIZE]);

void poly_shiftleft_custom(const uint8_t poly_in[], const uint16_t offset, const uint16_t poly_len, const uint8_t addi, uint8_t poly_out[]);