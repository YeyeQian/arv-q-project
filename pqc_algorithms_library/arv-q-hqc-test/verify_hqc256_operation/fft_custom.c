/**
 * @file fft.cpp
 * @brief Implementation of the additive FFT and its transpose.
 * This implementation is based on the paper from Gao and Mateer: <br>
 * Shuhong Gao and Todd Mateer, Additive Fast Fourier Transforms over Finite Fields,
 * IEEE Transactions on Information Theory 56 (2010), 6265--6272.
 * http://www.math.cpplemson.edu/~sgao/papers/GM10.pdf <br>
 * and includes improvements proposed by Bernstein, Chou and Schwabe here:
 * https://binary.cppr.yp.to/mcbits-20130616.pdf
 */

#include "fft.h"
#include "gf.h"
#include "parameters.h"
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

static void compute_fft_betas(uint16_t *betas);
static void compute_subset_sums(uint16_t *subset_sums, const uint16_t *set, uint16_t set_size);
static void radix(uint16_t *f0, uint16_t *f1, const uint16_t *f, uint32_t m_f);
static void radix_big(uint16_t *f0, uint16_t *f1, const uint16_t *f, uint32_t m_f);
static void fft_rec(uint16_t *w, uint16_t *f, size_t f_coeffs, uint8_t m, uint32_t m_f, const uint16_t *betas);


/**
 * @brief Computes the basis of betas (omitting 1) used in the additive FFT and its transpose
 *
 * @param[out] betas Array of size PARAM_M-1
 */
static void compute_fft_betas(uint16_t *betas) {
    size_t i;
    for (i = 0; i < PARAM_M - 1; ++i) {
        betas[i] = 1 << (PARAM_M - 1 - i);
    }
}



/**
 * @brief Computes the subset sums of the given set
 *
 * The array subset_sums is such that its ith element is
 * the subset sum of the set elements given by the binary form of i.
 *
 * @param[out] subset_sums Array of size 2^set_size receiving the subset sums
 * @param[in] set Array of set_size elements
 * @param[in] set_size Size of the array set
 */
static void compute_subset_sums(uint16_t *subset_sums, const uint16_t *set, uint16_t set_size) {
    uint16_t i, j;
    subset_sums[0] = 0;

    for (i = 0; i < set_size; ++i) {
        for (j = 0; j < (1 << i); ++j) {
            subset_sums[(1 << i) + j] = set[i] ^ subset_sums[j];
        }
    }
}



/**
 * @brief Computes the radix conversion of a polynomial f in GF(2^m)[x]
 *
 * Computes f0 and f1 such that f(x) = f0(x^2-x) + x.f1(x^2-x)
 * as proposed by Bernstein, Chou and Schwabe:
 * https://binary.cppr.yp.to/mcbits-20130616.pdf
 *
 * @param[out] f0 Array half the size of f
 * @param[out] f1 Array half the size of f
 * @param[in] f Array of size a power of 2
 * @param[in] m_f 2^{m_f} is the smallest power of 2 greater or equal to the number of coefficients of f
 */
static void radix(uint16_t *f0, uint16_t *f1, const uint16_t *f, uint32_t m_f) {
    switch (m_f) {
        case 4:
            f0[4] = f[8] ^ f[12];
            f0[6] = f[12] ^ f[14];
            f0[7] = f[14] ^ f[15];
            f1[5] = f[11] ^ f[13];
            f1[6] = f[13] ^ f[14];
            f1[7] = f[15];
            f0[5] = f[10] ^ f[12] ^ f1[5];
            f1[4] = f[9] ^ f[13] ^ f0[5];

            f0[0] = f[0];
            f1[3] = f[7] ^ f[11] ^ f[15];
            f0[3] = f[6] ^ f[10] ^ f[14] ^ f1[3];
            f0[2] = f[4] ^ f0[4] ^ f0[3] ^ f1[3];
            f1[1] = f[3] ^ f[5] ^ f[9] ^ f[13] ^ f1[3];
            f1[2] = f[3] ^ f1[1] ^ f0[3];
            f0[1] = f[2] ^ f0[2] ^ f1[1];
            f1[0] = f[1] ^ f0[1];
            break;

        case 3:
            f0[0] = f[0];
            f0[2] = f[4] ^ f[6];
            f0[3] = f[6] ^ f[7];
            f1[1] = f[3] ^ f[5] ^ f[7];
            f1[2] = f[5] ^ f[6];
            f1[3] = f[7];
            f0[1] = f[2] ^ f0[2] ^ f1[1];
            f1[0] = f[1] ^ f0[1];
            break;

        case 2:
            f0[0] = f[0];
            f0[1] = f[2] ^ f[3];
            f1[0] = f[1] ^ f0[1];
            f1[1] = f[3];
            break;

        case 1:
            f0[0] = f[0];
            f1[0] = f[1];
            break;

        default:
            radix_big(f0, f1, f, m_f);
            break;
    }
}

static void radix_big(uint16_t *f0, uint16_t *f1, const uint16_t *f, uint32_t m_f) {
    uint16_t Q[2 * (1 << (PARAM_FFT - 2))] = {0};
    uint16_t R[2 * (1 << (PARAM_FFT - 2))] = {0};

    uint16_t Q0[1 << (PARAM_FFT - 2)] = {0};
    uint16_t Q1[1 << (PARAM_FFT - 2)] = {0};
    uint16_t R0[1 << (PARAM_FFT - 2)] = {0};
    uint16_t R1[1 << (PARAM_FFT - 2)] = {0};

    size_t i, n;

    n = 1;
    n <<= (m_f - 2);
    memcpy(Q, f + 3 * n, 2 * n);
    memcpy(Q + n, f + 3 * n, 2 * n);
    memcpy(R, f, 4 * n);

    for (i = 0; i < n; ++i) {
        Q[i] ^= f[2 * n + i];
        R[n + i] ^= Q[i];
    }

    radix(Q0, Q1, Q, m_f - 1);
    radix(R0, R1, R, m_f - 1);

    memcpy(f0, R0, 2 * n);
    memcpy(f0 + n, Q0, 2 * n);
    memcpy(f1, R1, 2 * n);
    memcpy(f1 + n, Q1, 2 * n);
}



/**
 * @brief Evaluates f at all subset sums of a given set
 *
 * This function is a subroutine of the function fft.
 *
 * @param[out] w Array
 * @param[in] f Array
 * @param[in] f_coeffs Number of coefficients of f
 * @param[in] m Number of betas
 * @param[in] m_f Number of coefficients of f (one more than its degree)
 * @param[in] betas FFT constants
 */
static void fft_rec(uint16_t *w, uint16_t *f, size_t f_coeffs, uint8_t m, uint32_t m_f, const uint16_t *betas) {
    uint16_t f0[1 << (PARAM_FFT - 2)] = {0};
    uint16_t f1[1 << (PARAM_FFT - 2)] = {0};
    uint16_t gammas[PARAM_M - 2] = {0};
    uint16_t deltas[PARAM_M - 2] = {0};
    uint16_t gammas_sums[1 << (PARAM_M - 2)] = {0};
    uint16_t u[1 << (PARAM_M - 2)] = {0};
    uint16_t v[1 << (PARAM_M - 2)] = {0};
    uint16_t tmp[PARAM_M - (PARAM_FFT - 1)] = {0};

    uint16_t beta_m_pow;
    size_t i, j, k;
    size_t x;

    // Step 1
    if (m_f == 1) {
        for (i = 0; i < m; ++i) {
            tmp[i] = gf_mul(betas[i], f[1]);
        }

        w[0] = f[0];
        x = 1;
        for (j = 0; j < m; ++j) {
            for (k = 0; k < x; ++k) {
                w[x + k] = w[k] ^ tmp[j];
            }
            x <<= 1;
        }

        return;
    }

    // Step 2: compute g
    if (betas[m - 1] != 1) {
        beta_m_pow = 1;
        x = 1;
        x <<= m_f;
        for (i = 1; i < x; ++i) {
            beta_m_pow = gf_mul(beta_m_pow, betas[m - 1]);
            f[i] = gf_mul(beta_m_pow, f[i]);
        }
    }

    // Step 3
    radix(f0, f1, f, m_f);

    // Step 4: compute gammas and deltas
    for (i = 0; i + 1 < m; ++i) {
        gammas[i] = gf_mul(betas[i], gf_inverse(betas[m - 1]));
        deltas[i] = gf_square(gammas[i]) ^ gammas[i];
    }

    // Compute gammas sums
    compute_subset_sums(gammas_sums, gammas, m - 1);

    // Step 5
    fft_rec(u, f0, (f_coeffs + 1) / 2, m - 1, m_f - 1, deltas);

    k = 1;
    k <<= ((m - 1) & 0xf); // &0xf is to let the compiler know that m-1 is small.
    if (f_coeffs <= 3) { // 3-coefficient polynomial f case: f1 is constant
        w[0] = u[0];
        w[k] = u[0] ^ f1[0];
        for (i = 1; i < k; ++i) {
            w[i] = u[i] ^ gf_mul(gammas_sums[i], f1[0]);
            w[k + i] = w[i] ^ f1[0];
        }
    } else {
        fft_rec(v, f1, f_coeffs / 2, m - 1, m_f - 1, deltas);

        // Step 6
        memcpy(w + k, v, 2 * k);
        w[0] = u[0];
        w[k] ^= u[0];
        for (i = 1; i < k; ++i) {
            w[i] = u[i] ^ gf_mul(gammas_sums[i], v[i]);
            w[k + i] ^= w[i];
        }
    }
}


static void fft_rec_custom(uint16_t *w, uint16_t *f, size_t f_coeffs, uint8_t m, uint32_t m_f, const uint16_t *betas) {
    uint16_t f0[1 << (PARAM_FFT - 2)] = {0};
    uint16_t f1[1 << (PARAM_FFT - 2)] = {0};
    uint16_t gammas[PARAM_M - 2] = {0};
    uint16_t deltas[PARAM_M - 2] = {0};
    uint16_t gammas_sums[1 << (PARAM_M - 2)] = {0};
    uint16_t u[1 << (PARAM_M - 2)] = {0};
    uint16_t v[1 << (PARAM_M - 2)] = {0};
    uint16_t tmp[PARAM_M - (PARAM_FFT - 1)] = {0};

    uint16_t beta_m_pow;
    size_t i, j, k;
    size_t x,vl,avl;
    vuint16m1_t va,vb,vc,vd;
    csr_primpoly_rw(PARAM_GF_POLY);
    // Step 1
    if (m_f == 1) {
        gfmul_vx_custom_u16(tmp,betas,f[1],m);

        w[0] = f[0];
        x = 1;
        for (j = 0; j < m; ++j) {
            for (k = 0; k < x; ++k) {
                w[x + k] = w[k] ^ tmp[j];
            }
            x <<= 1;
        }

        return;
    }

    // Step 2: compute g
    if (betas[m - 1] != 1) {
        beta_m_pow = 1;
        x = 1;
        x <<= m_f;
        vl=vsetvl_e16m1(1);
        va=vmv_s_x_u16m1(va,beta_m_pow,vl);
        for (i = 1; i < x; ++i) {
            vsetvl_e16m1_wrapper(1);
            va=vgfmul_vx_u16m1(va,betas[m-1]);
            vb=vgfmul_vx_u16m1(va,f[i]);
            f[i]=vmv_x_s_u16m1_u16(vb);
        }
    }

    // Step 3
    radix(f0, f1, f, m_f);

    // Step 4: compute gammas and deltas
    uint16_t temp_inv=gfinv_u16(betas[m-1]);
    avl=m-1;
    uint16_t* gammas_addr=gammas;
    uint16_t* deltas_addr=deltas;
    const uint16_t* betas_addr=betas;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        va=vle16_v_u16m1(betas_addr,vl);
        va=vgfmul_vx_u16m1(va,temp_inv);
        vse16_v_u16m1(gammas_addr,va,vl);
        vb=vgfmul_vv_u16m1(va,va);
        vb=vxor_vv_u16m1(vb,va,vl);
        vse16_v_u16m1(deltas_addr,vb,vl);
        gammas_addr+=vl,deltas_addr+=vl,betas_addr+=vl,avl-=vl;
    }

    // Compute gammas sums
    compute_subset_sums(gammas_sums, gammas, m - 1);

    // Step 5
    fft_rec_custom(u, f0, (f_coeffs + 1) / 2, m - 1, m_f - 1, deltas);

    k = 1;
    k <<= ((m - 1) & 0xf); // &0xf is to let the compiler know that m-1 is small.
    if (f_coeffs <= 3) { // 3-coefficient polynomial f case: f1 is constant
        w[0] = u[0];
        w[k] = u[0] ^ f1[0];

        avl=k-1;
        uint16_t* gammas_sums_addr=gammas_sums+1;
        uint16_t* w_addr=w+1;
        uint16_t* u_addr=u+1;
        while(avl>0){
            vl=vsetvl_e16m1(avl);
            va=vle16_v_u16m1(gammas_sums_addr,vl);
            vb=vle16_v_u16m1(u_addr,vl);
            va=vgfmul_vx_u16m1(va,f1[0]);
            va=vxor_vv_u16m1(va,vb,vl);
            vse16_v_u16m1(w_addr,va,vl);
            va=vxor_vx_u16m1(va,f1[0],vl);
            vse16_v_u16m1(w_addr+k,va,vl);
            gammas_sums_addr+=vl,w_addr+=vl,u_addr+=vl,avl-=vl;
        }
    } else {
        fft_rec_custom(v, f1, f_coeffs / 2, m - 1, m_f - 1, deltas);

        // Step 6
        memcpy(w + k, v, 2 * k);
        w[0] = u[0];
        w[k] ^= u[0];
        
        avl=k-1;
        uint16_t* gammas_sums_addr=gammas_sums+1;
        uint16_t* w_addr=w+1;
        uint16_t* u_addr=u+1;
        uint16_t* v_addr=v+1;
        while(avl>0){
            vl=vsetvl_e16m1(avl);
            va=vle16_v_u16m1(gammas_sums_addr,vl);
            vb=vle16_v_u16m1(u_addr,vl);
            vc=vle16_v_u16m1(v_addr,vl);
            vd=vle16_v_u16m1(w_addr+k,vl);
            va=vgfmul_vv_u16m1(va,vc);
            va=vxor_vv_u16m1(va,vb,vl);
            vse16_v_u16m1(w_addr,va,vl);

            va=vxor_vv_u16m1(va,vd,vl);
            vse16_v_u16m1(w_addr+k,va,vl);
            gammas_sums_addr+=vl,w_addr+=vl,u_addr+=vl,v_addr+=vl,avl-=vl;
        }
    }
}



/**
 * @brief Evaluates f on all fields elements using an additive FFT algorithm
 *
 * f_coeffs is the number of coefficients of f (one less than its degree). <br>
 * The FFT proceeds recursively to evaluate f at all subset sums of a basis B. <br>
 * This implementation is based on the paper from Gao and Mateer: <br>
 * Shuhong Gao and Todd Mateer, Additive Fast Fourier Transforms over Finite Fields,
 * IEEE Transactions on Information Theory 56 (2010), 6265--6272.
 * http://www.math.cpplemson.edu/~sgao/papers/GM10.pdf <br>
 * and includes improvements proposed by Bernstein, Chou and Schwabe here:
 * https://binary.cppr.yp.to/mcbits-20130616.pdf <br>
 * Note that on this first call (as opposed to the recursive calls to fft_rec), gammas are equal to betas,
 * meaning the first gammas subset sums are actually the subset sums of betas (except 1). <br>
 * Also note that f is altered during computation (twisted at each level).
 *
 * @param[out] w Array
 * @param[in] f Array of 2^PARAM_FFT elements
 * @param[in] f_coeffs Number coefficients of f (i.e. deg(f)+1)
 */
void fft(uint16_t *w, const uint16_t *f, size_t f_coeffs) {
    uint16_t betas[PARAM_M - 1] = {0};
    uint16_t betas_sums[1 << (PARAM_M - 1)] = {0};
    uint16_t f0[1 << (PARAM_FFT - 1)] = {0};
    uint16_t f1[1 << (PARAM_FFT - 1)] = {0};
    uint16_t deltas[PARAM_M - 1] = {0};
    uint16_t u[1 << (PARAM_M - 1)] = {0};
    uint16_t v[1 << (PARAM_M - 1)] = {0};

    size_t i, k;

    // Follows Gao and Mateer algorithm
    compute_fft_betas(betas);

    // Step 1: PARAM_FFT > 1, nothing to do

    // Compute gammas sums
    compute_subset_sums(betas_sums, betas, PARAM_M - 1);

    // Step 2: beta_m = 1, nothing to do

    // Step 3
    radix(f0, f1, f, PARAM_FFT);

    // Step 4: Compute deltas
    for (i = 0; i < PARAM_M - 1; ++i) {
        deltas[i] = gf_square(betas[i]) ^ betas[i];
    }

    // Step 5
    fft_rec(u, f0, (f_coeffs + 1) / 2, PARAM_M - 1, PARAM_FFT - 1, deltas);
    fft_rec(v, f1, f_coeffs / 2, PARAM_M - 1, PARAM_FFT - 1, deltas);

    k = 1 << (PARAM_M - 1);
    // Step 6, 7 and error polynomial computation
    memcpy(w + k, v, 2 * k);

    // Check if 0 is root
    w[0] = u[0];

    // Check if 1 is root
    w[k] ^= u[0];

    // Find other roots
    for (i = 1; i < k; ++i) {
        w[i] = u[i] ^ gf_mul(betas_sums[i], v[i]);
        w[k + i] ^= w[i];
    }
}

void fft_custom(uint16_t *w, const uint16_t *f, size_t f_coeffs) {
    uint16_t betas[PARAM_M - 1] = {0};
    uint16_t betas_sums[1 << (PARAM_M - 1)] = {0};
    uint16_t f0[1 << (PARAM_FFT - 1)] = {0};
    uint16_t f1[1 << (PARAM_FFT - 1)] = {0};
    uint16_t deltas[PARAM_M - 1] = {0};
    uint16_t u[1 << (PARAM_M - 1)] = {0};
    uint16_t v[1 << (PARAM_M - 1)] = {0};
    vuint16m1_t va,vb,vc,vd;

    size_t i, k,vl,avl;
    csr_primpoly_rw(PARAM_GF_POLY);
    // Follows Gao and Mateer algorithm
    compute_fft_betas(betas);

    // Step 1: PARAM_FFT > 1, nothing to do

    // Compute gammas sums
    compute_subset_sums(betas_sums, betas, PARAM_M - 1);

    // Step 2: beta_m = 1, nothing to do

    // Step 3
    radix(f0, f1, f, PARAM_FFT);

    // Step 4: Compute deltas
    avl=PARAM_M-1;
    uint16_t* deltas_addr=deltas;
    uint16_t* betas_addr=betas;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        va=vle16_v_u16m1(betas_addr,vl);
        vb=vgfmul_vv_u16m1(va,va);
        va=vxor_vv_u16m1(va,vb,vl);
        vse16_v_u16m1(deltas_addr,va,vl);
        deltas_addr+=vl,betas_addr+=vl,avl-=vl;
    }

    // Step 5
    fft_rec_custom(u, f0, (f_coeffs + 1) / 2, PARAM_M - 1, PARAM_FFT - 1, deltas);
    fft_rec_custom(v, f1, f_coeffs / 2, PARAM_M - 1, PARAM_FFT - 1, deltas);

    k = 1 << (PARAM_M - 1);
    // Step 6, 7 and error polynomial computation
    memcpy(w + k, v, 2 * k);

    // Check if 0 is root
    w[0] = u[0];

    // Check if 1 is root
    w[k] ^= u[0];

    // Find other roots
    avl=k-1;
    uint16_t* betas_sums_addr=betas_sums+1;
    uint16_t* w_addr=w+1;
    uint16_t* u_addr=u+1;
    uint16_t* v_addr=v+1;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        va=vle16_v_u16m1(betas_sums_addr,vl);
        vb=vle16_v_u16m1(u_addr,vl);
        vc=vle16_v_u16m1(v_addr,vl);
        vd=vle16_v_u16m1(w_addr+k,vl);
        va=vgfmul_vv_u16m1(va,vc);
        va=vxor_vv_u16m1(va,vb,vl);
        vse16_v_u16m1(w_addr,va,vl);

        va=vxor_vv_u16m1(va,vd,vl);
        vse16_v_u16m1(w_addr+k,va,vl);
        betas_sums_addr+=vl,w_addr+=vl,u_addr+=vl,v_addr+=vl,avl-=vl;
    }
}

/**
 * @brief Retrieves the error polynomial error from the evaluations w of the ELP (Error Locator Polynomial) on all field elements.
 *
 * @param[out] error Array with the error
 * @param[out] error_compact Array with the error in a compact form
 * @param[in] w Array of size 2^PARAM_M
 */
void fft_retrieve_error_poly(uint8_t *error, const uint16_t *w) {
    uint16_t gammas[PARAM_M - 1] = {0};
    uint16_t gammas_sums[1 << (PARAM_M - 1)] = {0};
    uint16_t k;
    size_t i, index;

    compute_fft_betas(gammas);
    compute_subset_sums(gammas_sums, gammas, PARAM_M - 1);

    k = 1 << (PARAM_M - 1);
    error[0] ^= 1 ^ ((uint16_t) - w[0] >> 15);
    error[0] ^= 1 ^ ((uint16_t) - w[k] >> 15);

    for (i = 1; i < k; ++i) {
        index = PARAM_GF_MUL_ORDER - gf_log[gammas_sums[i]];
        error[index] ^= 1 ^ ((uint16_t) - w[i] >> 15);

        index = PARAM_GF_MUL_ORDER - gf_log[gammas_sums[i] ^ 1];
        error[index] ^= 1 ^ ((uint16_t) - w[k + i] >> 15);
    }
}

void fft_retrieve_error_poly_custom(uint8_t *error, const uint16_t *w) {
    uint16_t gammas[PARAM_M - 1] = {0};
    uint16_t gammas_sums[1 << (PARAM_M - 1)] = {0};
    uint16_t k;
    size_t i, index,vl,avl;
    vuint16m4_t va,vb,vd,vg;
    vuint8m2_t vc,ve,vf;

    compute_fft_betas(gammas);
    compute_subset_sums(gammas_sums, gammas, PARAM_M - 1);

    k = 1 << (PARAM_M - 1);
    error[0] ^= 1 ^ ((uint16_t) - w[0] >> 15);
    error[0] ^= 1 ^ ((uint16_t) - w[k] >> 15);

    avl=k-1;
    uint16_t* gammas_sums_addr=gammas_sums+1;
    const uint16_t* w_addr=w+1;

    while(avl>0){
        vl=vsetvl_e16m4(avl);
        va=vle16_v_u16m4(gammas_sums_addr,vl);
        vg=vxor_vx_u16m4(va,1,vl);
        va=vsll_vx_u16m4(va,1,vl);
        vb=vluxei16_v_u16m4(gf_log,va,vl);
        vb=vrsub_vx_u16m4(vb,PARAM_GF_MUL_ORDER,vl);
        vd=vle16_v_u16m4(w_addr,vl);
        vd=vrsub_vx_u16m4(vd,0,vl);
        vd=vsrl_vx_u16m4(vd,15,vl);
        vd=vxor_vx_u16m4(vd,1,vl);
        vf=vnsrl_wx_u8m2(vd,0,vl);
        ve=vnsrl_wx_u8m2(vb,0,vl);//index
        vl=vsetvl_e8m2(avl);
        vc=vluxei8_v_u8m2(error,ve,vl);//error[index]
        vc=vxor_vv_u8m2(vc,vf,vl);
        vsuxei8_v_u8m2(error,ve,vc,vl);


        vg=vsll_vx_u16m4(vg,1,vl);
        vb=vluxei16_v_u16m4(gf_log,vg,vl);
        vb=vrsub_vx_u16m4(vb,PARAM_GF_MUL_ORDER,vl);
        vd=vle16_v_u16m4(w_addr+k,vl);
        vd=vrsub_vx_u16m4(vd,0,vl);
        vd=vsrl_vx_u16m4(vd,15,vl);
        vd=vxor_vx_u16m4(vd,1,vl);
        vf=vnsrl_wx_u8m2(vd,0,vl);
        ve=vnsrl_wx_u8m2(vb,0,vl);//index
        vl=vsetvl_e8m2(avl);
        vc=vluxei8_v_u8m2(error,ve,vl);//error[index]
        vc=vxor_vv_u8m2(vc,vf,vl);
        vsuxei8_v_u8m2(error,ve,vc,vl);

        gammas_sums_addr+=vl,w_addr+=vl,avl-=vl;
    }
}

///////////////////////////////////////////////////////////////////////////
//                               Tests
///////////////////////////////////////////////////////////////////////////
bool test_fft_rec_custom(){
    uint16_t u_ref[1 << (PARAM_M - 1)] = {
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
    };
    uint16_t u_res[1 << (PARAM_M - 1)] = {0};
    uint16_t f0[1 << (PARAM_FFT - 1)] = {
        1,0,0,0,0,0,0,0
    };
    uint16_t deltas[PARAM_M - 1] = {
        147,141,84,13,72,20,6
    };

    fft_rec_custom(u_res, f0, 8, PARAM_M - 1, PARAM_FFT - 1, deltas);

    //CHECK
    bool flag=true;
    for(int i=0;i<(1 << (PARAM_M - 1));i++){
        if(u_ref[i]!=u_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_fft_rec_custom return with flag %d\n",flag);

    return flag;
}

bool test_fft_custom(){
    uint16_t sigma[1 << PARAM_FFT] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<(1<<PARAM_FFT);i++){
        sigma[i]=rand()&255;
    }

    uint16_t w_ref[1 << PARAM_M] = {0};
    uint16_t w_res[1 << PARAM_M] = {0};

    uint64_t start, end;

    start=read_cycle();
    fft(w_ref, sigma, PARAM_DELTA + 1);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    fft_custom(w_res, sigma, PARAM_DELTA + 1);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(1 << PARAM_M);i++){
        if(w_ref[i]!=w_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_fft_custom return with flag %d\n",flag);

    return flag;
}

bool test_fft_retrieve_error_poly_custom(){
    uint16_t w[1 << PARAM_M] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<(1<<PARAM_M);i++){
        w[i]=rand()&255;
    }

    uint8_t error_ref[1 << PARAM_M] = {0};
    uint8_t error_res[1 << PARAM_M] = {0};

    uint64_t start, end;

    start=read_cycle();
    fft_retrieve_error_poly(error_ref, w);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    fft_retrieve_error_poly_custom(error_res, w);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(1 << PARAM_M);i++){
        if(error_ref[i]!=error_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_fft_retrieve_error_poly_custom return with flag %d\n",flag);

    return flag;
}

// int main(){
//     // test_fft_rec_custom();
//     // test_fft_custom();
//     test_fft_retrieve_error_poly_custom();

//     return 0;
// }