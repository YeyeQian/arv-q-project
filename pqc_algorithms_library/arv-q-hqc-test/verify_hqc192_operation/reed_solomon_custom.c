/**
 * @file reed_solomon.cpp
 * @brief Constant time implementation of Reed-Solomon codes
 */

#include "fft.h"
#include "gf.h"
#include "reed_solomon.h"
#include "parameters.h"
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#ifdef VERBOSE
#include <stdbool.h>
#include <stdio.h>
#endif

static uint16_t mod(uint16_t i, uint16_t modulus);
static void compute_syndromes(uint16_t *syndromes, uint8_t *cdw);
static uint16_t compute_elp(uint16_t *sigma, const uint16_t *syndromes);
static void compute_roots(uint8_t *error, uint16_t *sigma);
static void compute_z_poly(uint16_t *z, const uint16_t *sigma, const uint16_t degree, const uint16_t *syndromes);
static void compute_error_values(uint16_t *error_values, const uint16_t *z, const uint8_t *error);
static void correct_errors(uint8_t *cdw, const uint16_t *error_values);


/**
 * Returns i modulo the given modulus.
 * i must be less than 2*modulus.
 * Therefore, the return value is either i or i-modulus.
 * @returns i mod (modulus)
 * @param[in] i The integer whose modulo is taken
 * @param[in] modulus The modulus
 */
static uint16_t mod(uint16_t i, uint16_t modulus) {
    uint16_t tmp = i - modulus;

    // mask = 0xffff if(i < PARAM_GF_MUL_ORDER)
    int16_t mask = -(tmp >> 15);

    return tmp + (mask & modulus);
}



/**
 * @brief Computes the generator polynomial of the primitive Reed-Solomon code with given parameters.
 *
 * Code length is 2^m-1. <br>
 * PARAM_DELTA is the targeted correction capacity of the code
 * and receives the real correction capacity (which is at least equal to the target). <br>
 * gf_exp and gf_log are arrays giving antilog and log of GF(2^m) elements.
 *
 * @param[out] poly Array of size (2*PARAM_DELTA + 1) receiving the coefficients of the generator polynomial
 */
void compute_generator_poly(uint16_t* poly) {
    poly[0] = 1;
    int tmp_degree = 0;

    for (uint16_t i = 1; i < (2 * PARAM_DELTA + 1); ++i) {

        for(size_t j = tmp_degree; j; --j) {
            poly[j] = gf_exp[mod(gf_log[poly[j]] + i, PARAM_GF_MUL_ORDER)] ^ poly[j - 1];
        }

        poly[0] = gf_exp[mod(gf_log[poly[0]] + i, PARAM_GF_MUL_ORDER)];
        poly[++tmp_degree] = 1;

    }

    printf("\n");
    for (int i = 0; i < (PARAM_G); ++i) {
        printf("%d, ", poly[i]);
    }
    printf("\n");
}



/**
 * @brief Encodes a message message of PARAM_K bits to a Reed-Solomon codeword codeword of PARAM_N1 bytes
 *
 * Following @cite lin1983error (Chapter 4 - Cyclic Codes),
 * We perform a systematic encoding using a linear (PARAM_N1 - PARAM_K)-stage shift register
 * with feedback connections based on the generator polynomial PARAM_RS_POLY of the Reed-Solomon code.
 *
 * @param[out] cdw Array of size VEC_N1_SIZE_64 receiving the encoded message
 * @param[in] msg Array of size VEC_K_SIZE_64 storing the message
 */
void reed_solomon_encode(uint64_t *cdw, const uint64_t *msg) {
    size_t i, j, k;
    uint8_t gate_value = 0;

    uint16_t tmp[PARAM_G] = {0};
    uint16_t PARAM_RS_POLY [] = {RS_POLY_COEFS};

    uint8_t msg_bytes[PARAM_K] = {0};
    uint8_t cdw_bytes[PARAM_N1] = {0};

    memcpy(msg_bytes, msg, PARAM_K);

    for (i = 0; i < PARAM_K; ++i) {
        gate_value = msg_bytes[PARAM_K - 1 - i] ^ cdw_bytes[PARAM_N1 - PARAM_K - 1];

        for (j = 0; j < PARAM_G; ++j) {
            tmp[j] = gf_mul(gate_value, PARAM_RS_POLY[j]);
        }

        for(k = PARAM_N1 - PARAM_K - 1; k; --k) {
            cdw_bytes[k] = cdw_bytes[k - 1] ^ tmp[k];
        }

        cdw_bytes[0] = tmp[0];
    }

    memcpy(cdw_bytes + PARAM_N1 - PARAM_K, msg_bytes, PARAM_K);
    memcpy(cdw, cdw_bytes, PARAM_N1);
}



/**
 * @brief Computes 2 * PARAM_DELTA syndromes
 *
 * @param[out] syndromes Array of size 2 * PARAM_DELTA receiving the computed syndromes
 * @param[in] cdw Array of size PARAM_N1 storing the received vector
 */
void compute_syndromes(uint16_t *syndromes, uint8_t *cdw) {	
    for (size_t i = 0; i < 2 * PARAM_DELTA; ++i) {
        for (size_t j = 1; j < PARAM_N1; ++j) {
            syndromes[i] ^= gf_mul(cdw[j], alpha_ij_pow[i][j-1]);
        }
        syndromes[i] ^= cdw[0];
    }
}



/**
 * @brief Computes the error locator polynomial (ELP) sigma
 *
 * This is a constant time implementation of Berlekamp's simplified algorithm (see @cite lin1983error (Chapter 6 - BCH Codes). <br>
 * We use the letter p for rho which is initialized at -1. <br>
 * The array X_sigma_p represents the polynomial X^(mu-rho)*sigma_p(X). <br>
 * Instead of maintaining a list of sigmas, we update in place both sigma and X_sigma_p. <br>
 * sigma_copy serves as a temporary save of sigma in case X_sigma_p needs to be updated. <br>
 * We can properly correct only if the degree of sigma does not exceed PARAM_DELTA.
 * This means only the first PARAM_DELTA + 1 coefficients of sigma are of value
 * and we only need to save its first PARAM_DELTA - 1 coefficients.
 *
 * @returns the degree of the ELP sigma
 * @param[out] sigma Array of size (at least) PARAM_DELTA receiving the ELP
 * @param[in] syndromes Array of size (at least) 2*PARAM_DELTA storing the syndromes
 */
static uint16_t compute_elp(uint16_t *sigma, const uint16_t *syndromes) {
    uint16_t deg_sigma = 0;
    uint16_t deg_sigma_p = 0;
    uint16_t deg_sigma_copy = 0;
    uint16_t sigma_copy[PARAM_DELTA + 1] = {0};
    uint16_t X_sigma_p[PARAM_DELTA + 1] = {0,1};
    uint16_t pp = (uint16_t) -1; // 2*rho
    uint16_t d_p = 1;
    uint16_t d = syndromes[0];

    uint16_t mask1, mask2, mask12;
    uint16_t deg_X, deg_X_sigma_p;
    uint16_t dd;
    uint16_t mu;

    uint16_t i;

    sigma[0] = 1;
    for (mu = 0; (mu < (2 * PARAM_DELTA)); ++mu) {
        // Save sigma in case we need it to update X_sigma_p
        memcpy(sigma_copy, sigma, 2 * (PARAM_DELTA));
        deg_sigma_copy = deg_sigma;

        dd = gf_mul(d, gf_inverse(d_p));

        for (i = 1; (i <= mu + 1) && (i <= PARAM_DELTA); ++i) {
            sigma[i] ^= gf_mul(dd, X_sigma_p[i]);
        }

        deg_X = mu - pp;
        deg_X_sigma_p = deg_X + deg_sigma_p;

        // mask1 = 0xffff if(d != 0) and 0 otherwise
        mask1 = -((uint16_t) - d >> 15);

        // mask2 = 0xffff if(deg_X_sigma_p > deg_sigma) and 0 otherwise
        mask2 = -((uint16_t) (deg_sigma - deg_X_sigma_p) >> 15);

        // mask12 = 0xffff if the deg_sigma increased and 0 otherwise
        mask12 = mask1 & mask2;
        deg_sigma ^= mask12 & (deg_X_sigma_p ^ deg_sigma);

        if (mu == (2 * PARAM_DELTA - 1)) {
            break;
        }

        pp ^= mask12 & (mu ^ pp);
        d_p ^= mask12 & (d ^ d_p);

        for (i = PARAM_DELTA; i ;--i) {
            X_sigma_p[i] = (mask12 & sigma_copy[i - 1]) ^ (~mask12 & X_sigma_p[i - 1]);
        }

        deg_sigma_p ^= mask12 & (deg_sigma_copy ^ deg_sigma_p);
        d = syndromes[mu + 1];

        for (i = 1; (i <= mu+1) && (i <= PARAM_DELTA); ++i) {
            d ^= gf_mul(sigma[i], syndromes[mu + 1 - i]);
        }
    }

    return deg_sigma;
}



/**
 * @brief Computes the error polynomial error from the error locator polynomial sigma
 *
 * See function fft for more details.
 *
 * @param[out] error Array of 2^PARAM_M elements receiving the error polynomial
 * @param[out] error_compact Array of PARAM_DELTA + PARAM_N1 elements receiving a compact representation of the vector error
 * @param[in] sigma Array of 2^PARAM_FFT elements storing the error locator polynomial
 */
static void compute_roots(uint8_t *error, uint16_t *sigma) {
    uint16_t w[1 << PARAM_M] = {0};

    fft(w, sigma, PARAM_DELTA + 1);
    fft_retrieve_error_poly(error, w);
}



/**
 * @brief Computes the polynomial z(x)
 *
 * See @cite lin1983error (Chapter 6 - BCH Codes) for more details.
 *
 * @param[out] z Array of PARAM_DELTA + 1 elements receiving the polynomial z(x)
 * @param[in] sigma Array of 2^PARAM_FFT elements storing the error locator polynomial
 * @param[in] degree Integer that is the degree of polynomial sigma
 * @param[in] syndromes Array of 2 * PARAM_DELTA storing the syndromes
 */
static void compute_z_poly(uint16_t *z, const uint16_t *sigma, const uint16_t degree, const uint16_t *syndromes) {
    size_t i, j;
    uint16_t mask;

    z[0] = 1;

    for (i = 1; i < PARAM_DELTA + 1; ++i) {
        mask = -((uint16_t) (i - degree - 1) >> 15);
        z[i] = mask & sigma[i];
    }

    z[1] ^= syndromes[0];

    for (i = 2; i <= PARAM_DELTA; ++i) {
        mask = -((uint16_t) (i - degree - 1) >> 15);
        z[i] ^= mask & syndromes[i - 1];

        for (j = 1; j < i; ++j) {
            z[i] ^= mask & gf_mul(sigma[j], syndromes[i - j - 1]);
        }
    }
}



/**
 * @brief Computes the error values
 *
 * See @cite lin1983error (Chapter 6 - BCH Codes) for more details.
 *
 * @param[out] error_values Array of PARAM_DELTA elements receiving the error values
 * @param[in] z Array of PARAM_DELTA + 1 elements storing the polynomial z(x)
 * @param[in] z_degree Integer that is the degree of polynomial z(x)
 * @param[in] error_compact Array of PARAM_DELTA + PARAM_N1 storing compact representation of the error
 */
static void compute_error_values(uint16_t *error_values, const uint16_t *z, const uint8_t *error) {
    uint16_t beta_j[PARAM_DELTA] = {0};
    uint16_t e_j[PARAM_DELTA] = {0};

    uint16_t delta_counter;
    uint16_t delta_real_value;
    uint16_t found;
    uint16_t mask1;
    uint16_t mask2;
    uint16_t tmp1;
    uint16_t tmp2;
    uint16_t inverse;
    uint16_t inverse_power_j;

    // Compute the beta_{j_i} page 31 of the documentation
    delta_counter = 0;

    for (size_t i = 0; i < PARAM_N1; i++) {
        found = 0;
        mask1 = (uint16_t) (-((int32_t)error[i]) >> 31); // error[i] != 0
        for (size_t j = 0; j < PARAM_DELTA; j++) {
            mask2 = ~((uint16_t) (-((int32_t) j ^ delta_counter) >> 31)); // j == delta_counter
            beta_j[j] += mask1 & mask2 & gf_exp[i];
            found += mask1 & mask2 & 1;
        }

        delta_counter += found;
    }
    delta_real_value = delta_counter;

    // Compute the e_{j_i} page 31 of the documentation
    for (size_t i = 0; i < PARAM_DELTA; ++i) {
        tmp1 = 1;
        tmp2 = 1;
        inverse = gf_inverse(beta_j[i]);
        inverse_power_j = 1;

        for (size_t j = 1; j <= PARAM_DELTA; ++j) {
            inverse_power_j = gf_mul(inverse_power_j, inverse);
            tmp1 ^= gf_mul(inverse_power_j, z[j]);
        }
        for (size_t k = 1; k < PARAM_DELTA; ++k) {
            tmp2 = gf_mul(tmp2, (1 ^ gf_mul(inverse, beta_j[(i + k) % PARAM_DELTA])));
        }
        mask1 = (uint16_t) (((int16_t) i - delta_real_value) >> 15); // i < delta_real_value
        e_j[i] = mask1 & gf_mul(tmp1, gf_inverse(tmp2));
    }

    // Place the delta e_{j_i} values at the right coordinates of the output vector
    delta_counter = 0;
    for (size_t i = 0; i < PARAM_N1; ++i) {
        found = 0;
        mask1 = (uint16_t) (-((int32_t)error[i]) >> 31); // error[i] != 0
        for (size_t j = 0; j < PARAM_DELTA; j++) {
            mask2 = ~((uint16_t) (-((int32_t) j ^ delta_counter) >> 31)); // j == delta_counter
            error_values[i] += mask1 & mask2 & e_j[j];
            found += mask1 & mask2 & 1;
        }
        delta_counter += found;
    }
}



/**
 * @brief Correct the errors
 *
 * @param[out] cdw Array of PARAM_N1 elements receiving the corrected vector
 * @param[in] error Array of the error vector
 * @param[in] error_values Array of PARAM_DELTA elements storing the error values
 */
static void correct_errors(uint8_t *cdw, const uint16_t *error_values) {
    for (size_t i = 0; i < PARAM_N1; ++i) {
        cdw[i] ^= error_values[i];
    }
}



/**
 * @brief Decodes the received word
 *
 * This function relies on six steps:
 *    <ol>
 *    <li> The first step, is the computation of the 2*PARAM_DELTA syndromes.
 *    <li> The second step is the computation of the error-locator polynomial sigma.
 *    <li> The third step, done by additive FFT, is finding the error-locator numbers by calculating the roots of the polynomial sigma and takings their inverses.
 *    <li> The fourth step, is the polynomial z(x).
 *    <li> The fifth step, is the computation of the error values.
 *    <li> The sixth step is the correction of the errors in the received polynomial.
 *    </ol>
 * For a more complete picture on Reed-Solomon decoding, see Shu. Lin and Daniel J. Costello in Error Control Coding: Fundamentals and Applications @cite lin1983error
 *
 * @param[out] msg Array of size VEC_K_SIZE_64 receiving the decoded message
 * @param[in] cdw Array of size VEC_N1_SIZE_64 storing the received word
 */
void reed_solomon_decode(uint64_t *msg, uint64_t *cdw) {
    uint8_t cdw_bytes[PARAM_N1] = {0};
    uint16_t syndromes[2 * PARAM_DELTA] = {0};
    uint16_t sigma[1 << PARAM_FFT] = {0};
    uint8_t error[1 << PARAM_M] = {0};
    uint16_t z[PARAM_N1] = {0};
    uint16_t error_values[PARAM_N1] = {0};
    uint16_t deg;

    // Copy the vector in an array of bytes
    memcpy(cdw_bytes, cdw, PARAM_N1);

    // Calculate the 2*PARAM_DELTA syndromes
    compute_syndromes(syndromes, cdw_bytes);

    // Compute the error locator polynomial sigma
    // Sigma's degree is at most PARAM_DELTA but the FFT requires the extra room
    deg = compute_elp(sigma, syndromes);

    // Compute the error polynomial error
    compute_roots(error, sigma);

    // Compute the polynomial z(x)
    compute_z_poly(z, sigma, deg, syndromes);

    // Compute the error values
    compute_error_values(error_values, z, error);

    // Correct the errors
    correct_errors(cdw_bytes, error_values);

    // Retrieve the message from the decoded codeword
    memcpy(msg, cdw_bytes + (PARAM_G - 1) , PARAM_K);

    #ifdef VERBOSE
        printf("\n\nThe syndromes: ");
        for (size_t i = 0 ; i < 2*PARAM_DELTA ; ++i) {
            printf("%u ", syndromes[i]);
        }
        printf("\n\nThe error locator polynomial: sigma(x) = ");
        bool first_coeff = true;
        if (sigma[0]) {
            printf("%u", sigma[0]);
            first_coeff = false;
        }
        for (size_t i = 1 ; i < (1 << PARAM_FFT) ; ++i) {
            if (sigma[i] == 0)
                continue;
            if (!first_coeff)
                printf(" + ");
            first_coeff = false;
            if(sigma[i] != 1)
                printf("%u ", sigma[i]);
            if (i == 1)
                printf("x");
            else
                printf("x^%zu", i);
        }
        if (first_coeff)
            printf("0");

        printf("\n\nThe polynomial: z(x) = ");
        bool first_coeff_1 = true;
        if (z[0]) {
            printf("%u", z[0]);
            first_coeff_1 = false;
        }
        for (size_t i = 1 ; i < (PARAM_DELTA + 1) ; ++i) {
            if (z[i] == 0)
                continue;
            if (!first_coeff_1)
                printf(" + ");
            first_coeff_1 = false;
            if(z[i] != 1)
                printf("%u ", z[i]);
            if (i == 1)
                printf("x");
            else
                printf("x^%zu", i);
        }
        if (first_coeff_1)
            printf("0");

        printf("\n\nThe pairs of (error locator numbers, error values): ");
        size_t j = 0;
        for (size_t i = 0 ; i < PARAM_N1 ; ++i) {
            if(error[i]){
                printf("(%zu, %d) ", i, error_values[j]);
                j++;
            }
        }
        printf("\n");
    #endif
}


///////////////////////////////////////////////////////////////////////////
//                      Customized Versions
///////////////////////////////////////////////////////////////////////////
void reed_solomon_encode_custom(uint64_t *cdw, const uint64_t *msg) {
    size_t i, j, k;
    uint8_t gate_value = 0;

    uint8_t tmp[PARAM_G] = {0};
    uint8_t PARAM_RS_POLY [] = {RS_POLY_COEFS};

    uint8_t msg_bytes[PARAM_K] = {0};
    uint8_t cdw_bytes[PARAM_N1] = {0};

    memcpy(msg_bytes, msg, PARAM_K);

    vuint8m1_t va,vb;

    for (i = 0; i < PARAM_K; ++i) {
        gate_value = msg_bytes[PARAM_K - 1 - i] ^ cdw_bytes[PARAM_N1 - PARAM_K - 1];

        gfmul_vx_custom_u8(tmp,PARAM_RS_POLY,gate_value,PARAM_G);

        size_t vl,avl;
        avl=PARAM_N1-PARAM_K-1;
        uint8_t* src1_addr=cdw_bytes;
        uint8_t* src2_addr=tmp+1;
        uint8_t* dst_addr=cdw_bytes+1;
        while(avl>0){
            vl=vsetvl_e8m1(avl);
            va=vle8_v_u8m1(src1_addr,vl);
            vb=vle8_v_u8m1(src2_addr,vl);
            va=vxor_vv_u8m1(va,vb,vl);
            vse8_v_u8m1(dst_addr,va,vl);
            src1_addr+=vl,src2_addr+=vl,dst_addr+=vl,avl-=vl;
        }

        cdw_bytes[0] = tmp[0];
    }

    memcpy(cdw_bytes + PARAM_N1 - PARAM_K, msg_bytes, PARAM_K);
    memcpy(cdw, cdw_bytes, PARAM_N1);
}

void compute_syndromes_custom(uint16_t *syndromes, uint8_t *cdw) {
    size_t vl,avl;
    vuint8m2_t va,vb;
    vuint8m1_t vc;
    csr_primpoly_rw(PARAM_GF_POLY);
    for (size_t i = 0; i < 2 * PARAM_DELTA; ++i) {
        avl=PARAM_N1-1;
        uint8_t* src1_addr=cdw+1;
        const uint8_t* src2_addr=&alpha_ij_pow[i][0];
        vl=vsetvl_e8m1(1);
        vc=vmv_v_x_u8m1(cdw[0],vl);
        while(avl>0){
            vl=vsetvl_e8m2(avl);
            va=vle8_v_u8m2(src1_addr,vl);
            vb=vle8_v_u8m2(src2_addr,vl);
            va=vgfmul_vv_u8m2(va,vb);
            vc=vredxor_vs_u8m2_u8m1(vc,va,vc,vl);
            src1_addr+=vl,src2_addr+=vl,avl-=vl;
        }
        syndromes[i]=vmv_x_s_u8m1_u8(vc);
    }
}

// static uint16_t compute_elp_custom(uint16_t *sigma, const uint16_t *syndromes) {
//     uint16_t deg_sigma = 0;
//     uint16_t deg_sigma_p = 0;
//     uint16_t deg_sigma_copy = 0;
//     uint16_t sigma_copy[PARAM_DELTA + 1] = {0};
//     uint16_t X_sigma_p[PARAM_DELTA + 1] = {0,1};
//     uint16_t pp = (uint16_t) -1; // 2*rho
//     uint16_t d_p = 1;
//     uint16_t d = syndromes[0];

//     uint16_t mask1, mask2, mask12;
//     uint16_t deg_X, deg_X_sigma_p;
//     uint16_t dd;
//     uint16_t mu;

//     uint16_t i;

//     sigma[0] = 1;
//     size_t vl,avl;
//     csr_primpoly_rw(PARAM_GF_POLY);

//     vuint16m1_t va;
//     vuint16m2_t vb,vc;
//     vuint16m2_t vidx;

//     for (mu = 0; (mu < (2 * PARAM_DELTA)); ++mu) {
//         // Save sigma in case we need it to update X_sigma_p
//         memcpy(sigma_copy, sigma, 2 * (PARAM_DELTA));
//         deg_sigma_copy = deg_sigma;

//         uint16_t inv_d_p=gfinv_u16(d_p);
//         vl=vsetvl_e16m1(1);
//         va=vmv_s_x_u16m1(va,d,vl);
//         va=vgfmul_vx_u16m1(va,inv_d_p);
//         dd=vmv_x_s_u16m1_u16(va);

//         uint32_t len=(PARAM_DELTA>=(mu+1))?mu+1:PARAM_DELTA;
//         avl=len;
//         uint16_t* X_sigma_p_addr=X_sigma_p+1;
//         uint16_t* sigma_addr=sigma+1;
//         while(avl>0){
//             vl=vsetvl_e16m2(avl);
//             vb=vle16_v_u16m2(X_sigma_p_addr,vl);
//             vc=vle16_v_u16m2(sigma_addr,vl);
//             vb=vgfmul_vx_u16m2(vb,dd);
//             vc=vxor_vv_u16m2(vb,vc,vl);
//             vse16_v_u16m2(sigma_addr,vc,vl);
//             sigma_addr+=vl,X_sigma_p_addr+=vl,avl-=vl;
//         }

//         deg_X = mu - pp;
//         deg_X_sigma_p = deg_X + deg_sigma_p;

//         // mask1 = 0xffff if(d != 0) and 0 otherwise
//         mask1 = -((uint16_t) - d >> 15);

//         // mask2 = 0xffff if(deg_X_sigma_p > deg_sigma) and 0 otherwise
//         mask2 = -((uint16_t) (deg_sigma - deg_X_sigma_p) >> 15);

//         // mask12 = 0xffff if the deg_sigma increased and 0 otherwise
//         mask12 = mask1 & mask2;
//         deg_sigma ^= mask12 & (deg_X_sigma_p ^ deg_sigma);

//         if (mu == (2 * PARAM_DELTA - 1)) {
//             break;
//         }

//         pp ^= mask12 & (mu ^ pp);
//         d_p ^= mask12 & (d ^ d_p);
        
//         avl=PARAM_DELTA;
//         uint16_t* src1_addr=sigma_copy;
//         uint16_t* src2_addr=X_sigma_p;
//         uint16_t* dst_addr=X_sigma_p+1;
//         while(avl>0){
//             vl=vsetvl_e16m2(avl);
//             vb=vle16_v_u16m2(src1_addr,vl);
//             vb=vand_vx_u16m2(vb,mask12,vl);
//             vc=vle16_v_u16m2(src2_addr,vl);
//             vc=vand_vx_u16m2(vc,~mask12,vl);
//             vb=vxor_vv_u16m2(vb,vc,vl);
//             vse16_v_u16m2(dst_addr,vb,vl);
//             src1_addr+=vl,src2_addr+=vl,dst_addr+=vl,avl-=vl;
//         }

//         deg_sigma_p ^= mask12 & (deg_sigma_copy ^ deg_sigma_p);
//         d = syndromes[mu + 1];

//         avl=len;
//         vl=vsetvl_e16m1(1);
//         va=vmv_s_x_u16m1(va,d,vl);
//         src1_addr=sigma+1;
//         while(avl>0){
//             vl=vsetvl_e16m2(avl);
//             vidx=vid_v_u16m2(vl);
//             vidx=vrsub_vx_u16m2(vidx,mu+avl-len,vl);
//             vidx=vsll_vx_u16m2(vidx,1,vl);
//             vb=vle16_v_u16m2(src1_addr,vl);
//             vc=vluxei16_v_u16m2(syndromes,vidx,vl);
//             vb=vgfmul_vv_u16m2(vb,vc);
//             va=vredxor_vs_u16m2_u16m1(va,vb,va,vl);
//             src1_addr+=vl,avl-=vl;
//         }
//         d=vmv_x_s_u16m1_u16(va);
//     }

//     return deg_sigma;
// }

static uint16_t compute_elp_custom(uint16_t *sigma, const uint16_t *syndromes) {
    uint16_t deg_sigma = 0;
    uint16_t deg_sigma_p = 0;
    uint16_t deg_sigma_copy = 0;
    uint16_t sigma_copy[PARAM_DELTA + 1] = {0};
    uint16_t X_sigma_p[PARAM_DELTA + 1] = {0,1};
    uint16_t pp = (uint16_t) -1; // 2*rho
    uint16_t d_p = 1;
    uint16_t d = syndromes[0];

    uint16_t mask1, mask2, mask12;
    uint16_t deg_X, deg_X_sigma_p;
    uint16_t dd;
    uint16_t mu;

    uint16_t i;

    sigma[0] = 1;
    size_t vl,avl;
    csr_primpoly_rw(PARAM_GF_POLY);

    vuint16m1_t va;
    vuint16m2_t vb,vc;
    vuint16m2_t vidx;

    for (mu = 0; (mu < (2 * PARAM_DELTA)); ++mu) {
        // Save sigma in case we need it to update X_sigma_p
        memcpy(sigma_copy, sigma, 2 * (PARAM_DELTA));
        deg_sigma_copy = deg_sigma;

        uint16_t inv_d_p=gfinv_u16(d_p);
        vl=vsetvl_e16m1(1);
        va=vmv_s_x_u16m1(va,d,vl);
        va=vgfmul_vx_u16m1(va,inv_d_p);
        dd=vmv_x_s_u16m1_u16(va);

        uint32_t len=(PARAM_DELTA>=(mu+1))?mu+1:PARAM_DELTA;
        avl=len;
        uint16_t* X_sigma_p_addr=X_sigma_p+1;
        uint16_t* sigma_addr=sigma+1;
        while(avl>0){
            vl=vsetvl_e16m2(avl);
            vb=vle16_v_u16m2(X_sigma_p_addr,vl);
            vc=vle16_v_u16m2(sigma_addr,vl);
            vb=vgfmul_vx_u16m2(vb,dd);
            vc=vxor_vv_u16m2(vb,vc,vl);
            vse16_v_u16m2(sigma_addr,vc,vl);
            sigma_addr+=vl,X_sigma_p_addr+=vl,avl-=vl;
        }

        deg_X = mu - pp;
        deg_X_sigma_p = deg_X + deg_sigma_p;

        // mask1 = 0xffff if(d != 0) and 0 otherwise
        mask1 = -((uint16_t) - d >> 15);

        // mask2 = 0xffff if(deg_X_sigma_p > deg_sigma) and 0 otherwise
        mask2 = -((uint16_t) (deg_sigma - deg_X_sigma_p) >> 15);

        // mask12 = 0xffff if the deg_sigma increased and 0 otherwise
        mask12 = mask1 & mask2;
        deg_sigma ^= mask12 & (deg_X_sigma_p ^ deg_sigma);

        if (mu == (2 * PARAM_DELTA - 1)) {
            break;
        }

        pp ^= mask12 & (mu ^ pp);
        d_p ^= mask12 & (d ^ d_p);
        
        avl=PARAM_DELTA;
        uint16_t* src1_addr=sigma_copy;
        uint16_t* src2_addr=X_sigma_p;
        uint16_t* dst_addr=X_sigma_p+1;
        while(avl>0){
            vl=vsetvl_e16m2(avl);
            vb=vle16_v_u16m2(src1_addr,vl);
            vb=vand_vx_u16m2(vb,mask12,vl);
            vc=vle16_v_u16m2(src2_addr,vl);
            vc=vand_vx_u16m2(vc,~mask12,vl);
            vb=vxor_vv_u16m2(vb,vc,vl);
            vse16_v_u16m2(dst_addr,vb,vl);
            src1_addr+=vl,src2_addr+=vl,dst_addr+=vl,avl-=vl;
        }

        deg_sigma_p ^= mask12 & (deg_sigma_copy ^ deg_sigma_p);
        d = syndromes[mu + 1];

        for (i = 1; (i <= mu+1) && (i <= PARAM_DELTA); ++i) {
            d ^= gfmul_u16(sigma[i], syndromes[mu + 1 - i]);
        }
    }

    return deg_sigma;
}

static void compute_roots_custom(uint8_t *error, uint16_t *sigma) {
    uint16_t w[1 << PARAM_M] = {0};

    fft_custom(w, sigma, PARAM_DELTA + 1);
    fft_retrieve_error_poly_custom(error, w);
}

static void compute_z_poly_custom(uint16_t *z, const uint16_t *sigma, const uint16_t degree, const uint16_t *syndromes) {
    size_t i, j,vl,avl;
    uint16_t mask;
    vuint16m1_t va,vb,vc;

    z[0] = 1;

    avl=PARAM_DELTA;
    uint16_t* z_addr=z+1;
    const uint16_t* sigma_addr=sigma+1;
    csr_primpoly_rw(PARAM_GF_POLY);
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        va=vid_v_u16m1(vl);
        va=vadd_vx_u16m1(va,PARAM_DELTA-avl,vl);
        va=vsub_vx_u16m1(va,degree,vl);
        va=vsrl_vx_u16m1(va,15,vl);
        va=vrsub_vx_u16m1(va,0,vl);
        vb=vle16_v_u16m1(sigma_addr,vl);
        vb=vand_vv_u16m1(vb,va,vl);
        vse16_v_u16m1(z_addr,vb,vl);
        z_addr+=vl,sigma_addr+=vl,avl-=vl;
    }

    z[1] ^= syndromes[0];

    for (i = 2; i <= PARAM_DELTA; ++i) {
        mask = -((uint16_t) (i - degree - 1) >> 15);
        z[i] ^= mask & syndromes[i - 1];

        avl=i-1;
        vuint16m1_t vidx;
        vl=vsetvl_e16m1(1);
        va=vmv_s_x_u16m1(va,z[i],vl);
        const uint16_t* sigma_addr=sigma+1;
        const uint16_t* syndromes_addr=syndromes;
        while(avl>0){
            vl=vsetvl_e16m1(avl);
            vidx=vid_v_u16m1(vl);
            vidx=vrsub_vx_u16m1(vidx,avl-1,vl);
            vidx=vsll_vx_u16m1(vidx,1,vl);
            vb=vle16_v_u16m1(sigma_addr,vl);
            vc=vluxei16_v_u16m1(syndromes_addr,vidx,vl);
            vb=vgfmul_vv_u16m1(vb,vc);
            vb=vand_vx_u16m1(vb,mask,vl);
            va=vredxor_vs_u16m1_u16m1(va,vb,va,vl);
            sigma_addr+=vl,avl-=vl;
        }
        z[i]=vmv_x_s_u16m1_u16(va);
    }
}


static void compute_error_values_custom(uint16_t *error_values, const uint16_t *z, const uint8_t *error) {
    uint16_t beta_j[PARAM_DELTA] = {0};
    uint16_t e_j[PARAM_DELTA] = {0};

    uint16_t delta_counter;
    uint16_t delta_real_value;
    uint16_t found;
    uint16_t mask1;
    uint16_t mask2;
    uint16_t tmp1;
    uint16_t tmp2;
    uint16_t inverse;
    uint16_t inverse_power_j;

    size_t vl,avl;
    vuint16m1_t va,vb,vc,vd;
    vint16m1_t ve;
    csr_primpoly_rw(PARAM_GF_POLY);
    // Compute the beta_{j_i} page 31 of the documentation
    delta_counter = 0;
    for (size_t i = 0; i < PARAM_N1; i++) {
        mask1 = (uint16_t) (-((int32_t)error[i]) >> 31); // error[i] != 0
        
        avl=PARAM_DELTA;
        vl=vsetvl_e16m1(1);
        vc=vmv_v_x_u16m1(0,vl);//found
        uint16_t* beta_j_addr=beta_j;

        while(avl>0){
            vl=vsetvl_e16m1(avl);
            va=vid_v_u16m1(vl);
            va=vadd_vx_u16m1(va,PARAM_DELTA-avl,vl);
            va=vxor_vx_u16m1(va,delta_counter,vl);
            va=vrsub_vx_u16m1(va,0,vl);
            ve=vreinterpret_v_u16m1_i16m1(va);
            ve=vsra_vx_i16m1(ve,31,vl);
            va=vreinterpret_v_i16m1_u16m1(ve);
            va=vnot_v_u16m1(va,vl);
            va=vand_vx_u16m1(va,mask1,vl);//mask1&mask2
            vb=vand_vx_u16m1(va,gf_exp[i],vl);
            vd=vle16_v_u16m1(beta_j_addr,vl);
            vd=vadd_vv_u16m1(vd,vb,vl);
            vse16_v_u16m1(beta_j_addr,vd,vl);
            va=vand_vx_u16m1(va,1,vl);
            vc=vredsum_vs_u16m1_u16m1(vc,va,vc,vl);
            beta_j_addr+=vl,avl-=vl;
        }

        delta_counter += vmv_x_s_u16m1_u16(vc);
    }
    delta_real_value = delta_counter;

    // Compute the e_{j_i} page 31 of the documentation
    for (size_t i = 0; i < PARAM_DELTA; ++i) {
        tmp1 = 1;
        tmp2 = 1;
        inverse = gfinv_u16(beta_j[i]);
        inverse_power_j = 1;

        vl=vsetvl_e16m1(1);
        va=vmv_s_x_u16m1(va,inverse_power_j,vl);
        vb=vmv_s_x_u16m1(vb,tmp1,vl);
        for (size_t j = 1; j <= PARAM_DELTA; ++j) {
            va=vgfmul_vx_u16m1(va,inverse);
            vc=vgfmul_vx_u16m1(va,z[j]);
            vb=vxor_vv_u16m1(vb,vc,vl);
        }
        tmp1=vmv_x_s_u16m1_u16(vb);

        va=vmv_s_x_u16m1(va,inverse,vl);
        vb=vmv_s_x_u16m1(vb,tmp2,vl);
        for (size_t k = 1; k < PARAM_DELTA; ++k) {
            vc=vgfmul_vx_u16m1(va,beta_j[(i + k) % PARAM_DELTA]);
            vc=vxor_vx_u16m1(vc,1,vl);
            vb=vgfmul_vv_u16m1(vb,vc);
        }
        tmp2=vmv_x_s_u16m1_u16(vb);

        mask1 = (uint16_t) (((int16_t) i - delta_real_value) >> 15); // i < delta_real_value

        uint16_t tmp2_inv=gfinv_u16(tmp2);
        va=vmv_s_x_u16m1(va,tmp1,vl);
        va=vgfmul_vx_u16m1(va,tmp2_inv);
        uint16_t tmp=vmv_x_s_u16m1_u16(va);
        e_j[i] = mask1 & tmp;
    }

    // Place the delta e_{j_i} values at the right coordinates of the output vector
    delta_counter = 0;
    for (size_t i = 0; i < PARAM_N1; ++i) {
        found = 0;
        mask1 = (uint16_t) (-((int32_t)error[i]) >> 31); // error[i] != 0

        avl=PARAM_DELTA;
        vl=vsetvl_e16m1(1);
        vc=vmv_v_x_u16m1(0,vl);//found
        uint16_t* e_j_addr=e_j;
        uint16_t* error_values_addr=error_values;
        while(avl>0){
            vl=vsetvl_e16m1(avl);
            va=vid_v_u16m1(vl);
            va=vadd_vx_u16m1(va,PARAM_DELTA-avl,vl);
            va=vxor_vx_u16m1(va,delta_counter,vl);
            va=vrsub_vx_u16m1(va,0,vl);
            ve=vreinterpret_v_u16m1_i16m1(va);
            ve=vsra_vx_i16m1(ve,31,vl);
            va=vreinterpret_v_i16m1_u16m1(ve);
            va=vnot_v_u16m1(va,vl);
            va=vand_vx_u16m1(va,mask1,vl);//mask1&mask2
            vb=vle16_v_u16m1(e_j_addr,vl);//e_j[j]
            vb=vand_vv_u16m1(vb,va,vl);
            vd=vle16_v_u16m1(error_values_addr,vl);
            vd=vadd_vv_u16m1(vd,vb,vl);
            vse16_v_u16m1(error_values_addr,vd,vl);
            va=vand_vx_u16m1(va,1,vl);
            vc=vredsum_vs_u16m1_u16m1(vc,va,vc,vl);
            e_j_addr+=vl,error_values_addr+=vl,avl-=vl;
        }

        delta_counter += vmv_x_s_u16m1_u16(vc);
    }
}

static void correct_errors_custom(uint8_t *cdw, const uint16_t *error_values) {
    size_t vl,avl;
    vuint16m4_t va;
    vuint8m2_t vb,vc;
    avl=PARAM_N1;
    while(avl>0){
        vl=vsetvl_e16m4(avl);
        va=vle16_v_u16m4(error_values,vl);
        vb=vnsrl_wx_u8m2(va,0,vl);
        vl=vsetvl_e8m2(avl);
        vc=vle8_v_u8m2(cdw,vl);
        vc=vxor_vv_u8m2(vc,vb,vl);
        vse8_v_u8m2(cdw,vc,vl);
        cdw+=vl,error_values+=vl,avl-=vl;
    }
}

void reed_solomon_decode_custom(uint64_t *msg, uint64_t *cdw) {
    uint8_t cdw_bytes[PARAM_N1] = {0};
    uint16_t syndromes[2 * PARAM_DELTA] = {0};
    uint16_t sigma[1 << PARAM_FFT] = {0};
    uint8_t error[1 << PARAM_M] = {0};
    uint16_t z[PARAM_N1] = {0};
    uint16_t error_values[PARAM_N1] = {0};
    uint16_t deg;

    // Copy the vector in an array of bytes
    memcpy(cdw_bytes, cdw, PARAM_N1);

    // Calculate the 2*PARAM_DELTA syndromes
    compute_syndromes_custom(syndromes, cdw_bytes);

    // Compute the error locator polynomial sigma
    // Sigma's degree is at most PARAM_DELTA but the FFT requires the extra room
    deg = compute_elp_custom(sigma, syndromes);

    // Compute the error polynomial error
    compute_roots_custom(error, sigma);

    // Compute the polynomial z(x)
    compute_z_poly_custom(z, sigma, deg, syndromes);

    // Compute the error values
    compute_error_values_custom(error_values, z, error);

    // Correct the errors
    correct_errors_custom(cdw_bytes, error_values);

    // Retrieve the message from the decoded codeword
    memcpy(msg, cdw_bytes + (PARAM_G - 1) , PARAM_K);

    #ifdef VERBOSE
        printf("\n\nThe syndromes: ");
        for (size_t i = 0 ; i < 2*PARAM_DELTA ; ++i) {
            printf("%u ", syndromes[i]);
        }
        printf("\n\nThe error locator polynomial: sigma(x) = ");
        bool first_coeff = true;
        if (sigma[0]) {
            printf("%u", sigma[0]);
            first_coeff = false;
        }
        for (size_t i = 1 ; i < (1 << PARAM_FFT) ; ++i) {
            if (sigma[i] == 0)
                continue;
            if (!first_coeff)
                printf(" + ");
            first_coeff = false;
            if(sigma[i] != 1)
                printf("%u ", sigma[i]);
            if (i == 1)
                printf("x");
            else
                printf("x^%zu", i);
        }
        if (first_coeff)
            printf("0");

        printf("\n\nThe polynomial: z(x) = ");
        bool first_coeff_1 = true;
        if (z[0]) {
            printf("%u", z[0]);
            first_coeff_1 = false;
        }
        for (size_t i = 1 ; i < (PARAM_DELTA + 1) ; ++i) {
            if (z[i] == 0)
                continue;
            if (!first_coeff_1)
                printf(" + ");
            first_coeff_1 = false;
            if(z[i] != 1)
                printf("%u ", z[i]);
            if (i == 1)
                printf("x");
            else
                printf("x^%zu", i);
        }
        if (first_coeff_1)
            printf("0");

        printf("\n\nThe pairs of (error locator numbers, error values): ");
        size_t j = 0;
        for (size_t i = 0 ; i < PARAM_N1 ; ++i) {
            if(error[i]){
                printf("(%zu, %d) ", i, error_values[j]);
                j++;
            }
        }
        printf("\n");
    #endif
}

///////////////////////////////////////////////////////////////////////////
//                               Tests
///////////////////////////////////////////////////////////////////////////
bool test_compute_syndromes_custom(){
    uint8_t cdw_bytes[PARAM_N1] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<PARAM_N1;i++){
        cdw_bytes[i]=rand()&255;
    }

    uint16_t syndromes_ref[2 * PARAM_DELTA] = {0};
    uint16_t syndromes_res[2 * PARAM_DELTA] = {0};

    uint64_t start, end;

    start=read_cycle();
    compute_syndromes(syndromes_ref,cdw_bytes);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    compute_syndromes_custom(syndromes_res,cdw_bytes);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(2*PARAM_DELTA);i++){
        if(syndromes_ref[i]!=syndromes_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_compute_syndromes_custom return with flag %d\n",flag);

    return flag;
}

bool test_compute_elp_custom(){
    uint16_t syndromes[2 * PARAM_DELTA] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<(2*PARAM_DELTA);i++){
        syndromes[i]=rand()&255;
    }

    uint16_t sigma_ref[1 << PARAM_FFT] = {0};
    uint16_t sigma_res[1 << PARAM_FFT] = {0};
    uint16_t deg_ref,deg_res;

    uint64_t start, end;

    start=read_cycle();
    deg_ref=compute_elp(sigma_ref,syndromes);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    deg_res=compute_elp_custom(sigma_res,syndromes);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(1 << PARAM_FFT);i++){
        if(sigma_ref[i]!=sigma_res[i]){
            flag=false;
            break;
        }
    }
    printf("%d\n",flag);

    if(deg_ref!=deg_res)flag=false;

    printf("test_compute_elp_custom return with flag %d\n",flag);

    return flag;
}

bool test_compute_roots_custom(){
    uint16_t sigma[1 << PARAM_FFT] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<(1<<PARAM_FFT);i++){
        sigma[i]=rand()&255;
    }

    uint8_t error_ref[1 << PARAM_M] = {0};
    uint8_t error_res[1 << PARAM_M] = {0};

    uint64_t start, end;

    start=read_cycle();
    compute_roots(error_ref,sigma);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    compute_roots_custom(error_res,sigma);
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

    printf("test_compute_roots_custom return with flag %d\n",flag);

    return flag;
}

bool test_compute_z_poly_custom(){
    uint16_t syndromes[2 * PARAM_DELTA]={
        62,115,221,54,53,107,23,241,86,58,80,134,168,87,251,143,
        47,62,115,221,54,53,107,23,241,86,58,80,134,168
    };
    uint16_t sigma[1 << PARAM_FFT] ={
        1,193,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    uint16_t deg=1;
    uint16_t z_ref[PARAM_N1] = {
        1,255,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
    uint16_t z_res[PARAM_N1] = {0};

    compute_z_poly_custom(z_res, sigma, deg, syndromes);

    //CHECK
    bool flag=true;
    for(int i=0;i<PARAM_N1;i++){
        if(z_ref[i]!=z_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_compute_z_poly_custom return with flag %d\n",flag);

    return flag;
}

bool test_compute_error_values_custom(){
    uint8_t error[1 << PARAM_M] = {0};
    uint16_t z[PARAM_N1] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<(1<<PARAM_M);i++){
        error[i]=rand()&255;
    }
    srand((unsigned)time(NULL));
    for(int i=0;i<(PARAM_N1);i++){
        z[i]=rand()&255;
    }

    uint16_t error_values_ref[PARAM_N1] = {0};
    uint16_t error_values_res[PARAM_N1] = {0};

    uint64_t start, end;

    start=read_cycle();
    compute_error_values(error_values_ref, z, error);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    compute_error_values_custom(error_values_res, z, error);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(PARAM_N1);i++){
        if(error_values_ref[i]!=error_values_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_compute_error_values_custom return with flag %d\n",flag);

    return flag;
}

bool test_correct_errors_custom(){
    uint8_t cdw_bytes_ref[PARAM_N1] = {0};
    uint8_t cdw_bytes_res[PARAM_N1] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<(PARAM_N1);i++){
        cdw_bytes_ref[i]=cdw_bytes_res[i]=rand()&255;
    }
    uint16_t error_values[PARAM_N1] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<(PARAM_N1);i++){
        error_values[i]=rand()&255;
    }

    uint64_t start, end;

    start=read_cycle();
    correct_errors(cdw_bytes_ref, error_values);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    correct_errors_custom(cdw_bytes_res, error_values);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(PARAM_N1);i++){
        if(cdw_bytes_ref[i]!=cdw_bytes_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_correct_errors_custom return with flag %d\n",flag);

    return flag;
}

// int main(){
//     // test_compute_syndromes_custom();
//     // test_compute_elp_custom();
//     // test_compute_roots_custom();
//     // test_compute_z_poly_custom();
//     // test_compute_error_values_custom();
//     test_correct_errors_custom();

//     return 0;
// }