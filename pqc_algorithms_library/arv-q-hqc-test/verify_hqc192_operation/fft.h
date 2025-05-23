#ifndef FFT_H
#define FFT_H

/**
 * @file fft.h
 * @brief Header file of fft.cpp
 */


#include <stddef.h>
#include <stdint.h>

void fft(uint16_t *w, const uint16_t *f, size_t f_coeffs);
void fft_retrieve_error_poly(uint8_t *error, const uint16_t *w);

//Customized Versions
void fft_custom(uint16_t *w, const uint16_t *f, size_t f_coeffs);
void fft_retrieve_error_poly_custom(uint8_t *error, const uint16_t *w);
#endif
