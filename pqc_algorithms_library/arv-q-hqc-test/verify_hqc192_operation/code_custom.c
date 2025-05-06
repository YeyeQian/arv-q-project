#include "code.h"
#include "reed_muller.h"
#include "reed_solomon.h"
#include "parameters.h"
#include <stdint.h>
#include <string.h>
#ifdef VERBOSE
#include <stdio.h>
#include "vector.h"
#endif

void code_encode_custom(uint64_t *em, const uint64_t *m) {
    uint64_t tmp[VEC_N1_SIZE_64] = {0};

    reed_solomon_encode_custom(tmp, m);
    reed_muller_encode_custom(em, tmp);

    #ifdef VERBOSE
        printf("\n\nReed-Solomon code word: "); vect_print(tmp, VEC_N1_SIZE_BYTES);
        printf("\n\nConcatenated code word: "); vect_print(em, VEC_N1N2_SIZE_BYTES);
    #endif
}

void code_decode_custom(uint64_t *m, const uint64_t *em) {
    uint64_t tmp[VEC_N1_SIZE_64] = {0};

    reed_muller_decode_custom(tmp, em);
    reed_solomon_decode_custom(m, tmp);


    #ifdef VERBOSE
        printf("\n\nReed-Muller decoding result (the input for the Reed-Solomon decoding algorithm): "); vect_print(tmp, VEC_N1_SIZE_BYTES);
    #endif
}