#include <stdlib.h>
#include <stdio.h>
#include "extGCD_polyinv64.h"

void clmul64(uint64_t a, uint64_t b, uint64_t res[2]) {
    uint64_t temp[64][2] = { 0 };
    for (int i = 0; i < 64; i++) {
        if ((b & 1) == 1) {
            temp[i][0] = i == 0 ? a : a << i;
            temp[i][1] = i==0?0: a >> (64 - i);
        }
        b >>= 1;
    }

    for (int j = 32; j > 0; j >>= 1) {
        for (int i = 0; i < j; i++) {
            temp[i][0] ^= temp[i + j][0];
            temp[i][1] ^= temp[i + j][1];
        }
    }
    res[0] = temp[0][0];
    res[1] = temp[0][1];
}

void clmul64_array(uint64_t* din, int len, uint64_t h, uint64_t* dout) {
    uint64_t temp[2] = { 0 };
    clmul64(din[0],h,temp);
    dout[0] = temp[0];
    uint64_t temp1 = temp[1];
    for (int i = 1; i < len; i++) {
        clmul64(din[i],h,temp);
        dout[i] = temp1 ^ temp[0];
        temp1 = temp[1];
    }
}


void gen_hmatrix64(int32_t* delta_ptr, uint64_t fm, uint64_t gm, uint32_t step_size, uint64_t h[4]) {
    for (int j = 1; j <= step_size; j++) {
        uint64_t g_msb = gm >> 63;
        if (g_msb == 0) {
            (*delta_ptr) += 1;
            gm <<= 1;
            h[2] <<= 1, h[3] <<= 1;
        }
        else {
            if ((*delta_ptr) > 0) {
                (*delta_ptr) = -(*delta_ptr) + 1;
                uint64_t temp = fm;
                fm = gm, gm = (temp ^ gm) << 1;
                uint64_t h0_temp = h[0];
                uint64_t h1_temp = h[1];
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

#define R_SIZE_PRIME_HALF 193
void extGCD_polyinv64(uint32_t gin[R_SIZE_PRIME], uint32_t gout[R_SIZE_PRIME]) {
    uint8_t shiftcnt = 63 - (R_BITS & 63);
    uint64_t* f = (uint64_t*)calloc(R_SIZE_PRIME_HALF, sizeof(uint64_t));
    f[0] = 1 << shiftcnt, f[R_SIZE_PRIME_HALF - 1] = 1ULL << 63;

    uint32_t* g_temp = (uint32_t*)calloc(R_SIZE_PRIME,sizeof(uint32_t));
    for (int i = R_SIZE_PRIME - 1; i > 0; i--) {
        g_temp[i] = (((uint64_t)gin[i] << SEW) | (uint64_t)gin[i - 1]) >> (1 + (R_BITS & (SEW - 1)));
    }
    g_temp[0] = gin[0] << (SEW - 1 - (R_BITS & (SEW - 1)));
    uint64_t* g = (uint64_t*)calloc(R_SIZE_PRIME_HALF, sizeof(uint64_t));
    for (int i = 0; i < R_SIZE_PRIME; i += 2)g[i >> 1] = ((uint64_t)g_temp[i + 1] << 32) | (uint64_t)g_temp[i];

    uint64_t w[R_SIZE_PRIME] = { 0 };
    w[0] = 1;
    uint64_t v[R_SIZE_PRIME] = { 0 };

    int32_t delta = 0;
    uint32_t Tau = 2 * R_BITS - 1;

    uint32_t step_size = 64 - 1;

    //define some temporary stack space
    uint64_t* h0fv = (uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));
    uint64_t* h1gw = (uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));
    uint64_t* h2fv = (uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));
    uint64_t* h3gw = (uint64_t*)calloc(R_SIZE_PRIME, sizeof(uint64_t));

    while (Tau >= step_size) {
        uint64_t h[4] = { 1,0,0,1 };
        uint64_t gm = g[R_SIZE_PRIME_HALF - 1];
        uint64_t fm = f[R_SIZE_PRIME_HALF - 1];
        gen_hmatrix64(&delta, fm, gm, step_size, h);
        //do carry-less multiplication
        clmul64_array(f,R_SIZE_PRIME_HALF,h[0],h0fv);
        clmul64_array(g,R_SIZE_PRIME_HALF,h[1],h1gw);
        clmul64_array(f,R_SIZE_PRIME_HALF,h[2],h2fv);
        clmul64_array(g,R_SIZE_PRIME_HALF,h[3],h3gw);

        //update
        for (int i = 0; i < R_SIZE_PRIME_HALF; i++) {
            f[i] = h0fv[i] ^ h1gw[i];
            g[i] = h2fv[i] ^ h3gw[i];
        }

        //do carry-less multiplication
        clmul64_array(v, R_SIZE_PRIME, h[0], h0fv);
        clmul64_array(w, R_SIZE_PRIME, h[1], h1gw);
        clmul64_array(v, R_SIZE_PRIME, h[2], h2fv);
        clmul64_array(w, R_SIZE_PRIME, h[3], h3gw);

        //update
        for (int i = 0; i < R_SIZE_PRIME; i++) {
            v[i] = h0fv[i] ^ h1gw[i];
            w[i] = h2fv[i] ^ h3gw[i];
        }

        Tau -= step_size;
    }
    if (Tau > 0) {
        uint64_t h[4] = { 1,0,0,1 };
        uint64_t gm = g[R_SIZE_PRIME_HALF - 1];
        uint64_t fm = f[R_SIZE_PRIME_HALF - 1];
        gen_hmatrix64(&delta, fm, gm, Tau, h);

        //do carry-less multiplication
        clmul64_array(v, R_SIZE_PRIME, h[0], h0fv);
        clmul64_array(w, R_SIZE_PRIME, h[1], h1gw);
        clmul64_array(v, R_SIZE_PRIME, h[2], h2fv);
        clmul64_array(w, R_SIZE_PRIME, h[3], h3gw);

        //update
        for (int i = 0; i < R_SIZE_PRIME; i++) {
            v[i] = h0fv[i] ^ h1gw[i];
            w[i] = h2fv[i] ^ h3gw[i];
        }
    }
    //shift v R_BITS right

    //temp print
    printf("v:\n");
    for (int i = 0; i < R_SIZE_PRIME; i++) {
        printf("%lluULL,",v[i]);
        if ((i + 1) % 32 == 0)printf("\n");
    }

    uint32_t* v_temp = (uint32_t*)calloc(R_SIZE_PRIME<<1, sizeof(uint32_t));
    for (int i = 0; i < R_SIZE_PRIME; i++) {
        v_temp[2 * i] = v[i] & CHUNK_MASK;
        v_temp[2 * i + 1] = (v[i] >> 32) & CHUNK_MASK;
    }
    for (int i = 0; i < R_SIZE_PRIME - 1; i++) {
        gout[i] = (((uint64_t)v_temp[i + R_SIZE_PRIME] << SEW) | ((uint64_t)v_temp[i + R_SIZE_PRIME - 1])) >> (R_BITS & (SEW - 1));
    }
    gout[R_SIZE_PRIME - 1] = v_temp[2 * R_SIZE_PRIME - 2] >> (R_BITS & (SEW - 1));

    //free stack spaces
    free(f),free(g), free(g_temp), free(v_temp);
    free(h0fv), free(h1gw), free(h2fv), free(h3gw);
}