#include <riscv_vector.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include "inner.h"
#include "params_custom.h"
#include "../apis/custom_inst_api.h"

#define Np 4

void gen_rand_fpr_array(fpr* a, fpr* b, size_t len) {
    srand((unsigned int)time(NULL));
    for (size_t i = 0; i < len; ++i)
    {
        double random = ((double) rand()) /(double) RAND_MAX;
        double range = 150.9428f - 10.5231f;
        double num = (random * range) + 10.5231f;
        a[i].v = num;
        b[i].v = num;
    }
    for (size_t i = 0; i < len; i += 2) {
        a[i].v = -a[i].v;
        b[i].v = -b[i].v;
    }
}

double comparefpr_1d(fpr* a, fpr* b, size_t len) {
	double error = 0;
    for (size_t i = 0; i < len; ++i)
        if (a[i].v != b[i].v) {
            double error_cur = a[i].v > b[i].v ? (a[i].v - b[i].v) : (b[i].v - a[i].v) ;
			if (error_cur > error)
				error = error_cur;
		}
    return error;
}

#define ERROR_LIMIT 0.001f

bool test_FFT_custom()
{
    fpr orig_a[FALCON_N], cust_a[FALCON_N];
    gen_rand_fpr_array(orig_a, cust_a, FALCON_N);

    Zf(FFT)(orig_a, LOGN);

    Zf(FFT_custom)(cust_a, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_a[%d] = %28.27f \n", i, cust_a[i].v);
        printf("orig_a[%d] = %28.27f \n", i, orig_a[i].v);
    }

    double error = comparefpr_1d(cust_a, orig_a, FALCON_N);

    printf("FFT execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_FFT_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_FFT_custom fail with result=%d!\n", result);
    }

    return result;
}

bool test_iFFT_custom()
{
    fpr orig_a[FALCON_N], cust_a[FALCON_N];
    gen_rand_fpr_array(orig_a, cust_a, FALCON_N);

    Zf(FFT)(orig_a, LOGN);

    Zf(FFT_custom)(cust_a, LOGN);

    Zf(iFFT)(orig_a, LOGN);
    Zf(iFFT_custom)(cust_a, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_a[%d] = %28.27f \n", i, cust_a[i].v);
        printf("orig_a[%d] = %28.27f \n", i, orig_a[i].v);
    }

    double error = comparefpr_1d(cust_a, orig_a, FALCON_N);

    printf("iFFT execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_iFFT_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_iFFT_custom fail with result=%d!\n", result);
    }

    return result;
}

bool test_poly_add_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);
    fpr cust_op2[FALCON_N], orig_op2[FALCON_N];
    gen_rand_fpr_array(cust_op2, orig_op2, FALCON_N);

    Zf(poly_add)(orig_op1, orig_op2, LOGN); 
    Zf(poly_add_custom)(cust_op1, cust_op2, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_add execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_add_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_add_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_sub_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);
    fpr cust_op2[FALCON_N], orig_op2[FALCON_N];
    gen_rand_fpr_array(cust_op2, orig_op2, FALCON_N);

    Zf(poly_sub)(orig_op1, orig_op2, LOGN); 
    Zf(poly_sub_custom)(cust_op1, cust_op2, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_sub execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_add_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_add_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_neg_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);

    Zf(poly_neg)(orig_op1, LOGN); 
    Zf(poly_neg_custom)(cust_op1, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_neg execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_neg_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_neg_custom fail with result=%d!\n", result);
    }

    return result;    
}


bool test_poly_adj_fft_custom()
{
    fpr orig_a[FALCON_N], cust_a[FALCON_N];
    gen_rand_fpr_array(orig_a, cust_a, FALCON_N);

    Zf(poly_adj_fft)(orig_a, LOGN);

    Zf(poly_adj_fft_custom)(cust_a, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_a[%d] = %28.27f \n", i, cust_a[i].v);
        printf("orig_a[%d] = %28.27f \n", i, orig_a[i].v);
    }

    double error = comparefpr_1d(cust_a, orig_a, FALCON_N);

    printf("poly_adj_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_adj_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_adj_fft_custom fail with result=%d!\n", result);
    }

    return result;
}

bool test_poly_mul_fft_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);
    fpr cust_op2[FALCON_N], orig_op2[FALCON_N];
    gen_rand_fpr_array(cust_op2, orig_op2, FALCON_N);

    Zf(poly_mul_fft)(orig_op1, orig_op2, LOGN); 
    Zf(poly_mul_fft_custom)(cust_op1, cust_op2, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_mul_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_mul_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_mul_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}


bool test_poly_muladj_fft_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);
    fpr cust_op2[FALCON_N], orig_op2[FALCON_N];
    gen_rand_fpr_array(cust_op2, orig_op2, FALCON_N);

    Zf(poly_muladj_fft)(orig_op1, orig_op2, LOGN); 
    Zf(poly_muladj_fft_custom)(cust_op1, cust_op2, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_muladj_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_muladj_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_muladj_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_mulconst_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);
    fpr x = orig_op1[0];

    Zf(poly_mulconst)(orig_op1, x, LOGN); 
    Zf(poly_mulconst_custom)(cust_op1, x, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_mulconst execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_mulconst_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_mulconst_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_div_fft_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);
    fpr cust_op2[FALCON_N], orig_op2[FALCON_N];
    gen_rand_fpr_array(cust_op2, orig_op2, FALCON_N);

    Zf(poly_div_fft)(orig_op1, orig_op2, LOGN); 
    Zf(poly_div_fft_custom)(cust_op1, cust_op2, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_div_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_div_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_div_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_invnorm2_fft_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);
    fpr cust_op2[FALCON_N], orig_op2[FALCON_N];
    gen_rand_fpr_array(cust_op2, orig_op2, FALCON_N);

    fpr cust_opd[FALCON_N], orig_opd[FALCON_N];

    Zf(poly_invnorm2_fft)(orig_opd, orig_op1, orig_op2, LOGN); 
    Zf(poly_invnorm2_fft_custom)(cust_opd, cust_op1, cust_op2, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_invnorm2_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_invnorm2_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_invnorm2_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}


bool test_poly_add_muladj_fft_custom()
{
    fpr cust_F[FALCON_N], orig_F[FALCON_N];
    gen_rand_fpr_array(cust_F, orig_F, FALCON_N);
    fpr cust_G[FALCON_N], orig_G[FALCON_N];
    gen_rand_fpr_array(cust_G, orig_G, FALCON_N);
    fpr cust_f[FALCON_N], orig_f[FALCON_N];
    gen_rand_fpr_array(cust_f, orig_f, FALCON_N);
    fpr cust_g[FALCON_N], orig_g[FALCON_N];
    gen_rand_fpr_array(cust_g, orig_g, FALCON_N);

    fpr cust_opd[FALCON_N], orig_opd[FALCON_N];

    Zf(poly_add_muladj_fft)(orig_opd, orig_F, orig_G, orig_f, orig_g, LOGN);

    Zf(poly_add_muladj_fft_custom)(cust_opd, cust_F, cust_G, cust_f, cust_g, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_opd[%d] = %28.27f \n", i, cust_opd[i].v);
        printf("orig_opd[%d] = %28.27f \n", i, orig_opd[i].v);
    }

    double error = comparefpr_1d(cust_opd, orig_opd, FALCON_N);

    printf("poly_add_muladj_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_add_muladj_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_add_muladj_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_mul_autoadj_fft_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);
    fpr cust_op2[FALCON_N], orig_op2[FALCON_N];
    gen_rand_fpr_array(cust_op2, orig_op2, FALCON_N);

    Zf(poly_mul_autoadj_fft)(orig_op1, orig_op2, LOGN); 
    Zf(poly_mul_autoadj_fft_custom)(cust_op1, cust_op2, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_mul_autoadj_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_mul_autoadj_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_mul_autoadj_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_div_autoadj_fft_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];
    gen_rand_fpr_array(cust_op1, orig_op1, FALCON_N);
    fpr cust_op2[FALCON_N], orig_op2[FALCON_N];
    gen_rand_fpr_array(cust_op2, orig_op2, FALCON_N);

    Zf(poly_div_autoadj_fft)(orig_op1, orig_op2, LOGN); 
    Zf(poly_div_autoadj_fft_custom)(cust_op1, cust_op2, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_op1[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_op1[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_div_autoadj_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_div_autoadj_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_div_autoadj_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_LDLmv_fft_custom()
{
    fpr cust_F[FALCON_N], orig_F[FALCON_N];
    gen_rand_fpr_array(cust_F, orig_F, FALCON_N);
    fpr cust_G[FALCON_N], orig_G[FALCON_N];
    gen_rand_fpr_array(cust_G, orig_G, FALCON_N);
    fpr cust_f[FALCON_N], orig_f[FALCON_N];
    gen_rand_fpr_array(cust_f, orig_f, FALCON_N);
    fpr cust_g[FALCON_N], orig_g[FALCON_N];
    gen_rand_fpr_array(cust_g, orig_g, FALCON_N);

    fpr cust_opd[FALCON_N], orig_opd[FALCON_N];

    Zf(poly_LDLmv_fft)(orig_opd, orig_F, orig_G, orig_f, orig_g, LOGN);

    Zf(poly_LDLmv_fft_custom)(cust_opd, cust_F, cust_G, cust_f, cust_g, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_opd[%d] = %28.27f \n", i, cust_opd[i].v);
        printf("orig_opd[%d] = %28.27f \n", i, orig_opd[i].v);
    }

    double error = comparefpr_1d(cust_opd, orig_opd, FALCON_N);

    printf("poly_LDLmv_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_LDLmv_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_LDLmv_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_LDL_fft_custom()
{
    fpr cust_F[FALCON_N], orig_F[FALCON_N];
    gen_rand_fpr_array(cust_F, orig_F, FALCON_N);
    fpr cust_G[FALCON_N], orig_G[FALCON_N];
    gen_rand_fpr_array(cust_G, orig_G, FALCON_N);
    fpr cust_f[FALCON_N], orig_f[FALCON_N];
    gen_rand_fpr_array(cust_f, orig_f, FALCON_N);
    fpr cust_g[FALCON_N], orig_g[FALCON_N];
    gen_rand_fpr_array(cust_g, orig_g, FALCON_N);


    Zf(poly_LDL_fft)(orig_G, orig_f, orig_g, LOGN);
    Zf(poly_LDL_fft_custom)(cust_G, cust_f, cust_g, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_g11[%d] = %28.27f \n", i, cust_g[i].v);
        printf("orig_g11[%d] = %28.27f \n", i, orig_g[i].v);
    }

    double error = comparefpr_1d(cust_g, orig_g, FALCON_N);

    printf("poly_LDL_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_LDL_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_LDL_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_split_fft_custom()
{
    fpr cust_G[FALCON_N], orig_G[FALCON_N];
    gen_rand_fpr_array(cust_G, orig_G, FALCON_N);

    fpr cust_f0[FALCON_N >> 1], orig_f0[FALCON_N >> 1];
    fpr cust_f1[FALCON_N >> 1], orig_f1[FALCON_N >> 1];


    Zf(poly_split_fft)(orig_f0, orig_f1, orig_G, LOGN);
    Zf(poly_split_fft_custom)(cust_f0, cust_f1, cust_G, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_g11[%d] = %28.27f \n", i, cust_f0[i].v);
        printf("orig_g11[%d] = %28.27f \n", i, orig_f0[i].v);
    }

    double error = comparefpr_1d(cust_f0, orig_f0, FALCON_N >> 1);

    printf("poly_split_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_split_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_split_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}

bool test_poly_merge_fft_custom()
{
    fpr cust_op1[FALCON_N], orig_op1[FALCON_N];

    fpr cust_f0[FALCON_N >> 1], orig_f0[FALCON_N >> 1];
    gen_rand_fpr_array(cust_f0, orig_f0, FALCON_N >> 1);
    fpr cust_f1[FALCON_N >> 1], orig_f1[FALCON_N >> 1];
    gen_rand_fpr_array(cust_f1, orig_f1, FALCON_N >> 1);

    Zf(poly_merge_fft)(orig_op1, orig_f0, orig_f1, LOGN);
    Zf(poly_merge_fft_custom)(cust_op1, cust_f0, cust_f1, LOGN);

    for(uint32_t i = 0; i < Np; i++) {
        printf("cust_g11[%d] = %28.27f \n", i, cust_op1[i].v);
        printf("orig_g11[%d] = %28.27f \n", i, orig_op1[i].v);
    }

    double error = comparefpr_1d(cust_op1, orig_op1, FALCON_N);

    printf("poly_merge_fft execution error = %28.27f \n", error);

    bool result = (error < ERROR_LIMIT);

    if(result) {
        printf("test_poly_merge_fft_custom pass with result=%d!\n", result);
    }
    else {
        printf("test_poly_merge_fft_custom fail with result=%d!\n", result);
    }

    return result;    
}

int main()
{
    // test_FFT_custom();
    // test_iFFT_custom();
    // test_poly_add_custom();
    // test_poly_neg_custom();
    // test_poly_adj_fft_custom();

    // test_poly_mul_fft_custom();
    // test_poly_muladj_fft_custom();
    // test_poly_mulconst_custom();

    // test_poly_div_fft_custom();
    // test_poly_invnorm2_fft_custom();

    // test_poly_add_muladj_fft_custom();
    // test_poly_mul_autoadj_fft_custom();
    // test_poly_div_autoadj_fft_custom();
    // test_poly_LDLmv_fft_custom();

    // test_poly_LDL_fft_custom();
    test_poly_split_fft_custom();
    // test_poly_merge_fft_custom();

    return 0;
}