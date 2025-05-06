#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "api.h"
#include "compiler.h"
#include "fips202.h"
#include "symmetric.h"
#include "poly.h"
#include "param.h"
#include "rounding.h"
#include "reduce.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "../common/debug.h"
#include <math.h>

int test_poly_reduce_rvv()
{   
    poly a;
    poly b;
    bool flag;
    int i;
    int t=pow(2,31)-pow(2,22);
    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++) {
        b.coeffs[i]=a.coeffs[i] = rand() % t;
    }    

    poly_reduce(&a);

    poly_reduce_rvv(&b);

    flag = true;
    for(i = 0; i < N; i++) {
        if(a.coeffs[i] != b.coeffs[i]) {
            log_e("a.coeff[%d]=%d, while b.coeff[%d]=%d", i, a.coeffs[i], i, b.coeffs[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("poly_reduce_rvv test pass!\n");
    }
    else {
        printf("poly_reduce_rvv test fail..\n");
    }

    return flag;
}

int test_poly_make_hint_rvv()
{   
    poly a0,a1,h,h_ref;
    uint32_t s,s_ref;
    bool flag;
    int i;
    
    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++) {
        int32_t a=rand()%Q;
        int32_t a0_val,a1_val;
        a1_val=power2round(&a0_val,a);
        a0.coeffs[i]=a0_val;
        a1.coeffs[i]=a1_val;
    }    

    s_ref=poly_make_hint(&h_ref,&a0,&a1);

    s=poly_make_hint_rvv(&h,&a0,&a1);

    flag = true;
    for(i = 0; i < N; i++) {
        if(h.coeffs[i] != h_ref.coeffs[i]) {
            log_e("h.coeff[%d]=%d, while h_ref.coeff[%d]=%d", i, h.coeffs[i], i, h_ref.coeffs[i]);
            flag = false;
        }
    }

    if(s!=s_ref) flag=false;

    if(flag) {
        printf("poly_make_hint_rvv test pass!\n");
    }
    else {
        printf("poly_make_hint_rvv test fail..\n");
    }

    return flag;
}

int test_poly_decompose_rvv()
{   
    poly a,a1,a0,a0_ref,a1_ref;
    bool flag;
    int i;
    
    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++) {
        a.coeffs[i]=rand()%Q;
    }    

    poly_decompose(&a1_ref,&a0_ref,&a);

    poly_decompose_rvv(&a1,&a0,&a);

    flag = true;
    for(i = 0; i < N; i++) {
        if(a1.coeffs[i] != a1_ref.coeffs[i]||a0.coeffs[i] != a0_ref.coeffs[i]) {
            printf("a1_ref[%d]=%d,a1[%d]=%d,a0_ref[%d]=%d,a0[%d]=%d\n",i,a1_ref.coeffs[i],i,a1.coeffs[i],i,a0_ref.coeffs[i],i,a0.coeffs[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("poly_decompose_rvv test pass!\n");
    }
    else {
        printf("poly_decompose_rvv test fail..\n");
    }

    return flag;
}

bool test_poly_use_hint_rvv()
{   
    poly a,h,b,b_ref;
    bool flag;
    int i;
    
    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++) {
        a.coeffs[i]=rand()%Q;
        h.coeffs[i]=rand()&1;
    }    

    poly_use_hint(&b_ref,&a,&h);

    poly_use_hint_rvv(&b,&a,&h);

    flag = true;
    for(i = 0; i < N; i++) {
        if(b.coeffs[i] != b_ref.coeffs[i]) {
            log_e("b.coeff[%d]=%d, while b_ref.coeff[%d]=%d", i, b.coeffs[i], i, b_ref.coeffs[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("poly_use_hint_rvv test pass!\n");
    }
    else {
        printf("poly_use_hint_rvv test fail..\n");
    }

    return flag;
}

bool test_poly_chknorm_rvv()
{   
    poly a;
    bool flag;
    int r,ref;
    int i;
    int t=pow(2,31)-pow(2,22);
    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++) {
        a.coeffs[i]=reduce32(rand()%t);
    }    

    ref=poly_chknorm(&a,GAMMA1 - BETA);

    r=poly_chknorm_rvv(&a,GAMMA1 - BETA);

    flag = true;
    if(r!=ref) flag=false;

    if(flag) {
        printf("poly_chknorm_rvv test pass!\n");
    }
    else {
        printf("poly_chknorm_rvv test fail..\n");
    }

    return flag;
}

bool test_polyeta_pack_custom()
{   
    poly a;
    uint8_t r[POLYETA_PACKEDBYTES];
    uint8_t r_ref[POLYETA_PACKEDBYTES];
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for (i = 0; i < N; i++) {
        a.coeffs[i]=rand()%(2*ETA)-ETA;
    }    

    polyeta_pack(r_ref,&a);

    polyeta_pack_custom(r,&a);

    flag = true;
    for(i = 0; i < POLYETA_PACKEDBYTES; i++) {
        if(r[i] != r_ref[i]) {
            log_e("r[%d]=%d, while r_ref[%d]=%d", i, r[i], i, r_ref[i]);
            flag = false;
        }
    }

    if(flag) {
        printf("polyeta_pack_custom test pass!\n");
    }
    else {
        printf("polyeta_pack_custom test fail..\n");
    }

    return flag;
}

int main(){
    test_polyeta_pack_custom();
    return 0;
}