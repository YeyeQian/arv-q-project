#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "parameters.h"
#include "vector.h"

bool test_vect_set_random_fixed_weight_xy_custom(){
    uint8_t sk_seed[SEED_BYTES];
    srand((unsigned)time(NULL));
    for(int i=0;i<SEED_BYTES;i++){
        sk_seed[i]=rand()&255;
    }

    uint64_t x_ref[VEC_N_SIZE_64] = {0};
    uint64_t y_ref[VEC_N_SIZE_64] = {0};
    uint32_t support_y_ref[PARAM_OMEGA_R] = { 0 };

    uint64_t x_res[VEC_N_SIZE_64] = {0};
    uint64_t y_res[VEC_N_SIZE_64] = {0};
    uint32_t support_y_res[PARAM_OMEGA_R] = { 0 };

    uint64_t start, end;

    start=read_cycle();
    vect_set_random_fixed_weight_xy(x_ref,y_ref,support_y_ref,sk_seed);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    vect_set_random_fixed_weight_xy_custom(x_res,y_res,support_y_res,sk_seed);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag_x=true;
    for(int i=0;i<VEC_N_SIZE_64;i++){
        if(x_ref[i]!=x_res[i]){
            flag_x=false;
            break;
        }
    }

    bool flag_y=true;
    for(int i=0;i<VEC_N_SIZE_64;i++){
        if(y_ref[i]!=y_res[i]){
            flag_y=false;
            break;
        }
    }

    bool flag_support=true;
    for(int i=0;i<PARAM_OMEGA_R;i++){
        if(support_y_ref[i]!=support_y_res[i]){
            flag_support=false;
            break;
        }
    }

    printf("test_vect_set_random_fixed_weight_xy_custom return with flag %d\n",flag_x&&flag_y&&flag_support);

    return flag_x&&flag_y&&flag_support;
}

bool test_vect_set_random_fixed_weight_r1r2e_custom(){
    uint8_t theta[SHAKE256_512_BYTES] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<SHAKE256_512_BYTES;i++){
        theta[i]=rand()&255;
    }

    uint64_t r1_ref[VEC_N_SIZE_64] = {0};
    uint64_t r2_ref[VEC_N_SIZE_64] = {0};
    uint64_t e_ref[VEC_N_SIZE_64] = {0};
    uint32_t support_r2_ref[PARAM_OMEGA_R] = { 0 };

    uint64_t r1_res[VEC_N_SIZE_64] = {0};
    uint64_t r2_res[VEC_N_SIZE_64] = {0};
    uint64_t e_res[VEC_N_SIZE_64] = {0};
    uint32_t support_r2_res[PARAM_OMEGA_R] = { 0 };

    uint64_t start, end;

    start=read_cycle();
    vect_set_random_fixed_weight_r1r2e(r1_ref,r2_ref,support_r2_ref,e_ref,theta);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    vect_set_random_fixed_weight_r1r2e_custom(r1_res,r2_res,support_r2_res,e_res,theta);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag_r1=true;
    for(int i=0;i<VEC_N_SIZE_64;i++){
        if(r1_ref[i]!=r1_res[i]){
            flag_r1=false;
            break;
        }
    }

    bool flag_r2=true;
    for(int i=0;i<VEC_N_SIZE_64;i++){
        if(r2_ref[i]!=r2_res[i]){
            flag_r2=false;
            break;
        }
    }

    bool flag_support=true;
    for(int i=0;i<PARAM_OMEGA_R;i++){
        if(support_r2_ref[i]!=support_r2_res[i]){
            flag_support=false;
            break;
        }
    }

    bool flag_e=true;
    for(int i=0;i<VEC_N_SIZE_64;i++){
        if(e_ref[i]!=e_res[i]){
            flag_e=false;
            break;
        }
    }

    printf("test_vect_set_random_fixed_weight_r1r2e_custom return with flag %d\n",flag_r1&&flag_r2&&flag_support&&flag_e);

    return flag_r1&&flag_r2&&flag_support&&flag_e;
}

bool test_vect_set_random_h_custom(){
    uint8_t pk_seed[SEED_BYTES];
    srand((unsigned)time(NULL));
    for(int i=0;i<SEED_BYTES;i++){
        pk_seed[i]=rand()&255;
    }

    uint64_t h_ref[VEC_N_SIZE_64] = {0};
    uint64_t h_res[VEC_N_SIZE_64] = {0};

    uint64_t start, end;

    start=read_cycle();
    vect_set_random_h(pk_seed,h_ref);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    vect_set_random_h_custom(pk_seed,h_res);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<VEC_N_SIZE_64;i++){
        if(h_ref[i]!=h_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_vect_set_random_h_custom return with flag %d\n",flag);

    return flag;
}

bool test_vect_add_custom(){
    uint64_t x[VEC_N_SIZE_64] = {0};
    uint64_t y[VEC_N_SIZE_64] = {0};
    uint64_t res[VEC_N_SIZE_64] = {0};
    uint64_t ref[VEC_N_SIZE_64] = {0};

    srand((unsigned)time(NULL));
    for(int i=0;i<VEC_N_SIZE_64;i++){
        x[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());
        y[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());
    }
    x[VEC_N_SIZE_64 - 1] &= BITMASK(PARAM_N, 64);
    y[VEC_N_SIZE_64 - 1] &= BITMASK(PARAM_N, 64);

    uint64_t start, end;

    start=read_cycle();
    vect_add(ref, x, y, VEC_N_SIZE_64);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    vect_add(res, x, y, VEC_N_SIZE_64);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<VEC_N_SIZE_64;i++){
        if(ref[i]!=res[i]){
            flag=false;
            break;
        }
    }

    printf("test_vect_add_custom return with flag %d\n",flag);

    return flag;
}

bool test_vect_compare_custom(){
    uint8_t result;
    uint64_t x[VEC_N_SIZE_64] = {0};
    uint64_t y[VEC_N_SIZE_64] = {0};
    srand((unsigned)time(NULL));
    for(int i=0;i<VEC_N_SIZE_64;i++){
        x[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());
        y[i]=x[i];
    }
    result=vect_compare_custom((uint8_t *)x, (uint8_t *)y, VEC_N_SIZE_BYTES);

    bool flag=result==0;

    printf("test_vect_compare_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    test_vect_set_random_fixed_weight_xy_custom();
    // test_vect_set_random_fixed_weight_r1r2e_custom();
    // test_vect_set_random_h_custom();
    // test_vect_add_custom();
    // test_vect_compare_custom();
    return 0;
}