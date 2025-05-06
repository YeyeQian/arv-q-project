#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "parameters.h"
#include "gf.h"

bool test_vect_mul_cycshift_custom(){
    uint64_t h[VEC_N_SIZE_64] = {0};
    uint32_t support_y[PARAM_OMEGA_R] = { 0 };
    uint64_t s_res[VEC_N_SIZE_64] = {0};
    uint64_t s_ref[VEC_N_SIZE_64] = {0};
    uint64_t start, end;

    srand((unsigned)time(NULL));
    for(int i=0;i<VEC_N_SIZE_64;i++){
        h[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());
    }
    h[VEC_N_SIZE_64 - 1] &= BITMASK(PARAM_N, 64);

    for(int i=0;i<PARAM_OMEGA;i++){
        support_y[i]=rand()%PARAM_N;
    }
    
    start=read_cycle();
    vect_mul_cycshift(s_ref,h,support_y,PARAM_OMEGA);
    //vect_mul_cycshift_in32(s_ref,h,support_y,PARAM_OMEGA);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    vect_mul_cycshift_custom(s_res,h,support_y,PARAM_OMEGA);
    end=read_cycle();
    printf("vect_mul_cycshift_custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<VEC_N_SIZE_64;i++){
        if(s_ref[i]!=s_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_vect_mul_cycshift_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    test_vect_mul_cycshift_custom();

    return 0;
}