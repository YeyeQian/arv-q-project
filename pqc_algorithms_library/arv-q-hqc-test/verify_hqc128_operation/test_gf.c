#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "parameters.h"
#include "gf.h"

#define LEN 1000
bool test_gfmulvx_u8(){
    uint8_t srcv[LEN];
    uint8_t scalar;
    uint8_t dstv_ref[LEN];
    uint8_t dstv_res[LEN];
    uint64_t start, end;

    srand((unsigned)time(NULL));
    for(int i=0;i<LEN;i++){
        srcv[i]=rand()&255;
        //srcv[i]=i;
    }
    scalar=rand()&255;
    //scalar=114;

    start=read_cycle();
    for(int i=0;i<LEN;i++){
        dstv_ref[i]=gf_mul(srcv[i],scalar);
    }
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    gfmul_vx_custom_u8(dstv_res,srcv,scalar,LEN);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<LEN;i++){
        if(dstv_ref[i]!=dstv_res[i]){
            flag=false;
            printf("dstv_ref[%d]=%d,dstv_res[%d]=%d\n",i,dstv_ref[i],i,dstv_res[i]);
            //break;
        }
    }

    printf("test_gfmulvx_u8 return with flag %d\n",flag);

    return flag;
}

bool test_gfmulvx_u8_scalarversion(){
    uint8_t srcv[LEN];
    uint8_t scalar;
    uint8_t dstv_ref[LEN];
    uint8_t dstv_res[LEN];
    uint64_t start, end;

    srand((unsigned)time(NULL));
    for(int i=0;i<LEN;i++){
        srcv[i]=rand()&255;
        //srcv[i]=i;
    }
    scalar=rand()&255;
    //scalar=114;

    start=read_cycle();
    for(int i=0;i<LEN;i++){
        dstv_ref[i]=gf_mul(srcv[i],scalar);
    }
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    csr_primpoly_rw(PARAM_GF_POLY);
    for(int i=0;i<LEN;i++){
        dstv_res[i]=gfmul_u8(srcv[i],scalar);
    }
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<LEN;i++){
        if(dstv_ref[i]!=dstv_res[i]){
            flag=false;
            printf("dstv_ref[%d]=%d,dstv_res[%d]=%d\n",i,dstv_ref[i],i,dstv_res[i]);
            //break;
        }
    }

    printf("test_gfmulvx_u8_scalarversion return with flag %d\n",flag);

    return flag;
}

bool test_gfmulvx_u16(){
    uint16_t srcv[LEN];
    uint16_t scalar;
    uint16_t dstv_ref[LEN];
    uint16_t dstv_res[LEN];
    uint64_t start, end;

    srand((unsigned)time(NULL));
    for(int i=0;i<LEN;i++){
        srcv[i]=rand()&255;
    }
    scalar=rand()&255;

    start=read_cycle();
    for(int i=0;i<LEN;i++){
        dstv_ref[i]=gf_mul(srcv[i],scalar);
    }
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    gfmul_vx_custom_u16(dstv_res,srcv,scalar,LEN);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<LEN;i++){
        if(dstv_ref[i]!=dstv_res[i]){
            flag=false;
            printf("dstv_ref[%d]=%d,dstv_res[%d]=%d\n",i,dstv_ref[i],i,dstv_res[i]);
            //break;
        }
    }

    printf("test_gfmulvx_u16 return with flag %d\n",flag);

    return flag;
}

bool test_gfinv_u8(){
    uint8_t srcv[LEN];
    uint8_t dstv_ref[LEN];
    uint8_t dstv_res[LEN];
    uint64_t start, end;

    srand((unsigned)time(NULL));
    for(int i=0;i<LEN;i++){
        srcv[i]=rand()&255;
    }

    start=read_cycle();
    for(int i=0;i<LEN;i++){
        dstv_ref[i]=gf_inverse(srcv[i]);
    }
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    csr_primpoly_rw(PARAM_GF_POLY);
    for(int i=0;i<LEN;i++){
        dstv_res[i]=gfinv_u8(srcv[i]);
    }
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<LEN;i++){
        if(dstv_ref[i]!=dstv_res[i]){
            flag=false;
            printf("dstv_ref[%d]=%d,dstv_res[%d]=%d\n",i,dstv_ref[i],i,dstv_res[i]);
            //break;
        }
    }

    printf("test_gfinv_u8 return with flag %d\n",flag);

    return flag;
}

bool test_gfinv_u16(){
    uint16_t srcv[LEN];
    uint16_t dstv_ref[LEN];
    uint16_t dstv_res[LEN];
    uint64_t start, end;

    srand((unsigned)time(NULL));
    for(int i=0;i<LEN;i++){
        srcv[i]=rand()&255;
    }

    start=read_cycle();
    for(int i=0;i<LEN;i++){
        dstv_ref[i]=gf_inverse(srcv[i]);
    }
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    csr_primpoly_rw(PARAM_GF_POLY);
    for(int i=0;i<LEN;i++){
        // dstv_res[i]=gfinv_u16(srcv[i]);
        dstv_res[i]=gf_inverse_custom(srcv[i]);
    }
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<LEN;i++){
        if(dstv_ref[i]!=dstv_res[i]){
            flag=false;
            printf("dstv_ref[%d]=%d,dstv_res[%d]=%d\n",i,dstv_ref[i],i,dstv_res[i]);
            //break;
        }
    }

    printf("test_gfinv_u16 return with flag %d\n",flag);

    return flag;
}

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
    //test_gfmulvx_u8();
    //test_gfmulvx_u8_scalarversion();
    //test_gfmulvx_u16();
    //test_gfinv_u8();
    //test_gfinv_u16();

    return 0;
}