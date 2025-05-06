#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "gf.h"

bool test_GF_mul_custom(){
    gf in0[SYS_T];
    gf in1[SYS_T];

    srand((unsigned)time(NULL));
    for(int i=0;i<SYS_T;i++){
        in0[i]=rand()&GFMASK;
        in1[i]=rand()&GFMASK;
    }

    gf out_ref[SYS_T];
    gf out_res[SYS_T];

    uint64_t start, end;

    start=read_cycle();
    GF_mul(out_ref, in0, in1);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    GF_mul_custom(out_res, in0, in1);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(SYS_T);i++){
        if(out_ref[i]!=out_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_GF_mul_custom return with flag %d\n",flag);

    return flag;
}

#define LEN 1000
bool test_gfmulvx_u16(){
    uint16_t srcv[LEN];
    uint16_t scalar;
    uint16_t dstv_ref[LEN];
    uint16_t dstv_res[LEN];
    uint64_t start, end;

    srand((unsigned)time(NULL));
    for(int i=0;i<LEN;i++){
        srcv[i]=rand()&GFMASK;
    }
    scalar=rand()&GFMASK;

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
            break;
        }
    }

    printf("test_gfmulvx_u16 return with flag %d\n",flag);

    return flag;
}

bool test_gfinv_u16(){
    uint16_t srcv[LEN];
    uint16_t dstv_ref[LEN];
    uint16_t dstv_res[LEN];
    uint64_t start, end;

    srand((unsigned)time(NULL));
    for(int i=0;i<LEN;i++){
        srcv[i]=rand()&GFMASK;
    }

    start=read_cycle();
    for(int i=0;i<LEN;i++){
        dstv_ref[i]=gf_inv(srcv[i]);
    }
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    csr_primpoly_rw(MC_GF_POLY);
    for(int i=0;i<LEN;i++){
        // dstv_res[i]=gfinv_u16(srcv[i]);
        dstv_res[i]=gf_inv_custom(srcv[i]);
    }
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<LEN;i++){
        if(dstv_ref[i]!=dstv_res[i]){
            flag=false;
            printf("dstv_ref[%d]=%d,dstv_res[%d]=%d\n",i,dstv_ref[i],i,dstv_res[i]);
            break;
        }
    }

    printf("test_gfinv_u16 return with flag %d\n",flag);

    return flag;
}

int main(){
    //test_GF_mul_custom();
    //test_gfmulvx_u16();
    test_gfinv_u16();

    return 0;
}