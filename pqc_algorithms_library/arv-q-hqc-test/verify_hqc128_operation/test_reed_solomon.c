#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "parameters.h"
#include "reed_solomon.h"

bool test_reed_solomon_encode_custom(){
    uint64_t m[VEC_K_SIZE_64]={0};
    srand((unsigned)time(NULL));
    for(int i=0;i<VEC_K_SIZE_64;i++){
        m[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());
    }
    uint64_t ref[VEC_N1_SIZE_64] = {0};
    uint64_t res[VEC_N1_SIZE_64] = {0};

    uint64_t start, end;

    start=read_cycle();
    reed_solomon_encode(ref,m);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    reed_solomon_encode_custom(res,m);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<VEC_N1_SIZE_64;i++){
        if(ref[i]!=res[i]){
            flag=false;
            break;
        }
    }

    printf("test_reed_solomon_encode_custom return with flag %d\n",flag);

    return flag;
}

bool test_reed_solomon_decode_custom(){
    uint64_t m_ref[VEC_K_SIZE_64] = {
        4164766533729176191ull,8371366596015814577ull
    };
    uint64_t m_res[VEC_K_SIZE_64] = {0};

    uint64_t tmp[VEC_N1_SIZE_64]={
        5380635517672513561ull,9442858707034731139ull,6699578827259551292ull,14591435790409135210ull,17848110347540019036ull,127736917053463ull
    };

    uint64_t start, end;

    start=read_cycle();
    reed_solomon_decode(m_ref,tmp);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    reed_solomon_decode_custom(m_res,tmp);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<VEC_K_SIZE_64;i++){
        if(m_ref[i]!=m_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_reed_solomon_decode_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    // test_reed_solomon_encode_custom();
    test_reed_solomon_decode_custom();

    return 0;
}