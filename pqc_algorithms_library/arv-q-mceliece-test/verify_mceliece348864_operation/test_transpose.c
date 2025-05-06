#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "transpose.h"

bool test_transpose_64x64_custom(){
    uint64_t bs[64];
    srand((unsigned)time(NULL));
    for(int i=0;i<64;i++){
        bs[i]=((uint64_t)rand()<<48)|((uint64_t)rand()<<32)|((uint64_t)rand()<<16)|((uint64_t)rand());;
    }
    uint64_t res[64], ref[64];

    uint64_t start, end;

    start=read_cycle();
    transpose_64x64(ref, bs);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    transpose_64x64_custom(res, bs);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<64;i++){
        if(ref[i]!=res[i]){
            flag=false;
            break;
        }
    }

    printf("test_transpose_64x64_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    test_transpose_64x64_custom();

    return 0;
}