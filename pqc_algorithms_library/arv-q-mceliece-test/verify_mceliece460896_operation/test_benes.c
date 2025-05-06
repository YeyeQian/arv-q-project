#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "benes.h"

bool test_support_gen_custom(){
    uint8_t s[COND_BYTES];
    srand((unsigned)time(NULL));
    for(int i=0;i<COND_BYTES;i++) s[i]=rand()&255;

    gf L_ref[ SYS_N ];
    gf L_res[ SYS_N ];

    uint64_t start, end;

    start=read_cycle();
    support_gen(L_ref, s);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    support_gen_custom(L_res, s);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<SYS_N;i++){
        if(L_ref[i]!=L_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_support_gen_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    test_support_gen_custom();

    return 0;
}