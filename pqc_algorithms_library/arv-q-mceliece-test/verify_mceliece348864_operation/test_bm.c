#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "gf.h"
#include "bm.h"

bool test_bm_custom(){
    gf s[ SYS_T*2 ];
    srand((unsigned)time(NULL));
    for(int i=0;i<SYS_T*2;i++) s[i]=rand()&GFMASK;

    gf locator_ref[ SYS_T+1 ];
    gf locator_res[ SYS_T+1 ];

    uint64_t start, end;

    start=read_cycle();
    bm(locator_ref, s);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    bm_custom(locator_res, s);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(SYS_T+1);i++){
        if(locator_ref[i]!=locator_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_bm_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    test_bm_custom();

    return 0;
}