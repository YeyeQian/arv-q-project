#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "synd.h"

bool test_synd_custom(){
    srand((unsigned)time(NULL));
    unsigned char r[ SYS_N/8 ];
    for (int i = 0; i < SYND_BYTES; i++)       r[i] = rand()&255;
    for (int i = SYND_BYTES; i < SYS_N/8; i++) r[i] = 0;
    gf g[ SYS_T+1 ];
    for(int i=0;i<SYS_T;i++)g[i]=rand()&GFMASK; 
    g[ SYS_T ] = 1;
	gf L[ SYS_N ];
    for(int i=0;i<SYS_N;i++)L[i]=rand()&GFMASK;

    gf s_ref[ SYS_T*2 ];
    gf s_res[ SYS_T*2 ];

    uint64_t start, end;

    start=read_cycle();
    synd(s_ref, g, L, r);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    synd_custom(s_res, g, L, r);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<SYS_T*2;i++){
        if(s_ref[i]!=s_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_synd_custom return with flag %d\n",flag);

    return flag;

}

int main(){
    test_synd_custom();

    return 0;
}