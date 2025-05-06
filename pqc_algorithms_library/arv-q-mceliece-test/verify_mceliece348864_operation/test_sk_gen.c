#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "sk_gen.h"

bool test_genpoly_gen_custom(){
    gf f[ SYS_T ];
    srand((unsigned)time(NULL));
    for(int i=0;i<SYS_T;i++){
        f[i]=rand()&GFMASK;
    }
	gf irr_ref[ SYS_T ];
    gf irr_res[ SYS_T ];

    uint64_t start, end;

    start=read_cycle();
    genpoly_gen(irr_ref, f);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    genpoly_gen_custom(irr_res, f);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(SYS_T);i++){
        if(irr_ref[i]!=irr_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_genpoly_gen_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    test_genpoly_gen_custom();

    return 0;
}