#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "root.h"
#include "crypto_hash.h"

bool test_root_custom(){
    gf locator[ SYS_T+1 ];
    gf L[ SYS_N ];
    srand((unsigned)time(NULL));
    for(int i=0;i<SYS_T+1;i++){
        locator[i]=rand()&GFMASK;
    }
    for(int i=0;i<SYS_N;i++){
        L[i]=rand()&GFMASK;
    }

    gf images_ref[ SYS_N ];
    gf images_res[ SYS_N ];

    uint64_t start, end;

    start=read_cycle();
    root(images_ref, locator, L);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    root_custom(images_res, locator, L);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(SYS_N);i++){
        if(images_ref[i]!=images_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_root_custom return with flag %d\n",flag);

    return flag;
}

bool test_root_shuffling_custom(){
    gf locator[ SYS_T+1 ];
    gf L[ SYS_N ];
    srand((unsigned)time(NULL));
    for(int i=0;i<SYS_T+1;i++){
        locator[i]=rand()&GFMASK;
    }
    for(int i=0;i<SYS_N;i++){
        L[i]=rand()&GFMASK;
    }

    gf images_ref[ SYS_N ];
    gf images_res[ SYS_N ];

    uint64_t start, end;

    start=read_cycle();
    root_custom(images_ref, locator, L);
    end=read_cycle();
    printf("root_custom finished with %lu cycles\n",end-start);

    //prepare shuffling index
    uint8_t seed[32] = {
    218, 26, 147, 106, 147, 243, 45, 212, 161, 50, 181, 182, 182, 128, 38, 156, 177, 46, 48, 128, 207, 132, 192, 50, 121, 195, 126, 255, 221, 222, 188, 83
    };
    uint8_t buffer[SYS_N<<1];
    shake256_custom(buffer,SYS_N<<1,seed,32);
    uint16_t* buffer_tmp=(uint16_t*)buffer;
    uint16_t shuffle_idx[SYS_N];
    for(int i=0;i<SYS_N;i++) {
        shuffle_idx[i]=i;
    }
    for(int i=SYS_N-1; i>=0; i--) {             //shuffling in pairs
        int selection=buffer_tmp[i]%(i+1);

        //shuffling zeta index
        uint16_t temp=shuffle_idx[i];
        shuffle_idx[i]=shuffle_idx[selection];
        shuffle_idx[selection]=temp;
    }
    for(int i=0;i<SYS_N;i++) {
        shuffle_idx[i]<<=1;
    }

    start=read_cycle();
    root_shuffling_custom(images_res, locator, L, shuffle_idx);
    end=read_cycle();
    printf("root_shuffling_custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(SYS_N);i++){
        if(images_ref[i]!=images_res[i]){
            flag=false;
            break;
        }
    }

    printf("test_root_custom return with flag %d\n",flag);

    return flag;
}

int main(){
    test_root_shuffling_custom();

    return 0;
}