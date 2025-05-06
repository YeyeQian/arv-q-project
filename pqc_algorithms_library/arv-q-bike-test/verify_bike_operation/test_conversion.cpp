#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include "conversions.h"
#include "defs.h"

void test_convertByteToBinary_custom(){
    uint8_t in[R_SIZE]={0};
    uint8_t out_ref[R_BITS]={0};
    uint8_t out_res[R_BITS]={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<R_SIZE;i++){
        in[i]=rand()&255;
    }

    convertByteToBinary(out_ref,in,R_BITS);
    convertByteToBinary_custom(out_res,in,R_BITS);

    bool flag=true;
    for(int i=0;i<R_BITS;i++){
        if(out_ref[i]!=out_res[i]){
            printf("ref[%d]=%d,res[%d]=%d\n",i,out_ref[i],i,out_res[i]);
            flag=false;
            break;
        }
    }
    if(flag)printf("success in convertByteToBinary\n");
    else printf("fail in convertByteToBinary\n");
}

void test_convertBinaryToByte_custom(){
    uint8_t in[R_BITS]={0};
    uint8_t out_ref[R_SIZE]={0};
    uint8_t out_res[R_SIZE]={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<R_BITS;i++){
        in[i]=rand()&1;
    }

    convertBinaryToByte(out_ref,in,R_BITS);
    convertBinaryToByte_custom(out_res,in,R_BITS);

    bool flag=true;
    for(int i=0;i<R_SIZE;i++){
        if(out_ref[i]!=out_res[i]){
            printf("ref[%d]=%d,res[%d]=%d\n",i,out_ref[i],i,out_res[i]);
            flag=false;
            break;
        }
    }
    if(flag)printf("success in convertBinaryToByte\n");
    else printf("fail in convertBinaryToByte\n");
}

int main(){
    // test_convertByteToBinary_custom();
    test_convertBinaryToByte_custom();
}