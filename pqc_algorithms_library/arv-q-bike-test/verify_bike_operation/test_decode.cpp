#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include "decode.h"
#include "defs.h"

void test_getHammingWeight_custom(){
    uint8_t in[R_BITS]={0};
    uint32_t count_ref,count_res;

    srand((unsigned)time(NULL));
    for(int i=0;i<R_BITS;i++){
        in[i]=rand()&1;
    }

    count_ref=getHammingWeight(in,R_BITS);
    count_res=getHammingWeight_custom(in,R_BITS);

    if(count_ref==count_res)printf("success in test_getHammingWeight_custom");
    else printf("fail in test_getHammingWeight_custom");
    return;
}

void test_recompute_syndrome_custom(){
    uint32_t h0_compact[DV]={0};
    uint32_t h1_compact[DV]={0};
    uint32_t pos=0;
    uint8_t s_ref[R_BITS]={0};
    uint8_t s_res[R_BITS]={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<DV;i++){
        h0_compact[i]=((((uint32_t)rand())<<16)|((uint32_t)rand()))%R_BITS;
        h1_compact[i]=((((uint32_t)rand())<<16)|((uint32_t)rand()))%R_BITS;
    }
    pos=((((uint32_t)rand())<<16)|((uint32_t)rand()))%(R_BITS<<1);

    recompute_syndrome(s_ref,pos,h0_compact,h1_compact);
    recompute_syndrome_custom(s_res,pos,h0_compact,h1_compact);

    bool flag=true;
    for(int i=0;i<R_BITS;i++){
        if(s_ref[i]!=s_res[i]){
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_recompute_syndrome_custom");
    else printf("fail in test_recompute_syndrome_custom");
    return;
}

void test_ctr_custom(){
    uint32_t h_compact_col[DV]={0};
    int position;
    uint8_t s[R_BITS]={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<DV;i++){
        h_compact_col[i]=((((uint32_t)rand())<<16)|((uint32_t)rand()))%R_BITS;
    }
    position=((((uint32_t)rand())<<16)|((uint32_t)rand()))%R_BITS;
    for(int i=0;i<R_BITS;i++){
        s[i]=rand()&1;
    }
    uint32_t ctr_ref,ctr_res;
    ctr_ref=ctr(h_compact_col,position,s);
    ctr_res=ctr_custom(h_compact_col,position,s);

    if(ctr_ref==ctr_res)printf("success in test_ctr_custom");
    else printf("fail in test_ctr_custom");
    return;
}

void test_BGF_decoder_custom(){
    uint8_t e[R_BITS*2]={0};
    uint8_t s[R_BITS]={0};
    uint32_t h0_compact[DV]={0};
    uint32_t h1_compact[DV]={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<DV;i++){
        h0_compact[i]=((((uint32_t)rand())<<16)|((uint32_t)rand()))%R_BITS;
        h1_compact[i]=((((uint32_t)rand())<<16)|((uint32_t)rand()))%R_BITS;
    }
    for(int i=0;i<R_BITS;i++){
        s[i]=rand()&1;
    }
    for(int i=0;i<2*R_BITS;i++){
        e[i]=rand()&1;
    }
    printf("Data init finished\n");

    int ref,res;
    ref=BGF_decoder(e,s,h0_compact,h1_compact);
    printf("BGF_decoder finished\n");
    res=BGF_decoder_custom(e,s,h0_compact,h1_compact);

    if(ref==res)printf("success in test_BGF_decoder_custom");
    else printf("fail in test_BGF_decoder_custom");
    return;
}

void test_transpose_custom(){
    uint8_t row[R_BITS]={0};
    
    uint8_t col_ref[R_BITS]={0};
    uint8_t col_res[R_BITS]={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<R_BITS;i++){
        row[i]=rand()&1;
    }

    transpose(col_ref,row);
    transpose_custom(col_res,row);

    bool flag=true;
    for(int i=0;i<R_BITS;i++){
        if(col_ref[i]!=col_res[i]){
            flag=false;
            break;
        }
    }
    if(flag)printf("success in test_transpose_custom");
    else printf("fail in test_transpose_custom");
    return;
}

int main(){
    // test_getHammingWeight_custom();
    // test_ctr_custom();
    // test_recompute_syndrome_custom();
    // test_BGF_decoder_custom();
    test_transpose_custom();
}