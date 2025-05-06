#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdbool.h>
#include "gf2x_port.h"

void test_gf2x_sqr_custom(){
    pad_r_t a={0};
    srand((unsigned)time(NULL));
    for(int i=0;i<R_SIZE;i++){
        a.val.raw[i]=rand()&255;
    }
    dbl_pad_r_t c_ref={0};
    dbl_pad_r_t c_res={0};
    gf2x_sqr_port(&c_ref,&a);
    gf2x_sqr_custom(&c_res,&a);
    
    bool flag=true;
    for(int i=0;i<sizeof(c_ref);i++){
        if(c_ref.raw[i]!=c_res.raw[i]){
            flag=false;
            break;
        }
    }

    printf("test_gf2x_sqr_custom return with flag=%d\n",flag);
    printf("test_gf2x_sqr_custom return with flag=%d\n",flag);
    return;
}

void test_gf2x_red_custom(){
    pad_r_t a={0};
    srand((unsigned)time(NULL));
    for(int i=0;i<R_SIZE;i++){
        a.val.raw[i]=rand()&255;
    }
    dbl_pad_r_t c={0};
    gf2x_sqr_custom(&c,&a);
    pad_r_t ref={0};
    pad_r_t res={0};
    gf2x_red_port(&ref,&c);
    gf2x_red_custom(&res,&c);

    bool flag=true;
    for(int i=0;i<R_SIZE;i++){
        if(ref.val.raw[i]!=res.val.raw[i]){
            flag=false;
            printf("ref.val.raw[%d]=%d,res.val,raw[%d]=%d\n",i,ref.val.raw[i],i,res.val.raw[i]);
            // break;
        }
    }

    for(int i=0;i<(R_PADDED_BYTES-R_SIZE);i++){
        if(ref.pad[i]!=res.pad[i]){
            flag=false;
            printf("ref.pad[%d]=%d,res.pad[%d]=%d\n",i,ref.pad[i],i,res.pad[i]);
            // break;
        }
    }

    printf("test_gf2x_red_custom return with flag=%d\n",flag);
    printf("test_gf2x_red_custom return with flag=%d\n",flag);
    return;
}

void test_gf2x_mod_mul_with_ctx_custom(){
    pad_r_t a={0};
    pad_r_t b={0};
    pad_r_t c_ref={0};
    pad_r_t c_res={0};
    srand((unsigned)time(NULL));
    for(int i=0;i<R_SIZE;i++){
        a.val.raw[i]=rand()&255;
        b.val.raw[i]=rand()&255;
    }
    uint64_t start,end;
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    gf2x_mod_mul_with_ctx(&c_ref,&a,&b);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("gf2x_mod_mul_with_ctx cost %lu cycles\n",end-start);

    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );
    gf2x_mod_mul_with_ctx_custom(&c_res,&a,&b);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );
    printf("gf2x_mod_mul_with_ctx_custom cost %lu cycles\n",end-start);

    bool flag=true;
    for(int i=0;i<R_SIZE;i++){
        if(c_ref.val.raw[i]!=c_res.val.raw[i]){
            flag=false;
            printf("ref.val.raw[%d]=%d,res.val,raw[%d]=%d\n",i,c_ref.val.raw[i],i,c_res.val.raw[i]);
            // break;
        }
    }

    for(int i=0;i<(R_PADDED_BYTES-R_SIZE);i++){
        if(c_ref.pad[i]!=c_res.pad[i]){
            flag=false;
            printf("ref.pad[%d]=%d,res.pad[%d]=%d\n",i,c_ref.pad[i],i,c_res.pad[i]);
            // break;
        }
    }

    printf("test_gf2x_mod_mul_with_ctx_custom return with flag=%d\n",flag);
    printf("test_gf2x_mod_mul_with_ctx_custom return with flag=%d\n",flag);
    return;
}

void test_mul2_512_custom(){
    uint64_t a[8]={0};
    uint64_t b[8]={0};
    srand((unsigned)time(NULL));
    for(int i=0;i<8;i++){
        a[i]=(uint64_t)rand()<<48 | (uint64_t)rand()<<32 | (uint64_t)rand()<<16 | (uint64_t)rand();
        b[i]=(uint64_t)rand()<<48 | (uint64_t)rand()<<32 | (uint64_t)rand()<<16 | (uint64_t)rand();
    }

    uint64_t ref[16]={0};
    uint64_t res[16]={0};

    vuint64m1_t va,vb;
    vuint64m2_t vab;
    va=vle64_v_u64m1(a,8);
    vb=vle64_v_u64m1(b,8);
    vab=mul2_512_custom(va,vb);
    vse64_v_u64m2(res,vab,16);

    mul2_512(&ref[8],&ref[0],a,b);

    bool flag=true;
    for(int i=0;i<16;i++){
        if(ref[i]!=res[i]){
            flag=false;
            printf("ref[%d]=%llx,res[%d]=%llx\n",i,ref[i],i,res[i]);
            // break;
        }
    }

    printf("test_mul2_512_custom return with flag=%d\n",flag);
    printf("test_mul2_512_custom return with flag=%d\n",flag);
    return;
}

void test_gf2x_mul8_512_int_custom(){
    uint64_t a[8]={0};
    uint64_t b[8]={0};
    srand((unsigned)time(NULL));
    for(int i=0;i<8;i++){
        a[i]=(uint64_t)rand()<<48 | (uint64_t)rand()<<32 | (uint64_t)rand()<<16 | (uint64_t)rand();
        b[i]=(uint64_t)rand()<<48 | (uint64_t)rand()<<32 | (uint64_t)rand()<<16 | (uint64_t)rand();
    }

    uint64_t ref[16]={0};
    uint64_t res[16]={0};

    vuint64m1_t va,vb;
    vuint64m2_t vab;
    va=vle64_v_u64m1(a,8);
    vb=vle64_v_u64m1(b,8);
    vab=gf2x_mul8_512_int_custom(va,vb);
    vse64_v_u64m2(res,vab,16);

    gf2x_mul8_512_int(&ref[8],&ref[0],a,b);

    bool flag=true;
    for(int i=0;i<16;i++){
        if(ref[i]!=res[i]){
            flag=false;
            printf("ref[%d]=%llx,res[%d]=%llx\n",i,ref[i],i,res[i]);
            // break;
        }
    }

    printf("test_gf2x_mul8_512_int_custom return with flag=%d\n",flag);
    printf("test_gf2x_mul8_512_int_custom return with flag=%d\n",flag);
    return;
}

void test_gf2x_mul_base_vpclmul_custom(){
    uint64_t a[16]={0};
    uint64_t b[16]={0};
    srand((unsigned)time(NULL));
    for(int i=0;i<16;i++){
        a[i]=(uint64_t)rand()<<48 | (uint64_t)rand()<<32 | (uint64_t)rand()<<16 | (uint64_t)rand();
        b[i]=(uint64_t)rand()<<48 | (uint64_t)rand()<<32 | (uint64_t)rand()<<16 | (uint64_t)rand();
    }

    uint64_t ref[32]={0};
    uint64_t res[32]={0};

    gf2x_mul_base_vpclmul_custom(res,a,b);

    gf2x_mul_base_vpclmul(ref,a,b);

    bool flag=true;
    for(int i=0;i<32;i++){
        if(ref[i]!=res[i]){
            flag=false;
            printf("ref[%d]=%llx,res[%d]=%llx\n",i,ref[i],i,res[i]);
            // break;
        }
    }

    printf("test_gf2x_mul_base_vpclmul_custom return with flag=%d\n",flag);
    printf("test_gf2x_mul_base_vpclmul_custom return with flag=%d\n",flag);
    return;
}

void test_k_sqr_custom(){
    pad_r_t a={0};
    pad_r_t c_ref={0};
    pad_r_t c_res={0};

    srand((unsigned)time(NULL));
    for(int i=0;i<R_SIZE;i++){
        a.val.raw[i]=rand()&255;
    }

    uint16_t l_param=6162;

    k_sqr_port(&c_ref,&a,l_param);

    k_sqr_custom(&c_res,&a,l_param);

    bool flag=true;
    for(int i=0;i<sizeof(c_ref.val);i++){
        if(c_ref.val.raw[i]!=c_res.val.raw[i]){
            flag=false;
            break;
        }
    }
    for(int i=0;i<sizeof(c_ref.pad);i++){
        if(c_ref.pad[i]!=c_res.pad[i]){
            flag=false;
            break;
        }
    }

    printf("test_k_sqr_custom return with flag=%d\n",flag);
    printf("test_k_sqr_custom return with flag=%d\n",flag);
    return;
}

int main(){
    // test_gf2x_sqr_custom();
    // test_gf2x_red_custom();
    // test_gf2x_mod_mul_with_ctx_custom();
    // test_mul2_512_custom();
    // test_gf2x_mul8_512_int_custom();
    test_gf2x_mul_base_vpclmul_custom();
    // test_k_sqr_custom();
}