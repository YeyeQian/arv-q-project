#include <stddef.h>
#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "params.h"
#include "reduce.h"
#include "poly.h"
#include "ntt.h"
#include "fips202.h"

#if (VLEN == 256)
#elif (VLEN == 512)
//#define FORTEST
#ifdef FORTEST
void print_data(uint16_t* data, uint16_t len) {
	for (int i = 0; i < len; i++) {
		if(i!=(len-1))printf("%d,", data[i]);
		else printf("%d", data[i]);
		if ((i + 1) % 16 == 0) {
			printf("\n");
		}
	}
}
#endif
void generate_intt_shuffling_idx(uint16_t srcpoly_idx[128*7],uint16_t dstpoly_idx[256*7],uint8_t prng_buf[128*7]){
    //intialize scrpoly_idx
    for(int i=0;i<7;i++){
        for(int j=0;j<128;j++){
            srcpoly_idx[i*128+j]=j;
        }
    }
    //shuffling srcpoly_idx
    for(int i=0;i<7;i++){
        for(int j=127; j>=0; j--) {
            int selection=prng_buf[i*128+j]%(j+1);
            uint16_t temp=srcpoly_idx[i*128+j];
            srcpoly_idx[i*128+j]=srcpoly_idx[i*128+selection];
            srcpoly_idx[i*128+selection]=temp;
        }
    }
    //compute dstpoly_idx
    for(int i=0;i<6;i++){
        for(int j=0;j<128;j++){
            dstpoly_idx[i*256+2*j+0]=2*srcpoly_idx[i*128+j];
            dstpoly_idx[i*256+2*j+1]=2*srcpoly_idx[i*128+j]+1;
        }
    }
    //for last stage, consider bit-reverse permutaiton
    uint16_t br_tree[256];
    for(int i=0;i<128;i++){
        uint16_t tmp=byteoffset_even[i]>>1;
        br_tree[i]=tmp;
        br_tree[i+128]=tmp+1;
    }
    uint16_t tmp_dstpoly_idx_stage7[256];
    for(int j=0;j<128;j++){
        tmp_dstpoly_idx_stage7[2*j+0]=2*srcpoly_idx[6*128+j];
        tmp_dstpoly_idx_stage7[2*j+1]=2*srcpoly_idx[6*128+j]+1;
    }
    for(int i=0;i<256;i++){
        dstpoly_idx[6*256+i]=br_tree[tmp_dstpoly_idx_stage7[i]];
    }
    //for test: print all the shuffling results for check
#ifdef FORTEST
    printf("srcpoly_idx:\n");
    for(int i=0;i<7;i++){
        printf("stage%d:\n",i+1);
        print_data(&srcpoly_idx[i*128],128);
    }
    printf("dstpoly_idx:\n");
    for(int i=0;i<7;i++){
        printf("stage%d:\n",i+1);
        print_data(&dstpoly_idx[i*256],256);
    }
#endif
    //turn idx into byte offset for indexed-load-store usage
    for(int i=0;i<7;i++){
        for(int j=0;j<128;j++){
            srcpoly_idx[i*128+j]<<=1;
        }
        for(int j=0;j<256;j++){
            dstpoly_idx[i*256+j]<<=1;
        }
    }
}

void invntt_cg_custom_reorder_shuffling(int16_t r[KYBER_N],const uint8_t seed[32]) {
    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    uint8_t buffer[24];
    uint8_t idx[24];
    uint8_t idx_offset; 
    for(int i=0;i<5;i++){
        for(int j=0;j<4;j++){
            idx[i*4+j]=j;
        }
    }
    idx[20]=0,idx[21]=1,idx[22]=0,idx[23]=1;
    shake128_custom(buffer,24,seed,32);
    for(int i=0;i<5;i++){
        for(int j=3;j>=0;j--){
            int selection=buffer[i*4+j]%(j+1);
            uint8_t temp=idx[i*4+j];
            idx[i*4+j]=idx[i*4+selection];
            idx[i*4+selection]=temp;
        }
    }
    for(int i=5;i<7;i++){
        for(int j=1;j>=0;j--){
            int selection=buffer[i*4+j]%(j+1);
            uint8_t temp=idx[i*4+j];
            idx[i*4+j]=idx[i*4+selection];
            idx[i*4+selection]=temp;
        }
    }

    int16_t r_temp0[256];
    int16_t r_temp1[256];
    unsigned int i, j;
    int16_t* r_ptr0 = NULL;
    int16_t* r_ptr1 = NULL;
    int16_t* r_ptr2 = NULL;
    const int16_t* zeta_ptr = zetas_inv_in_order;

    vint16m1_t v_zetas;
    vint16m1_t v_in0, v_in1;
    vint16m2_t v_out0;
    vint16m2_t v_in2, v_in3;
    vint16m4_t v_out1;

    // stage 1
    idx_offset=0;
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    // same_num = 2
    csr_zeta_selMode_rw(GET_ZETA_MODE(1, 5));
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/2=16 zetas
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx<<4;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 16);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=64;

    // stage 2
    idx_offset=4;
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r;
    // same_num = 4
    csr_zeta_selMode_rw(GET_ZETA_MODE(2, 5));
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/4=8 zetas
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx<<3;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 8);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=32;

    // stage 3
    idx_offset=8;
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    // same_num = 8
    csr_zeta_selMode_rw(GET_ZETA_MODE(3, 5));  
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/8=4 zetas
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx<<2;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 4);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=16;

    // stage 4
    idx_offset=12;
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r;
    // same_num = 16
    csr_zeta_selMode_rw(GET_ZETA_MODE(4, 5));  
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/16=2 zetas
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx<<1;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 2);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=8;

    // stage 5
    idx_offset=16;
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    // same_num = 32
    csr_zeta_selMode_rw(GET_ZETA_MODE(5, 5));
    for(i = 0; i < 4; i++) {                // each 32* butterfly, all same zeta
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 1);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=4;

    // stage 6
    idx_offset=20;
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r_temp1;
    // same_num = 64
    csr_zeta_selMode_rw(GET_ZETA_MODE(6, 5));
    for(i = 0; i < 2; i++) {              // each 64* butterfly, all same zeta
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx;
        uint16_t src_offset=shuffled_idx<<6;
        uint16_t dst_offset=shuffled_idx<<7;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 1);
        v_in2 = vle16_v_i16m2(r_ptr0+src_offset, 64);
        v_in3 = vle16_v_i16m2(r_ptr1+src_offset, 64);
        v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
        vse16_v_i16m4(r_ptr2+dst_offset, v_out1, 128);
    }
    zeta_ptr+=2;

    // indexed store: used as offsets
    vuint16m4_t v_offset;
    vuint16m4_t v_offset_act;
    // tree pointer
    const uint16_t* tree_prt = NULL;

    // stage 7
    // same_num = 128
    csr_zeta_selMode_rw(GET_ZETA_MODE(7, 5));
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    idx_offset=22;
    r_ptr0 = r_temp1;
    r_ptr1 = &r_temp1[128];
    tree_prt = byteoffset_even;
    v_offset = vle16_v_u16m4(tree_prt, 128);
    for(j = 0; j < 2; j++) {                // each 64* butterfly, all same zeta
        uint8_t shuffled_idx=idx[idx_offset+j];
        uint16_t src_offset=shuffled_idx<<6;
        uint16_t dst_offset=shuffled_idx<<7;
        v_in2 = vle16_v_i16m2(r_ptr0+src_offset, 64);
        v_in3 = vle16_v_i16m2(r_ptr1+src_offset, 64);
        v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
        // reorder, from special bit-reversed order to standard order
        v_offset_act = vadd_vx_u16m4(v_offset, shuffled_idx<<1, 128);
        vsuxei16_v_i16m4(r, v_offset_act, v_out1, 128);
    }    
}

void invntt_cg_custom_reorder_shuffling_fine(int16_t r[KYBER_N],const uint8_t seed[32]) {
    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    uint8_t buffer[28];
    uint8_t idx[28];
    uint8_t idx_offset;
    for(int i=0;i<7;i++){
        for(int j=0;j<4;j++){
            idx[i*4+j]=j;
        }
    }
    shake128_custom(buffer,28,seed,32);
    for(int i=0;i<7;i++){
        for(int j=3;j>=0;j--){
            int selection=buffer[i*4+j]%(j+1);
            uint8_t temp=idx[i*4+j];
            idx[i*4+j]=idx[i*4+selection];
            idx[i*4+selection]=temp;
        }
    }

    int16_t r_temp0[256];
    int16_t r_temp1[256];
    unsigned int i, j;
    int16_t* r_ptr0 = NULL;
    int16_t* r_ptr1 = NULL;
    int16_t* r_ptr2 = NULL;
    const int16_t* zeta_ptr = zetas_inv_in_order;

    vint16m1_t v_zetas;
    vint16m1_t v_in0, v_in1;
    vint16m2_t v_out0;

    // stage 1
    idx_offset=0;
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    // same_num = 2
    csr_zeta_selMode_rw(GET_ZETA_MODE(1, 5));
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/2=16 zetas
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx<<4;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 16);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=64;

    // stage 2
    idx_offset=4;
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r;
    // same_num = 4
    csr_zeta_selMode_rw(GET_ZETA_MODE(2, 5));
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/4=8 zetas
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx<<3;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 8);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=32;

    // stage 3
    idx_offset=8;
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    // same_num = 8
    csr_zeta_selMode_rw(GET_ZETA_MODE(3, 5));  
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/8=4 zetas
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx<<2;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 4);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=16;

    // stage 4
    idx_offset=12;
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r;
    // same_num = 16
    csr_zeta_selMode_rw(GET_ZETA_MODE(4, 5));  
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/16=2 zetas
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx<<1;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 2);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=8;

    // stage 5
    idx_offset=16;
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    // same_num = 32
    csr_zeta_selMode_rw(GET_ZETA_MODE(5, 5));
    for(i = 0; i < 4; i++) {                // each 32* butterfly, all same zeta
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 1);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=4;

    // stage 6
    idx_offset=20;
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r_temp1;
    // same_num = 64
    csr_zeta_selMode_rw(GET_ZETA_MODE(6, 5));
    for(i = 0; i < 4; i++) {              // each 64* butterfly, all same zeta
        uint8_t shuffled_idx=idx[idx_offset+i];
        uint16_t zeta_offset=shuffled_idx>>1;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_zetas = vle16_v_i16m1(zeta_ptr+zeta_offset, 1);
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vse16_v_i16m2(r_ptr2+dst_offset, v_out0, 64);
    }
    zeta_ptr+=2;

    // indexed store: used as offsets
    vuint16m2_t v_offset;
    // tree pointer
    const uint16_t* tree_prt = NULL;

    // stage 7
    // same_num = 128
    csr_zeta_selMode_rw(GET_ZETA_MODE(7, 5));
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    idx_offset=24;
    r_ptr0 = r_temp1;
    r_ptr1 = &r_temp1[128];
    tree_prt = byteoffset_even;
    for(j = 0; j < 4; j++) {                // each 64* butterfly, all same zeta
        uint8_t shuffled_idx=idx[idx_offset+j];
        uint16_t tree_ptr_offset=(shuffled_idx&1)<<6;
        uint16_t src_offset=shuffled_idx<<5;
        uint16_t dst_offset=shuffled_idx<<6;
        v_in0 = vle16_v_i16m1(r_ptr0+src_offset, 32);
        v_in1 = vle16_v_i16m1(r_ptr1+src_offset, 32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        // reorder, from special bit-reversed order to standard order
        v_offset = vle16_v_u16m2(tree_prt+tree_ptr_offset, 64);
        v_offset = vadd_vx_u16m2(v_offset, (shuffled_idx>>1)<<1, 64);
        vsuxei16_v_i16m2(r, v_offset, v_out0, 64);
    }
}

void invntt_cg_custom_reorder_redundant(int16_t r[KYBER_N], const uint8_t seed[32]) {
    uint8_t buffer[128];
    shake128_custom(buffer,128,seed,32);

    int16_t r_temp0[256];
    int16_t r_temp1[256];
    unsigned int i, j;
    int16_t* r_ptr0 = NULL;
    int16_t* r_ptr1 = NULL;
    int16_t* r_ptr2 = NULL;
    const int16_t* zeta_ptr = zetas_inv_in_order;

    vint16m1_t v_zetas;
    vint16m1_t v_in0, v_in1;
    vint16m2_t v_out0;
    vint16m2_t v_in2, v_in3;
    vint16m4_t v_out1;

    // stage 1
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    // same_num = 2
    csr_zeta_selMode_rw(GET_ZETA_MODE(1, 5));
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/2=16 zetas
        v_zetas = vle16_v_i16m1(zeta_ptr, 16);
        zeta_ptr += 16;
        v_in0 = vle16_v_i16m1(r_ptr0, 32);
        v_in1 = vle16_v_i16m1(r_ptr1, 32);
        r_ptr0 += 32;
        r_ptr1 += 32;
        keccak_squeeze();
        keccak_squeeze();
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        keccak_squeeze();
        keccak_squeeze();
        vse16_v_i16m2(r_ptr2, v_out0, 64);
        r_ptr2 += 64;
    }

    // stage 2
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r;
    // same_num = 4
    csr_zeta_selMode_rw(GET_ZETA_MODE(2, 5));
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/4=8 zetas
        v_zetas = vle16_v_i16m1(zeta_ptr, 8);
        zeta_ptr += 8;
        v_in0 = vle16_v_i16m1(r_ptr0, 32);
        v_in1 = vle16_v_i16m1(r_ptr1, 32);
        r_ptr0 += 32;
        r_ptr1 += 32;
        keccak_squeeze();
        keccak_squeeze();
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        keccak_squeeze();
        keccak_squeeze();
        vse16_v_i16m2(r_ptr2, v_out0, 64);
        r_ptr2 += 64;
    }

    // stage 3
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    // same_num = 8
    csr_zeta_selMode_rw(GET_ZETA_MODE(3, 5));  
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/8=4 zetas
        v_zetas = vle16_v_i16m1(zeta_ptr, 4);
        zeta_ptr += 4;
        v_in0 = vle16_v_i16m1(r_ptr0, 32);
        v_in1 = vle16_v_i16m1(r_ptr1, 32);
        r_ptr0 += 32;
        r_ptr1 += 32;
        keccak_squeeze();
        keccak_squeeze();
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        keccak_squeeze();
        keccak_squeeze();
        vse16_v_i16m2(r_ptr2, v_out0, 64);
        r_ptr2 += 64;
    }

    // stage 4
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r;
    // same_num = 16
    csr_zeta_selMode_rw(GET_ZETA_MODE(4, 5));  
    for(i = 0; i < 4; i++) {                // each 32* butterfly, 32/16=2 zetas
        v_zetas = vle16_v_i16m1(zeta_ptr, 2);
        zeta_ptr += 2;
        v_in0 = vle16_v_i16m1(r_ptr0, 32);
        v_in1 = vle16_v_i16m1(r_ptr1, 32);
        r_ptr0 += 32;
        r_ptr1 += 32;
        keccak_squeeze();
        keccak_squeeze();
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        keccak_squeeze();
        keccak_squeeze();
        vse16_v_i16m2(r_ptr2, v_out0, 64);
        r_ptr2 += 64;
    }

    // stage 5
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    // same_num = 32
    csr_zeta_selMode_rw(GET_ZETA_MODE(5, 5));
    for(i = 0; i < 4; i++) {                // each 32* butterfly, all same zeta
        v_zetas = vle16_v_i16m1(zeta_ptr, 1);
        zeta_ptr += 1;
        v_in0 = vle16_v_i16m1(r_ptr0, 32);
        v_in1 = vle16_v_i16m1(r_ptr1, 32);
        r_ptr0 += 32;
        r_ptr1 += 32;
        keccak_squeeze();
        keccak_squeeze();
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        keccak_squeeze();
        keccak_squeeze();
        vse16_v_i16m2(r_ptr2, v_out0, 64);
        r_ptr2 += 64;
    }

    // stage 6
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r_temp1;
    // same_num = 64
    csr_zeta_selMode_rw(GET_ZETA_MODE(6, 5));
    for(i = 0; i < 2; i++) {              // each 64* butterfly, all same zeta
        v_zetas = vle16_v_i16m1(zeta_ptr, 1);
        zeta_ptr += 1;
        
        v_in2 = vle16_v_i16m2(r_ptr0, 64);
        v_in3 = vle16_v_i16m2(r_ptr1, 64);
        r_ptr0 += 64;
        r_ptr1 += 64;
        keccak_squeeze();
        keccak_squeeze();
        v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
        keccak_squeeze();
        keccak_squeeze();
        vse16_v_i16m4(r_ptr2, v_out1, 128);
        r_ptr2 += 128;
    }

    // indexed store: used as offsets
    vuint16m4_t v_offset;
    // tree pointer
    const uint16_t* tree_prt = NULL;

    // stage 7
    // same_num = 128
    csr_zeta_selMode_rw(GET_ZETA_MODE(7, 5));
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    r_ptr0 = r_temp1;
    r_ptr1 = &r_temp1[128];
    tree_prt = byteoffset_even;
    v_offset = vle16_v_u16m4(tree_prt, 128);
    for(j = 0; j < 2; j++) {                // each 64* butterfly, all same zeta
        v_in2 = vle16_v_i16m2(r_ptr0, 64);
        v_in3 = vle16_v_i16m2(r_ptr1, 64);
        r_ptr0 += 64;
        r_ptr1 += 64;
        keccak_squeeze();
        keccak_squeeze();
        v_out1 = vbutterfly_gs_vvm_i16m4(v_in2, v_in3, v_zetas);
        keccak_squeeze();
        keccak_squeeze();
        // reorder, from special bit-reversed order to standard order
        vsuxei16_v_i16m4(r, v_offset, v_out1, 128);
        v_offset = vadd_vx_u16m4(v_offset, 2, 128);
    }    
}

void invntt_cg_custom_reorder_shuffling_ultrafine(int16_t r[KYBER_N],const uint16_t srcpoly_idx[128*7], const uint16_t dstpoly_idx[256*7]) {
    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);
    csr_zeta_selMode_rw(GET_ZETA_MODE(0, 5));

    int16_t r_temp0[256];
    int16_t r_temp1[256];
    unsigned int i, j;
    int16_t* r_ptr0 = NULL;
    int16_t* r_ptr1 = NULL;
    int16_t* r_ptr2 = NULL;
    const int16_t* zeta_ptr=NULL;
    const uint16_t* src_idx_ptr=NULL;
    const uint16_t* dst_idx_ptr=NULL;

    vint16m1_t v_zetas;
    vint16m1_t v_in0, v_in1;
    vint16m2_t v_out0;
    vuint16m1_t vsrcidx;
    vuint16m2_t vdstidx;

    // stage 1
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    zeta_ptr = zetas_inv_in_order_stage1;
    src_idx_ptr=srcpoly_idx;
    dst_idx_ptr=dstpoly_idx;
    for(i = 0; i < 4; i++) {
        vsrcidx=vle16_v_u16m1(src_idx_ptr,32);
        vdstidx=vle16_v_u16m2(dst_idx_ptr,64);
        v_zetas = vluxei16_v_i16m1(zeta_ptr,vsrcidx,32);
        v_in0 = vluxei16_v_i16m1(r_ptr0,vsrcidx,32);
        v_in1 = vluxei16_v_i16m1(r_ptr1,vsrcidx,32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vsuxei16_v_i16m2(r_ptr2,vdstidx,v_out0,64);
        src_idx_ptr+=32;
        dst_idx_ptr+=64;
    }

    // stage 2
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r;
    zeta_ptr = zetas_inv_in_order_stage2;
    src_idx_ptr=&srcpoly_idx[128*1];
    dst_idx_ptr=&dstpoly_idx[256*1];
    for(i = 0; i < 4; i++) {
        vsrcidx=vle16_v_u16m1(src_idx_ptr,32);
        vdstidx=vle16_v_u16m2(dst_idx_ptr,64);
        v_zetas = vluxei16_v_i16m1(zeta_ptr,vsrcidx,32);
        v_in0 = vluxei16_v_i16m1(r_ptr0,vsrcidx,32);
        v_in1 = vluxei16_v_i16m1(r_ptr1,vsrcidx,32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vsuxei16_v_i16m2(r_ptr2,vdstidx,v_out0,64);
        src_idx_ptr+=32;
        dst_idx_ptr+=64;
    }

    // stage 3
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    zeta_ptr = zetas_inv_in_order_stage3;
    src_idx_ptr=&srcpoly_idx[128*2];
    dst_idx_ptr=&dstpoly_idx[256*2];
    for(i = 0; i < 4; i++) {
        vsrcidx=vle16_v_u16m1(src_idx_ptr,32);
        vdstidx=vle16_v_u16m2(dst_idx_ptr,64);
        v_zetas = vluxei16_v_i16m1(zeta_ptr,vsrcidx,32);
        v_in0 = vluxei16_v_i16m1(r_ptr0,vsrcidx,32);
        v_in1 = vluxei16_v_i16m1(r_ptr1,vsrcidx,32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vsuxei16_v_i16m2(r_ptr2,vdstidx,v_out0,64);
        src_idx_ptr+=32;
        dst_idx_ptr+=64;
    }

    // stage 4
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r;
    zeta_ptr = zetas_inv_in_order_stage4;
    src_idx_ptr=&srcpoly_idx[128*3];
    dst_idx_ptr=&dstpoly_idx[256*3];
    for(i = 0; i < 4; i++) {
        vsrcidx=vle16_v_u16m1(src_idx_ptr,32);
        vdstidx=vle16_v_u16m2(dst_idx_ptr,64);
        v_zetas = vluxei16_v_i16m1(zeta_ptr,vsrcidx,32);
        v_in0 = vluxei16_v_i16m1(r_ptr0,vsrcidx,32);
        v_in1 = vluxei16_v_i16m1(r_ptr1,vsrcidx,32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vsuxei16_v_i16m2(r_ptr2,vdstidx,v_out0,64);
        src_idx_ptr+=32;
        dst_idx_ptr+=64;
    }

    // stage 5
    r_ptr0 = r;
    r_ptr1 = &r[128];
    r_ptr2 = r_temp0;
    zeta_ptr = zetas_inv_in_order_stage5;
    src_idx_ptr=&srcpoly_idx[128*4];
    dst_idx_ptr=&dstpoly_idx[256*4];
    for(i = 0; i < 4; i++) {
        vsrcidx=vle16_v_u16m1(src_idx_ptr,32);
        vdstidx=vle16_v_u16m2(dst_idx_ptr,64);
        v_zetas = vluxei16_v_i16m1(zeta_ptr,vsrcidx,32);
        v_in0 = vluxei16_v_i16m1(r_ptr0,vsrcidx,32);
        v_in1 = vluxei16_v_i16m1(r_ptr1,vsrcidx,32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vsuxei16_v_i16m2(r_ptr2,vdstidx,v_out0,64);
        src_idx_ptr+=32;
        dst_idx_ptr+=64;
    }

    // stage 6
    r_ptr0 = r_temp0;
    r_ptr1 = &r_temp0[128];
    r_ptr2 = r_temp1;
    zeta_ptr = zetas_inv_in_order_stage6;
    src_idx_ptr=&srcpoly_idx[128*5];
    dst_idx_ptr=&dstpoly_idx[256*5];
    for(i = 0; i < 4; i++) {
        vsrcidx=vle16_v_u16m1(src_idx_ptr,32);
        vdstidx=vle16_v_u16m2(dst_idx_ptr,64);
        v_zetas = vluxei16_v_i16m1(zeta_ptr,vsrcidx,32);
        v_in0 = vluxei16_v_i16m1(r_ptr0,vsrcidx,32);
        v_in1 = vluxei16_v_i16m1(r_ptr1,vsrcidx,32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vsuxei16_v_i16m2(r_ptr2,vdstidx,v_out0,64);
        src_idx_ptr+=32;
        dst_idx_ptr+=64;
    }

    // stage 7
    zeta_ptr=zetas_inv_in_order_stage7;
    src_idx_ptr=&srcpoly_idx[128*6];
    dst_idx_ptr=&dstpoly_idx[256*6];
    csr_zeta_selMode_rw(GET_ZETA_MODE(7, 5));
    v_zetas = vle16_v_i16m1(zeta_ptr, 1);
    r_ptr0 = r_temp1;
    r_ptr1 = &r_temp1[128];
    for(j = 0; j < 4; j++) {
        vsrcidx=vle16_v_u16m1(src_idx_ptr,32);
        vdstidx=vle16_v_u16m2(dst_idx_ptr,64);
        v_in0 = vluxei16_v_i16m1(r_ptr0,vsrcidx,32);
        v_in1 = vluxei16_v_i16m1(r_ptr1,vsrcidx,32);
        v_out0 = vbutterfly_gs_vvm_i16m2(v_in0, v_in1, v_zetas);
        vsuxei16_v_i16m2(r,vdstidx,v_out0,64);
        src_idx_ptr+=32;
        dst_idx_ptr+=64;
    }
}

#elif (VLEN == 1024)
# else
#error "VLEN must be 256/512/1024"
#endif

void test_invntt_cg_custom_reorder_sca(){
    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);
    srand(time(NULL));
    uint8_t seed[32];
    for(int i=0;i<32;i++)seed[i]=rand()+2;

    int16_t r_ref[KYBER_N];
    int16_t r_res[KYBER_N];
    for(int i=0;i<KYBER_N;i++){
        r_ref[i]=r_res[i]=rand()%KYBER_Q;
    }

    invntt_cg_custom_reorder(r_ref);
    //invntt_cg_custom_reorder_shuffling(r_res,seed);
    //invntt_cg_custom_reorder_shuffling_fine(r_res,seed);
    invntt_cg_custom_reorder_redundant(r_res,seed);

    //check
    bool flag=true;
    for(int i=0;i<KYBER_N;i++){
        if(r_ref[i]!=r_res[i]){
            flag=false;
            break;
            //printf("r_ref[%d]=%d,r_res[%d]=%d\n",i,r_ref[i],i,r_res[i]);
        }
    }
    if(flag) {
        printf("test_invntt_cg_custom_reorder_shuffling test pass!!!\n");
    }
    else {
        printf("test_invntt_cg_custom_reorder_shuffling test fail...\n");
    }
    return;
}

void try_generate_intt_shuffling_idx(){
    uint16_t srcpoly_idx[128*7];
    uint16_t dstpoly_idx[256*7];
    uint8_t prng_buf[128*7];
    uint8_t seed[32];
    for(int i=0;i<32;i++)seed[i]=rand();
    shake128_custom(prng_buf,128*7,seed,32);
    generate_intt_shuffling_idx(srcpoly_idx,dstpoly_idx,prng_buf);
}

void test_invntt_cg_custom_reorder_shuffling_ultrafine(){
    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);
    uint16_t srcpoly_idx[128*7];
    uint16_t dstpoly_idx[256*7];
    uint8_t prng_buf[128*7];
    srand(time(NULL));
    uint8_t seed[32];
    for(int i=0;i<32;i++)seed[i]=rand()+2;
    shake128_custom(prng_buf,128*7,seed,32);
    generate_intt_shuffling_idx(srcpoly_idx,dstpoly_idx,prng_buf);

    int16_t r_ref[KYBER_N];
    int16_t r_res[KYBER_N];
    for(int i=0;i<KYBER_N;i++){
        r_ref[i]=r_res[i]=rand()%KYBER_Q;
    }

    invntt_cg_custom_reorder(r_ref);
    invntt_cg_custom_reorder_shuffling_ultrafine(r_res,srcpoly_idx,dstpoly_idx);

    //check
    bool flag=true;
    for(int i=0;i<KYBER_N;i++){
        if(r_ref[i]!=r_res[i]){
            flag=false;
            break;
            //printf("r_ref[%d]=%d,r_res[%d]=%d\n",i,r_ref[i],i,r_res[i]);
        }
    }
    if(flag) {
        printf("test_invntt_cg_custom_reorder_shuffling_ultrafine test pass!!!\n");
    }
    else {
        printf("test_invntt_cg_custom_reorder_shuffling_ultrafine test fail...\n");
    }
    return;
}

int main(){
    //test_invntt_cg_custom_reorder_sca();
    //try_generate_intt_shuffling_idx();
    test_invntt_cg_custom_reorder_shuffling_ultrafine();
}