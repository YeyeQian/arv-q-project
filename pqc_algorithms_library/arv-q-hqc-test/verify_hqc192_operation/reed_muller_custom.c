/**
 * @file reed_muller.cpp
 * @brief Constant time implementation of Reed-Muller code RM(1,7)
 */

#include "reed_muller.h"
#include "parameters.h"
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

// number of repeated code words
#define MULTIPLICITY                   CEIL_DIVIDE(PARAM_N2, 128)

// codeword is 128 bits, seen multiple ways
typedef union {
    uint8_t u8[16];
    uint32_t u32[4];
} codeword;

// Expanded codeword has a short for every bit, for internal calculations
typedef int16_t expandedCodeword[128];

// copy bit 0 into all bits of a 32 bit value
#define BIT0MASK(x) (int32_t)(-((x) & 1))


void encode(codeword *word, int32_t message);
void hadamard(expandedCodeword *src, expandedCodeword *dst);
void expand_and_sum(expandedCodeword *dest, codeword src[]);
int32_t find_peaks(expandedCodeword *transform);



/**
 * @brief Encode a single byte into a single codeword using RM(1,7)
 *
 * Encoding matrix of this code:
 * bit pattern (note that bits are numbered big endian)
 * 0   aaaaaaaa aaaaaaaa aaaaaaaa aaaaaaaa
 * 1   cccccccc cccccccc cccccccc cccccccc
 * 2   f0f0f0f0 f0f0f0f0 f0f0f0f0 f0f0f0f0
 * 3   ff00ff00 ff00ff00 ff00ff00 ff00ff00
 * 4   ffff0000 ffff0000 ffff0000 ffff0000
 * 5   ffffffff 00000000 ffffffff 00000000
 * 6   ffffffff ffffffff 00000000 00000000
 * 7   ffffffff ffffffff ffffffff ffffffff
 *
 * @param[out] word An RM(1,7) codeword
 * @param[in] message A message
 */
void encode(codeword *word, int32_t message) {
    int32_t first_word;

    first_word = BIT0MASK(message >> 7);

    first_word ^= BIT0MASK(message >> 0) & 0xaaaaaaaa;
    first_word ^= BIT0MASK(message >> 1) & 0xcccccccc;
    first_word ^= BIT0MASK(message >> 2) & 0xf0f0f0f0;
    first_word ^= BIT0MASK(message >> 3) & 0xff00ff00;
    first_word ^= BIT0MASK(message >> 4) & 0xffff0000;

    word->u32[0] = first_word;

    first_word ^= BIT0MASK(message >> 5);
    word->u32[1] = first_word;
    first_word ^= BIT0MASK(message >> 6);
    word->u32[3] = first_word;
    first_word ^= BIT0MASK(message >> 5);
    word->u32[2] = first_word;
    return;
}



/**
 * @brief Hadamard transform
 *
 * Perform hadamard transform of src and store result in dst
 * src is overwritten
 *
 * @param[out] src Structure that contain the expanded codeword
 * @param[out] dst Structure that contain the expanded codeword
 */
void hadamard(expandedCodeword *src, expandedCodeword *dst) {
    // the passes move data:
    // src -> dst -> src -> dst -> src -> dst -> src -> dst
    // using p1 and p2 alternately
    expandedCodeword *p1 = src;
    expandedCodeword *p2 = dst;
    for (int32_t pass = 0; pass < 7; pass++) {
        for (int32_t i = 0; i < 64; i++) {
            (*p2)[i] = (*p1)[2 * i] + (*p1)[2 * i + 1];
            (*p2)[i + 64] = (*p1)[2 * i] - (*p1)[2 * i + 1];
        }
        // swap p1, p2 for next round
        expandedCodeword *p3 = p1;
        p1 = p2;
        p2 = p3;
    }
}



/**
 * @brief Add multiple codewords into expanded codeword
 *
 * Accesses memory in order
 * Note: this does not write the codewords as -1 or +1 as the green machine does
 * instead, just 0 and 1 is used.
 * The resulting hadamard transform has:
 * all values are halved
 * the first entry is 64 too high
 *
 * @param[out] dest Structure that contain the expanded codeword
 * @param[in] src Structure that contain the codeword
 */
void expand_and_sum(expandedCodeword *dest, codeword src[]) {
    // start with the first copy
    for (int32_t part = 0; part < 4; part++) {
        for (int32_t bit = 0; bit < 32; bit++) {
            (*dest)[part * 32 + bit] = src[0].u32[part] >> bit & 1;
        }
    }
    // sum the rest of the copies
    for (int32_t copy = 1; copy < MULTIPLICITY; copy++) {
        for (int32_t part = 0; part < 4; part++) {
            for (int32_t bit = 0; bit < 32; bit++) {
                (*dest)[part * 32 + bit] += src[copy].u32[part] >> bit & 1;
            }
        }
    }
}



/**
 * @brief Finding the location of the highest value
 *
 * This is the final step of the green machine: find the location of the highest value,
 * and add 128 if the peak is positive
 * if there are two identical peaks, the peak with smallest value
 * in the lowest 7 bits it taken
 * @param[in] transform Structure that contain the expanded codeword
 */
int32_t find_peaks(expandedCodeword *transform) {
    int32_t peak_abs_value = 0;
    int32_t peak_value = 0;
    int32_t peak_pos = 0;
    for (int32_t i = 0; i < 128; i++) {
        // get absolute value
        int32_t t = (*transform)[i];
        int32_t pos_mask = -(t > 0);
        int32_t absolute = (pos_mask & t) | (~pos_mask & -t);
        // all compilers nowadays compile with a conditional move
        peak_value = absolute > peak_abs_value ? t : peak_value;
        peak_pos = absolute > peak_abs_value ? i : peak_pos;
        peak_abs_value = absolute > peak_abs_value ? absolute : peak_abs_value;
    }
    // set bit 7
    peak_pos |= 128 * (peak_value > 0);
    return peak_pos;
}



/**
 * @brief Encodes the received word
 *
 * The message consists of N1 bytes each byte is encoded into PARAM_N2 bits,
 * or MULTIPLICITY repeats of 128 bits
 *
 * @param[out] cdw Array of size VEC_N1N2_SIZE_64 receiving the encoded message
 * @param[in] msg Array of size VEC_N1_SIZE_64 storing the message
 */
void reed_muller_encode(uint64_t *cdw, const uint64_t *msg) {
    uint8_t *message_array = (uint8_t *) msg;
    codeword *codeArray = (codeword *) cdw;
    for (size_t i = 0; i < VEC_N1_SIZE_BYTES; i++) {
        // fill entries i * MULTIPLICITY to (i+1) * MULTIPLICITY
        int32_t pos = i * MULTIPLICITY;
        // encode first word
        encode(&codeArray[pos], message_array[i]);
        // copy to other identical codewords
        for (size_t copy = 1; copy < MULTIPLICITY; copy++) {
            memcpy(&codeArray[pos + copy], &codeArray[pos], sizeof(codeword));
        }
    }
    return;
}



/**
 * @brief Decodes the received word
 *
 * Decoding uses fast hadamard transform, for a more complete picture on Reed-Muller decoding, see MacWilliams, Florence Jessie, and Neil James Alexander Sloane.
 * The theory of error-correcting codes codes @cite macwilliams1977theory
 *
 * @param[out] msg Array of size VEC_N1_SIZE_64 receiving the decoded message
 * @param[in] cdw Array of size VEC_N1N2_SIZE_64 storing the received word
 */
void reed_muller_decode(uint64_t *msg, const uint64_t *cdw) {
    uint8_t *message_array = (uint8_t *) msg;
    codeword *codeArray = (codeword *) cdw;
    expandedCodeword expanded;
    for (size_t i = 0; i < VEC_N1_SIZE_BYTES; i++) {
        // collect the codewords
        expand_and_sum(&expanded, &codeArray[i * MULTIPLICITY]);
        // apply hadamard transform
        expandedCodeword transform;
        hadamard(&expanded, &transform);
        // fix the first entry to get the half Hadamard transform
        transform[0] -= 64 * MULTIPLICITY;
        // finish the decoding
        message_array[i] = find_peaks(&transform);
    }
}

///////////////////////////////////////////////////////////////////////////
//                      Customized Versions
///////////////////////////////////////////////////////////////////////////

void encode_custom(codeword *word, int32_t message) {//VLEN Must Be Greater Than 256bits
    uint32_t column0[8]={0xaaaaaaaa,0xcccccccc,0x0f0f0f0f,0x00ff00ff,0x0000ffff,0xffffffff,0xffffffff,0xffffffff};
    uint32_t column1[8]={0xaaaaaaaa,0xcccccccc,0x0f0f0f0f,0x00ff00ff,0x0000ffff,0x00000000,0xffffffff,0xffffffff};
    uint32_t column2[8]={0xaaaaaaaa,0xcccccccc,0x0f0f0f0f,0x00ff00ff,0x0000ffff,0xffffffff,0x00000000,0xffffffff};
    uint32_t column3[8]={0xaaaaaaaa,0xcccccccc,0x0f0f0f0f,0x00ff00ff,0x0000ffff,0x00000000,0x00000000,0xffffffff};

    size_t vl,avl;
    avl=8;
    vuint32m1_t va,vb;
    vbool32_t vc;
    uint8_t* msg_addr=NULL;
    msg_addr=(uint8_t*)&message;
    vl=vsetvl_e8m1(1);
    vc=vlm_v_b32(msg_addr,vl);

    vl=vsetvl_e32m1(8);
    //first column process
    va=vmv_s_x_u32m1(va,0,vl);
    vb=vle32_v_u32m1(column0,vl);
    va=vredxor_vs_u32m1_u32m1_m(vc,va,vb,va,vl);
    word->u32[0]=vmv_x_s_u32m1_u32(va);

    //second column process
    va=vmv_s_x_u32m1(va,0,vl);
    vb=vle32_v_u32m1(column1,vl);
    va=vredxor_vs_u32m1_u32m1_m(vc,va,vb,va,vl);
    word->u32[1]=vmv_x_s_u32m1_u32(va);

    //Third column process
    va=vmv_s_x_u32m1(va,0,vl);
    vb=vle32_v_u32m1(column2,vl);
    va=vredxor_vs_u32m1_u32m1_m(vc,va,vb,va,vl);
    word->u32[2]=vmv_x_s_u32m1_u32(va);

    //Forth column process
    va=vmv_s_x_u32m1(va,0,vl);
    vb=vle32_v_u32m1(column3,vl);
    va=vredxor_vs_u32m1_u32m1_m(vc,va,vb,va,vl);
    word->u32[3]=vmv_x_s_u32m1_u32(va);

    return;
}

void hadamard_custom(expandedCodeword *src, expandedCodeword *dst) {
    // the passes move data:
    // src -> dst -> src -> dst -> src -> dst -> src -> dst
    // using p1 and p2 alternately
    expandedCodeword *p1 = src;
    expandedCodeword *p2 = dst;
    size_t vl,avl;
    vint16m4_t va,vb,vc;
    vuint16m4_t vidx1,vidx2;
    for (int32_t pass = 0; pass < 7; pass++) {
        avl=64;
        int16_t* p1_addr=&(*p1)[0];
        int16_t* p2_addr=&(*p2)[0];
        while(avl>0){
            vl=vsetvl_e16m4(avl);
            vidx1=vid_v_u16m4(vl);
            vidx1=vadd_vx_u16m4(vidx1,64-avl,vl);
            vidx1=vmul_vx_u16m4(vidx1,2,vl);//2*i
            vidx2=vadd_vx_u16m4(vidx1,1,vl);//2*i+1
            vidx1=vsll_vx_u16m4(vidx1,1,vl);
            vidx2=vsll_vx_u16m4(vidx2,1,vl);
            va=vluxei16_v_i16m4(p1_addr,vidx1,vl);//(*p1)[2 * i]
            vb=vluxei16_v_i16m4(p1_addr,vidx2,vl);//(*p1)[2 * i + 1]
            vc=vadd_vv_i16m4(va,vb,vl);
            vse16_v_i16m4(p2_addr,vc,vl);
            vc=vsub_vv_i16m4(va,vb,vl);
            vse16_v_i16m4(p2_addr+64,vc,vl);

            p2_addr+=vl,avl-=vl;
        }
        // swap p1, p2 for next round
        expandedCodeword *p3 = p1;
        p1 = p2;
        p2 = p3;
    }
}

void expand_and_sum_custom(expandedCodeword *dest, codeword src[]) {//VLEN Must Be Greater Than 256bits
    // size_t vl,avl;
    // vint16m2_t va,vb,vres;
    // vbool8_t vflag;
    // int16_t* dest_addr=&(*dest)[0];

    // int16_t buf[32]={0};//for test

    // for(int part=0;part<4;part++){
    //     vl=vsetvl_e16m2(32);
    //     vres=vmv_s_x_i16m2(vres,0,vl);
    //     printf("part%d:\n",part);//for test
    //     for(int i=0;i<MULTIPLICITY;i++){
    //         vl=vsetvl_e16m2(32);
    //         va=vand_vx_i16m2(va,0,vl);
    //         vflag=vlm_v_b8((uint8_t*)&src[i].u32[part],vl);
    //         va=vor_vx_i16m2_m(vflag,va,va,1,vl);
    //             printf("src[%d].u32[%d]=%x\n",i,part,src[i].u32[part]);// for test
    //             vse16_v_i16m2(buf,va,vl);//for test
    //             for(int i=31;i>=0;i--)printf("%d,",buf[i]);//for test
    //             printf("\n");//for test
    //         vres=vadd_vv_i16m2(vres,va,vl);
    //     }
    //     vse16_v_i16m2(dest_addr+(part<<5),vres,vl);
    // }
    // start with the first copy
    for (int32_t part = 0; part < 4; part++) {
        for (int32_t bit = 0; bit < 32; bit++) {
            (*dest)[part * 32 + bit] = src[0].u32[part] >> bit & 1;
        }
    }
    // sum the rest of the copies
    for (int32_t copy = 1; copy < MULTIPLICITY; copy++) {
        for (int32_t part = 0; part < 4; part++) {
            for (int32_t bit = 0; bit < 32; bit++) {
                (*dest)[part * 32 + bit] += src[copy].u32[part] >> bit & 1;
            }
        }
    }
}

void reed_muller_encode_custom(uint64_t *cdw, const uint64_t *msg) {
    uint8_t *message_array = (uint8_t *) msg;
    codeword *codeArray = (codeword *) cdw;
    for (size_t i = 0; i < VEC_N1_SIZE_BYTES; i++) {
        // fill entries i * MULTIPLICITY to (i+1) * MULTIPLICITY
        int32_t pos = i * MULTIPLICITY;
        // encode first word
        encode(&codeArray[pos], message_array[i]);
        // copy to other identical codewords
        for (size_t copy = 1; copy < MULTIPLICITY; copy++) {
            memcpy(&codeArray[pos + copy], &codeArray[pos], sizeof(codeword));
        }
    }
    return;
}

void reed_muller_decode_custom(uint64_t *msg, const uint64_t *cdw) {
    uint8_t *message_array = (uint8_t *) msg;
    codeword *codeArray = (codeword *) cdw;
    expandedCodeword expanded;
    for (size_t i = 0; i < VEC_N1_SIZE_BYTES; i++) {
        // collect the codewords
        expand_and_sum(&expanded, &codeArray[i * MULTIPLICITY]);
        // apply hadamard transform
        expandedCodeword transform;
        hadamard_custom(&expanded, &transform);
        // fix the first entry to get the half Hadamard transform
        transform[0] -= 64 * MULTIPLICITY;
        // finish the decoding
        message_array[i] = find_peaks(&transform);
    }
}

///////////////////////////////////////////////////////////////////////////
//                               Tests
///////////////////////////////////////////////////////////////////////////
bool test_encode_custom(){
    codeword ref,res;
    uint8_t msg;
    srand((unsigned)time(NULL));
    msg=rand()&255;

    uint64_t start, end;

    start=read_cycle();
    encode(&ref,msg);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    encode_custom(&res,msg);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<(4);i++){
        if(ref.u32[i]!=res.u32[i]){
            flag=false;
            //break;
            printf("ref[%d]=%lu,res[%d]=%lu\n",i,ref.u32[i],i,res.u32[i]);
        }
    }

    printf("test_encode_custom return with flag %d\n",flag);

    return flag;
}


bool test_expand_and_sum_custom(){
    expandedCodeword expanded_ref;
    expandedCodeword expanded_res;

    codeword codeArray[MULTIPLICITY];
    srand((unsigned)time(NULL));
    for(int i=0;i<MULTIPLICITY;i++){
        for(int j=0;j<4;j++){
            codeArray[i].u32[j]=((uint32_t)rand()<<16)|((uint32_t)rand());
        }
    }

    uint64_t start, end;

    start=read_cycle();
    expand_and_sum(&expanded_ref, &codeArray[0]);
    end=read_cycle();
    printf("Reference finished with %lu cycles\n",end-start);

    start=read_cycle();
    expand_and_sum_custom(&expanded_res, &codeArray[0]);
    end=read_cycle();
    printf("Custom finished with %lu cycles\n",end-start);

    //CHECK
    bool flag=true;
    for(int i=0;i<128;i++){
        if(expanded_ref[i]!=expanded_res[i]){
            flag=false;
            break;
            // printf("ref[%d]=%lu,res[%d]=%lu\n",i,expanded_ref[i],i,expanded_res[i]);
        }
    }

    printf("test_expand_and_sum_custom return with flag %d\n",flag);

    return flag;
}

// int main(){
//     // test_encode_custom();
//     test_expand_and_sum_custom();

//     return 0;
// }