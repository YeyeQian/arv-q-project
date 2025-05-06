#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include "gf.h"
#include "parameters.h"
#include <stdint.h>

void gfmul_vx_custom_u8(uint8_t* dst_v, const uint8_t* src_v, uint8_t s, uint32_t len){
    csr_primpoly_rw(PARAM_GF_POLY);

    size_t vl,avl;
    avl=len;

    vuint8m1_t va,vb;

    while(avl>0){
        vl=vsetvl_e8m1(avl);
        va=vle8_v_u8m1(src_v,vl);
        vb=vgfmul_vx_u8m1(va,s);
        vse8_v_u8m1(dst_v,vb,vl);

        dst_v+=vl,src_v+=vl,avl-=vl;
    }
}

void gfmul_vx_custom_u16(uint16_t* dst_v, const uint16_t* src_v, uint16_t s, uint32_t len){
    csr_primpoly_rw(PARAM_GF_POLY);

    size_t vl,avl;
    avl=len;

    vuint16m2_t va,vb;

    while(avl>0){
        vl=vsetvl_e16m2(avl);
        va=vle16_v_u16m2(src_v,vl);
        vb=vgfmul_vx_u16m2(va,s);
        vse16_v_u16m2(dst_v,vb,vl);

        dst_v+=vl,src_v+=vl,avl-=vl;
    }
}

/**********************************************************************
 *                Binary Polynomail Multiplication
 *                      (Only Support SEW=8)
**********************************************************************/
void cyclic_shift_left_custom(uint8_t poly_res[VEC_N_SIZE_BYTES], const uint8_t poly_dense[VEC_N_SIZE_BYTES], const uint8_t addi,uint16_t offset){
    size_t vl;
    
    uint8_t addi_data;
    int32_t addi_idx=-1;

    const uint8_t* inspect_addr=poly_dense;

    uint32_t chunk_shift=offset>>3;
    uint8_t* res_chunk_addr=poly_res+chunk_shift;

    int32_t remains_num=VEC_N_SIZE_BYTES-chunk_shift;

    uint8_t bit_shift1=offset&7;
    uint8_t bit_shift2=bit_shift1+8-(PARAM_N&7);
    bool flag=bit_shift2>=8;
    bit_shift2=flag?bit_shift2-8:bit_shift2;

    //configure the BSeqSftOffset CSR
    uint8_t old_offset;//no useful meaning
    old_offset=csr_bseqoffset_rw(bit_shift1);

    while(remains_num>0){
        vl=vsetvl_e8m1(remains_num);
        addi_data=addi_idx==-1?addi:poly_dense[addi_idx];
        vuint8m1_t v_inspect=vle8_v_u8m1(inspect_addr,vl);
        vuint8m1_t v_cycres;
        v_cycres=bitseqsll_vx_u8m1(v_inspect,addi_data);
        vuint8m1_t v_res=vle8_v_u8m1(res_chunk_addr,vl);
        v_res=vxor_vv_u8m1(v_res,v_cycres,vl);
        vse8_v_u8m1(res_chunk_addr,v_res,vl);

        addi_idx+=vl;
        inspect_addr+=vl;
        res_chunk_addr+=vl;
        remains_num-=vl;
    }

    remains_num=chunk_shift;
    res_chunk_addr=poly_res;
    if(flag){
        addi_idx-=1;
        inspect_addr-=1;
    }
    old_offset=csr_bseqoffset_rw(bit_shift2);
    while(remains_num>0){
        vl=vsetvl_e8m1(remains_num);
        addi_data=poly_dense[addi_idx];
        vuint8m1_t v_inspect=vle8_v_u8m1(inspect_addr,vl);
        vuint8m1_t v_cycres;
        v_cycres=bitseqsll_vx_u8m1(v_inspect,addi_data);
        vuint8m1_t v_res=vle8_v_u8m1(res_chunk_addr,vl);
        v_res=vxor_vv_u8m1(v_res,v_cycres,vl);
        vse8_v_u8m1(res_chunk_addr,v_res,vl);

        addi_idx+=vl;
        inspect_addr+=vl;
        res_chunk_addr+=vl;
        remains_num-=vl;
    }
}

void vect_mul_cycshift_custom(uint64_t* out, const uint64_t* dense, const uint32_t* support, const uint32_t weight) {
    uint8_t last_chunk_bits = PARAM_N & (8 - 1);
    const uint8_t* dense_byte=NULL;
    dense_byte=(uint8_t*)dense;
    uint8_t* out_byte=NULL;
    out_byte=(uint8_t*)out;
    uint8_t addi = (((uint16_t)dense_byte[VEC_N_SIZE_BYTES-1] << 8 | (uint16_t)dense_byte[VEC_N_SIZE_BYTES - 2]) >> last_chunk_bits) & 255;

    for (int i = 0; i < weight; i++) {
        cyclic_shift_left_custom(out_byte,dense_byte,addi,support[i]);
    }

    uint8_t final_mask=((1 << (PARAM_N & 7)) - 1);
    out_byte[VEC_N_SIZE_BYTES-1]&=final_mask;
}
