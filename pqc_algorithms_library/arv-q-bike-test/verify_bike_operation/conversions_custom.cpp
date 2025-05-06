#include "types.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include <stdio.h>


// convert a sequence of uint8_t elements which fully uses all 8-bits of an uint8_t element to
// a sequence of uint8_t which uses just a single bit per byte (either 0 or 1).
int convertByteToBinary_custom(uint8_t* out, const uint8_t * in, uint32_t length)
{

    int32_t avl=length;
    size_t vl;
    const uint8_t* vsrc_addr=in;
    uint8_t* vdst_addr=out;
    while(avl>0){
        vl=vsetvl_e8m1(avl);
        vbool8_t vmask=vlm_v_b8(vsrc_addr,vl);
        vuint8m1_t vres=vmv_v_x_u8m1(0,vl);
        vres=vadd_vx_u8m1_m(vmask,vres,vres,1,vl);
        vse8_v_u8m1(vdst_addr,vres,vl);

        vsrc_addr+=(vl>>3),vdst_addr+=vl,avl-=vl;
    }
    return 0;
}

// convert a sequence of uint8_t elements which uses just a single bit per byte (either 0 or 1) to
// a sequence of uint8_t which fully uses all 8-bits of an uint8_t element.
int convertBinaryToByte_custom(uint8_t * out, const uint8_t* in, uint32_t length)
{

    int32_t avl=length;
    size_t vl;
    const uint8_t* vsrc_addr=in;
    uint8_t* vdst_addr=out;
    while(avl>0){
        vl=vsetvl_e8m1(avl);
        vuint8m1_t vsrc=vle8_v_u8m1(vsrc_addr,vl);
        vbool8_t vmask=vmseq_vx_u8m1_b8(vsrc,1,vl);
        vsm_v_b8(vdst_addr,vmask,vl);

        vsrc_addr+=vl,vdst_addr+=(vl>>3),avl-=vl;
    }

    return 0;
}