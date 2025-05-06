#include "sampling.h"
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"

void setZero_custom(uint8_t * r, uint32_t length)
{
    int32_t avl=length;
    size_t vl;
    uint8_t* vdst_addr=r;

    while(avl>0){
        vl=vsetvl_e8m8(avl);
        vuint8m8_t vres=vmv_v_x_u8m8(0,vl);
        vse8_v_u8m8(vdst_addr,vres,vl);

        vdst_addr+=vl,avl-=vl;
    }
}