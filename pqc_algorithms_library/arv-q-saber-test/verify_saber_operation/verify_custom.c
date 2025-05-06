#include <string.h>
#include <stdint.h>
#include "verify.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

/* returns 0 for equal strings, 1 for non-equal strings */
int verify_custom(const unsigned char *a, const unsigned char *b, size_t len)
{
    uint64_t r;
    size_t i;
    r = 0;
    size_t vl;
    size_t avl=len;
    
    vuint8m8_t va,vb;
    uint8_t* a_addr=a;
    uint8_t* b_addr=b;

    uint8_t r_tmp=0;
    vuint8m1_t vr;
    vl=vsetvl_e8m1(1);
    vr=vmv_v_x_u8m1(r_tmp,vl);

    while(avl>0){
        vl=vsetvl_e8m8(avl);
        va=vle8_v_u8m8(a_addr,vl);
        vb=vle8_v_u8m8(b_addr,vl);
        va=vxor_vv_u8m8(va,vb,vl);
        vr=vredor_vs_u8m8_u8m1(vr,va,vr,vl);
        a_addr+=vl,b_addr+=vl,avl-=vl;
    }

    r_tmp=vmv_x_s_u8m1_u8(vr);
    r=(uint64_t)r_tmp;

    r = (-r) >> 63;
    return r;
}

/* b = 1 means mov, b = 0 means don't mov*/
void cmov_custom(unsigned char *r, const unsigned char *x, size_t len, unsigned char b)
{
    size_t vl;
    size_t avl=len;
    vuint8m8_t vr,vx;
    uint8_t* r_addr=r;
    uint8_t* x_addr=x;

    b = -b;

    while(avl>0){
        vl=vsetvl_e8m8(avl);
        vr=vle8_v_u8m8(r_addr,vl);
        vx=vle8_v_u8m8(x_addr,vl);
        vx=vxor_vv_u8m8(vx,vr,vl);
        vx=vand_vx_u8m8(vx,b,vl);
        vr=vxor_vv_u8m8(vr,vx,vl);
        vse8_v_u8m8(r_addr,vr,vl);
        r_addr+=vl,x_addr+=vl,avl-=vl;
    }
}