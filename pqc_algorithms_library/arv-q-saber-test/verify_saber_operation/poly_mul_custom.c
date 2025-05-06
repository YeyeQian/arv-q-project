#include "poly_mul.h"
#include <stdint.h>
#include <string.h>
#include "ntt.h"
#include "reduce.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

void poly_mul_acc_bf_custom(uint32_t a[SABER_N], uint32_t b[SABER_N], uint32_t res[SABER_N]) {
    size_t vl;
    size_t avl=SABER_N;
    uint32_t* res_ptr=res;
    uint32_t* a_ptr=a;
    uint32_t* b_ptr=b;
    vuint32m8_t va;
    vuint32m8_t vb;
    vuint32m8_t vc;
    while(avl>0){
        vl=vsetvl_e32m8(avl);
        va=vle32_v_u32m8(a_ptr,vl);
        vb=vle32_v_u32m8(b_ptr,vl);
        vc=vle32_v_u32m8(res_ptr,vl);
        va=vmod_mul_vv_u32m8(va,vb);
        vc=vmod_add_vv_u32m8(va,vc);
        vse32_v_u32m8(res_ptr,vc,vl);
        a_ptr+=vl;
        b_ptr+=vl;
        res_ptr+=vl;
        avl-=vl;
    }
}