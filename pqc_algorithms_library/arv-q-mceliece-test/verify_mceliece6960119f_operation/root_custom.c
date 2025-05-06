#include "root.h"
#include "params.h"
#include "gf.h"

#include <stdio.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

void root_custom(gf *out, gf *f, gf *L)
{
	int i; 
    size_t vl,avl;
    vuint16m4_t va,vb;
    gf* dst_addr=out;
    gf* src_addr=L;
    avl=SYS_N;
    csr_primpoly_rw(MC_GF_POLY);
    while(avl>0){
        vl=vsetvl_e16m4(avl);
        va=vmv_v_x_u16m4(f[SYS_T],vl);
        vb=vle16_v_u16m4(src_addr,vl);
        for(int i=SYS_T-1;i>=0;i--){
            va=vgfmul_vv_u16m4(va,vb);
            va=vxor_vx_u16m4(va,f[i],vl);
        }
        vse16_v_u16m4(dst_addr,va,vl);
        dst_addr+=vl,src_addr+=vl,avl-=vl;
    }
	
}