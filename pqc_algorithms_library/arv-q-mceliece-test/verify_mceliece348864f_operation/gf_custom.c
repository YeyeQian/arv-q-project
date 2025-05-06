#include "gf.h"
#include "params.h"
#include <stdio.h>

void GF_mul_custom(gf *out, gf *in0, gf *in1)   //VLEN Must be greater equal than 256bits
{
	int i, j;
    size_t vl;
    vuint16m4_t va,vb;
    csr_primpoly_rw(MC_GF_POLY);
	gf prod[ SYS_T*2-1 ]={0};

	for (i = 0; i < SYS_T; i++){
        vl=vsetvl_e16m4(SYS_T);
        va=vle16_v_u16m4(in1,vl);
        vb=vle16_v_u16m4(&prod[i],vl);
        va=vgfmul_vx_u16m4(va,in0[i]);
        vb=vxor_vv_u16m4(vb,va,vl);
        vse16_v_u16m4(&prod[i],vb,vl);
    }

    vl=vsetvl_e16m4(SYS_T-1);
    va=vle16_v_u16m4(prod+SYS_T,vl);
    vb=vle16_v_u16m4(prod+3,vl);
    vb=vxor_vv_u16m4(vb,va,vl);
    vse16_v_u16m4(prod+3,vb,vl);
    vb=vle16_v_u16m4(prod+1,vl);
    vb=vxor_vv_u16m4(vb,va,vl);
    vse16_v_u16m4(prod+1,vb,vl);
    vb=vle16_v_u16m4(prod,vl);
    va=vgfmul_vx_u16m4(va,2);
    vb=vxor_vv_u16m4(vb,va,vl);
    vse16_v_u16m4(prod,vb,vl);

	memcpy(out,prod,SYS_T*sizeof(gf));
}