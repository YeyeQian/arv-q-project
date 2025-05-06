#include "gf.h"
#include "params.h"
#include <stdio.h>

void GF_mul_custom(gf *out, gf *in0, gf *in1)   //VLEN Must be greater equal than 512bits
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

    vb=vle16_v_u16m4(prod+10,vl);
    vb=vxor_vv_u16m4(vb,va,vl);
    vse16_v_u16m4(prod+10,vb,vl);

    vb=vle16_v_u16m4(prod+9,vl);
    vb=vxor_vv_u16m4(vb,va,vl);
    vse16_v_u16m4(prod+9,vb,vl);

    vb=vle16_v_u16m4(prod+6,vl);
    vb=vxor_vv_u16m4(vb,va,vl);
    vse16_v_u16m4(prod+6,vb,vl);

    vb=vle16_v_u16m4(prod+0,vl);
    vb=vxor_vv_u16m4(vb,va,vl);
    vse16_v_u16m4(prod+0,vb,vl);

	memcpy(out,prod,SYS_T*sizeof(gf));
}

gf gf_inv_custom(gf den)
{
    gf num=1;
	gf tmp_11;
	gf tmp_1111;
	gf out;

    tmp_11=gfmul_u16(den,den);
    tmp_11=gfmul_u16(tmp_11,den);//^3
    tmp_1111=gfmul_u16(tmp_11,tmp_11);
    tmp_1111=gfmul_u16(tmp_1111,tmp_1111);
    tmp_1111=gfmul_u16(tmp_1111,tmp_11);//^15
    out=gfmul_u16(tmp_1111,tmp_1111);//^30
    out=gfmul_u16(out,out);//^60
    out=gfmul_u16(out,out);//^120
    out=gfmul_u16(out,out);//^240
    out=gfmul_u16(out,tmp_1111);//^255

    out=gfmul_u16(out,out);//^510
    out=gfmul_u16(out,out);//^1020
    out=gfmul_u16(out,out);//^2040
    out=gfmul_u16(out,out);//^4080
    out=gfmul_u16(out,tmp_1111);//^4095

    out=gfmul_u16(out,out);//^8190
    //out=gfmul_u16(out,1);

    return out;
}