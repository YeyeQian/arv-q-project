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

gf gf_inv_custom(gf in)
{
    size_t vl;
    size_t avl=1;
    vuint16m1_t vtmp_11,vtmp_1111,vout;
    csr_primpoly_rw(MC_GF_POLY);

    vl=vsetvl_e16m1(avl);
    vout=vmv_s_x_u16m1(vout,in,vl);

	vout=vgfmul_vv_u16m1(vout,vout);//out = gf_sq(out);
	vtmp_11=vgfmul_vx_u16m1(vout,in);//tmp_11 = gf_mul(out, in); // 11

	vout=vgfmul_vv_u16m1(vtmp_11,vtmp_11);//out = gf_sq(tmp_11);
	vout=vgfmul_vv_u16m1(vout,vout);//out = gf_sq(out);
	vtmp_1111=vgfmul_vv_u16m1(vout,vtmp_11);//tmp_1111 = gf_mul(out, tmp_11); // 1111

	vout=vgfmul_vv_u16m1(vtmp_1111,vtmp_1111);//out = gf_sq(tmp_1111);
	vout=vgfmul_vv_u16m1(vout,vout);//out = gf_sq(out);
	vout=vgfmul_vv_u16m1(vout,vout);//out = gf_sq(out);
	vout=vgfmul_vv_u16m1(vout,vout);//out = gf_sq(out);
	vout=vgfmul_vv_u16m1(vout,vtmp_1111);//out = gf_mul(out, tmp_1111); // 11111111

	vout=vgfmul_vv_u16m1(vout,vout);//out = gf_sq(out);
	vout=vgfmul_vv_u16m1(vout,vout);//out = gf_sq(out);
	vout=vgfmul_vv_u16m1(vout,vtmp_11);//out = gf_mul(out, tmp_11); // 1111111111

	vout=vgfmul_vv_u16m1(vout,vout);//out = gf_sq(out);
	vout=vgfmul_vx_u16m1(vout,in);//out = gf_mul(out, in); // 11111111111
    vout=vgfmul_vv_u16m1(vout,vout);

    uint16_t res;
    res=vmv_x_s_u16m1_u16(vout);

	return res;
}