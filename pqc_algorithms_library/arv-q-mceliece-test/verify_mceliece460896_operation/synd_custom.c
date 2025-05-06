#include "synd.h"

#include "params.h"
#include "root.h"

#include <stdio.h>
#include <string.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

// void synd_custom(gf *out, gf *f, gf *L, unsigned char *r)
// {
// 	int i, j;
// 	gf e, e_inv, c;
//     size_t vl;
//     vuint16m1_t va,vb,vc;
//     csr_primpoly_rw(MC_GF_POLY);
// 	memset(out,0,SYS_T<<2);

//     vl=vsetvl_e16m1(1);
// 	for (i = 0; i < SYS_N; i++)
// 	{
// 		c = (r[i/8] >> (i%8)) & 1;

//         va=vmv_s_x_u16m1(va,f[SYS_T],vl);
//         for(int j=SYS_T-1;j>=0;j--){
//             va=vgfmul_vx_u16m1(va,L[i]);
//             va=vxor_vx_u16m1(va,f[j],vl);
//         }
        
//         va=vgfmul_vv_u16m1(va,va);
//         e_inv=vmv_x_s_u16m1_u16(va);
//         e_inv=gfinv_u16(e_inv);
//         va=vmv_s_x_u16m1(va,e_inv,vl);

// 		for (j = 0; j < 2*SYS_T; j++)
// 		{
// 			vb=vgfmul_vx_u16m1(va,c);
//             vb=vxor_vx_u16m1(vb,out[j],vl);
//             out[j]=vmv_x_s_u16m1_u16(vb);
//             va=vgfmul_vx_u16m1(va,L[i]);
// 		}
// 	}
// }

void synd_custom(gf *out, gf *f, gf *L, unsigned char *r)
{
	int i, j;
	gf e, e_inv, c;
    size_t vl,avl;
    vuint16m4_t va,vb,vc,vf;
    vuint8m2_t vd,vidx;
    vuint16m1_t ve;
    csr_primpoly_rw(MC_GF_POLY);
	memset(out,0,SYS_T<<2);

    gf* src_addr=L;

    gf temp_buf[256]={0};

    avl=SYS_N;
    while(avl>0){
        vl=vsetvl_e8m2(avl);
        vidx=vid_v_u8m2(vl);
        vidx=vadd_vx_u8m2(vidx,SYS_N-avl,vl);   // i
        vidx=vsrl_vx_u8m2(vidx,3,vl);           // i/8
        vd=vluxei8_v_u8m2(r,vidx,vl);           // r[i/8]
        vidx=vand_vx_u8m2(vidx,7,vl);           // i%8
        vd=vsrl_vv_u8m2(vd,vidx,vl);            // r[i/8] >> (i%8)
        vd=vand_vx_u8m2(vd,1,vl);
        vc=vwaddu_vx_u16m4(vd,0,vl);            // c
        vl=vsetvl_e16m4(avl);

        va=vmv_v_x_u16m4(f[SYS_T],vl);
        vb=vle16_v_u16m4(src_addr,vl);          // L[i]
        for(int i=SYS_T-1;i>=0;i--){
            va=vgfmul_vv_u16m4(va,vb);
            va=vxor_vx_u16m4(va,f[i],vl);
        }                                       // va stores e
        va=vgfmul_vv_u16m4(va,va);              // gfmul(e,e)
        vse16_v_u16m4(temp_buf,va,vl);
        for(int i=0;i<vl;i++) temp_buf[i]=gfinv_u16(temp_buf[i]);
        va=vle16_v_u16m4(temp_buf,vl);          // einv

        for(j=0;j<2*SYS_T;j++){
            vl=vsetvl_e16m1(1);
            ve=vmv_s_x_u16m1(ve,out[j],vl);
            vl=vsetvl_e16m4(avl);
            vf=vgfmul_vv_u16m4(va,vc);          // gf_mul(e_inv, c)
            ve=vredxor_vs_u16m4_u16m1(ve,vf,ve,vl);
            out[j]=vmv_x_s_u16m1_u16(ve);
            va=vgfmul_vv_u16m4(va,vb);
        }
        src_addr+=vl,avl-=vl;
    }
}