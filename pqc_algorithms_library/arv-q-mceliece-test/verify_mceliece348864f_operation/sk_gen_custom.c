#include "sk_gen.h"

#include "controlbits.h"
#include "params.h"
#include "util.h"
#include "gf.h"
#include "crypto_declassify.h"
#include "crypto_uint16.h"
#include <string.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

static inline crypto_uint16 gf_is_zero_declassify(gf t)
{
  crypto_uint16 mask = crypto_uint16_zero_mask(t);
  crypto_declassify(&mask,sizeof mask);
  return mask;
}

int genpoly_gen_custom(gf *out, gf *f)  //VLEN Must be greater equal than 256
{
	int i, j, k, c;
    size_t vl,avl;
    vuint16m4_t va,vb,vidx;
    vbool4_t vflag;
    vuint16m1_t vc;
    csr_primpoly_rw(MC_GF_POLY);

	gf mat[ SYS_T+1 ][ SYS_T ];
	gf mask, inv, t;

	// fill matrix

	mat[0][0] = 1;

	for (i = 1; i < SYS_T; i++)
		mat[0][i] = 0;

	for (i = 0; i < SYS_T; i++)
		mat[1][i] = f[i];

	for (j = 2; j <= SYS_T; j++)
		GF_mul_custom(mat[j], mat[j-1], f);

	// gaussian //Gauss-Jordan Elimination (Column Major)

	for (j = 0; j < SYS_T; j++)
	{
        // mask = gf_iszero(mat[ j ][ j ]);
		// for (c = j; c < SYS_T + 1; c++)
		// {   
        //     vl=vsetvl_e16m1(1);
        //     vc=vmv_s_x_u16m1(vc,mat[c][j],vl);
		// 	avl=SYS_T-j-1;
        //     vl=vsetvl_e16m4(avl);
        //     va=vle16_v_u16m4(&mat[c][j+1],vl);
        //     va=vand_vx_u16m4(va,mask,vl);
        //     vc=vredxor_vs_u16m4_u16m1(vc,va,vc,vl);
        //     mat[c][j]=vmv_x_s_u16m1_u16(vc);
		// }
        for (k = j + 1; k < SYS_T; k++)
		{
			mask = gf_iszero(mat[ j ][ j ]);

			for (c = j; c < SYS_T + 1; c++)
				mat[ c ][ j ] ^= mat[ c ][ k ] & mask;

		}

		if ( gf_is_zero_declassify(mat[ j ][ j ]) ) // return if not systematic
		{
			return -1;
		}

		inv = gfinv_u16(mat[j][j]);

        avl=SYS_T+1-j;
        vl=vsetvl_e16m4(avl);
        va=vlse16_v_u16m4(&mat[j][j],SYS_T<<1,vl);
        va=vgfmul_vx_u16m4(va,inv);
        vsse16_v_u16m4(&mat[j][j],SYS_T<<1,va,vl);

		for (c = j; c < SYS_T + 1; c++) {
            vl=vsetvl_e16m4(SYS_T);
            va=vle16_v_u16m4(&mat[c][0],vl);
            vb=vle16_v_u16m4(&mat[j][0],vl);
            vidx=vid_v_u16m4(vl);
            vflag=vmsne_vx_u16m4_b4(vidx,j,vl);
            vb=vgfmul_vx_u16m4(vb,mat[c][j]);
            va=vxor_vv_u16m4_m(vflag,va,va,vb,vl);
            vse16_v_u16m4(&mat[c][0],va,vl);
		}
	}

	memcpy(out,&mat[SYS_T][0],SYS_T*sizeof(gf));

	return 0;
}