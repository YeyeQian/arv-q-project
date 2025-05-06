#include "params.h"
#include "gf.h"
#include "bm.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include <stdio.h>

#define min(a, b) ((a < b) ? a : b)

void bm_custom(gf *out, gf *s) //VLEN Must be greater equal than 256bit
{
	int i;
    size_t vl,avl;
    vuint16m8_t va,vb,vidx;
    vuint16m1_t vc;
    csr_primpoly_rw(MC_GF_POLY);

	uint16_t N = 0;
	uint16_t L = 0;
	uint16_t mle;
	uint16_t mne;

	gf T[ SYS_T+1  ];
	gf C[ SYS_T+1 ]={0};
	gf B[ SYS_T+1 ]={0};

	gf b = 1, d, f;

	B[1] = C[0] = 1;

	for (N = 0; N < 2 * SYS_T; N++)
	{
		d = 0;

        vl=vsetvl_e16m1(1);
        vc=vmv_s_x_u16m1(vc,0,vl);
		vl=vsetvl_e16m8(min(N+1,SYS_T+1));
        va=vle16_v_u16m8(C,vl);
        vidx=vid_v_u16m8(vl);
        vidx=vrsub_vx_u16m8(vidx,N,vl);
		vidx=vsll_vx_u16m8(vidx,1,vl);
        vb=vluxei16_v_u16m8(s,vidx,vl);
        va=vgfmul_vv_u16m8(va,vb);
        vc=vredxor_vs_u16m8_u16m1(vc,va,vc,vl);
        d=vmv_x_s_u16m1_u16(vc);
		
	
		mne = d; mne -= 1;   mne >>= 15; mne -= 1;
		mle = N; mle -= 2*L; mle >>= 15; mle -= 1;
		mle &= mne;

		memcpy(T,C,(SYS_T+1)<<1);

		f=gfinv_u16(b);
        vl=vsetvl_e16m1(1);
        vc=vgfmul_vx_u16m1(vc,f);
        f=vmv_x_s_u16m1_u16(vc);


		vl=vsetvl_e16m8(1+SYS_T);
        va=vle16_v_u16m8(B,vl);
        vb=vle16_v_u16m8(C,vl);
        va=vgfmul_vx_u16m8(va,f);
        va=vand_vx_u16m8(va,mne,vl);
        vb=vxor_vv_u16m8(vb,va,vl);
        vse16_v_u16m8(C,vb,vl);

		L = (L & ~mle) | ((N+1-L) & mle);

		va=vle16_v_u16m8(B,vl);
        vb=vle16_v_u16m8(T,vl);
        va=vand_vx_u16m8(va,~mle,vl);
        vb=vand_vx_u16m8(vb,mle,vl);
        va=vor_vv_u16m8(va,vb,vl);
        vse16_v_u16m8(B,va,vl);

		b = (b & ~mle) | (d & mle);

		for (i = SYS_T; i >= 1; i--) B[i] = B[i-1];
		B[0] = 0;
	}

	for (i = 0; i <= SYS_T; i++)
		out[i] = C[ SYS_T-i ];
}