#include "transpose.h"

#include <stdint.h>
#include <string.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

/* input: in, a 64x64 matrix over GF(2) */
/* output: out, transpose of in */
void transpose_64x64_custom(uint64_t * out, uint64_t * in)
{
	int i, j, s, d;
    size_t vl,avl;
    vuint64m4_t va,vb,vc,vd;

	uint64_t x, y;
	uint64_t masks[6][2] = {
	                        {0x5555555555555555, 0xAAAAAAAAAAAAAAAA},
	                        {0x3333333333333333, 0xCCCCCCCCCCCCCCCC},
	                        {0x0F0F0F0F0F0F0F0F, 0xF0F0F0F0F0F0F0F0},
	                        {0x00FF00FF00FF00FF, 0xFF00FF00FF00FF00},
	                        {0x0000FFFF0000FFFF, 0xFFFF0000FFFF0000},
	                        {0x00000000FFFFFFFF, 0xFFFFFFFF00000000}
	                       };

	memcpy(out,in,64<<3);

	for (d = 5; d >= 0; d--)
	{
		s = 1 << d;

		for (i = 0; i < 64; i += s*2){
            avl=s;
            uint64_t* src1_addr=&out[i];
            uint64_t* src2_addr=&out[i+s];
            while(avl>0){
                vl=vsetvl_e64m4(avl);
                va=vle64_v_u64m4(src1_addr,vl);
                vb=vle64_v_u64m4(src2_addr,vl);
                vc=vand_vx_u64m4(va,masks[d][0],vl);
                vd=vand_vx_u64m4(vb,masks[d][0],vl);
                vd=vsll_vx_u64m4(vd,s,vl);
                vd=vor_vv_u64m4(vd,vc,vl);
                vse64_v_u64m4(src1_addr,vd,vl);
                vc=vand_vx_u64m4(va,masks[d][1],vl);
                vc=vsrl_vx_u64m4(vc,s,vl);
                vd=vand_vx_u64m4(vb,masks[d][1],vl);
                vd=vor_vv_u64m4(vd,vc,vl);
                vse64_v_u64m4(src2_addr,vd,vl);
                src1_addr+=vl,src2_addr+=vl,avl-=vl;
            }
        }
	}
}