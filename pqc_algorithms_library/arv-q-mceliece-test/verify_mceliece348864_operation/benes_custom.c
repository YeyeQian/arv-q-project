#include "util.h"
#include "transpose.h"
#include "params.h"
#include "benes.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

static void layer_custom(uint64_t * data, uint64_t * bits, int lgs) 
{
	int i, j, s;
    size_t vl,avl;
    vuint64m4_t va,vb,vc,vd;

	uint64_t d;

	s = 1 << lgs;
    uint64_t* src3_addr=bits;
	for (i = 0; i < 64; i += s<<1){
        avl=s;
        uint64_t *src1_addr, *src2_addr;
        src1_addr=&data[i];
        src2_addr=&data[i+s];
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            va=vle64_v_u64m4(src1_addr,vl);
            vb=vle64_v_u64m4(src2_addr,vl);
            vc=vxor_vv_u64m4(va,vb,vl);
            vd=vle64_v_u64m4(src3_addr,vl);
            vd=vand_vv_u64m4(vd,vc,vl);
            va=vxor_vv_u64m4(va,vd,vl);
            vb=vxor_vv_u64m4(vb,vd,vl);
            vse64_v_u64m4(src1_addr,va,vl);
            vse64_v_u64m4(src2_addr,vb,vl);
            src1_addr+=vl,src2_addr+=vl,src3_addr+=vl,avl-=vl;
        }
    }
}

void apply_benes_custom(unsigned char * r, const unsigned char * bits, int rev)
{
	int i;

	const unsigned char *cond_ptr; 
	int inc, low;

	uint64_t bs[64];
	uint64_t cond[64];

	//

	for (i = 0; i < 64; i++)
	{
		bs[i] = load8(r + (i<<3));
	}

	if (rev == 0) 
	{
		inc = 256;
		cond_ptr = bits;
	}
	else
	{
		inc = -256;
		cond_ptr = bits + ((2*GFBITS-2)<<8);
	}

	//

	transpose_64x64_custom(bs, bs);

	for (low = 0; low <= 5; low++) 
	{ 
		for (i = 0; i < 64; i++) cond[i] = load4(cond_ptr + (i<<2));
		transpose_64x64_custom(cond, cond);
		layer_custom(bs, cond, low); 
		cond_ptr += inc; 
	}
	
	transpose_64x64_custom(bs, bs);
	
	for (low = 0; low <= 5; low++) 
	{ 
		for (i = 0; i < 32; i++) cond[i] = load8(cond_ptr + (i<<3));
		layer_custom(bs, cond, low); 
		cond_ptr += inc; 
	}
	for (low = 4; low >= 0; low--) 
	{ 
		for (i = 0; i < 32; i++) cond[i] = load8(cond_ptr + (i<<3));
		layer_custom(bs, cond, low); 
		cond_ptr += inc; 
	}

	transpose_64x64_custom(bs, bs);
	
	for (low = 5; low >= 0; low--) 
	{ 
		for (i = 0; i < 64; i++) cond[i] = load4(cond_ptr + (i<<2));
		transpose_64x64_custom(cond, cond);
		layer_custom(bs, cond, low); 
		cond_ptr += inc; 
	}

	transpose_64x64_custom(bs, bs);

	//

	for (i = 0; i < 64; i++)
	{
		store8(r + (i<<3), bs[i]);
	}
}

void support_gen_custom(gf * s, const unsigned char *c)
{
	gf a;
	int i, j;
	unsigned char L[ GFBITS ][ (1 << GFBITS)>>3 ]={0};

	for (i = 0; i < (1 << GFBITS); i++)
	{
		a = bitrev((gf) i);

		for (j = 0; j < GFBITS; j++)
			L[j][ i>>3 ] |= ((a >> j) & 1) << (i&7);
	}
			
	for (j = 0; j < GFBITS; j++)
		apply_benes_custom(L[j], c, 0);

	for (i = 0; i < SYS_N; i++)
	{
		s[i] = 0;
		for (j = GFBITS-1; j >= 0; j--)
		{
			s[i] <<= 1;
			s[i] |= (L[j][i>>3] >> (i&7)) & 1;
		}
	}
}