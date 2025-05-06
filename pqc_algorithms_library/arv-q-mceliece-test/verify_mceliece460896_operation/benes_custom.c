#include "util.h"
#include "transpose.h"
#include "params.h"
#include "benes.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

static void layer_in_custom(uint64_t data[2][64], uint64_t * bits, int lgs)
{
	int i, j, s;
    size_t vl,avl;
    vuint64m4_t va,vb,vc,vd;

	uint64_t d;
    uint64_t* src3_addr=bits;
	s = 1 << lgs;
	for (i = 0; i < 64; i += s<<1){
        avl=s;
        uint64_t *src1_addr, *src2_addr;
        uint64_t* src3_addr_tmp=src3_addr;
        src1_addr=&data[0][i];
        src2_addr=&data[0][i+s];
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            va=vle64_v_u64m4(src1_addr,vl);
            vb=vle64_v_u64m4(src2_addr,vl);
            vc=vxor_vv_u64m4(va,vb,vl);
            vd=vlse64_v_u64m4(src3_addr_tmp,8,vl);
            vd=vand_vv_u64m4(vd,vc,vl);
            va=vxor_vv_u64m4(va,vd,vl);
            vb=vxor_vv_u64m4(vb,vd,vl);
            vse64_v_u64m4(src1_addr,va,vl);
            vse64_v_u64m4(src2_addr,vb,vl);
            src1_addr+=vl,src2_addr+=vl,src3_addr_tmp+=(vl<<1),avl-=vl;
        }

        avl=s;
        src3_addr_tmp=src3_addr+1;
        src1_addr=&data[1][i];
        src2_addr=&data[1][i+s];
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            va=vle64_v_u64m4(src1_addr,vl);
            vb=vle64_v_u64m4(src2_addr,vl);
            vc=vxor_vv_u64m4(va,vb,vl);
            vd=vlse64_v_u64m4(src3_addr_tmp,8,vl);
            vd=vand_vv_u64m4(vd,vc,vl);
            va=vxor_vv_u64m4(va,vd,vl);
            vb=vxor_vv_u64m4(vb,vd,vl);
            vse64_v_u64m4(src1_addr,va,vl);
            vse64_v_u64m4(src2_addr,vb,vl);
            src1_addr+=vl,src2_addr+=vl,src3_addr_tmp+=(vl<<1),avl-=vl;
        }

        src3_addr+=(s<<1);
    }
}

static void layer_ex_custom(uint64_t * data, uint64_t * bits, int lgs)
{
	int i, j, s;
    size_t vl,avl;
    vuint64m4_t va,vb,vc,vd;

	uint64_t d;
    uint64_t* src3_addr=bits;
	s = 1 << lgs;

	for (i = 0; i < 128; i += s*2){
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
	int i, iter, inc;

	unsigned char *r_ptr = r;
	const unsigned char *bits_ptr;

	uint64_t r_int_v[2][64];
	uint64_t r_int_h[2][64];
	uint64_t b_int_v[64];
	uint64_t b_int_h[64];

	//

	if (rev) { bits_ptr = bits + 12288; inc = -1024; }
	else     { bits_ptr = bits;         inc = 0;    }
		
	for (i = 0; i < 64; i++)
	{
		r_int_v[0][i] = load8(r_ptr + i*16 + 0);
		r_int_v[1][i] = load8(r_ptr + i*16 + 8);
	}

	transpose_64x64_custom(r_int_h[0], r_int_v[0]);
	transpose_64x64_custom(r_int_h[1], r_int_v[1]);

	for (iter = 0; iter <= 6; iter++)
	{
		for (i = 0; i < 64; i++)
		{
			b_int_v[i] = load8(bits_ptr); bits_ptr += 8;
		}

		bits_ptr += inc;

		transpose_64x64_custom(b_int_h, b_int_v);
	
		layer_ex_custom(r_int_h[0], b_int_h, iter);
	}

	transpose_64x64_custom(r_int_v[0], r_int_h[0]);
	transpose_64x64_custom(r_int_v[1], r_int_h[1]);

	for (iter = 0; iter <= 5; iter++) 
	{
		for (i = 0; i < 64; i++) { b_int_v[i] = load8(bits_ptr); bits_ptr += 8; }
		bits_ptr += inc;

		layer_in_custom(r_int_v, b_int_v, iter);
	}

	for (iter = 4; iter >= 0; iter--) 
	{
		for (i = 0; i < 64; i++) { b_int_v[i] = load8(bits_ptr); bits_ptr += 8; }
		bits_ptr += inc;

		layer_in_custom(r_int_v, b_int_v, iter);
	}

	transpose_64x64_custom(r_int_h[0], r_int_v[0]);
	transpose_64x64_custom(r_int_h[1], r_int_v[1]);

	for (iter = 6; iter >= 0; iter--)
	{
		for (i = 0; i < 64; i++)
		{
			b_int_v[i] = load8(bits_ptr); bits_ptr += 8;
		}

		bits_ptr += inc;

		transpose_64x64_custom(b_int_h, b_int_v);

		layer_ex_custom(r_int_h[0], b_int_h, iter);
	}

	transpose_64x64_custom(r_int_v[0], r_int_h[0]);
	transpose_64x64_custom(r_int_v[1], r_int_h[1]);

	for (i = 0; i < 64; i++)
	{
		store8(r_ptr + i*16 + 0, r_int_v[0][i]);
		store8(r_ptr + i*16 + 8, r_int_v[1][i]);
	}
}

void support_gen_custom(gf * s, const unsigned char *c)
{
	gf a;
	int i, j;
	unsigned char L[ GFBITS ][ (1 << GFBITS)/8 ];

	for (i = 0; i < GFBITS; i++)
		for (j = 0; j < (1 << GFBITS)/8; j++)
			L[i][j] = 0;

	for (i = 0; i < (1 << GFBITS); i++)
	{
		a = bitrev((gf) i);

		for (j = 0; j < GFBITS; j++)
			L[j][ i/8 ] |= ((a >> j) & 1) << (i%8);
	}
			
	for (j = 0; j < GFBITS; j++)
		apply_benes_custom(L[j], c, 0);

	for (i = 0; i < SYS_N; i++)
	{
		s[i] = 0;
		for (j = GFBITS-1; j >= 0; j--)
		{
			s[i] <<= 1;
			s[i] |= (L[j][i/8] >> (i%8)) & 1;
		}
	}
}