#include <stdint.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "controlbits.h"
#include "uint64_sort.h"
#include "pk_gen.h"
#include "params.h"
#include "benes.h"
#include "root.h"
#include "util.h"
#include "crypto_declassify.h"
#include "crypto_uint64.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

static crypto_uint64 uint64_is_equal_declassify(uint64_t t,uint64_t u)
{
  crypto_uint64 mask = crypto_uint64_equal_mask(t,u);
  crypto_declassify(&mask,sizeof mask);
  return mask;
}

static crypto_uint64 uint64_is_zero_declassify(uint64_t t)
{
  crypto_uint64 mask = crypto_uint64_zero_mask(t);
  crypto_declassify(&mask,sizeof mask);
  return mask;
}

int pk_gen_custom(unsigned char * pk, unsigned char * sk, uint32_t * perm, int16_t * pi)    //VLEN Must be greater equal then 256bit
{
	int i, j, k;
	int row, c;
    size_t vl,avl;
    vuint16m2_t va,vb;
    csr_primpoly_rw(MC_GF_POLY);
	uint64_t buf[ 1 << GFBITS ];

	unsigned char mat[ PK_NROWS ][ SYS_N/8 ]={0};
	unsigned char mask;
	unsigned char b;

	gf g[ SYS_T+1 ]; // Goppa polynomial
	gf L[ SYS_N ]; // support
	gf inv[ SYS_N ];

	//

	g[ SYS_T ] = 1;

    vuint16m4_t vd;
	vl=vsetvl_e16m4(SYS_T);
    vd=vle16_v_u16m4((uint16_t*)sk,vl);
    vd=vand_vx_u16m4(vd,GFMASK,vl);
    vse16_v_u16m4(g,vd,vl);

	avl=1<<GFBITS;
    uint32_t* perm_addr=perm;
    uint64_t* buf_addr=buf;
    vuint32m4_t vc;
    vuint64m8_t vbuf,vidx;
    while(avl>0){
        vl=vsetvl_e32m4(avl);
        vc=vle32_v_u32m4(perm_addr,vl);
        vbuf=vwaddu_vx_u64m8(vc,0,vl);
        vl=vsetvl_e64m8(avl);
        vidx=vid_v_u64m8(vl);
        vidx=vadd_vx_u64m8(vidx,(1<<GFBITS)-avl,vl);
        vbuf=vsll_vx_u64m8(vbuf,31,vl);
        vbuf=vor_vv_u64m8(vbuf,vidx,vl);
        vse64_v_u64m8(buf_addr,vbuf,vl);
        perm_addr+=vl,buf_addr+=vl,avl-=vl;
    }

	uint64_sort(buf, 1 << GFBITS);

	for (i = 1; i < (1 << GFBITS); i++)
		if (uint64_is_equal_declassify(buf[i-1] >> 31,buf[i] >> 31))
			return -1;

	avl=SYS_N;
    buf_addr=buf;
    uint16_t* pi_addr=NULL;
    pi_addr=(uint16_t*) pi;
    uint16_t* L_addr=L;
    while(avl>0){
        vl=vsetvl_e64m8(avl);
        vbuf=vle64_v_u64m8(buf_addr,vl);
        vbuf=vand_vx_u64m8(vbuf,GFMASK,vl);
        vc=vnsrl_wx_u32m4(vbuf,0,vl);
        vl=vsetvl_e32m4(avl);
        va=vnsrl_wx_u16m2(vc,0,vl);
        vl=vsetvl_e16m2(avl);
        vse16_v_u16m2(pi_addr,va,vl);
        vb=vmv_v_v_u16m2(va,vl);

        va=vand_vx_u16m2(va,0x00FF,vl);
        va=vsll_vx_u16m2(va,8,vl);
        vb=vand_vx_u16m2(vb,0xFF00,vl);
        vb=vsrl_vx_u16m2(vb,8,vl);
        va=vor_vv_u16m2(va,vb,vl);
        vb=vmv_v_v_u16m2(va,vl);

        va=vand_vx_u16m2(va,0x0F0F,vl);
        va=vsll_vx_u16m2(va,4,vl);
        vb=vand_vx_u16m2(vb,0xF0F0,vl);
        vb=vsrl_vx_u16m2(vb,4,vl);
        va=vor_vv_u16m2(va,vb,vl);
        vb=vmv_v_v_u16m2(va,vl);

        va=vand_vx_u16m2(va,0x3333,vl);
        va=vsll_vx_u16m2(va,2,vl);
        vb=vand_vx_u16m2(vb,0xCCCC,vl);
        vb=vsrl_vx_u16m2(vb,2,vl);
        va=vor_vv_u16m2(va,vb,vl);
        vb=vmv_v_v_u16m2(va,vl);

        va=vand_vx_u16m2(va,0x5555,vl);
        va=vsll_vx_u16m2(va,1,vl);
        vb=vand_vx_u16m2(vb,0xAAAA,vl);
        vb=vsrl_vx_u16m2(vb,1,vl);
        va=vor_vv_u16m2(va,vb,vl);
        
        va=vsrl_vx_u16m2(va,4,vl);
        vse16_v_u16m2(L_addr,va,vl);

        buf_addr+=vl,L_addr+=vl,pi_addr+=vl,avl-=vl;
    }
    for(int i=SYS_N;i<(1 << GFBITS);i++)pi[i] = buf[i] & GFMASK;

	// filling the matrix

	root_custom(inv, g, L);
		
	for (i = 0; i < SYS_N; i++)
		inv[i] = gfinv_u16(inv[i]);


	for (i = 0; i < SYS_T; i++)
	{
		for (j = 0; j < SYS_N; j+=8)
		for (k = 0; k < GFBITS;  k++)
		{
			b  = (inv[j+7] >> k) & 1; b <<= 1;
			b |= (inv[j+6] >> k) & 1; b <<= 1;
			b |= (inv[j+5] >> k) & 1; b <<= 1;
			b |= (inv[j+4] >> k) & 1; b <<= 1;
			b |= (inv[j+3] >> k) & 1; b <<= 1;
			b |= (inv[j+2] >> k) & 1; b <<= 1;
			b |= (inv[j+1] >> k) & 1; b <<= 1;
			b |= (inv[j+0] >> k) & 1;

			mat[ i*GFBITS + k ][ j/8 ] = b;
		}

		avl=SYS_N;
        L_addr=L;
        gf* inv_addr=inv;
        while(avl>0){
            vl=vsetvl_e16m2(avl);
            va=vle16_v_u16m2(inv_addr,vl);
            vb=vle16_v_u16m2(L_addr,vl);
            va=vgfmul_vv_u16m2(va,vb);
            vse16_v_u16m2(inv_addr,va,vl);
            inv_addr+=vl,L_addr+=vl,avl-=vl;
        }
	}

	// gaussian elimination
    vuint8m4_t ve,vf;
	for (i = 0; i < (PK_NROWS + 7) / 8; i++)
	for (j = 0; j < 8; j++)
	{
		row = i*8 + j;			

		if (row >= PK_NROWS)
			break;

		for (k = row + 1; k < PK_NROWS; k++)
		{
			mask = mat[ row ][ i ] ^ mat[ k ][ i ];
			mask >>= j;
			mask &= 1;
			mask = -mask;

            avl=SYS_N>>3;
            uint8_t* src_addr=&mat[k][0];
            uint8_t* dst_addr=&mat[row][0];
			while(avl>0){
                vl=vsetvl_e8m4(avl);
                ve=vle8_v_u8m4(src_addr,vl);
                vf=vle8_v_u8m4(dst_addr,vl);
                ve=vand_vx_u8m4(ve,mask,vl);
                vf=vxor_vv_u8m4(vf,ve,vl);
                vse8_v_u8m4(dst_addr,vf,vl);
                dst_addr+=vl,src_addr+=vl,avl-=vl;
            }
		}

                if ( uint64_is_zero_declassify((mat[ row ][ i ] >> j) & 1) ) // return if not systematic
		{
			return -1;
		}

		for (k = 0; k < PK_NROWS; k++)
		{
			if (k != row)
			{
				mask = mat[ k ][ i ] >> j;
				mask &= 1;
				mask = -mask;

				avl=SYS_N>>3;
                uint8_t* src_addr=&mat[row][0];
                uint8_t* dst_addr=&mat[k][0];
                while(avl>0){
                    vl=vsetvl_e8m4(avl);
                    ve=vle8_v_u8m4(src_addr,vl);
                    vf=vle8_v_u8m4(dst_addr,vl);
                    ve=vand_vx_u8m4(ve,mask,vl);
                    vf=vxor_vv_u8m4(vf,ve,vl);
                    vse8_v_u8m4(dst_addr,vf,vl);
                    dst_addr+=vl,src_addr+=vl,avl-=vl;
                }
			}
		}
	}

	for (i = 0; i < PK_NROWS; i++)
		memcpy(pk + i*PK_ROW_BYTES, mat[i] + PK_NROWS/8, PK_ROW_BYTES);

	return 0;
}