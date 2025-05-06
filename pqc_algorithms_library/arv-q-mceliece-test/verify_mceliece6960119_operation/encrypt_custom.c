#include "encrypt.h"

#include "util.h"
#include "params.h"

#include <stdint.h>
#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "gf.h"
#include "crypto_declassify.h"
#include "crypto_uint16.h"
#include "crypto_uint32.h"
#include "crypto_hash.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

#include "crypto_kem.h"
#include "crypto_kem_mceliece6960119.h"
#include "pk_sk_cons.h"

static inline crypto_uint16 uint16_is_smaller_declassify(uint16_t t,uint16_t u)
{
  crypto_uint16 mask = crypto_uint16_smaller_mask(t,u);
  crypto_declassify(&mask,sizeof mask);
  return mask;
}

static inline crypto_uint32 uint32_is_equal_declassify(uint32_t t,uint32_t u)
{
  crypto_uint32 mask = crypto_uint32_equal_mask(t,u);
  crypto_declassify(&mask,sizeof mask);
  return mask;
}

static inline unsigned char same_mask(uint16_t x, uint16_t y)
{
	uint32_t mask;

	mask = x ^ y;
	mask -= 1;
	mask >>= 31;
	mask = -mask;

	return mask & 0xFF;
}

void gen_e_custom(unsigned char *e)  //VLEN Must be greater equal than 256bit
{
	int i, j, eq;
    size_t ctr, avl, vl, valid_num;
    vuint16m1_t va;
    csr_validnum_rw();
    ctr = 0;

	union 
	{
		uint16_t nums[ SYS_T*2 ];
		unsigned char bytes[ SYS_T*2 * sizeof(uint16_t) ];
	} buf;

	srand(time(NULL));
	for (int i = 0; i < 2 * SYS_T; i++){
    // buf.nums[i] = rand();
    buf.nums[i]=i;//for test
  }

	uint16_t ind[ SYS_T ];
	unsigned char mask;	
	unsigned char val[ SYS_T ];	

	while (1)
	{
		shake256_custom(buf.bytes,sizeof(buf),buf.bytes,sizeof(buf));

    avl=SYS_T<<1;
    uint16_t* buf_addr=buf.nums;
    while(avl>0){
      vl=vsetvl_e16m1(avl);
      va=vle16_v_u16m1(buf_addr,vl);
      va=vand_vx_u16m1(va,GFMASK,vl);
      vse16_v_u16m1(buf_addr,va,vl);
      buf_addr+=vl;
      avl-=vl;
    }

    avl=SYS_T<<1;
    buf_addr=buf.nums;
    ctr=0;
    valid_num=0;
    while(avl>0 && ctr<SYS_T){
        vl=vsetvl_e16m1(avl);
        va=vle16_v_u16m1(buf_addr,vl);
        va=sample_rej_vx_u16m1(va,SYS_N);
        valid_num=csr_validnum_rw();
        if(ctr+valid_num>SYS_T){
            valid_num=SYS_T-ctr;
        }
        buf_addr+=vl;
        avl-=vl;
        vl=vsetvl_e16m1(valid_num);
        vse16_v_u16m1(ind+ctr,va,valid_num);
        ctr+=valid_num;
    }
		
		if (ctr < SYS_T) continue;

		// check for repetition

		eq = 0;

		for (i = 1; i < SYS_T; i++) 
			for (j = 0; j < i; j++)
			        if (uint32_is_equal_declassify(ind[i],ind[j]))
					eq = 1;

		if (eq == 0)
			break;
	}

	for(int i=0;i<SYS_T;i++){
        uint16_t index=ind[i]>>3;
        uint8_t pos=ind[i]&7;
        e[index]|=1<<pos;
    }
}

void syndrome_custom(unsigned char *s, const unsigned char *pk, unsigned char *e)
{
	unsigned char b;
  unsigned char pk_tmp[PK_ROW_BYTES];
  const unsigned char *pk_ptr=pk;

	int i, j, tail = PK_NROWS & 7;
  size_t vl,avl;
  vuint8m8_t va,vb;
  vuint8m1_t vc;

	memset(s,0,SYND_BYTES);
  uint8_t old_offset;//no useful meaning
  old_offset=csr_bseqoffset_rw(tail);

	for (i = 0; i < PK_NROWS; i++)	
	{
		//shift left pk
    int32_t remains_num=PK_ROW_BYTES;
    int32_t addi_idx=-1;
    const uint8_t* src_addr=pk_ptr;
    uint8_t* dst_addr=pk_tmp;
    while(remains_num>0){
      vl=vsetvl_e8m1(remains_num);
      uint8_t addi_data=addi_idx==-1?0:pk_ptr[addi_idx];
      vuint8m1_t vsrc=vle8_v_u8m1(src_addr,vl);
      vuint8m1_t vres;
      vres=bitseqsll_vx_u8m1(vsrc,addi_data);
      vse8_v_u8m1(dst_addr,vres,vl);
      addi_idx+=vl,src_addr+=vl,dst_addr+=vl,remains_num-=vl;
    }

		b = 0;
    vl=vsetvl_e8m1(1);
    vc=vmv_s_x_u8m1(vc,0,vl);
    avl=PK_ROW_BYTES;
    const uint8_t* row_addr=pk_tmp;
    uint8_t* e_addr=&e[(SYS_N>>3)-PK_ROW_BYTES];
		while(avl>0){
      vl=vsetvl_e8m8(avl);
      va=vle8_v_u8m8(row_addr,vl);
      vb=vle8_v_u8m8(e_addr,vl);
      va=vand_vv_u8m8(va,vb,vl);
      vc=vredxor_vs_u8m8_u8m1(vc,va,vc,vl);
      row_addr+=vl,e_addr+=vl,avl-=vl;
    }
    b=vmv_x_s_u8m1_u8(vc);
    b^=(1<<(i&7)) & e[i>>3];

		b ^= b >> 4;
		b ^= b >> 2;
		b ^= b >> 1;
		b &= 1;

		s[ i>>3 ] |= (b << (i&7));

		pk_ptr += PK_ROW_BYTES;
	}
}

void encrypt_custom(unsigned char *s, const unsigned char *pk, unsigned char *e)
{
	gen_e_custom(e);

#ifdef KAT
  {
    int k;
    printf("encrypt e: positions");
    for (k = 0;k < SYS_N;++k)
      if (e[k/8] & (1 << (k&7)))
        printf(" %d",k);
    printf("\n");
  }
#endif

	syndrome_custom(s, pk, e);
}