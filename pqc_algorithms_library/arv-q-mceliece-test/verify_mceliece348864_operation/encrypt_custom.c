#include "encrypt.h"

#include "util.h"
#include "params.h"

#include <stdint.h>
#include <stdio.h>
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

static void gen_e_custom(unsigned char *e)  //VLEN Must be greater equal than 256bit
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

static void gen_e_custom_fortest(unsigned char *e)  //VLEN Must be greater equal than 256bit
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
  

  // check for repetition

  eq = 0;

  for (i = 1; i < SYS_T; i++) 
    for (j = 0; j < i; j++)
            if (uint32_is_equal_declassify(ind[i],ind[j]))
            eq = 1;

	for(int i=0;i<SYS_T;i++){
    uint16_t index=ind[i]>>3;
    uint8_t pos=ind[i]&7;
    e[index]|=1<<pos;
  }
}

static void syndrome_custom(unsigned char *s, const unsigned char *pk, unsigned char *e)
{
	unsigned char b;
	const unsigned char *pk_ptr = pk;

	int i, j;
    size_t vl,avl;
    vuint8m8_t va,vb;
    vuint8m1_t vc;

	memset(s,0,SYND_BYTES);

	for (i = 0; i < PK_NROWS; i++)	
	{
		
		b = 0;
    vl=vsetvl_e8m1(1);
    vc=vmv_s_x_u8m1(vc,0,vl);
    avl=PK_ROW_BYTES;
    const uint8_t* row_addr=pk_ptr;
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

#if (VLEN == 256)
#elif (VLEN == 512)
static void syndrome_custom_asm(unsigned char *s, const unsigned char *pk, unsigned char *e)
{
	memset(s,0,SYND_BYTES);
  size_t avl1=8;
  size_t avl2=64;
  size_t avl3=128;
  size_t pack_unpack_ewidth=1;
	//firstly, load first SYND_BYTES of e into accumulator vregs as follows
  //(v2,v3),(v4,v5),(v6,v7),(v8,v9),(v10,11),(v12,v13), each pair stores 128bytes, total 6 pairs stores 768bytes=PK_NROWS bytes
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl1],e8,m1,ta,mu\n"
    "vle8.v v2, (%[addr0])\n"
    "vle8.v v3, (%[addr1])\n"
    "vle8.v v4, (%[addr2])\n"
    "vle8.v v5, (%[addr3])\n"
    "vle8.v v6, (%[addr4])\n"
    "vle8.v v7, (%[addr5])\n"
    "vle8.v v8, (%[addr6])\n"
    "vle8.v v9, (%[addr7])\n"
    "vle8.v v10,(%[addr8])\n"
    "vle8.v v11, (%[addr9])\n"
    "vle8.v v12, (%[addr10])\n"
    "vle8.v v13, (%[addr11])\n"
    "vsetvli	zero,%[avl2],e8,m1,ta,mu\n"
    "unpack.vx	v2,v2,%[pack_unpack_ewidth]\n"
    "unpack.vx	v3,v3,%[pack_unpack_ewidth]\n"
    "unpack.vx	v4,v4,%[pack_unpack_ewidth]\n"
    "unpack.vx	v5,v5,%[pack_unpack_ewidth]\n"
    "unpack.vx	v6,v6,%[pack_unpack_ewidth]\n"
    "unpack.vx	v7,v7,%[pack_unpack_ewidth]\n"
    "unpack.vx	v8,v8,%[pack_unpack_ewidth]\n"
    "unpack.vx	v9,v9,%[pack_unpack_ewidth]\n"
    "unpack.vx	v10,v10,%[pack_unpack_ewidth]\n"
    "unpack.vx	v11,v11,%[pack_unpack_ewidth]\n"
    "unpack.vx	v12,v12,%[pack_unpack_ewidth]\n"
    "unpack.vx	v13,v13,%[pack_unpack_ewidth]\n"
    :
    : [avl1]"r"(avl1), [addr0]"r"(&e[0]), [addr1]"r"(&e[8]), [addr2]"r"(&e[16]), [addr3]"r"(&e[24]), [addr4]"r"(&e[32]),
      [addr5]"r"(&e[40]), [addr6]"r"(&e[48]), [addr7]"r"(&e[56]), [addr8]"r"(&e[64]), [addr9]"r"(&e[72]), [addr10]"r"(&e[80]),
      [addr11]"r"(&e[88]), [avl2]"r"(avl2), [pack_unpack_ewidth]"r"(pack_unpack_ewidth)
  );

  //Main loop, process PK Matrix columnwise
  unsigned char * e_ptr=&e[(SYS_N>>3)-PK_ROW_BYTES];
  size_t byte_stride=PK_ROW_BYTES;
  size_t total_byte_stride=PK_ROW_BYTES<<7;

  for(int i=0;i<PK_ROW_BYTES;i++){
    const unsigned char *pk_ptr = pk+i;
    unsigned char e_coeff=e_ptr[i];

    //1st pair
    __asm__ __volatile__ (
      "vsetvli	zero,%[avl3],e8,m2,ta,mu\n"
      "vlse8.v v14,(%[pk_ptr]),%[byte_stride]\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v2,v2,v14\n"
      :
      : [avl3]"r"(avl3), [pk_ptr]"r"(pk_ptr), [byte_stride]"r"(byte_stride), [e_coeff]"r"(e_coeff)
    );

    //2nd pair
    pk_ptr+=total_byte_stride;
    __asm__ __volatile__ (
      "vlse8.v v14,(%[pk_ptr]),%[byte_stride]\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v4,v4,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [byte_stride]"r"(byte_stride), [e_coeff]"r"(e_coeff)
    );

    //3rd pair
    pk_ptr+=total_byte_stride;
    __asm__ __volatile__ (
      "vlse8.v v14,(%[pk_ptr]),%[byte_stride]\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v6,v6,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [byte_stride]"r"(byte_stride), [e_coeff]"r"(e_coeff)
    );

    //4th pair
    pk_ptr+=total_byte_stride;
    __asm__ __volatile__ (
      "vlse8.v v14,(%[pk_ptr]),%[byte_stride]\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v8,v8,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [byte_stride]"r"(byte_stride), [e_coeff]"r"(e_coeff)
    );

    //5th pair
    pk_ptr+=total_byte_stride;
    __asm__ __volatile__ (
      "vlse8.v v14,(%[pk_ptr]),%[byte_stride]\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v10,v10,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [byte_stride]"r"(byte_stride), [e_coeff]"r"(e_coeff)
    );

    //6th pair
    pk_ptr+=total_byte_stride;
    __asm__ __volatile__ (
      "vlse8.v v14,(%[pk_ptr]),%[byte_stride]\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v12,v12,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [byte_stride]"r"(byte_stride), [e_coeff]"r"(e_coeff)
    );
  }

  //Byte reduction, for 768bytes in vregs, turn each byte from 8bit to 1bit
    //1st pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v2,4\n"
    "vxor.vv v2,v2,v14\n"
    "vsrl.vi v14,v2,2\n"
    "vxor.vv v2,v2,v14\n"
    "vsrl.vi v14,v2,1\n"
    "vxor.vv v2,v2,v14\n"
    "vand.vi v2,v2,1\n"
    :
    :
  );

    //2nd pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v4,4\n"
    "vxor.vv v4,v4,v14\n"
    "vsrl.vi v14,v4,2\n"
    "vxor.vv v4,v4,v14\n"
    "vsrl.vi v14,v4,1\n"
    "vxor.vv v4,v4,v14\n"
    "vand.vi v4,v4,1\n"
    :
    :
  );

    //3rd pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v6,4\n"
    "vxor.vv v6,v6,v14\n"
    "vsrl.vi v14,v6,2\n"
    "vxor.vv v6,v6,v14\n"
    "vsrl.vi v14,v6,1\n"
    "vxor.vv v6,v6,v14\n"
    "vand.vi v6,v6,1\n"
    :
    :
  );

    //4th pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v8,4\n"
    "vxor.vv v8,v8,v14\n"
    "vsrl.vi v14,v8,2\n"
    "vxor.vv v8,v8,v14\n"
    "vsrl.vi v14,v8,1\n"
    "vxor.vv v8,v8,v14\n"
    "vand.vi v8,v8,1\n"
    :
    :
  );

    //5th pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v10,4\n"
    "vxor.vv v10,v10,v14\n"
    "vsrl.vi v14,v10,2\n"
    "vxor.vv v10,v10,v14\n"
    "vsrl.vi v14,v10,1\n"
    "vxor.vv v10,v10,v14\n"
    "vand.vi v10,v10,1\n"
    :
    :
  );

    //6th pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v12,4\n"
    "vxor.vv v12,v12,v14\n"
    "vsrl.vi v14,v12,2\n"
    "vxor.vv v12,v12,v14\n"
    "vsrl.vi v14,v12,1\n"
    "vxor.vv v12,v12,v14\n"
    "vand.vi v12,v12,1\n"
    :
    :
  );

  //Pack 768 1bit into 96 bytes, now each vreg contains 8bytes
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl2],e8,m1,ta,mu\n"
    "pack.vx  v2,v2,%[pack_unpack_ewidth]\n"
    "pack.vx  v3,v3,%[pack_unpack_ewidth]\n"
    "pack.vx	v4,v4,%[pack_unpack_ewidth]\n"
    "pack.vx	v5,v5,%[pack_unpack_ewidth]\n"
    "pack.vx	v6,v6,%[pack_unpack_ewidth]\n"
    "pack.vx	v7,v7,%[pack_unpack_ewidth]\n"
    "pack.vx	v8,v8,%[pack_unpack_ewidth]\n"
    "pack.vx	v9,v9,%[pack_unpack_ewidth]\n"
    "pack.vx	v10,v10,%[pack_unpack_ewidth]\n"
    "pack.vx	v11,v11,%[pack_unpack_ewidth]\n"
    "pack.vx	v12,v12,%[pack_unpack_ewidth]\n"
    "pack.vx	v13,v13,%[pack_unpack_ewidth]\n"
    :
    : [avl2]"r"(avl2), [pack_unpack_ewidth]"r"(pack_unpack_ewidth)
  );

  //Store back to s
  unsigned char* s_ptr=s;
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl1],e8,m1,ta,mu\n"
    "vse8.v v2,(%[addr0])\n"
    "vse8.v v3,(%[addr1])\n"
    "vse8.v v4,(%[addr2])\n"
    "vse8.v v5,(%[addr3])\n"
    "vse8.v v6,(%[addr4])\n"
    "vse8.v v7,(%[addr5])\n"
    "vse8.v v8,(%[addr6])\n"
    "vse8.v v9,(%[addr7])\n"
    "vse8.v v10,(%[addr8])\n"
    "vse8.v v11,(%[addr9])\n"
    "vse8.v v12,(%[addr10])\n"
    "vse8.v v13,(%[addr11])\n"
    :
    : [avl1]"r"(avl1), [addr0]"r"(&s_ptr[0]), [addr1]"r"(&s_ptr[8]), [addr2]"r"(&s_ptr[16]), [addr3]"r"(&s_ptr[24]), [addr4]"r"(&s_ptr[32]),
      [addr5]"r"(&s_ptr[40]), [addr6]"r"(&s_ptr[48]), [addr7]"r"(&s_ptr[56]), [addr8]"r"(&s_ptr[64]), [addr9]"r"(&s_ptr[72]), [addr10]"r"(&s_ptr[80]),
      [addr11]"r"(&s_ptr[88])
  );
}

static void syndrome_custom_asm_trans(unsigned char *s, const unsigned char *pk_trans, unsigned char *e) //input pk is transposed
{
	memset(s,0,SYND_BYTES);
  size_t avl1=8;
  size_t avl2=64;
  size_t avl3=128;
  size_t pack_unpack_ewidth=1;
	//firstly, load first SYND_BYTES of e into accumulator vregs as follows
  //(v2,v3),(v4,v5),(v6,v7),(v8,v9),(v10,11),(v12,v13), each pair stores 128bytes, total 6 pairs stores 768bytes=PK_NROWS bytes
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl1],e8,m1,ta,mu\n"
    "vle8.v v2, (%[addr0])\n"
    "vle8.v v3, (%[addr1])\n"
    "vle8.v v4, (%[addr2])\n"
    "vle8.v v5, (%[addr3])\n"
    "vle8.v v6, (%[addr4])\n"
    "vle8.v v7, (%[addr5])\n"
    "vle8.v v8, (%[addr6])\n"
    "vle8.v v9, (%[addr7])\n"
    "vle8.v v10,(%[addr8])\n"
    "vle8.v v11, (%[addr9])\n"
    "vle8.v v12, (%[addr10])\n"
    "vle8.v v13, (%[addr11])\n"
    "vsetvli	zero,%[avl2],e8,m1,ta,mu\n"
    "unpack.vx	v2,v2,%[pack_unpack_ewidth]\n"
    "unpack.vx	v3,v3,%[pack_unpack_ewidth]\n"
    "unpack.vx	v4,v4,%[pack_unpack_ewidth]\n"
    "unpack.vx	v5,v5,%[pack_unpack_ewidth]\n"
    "unpack.vx	v6,v6,%[pack_unpack_ewidth]\n"
    "unpack.vx	v7,v7,%[pack_unpack_ewidth]\n"
    "unpack.vx	v8,v8,%[pack_unpack_ewidth]\n"
    "unpack.vx	v9,v9,%[pack_unpack_ewidth]\n"
    "unpack.vx	v10,v10,%[pack_unpack_ewidth]\n"
    "unpack.vx	v11,v11,%[pack_unpack_ewidth]\n"
    "unpack.vx	v12,v12,%[pack_unpack_ewidth]\n"
    "unpack.vx	v13,v13,%[pack_unpack_ewidth]\n"
    :
    : [avl1]"r"(avl1), [addr0]"r"(&e[0]), [addr1]"r"(&e[8]), [addr2]"r"(&e[16]), [addr3]"r"(&e[24]), [addr4]"r"(&e[32]),
      [addr5]"r"(&e[40]), [addr6]"r"(&e[48]), [addr7]"r"(&e[56]), [addr8]"r"(&e[64]), [addr9]"r"(&e[72]), [addr10]"r"(&e[80]),
      [addr11]"r"(&e[88]), [avl2]"r"(avl2), [pack_unpack_ewidth]"r"(pack_unpack_ewidth)
  );

  //Main loop, process PK Matrix rowwise
  unsigned char * e_ptr=&e[(SYS_N>>3)-PK_ROW_BYTES];
  const unsigned char *pk_ptr = pk_trans;

  for(int i=0;i<PK_ROW_BYTES;i++){
    unsigned char e_coeff=e_ptr[i];

    //1st pair
    __asm__ __volatile__ (
      "vsetvli	zero,%[avl3],e8,m2,ta,mu\n"
      "vle8.v v14,(%[pk_ptr])\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v2,v2,v14\n"
      :
      : [avl3]"r"(avl3), [pk_ptr]"r"(pk_ptr), [e_coeff]"r"(e_coeff)
    );

    //2nd pair
    pk_ptr+=128;
    __asm__ __volatile__ (
      "vle8.v v14,(%[pk_ptr])\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v4,v4,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [e_coeff]"r"(e_coeff)
    );

    //3rd pair
    pk_ptr+=128;
    __asm__ __volatile__ (
      "vle8.v v14,(%[pk_ptr])\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v6,v6,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [e_coeff]"r"(e_coeff)
    );

    //4th pair
    pk_ptr+=128;
    __asm__ __volatile__ (
      "vle8.v v14,(%[pk_ptr])\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v8,v8,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [e_coeff]"r"(e_coeff)
    );

    //5th pair
    pk_ptr+=128;
    __asm__ __volatile__ (
      "vle8.v v14,(%[pk_ptr])\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v10,v10,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [e_coeff]"r"(e_coeff)
    );

    //6th pair
    pk_ptr+=128;
    __asm__ __volatile__ (
      "vle8.v v14,(%[pk_ptr])\n"
      "vand.vx v14,v14,%[e_coeff]\n"
      "vxor.vv v12,v12,v14\n"
      :
      : [pk_ptr]"r"(pk_ptr), [e_coeff]"r"(e_coeff)
    );

    pk_ptr+=128;
  }

  //Byte reduction, for 768bytes in vregs, turn each byte from 8bit to 1bit
    //1st pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v2,4\n"
    "vxor.vv v2,v2,v14\n"
    "vsrl.vi v14,v2,2\n"
    "vxor.vv v2,v2,v14\n"
    "vsrl.vi v14,v2,1\n"
    "vxor.vv v2,v2,v14\n"
    "vand.vi v2,v2,1\n"
    :
    :
  );

    //2nd pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v4,4\n"
    "vxor.vv v4,v4,v14\n"
    "vsrl.vi v14,v4,2\n"
    "vxor.vv v4,v4,v14\n"
    "vsrl.vi v14,v4,1\n"
    "vxor.vv v4,v4,v14\n"
    "vand.vi v4,v4,1\n"
    :
    :
  );

    //3rd pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v6,4\n"
    "vxor.vv v6,v6,v14\n"
    "vsrl.vi v14,v6,2\n"
    "vxor.vv v6,v6,v14\n"
    "vsrl.vi v14,v6,1\n"
    "vxor.vv v6,v6,v14\n"
    "vand.vi v6,v6,1\n"
    :
    :
  );

    //4th pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v8,4\n"
    "vxor.vv v8,v8,v14\n"
    "vsrl.vi v14,v8,2\n"
    "vxor.vv v8,v8,v14\n"
    "vsrl.vi v14,v8,1\n"
    "vxor.vv v8,v8,v14\n"
    "vand.vi v8,v8,1\n"
    :
    :
  );

    //5th pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v10,4\n"
    "vxor.vv v10,v10,v14\n"
    "vsrl.vi v14,v10,2\n"
    "vxor.vv v10,v10,v14\n"
    "vsrl.vi v14,v10,1\n"
    "vxor.vv v10,v10,v14\n"
    "vand.vi v10,v10,1\n"
    :
    :
  );

    //6th pair
  __asm__ __volatile__ (
    "vsrl.vi v14,v12,4\n"
    "vxor.vv v12,v12,v14\n"
    "vsrl.vi v14,v12,2\n"
    "vxor.vv v12,v12,v14\n"
    "vsrl.vi v14,v12,1\n"
    "vxor.vv v12,v12,v14\n"
    "vand.vi v12,v12,1\n"
    :
    :
  );

  //Pack 768 1bit into 96 bytes, now each vreg contains 8bytes
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl2],e8,m1,ta,mu\n"
    "pack.vx  v2,v2,%[pack_unpack_ewidth]\n"
    "pack.vx  v3,v3,%[pack_unpack_ewidth]\n"
    "pack.vx	v4,v4,%[pack_unpack_ewidth]\n"
    "pack.vx	v5,v5,%[pack_unpack_ewidth]\n"
    "pack.vx	v6,v6,%[pack_unpack_ewidth]\n"
    "pack.vx	v7,v7,%[pack_unpack_ewidth]\n"
    "pack.vx	v8,v8,%[pack_unpack_ewidth]\n"
    "pack.vx	v9,v9,%[pack_unpack_ewidth]\n"
    "pack.vx	v10,v10,%[pack_unpack_ewidth]\n"
    "pack.vx	v11,v11,%[pack_unpack_ewidth]\n"
    "pack.vx	v12,v12,%[pack_unpack_ewidth]\n"
    "pack.vx	v13,v13,%[pack_unpack_ewidth]\n"
    :
    : [avl2]"r"(avl2), [pack_unpack_ewidth]"r"(pack_unpack_ewidth)
  );

  //Store back to s
  unsigned char* s_ptr=s;
  __asm__ __volatile__ (
    "vsetvli	zero,%[avl1],e8,m1,ta,mu\n"
    "vse8.v v2,(%[addr0])\n"
    "vse8.v v3,(%[addr1])\n"
    "vse8.v v4,(%[addr2])\n"
    "vse8.v v5,(%[addr3])\n"
    "vse8.v v6,(%[addr4])\n"
    "vse8.v v7,(%[addr5])\n"
    "vse8.v v8,(%[addr6])\n"
    "vse8.v v9,(%[addr7])\n"
    "vse8.v v10,(%[addr8])\n"
    "vse8.v v11,(%[addr9])\n"
    "vse8.v v12,(%[addr10])\n"
    "vse8.v v13,(%[addr11])\n"
    :
    : [avl1]"r"(avl1), [addr0]"r"(&s_ptr[0]), [addr1]"r"(&s_ptr[8]), [addr2]"r"(&s_ptr[16]), [addr3]"r"(&s_ptr[24]), [addr4]"r"(&s_ptr[32]),
      [addr5]"r"(&s_ptr[40]), [addr6]"r"(&s_ptr[48]), [addr7]"r"(&s_ptr[56]), [addr8]"r"(&s_ptr[64]), [addr9]"r"(&s_ptr[72]), [addr10]"r"(&s_ptr[80]),
      [addr11]"r"(&s_ptr[88])
  );
}
#elif (VLEN == 1024)
# else
#error "VLEN must be 256/512/1024"
#endif

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

void encrypt_custom_asm(unsigned char *s, const unsigned char *pk, unsigned char *e)
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

	syndrome_custom_asm(s, pk, e);
}

void encrypt_custom_asm_trans(unsigned char *s, const unsigned char *pk_trans, unsigned char *e)
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

	syndrome_custom_asm_trans(s, pk_trans, e);
}

void encrypt_custom_fortest(unsigned char *s, const unsigned char *pk, unsigned char *e)
{
  uint64_t start, end;

  start=read_cycle();
	gen_e_custom_fortest(e);
  end=read_cycle();
  printf("gen_e_custom_fortest finished with %lu cycles\n",end-start);

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
  start=read_cycle();
	syndrome_custom_asm_trans(s, pk, e);
  end=read_cycle();
  printf("syndrome_custom_asm finished with %lu cycles\n",end-start);
}