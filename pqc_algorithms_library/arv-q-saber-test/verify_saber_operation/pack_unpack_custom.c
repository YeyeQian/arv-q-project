#include "pack_unpack.h"
#include "api.h"
#include <string.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

void POLT2BS_custom(uint8_t bytes[SABER_SCALEBYTES_KEM], const uint16_t data[SABER_N])
{
#if SABER_ET == 3
	size_t vl;
    size_t avl=SABER_N;
    uint16_t* a_ptr=data;
    uint8_t* r_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        vt=vle16_v_u16m1(a_ptr,vl);
        size_t vl_packed=((vl<<1)+vl)>>3;
        vx=pack_vx_u16m1(vt,3);
        a_ptr+=vl;
        avl-=vl;
        vse8_v_u8m1(r_ptr,vx,vl_packed);
        r_ptr+=vl_packed;
    }
#elif SABER_ET == 4
	size_t vl;
    size_t avl=SABER_N;
    uint16_t* a_ptr=data;
    uint8_t* r_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        vt=vle16_v_u16m1(a_ptr,vl);
        size_t vl_packed=vl>>1;
        vx=pack_vx_u16m1(vt,4);
        a_ptr+=vl;
        avl-=vl;
        vse8_v_u8m1(r_ptr,vx,vl_packed);
        r_ptr+=vl_packed;
    }
#elif SABER_ET == 6
	size_t vl;
    size_t avl=SABER_N;
    uint16_t* a_ptr=data;
    uint8_t* r_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        vt=vle16_v_u16m1(a_ptr,vl);
        size_t vl_packed=((vl<<2)+(vl<<1))>>3;
        vx=pack_vx_u16m1(vt,6);
        a_ptr+=vl;
        avl-=vl;
        vse8_v_u8m1(r_ptr,vx,vl_packed);
        r_ptr+=vl_packed;
    }
#else
#error "Unsupported SABER parameter."
#endif
}

void BS2POLT_custom(const uint8_t bytes[SABER_SCALEBYTES_KEM], uint16_t data[SABER_N])
{
#if SABER_ET == 3
	size_t vl;
    size_t avl=SABER_N;
    uint16_t* r_ptr=data;
    uint8_t* a_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        size_t vl_packed=((vl<<1)+vl)>>3;
        vx=vle8_v_u8m1(a_ptr,vl_packed);
        vsetvl_e16m1_wrapper(avl);
        vt=unpack_vx_u16m1(vx,3);
        vse16_v_u16m1(r_ptr,vt,vl);
        r_ptr+=vl;
        a_ptr+=vl_packed;
        avl-=vl;
    }
#elif SABER_ET == 4
	size_t vl;
    size_t avl=SABER_N;
    uint16_t* r_ptr=data;
    uint8_t* a_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        size_t vl_packed=vl>>1;
        vx=vle8_v_u8m1(a_ptr,vl_packed);
        vsetvl_e16m1_wrapper(avl);
        vt=unpack_vx_u16m1(vx,4);
        vse16_v_u16m1(r_ptr,vt,vl);
        r_ptr+=vl;
        a_ptr+=vl_packed;
        avl-=vl;
    }
#elif SABER_ET == 6
	size_t vl;
    size_t avl=SABER_N;
    uint16_t* r_ptr=data;
    uint8_t* a_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        size_t vl_packed=((vl<<2)+(vl<<1))>>3;
        vx=vle8_v_u8m1(a_ptr,vl_packed);
        vsetvl_e16m1_wrapper(avl);
        vt=unpack_vx_u16m1(vx,6);
        vse16_v_u16m1(r_ptr,vt,vl);
        r_ptr+=vl;
        a_ptr+=vl_packed;
        avl-=vl;
    }
#else
#error "Unsupported SABER parameter."
#endif
}

static void POLq2BS_custom(uint8_t bytes[SABER_POLYBYTES], const uint16_t data[SABER_N])
{
	size_t vl;
    size_t avl=SABER_N;
    uint16_t* a_ptr=data;
    uint8_t* r_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        vt=vle16_v_u16m1(a_ptr,vl);
        size_t vl_packed=((vl<<3)+(vl<<2)+vl)>>3;
        vx=pack_vx_u16m1(vt,13);
        a_ptr+=vl;
        avl-=vl;
        vse8_v_u8m1(r_ptr,vx,vl_packed);
        r_ptr+=vl_packed;
    }
}

static void BS2POLq_custom(const uint8_t bytes[SABER_POLYBYTES], uint16_t data[SABER_N])
{
	size_t vl;
    size_t avl=SABER_N;
    uint16_t* r_ptr=data;
    uint8_t* a_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        size_t vl_packed=((vl<<3)+(vl<<2)+vl)>>3;
        vx=vle8_v_u8m1(a_ptr,vl_packed);
        vsetvl_e16m1_wrapper(avl);
        vt=unpack_vx_u16m1(vx,13);
        vse16_v_u16m1(r_ptr,vt,vl);
        r_ptr+=vl;
        a_ptr+=vl_packed;
        avl-=vl;
    }
}

static void POLp2BS_custom(uint8_t bytes[SABER_POLYCOMPRESSEDBYTES], const uint16_t data[SABER_N])
{
	size_t vl;
    size_t avl=SABER_N;
    uint16_t* a_ptr=data;
    uint8_t* r_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        vt=vle16_v_u16m1(a_ptr,vl);
        size_t vl_packed=((vl<<3)+(vl<<1))>>3;
        vx=pack_vx_u16m1(vt,10);
        a_ptr+=vl;
        avl-=vl;
        vse8_v_u8m1(r_ptr,vx,vl_packed);
        r_ptr+=vl_packed;
    }
}

static void BS2POLp_custom(const uint8_t bytes[SABER_POLYCOMPRESSEDBYTES], uint32_t data[SABER_N])
{
	size_t vl;
    size_t avl=SABER_N;
    uint32_t* r_ptr=data;
    uint8_t* a_ptr=bytes;
    vuint32m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e32m1(avl);
        size_t vl_packed=((vl<<3)+(vl<<1))>>3;
        vx=vle8_v_u8m1(a_ptr,vl_packed);
        vsetvl_e32m1_wrapper(avl);
        vt=unpack_vx_u32m1(vx,10);
        vse32_v_u32m1(r_ptr,vt,vl);
        r_ptr+=vl;
        a_ptr+=vl_packed;
        avl-=vl;
    }
}

void POLVECq2BS_custom(uint8_t bytes[SABER_POLYVECBYTES], const uint16_t data[SABER_L][SABER_N])
{
	size_t i;
	for (i = 0; i < SABER_L; i++)
	{
		POLq2BS_custom(bytes + i * SABER_POLYBYTES, data[i]);
	}
}

void BS2POLVECq_custom(const uint8_t bytes[SABER_POLYVECBYTES], uint16_t data[SABER_L][SABER_N])
{
	size_t i;
	for (i = 0; i < SABER_L; i++)
	{
		BS2POLq_custom(bytes + i * SABER_POLYBYTES, data[i]);
	}
}

void POLVECp2BS_custom(uint8_t bytes[SABER_POLYVECCOMPRESSEDBYTES], const uint16_t data[SABER_L][SABER_N])
{
	size_t i;
	for (i = 0; i < SABER_L; i++)
	{
		POLp2BS_custom(bytes + i * (SABER_EP * SABER_N / 8), data[i]);
	}
}

void BS2POLVECp_custom(const uint8_t bytes[SABER_POLYVECCOMPRESSEDBYTES], uint32_t data[SABER_L][SABER_N])
{
	size_t i;
	for (i = 0; i < SABER_L; i++)
	{
		BS2POLp_custom(bytes + i * (SABER_EP * SABER_N / 8), data[i]);
	}
}

void POLmsg2BS_custom(uint8_t bytes[SABER_KEYBYTES], const uint16_t data[SABER_N])
{
	memset(bytes, 0, SABER_KEYBYTES);

	size_t vl;
    size_t avl=SABER_N;
    uint16_t* a_ptr=data;
    uint8_t* r_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        vt=vle16_v_u16m1(a_ptr,vl);
        size_t vl_packed=vl>>3;
        vx=pack_vx_u16m1(vt,1);
        a_ptr+=vl;
        avl-=vl;
        vse8_v_u8m1(r_ptr,vx,vl_packed);
        r_ptr+=vl_packed;
    }
}

void BS2POLmsg_custom(const uint8_t bytes[SABER_KEYBYTES], uint16_t data[SABER_N])
{
    size_t vl;
    size_t avl=SABER_N;
    uint16_t* r_ptr=data;
    uint8_t* a_ptr=bytes;
    vuint16m1_t vt;
    vuint8m1_t vx;
    while(avl>0){
        vl=vsetvl_e16m1(avl);
        size_t vl_packed=vl>>3;
        vx=vle8_v_u8m1(a_ptr,vl_packed);
        vsetvl_e16m1_wrapper(avl);
        vt=unpack_vx_u16m1(vx,1);
        vse16_v_u16m1(r_ptr,vt,vl);
        r_ptr+=vl;
        a_ptr+=vl_packed;
        avl-=vl;
    }
}