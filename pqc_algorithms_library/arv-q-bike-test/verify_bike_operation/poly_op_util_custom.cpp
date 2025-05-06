#include <stdint.h>
#include <stdio.h>
#include "poly_op_util.h"
#include "conversions.h"
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"

void poly_add_custom(OUT uint8_t res_bin[],
    IN const uint8_t a_bin[],
    IN const uint8_t b_bin[], IN const uint32_t len)
{
    int32_t avl=len;
    size_t vl;
    const uint8_t* a_bin_addr=a_bin;
    const uint8_t* b_bin_addr=b_bin;
    uint8_t* res_bin_addr=res_bin;
    while(avl>0){
        vl=vsetvl_e8m8(avl);
        vuint8m8_t va=vle8_v_u8m8(a_bin_addr,vl);
        vuint8m8_t vb=vle8_v_u8m8(b_bin_addr,vl);
        vuint8m8_t vres=vxor_vv_u8m8(va,vb,vl);
        vse8_v_u8m8(res_bin_addr,vres,vl);

        a_bin_addr+=vl,b_bin_addr+=vl,res_bin_addr+=vl,avl-=vl;
    }
}

//only shiftleft a Bit Sequence, no wrap around
void poly_shiftleft_custom(const uint8_t poly_in[], const uint16_t offset, const uint16_t poly_len, const uint8_t addi, uint8_t poly_out[]){
    size_t vl;
    
    uint8_t addi_data;
    int32_t addi_idx=-1;

    const uint8_t* inspect_addr=poly_in;

    uint32_t chunk_shift=offset>>3;
    uint8_t* res_chunk_addr=poly_out+chunk_shift;

    int32_t remains_num=poly_len-chunk_shift;

    uint8_t bit_shift=offset&7;

    //configure the BSeqSftOffset CSR
    uint32_t old_offset;//no useful meaning
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]": [rd] "=r"(old_offset): [rt] "r"(bit_shift), [csrnum] "i"(6));

    while(remains_num>0){
        vl=vsetvl_e8m1(remains_num);
        addi_data=addi_idx==-1?addi:poly_in[addi_idx];
        vuint8m1_t v_inspect=vle8_v_u8m1(inspect_addr,vl);
        vuint8m1_t v_cycres;
        __asm__ __volatile__ ( "bitseqsll.vx %[vd], %[vt], %[op2]": [vd] "=vr"(v_cycres) : [vt] "vr"(v_inspect), [op2] "r"(addi_data));
        vse8_v_u8m1(res_chunk_addr,v_cycres,vl);

        addi_idx+=vl;
        inspect_addr+=vl;
        res_chunk_addr+=vl;
        remains_num-=vl;
    }
    for(int i=0;i<chunk_shift;i++)poly_out[i]=0;
    return;
}

void poly_split_polynomial_custom(OUT uint8_t e0[R_SIZE],
    OUT uint8_t e1[R_SIZE],
    IN const uint8_t e[2 * R_SIZE])
{
    uint8_t bits_remained = R_BITS & 7;
    uint8_t mask = (1 << bits_remained) - 1;

    int32_t avl=R_SIZE;
    size_t vl;
    const uint8_t* vsrc_addr=e;
    uint8_t* vdst_addr=e0;
    while(avl>0){
        vl=vsetvl_e8m8(avl);
        vuint8m8_t vsrc=vle8_v_u8m8(vsrc_addr,vl);
        vse8_v_u8m8(vdst_addr,vsrc,vl);

        vsrc_addr+=vl,vdst_addr+=vl,avl-=vl;
    }
    e0[R_SIZE - 1] &= mask;

    uint16_t shift_offset=8-bits_remained;
    poly_shiftleft_custom(&e[R_SIZE],shift_offset,R_SIZE,e[R_SIZE-1],e1);
    return;
}


/**********************************************************************
 *                Binary Polynomail Multiplication
 *                      (Only Support SEW=8)
**********************************************************************/
void poly_cyclic_shiftleft(const uint8_t poly_dense[R_SIZE], const uint8_t addi,uint16_t offset, uint8_t poly_res[R_SIZE]){
    size_t vl;
    
    uint8_t addi_data;
    int32_t addi_idx=-1;

    const uint8_t* inspect_addr=poly_dense;

    uint32_t chunk_shift=offset>>3;
    uint8_t* res_chunk_addr=poly_res+chunk_shift;

    int32_t remains_num=R_SIZE-chunk_shift;

    uint8_t bit_shift1=offset&7;
    uint8_t bit_shift2=bit_shift1+8-(R_BITS&7);
    bool flag=bit_shift2>=8;
    bit_shift2=flag?bit_shift2-8:bit_shift2;

    //configure the BSeqSftOffset CSR
    uint32_t old_offset;//no useful meaning
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]": [rd] "=r"(old_offset): [rt] "r"(bit_shift1), [csrnum] "i"(6));

    while(remains_num>0){
        vl=vsetvl_e8m1(remains_num);
        addi_data=addi_idx==-1?addi:poly_dense[addi_idx];
        vuint8m1_t v_inspect=vle8_v_u8m1(inspect_addr,vl);
        vuint8m1_t v_cycres;
        __asm__ __volatile__ ( "bitseqsll.vx %[vd], %[vt], %[op2]": [vd] "=vr"(v_cycres) : [vt] "vr"(v_inspect), [op2] "r"(addi_data));
        vuint8m1_t v_res=vle8_v_u8m1(res_chunk_addr,vl);
        v_res=vxor_vv_u8m1(v_res,v_cycres,vl);
        vse8_v_u8m1(res_chunk_addr,v_res,vl);

        addi_idx+=vl;
        inspect_addr+=vl;
        res_chunk_addr+=vl;
        remains_num-=vl;
    }

    remains_num=chunk_shift;
    res_chunk_addr=poly_res;
    if(flag){
        addi_idx-=1;
        inspect_addr-=1;
    }
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]": [rd] "=r"(old_offset): [rt] "r"(bit_shift2), [csrnum] "i"(6));
    while(remains_num>0){
        vl=vsetvl_e8m1(remains_num);
        addi_data=poly_dense[addi_idx];
        vuint8m1_t v_inspect=vle8_v_u8m1(inspect_addr,vl);
        vuint8m1_t v_cycres;
        __asm__ __volatile__ ( "bitseqsll.vx %[vd], %[vt], %[op2]": [vd] "=vr"(v_cycres) : [vt] "vr"(v_inspect), [op2] "r"(addi_data));
        vuint8m1_t v_res=vle8_v_u8m1(res_chunk_addr,vl);
        v_res=vxor_vv_u8m1(v_res,v_cycres,vl);
        vse8_v_u8m1(res_chunk_addr,v_res,vl);

        addi_idx+=vl;
        inspect_addr+=vl;
        res_chunk_addr+=vl;
        remains_num-=vl;
    }
}

void poly_mul_binary_custom(uint8_t poly_res[R_SIZE],const uint8_t poly_dense[R_SIZE], const uint8_t poly_sparse[R_SIZE]){
    //Convert poly_sparse to poly_sparse_compact
    uint32_t poly_sparse_compact[DV]={0};
    uint64_t start,end;//for test
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(start): );//for test
    convert2compact(poly_sparse_compact,poly_sparse);
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(end): );//for test
    printf("convert2compact cost %lu cycles\n",end-start);//for test
    //Compute the additional byte
    uint8_t addi=((((uint16_t)poly_dense[R_SIZE-1])<<8)|((uint16_t)poly_dense[R_SIZE-2]))>>(R_BITS&7);
    for(int i=0;i<DV;i++){
        poly_cyclic_shiftleft(poly_dense,addi,poly_sparse_compact[i],poly_res);
    }
    uint8_t final_mask=((1 << (R_BITS & 7)) - 1);
    poly_res[R_SIZE-1]&=final_mask;
}

void poly_mul_binary_shuffling_custom(uint8_t poly_res[R_SIZE],const uint8_t poly_dense[R_SIZE], const uint8_t poly_sparse[R_SIZE], const uint8_t shuffle_idx[DV]){
    //Convert poly_sparse to poly_sparse_compact
    uint32_t poly_sparse_compact[DV]={0};
    convert2compact(poly_sparse_compact,poly_sparse);
    //Compute the additional byte
    uint8_t addi=((((uint16_t)poly_dense[R_SIZE-1])<<8)|((uint16_t)poly_dense[R_SIZE-2]))>>(R_BITS&7);
    for(int i=0;i<DV;i++){
        poly_cyclic_shiftleft(poly_dense,addi,poly_sparse_compact[shuffle_idx[i]],poly_res);
    }
    uint8_t final_mask=((1 << (R_BITS & 7)) - 1);
    poly_res[R_SIZE-1]&=final_mask;
}


/**********************************************************************
 *                Binary Polynomail Inversion
 *                      (SEW=8)
**********************************************************************/
void gen_hmatrix_sew8(int32_t* delta_ptr, uint8_t fm, uint8_t gm,uint32_t step_size,uint8_t h[4]) {
	for (int j = 1; j <= step_size; j++) {
		uint8_t g_msb = gm >> 7;
		if (g_msb == 0) {
			(*delta_ptr) += 1;
			gm <<= 1;
			h[2] <<= 1, h[3] <<= 1;
		}
		else {
			if ((*delta_ptr) > 0) {
				(*delta_ptr) = -(*delta_ptr) + 1;
				uint8_t temp = fm;
				fm = gm, gm = (temp ^ gm) << 1;
				uint8_t h0_temp = h[0];
				uint8_t h1_temp = h[1];
				h[0] = h[2], h[1] = h[3], h[2] = (h[2] ^ h0_temp) << 1, h[3] = (h[3] ^ h1_temp) << 1;
			}
			else {
				(*delta_ptr) += 1;
				gm = (gm ^ fm) << 1;
				h[2] = (h[2] ^ h[0]) << 1, h[3] = (h[3] ^ h[1]) << 1;
			}
		}
	}
}

static uint8_t f[R_SIZE]={0};
static uint8_t g[R_SIZE]={0};
static uint8_t w[R_SIZE<<1]={0};
static uint8_t v[R_SIZE<<1]={0};
static uint8_t h0fv[R_SIZE<<1]={0};
static uint8_t h1gw[R_SIZE<<1]={0};
static uint8_t h2fv[R_SIZE<<1]={0};
static uint8_t h3gw[R_SIZE<<1]={0};
static uint8_t buffer[R_SIZE<<2]={0};

void poly_inv_binary_custom8(uint8_t gout[R_SIZE],uint8_t gin[R_SIZE]){
    int32_t avl;
    size_t vl;
    size_t tmp_vl;
    
    uint8_t shiftcnt=7 - (R_BITS & 7);
    //initial some work variables
    f[0] = 1 << shiftcnt, f[R_SIZE - 1] = 1ULL << 7;
    
    poly_shiftleft_custom(gin,shiftcnt,R_SIZE,0,g);

    w[0]=1;

    int32_t delta = 0;
	uint32_t Tau = 2 * R_BITS - 1;

	uint32_t step_size = 8-1;

    uint8_t* vsrc1_addr;
    uint8_t* vsrc2_addr;
    uint8_t* vdst_addr;

    //Compute the even and odd index (Using VLEN>>5, fixed for SEW=8 and lmul=4)
    uint8_t odd_idx[VLEN>>2];
    uint8_t even_idx[VLEN>>2];
    for(int i=0;i<(VLEN>>2);i++){
        even_idx[i]=i<<1;
        odd_idx[i]=even_idx[i]+1;
    }
    vuint8m4_t v_odd_idx=vle8_v_u8m4(odd_idx,VLEN>>2);
    vuint8m4_t v_even_idx=vle8_v_u8m4(even_idx,VLEN>>2);

    while(Tau>=step_size){
        uint8_t h[4]={1,0,0,1};
        uint8_t gm=g[R_SIZE-1];
        uint8_t fm=f[R_SIZE-1];
        gen_hmatrix_sew8(&delta,fm,gm,step_size,h);
        
        //do carry-less multiplication
        //1.h0 multiply f
        vsrc2_addr=f;
        vdst_addr=buffer;
        avl=R_SIZE;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[0]));
            tmp_vl=vsetvl_e8m8(vl<<1);
            vse8_v_u8m8(vdst_addr,vdst,tmp_vl);
            //address update
            vsrc2_addr+=vl,vdst_addr+=tmp_vl;
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h0fv[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h0fv[1];
        avl=(R_SIZE<<1)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }
        
        //2.h2 multiply f
        vsrc2_addr=f;
        vdst_addr=buffer;
        avl=R_SIZE;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[2]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h2fv[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h2fv[1];
        avl=(R_SIZE<<1)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //3.h1 multiply g
        vsrc2_addr=g;
        vdst_addr=buffer;
        avl=R_SIZE;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[1]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h1gw[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h1gw[1];
        avl=(R_SIZE<<1)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //4.h3 multiply g
        vsrc2_addr=g;
        vdst_addr=buffer;
        avl=R_SIZE;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[3]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h3gw[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h3gw[1];
        avl=(R_SIZE<<1)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //update f, xor h0fv and h1gw
        vsrc1_addr=h0fv;
        vsrc2_addr=h1gw;
        vdst_addr=f;
        avl=R_SIZE;
        while(avl>0){
            vl=vsetvl_e8m8(avl);
            vuint8m8_t vsrc1=vle8_v_u8m8(vsrc1_addr,vl);
            vuint8m8_t vsrc2=vle8_v_u8m8(vsrc2_addr,vl);
            vuint8m8_t vdst=vxor_vv_u8m8(vsrc1,vsrc2,vl);
            vse8_v_u8m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }

        //update g, xor h2fv and h3gw
        vsrc1_addr=h2fv;
        vsrc2_addr=h3gw;
        vdst_addr=g;
        avl=R_SIZE;
        while(avl>0){
            vl=vsetvl_e8m8(avl);
            vuint8m8_t vsrc1=vle8_v_u8m8(vsrc1_addr,vl);
            vuint8m8_t vsrc2=vle8_v_u8m8(vsrc2_addr,vl);
            vuint8m8_t vdst=vxor_vv_u8m8(vsrc1,vsrc2,vl);
            vse8_v_u8m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }

        //do carry-less multiplication
        //1.h0 multiply v
        vsrc2_addr=v;
        vdst_addr=buffer;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[0]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h0fv[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h0fv[1];
        avl=(R_SIZE<<2)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //2.h2 multiply v
        vsrc2_addr=v;
        vdst_addr=buffer;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[2]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h2fv[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h2fv[1];
        avl=(R_SIZE<<2)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //3.h1 multiply w
        vsrc2_addr=w;
        vdst_addr=buffer;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[1]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h1gw[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h1gw[1];
        avl=(R_SIZE<<2)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //4.h3 multiply w
        vsrc2_addr=w;
        vdst_addr=buffer;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[3]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h3gw[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h3gw[1];
        avl=(R_SIZE<<2)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //update v, xor h0fv and h1gw
        vsrc1_addr=h0fv;
        vsrc2_addr=h1gw;
        vdst_addr=v;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m8(avl);
            vuint8m8_t vsrc1=vle8_v_u8m8(vsrc1_addr,vl);
            vuint8m8_t vsrc2=vle8_v_u8m8(vsrc2_addr,vl);
            vuint8m8_t vdst=vxor_vv_u8m8(vsrc1,vsrc2,vl);
            vse8_v_u8m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }
        
        //update w, xor h2fv and h3gw
        vsrc1_addr=h2fv;
        vsrc2_addr=h3gw;
        vdst_addr=w;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m8(avl);
            vuint8m8_t vsrc1=vle8_v_u8m8(vsrc1_addr,vl);
            vuint8m8_t vsrc2=vle8_v_u8m8(vsrc2_addr,vl);
            vuint8m8_t vdst=vxor_vv_u8m8(vsrc1,vsrc2,vl);
            vse8_v_u8m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }

        Tau-=step_size;
    }
    if(Tau>0){
        uint8_t h[4]={1,0,0,1};
        uint8_t gm=g[R_SIZE-1];
        uint8_t fm=f[R_SIZE-1];
        gen_hmatrix_sew8(&delta,fm,gm,Tau,h);

        //do carry-less multiplication
        //1.h0 multiply v
        vsrc2_addr=v;
        vdst_addr=buffer;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[0]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h0fv[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h0fv[1];
        avl=(R_SIZE<<2)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //2.h2 multiply v
        vsrc2_addr=v;
        vdst_addr=buffer;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[2]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h2fv[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h2fv[1];
        avl=(R_SIZE<<2)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //3.h1 multiply w
        vsrc2_addr=w;
        vdst_addr=buffer;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[1]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h1gw[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h1gw[1];
        avl=(R_SIZE<<2)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //4.h3 multiply w
        vsrc2_addr=w;
        vdst_addr=buffer;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[3]));
            vse8_v_u8m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer, should be done with vld2q and xor in ARM NEON
        h3gw[0]=buffer[0];
        vsrc2_addr=&buffer[1];
        vdst_addr=&h3gw[1];
        avl=(R_SIZE<<2)-2;
        while(avl>0){
            vl=vsetvl_e8m4(avl);
            size_t temp_vl=vl>>1;
            vuint8m4_t vsrc2=vle8_v_u8m4(vsrc2_addr,vl);
            vuint8m4_t vsrc2_odd=vrgather_vv_u8m4(vsrc2,v_odd_idx,temp_vl);
            vuint8m4_t vsrc2_even=vrgather_vv_u8m4(vsrc2,v_even_idx,temp_vl);
            vuint8m4_t vdst=vxor_vv_u8m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse8_v_u8m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //update v, xor h0fv and h1gw
        vsrc1_addr=h0fv;
        vsrc2_addr=h1gw;
        vdst_addr=v;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m8(avl);
            vuint8m8_t vsrc1=vle8_v_u8m8(vsrc1_addr,vl);
            vuint8m8_t vsrc2=vle8_v_u8m8(vsrc2_addr,vl);
            vuint8m8_t vdst=vxor_vv_u8m8(vsrc1,vsrc2,vl);
            vse8_v_u8m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }
        
        //update w, xor h2fv and h3gw
        vsrc1_addr=h2fv;
        vsrc2_addr=h3gw;
        vdst_addr=w;
        avl=R_SIZE<<1;
        while(avl>0){
            vl=vsetvl_e8m8(avl);
            vuint8m8_t vsrc1=vle8_v_u8m8(vsrc1_addr,vl);
            vuint8m8_t vsrc2=vle8_v_u8m8(vsrc2_addr,vl);
            vuint8m8_t vdst=vxor_vv_u8m8(vsrc1,vsrc2,vl);
            vse8_v_u8m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }
    }

    //shift v R_BITS right, process in equivalent left shift
    poly_shiftleft_custom(&v[R_SIZE],(8-(R_BITS&7)),R_SIZE,v[R_SIZE-1],gout);
}

/**********************************************************************
 *                Binary Polynomail Inversion
 *                      (SEW=64)
**********************************************************************/
void gen_hmatrix_64(int32_t* delta_ptr, uint64_t fm, uint64_t gm,uint32_t step_size,uint64_t h[4]) {
	for (int j = 1; j <= step_size; j++) {
		uint64_t g_msb = gm >> 63;
		if (g_msb == 0) {
			(*delta_ptr) += 1;
			gm <<= 1;
			h[2] <<= 1, h[3] <<= 1;
		}
		else {
			if ((*delta_ptr) > 0) {
				(*delta_ptr) = -(*delta_ptr) + 1;
				uint64_t temp = fm;
				fm = gm, gm = (temp ^ gm) << 1;
				uint64_t h0_temp = h[0];
				uint64_t h1_temp = h[1];
				h[0] = h[2], h[1] = h[3], h[2] = (h[2] ^ h0_temp) << 1, h[3] = (h[3] ^ h1_temp) << 1;
			}
			else {
				(*delta_ptr) += 1;
				gm = (gm ^ fm) << 1;
				h[2] = (h[2] ^ h[0]) << 1, h[3] = (h[3] ^ h[1]) << 1;
			}
		}
	}
}

static uint64_t f_64[R_SIZE_64]={0};
static uint64_t g_64[R_SIZE_64]={0};
static uint64_t w_64[R_SIZE_64<<1]={0};
static uint64_t v_64[R_SIZE_64<<1]={0};
static uint64_t h0fv_64[R_SIZE_64<<1]={0};
static uint64_t h1gw_64[R_SIZE_64<<1]={0};
static uint64_t h2fv_64[R_SIZE_64<<1]={0};
static uint64_t h3gw_64[R_SIZE_64<<1]={0};
static uint64_t buffer_64[R_SIZE_64<<2]={0};

void poly_inv_binary_custom64(uint64_t gout[R_SIZE_64],uint64_t gin[R_SIZE_64]){
    int32_t avl;
    size_t vl;
    size_t tmp_vl;
    
    uint8_t shiftcnt=63 - (R_BITS & 63);
    //initial some work variables
    f_64[0] = 1 << shiftcnt, f_64[R_SIZE_64 - 1] = 1ULL << 63;
    
    poly_shiftleft_custom((uint8_t*)gin,shiftcnt,R_SIZE_64<<3,0,(uint8_t*)g_64);

    w_64[0]=1;

    int32_t delta = 0;
	uint32_t Tau = 2 * R_BITS - 1;

	uint32_t step_size = 64-1;

    uint64_t* vsrc1_addr;
    uint64_t* vsrc2_addr;
    uint64_t* vdst_addr;

    //Compute the even and odd index (Using VLEN>>5, fixed for SEW=64 and lmul=4)
    uint64_t odd_idx[VLEN>>5];
    uint64_t even_idx[VLEN>>5];
    for(int i=0;i<(VLEN>>5);i++){
        even_idx[i]=i<<1;
        odd_idx[i]=even_idx[i]+1;
    }
    vuint64m4_t v_odd_idx=vle64_v_u64m4(odd_idx,VLEN>>5);
    vuint64m4_t v_even_idx=vle64_v_u64m4(even_idx,VLEN>>5);

    while(Tau>=step_size){
        uint64_t h[4]={1,0,0,1};
        uint64_t gm=g_64[R_SIZE_64-1];
        uint64_t fm=f_64[R_SIZE_64-1];
        gen_hmatrix_64(&delta,fm,gm,step_size,h);
        
        //do carry-less multiplication
        //1.h0 multiply f_64
        vsrc2_addr=f_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[0]));
            tmp_vl=vsetvl_e64m8(vl<<1);
            vse64_v_u64m8(vdst_addr,vdst,tmp_vl);
            //address update
            vsrc2_addr+=vl,vdst_addr+=tmp_vl;
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h0fv_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h0fv_64[1];
        avl=(R_SIZE_64<<1)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }
        
        //2.h2 multiply f_64
        vsrc2_addr=f_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[2]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h2fv_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h2fv_64[1];
        avl=(R_SIZE_64<<1)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //3.h1 multiply g_64
        vsrc2_addr=g_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[1]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h1gw_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h1gw_64[1];
        avl=(R_SIZE_64<<1)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //4.h3 multiply g_64
        vsrc2_addr=g_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[3]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h3gw_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h3gw_64[1];
        avl=(R_SIZE_64<<1)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //update f_64, xor h0fv_64 and h1gw_64
        vsrc1_addr=h0fv_64;
        vsrc2_addr=h1gw_64;
        vdst_addr=f_64;
        avl=R_SIZE_64;
        while(avl>0){
            vl=vsetvl_e64m8(avl);
            vuint64m8_t vsrc1=vle64_v_u64m8(vsrc1_addr,vl);
            vuint64m8_t vsrc2=vle64_v_u64m8(vsrc2_addr,vl);
            vuint64m8_t vdst=vxor_vv_u64m8(vsrc1,vsrc2,vl);
            vse64_v_u64m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }

        //update g_64, xor h2fv_64 and h3gw_64
        vsrc1_addr=h2fv_64;
        vsrc2_addr=h3gw_64;
        vdst_addr=g_64;
        avl=R_SIZE_64;
        while(avl>0){
            vl=vsetvl_e64m8(avl);
            vuint64m8_t vsrc1=vle64_v_u64m8(vsrc1_addr,vl);
            vuint64m8_t vsrc2=vle64_v_u64m8(vsrc2_addr,vl);
            vuint64m8_t vdst=vxor_vv_u64m8(vsrc1,vsrc2,vl);
            vse64_v_u64m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }

        //do carry-less multiplication
        //1.h0 multiply v_64
        vsrc2_addr=v_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[0]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h0fv_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h0fv_64[1];
        avl=(R_SIZE_64<<2)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //2.h2 multiply v_64
        vsrc2_addr=v_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[2]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h2fv_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h2fv_64[1];
        avl=(R_SIZE_64<<2)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //3.h1 multiply w_64
        vsrc2_addr=w_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[1]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h1gw_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h1gw_64[1];
        avl=(R_SIZE_64<<2)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //4.h3 multiply w_64
        vsrc2_addr=w_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[3]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h3gw_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h3gw_64[1];
        avl=(R_SIZE_64<<2)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //update v_64, xor h0fv_64 and h1gw_64
        vsrc1_addr=h0fv_64;
        vsrc2_addr=h1gw_64;
        vdst_addr=v_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m8(avl);
            vuint64m8_t vsrc1=vle64_v_u64m8(vsrc1_addr,vl);
            vuint64m8_t vsrc2=vle64_v_u64m8(vsrc2_addr,vl);
            vuint64m8_t vdst=vxor_vv_u64m8(vsrc1,vsrc2,vl);
            vse64_v_u64m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }
        
        //update w_64, xor h2fv_64 and h3gw_64
        vsrc1_addr=h2fv_64;
        vsrc2_addr=h3gw_64;
        vdst_addr=w_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m8(avl);
            vuint64m8_t vsrc1=vle64_v_u64m8(vsrc1_addr,vl);
            vuint64m8_t vsrc2=vle64_v_u64m8(vsrc2_addr,vl);
            vuint64m8_t vdst=vxor_vv_u64m8(vsrc1,vsrc2,vl);
            vse64_v_u64m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }

        Tau-=step_size;
    }
    if(Tau>0){
        uint64_t h[4]={1,0,0,1};
        uint64_t gm=g_64[R_SIZE_64-1];
        uint64_t fm=f_64[R_SIZE_64-1];
        gen_hmatrix_64(&delta,fm,gm,Tau,h);

        //do carry-less multiplication
        //1.h0 multiply v_64
        vsrc2_addr=v_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[0]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h0fv_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h0fv_64[1];
        avl=(R_SIZE_64<<2)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //2.h2 multiply v_64
        vsrc2_addr=v_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[2]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h2fv_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h2fv_64[1];
        avl=(R_SIZE_64<<2)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //3.h1 multiply w_64
        vsrc2_addr=w_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[1]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h1gw_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h1gw_64[1];
        avl=(R_SIZE_64<<2)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //4.h3 multiply w_64
        vsrc2_addr=w_64;
        vdst_addr=buffer_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m8_t vdst;
            __asm__ __volatile__ ("clmul.vx %[vd],%[vt],%[rs]":[vd]"=&vr"(vdst):[vt]"vr"(vsrc2),[rs]"r"(h[3]));
            vse64_v_u64m8(vdst_addr,vdst,vl<<1);
            //address update
            vsrc2_addr+=vl,vdst_addr+=(vl<<1);
            avl-=vl;
        }
            //reduce value in the buffer_64, should be done with vld2q and xor in ARM NEON
        h3gw_64[0]=buffer_64[0];
        vsrc2_addr=&buffer_64[1];
        vdst_addr=&h3gw_64[1];
        avl=(R_SIZE_64<<2)-2;
        while(avl>0){
            vl=vsetvl_e64m4(avl);
            size_t temp_vl=vl>>1;
            vuint64m4_t vsrc2=vle64_v_u64m4(vsrc2_addr,vl);
            vuint64m4_t vsrc2_odd=vrgather_vv_u64m4(vsrc2,v_odd_idx,temp_vl);
            vuint64m4_t vsrc2_even=vrgather_vv_u64m4(vsrc2,v_even_idx,temp_vl);
            vuint64m4_t vdst=vxor_vv_u64m4(vsrc2_odd,vsrc2_even,temp_vl);
            vse64_v_u64m4(vdst_addr,vdst,temp_vl);
            vdst_addr+=temp_vl,vsrc2_addr+=vl,avl-=vl;
        }

        //update v_64, xor h0fv_64 and h1gw_64
        vsrc1_addr=h0fv_64;
        vsrc2_addr=h1gw_64;
        vdst_addr=v_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m8(avl);
            vuint64m8_t vsrc1=vle64_v_u64m8(vsrc1_addr,vl);
            vuint64m8_t vsrc2=vle64_v_u64m8(vsrc2_addr,vl);
            vuint64m8_t vdst=vxor_vv_u64m8(vsrc1,vsrc2,vl);
            vse64_v_u64m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }
        
        //update w_64, xor h2fv_64 and h3gw_64
        vsrc1_addr=h2fv_64;
        vsrc2_addr=h3gw_64;
        vdst_addr=w_64;
        avl=R_SIZE_64<<1;
        while(avl>0){
            vl=vsetvl_e64m8(avl);
            vuint64m8_t vsrc1=vle64_v_u64m8(vsrc1_addr,vl);
            vuint64m8_t vsrc2=vle64_v_u64m8(vsrc2_addr,vl);
            vuint64m8_t vdst=vxor_vv_u64m8(vsrc1,vsrc2,vl);
            vse64_v_u64m8(vdst_addr,vdst,vl);
            //address update
            vsrc1_addr+=vl,vsrc2_addr+=vl,vdst_addr+=vl;
            avl-=vl;
        }
    }

    //shift v_64 R_BITS right, process in equivalent left shift
    uint8_t byte_idx=R_SIZE-((R_SIZE_64-1)<<3);
    poly_shiftleft_custom((uint8_t*)(&v_64[R_SIZE_64-1])+byte_idx,(8-(R_BITS&7)),R_SIZE,*((uint8_t*)(&v_64[R_SIZE_64-1])+byte_idx-1),(uint8_t*)gout);
}

void poly_inv_binary_custom64_wrapper(uint8_t gout[R_SIZE],uint8_t gin[R_SIZE]){
    uint64_t gin_temp[R_SIZE_64]={0};
    uint64_t gout_temp[R_SIZE_64]={0};

    uint8_t* gin_addr=gin;
    uint64_t* gin_temp_addr=gin_temp;
    int32_t avl=R_SIZE;
    size_t vl;
    while(avl>0){
        vl=vsetvl_e8m8(avl);
        vuint8m8_t vgin=vle8_v_u8m8(gin_addr,vl);
        vuint64m8_t vgin_temp=vreinterpret_v_u8m8_u64m8(vgin);
        size_t temp_vl=(vl+7)>>3;
        vse64_v_u64m8(gin_temp_addr,vgin_temp,temp_vl);

        gin_addr+=vl,gin_temp_addr+=temp_vl,avl-=vl;
    }

    poly_inv_binary_custom64(gout_temp,gin_temp);

    uint8_t* gout_addr=gout;
    uint64_t* gout_temp_addr=gout_temp;
    avl=R_SIZE;
    while(avl>0){
        vl=vsetvl_e8m8(avl);
        size_t temp_vl=(vl+7)>>3;
        vuint64m8_t vgout_temp=vle64_v_u64m8(gout_temp_addr,temp_vl);
        vuint8m8_t vgout=vreinterpret_v_u64m8_u8m8(vgout_temp);
        vse8_v_u8m8(gout_addr,vgout,vl);

        gout_addr+=vl,gout_temp_addr+=temp_vl,avl-=vl;
    }
}