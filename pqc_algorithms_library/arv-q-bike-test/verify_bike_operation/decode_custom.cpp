#include "decode.h"

#include "kem.h"
#include "sampling.h"


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <riscv_vector.h>

// count number of 1's in tmp:
uint32_t getHammingWeight_custom(const uint8_t tmp[R_BITS], const uint32_t length)
{
    uint32_t count = 0;
    int32_t avl=R_BITS;
    size_t vl;
    const uint8_t* vsrc_addr=tmp;
    while(avl>0){
        vl=vsetvl_e8m8(avl);
        vuint8m8_t vsrc=vle8_v_u8m8(vsrc_addr,vl);
        vbool1_t vmask=vmseq_vx_u8m8_b1(vsrc,1,vl);
        count+=vcpop_m_b1(vmask,vl);
        
        vsrc_addr+=vl,avl-=vl;
    }

    return count;
}

void recompute_syndrome_custom(uint8_t s[R_BITS],
        const uint32_t pos,
        const uint32_t h0_compact[DV],
        const uint32_t h1_compact[DV])
{
    size_t avl=DV;
    size_t vl=0;

    const uint32_t* addr_h0=h0_compact;
    const uint32_t* addr_h1=h1_compact;

    uint32_t temp=pos<R_BITS?pos+R_BITS:pos-R_BITS;

    vuint32m8_t vx,vy;
    vuint8m2_t vz;
    vbool4_t vf;
    
    if(pos<R_BITS){
        while(avl>0){
            vl=vsetvl_e32m8(avl);
            vx=vle32_v_u32m8(addr_h0,vl);
            vf=vmsleu_vx_u32m8_b4(vx,pos,vl);
            vy=vrsub_vx_u32m8_m(vf,vundefined_u32m8(),vx,pos,vl);
            vx=vrsub_vx_u32m8(vx,temp,vl);
            vx=vmerge_vvm_u32m8(vf,vx,vy,vl);

            vl=vsetvl_e8m2(avl);
            vz=vluxei32_v_u8m2(s,vx,vl);
            vz=vxor_vx_u8m2(vz,1,vl);
            vsuxei32_v_u8m2(s,vx,vz,vl);

            avl-=vl;
            addr_h0+=vl;
        }
    }else{
        while(avl>0){
            vl=vsetvl_e32m8(avl);
            vx=vle32_v_u32m8(addr_h1,vl);
            vf=vmsleu_vx_u32m8_b4(vx,temp,vl);
            vy=vrsub_vx_u32m8_m(vf,vundefined_u32m8(),vx,temp,vl);
            vx=vrsub_vx_u32m8(vx,pos,vl);
            vx=vmerge_vvm_u32m8(vf,vx,vy,vl);

            vl=vsetvl_e8m2(avl);
            vz=vluxei32_v_u8m2(s,vx,vl);
            vz=vxor_vx_u8m2(vz,1,vl);
            vsuxei32_v_u8m2(s,vx,vz,vl);

            avl-=vl;
            addr_h1+=vl;
        }
    }
}

uint32_t ctr_custom(
        uint32_t h_compact_col[DV],
        int position,
        uint8_t s[R_BITS])
{
    uint32_t count = 0;
    
    int32_t avl=DV;
    size_t vl;
    uint32_t* addr_h=h_compact_col;

    while(avl>0){
        vl=vsetvl_e32m8(avl);
        vuint32m8_t vx=vle32_v_u32m8(addr_h,vl);
        vx=vadd_vx_u32m8(vx,(uint32_t)position,vl);
        vbool4_t vf=vmsgeu_vx_u32m8_b4(vx,R_BITS,vl);
        vx=vsub_vx_u32m8_m(vf,vx,vx,R_BITS,vl);
        vl=vsetvl_e8m2(avl);
        vuint8m2_t vy=vluxei32_v_u8m2(s,vx,vl);
        vf=vmseq_vx_u8m2_b4(vy,1,vl);
        count+=vcpop_m_b4(vf,vl);

        addr_h+=vl,avl-=vl;
    }

    return count;
}

void BFMaskedIter_custom(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint8_t mask[R_BITS*2],
    uint32_t T,
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV],
    uint32_t h0_compact_col[DV],
    uint32_t h1_compact_col[DV])
{
    uint8_t pos[R_BITS*2] = {0};

    for (uint32_t j = 0; j < R_BITS; j++)
    {
        uint32_t counter = ctr_custom(h0_compact_col, j, s);
        if (counter >= T && mask[j])
        {
            //e[j] ^= 1;
            flipAdjustedErrorPosition(e, j);
            pos[j] = 1;
        }
    }

    for (uint32_t j = 0; j < R_BITS; j++)
    {
        uint32_t counter = ctr_custom(h1_compact_col, j, s);
        if (counter >= T && mask[R_BITS+j])
        {
            //e[R_BITS+j] ^= 1;
            flipAdjustedErrorPosition(e, R_BITS+j);
            pos[R_BITS+j] = 1;
        }
    }

    // flip bits at the end - as defined in the BGF decoder
    for(uint32_t j=0; j < 2*R_BITS; j++){
        if(pos[j] == 1){
            recompute_syndrome_custom(s, j, h0_compact, h1_compact);
        }
    }
}

void BFIter_custom(uint8_t e[R_BITS*2],
    uint8_t black[R_BITS*2],
    uint8_t gray[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t T,
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV],
    uint32_t h0_compact_col[DV],
    uint32_t h1_compact_col[DV])
{
    uint8_t pos[R_BITS*2] = {0};

    for (uint32_t j = 0; j < R_BITS; j++)
    {
        uint32_t counter = ctr_custom(h0_compact_col, j, s);
        if (counter >= T)
        {
            //e[j] ^= 1;
            flipAdjustedErrorPosition(e, j);
            pos[j] = 1;
            black[j] = 1;
        } else if(counter >= T - tau)
        {
            gray[j] = 1;
        }
    }
    for (uint32_t j = 0; j < R_BITS; j++)
    {
        uint32_t counter = ctr_custom(h1_compact_col, j, s);

        if (counter >= T)
        {
            //e[R_BITS+j] ^= 1;
            flipAdjustedErrorPosition(e, R_BITS+j);
            pos[R_BITS+j] = 1;
            black[R_BITS+j] = 1;
        } else if(counter >= T - tau)
            {
                gray[R_BITS+j] = 1;
            }
    }

    // flip bits at the end
    for(uint32_t j=0; j < 2*R_BITS; j++){
        if(pos[j] == 1){
            recompute_syndrome_custom(s, j, h0_compact, h1_compact);
        }
    }
}

int BGF_decoder_custom(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV])
{
    memset(e, 0, R_BITS*2);

    // computing the first column of each parity-check block:
    uint32_t h0_compact_col[DV] = {0};
    uint32_t h1_compact_col[DV] = {0};
    getCol(h0_compact_col, h0_compact);
    getCol(h1_compact_col, h1_compact);

    uint8_t black[R_BITS*2] = {0};
    uint8_t gray[R_BITS*2] = {0};

    for (int i = 1; i <= NbIter; i++)
    {
        memset(black, 0, R_BITS*2);
        memset(gray, 0, R_BITS*2);

        uint32_t T = floor(VAR_TH_FCT(getHammingWeight_custom(s, R_BITS)));

        BFIter_custom(e, black, gray, s, T, h0_compact, h1_compact, h0_compact_col, h1_compact_col);

        if (i == 1)
        {
            BFMaskedIter_custom(e, s, black, (DV+1)/2 + 1, h0_compact, h1_compact, h0_compact_col, h1_compact_col);
            BFMaskedIter_custom(e, s, gray, (DV+1)/2 + 1, h0_compact, h1_compact, h0_compact_col, h1_compact_col);
        }
    }
    if (getHammingWeight_custom(s, R_BITS) == 0)
        return 0; // SUCCESS
    else
        return 1; // FAILURE
}