/******************************************************************************
 * BIKE -- Bit Flipping Key Encapsulation
 *
 * Copyright (c) 2021 Nir Drucker, Shay Gueron, Rafael Misoczki, Tobias Oder,
 * Tim Gueneysu, Jan Richter-Brockmann.
 * Contact: drucker.nir@gmail.com, shay.gueron@gmail.com,
 * rafaelmisoczki@google.com, tobias.oder@rub.de, tim.gueneysu@rub.de,
 * jan.richter-brockmann@rub.de.
 *
 * Permission to use this code for BIKE is granted.
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice,
 *   this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 *
 * * The names of the contributors may not be used to endorse or promote
 *   products derived from this software without specific prior written
 *   permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ""AS IS"" AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS CORPORATION OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ******************************************************************************/


#ifndef _R_DECAPS_H_
#define _R_DECAPS_H_

#include "types.h"
#include "conversions.h"
#include <riscv_vector.h>
#include <stdio.h>

// transpose a row into a column:
_INLINE_ void transpose(uint8_t col[R_BITS], uint8_t row[R_BITS])
{
    col[0] = row[0];
    for (uint64_t i = 1; i < R_BITS ; ++i)
    {
        col[i] = row[(R_BITS) - i];
    }
}

_INLINE_ void transpose_custom(uint8_t col[R_BITS], uint8_t row[R_BITS])
{
    col[0] = row[0];
    size_t vl=vsetvl_e8m4(R_BITS);
    size_t vlmax=vl;
    vuint8m4_t vidx=vid_v_u8m4(vlmax);
    vidx=vrsub_vx_u8m4(vidx,vlmax-1,vlmax);
    int32_t avl=R_BITS-1-vl;
    uint8_t* row_addr=row+R_BITS-vlmax;
    uint8_t* col_addr=col+1;
    while(vl==vlmax){
        vuint8m4_t vrow=vle8_v_u8m4(row_addr,vl);
        vuint8m4_t vrow_permu=vrgather_vv_u8m4(vrow,vidx,vl);
        vse8_v_u8m4(col_addr,vrow_permu,vl);

        col_addr+=vl;
        vl=vsetvl_e8m4(avl);
        row_addr-=vl,avl-=vl;
    }
    for(int i=0;i<vl;i++){
        *(col_addr+i)=*(row_addr+vl-1-i);
    }
}

// Count number of 1's in tmp:
uint32_t getHammingWeight(const uint8_t tmp[R_BITS], const uint32_t length);

int BGF_decoder(uint8_t e[R_BITS*2],
        uint8_t s[R_BITS],
        uint32_t h0_compact[DV],
        uint32_t h1_compact[DV]);

uint32_t ctr(
        uint32_t h_compact_col[DV],
        int position,
        uint8_t s[R_BITS]);

void recompute_syndrome(uint8_t s[R_BITS],
        const uint32_t pos,
        const uint32_t h0_compact[DV],
        const uint32_t h1_compact[DV]);

void BFMaskedIter(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint8_t mask[R_BITS*2],
    uint32_t T,
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV],
    uint32_t h0_compact_col[DV],
    uint32_t h1_compact_col[DV]);

void BFIter(uint8_t e[R_BITS*2],
    uint8_t black[R_BITS*2],
    uint8_t gray[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t T,
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV],
    uint32_t h0_compact_col[DV],
    uint32_t h1_compact_col[DV]);

void flipAdjustedErrorPosition(uint8_t e[R_BITS*2], uint32_t position);

void getCol(
        uint32_t h_compact_col[DV],
        uint32_t h_compact_row[DV]);

//=================================================================================

uint32_t getHammingWeight_custom(const uint8_t tmp[R_BITS], const uint32_t length);

uint32_t ctr_custom(
        uint32_t h_compact_col[DV],
        int position,
        uint8_t s[R_BITS]);

void recompute_syndrome_custom(uint8_t s[R_BITS],
        const uint32_t pos,
        const uint32_t h0_compact[DV],
        const uint32_t h1_compact[DV]);

void BFMaskedIter_custom(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint8_t mask[R_BITS*2],
    uint32_t T,
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV],
    uint32_t h0_compact_col[DV],
    uint32_t h1_compact_col[DV]);

void BFIter_custom(uint8_t e[R_BITS*2],
    uint8_t black[R_BITS*2],
    uint8_t gray[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t T,
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV],
    uint32_t h0_compact_col[DV],
    uint32_t h1_compact_col[DV]);

int BGF_decoder_custom(uint8_t e[R_BITS*2],
    uint8_t s[R_BITS],
    uint32_t h0_compact[DV],
    uint32_t h1_compact[DV]);

#endif //_R_DECAPS_H_
