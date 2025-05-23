/*
 * Copyright (c) 2012 Google
 * Copyright (c) The University of Virginia
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer;
 * redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution;
 * neither the name of the copyright holders nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <string>
#include <assert.h>

#include "arch/riscv/decoder.hh"
#include "arch/riscv/types.hh"
#include "base/bitfield.hh"
#include "debug/Decode.hh"
#include "arch/riscv/faults.hh"

namespace gem5
{

namespace RiscvISA
{

GenericISA::BasicDecodeCache<Decoder, ExtMachInst> Decoder::defaultCache;

void Decoder::reset()
{
    aligned = true;
    mid = false;
    vConfigDone = true;
    machInst = 0;
    emi = 0;
    pqcConfigDone = true;
}

void
Decoder::moreBytes(const PCStateBase &pc, Addr fetchPC)
{
    if (GEM5_UNLIKELY(!this->vConfigDone)) {
        DPRINTF(Decode, "Waiting for vset*vl* to be executed\n");
        instDone = false;
        outOfBytes = false;
        stall = true;
        return;
    }
    else if(GEM5_UNLIKELY(!this->pqcConfigDone)) {
        DPRINTF(Decode, "Waiting for pqc configuration instruction to be executed\n");
        instDone = false;
        outOfBytes = false;
        stall = true;
        return;        
    }
    stall = false;

    // The MSB of the upper and lower halves of a machine instruction.
    constexpr size_t max_bit = sizeof(machInst) * 8 - 1;
    constexpr size_t mid_bit = sizeof(machInst) * 4 - 1;

    auto inst = letoh(machInst);
    DPRINTF(Decode, "Requesting bytes 0x%08x from address %#x\n", inst,
            fetchPC);

    bool aligned = pc.instAddr() % sizeof(machInst) == 0;
    if (aligned) {
        emi.instBits = inst;
        if (compressed(inst))
            emi.instBits = bits(inst, mid_bit, 0);
        outOfBytes = !compressed(emi);
        instDone = true;
    } else {
        if (mid) {
            assert(bits(emi.instBits, max_bit, mid_bit + 1) == 0);
            replaceBits(emi.instBits, max_bit, mid_bit + 1, inst);
            mid = false;
            outOfBytes = false;
            instDone = true;
        } else {
            emi.instBits = bits(inst, max_bit, mid_bit + 1);
            mid = !compressed(emi);
            outOfBytes = true;
            instDone = compressed(emi);
        }
    }
    if (instDone) {
        emi.vl      = this->machVl;
        emi.vtype8   = this->machVtype & 0xff;
        emi.vill    = this->machVtype.vill;
        if (vconf(emi)) {
            this->vConfigDone = false; // set true when vconfig inst execute
        }
        else if(is_pqc_conf(emi)) {
            this->pqcConfigDone = false; // set true when pqc configuration inst execute
        }
    }
}

StaticInstPtr
Decoder::decode(ExtMachInst mach_inst, Addr addr)
{
    DPRINTF(Decode, "Decoding instruction 0x%08x at address %#x\n",
            mach_inst.instBits, addr);

    StaticInstPtr si = defaultCache.decode(this, mach_inst, addr);

    DPRINTF(Decode, "Decode: Decoded %s instruction: %#x\n",
            si->getName(), mach_inst);
    return si;
}

StaticInstPtr
Decoder::decode(PCStateBase &_next_pc)
{
    if (!instDone)
        return nullptr;
    instDone = false;

    auto &next_pc = _next_pc.as<PCState>();

    if (compressed(emi)) {
        next_pc.npc(next_pc.instAddr() + sizeof(machInst) / 2);
        next_pc.compressed(true);
    } else {
        next_pc.npc(next_pc.instAddr() + sizeof(machInst));
        next_pc.compressed(false);
    }

    return decode(emi, next_pc.instAddr());
}

void
Decoder::setVlAndVtype(uint32_t vl, VTYPE vtype)
{
    this->machVtype = vtype;
    this->machVl = vl;

    this->vConfigDone = true;
}

void
Decoder::setPqcConfig(RegIndex pqc_csr_sel)
{
    bool pqc_csr_sel_valid = pqc_csr_sel == 0 || pqc_csr_sel == 1 || pqc_csr_sel == 2 || pqc_csr_sel == 3 || pqc_csr_sel == 4 || pqc_csr_sel==5 || pqc_csr_sel==6 || pqc_csr_sel==7;
    assert(pqc_csr_sel_valid);

    this->pqcConfigDone = true;
}

} // namespace RiscvISA
} // namespace gem5