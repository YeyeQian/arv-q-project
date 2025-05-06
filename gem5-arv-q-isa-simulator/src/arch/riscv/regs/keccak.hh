/*
 * Copyright (c) 2022 PLCT Lab
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


#ifndef __ARCH_RISCV_REGS_KECCAK_HH__
#define __ARCH_RISCV_REGS_KECCAK_HH__

#include <cstdint>
#include <string>
#include <vector>
#include <cmath>

#include "arch/generic/vec_pred_reg.hh"
#include "arch/generic/vec_reg.hh"
#include "cpu/reg_class.hh"
#include "base/bitunion.hh"
#include "debug/KeccakRegs.hh"

namespace gem5
{

namespace RiscvISA
{
constexpr size_t KECCAK_REG_BSIZE = 200;
using KeccakRegContainer =
    gem5::VecRegContainer<KECCAK_REG_BSIZE>;

const int NumKeccakRegs = 2;

const std::vector<std::string> KeccakRegNames = {
    "message_reg", "state_reg"
};

static inline TypedRegClassOps<RiscvISA::KeccakRegContainer> keccakRegClassOps;

inline constexpr RegClass keccakRegClass =
    RegClass(KeccakRegClass, KeccakRegClassName, NumKeccakRegs, debug::KeccakRegs).
        ops(keccakRegClassOps).
        regType<KeccakRegContainer>();

// Keccak Related
constexpr unsigned B_WIDTH_KECCAK_PERMU = 1600;
constexpr unsigned B_WIDTH_KECCAK_PERMU_BYTE = B_WIDTH_KECCAK_PERMU / 8;
constexpr unsigned SHAKE128_RATE = 168;
constexpr unsigned SHAKE256_RATE = 136;
constexpr unsigned SHA3_256_RATE = 136;
constexpr unsigned SHA3_512_RATE = 72;
constexpr unsigned SHA3_384_RATE = 104;
constexpr uint8_t DOMAIN_SEPA_BYTE_SHAKE128 = 0x1F;
constexpr uint8_t DOMAIN_SEPA_BYTE_SHAKE256 = 0x1F;
constexpr uint8_t DOMAIN_SEPA_BYTE_SHA3256 = 0x06;
constexpr uint8_t DOMAIN_SEPA_BYTE_SHA3512 = 0x06;
constexpr uint8_t DOMAIN_SEPA_BYTE_SHA3384 = 0x06;

} // namespace RiscvISA
} // namespace gem5

#endif // __ARCH_RISCV_REGS_KECCAK_HH__