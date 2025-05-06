/*
 * Copyright (c) 2023 Xinglong Yu
 */

#ifndef __ARCH_RISCV_INSTS_PQC_FALCON_HH__
#define __ARCH_RISCV_INSTS_PQC_FALCON_HH__

#include <string>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/regs/misc.hh"
#include "arch/riscv/utility.hh"
#include "arch/riscv/types.hh"
#include "cpu/exec_context.hh"
#include "cpu/static_inst.hh"

#include "arch/riscv/insts/vector.hh"
#include "arch/riscv/insts/pqc_base.hh"

namespace gem5
{

namespace RiscvISA
{

/**
 * Base class for ChaCha20Init macro instructions
 */
class ChaCha20InitMacroInst : public RiscvMacroInst
{
protected:
    uint32_t vl;
    uint8_t vtype;
    
    ChaCha20InitMacroInst(const char* mnem, ExtMachInst _machInst, OpClass __opClass)
        : RiscvMacroInst(mnem, _machInst, __opClass),
        vl(_machInst.vl),
        vtype(checked_vtype(_machInst.vill, _machInst.vtype8))
    {
        this->flags[IsVector] = true;
    }

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;    
};

/**
 * Base class for ChaCha20Init micro instructions
 */
class ChaCha20InitMicroInst : public RiscvMicroInst
{
protected:
    uint16_t microIdx;          // the index of this micro inst
    uint16_t microVl;			// number of vector elements that this micro inst process
    uint16_t microSrcRegOffset;	// offset of vector source register index from the base index
    bool micro_last;            // denote if the last micro instruction
    uint8_t vtype;
    ChaCha20InitMicroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass,
                    uint16_t _microIdx,
                    uint16_t _microVl,
                    uint16_t _microSrcRegOffset,
                    bool _micro_last)
        : RiscvMicroInst(mnem, _machInst, __opClass),
        microIdx(_microIdx),
        microVl(_microVl),
        microSrcRegOffset(_microSrcRegOffset),
        micro_last(_micro_last),
        vtype(_machInst.vtype8)
    {
        this->flags[IsVector] = true;
    }            

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;    
};

class FalconSamplerZInst : public RiscvStaticInst
{
  protected:
    FalconSamplerZInst(const char* mnem, ExtMachInst _machInst,
                   OpClass __opClass)
        : RiscvStaticInst(mnem, _machInst, __opClass)
    {
        flags[IsFloating] = true;
    }

    std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};

}

}

#endif // __ARCH_RISCV_INSTS_PQC_FALCON_HH__