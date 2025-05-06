#ifndef __ARCH_RISCV_INSTS_PQC_BASE_HH__
#define __ARCH_RISCV_INSTS_PQC_BASE_HH__

#include <string>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/regs/misc.hh"
#include "arch/riscv/utility.hh"
#include "arch/riscv/types.hh"
#include "cpu/exec_context.hh"
#include "cpu/static_inst.hh"

#include <assert.h>
#include "arch/riscv/insts/vector.hh"


namespace gem5
{

namespace RiscvISA
{

class VectorPqcMacroInst : public RiscvMacroInst
{
  protected:
    uint32_t vl;
    uint8_t vtype;
    
    VectorPqcMacroInst(const char* mnem, ExtMachInst _machInst,
                   OpClass __opClass)
        : RiscvMacroInst(mnem, _machInst, __opClass),
        vl(_machInst.vl),
        vtype(checked_vtype(_machInst.vill, _machInst.vtype8))
    {
        this->flags[IsVector] = true;
    }
};

class VectorPqcMicroInst : public RiscvMicroInst
{
protected:
    uint16_t microIdx;          // the index of this micro inst
    uint16_t microVl;			// number of vector elements that this micro inst process
    uint16_t microSrcRegOffset;	// offset of vector source register index from the base index
    uint16_t microDstRegOffset; // offset of vector destination register index from the base index
    uint8_t vtype;
    VectorPqcMicroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass,
                    uint16_t _microIdx,
                    uint16_t _microVl,
                    uint16_t _microSrcRegOffset,
                    uint16_t _microDstRegOffset)
        : RiscvMicroInst(mnem, _machInst, __opClass),
        microIdx(_microIdx),
        microVl(_microVl),
        microSrcRegOffset(_microSrcRegOffset),
        microDstRegOffset(_microDstRegOffset),
        vtype(_machInst.vtype8)
    {
        this->flags[IsVector] = true;
    }
};

}

}

#endif // __ARCH_RISCV_INSTS_PQC_BASE_HH__

