#ifndef __ARCH_RISCV_INSTS_PQC_BUTTERFLY_HH__
#define __ARCH_RISCV_INSTS_PQC_BUTTERFLY_HH__

#include <string>
#include <assert.h>

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
 * Base class for pqc butterfly macro instruction.
 */
class ButterflyMacroInst : public VectorPqcMacroInst
{
  protected:
    uint32_t elemPerVecReg;     // the number of register elements of current data type in one vector register
    uint16_t vlane;             // local parallism vlane, for 32-bit data width, vlane=VLANE_BASE; for 16-bit data width, vlane=2*VLANE_BASE

    ButterflyMacroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass)
        : VectorPqcMacroInst(mnem, _machInst, __opClass)
    {
        if(_machInst.vtype8.vsew == 0b001) {    // data width is 16bit
            vlane = 2 * VLANE_BASE;
        }
        else {                                  // data width is 32bit
            vlane = VLANE_BASE;
        }
        elemPerVecReg = vtype_VLMAX(vtype, true);
    }

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;    
};

/**
 * Base class for pqc butterfly micro instructions 
 */
class ButterflyMicroInst : public VectorPqcMicroInst
{
protected:
    uint32_t elemPerVecReg;     // the number of register elements of current data type in one vector register
    uint16_t vlane;             // local parallism vlane, for 32-bit data width, vlane=VLANE_BASE; for 16-bit data width, vlane=2*VLANE_BASE

    ButterflyMicroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass,
                    uint16_t _microIdx,
                    uint16_t _microVl,
                    uint16_t _microSrcRegOffset,
                    uint16_t _microDstRegOffset)
        : VectorPqcMicroInst(mnem, _machInst, __opClass, _microIdx,
                             _microVl, _microSrcRegOffset, _microDstRegOffset)
    {
        if(_machInst.vtype8.vsew == 0b001) {    // data width is 16bit
            vlane = 2 * VLANE_BASE;
        }
        else {                                  // data width is 32bit
            vlane = VLANE_BASE;
        }
        elemPerVecReg = vtype_VLMAX(vtype, true);        
    }                    

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;    
};

}

}

#endif // __ARCH_RISCV_INSTS_PQC_BUTTERFLY_HH__