#ifndef __ARCH_RISCV_INSTS_PQC_ARITH_HH__
#define __ARCH_RISCV_INSTS_PQC_ARITH_HH__

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
 * Base class for pqc modular arithmetic macro instructions
 */
class VectorModArithMacroInst : public VectorPqcMacroInst
{
  protected:
    uint32_t elemPerVecReg;     // the number of register elements of current data type in one vector register
    uint16_t vlane;             // local parallism vlane, for 32-bit data width, vlane=VLANE_BASE; for 16-bit data width, vlane=2*VLANE_BASE

    VectorModArithMacroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass)
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
 * Base class for pqc modular arithmetic micro instructions
 */
class VectorModArithMicroInst : public VectorPqcMicroInst
{
protected:
    uint32_t elemPerVecReg;     // the number of register elements of current data type in one vector register
    uint16_t vlane;             // local parallism vlane, for 32-bit data width, vlane=VLANE_BASE; for 16-bit data width, vlane=2*VLANE_BASE

    VectorModArithMicroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass,
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

/**
 * Base class for pqc modular base multiplication 2x2 macro inst
 */
class VectorModBaseMulMacroInst : public VectorPqcMacroInst
{
  protected:
    uint32_t elemPerVecReg;     // the number of register elements of current data type in one vector register
    uint16_t vlane;             // local parallism vlane, for 32-bit data width, vlane=VLANE_BASE; for 16-bit data width, vlane=2*VLANE_BASE

    VectorModBaseMulMacroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass)
        : VectorPqcMacroInst(mnem, _machInst, __opClass)
    {
        if(_machInst.vtype8.vsew == 0b001) {    // data width is 16bit
            if(VLANE_BASE >= 2) {
                vlane = VLANE_BASE / 2;
            }
            else {
                vlane = 1;
            }
        }
        else {                                  // data width is 32bit
            if(VLANE_BASE >= 4) {
                vlane = VLANE_BASE / 4;
            }
            else {
                vlane = 1;
            }
        }
        elemPerVecReg = vtype_VLMAX(vtype, true);
    }

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;    
};

/**
 * Base class for pqc modular base multiplication 2x2 micro inst
 */
class VectorModBaseMulMicroInst : public VectorPqcMicroInst
{
protected:
    uint32_t elemPerVecReg;     // the number of register elements of current data type in one vector register
    uint16_t vlane;             // local parallism vlane, for 32-bit data width, vlane=VLANE_BASE; for 16-bit data width, vlane=2*VLANE_BASE

    VectorModBaseMulMicroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass,
                    uint16_t _microIdx,
                    uint16_t _microVl,
                    uint16_t _microSrcRegOffset,
                    uint16_t _microDstRegOffset)
        : VectorPqcMicroInst(mnem, _machInst, __opClass, _microIdx,
                             _microVl, _microSrcRegOffset, _microDstRegOffset)
    {
        if(_machInst.vtype8.vsew == 0b001) {    // data width is 16bit
            if(VLANE_BASE >= 2) {
                vlane = VLANE_BASE / 2;
            }
            else {
                vlane = 1;
            }
        }
        else {                                  // data width is 32bit
            if(VLANE_BASE >= 4) {
                vlane = VLANE_BASE / 4;
            }
            else {
                vlane = 1;
            }
        }
        elemPerVecReg = vtype_VLMAX(vtype, true);        
    }                    

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;    
};


class PqcNonSplitVectorInst : public RiscvStaticInst
{
  protected:
    uint32_t vl;
    uint8_t vtype;
    PqcNonSplitVectorInst(const char* mnem, ExtMachInst _machInst,
                   OpClass __opClass)
        : RiscvStaticInst(mnem, _machInst, __opClass),
        vl(_machInst.vl),
        vtype(checked_vtype(_machInst.vill, _machInst.vtype8))
    {
        this->flags[IsVector] = true;
    }

    std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};

}

}

#endif // __ARCH_RISCV_INSTS_PQC_ARITH_HH__