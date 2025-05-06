#ifndef __ARCH_RISCV_INSTS_PQC_KECCAK_HH__
#define __ARCH_RISCV_INSTS_PQC_KECCAK_HH__

#include <string>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/insts/pqc_base.hh"
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
class Keccak_initBase : public RiscvStaticInst
{
  protected:
    Keccak_initBase(const char* mnem, ExtMachInst _machInst,
                   OpClass __opClass)
        : RiscvStaticInst(mnem, _machInst, __opClass)
    {
        this->flags[IsVector] = true;
    }

    std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};

class Keccak_squeezeBase : public RiscvStaticInst
{
  protected:
    Keccak_squeezeBase(const char* mnem, ExtMachInst _machInst,
                   OpClass __opClass)
        : RiscvStaticInst(mnem, _machInst, __opClass)
    {
        this->flags[IsVector] = true;
    }

    std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};

class Keccak_absorbBase : public RiscvMacroInst      
{
  protected:
    Keccak_absorbBase(const char* mnem, ExtMachInst _machInst,
                   OpClass __opClass)
        : RiscvMacroInst(mnem, _machInst, __opClass)
    {
        this->flags[IsVector] = true;
    }

    std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};

class Keccak_permuBase : public RiscvMicroInst      
{
  protected:
    Keccak_permuBase(const char* mnem, ExtMachInst _machInst,
                   OpClass __opClass)
        : RiscvMicroInst(mnem, _machInst, __opClass)
    {
        this->flags[IsVector] = true;
    }

    std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};

class Keccak_feedbackBase : public RiscvMicroInst      
{
  protected:
    Keccak_feedbackBase(const char* mnem, ExtMachInst _machInst,
                   OpClass __opClass)
        : RiscvMicroInst(mnem, _machInst, __opClass)
    {
        this->flags[IsInteger] = true;
    }

    std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};

/**
 * Base class for pqc keccak load and store macro instruction
 */
class VectorKeccakMacroInst : public VectorPqcMacroInst
{
  protected:

    VectorKeccakMacroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass)
        : VectorPqcMacroInst(mnem, _machInst, __opClass)
    {}

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;    
};

/**
 * Base class for pqc keccak load and store micro instruction
 */
class VectorKeccakMicroInst : public VectorPqcMicroInst
{
protected:

    VectorKeccakMicroInst(const char *mnem, ExtMachInst _machInst, OpClass __opClass,
                    uint16_t _microIdx,
                    uint16_t _microVl,
                    uint16_t _microSrcRegOffset,
                    uint16_t _microDstRegOffset)
        : VectorPqcMicroInst(mnem, _machInst, __opClass, _microIdx,
                             _microVl, _microSrcRegOffset, _microDstRegOffset)
    {}                    

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;    
};

}

}

#endif // __ARCH_RISCV_REGS_PQC_KECCAK_HH__

