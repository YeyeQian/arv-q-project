#ifndef __ARCH_RISCV_INSTS_VECTOR_CRYPTO_HH__
#define __ARCH_RISCV_INSTS_VECTOR_CRYPTO_HH__

#include <string>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/regs/misc.hh"
#include "arch/riscv/utility.hh"
#include "arch/riscv/types.hh"
#include "cpu/exec_context.hh"
#include "cpu/static_inst.hh"

#include "arch/riscv/insts/vector.hh"

namespace gem5
{

namespace RiscvISA
{

class VectorCryptoArithMacroInst : public VectorMacroInst
{
  protected:
    VectorCryptoArithMacroInst(const char* mnem, ExtMachInst _machInst,
                         OpClass __opClass)
        : VectorMacroInst(mnem, _machInst, __opClass)
    {}

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;
};

class VectorCryptoArithMicroInst : public VectorMicroInst
{
protected:
    VectorCryptoArithMicroInst(const char *mnem, ExtMachInst _machInst,
                         OpClass __opClass, uint16_t _microVl,
                         uint8_t _microIdx)
        : VectorMicroInst(mnem, _machInst, __opClass, _microVl, _microIdx)
    {}

    std::string generateDisassembly(
            Addr pc, const loader::SymbolTable *symtab) const override;
};

}

}

#endif // __ARCH_RISCV_INSTS_VECTOR_CRYPTO_HH__

