#include <sstream>
#include <string>
#include <algorithm>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/utility.hh"
#include "cpu/static_inst.hh"
#include "arch/riscv/pqc_utility.hh"
#include "arch/riscv/regs/vector.hh"
#include "arch/riscv/insts/pqc_arith.hh"

namespace gem5
{

namespace RiscvISA
{

std::string VectorModArithMacroInst::generateDisassembly(Addr pc,
        const Loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss  << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    
    if (machInst.vm == 0) ss << ", v0.t";
    return ss.str();
}

std::string VectorModArithMicroInst::generateDisassembly(Addr pc,
        const Loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss  << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    
    if (machInst.vm == 0) ss << ", v0.t";
    return ss.str();
}

std::string VectorModBaseMulMacroInst::generateDisassembly(Addr pc,
        const Loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss  << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    ss << ", v0";

    return ss.str();
}

std::string VectorModBaseMulMicroInst::generateDisassembly(Addr pc,
        const loader::SymbolTable *symtab) const
{
    std::stringstream ss;

    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss  << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    ss << ", v0";

    return ss.str();
}

std::string
PqcNonSplitVectorInst::generateDisassembly(Addr pc,
        const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss  << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    
    return ss.str();
}

}

}