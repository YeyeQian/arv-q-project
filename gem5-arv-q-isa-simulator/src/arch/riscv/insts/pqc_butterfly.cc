#include <sstream>
#include <string>
#include <algorithm>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/utility.hh"
#include "cpu/static_inst.hh"
#include "arch/riscv/pqc_utility.hh"
#include "arch/riscv/regs/vector.hh"
#include "arch/riscv/insts/pqc_butterfly.hh"

namespace gem5
{

namespace RiscvISA
{

std::string ButterflyMacroInst::generateDisassembly(Addr pc,
        const Loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss  << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    ss << ", v0";

    return ss.str();
}

std::string ButterflyMicroInst::generateDisassembly(Addr pc,
        const loader::SymbolTable *symtab) const
{
    std::stringstream ss;

    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss  << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    ss << ", v0";

    return ss.str();
}

}

}