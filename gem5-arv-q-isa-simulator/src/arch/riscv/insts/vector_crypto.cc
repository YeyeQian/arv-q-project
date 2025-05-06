#include "arch/riscv/insts/vector_crypto.hh"
#include <sstream>
#include <string>
#include <algorithm>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/utility.hh"
#include "cpu/static_inst.hh"

#include "arch/riscv/pqc_utility.hh"

#include "arch/riscv/regs/vector.hh"

namespace gem5
{

namespace RiscvISA
{

std::string VectorCryptoArithMacroInst::generateDisassembly(Addr pc,
        const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    if (machInst.funct3 == 0b010) {
        // OPMVV
        ss  << registerName(srcRegIdx(0));
    } 
    else {    // OPIVV or OPIVX
        ss  << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    }
    if (machInst.vm == 0) 
        ss << ", v0.t";
    return ss.str();
}    

std::string VectorCryptoArithMicroInst::generateDisassembly(Addr pc,
        const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    if (machInst.funct3 == 0b010) { // OPMVV
        
        ss  << registerName(srcRegIdx(0));
    } 
    else {                          // OPIVV or OPIVX
        ss  << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    }
    if (machInst.vm == 0) 
        ss << ", v0.t";
    return ss.str();
}

}

}