#include "arch/riscv/insts/pqc_keccak.hh"
#include <sstream>
#include <string>
#include <algorithm>
#include <cstring>

#include "arch/riscv/insts/static_inst.hh"
#include "arch/riscv/utility.hh"
#include "cpu/static_inst.hh"

#include "arch/riscv/pqc_utility.hh"
#include "arch/riscv/regs/vector.hh"

namespace gem5
{

namespace RiscvISA
{

std::string 
Keccak_initBase::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", " << "0";
    
    return ss.str();
}


std::string 
Keccak_squeezeBase::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", " << registerName(srcRegIdx(0));
    
    return ss.str();
}


std::string 
Keccak_absorbBase::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", " << registerName(srcRegIdx(0));
    
    return ss.str();
}

std::string 
Keccak_permuBase::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss << registerName(srcRegIdx(0)) << ", " << registerName(srcRegIdx(1));
    
    return ss.str();
}

std::string 
Keccak_feedbackBase::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss << registerName(srcRegIdx(0)) << ", " << registerName(srcRegIdx(1));
    
    return ss.str();
}

std::string 
VectorKeccakMacroInst::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    if(machInst.customfunct6 == 0b000001) {         // keccak load message inst
        ss << mnemonic << ' ' << "message_reg" << ", " << registerName(srcRegIdx(0));
    }
    else if(machInst.customfunct6 == 0b000010) {    // keccak store state inst
        ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", " << "state_reg";
    }
    
    return ss.str();
}

std::string 
VectorKeccakMicroInst::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", " << registerName(srcRegIdx(0));
    
    return ss.str();
}


}

}