/*
 * Copyright (c) 2023 Xinglong Yu
 */

#include "arch/riscv/insts/pqc_falcon.hh"
#include <sstream>

namespace gem5
{

namespace RiscvISA
{


std::string ChaCha20InitMacroInst::generateDisassembly(Addr pc,
        const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss << registerName(srcRegIdx(0));
    return ss.str();
}

std::string ChaCha20InitMicroInst::generateDisassembly(Addr pc,
        const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss << registerName(srcRegIdx(0));
    return ss.str();
}


std::string FalconSamplerZInst::generateDisassembly(Addr pc, 
    const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    ss << registerName(srcRegIdx(1)) << ", " << registerName(srcRegIdx(0));
    return ss.str();
}

}

}