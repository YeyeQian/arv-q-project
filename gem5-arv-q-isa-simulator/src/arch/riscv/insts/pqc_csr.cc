#include "arch/riscv/insts/pqc_csr.hh"
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

std::string
PqcCSROp::generateDisassembly(Addr pc, const loader::SymbolTable *symtab) const
{
    std::stringstream ss;
    ss << mnemonic << ' ' << registerName(destRegIdx(0)) << ", ";
    
    auto data = CSRData.find(pqc_csr_index);
    if (data != CSRData.end())
        ss << data->second.name;
    else
        ss << "?? (" << std::hex << "0x" << pqc_csr_index << std::dec << ")";
    
    if (_numSrcRegs == 2)
        ss << ", " << registerName(srcRegIdx(0)) << ", " << registerName(srcRegIdx(1));
    else if(_numSrcRegs == 1)
        ss << ", " << registerName(srcRegIdx(0));

    return ss.str();
}

}

}