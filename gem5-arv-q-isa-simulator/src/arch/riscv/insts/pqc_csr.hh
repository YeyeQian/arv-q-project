#ifndef __ARCH_RISCV_INSTS_PQC_CSR_HH__
#define __ARCH_RISCV_INSTS_PQC_CSR_HH__

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

/**
 * Base class for Custom CSR Config Op
 */
class PqcCSROp : public RiscvStaticInst
{
  protected:
  	RegIndex pqc_csr_offset;   // pqc csr index offset
	RegIndex pqc_csr_index;    // pqc csr index
  	RegIndex midx = 0;
  	std::string csrName;
  	bool valid = false;
  	RegIndex pqc_csr_base_index = 0xC23;

	PqcCSROp(const char *mnem, ExtMachInst _extMachInst, OpClass __opClass)
    	: RiscvStaticInst(mnem, _extMachInst, __opClass),
	      pqc_csr_offset(_extMachInst.rs1)
	{
    	pqc_csr_index = pqc_csr_offset + pqc_csr_base_index;

    	auto csr_data_it = CSRData.find(pqc_csr_index);
    	if (csr_data_it == CSRData.end()) {
        	valid = false;
    	} else {
        	valid = true;
        	midx = csr_data_it->second.physIndex;
        	csrName = csr_data_it->second.name;
    	}    
	}

    std::string generateDisassembly(
        Addr pc, const loader::SymbolTable *symtab) const override;
};


}

}

#endif // __ARCH_RISCV_REGS_PQC_CSR_HH__

