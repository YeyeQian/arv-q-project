/* riscv-opc.h.  RISC-V instruction opcode and CSR macros.
   Copyright (C) 2011-2020 Free Software Foundation, Inc.
   Contributed by Andrew Waterman

   This file is part of GDB, GAS, and the GNU binutils.

   GDB, GAS, and the GNU binutils are free software; you can redistribute
   them and/or modify them under the terms of the GNU General Public
   License as published by the Free Software Foundation; either version
   3, or (at your option) any later version.

   GDB, GAS, and the GNU binutils are distributed in the hope that they
   will be useful, but WITHOUT ANY WARRANTY; without even the implied
   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See
   the GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING3. If not,
   see <http://www.gnu.org/licenses/>.  */

/*This file is in directory of "xuantie-gnu-toolchain/riscv-binutils/include/opcode"*/

#ifndef RISCV_ENCODING_H
#define RISCV_ENCODING_H
/* Instruction opcode macros.  */

// start custom instructions
#define MATCH_PQCCSRRW 0x505b
#define MASK_PQCCSRRW  0xfe00707f
#define MATCH_BFCT 0x2b
#define MASK_BFCT  0xfe00707f
#define MATCH_BFGS 0x400002b
#define MASK_BFGS  0xfe00707f    
#define MATCH_MODADDVV 0x800002b
#define MASK_MODADDVV  0xfc00707f
#define MATCH_MODSUBVV 0xc00002b
#define MASK_MODSUBVV  0xfc00707f
#define MATCH_MODMULVV 0x1000002b
#define MASK_MODMULVV  0xfc00707f
#define MATCH_MODADDVX 0x800202b
#define MASK_MODADDVX  0xfc00707f
#define MATCH_MODSUBVX 0xc00202b
#define MASK_MODSUBVX  0xfc00707f
#define MATCH_MODMULVX 0x1000202b
#define MASK_MODMULVX  0xfc00707f
#define MATCH_MODBASEMUL 0x1400002b
#define MASK_MODBASEMUL  0xfe00707f
#define MATCH_REJSAMPLE 0x1a00202b
#define MASK_REJSAMPLE  0xfe00707f
#define MATCH_PACKVX 0x1e00202b
#define MASK_PACKVX  0xfe00707f
#define MATCH_UNPACKVX 0x2200202b
#define MASK_UNPACKVX  0xfe00707f
#define MATCH_KECCAKLOADM 0x600105b
#define MASK_KECCAKLOADM  0xfe0fffff
#define MATCH_KECCAKSTORES 0xa00405b
#define MASK_KECCAKSTORES  0xfffff07f
#define MATCH_KECCAKINIT 0xc00005b
#define MASK_KECCAKINIT  0xfff07fff
#define MATCH_KECCAKABSORB 0x1000605b
#define MASK_KECCAKABSORB  0xfff0707f
#define MATCH_KECCAKSQUEEZE 0x1400005b
#define MASK_KECCAKSQUEEZE  0xfff07fff
// #define MATCH_VCPOPV 0x48072057 // this is used in vector cryptography spec
// #define MASK_VCPOPV 0xfc0ff07f
#define MATCH_VCPOPV 0x1a00505b  // this is customized
#define MASK_VCPOPV 0xfe0ff07f
#define MATCH_UNPACKSVX 0x2e00202b
#define MASK_UNPACKSVX  0xfe00707f
#define MATCH_CHACHA20INITV 0x1e00105b
#define MASK_CHACHA20INITV  0xfe0fffff
#define MATCH_SAMPLERZ 0x700b
#define MASK_SAMPLERZ  0xfe00707f
#define MATCH_BITSEQSLL 0x2600202b
#define MASK_BITSEQSLL  0xfe00707f
#define MATCH_VXCLMUL 0x2a00202b
#define MASK_VXCLMUL  0xfe00707f
#define MATCH_VVCLMUL 0x3200002b
#define MASK_VVCLMUL  0xfe00707f
// end custom instructions

#define MATCH_SLLI_RV32 0x1013
#define MASK_SLLI_RV32  0xfe00707f
#define MATCH_SRLI_RV32 0x5013
#define MASK_SRLI_RV32  0xfe00707f
#define MATCH_SRAI_RV32 0x40005013
#define MASK_SRAI_RV32  0xfe00707f
#define MATCH_FRFLAGS 0x102073
#define MASK_FRFLAGS  0xfffff07f
#define MATCH_FSFLAGS 0x101073
#define MASK_FSFLAGS  0xfff0707f
#define MATCH_FSFLAGSI 0x105073
#define MASK_FSFLAGSI  0xfff0707f
#define MATCH_FRRM 0x202073
#define MASK_FRRM  0xfffff07f
#define MATCH_FSRM 0x201073
#define MASK_FSRM  0xfff0707f
#define MATCH_FSRMI 0x205073
#define MASK_FSRMI  0xfff0707f
#define MATCH_FSCSR 0x301073
#define MASK_FSCSR  0xfff0707f
#define MATCH_FRCSR 0x302073
#define MASK_FRCSR  0xfffff07f
#define MATCH_RDCYCLE 0xc0002073
#define MASK_RDCYCLE  0xfffff07f
#define MATCH_RDTIME 0xc0102073
#define MASK_RDTIME  0xfffff07f
#define MATCH_RDINSTRET 0xc0202073
#define MASK_RDINSTRET  0xfffff07f
#define MATCH_RDCYCLEH 0xc8002073
#define MASK_RDCYCLEH  0xfffff07f
#define MATCH_RDTIMEH 0xc8102073
#define MASK_RDTIMEH  0xfffff07f
#define MATCH_RDINSTRETH 0xc8202073
#define MASK_RDINSTRETH  0xfffff07f
#define MATCH_SCALL 0x73
#define MASK_SCALL  0xffffffff
#define MATCH_SBREAK 0x100073
#define MASK_SBREAK  0xffffffff
#define MATCH_BEQ 0x63
#define MASK_BEQ  0x707f
#define MATCH_BNE 0x1063
#define MASK_BNE  0x707f
#define MATCH_BLT 0x4063
#define MASK_BLT  0x707f
#define MATCH_BGE 0x5063
#define MASK_BGE  0x707f
#define MATCH_BLTU 0x6063
#define MASK_BLTU  0x707f
#define MATCH_BGEU 0x7063
#define MASK_BGEU  0x707f
#define MATCH_JALR 0x67
#define MASK_JALR  0x707f
#define MATCH_JAL 0x6f
#define MASK_JAL  0x7f
#define MATCH_LUI 0x37
#define MASK_LUI  0x7f
#define MATCH_AUIPC 0x17
#define MASK_AUIPC  0x7f
#define MATCH_ADDI 0x13
#define MASK_ADDI  0x707f
#define MATCH_SLLI 0x1013
#define MASK_SLLI  0xfc00707f
#define MATCH_SLTI 0x2013
#define MASK_SLTI  0x707f
#define MATCH_SLTIU 0x3013
#define MASK_SLTIU  0x707f
#define MATCH_XORI 0x4013
#define MASK_XORI  0x707f
#define MATCH_SRLI 0x5013
#define MASK_SRLI  0xfc00707f
#define MATCH_SRAI 0x40005013
#define MASK_SRAI  0xfc00707f
#define MATCH_ORI 0x6013
#define MASK_ORI  0x707f
#define MATCH_ANDI 0x7013
#define MASK_ANDI  0x707f
#define MATCH_ADD 0x33
#define MASK_ADD  0xfe00707f
#define MATCH_SUB 0x40000033
#define MASK_SUB  0xfe00707f
#define MATCH_SLL 0x1033
#define MASK_SLL  0xfe00707f
#define MATCH_SLT 0x2033
#define MASK_SLT  0xfe00707f
#define MATCH_SLTU 0x3033
#define MASK_SLTU  0xfe00707f
#define MATCH_XOR 0x4033
#define MASK_XOR  0xfe00707f
#define MATCH_SRL 0x5033
#define MASK_SRL  0xfe00707f
#define MATCH_SRA 0x40005033
#define MASK_SRA  0xfe00707f
#define MATCH_OR 0x6033
#define MASK_OR  0xfe00707f
#define MATCH_AND 0x7033
#define MASK_AND  0xfe00707f
#define MATCH_ADDIW 0x1b
#define MASK_ADDIW  0x707f
#define MATCH_SLLIW 0x101b
#define MASK_SLLIW  0xfe00707f
#define MATCH_SRLIW 0x501b
#define MASK_SRLIW  0xfe00707f
#define MATCH_SRAIW 0x4000501b
#define MASK_SRAIW  0xfe00707f
#define MATCH_ADDW 0x3b
#define MASK_ADDW  0xfe00707f
#define MATCH_SUBW 0x4000003b
#define MASK_SUBW  0xfe00707f
#define MATCH_SLLW 0x103b
#define MASK_SLLW  0xfe00707f
#define MATCH_SRLW 0x503b
#define MASK_SRLW  0xfe00707f
#define MATCH_SRAW 0x4000503b
#define MASK_SRAW  0xfe00707f
#define MATCH_LB 0x3
#define MASK_LB  0x707f
#define MATCH_LH 0x1003
#define MASK_LH  0x707f
#define MATCH_LW 0x2003
#define MASK_LW  0x707f
#define MATCH_LD 0x3003
#define MASK_LD  0x707f
#define MATCH_LBU 0x4003
#define MASK_LBU  0x707f
#define MATCH_LHU 0x5003
#define MASK_LHU  0x707f
#define MATCH_LWU 0x6003
#define MASK_LWU  0x707f
#define MATCH_SB 0x23
#define MASK_SB  0x707f
#define MATCH_SH 0x1023
#define MASK_SH  0x707f
#define MATCH_SW 0x2023
#define MASK_SW  0x707f
#define MATCH_SD 0x3023
#define MASK_SD  0x707f
#define MATCH_PAUSE 0x0100000f
#define MASK_PAUSE  0xffffffff
#define MATCH_FENCE 0xf
#define MASK_FENCE  0x707f
#define MATCH_FENCE_I 0x100f
#define MASK_FENCE_I  0x707f
#define MATCH_FENCE_TSO 0x8330000f
#define MASK_FENCE_TSO  0xfff0707f
#define MATCH_MUL 0x2000033
#define MASK_MUL  0xfe00707f
#define MATCH_MULH 0x2001033
#define MASK_MULH  0xfe00707f
#define MATCH_MULHSU 0x2002033
#define MASK_MULHSU  0xfe00707f
#define MATCH_MULHU 0x2003033
#define MASK_MULHU  0xfe00707f
#define MATCH_DIV 0x2004033
#define MASK_DIV  0xfe00707f
#define MATCH_DIVU 0x2005033
#define MASK_DIVU  0xfe00707f
#define MATCH_REM 0x2006033
#define MASK_REM  0xfe00707f
#define MATCH_REMU 0x2007033
#define MASK_REMU  0xfe00707f
#define MATCH_MULW 0x200003b
#define MASK_MULW  0xfe00707f
#define MATCH_DIVW 0x200403b
#define MASK_DIVW  0xfe00707f
#define MATCH_DIVUW 0x200503b
#define MASK_DIVUW  0xfe00707f
#define MATCH_REMW 0x200603b
#define MASK_REMW  0xfe00707f
#define MATCH_REMUW 0x200703b
#define MASK_REMUW  0xfe00707f
#define MATCH_AMOADD_W 0x202f
#define MASK_AMOADD_W  0xf800707f
#define MATCH_AMOXOR_W 0x2000202f
#define MASK_AMOXOR_W  0xf800707f
#define MATCH_AMOOR_W 0x4000202f
#define MASK_AMOOR_W  0xf800707f
#define MATCH_AMOAND_W 0x6000202f
#define MASK_AMOAND_W  0xf800707f
#define MATCH_AMOMIN_W 0x8000202f
#define MASK_AMOMIN_W  0xf800707f
#define MATCH_AMOMAX_W 0xa000202f
#define MASK_AMOMAX_W  0xf800707f
#define MATCH_AMOMINU_W 0xc000202f
#define MASK_AMOMINU_W  0xf800707f
#define MATCH_AMOMAXU_W 0xe000202f
#define MASK_AMOMAXU_W  0xf800707f
#define MATCH_AMOSWAP_W 0x800202f
#define MASK_AMOSWAP_W  0xf800707f
#define MATCH_LR_W 0x1000202f
#define MASK_LR_W  0xf9f0707f
#define MATCH_SC_W 0x1800202f
#define MASK_SC_W  0xf800707f
#define MATCH_AMOADD_D 0x302f
#define MASK_AMOADD_D  0xf800707f
#define MATCH_AMOXOR_D 0x2000302f
#define MASK_AMOXOR_D  0xf800707f
#define MATCH_AMOOR_D 0x4000302f
#define MASK_AMOOR_D  0xf800707f
#define MATCH_AMOAND_D 0x6000302f
#define MASK_AMOAND_D  0xf800707f
#define MATCH_AMOMIN_D 0x8000302f
#define MASK_AMOMIN_D  0xf800707f
#define MATCH_AMOMAX_D 0xa000302f
#define MASK_AMOMAX_D  0xf800707f
#define MATCH_AMOMINU_D 0xc000302f
#define MASK_AMOMINU_D  0xf800707f
#define MATCH_AMOMAXU_D 0xe000302f
#define MASK_AMOMAXU_D  0xf800707f
#define MATCH_AMOSWAP_D 0x800302f
#define MASK_AMOSWAP_D  0xf800707f
#define MATCH_LR_D 0x1000302f
#define MASK_LR_D  0xf9f0707f
#define MATCH_SC_D 0x1800302f
#define MASK_SC_D  0xf800707f
#define MATCH_ECALL 0x73
#define MASK_ECALL  0xffffffff
#define MATCH_EBREAK 0x100073
#define MASK_EBREAK  0xffffffff
#define MATCH_URET 0x200073
#define MASK_URET  0xffffffff
#define MATCH_SRET 0x10200073
#define MASK_SRET  0xffffffff
#define MATCH_HRET 0x20200073
#define MASK_HRET  0xffffffff
#define MATCH_MRET 0x30200073
#define MASK_MRET  0xffffffff
#define MATCH_DRET 0x7b200073
#define MASK_DRET  0xffffffff
#define MATCH_SFENCE_VM 0x10400073
#define MASK_SFENCE_VM  0xfff07fff
#define MATCH_SFENCE_VMA 0x12000073
#define MASK_SFENCE_VMA  0xfe007fff
#define MATCH_WFI 0x10500073
#define MASK_WFI  0xffffffff
#define MATCH_CSRRW 0x1073
#define MASK_CSRRW  0x707f
#define MATCH_CSRRS 0x2073
#define MASK_CSRRS  0x707f
#define MATCH_CSRRC 0x3073
#define MASK_CSRRC  0x707f
#define MATCH_CSRRWI 0x5073
#define MASK_CSRRWI  0x707f
#define MATCH_CSRRSI 0x6073
#define MASK_CSRRSI  0x707f
#define MATCH_CSRRCI 0x7073
#define MASK_CSRRCI  0x707f
#define MATCH_FADD_S 0x53
#define MASK_FADD_S  0xfe00007f
#define MATCH_FSUB_S 0x8000053
#define MASK_FSUB_S  0xfe00007f
#define MATCH_FMUL_S 0x10000053
#define MASK_FMUL_S  0xfe00007f
#define MATCH_FDIV_S 0x18000053
#define MASK_FDIV_S  0xfe00007f
#define MATCH_FSGNJ_S 0x20000053
#define MASK_FSGNJ_S  0xfe00707f
#define MATCH_FSGNJN_S 0x20001053
#define MASK_FSGNJN_S  0xfe00707f
#define MATCH_FSGNJX_S 0x20002053
#define MASK_FSGNJX_S  0xfe00707f
#define MATCH_FMIN_S 0x28000053
#define MASK_FMIN_S  0xfe00707f
#define MATCH_FMAX_S 0x28001053
#define MASK_FMAX_S  0xfe00707f
#define MATCH_FSQRT_S 0x58000053
#define MASK_FSQRT_S  0xfff0007f
#define MATCH_FADD_D 0x2000053
#define MASK_FADD_D  0xfe00007f
#define MATCH_FSUB_D 0xa000053
#define MASK_FSUB_D  0xfe00007f
#define MATCH_FMUL_D 0x12000053
#define MASK_FMUL_D  0xfe00007f
#define MATCH_FDIV_D 0x1a000053
#define MASK_FDIV_D  0xfe00007f
#define MATCH_FSGNJ_D 0x22000053
#define MASK_FSGNJ_D  0xfe00707f
#define MATCH_FSGNJN_D 0x22001053
#define MASK_FSGNJN_D  0xfe00707f
#define MATCH_FSGNJX_D 0x22002053
#define MASK_FSGNJX_D  0xfe00707f
#define MATCH_FMIN_D 0x2a000053
#define MASK_FMIN_D  0xfe00707f
#define MATCH_FMAX_D 0x2a001053
#define MASK_FMAX_D  0xfe00707f
#define MATCH_FCVT_S_D 0x40100053
#define MASK_FCVT_S_D  0xfff0007f
#define MATCH_FCVT_D_S 0x42000053
#define MASK_FCVT_D_S  0xfff0007f
#define MATCH_FSQRT_D 0x5a000053
#define MASK_FSQRT_D  0xfff0007f
#define MATCH_FADD_Q 0x6000053
#define MASK_FADD_Q  0xfe00007f
#define MATCH_FSUB_Q 0xe000053
#define MASK_FSUB_Q  0xfe00007f
#define MATCH_FMUL_Q 0x16000053
#define MASK_FMUL_Q  0xfe00007f
#define MATCH_FDIV_Q 0x1e000053
#define MASK_FDIV_Q  0xfe00007f
#define MATCH_FSGNJ_Q 0x26000053
#define MASK_FSGNJ_Q  0xfe00707f
#define MATCH_FSGNJN_Q 0x26001053
#define MASK_FSGNJN_Q  0xfe00707f
#define MATCH_FSGNJX_Q 0x26002053
#define MASK_FSGNJX_Q  0xfe00707f
#define MATCH_FMIN_Q 0x2e000053
#define MASK_FMIN_Q  0xfe00707f
#define MATCH_FMAX_Q 0x2e001053
#define MASK_FMAX_Q  0xfe00707f
#define MATCH_FCVT_S_Q 0x40300053
#define MASK_FCVT_S_Q  0xfff0007f
#define MATCH_FCVT_Q_S 0x46000053
#define MASK_FCVT_Q_S  0xfff0007f
#define MATCH_FCVT_D_Q 0x42300053
#define MASK_FCVT_D_Q  0xfff0007f
#define MATCH_FCVT_Q_D 0x46100053
#define MASK_FCVT_Q_D  0xfff0007f
#define MATCH_FSQRT_Q 0x5e000053
#define MASK_FSQRT_Q  0xfff0007f
#define MATCH_FLE_S 0xa0000053
#define MASK_FLE_S  0xfe00707f
#define MATCH_FLT_S 0xa0001053
#define MASK_FLT_S  0xfe00707f
#define MATCH_FEQ_S 0xa0002053
#define MASK_FEQ_S  0xfe00707f
#define MATCH_FLE_D 0xa2000053
#define MASK_FLE_D  0xfe00707f
#define MATCH_FLT_D 0xa2001053
#define MASK_FLT_D  0xfe00707f
#define MATCH_FEQ_D 0xa2002053
#define MASK_FEQ_D  0xfe00707f
#define MATCH_FLE_Q 0xa6000053
#define MASK_FLE_Q  0xfe00707f
#define MATCH_FLT_Q 0xa6001053
#define MASK_FLT_Q  0xfe00707f
#define MATCH_FEQ_Q 0xa6002053
#define MASK_FEQ_Q  0xfe00707f
#define MATCH_FCVT_W_S 0xc0000053
#define MASK_FCVT_W_S  0xfff0007f
#define MATCH_FCVT_WU_S 0xc0100053
#define MASK_FCVT_WU_S  0xfff0007f
#define MATCH_FCVT_L_S 0xc0200053
#define MASK_FCVT_L_S  0xfff0007f
#define MATCH_FCVT_LU_S 0xc0300053
#define MASK_FCVT_LU_S  0xfff0007f
#define MATCH_FMV_X_S 0xe0000053
#define MASK_FMV_X_S  0xfff0707f
#define MATCH_FCLASS_S 0xe0001053
#define MASK_FCLASS_S  0xfff0707f
#define MATCH_FCVT_W_D 0xc2000053
#define MASK_FCVT_W_D  0xfff0007f
#define MATCH_FCVT_WU_D 0xc2100053
#define MASK_FCVT_WU_D  0xfff0007f
#define MATCH_FCVT_L_D 0xc2200053
#define MASK_FCVT_L_D  0xfff0007f
#define MATCH_FCVT_LU_D 0xc2300053
#define MASK_FCVT_LU_D  0xfff0007f
#define MATCH_FMV_X_D 0xe2000053
#define MASK_FMV_X_D  0xfff0707f
#define MATCH_FCLASS_D 0xe2001053
#define MASK_FCLASS_D  0xfff0707f
#define MATCH_FCVT_W_Q 0xc6000053
#define MASK_FCVT_W_Q  0xfff0007f
#define MATCH_FCVT_WU_Q 0xc6100053
#define MASK_FCVT_WU_Q  0xfff0007f
#define MATCH_FCVT_L_Q 0xc6200053
#define MASK_FCVT_L_Q  0xfff0007f
#define MATCH_FCVT_LU_Q 0xc6300053
#define MASK_FCVT_LU_Q  0xfff0007f
#define MATCH_FMV_X_Q 0xe6000053
#define MASK_FMV_X_Q  0xfff0707f
#define MATCH_FCLASS_Q 0xe6001053
#define MASK_FCLASS_Q  0xfff0707f
#define MATCH_FCVT_S_W 0xd0000053
#define MASK_FCVT_S_W  0xfff0007f
#define MATCH_FCVT_S_WU 0xd0100053
#define MASK_FCVT_S_WU  0xfff0007f
#define MATCH_FCVT_S_L 0xd0200053
#define MASK_FCVT_S_L  0xfff0007f
#define MATCH_FCVT_S_LU 0xd0300053
#define MASK_FCVT_S_LU  0xfff0007f
#define MATCH_FMV_S_X 0xf0000053
#define MASK_FMV_S_X  0xfff0707f
#define MATCH_FCVT_D_W 0xd2000053
#define MASK_FCVT_D_W  0xfff0007f
#define MATCH_FCVT_D_WU 0xd2100053
#define MASK_FCVT_D_WU  0xfff0007f
#define MATCH_FCVT_D_L 0xd2200053
#define MASK_FCVT_D_L  0xfff0007f
#define MATCH_FCVT_D_LU 0xd2300053
#define MASK_FCVT_D_LU  0xfff0007f
#define MATCH_FMV_D_X 0xf2000053
#define MASK_FMV_D_X  0xfff0707f
#define MATCH_FCVT_Q_W 0xd6000053
#define MASK_FCVT_Q_W  0xfff0007f
#define MATCH_FCVT_Q_WU 0xd6100053
#define MASK_FCVT_Q_WU  0xfff0007f
#define MATCH_FCVT_Q_L 0xd6200053
#define MASK_FCVT_Q_L  0xfff0007f
#define MATCH_FCVT_Q_LU 0xd6300053
#define MASK_FCVT_Q_LU  0xfff0007f
#define MATCH_FMV_Q_X 0xf6000053
#define MASK_FMV_Q_X  0xfff0707f
#define MATCH_CLZ 0x60001013
#define MASK_CLZ  0xfff0707f
#define MATCH_CTZ 0x60101013
#define MASK_CTZ  0xfff0707f
#define MATCH_CPOP 0x60201013
#define MASK_CPOP  0xfff0707f
#define MATCH_MIN 0xa004033
#define MASK_MIN  0xfe00707f
#define MATCH_MINU 0xa005033
#define MASK_MINU  0xfe00707f
#define MATCH_MAX 0xa006033
#define MASK_MAX  0xfe00707f
#define MATCH_MAXU 0xa007033
#define MASK_MAXU  0xfe00707f
#define MATCH_SEXT_B 0x60401013
#define MASK_SEXT_B  0xfff0707f
#define MATCH_SEXT_H 0x60501013
#define MASK_SEXT_H  0xfff0707f
#define MATCH_PACK 0x8004033
#define MASK_PACK  0xfe00707f
#define MATCH_PACKW 0x800403b
#define MASK_PACKW  0xfe00707f
#define MATCH_ANDN 0x40007033
#define MASK_ANDN  0xfe00707f
#define MATCH_ORN 0x40006033
#define MASK_ORN  0xfe00707f
#define MATCH_XNOR 0x40004033
#define MASK_XNOR  0xfe00707f
#define MATCH_ROL 0x60001033
#define MASK_ROL  0xfe00707f
#define MATCH_ROR 0x60005033
#define MASK_ROR  0xfe00707f
#define MATCH_RORI 0x60005013
#define MASK_RORI  0xfc00707f
#define MATCH_GREVI 0x68005013
#define MASK_GREVI  0xfc00707f
#define MATCH_GORCI 0x28005013
#define MASK_GORCI  0xfc00707f
#define MATCH_CLZW 0x6000101b
#define MASK_CLZW  0xfff0707f
#define MATCH_CTZW 0x6010101b
#define MASK_CTZW  0xfff0707f
#define MATCH_CPOPW 0x6020101b
#define MASK_CPOPW  0xfff0707f
#define MATCH_ROLW 0x6000103b
#define MASK_ROLW  0xfe00707f
#define MATCH_RORW 0x6000503b
#define MASK_RORW  0xfe00707f
#define MATCH_RORIW 0x6000501b
#define MASK_RORIW  0xfe00707f
#define MATCH_SH1ADD 0x20002033
#define MASK_SH1ADD  0xfe00707f
#define MATCH_SH2ADD 0x20004033
#define MASK_SH2ADD  0xfe00707f
#define MATCH_SH3ADD 0x20006033
#define MASK_SH3ADD  0xfe00707f
#define MATCH_SH1ADD_UW 0x2000203b
#define MASK_SH1ADD_UW  0xfe00707f
#define MATCH_SH2ADD_UW 0x2000403b
#define MASK_SH2ADD_UW  0xfe00707f
#define MATCH_SH3ADD_UW 0x2000603b
#define MASK_SH3ADD_UW  0xfe00707f
#define MATCH_ADD_UW 0x800003b
#define MASK_ADD_UW  0xfe00707f
#define MATCH_SLLI_UW 0x800101b
#define MASK_SLLI_UW  0xfc00707f
#define MATCH_CLMUL 0xa001033
#define MASK_CLMUL  0xfe00707f
#define MATCH_CLMULH 0xa003033
#define MASK_CLMULH  0xfe00707f
#define MATCH_CLMULR 0xa002033
#define MASK_CLMULR  0xfe00707f
#define MATCH_BCLR 0x48001033
#define MASK_BCLR  0xfe00707f
#define MATCH_BCLRI 0x48001013
#define MASK_BCLRI  0xfc00707f
#define MATCH_BEXT 0x48005033
#define MASK_BEXT  0xfe00707f
#define MATCH_BEXTI 0x48005013
#define MASK_BEXTI  0xfc00707f
#define MATCH_BINV 0x68001033
#define MASK_BINV  0xfe00707f
#define MATCH_BINVI 0x68001013
#define MASK_BINVI  0xfc00707f
#define MATCH_BSET 0x28001033
#define MASK_BSET  0xfe00707f
#define MATCH_BSETI 0x28001013
#define MASK_BSETI  0xfc00707f
#define MATCH_FLW 0x2007
#define MASK_FLW  0x707f
#define MATCH_FLD 0x3007
#define MASK_FLD  0x707f
#define MATCH_FLQ 0x4007
#define MASK_FLQ  0x707f
#define MATCH_FSW 0x2027
#define MASK_FSW  0x707f
#define MATCH_FSD 0x3027
#define MASK_FSD  0x707f
#define MATCH_FSQ 0x4027
#define MASK_FSQ  0x707f
#define MATCH_FMADD_S 0x43
#define MASK_FMADD_S  0x600007f
#define MATCH_FMSUB_S 0x47
#define MASK_FMSUB_S  0x600007f
#define MATCH_FNMSUB_S 0x4b
#define MASK_FNMSUB_S  0x600007f
#define MATCH_FNMADD_S 0x4f
#define MASK_FNMADD_S  0x600007f
#define MATCH_FMADD_D 0x2000043
#define MASK_FMADD_D  0x600007f
#define MATCH_FMSUB_D 0x2000047
#define MASK_FMSUB_D  0x600007f
#define MATCH_FNMSUB_D 0x200004b
#define MASK_FNMSUB_D  0x600007f
#define MATCH_FNMADD_D 0x200004f
#define MASK_FNMADD_D  0x600007f
#define MATCH_FMADD_Q 0x6000043
#define MASK_FMADD_Q  0x600007f
#define MATCH_FMSUB_Q 0x6000047
#define MASK_FMSUB_Q  0x600007f
#define MATCH_FNMSUB_Q 0x600004b
#define MASK_FNMSUB_Q  0x600007f
#define MATCH_FNMADD_Q 0x600004f
#define MASK_FNMADD_Q  0x600007f
#define MATCH_C_ADDI4SPN 0x0
#define MASK_C_ADDI4SPN  0xe003
#define MATCH_C_FLD 0x2000
#define MASK_C_FLD  0xe003
#define MATCH_C_LW 0x4000
#define MASK_C_LW  0xe003
#define MATCH_C_FLW 0x6000
#define MASK_C_FLW  0xe003
#define MATCH_C_FSD 0xa000
#define MASK_C_FSD  0xe003
#define MATCH_C_SW 0xc000
#define MASK_C_SW  0xe003
#define MATCH_C_FSW 0xe000
#define MASK_C_FSW  0xe003
#define MATCH_C_ADDI 0x1
#define MASK_C_ADDI  0xe003
#define MATCH_C_JAL 0x2001
#define MASK_C_JAL  0xe003
#define MATCH_C_LI 0x4001
#define MASK_C_LI  0xe003
#define MATCH_C_LUI 0x6001
#define MASK_C_LUI  0xe003
#define MATCH_C_SRLI 0x8001
#define MASK_C_SRLI  0xec03
#define MATCH_C_SRLI64 0x8001
#define MASK_C_SRLI64  0xfc7f
#define MATCH_C_SRAI 0x8401
#define MASK_C_SRAI  0xec03
#define MATCH_C_SRAI64 0x8401
#define MASK_C_SRAI64  0xfc7f
#define MATCH_C_ANDI 0x8801
#define MASK_C_ANDI  0xec03
#define MATCH_C_SUB 0x8c01
#define MASK_C_SUB  0xfc63
#define MATCH_C_XOR 0x8c21
#define MASK_C_XOR  0xfc63
#define MATCH_C_OR 0x8c41
#define MASK_C_OR  0xfc63
#define MATCH_C_AND 0x8c61
#define MASK_C_AND  0xfc63
#define MATCH_C_SUBW 0x9c01
#define MASK_C_SUBW  0xfc63
#define MATCH_C_ADDW 0x9c21
#define MASK_C_ADDW  0xfc63
#define MATCH_C_J 0xa001
#define MASK_C_J  0xe003
#define MATCH_C_BEQZ 0xc001
#define MASK_C_BEQZ  0xe003
#define MATCH_C_BNEZ 0xe001
#define MASK_C_BNEZ  0xe003
#define MATCH_C_SLLI 0x2
#define MASK_C_SLLI  0xe003
#define MATCH_C_SLLI64 0x2
#define MASK_C_SLLI64 0xf07f
#define MATCH_C_FLDSP 0x2002
#define MASK_C_FLDSP  0xe003
#define MATCH_C_LWSP 0x4002
#define MASK_C_LWSP  0xe003
#define MATCH_C_FLWSP 0x6002
#define MASK_C_FLWSP  0xe003
#define MATCH_C_MV 0x8002
#define MASK_C_MV  0xf003
#define MATCH_C_ADD 0x9002
#define MASK_C_ADD  0xf003
#define MATCH_C_FSDSP 0xa002
#define MASK_C_FSDSP  0xe003
#define MATCH_C_SWSP 0xc002
#define MASK_C_SWSP  0xe003
#define MATCH_C_FSWSP 0xe002
#define MASK_C_FSWSP  0xe003
#define MATCH_C_NOP 0x1
#define MASK_C_NOP  0xffff
#define MATCH_C_ADDI16SP 0x6101
#define MASK_C_ADDI16SP  0xef83
#define MATCH_C_JR 0x8002
#define MASK_C_JR  0xf07f
#define MATCH_C_JALR 0x9002
#define MASK_C_JALR  0xf07f
#define MATCH_C_EBREAK 0x9002
#define MASK_C_EBREAK  0xffff
#define MATCH_C_LD 0x6000
#define MASK_C_LD  0xe003
#define MATCH_C_SD 0xe000
#define MASK_C_SD  0xe003
#define MATCH_C_ADDIW 0x2001
#define MASK_C_ADDIW  0xe003
#define MATCH_C_LDSP 0x6002
#define MASK_C_LDSP  0xe003
#define MATCH_C_SDSP 0xe002
#define MASK_C_SDSP  0xe003

/* RVV */
/* Version 1.0-draft-20210130.  */

/* Temporary configuration-setting encoding info

`-` means zimm

31 30 zimm  RS2   RS1/uimm funct3 RD    opcode
1  0  00000 xxxxx xxxxx    111    xxxxx 1010111 vsetvl
1  1  ----- ----- xxxxx    111    xxxxx 1010111 vsetivli
0  -  ----- ----- xxxxx    111    xxxxx 1010111 vsetvli
*/

#define MATCH_VSETVL   0x80007057
#define MASK_VSETVL    0xfe00707f
#define MATCH_VSETIVLI 0xc0007057
#define MASK_VSETIVLI  0xc000707f
#define MATCH_VSETVLI  0x00007057
#define MASK_VSETVLI   0x8000707f

/* Temporary Load/store encoding info

MOP load
00 unit-stride		LE<EEW>, VLE<EEW>FF, VL<nf>RE<EEW> (nf = 1, 2, 4, 8)
01 indexed-unordered	VLUXEI<EEW>
10 strided		VLSE<EEW>
11 indexed-ordered	VLOXEI<EEW>

MOP store
00 unit-stride		VSE<EEW>, VS<nf>R (nf = 1, 2, 4, 8)
01 indexed-unordered	VSUXEI<EEW>
10 strided		VSSE<EEW>
11 indexed-ordered	VSOXEI<EEW>

VM 0 masked
VM 1 unmasked

LUMOP
00000 unit-stride load
01000 unit-stride, whole registers load
01011 unit-stride, mask load, EEW = 1
10000 unit-stride first-fault
xxxxx other encodings reserved, x != 0

SUMOP
00000 unit-stride store
01000 unit-stride, whole registers store
01011 unit-stride, mask store, EEW = 1
0xxxx other encodings reserved, x != 0

`-` means EEW =
MEW WIDTH
x   001   FLH/FSH
x   010   FLW/FSW
x   011   FLD/FSW
x   100   FLQ/FSQ
0   000   VLxE8/VSxE8, VLxEI8/VSxEI8, VL<nf>RE8, VS<nf>R
0   101   VLxE16/VSxE16, VLxEI16/VSxEI16, VL<nf>RE16
0   110   VLxE32/VSxE32, VLxEI32/VSxEI32, VL<nf>RE32
0   111   VLxE64/VSxE64, VLxEI64/VSxEI64, VL<nf>RE64
1   000   Reserved (VLxE128/VSxE128, VL<nf>RE128)
1   101   Reserved (VLxE256/VSxE256, VL<nf>RE256)
1   110   Reserved (VLxE512/VSxE512, VL<nf>RE512)
1   111   Reserved (VLxE1024/VSxE1024, VL<nf>RE1024)

NF  MEW MOP VM LUMOP/RS2 RS1   WIDTH VD    opcode
000 -   00  x  00000     xxxxx ---   xxxxx 0000111 VLE<EEW>
000 -   00  x  00000     xxxxx ---   xxxxx 0100111 VSE<EEW>
000 -   00  1  01011     xxxxx ---   xxxxx 0000111 VLE, EEW = 1
000 -   00  1  01011     xxxxx ---   xxxxx 0100111 VSE, EEW = 1
000 -   10  x  xxxxx     xxxxx ---   xxxxx 0000111 VLSE<EEW>
000 -   10  x  xxxxx     xxxxx ---   xxxxx 0100111 VSSE<EEW>
000 0   11  x  xxxxx     xxxxx ---   xxxxx 0000111 VLOXE<EEW>I
000 0   11  x  xxxxx     xxxxx ---   xxxxx 0100111 VSOXE<EEW>I
000 0   01  x  xxxxx     xxxxx ---   xxxxx 0000111 VLUXE<EEW>I
000 0   01  x  xxxxx     xxxxx ---   xxxxx 0100111 VSUXE<EEW>I
000 -   00  x  10000     xxxxx ---   xxxxx 0000111 VLE<EEW>FF
xxx -   00  1  01000     xxxxx ---   xxxxx 0000111 VL<nf>RE<EEW>, nf = 1,2,4,8
xxx 0   00  1  01000     xxxxx 000   xxxxx 0100111 VS<nf>R, nf = 1,2,4,8

xxx -   00  x  00000     xxxxx ---   xxxxx 0000111 VLSEG<nf>E<EEW>
xxx -   00  x  00000     xxxxx ---   xxxxx 0100111 VSSEG<nf>E<EEW>
xxx -   10  x  00000     xxxxx ---   xxxxx 0000111 VLSSEG<nf>E<EEW>
xxx -   10  x  00000     xxxxx ---   xxxxx 0100111 VSSSEG<nf>E<EEW>
xxx -   11  x  00000     xxxxx ---   xxxxx 0000111 VLOXSEG<nf>E<EEW>I
xxx -   11  x  00000     xxxxx ---   xxxxx 0100111 VSOXSEG<nf>E<EEW>I
xxx -   01  x  00000     xxxxx ---   xxxxx 0000111 VLUXSEG<nf>E<EEW>I
xxx -   01  x  00000     xxxxx ---   xxxxx 0100111 VSUXSEG<nf>E<EEW>I
xxx -   00  x  10000     xxxxx ---   xxxxx 0000111 VLSEG<nf>E<EEW>FF
*/

#define MATCH_VLMV     0x02b00007
#define MASK_VLMV      0xfff0707f
#define MATCH_VSMV     0x02b00027
#define MASK_VSMV      0xfff0707f

#define MATCH_VLE8V    0x00000007
#define MASK_VLE8V     0xfdf0707f
#define MATCH_VLE16V   0x00005007
#define MASK_VLE16V    0xfdf0707f
#define MATCH_VLE32V   0x00006007
#define MASK_VLE32V    0xfdf0707f
#define MATCH_VLE64V   0x00007007
#define MASK_VLE64V    0xfdf0707f

#define MATCH_VSE8V    0x00000027
#define MASK_VSE8V     0xfdf0707f
#define MATCH_VSE16V   0x00005027
#define MASK_VSE16V    0xfdf0707f
#define MATCH_VSE32V   0x00006027
#define MASK_VSE32V    0xfdf0707f
#define MATCH_VSE64V   0x00007027
#define MASK_VSE64V    0xfdf0707f

#define MATCH_VLSE8V    0x08000007
#define MASK_VLSE8V     0xfc00707f
#define MATCH_VLSE16V   0x08005007
#define MASK_VLSE16V    0xfc00707f
#define MATCH_VLSE32V   0x08006007
#define MASK_VLSE32V    0xfc00707f
#define MATCH_VLSE64V   0x08007007
#define MASK_VLSE64V    0xfc00707f

#define MATCH_VSSE8V    0x08000027
#define MASK_VSSE8V     0xfc00707f
#define MATCH_VSSE16V   0x08005027
#define MASK_VSSE16V    0xfc00707f
#define MATCH_VSSE32V   0x08006027
#define MASK_VSSE32V    0xfc00707f
#define MATCH_VSSE64V   0x08007027
#define MASK_VSSE64V    0xfc00707f

#define MATCH_VLOXEI8V    0x0c000007
#define MASK_VLOXEI8V     0xfc00707f
#define MATCH_VLOXEI16V   0x0c005007
#define MASK_VLOXEI16V    0xfc00707f
#define MATCH_VLOXEI32V   0x0c006007
#define MASK_VLOXEI32V    0xfc00707f
#define MATCH_VLOXEI64V   0x0c007007
#define MASK_VLOXEI64V    0xfc00707f

#define MATCH_VSOXEI8V    0x0c000027
#define MASK_VSOXEI8V     0xfc00707f
#define MATCH_VSOXEI16V   0x0c005027
#define MASK_VSOXEI16V    0xfc00707f
#define MATCH_VSOXEI32V   0x0c006027
#define MASK_VSOXEI32V    0xfc00707f
#define MATCH_VSOXEI64V   0x0c007027
#define MASK_VSOXEI64V    0xfc00707f

#define MATCH_VLUXEI8V    0x04000007
#define MASK_VLUXEI8V     0xfc00707f
#define MATCH_VLUXEI16V   0x04005007
#define MASK_VLUXEI16V    0xfc00707f
#define MATCH_VLUXEI32V   0x04006007
#define MASK_VLUXEI32V    0xfc00707f
#define MATCH_VLUXEI64V   0x04007007
#define MASK_VLUXEI64V    0xfc00707f

#define MATCH_VSUXEI8V    0x04000027
#define MASK_VSUXEI8V     0xfc00707f
#define MATCH_VSUXEI16V   0x04005027
#define MASK_VSUXEI16V    0xfc00707f
#define MATCH_VSUXEI32V   0x04006027
#define MASK_VSUXEI32V    0xfc00707f
#define MATCH_VSUXEI64V   0x04007027
#define MASK_VSUXEI64V    0xfc00707f

#define MATCH_VLE8FFV    0x01000007
#define MASK_VLE8FFV     0xfdf0707f
#define MATCH_VLE16FFV   0x01005007
#define MASK_VLE16FFV    0xfdf0707f
#define MATCH_VLE32FFV   0x01006007
#define MASK_VLE32FFV    0xfdf0707f
#define MATCH_VLE64FFV   0x01007007
#define MASK_VLE64FFV    0xfdf0707f

#define MATCH_VLSEG2E8V  0x20000007
#define MASK_VLSEG2E8V   0xfdf0707f
#define MATCH_VSSEG2E8V  0x20000027
#define MASK_VSSEG2E8V   0xfdf0707f
#define MATCH_VLSEG3E8V  0x40000007
#define MASK_VLSEG3E8V   0xfdf0707f
#define MATCH_VSSEG3E8V  0x40000027
#define MASK_VSSEG3E8V   0xfdf0707f
#define MATCH_VLSEG4E8V  0x60000007
#define MASK_VLSEG4E8V   0xfdf0707f
#define MATCH_VSSEG4E8V  0x60000027
#define MASK_VSSEG4E8V   0xfdf0707f
#define MATCH_VLSEG5E8V  0x80000007
#define MASK_VLSEG5E8V   0xfdf0707f
#define MATCH_VSSEG5E8V  0x80000027
#define MASK_VSSEG5E8V   0xfdf0707f
#define MATCH_VLSEG6E8V  0xa0000007
#define MASK_VLSEG6E8V   0xfdf0707f
#define MATCH_VSSEG6E8V  0xa0000027
#define MASK_VSSEG6E8V   0xfdf0707f
#define MATCH_VLSEG7E8V  0xc0000007
#define MASK_VLSEG7E8V   0xfdf0707f
#define MATCH_VSSEG7E8V  0xc0000027
#define MASK_VSSEG7E8V   0xfdf0707f
#define MATCH_VLSEG8E8V  0xe0000007
#define MASK_VLSEG8E8V   0xfdf0707f
#define MATCH_VSSEG8E8V  0xe0000027
#define MASK_VSSEG8E8V   0xfdf0707f

#define MATCH_VLSEG2E16V  0x20005007
#define MASK_VLSEG2E16V   0xfdf0707f
#define MATCH_VSSEG2E16V  0x20005027
#define MASK_VSSEG2E16V   0xfdf0707f
#define MATCH_VLSEG3E16V  0x40005007
#define MASK_VLSEG3E16V   0xfdf0707f
#define MATCH_VSSEG3E16V  0x40005027
#define MASK_VSSEG3E16V   0xfdf0707f
#define MATCH_VLSEG4E16V  0x60005007
#define MASK_VLSEG4E16V   0xfdf0707f
#define MATCH_VSSEG4E16V  0x60005027
#define MASK_VSSEG4E16V   0xfdf0707f
#define MATCH_VLSEG5E16V  0x80005007
#define MASK_VLSEG5E16V   0xfdf0707f
#define MATCH_VSSEG5E16V  0x80005027
#define MASK_VSSEG5E16V   0xfdf0707f
#define MATCH_VLSEG6E16V  0xa0005007
#define MASK_VLSEG6E16V   0xfdf0707f
#define MATCH_VSSEG6E16V  0xa0005027
#define MASK_VSSEG6E16V   0xfdf0707f
#define MATCH_VLSEG7E16V  0xc0005007
#define MASK_VLSEG7E16V   0xfdf0707f
#define MATCH_VSSEG7E16V  0xc0005027
#define MASK_VSSEG7E16V   0xfdf0707f
#define MATCH_VLSEG8E16V  0xe0005007
#define MASK_VLSEG8E16V   0xfdf0707f
#define MATCH_VSSEG8E16V  0xe0005027
#define MASK_VSSEG8E16V   0xfdf0707f

#define MATCH_VLSEG2E32V  0x20006007
#define MASK_VLSEG2E32V   0xfdf0707f
#define MATCH_VSSEG2E32V  0x20006027
#define MASK_VSSEG2E32V   0xfdf0707f
#define MATCH_VLSEG3E32V  0x40006007
#define MASK_VLSEG3E32V   0xfdf0707f
#define MATCH_VSSEG3E32V  0x40006027
#define MASK_VSSEG3E32V   0xfdf0707f
#define MATCH_VLSEG4E32V  0x60006007
#define MASK_VLSEG4E32V   0xfdf0707f
#define MATCH_VSSEG4E32V  0x60006027
#define MASK_VSSEG4E32V   0xfdf0707f
#define MATCH_VLSEG5E32V  0x80006007
#define MASK_VLSEG5E32V   0xfdf0707f
#define MATCH_VSSEG5E32V  0x80006027
#define MASK_VSSEG5E32V   0xfdf0707f
#define MATCH_VLSEG6E32V  0xa0006007
#define MASK_VLSEG6E32V   0xfdf0707f
#define MATCH_VSSEG6E32V  0xa0006027
#define MASK_VSSEG6E32V   0xfdf0707f
#define MATCH_VLSEG7E32V  0xc0006007
#define MASK_VLSEG7E32V   0xfdf0707f
#define MATCH_VSSEG7E32V  0xc0006027
#define MASK_VSSEG7E32V   0xfdf0707f
#define MATCH_VLSEG8E32V  0xe0006007
#define MASK_VLSEG8E32V   0xfdf0707f
#define MATCH_VSSEG8E32V  0xe0006027
#define MASK_VSSEG8E32V   0xfdf0707f

#define MATCH_VLSEG2E64V  0x20007007
#define MASK_VLSEG2E64V   0xfdf0707f
#define MATCH_VSSEG2E64V  0x20007027
#define MASK_VSSEG2E64V   0xfdf0707f
#define MATCH_VLSEG3E64V  0x40007007
#define MASK_VLSEG3E64V   0xfdf0707f
#define MATCH_VSSEG3E64V  0x40007027
#define MASK_VSSEG3E64V   0xfdf0707f
#define MATCH_VLSEG4E64V  0x60007007
#define MASK_VLSEG4E64V   0xfdf0707f
#define MATCH_VSSEG4E64V  0x60007027
#define MASK_VSSEG4E64V   0xfdf0707f
#define MATCH_VLSEG5E64V  0x80007007
#define MASK_VLSEG5E64V   0xfdf0707f
#define MATCH_VSSEG5E64V  0x80007027
#define MASK_VSSEG5E64V   0xfdf0707f
#define MATCH_VLSEG6E64V  0xa0007007
#define MASK_VLSEG6E64V   0xfdf0707f
#define MATCH_VSSEG6E64V  0xa0007027
#define MASK_VSSEG6E64V   0xfdf0707f
#define MATCH_VLSEG7E64V  0xc0007007
#define MASK_VLSEG7E64V   0xfdf0707f
#define MATCH_VSSEG7E64V  0xc0007027
#define MASK_VSSEG7E64V   0xfdf0707f
#define MATCH_VLSEG8E64V  0xe0007007
#define MASK_VLSEG8E64V   0xfdf0707f
#define MATCH_VSSEG8E64V  0xe0007027
#define MASK_VSSEG8E64V   0xfdf0707f

#define MATCH_VLSSEG2E8V  0x28000007
#define MASK_VLSSEG2E8V   0xfc00707f
#define MATCH_VSSSEG2E8V  0x28000027
#define MASK_VSSSEG2E8V   0xfc00707f
#define MATCH_VLSSEG3E8V  0x48000007
#define MASK_VLSSEG3E8V   0xfc00707f
#define MATCH_VSSSEG3E8V  0x48000027
#define MASK_VSSSEG3E8V   0xfc00707f
#define MATCH_VLSSEG4E8V  0x68000007
#define MASK_VLSSEG4E8V   0xfc00707f
#define MATCH_VSSSEG4E8V  0x68000027
#define MASK_VSSSEG4E8V   0xfc00707f
#define MATCH_VLSSEG5E8V  0x88000007
#define MASK_VLSSEG5E8V   0xfc00707f
#define MATCH_VSSSEG5E8V  0x88000027
#define MASK_VSSSEG5E8V   0xfc00707f
#define MATCH_VLSSEG6E8V  0xa8000007
#define MASK_VLSSEG6E8V   0xfc00707f
#define MATCH_VSSSEG6E8V  0xa8000027
#define MASK_VSSSEG6E8V   0xfc00707f
#define MATCH_VLSSEG7E8V  0xc8000007
#define MASK_VLSSEG7E8V   0xfc00707f
#define MATCH_VSSSEG7E8V  0xc8000027
#define MASK_VSSSEG7E8V   0xfc00707f
#define MATCH_VLSSEG8E8V  0xe8000007
#define MASK_VLSSEG8E8V   0xfc00707f
#define MATCH_VSSSEG8E8V  0xe8000027
#define MASK_VSSSEG8E8V   0xfc00707f

#define MATCH_VLSSEG2E16V  0x28005007
#define MASK_VLSSEG2E16V   0xfc00707f
#define MATCH_VSSSEG2E16V  0x28005027
#define MASK_VSSSEG2E16V   0xfc00707f
#define MATCH_VLSSEG3E16V  0x48005007
#define MASK_VLSSEG3E16V   0xfc00707f
#define MATCH_VSSSEG3E16V  0x48005027
#define MASK_VSSSEG3E16V   0xfc00707f
#define MATCH_VLSSEG4E16V  0x68005007
#define MASK_VLSSEG4E16V   0xfc00707f
#define MATCH_VSSSEG4E16V  0x68005027
#define MASK_VSSSEG4E16V   0xfc00707f
#define MATCH_VLSSEG5E16V  0x88005007
#define MASK_VLSSEG5E16V   0xfc00707f
#define MATCH_VSSSEG5E16V  0x88005027
#define MASK_VSSSEG5E16V   0xfc00707f
#define MATCH_VLSSEG6E16V  0xa8005007
#define MASK_VLSSEG6E16V   0xfc00707f
#define MATCH_VSSSEG6E16V  0xa8005027
#define MASK_VSSSEG6E16V   0xfc00707f
#define MATCH_VLSSEG7E16V  0xc8005007
#define MASK_VLSSEG7E16V   0xfc00707f
#define MATCH_VSSSEG7E16V  0xc8005027
#define MASK_VSSSEG7E16V   0xfc00707f
#define MATCH_VLSSEG8E16V  0xe8005007
#define MASK_VLSSEG8E16V   0xfc00707f
#define MATCH_VSSSEG8E16V  0xe8005027
#define MASK_VSSSEG8E16V   0xfc00707f

#define MATCH_VLSSEG2E32V  0x28006007
#define MASK_VLSSEG2E32V   0xfc00707f
#define MATCH_VSSSEG2E32V  0x28006027
#define MASK_VSSSEG2E32V   0xfc00707f
#define MATCH_VLSSEG3E32V  0x48006007
#define MASK_VLSSEG3E32V   0xfc00707f
#define MATCH_VSSSEG3E32V  0x48006027
#define MASK_VSSSEG3E32V   0xfc00707f
#define MATCH_VLSSEG4E32V  0x68006007
#define MASK_VLSSEG4E32V   0xfc00707f
#define MATCH_VSSSEG4E32V  0x68006027
#define MASK_VSSSEG4E32V   0xfc00707f
#define MATCH_VLSSEG5E32V  0x88006007
#define MASK_VLSSEG5E32V   0xfc00707f
#define MATCH_VSSSEG5E32V  0x88006027
#define MASK_VSSSEG5E32V   0xfc00707f
#define MATCH_VLSSEG6E32V  0xa8006007
#define MASK_VLSSEG6E32V   0xfc00707f
#define MATCH_VSSSEG6E32V  0xa8006027
#define MASK_VSSSEG6E32V   0xfc00707f
#define MATCH_VLSSEG7E32V  0xc8006007
#define MASK_VLSSEG7E32V   0xfc00707f
#define MATCH_VSSSEG7E32V  0xc8006027
#define MASK_VSSSEG7E32V   0xfc00707f
#define MATCH_VLSSEG8E32V  0xe8006007
#define MASK_VLSSEG8E32V   0xfc00707f
#define MATCH_VSSSEG8E32V  0xe8006027
#define MASK_VSSSEG8E32V   0xfc00707f

#define MATCH_VLSSEG2E64V  0x28007007
#define MASK_VLSSEG2E64V   0xfc00707f
#define MATCH_VSSSEG2E64V  0x28007027
#define MASK_VSSSEG2E64V   0xfc00707f
#define MATCH_VLSSEG3E64V  0x48007007
#define MASK_VLSSEG3E64V   0xfc00707f
#define MATCH_VSSSEG3E64V  0x48007027
#define MASK_VSSSEG3E64V   0xfc00707f
#define MATCH_VLSSEG4E64V  0x68007007
#define MASK_VLSSEG4E64V   0xfc00707f
#define MATCH_VSSSEG4E64V  0x68007027
#define MASK_VSSSEG4E64V   0xfc00707f
#define MATCH_VLSSEG5E64V  0x88007007
#define MASK_VLSSEG5E64V   0xfc00707f
#define MATCH_VSSSEG5E64V  0x88007027
#define MASK_VSSSEG5E64V   0xfc00707f
#define MATCH_VLSSEG6E64V  0xa8007007
#define MASK_VLSSEG6E64V   0xfc00707f
#define MATCH_VSSSEG6E64V  0xa8007027
#define MASK_VSSSEG6E64V   0xfc00707f
#define MATCH_VLSSEG7E64V  0xc8007007
#define MASK_VLSSEG7E64V   0xfc00707f
#define MATCH_VSSSEG7E64V  0xc8007027
#define MASK_VSSSEG7E64V   0xfc00707f
#define MATCH_VLSSEG8E64V  0xe8007007
#define MASK_VLSSEG8E64V   0xfc00707f
#define MATCH_VSSSEG8E64V  0xe8007027
#define MASK_VSSSEG8E64V   0xfc00707f

#define MATCH_VLOXSEG2EI8V  0x2c000007
#define MASK_VLOXSEG2EI8V   0xfc00707f
#define MATCH_VSOXSEG2EI8V  0x2c000027
#define MASK_VSOXSEG2EI8V   0xfc00707f
#define MATCH_VLOXSEG3EI8V  0x4c000007
#define MASK_VLOXSEG3EI8V   0xfc00707f
#define MATCH_VSOXSEG3EI8V  0x4c000027
#define MASK_VSOXSEG3EI8V   0xfc00707f
#define MATCH_VLOXSEG4EI8V  0x6c000007
#define MASK_VLOXSEG4EI8V   0xfc00707f
#define MATCH_VSOXSEG4EI8V  0x6c000027
#define MASK_VSOXSEG4EI8V   0xfc00707f
#define MATCH_VLOXSEG5EI8V  0x8c000007
#define MASK_VLOXSEG5EI8V   0xfc00707f
#define MATCH_VSOXSEG5EI8V  0x8c000027
#define MASK_VSOXSEG5EI8V   0xfc00707f
#define MATCH_VLOXSEG6EI8V  0xac000007
#define MASK_VLOXSEG6EI8V   0xfc00707f
#define MATCH_VSOXSEG6EI8V  0xac000027
#define MASK_VSOXSEG6EI8V   0xfc00707f
#define MATCH_VLOXSEG7EI8V  0xcc000007
#define MASK_VLOXSEG7EI8V   0xfc00707f
#define MATCH_VSOXSEG7EI8V  0xcc000027
#define MASK_VSOXSEG7EI8V   0xfc00707f
#define MATCH_VLOXSEG8EI8V  0xec000007
#define MASK_VLOXSEG8EI8V   0xfc00707f
#define MATCH_VSOXSEG8EI8V  0xec000027
#define MASK_VSOXSEG8EI8V   0xfc00707f

#define MATCH_VLUXSEG2EI8V  0x24000007
#define MASK_VLUXSEG2EI8V   0xfc00707f
#define MATCH_VSUXSEG2EI8V  0x24000027
#define MASK_VSUXSEG2EI8V   0xfc00707f
#define MATCH_VLUXSEG3EI8V  0x44000007
#define MASK_VLUXSEG3EI8V   0xfc00707f
#define MATCH_VSUXSEG3EI8V  0x44000027
#define MASK_VSUXSEG3EI8V   0xfc00707f
#define MATCH_VLUXSEG4EI8V  0x64000007
#define MASK_VLUXSEG4EI8V   0xfc00707f
#define MATCH_VSUXSEG4EI8V  0x64000027
#define MASK_VSUXSEG4EI8V   0xfc00707f
#define MATCH_VLUXSEG5EI8V  0x84000007
#define MASK_VLUXSEG5EI8V   0xfc00707f
#define MATCH_VSUXSEG5EI8V  0x84000027
#define MASK_VSUXSEG5EI8V   0xfc00707f
#define MATCH_VLUXSEG6EI8V  0xa4000007
#define MASK_VLUXSEG6EI8V   0xfc00707f
#define MATCH_VSUXSEG6EI8V  0xa4000027
#define MASK_VSUXSEG6EI8V   0xfc00707f
#define MATCH_VLUXSEG7EI8V  0xc4000007
#define MASK_VLUXSEG7EI8V   0xfc00707f
#define MATCH_VSUXSEG7EI8V  0xc4000027
#define MASK_VSUXSEG7EI8V   0xfc00707f
#define MATCH_VLUXSEG8EI8V  0xe4000007
#define MASK_VLUXSEG8EI8V   0xfc00707f
#define MATCH_VSUXSEG8EI8V  0xe4000027
#define MASK_VSUXSEG8EI8V   0xfc00707f

#define MATCH_VLOXSEG2EI16V  0x2c005007
#define MASK_VLOXSEG2EI16V   0xfc00707f
#define MATCH_VSOXSEG2EI16V  0x2c005027
#define MASK_VSOXSEG2EI16V   0xfc00707f
#define MATCH_VLOXSEG3EI16V  0x4c005007
#define MASK_VLOXSEG3EI16V   0xfc00707f
#define MATCH_VSOXSEG3EI16V  0x4c005027
#define MASK_VSOXSEG3EI16V   0xfc00707f
#define MATCH_VLOXSEG4EI16V  0x6c005007
#define MASK_VLOXSEG4EI16V   0xfc00707f
#define MATCH_VSOXSEG4EI16V  0x6c005027
#define MASK_VSOXSEG4EI16V   0xfc00707f
#define MATCH_VLOXSEG5EI16V  0x8c005007
#define MASK_VLOXSEG5EI16V   0xfc00707f
#define MATCH_VSOXSEG5EI16V  0x8c005027
#define MASK_VSOXSEG5EI16V   0xfc00707f
#define MATCH_VLOXSEG6EI16V  0xac005007
#define MASK_VLOXSEG6EI16V   0xfc00707f
#define MATCH_VSOXSEG6EI16V  0xac005027
#define MASK_VSOXSEG6EI16V   0xfc00707f
#define MATCH_VLOXSEG7EI16V  0xcc005007
#define MASK_VLOXSEG7EI16V   0xfc00707f
#define MATCH_VSOXSEG7EI16V  0xcc005027
#define MASK_VSOXSEG7EI16V   0xfc00707f
#define MATCH_VLOXSEG8EI16V  0xec005007
#define MASK_VLOXSEG8EI16V   0xfc00707f
#define MATCH_VSOXSEG8EI16V  0xec005027
#define MASK_VSOXSEG8EI16V   0xfc00707f

#define MATCH_VLUXSEG2EI16V  0x24005007
#define MASK_VLUXSEG2EI16V   0xfc00707f
#define MATCH_VSUXSEG2EI16V  0x24005027
#define MASK_VSUXSEG2EI16V   0xfc00707f
#define MATCH_VLUXSEG3EI16V  0x44005007
#define MASK_VLUXSEG3EI16V   0xfc00707f
#define MATCH_VSUXSEG3EI16V  0x44005027
#define MASK_VSUXSEG3EI16V   0xfc00707f
#define MATCH_VLUXSEG4EI16V  0x64005007
#define MASK_VLUXSEG4EI16V   0xfc00707f
#define MATCH_VSUXSEG4EI16V  0x64005027
#define MASK_VSUXSEG4EI16V   0xfc00707f
#define MATCH_VLUXSEG5EI16V  0x84005007
#define MASK_VLUXSEG5EI16V   0xfc00707f
#define MATCH_VSUXSEG5EI16V  0x84005027
#define MASK_VSUXSEG5EI16V   0xfc00707f
#define MATCH_VLUXSEG6EI16V  0xa4005007
#define MASK_VLUXSEG6EI16V   0xfc00707f
#define MATCH_VSUXSEG6EI16V  0xa4005027
#define MASK_VSUXSEG6EI16V   0xfc00707f
#define MATCH_VLUXSEG7EI16V  0xc4005007
#define MASK_VLUXSEG7EI16V   0xfc00707f
#define MATCH_VSUXSEG7EI16V  0xc4005027
#define MASK_VSUXSEG7EI16V   0xfc00707f
#define MATCH_VLUXSEG8EI16V  0xe4005007
#define MASK_VLUXSEG8EI16V   0xfc00707f
#define MATCH_VSUXSEG8EI16V  0xe4005027
#define MASK_VSUXSEG8EI16V   0xfc00707f

#define MATCH_VLOXSEG2EI32V  0x2c006007
#define MASK_VLOXSEG2EI32V   0xfc00707f
#define MATCH_VSOXSEG2EI32V  0x2c006027
#define MASK_VSOXSEG2EI32V   0xfc00707f
#define MATCH_VLOXSEG3EI32V  0x4c006007
#define MASK_VLOXSEG3EI32V   0xfc00707f
#define MATCH_VSOXSEG3EI32V  0x4c006027
#define MASK_VSOXSEG3EI32V   0xfc00707f
#define MATCH_VLOXSEG4EI32V  0x6c006007
#define MASK_VLOXSEG4EI32V   0xfc00707f
#define MATCH_VSOXSEG4EI32V  0x6c006027
#define MASK_VSOXSEG4EI32V   0xfc00707f
#define MATCH_VLOXSEG5EI32V  0x8c006007
#define MASK_VLOXSEG5EI32V   0xfc00707f
#define MATCH_VSOXSEG5EI32V  0x8c006027
#define MASK_VSOXSEG5EI32V   0xfc00707f
#define MATCH_VLOXSEG6EI32V  0xac006007
#define MASK_VLOXSEG6EI32V   0xfc00707f
#define MATCH_VSOXSEG6EI32V  0xac006027
#define MASK_VSOXSEG6EI32V   0xfc00707f
#define MATCH_VLOXSEG7EI32V  0xcc006007
#define MASK_VLOXSEG7EI32V   0xfc00707f
#define MATCH_VSOXSEG7EI32V  0xcc006027
#define MASK_VSOXSEG7EI32V   0xfc00707f
#define MATCH_VLOXSEG8EI32V  0xec006007
#define MASK_VLOXSEG8EI32V   0xfc00707f
#define MATCH_VSOXSEG8EI32V  0xec006027
#define MASK_VSOXSEG8EI32V   0xfc00707f

#define MATCH_VLUXSEG2EI32V  0x24006007
#define MASK_VLUXSEG2EI32V   0xfc00707f
#define MATCH_VSUXSEG2EI32V  0x24006027
#define MASK_VSUXSEG2EI32V   0xfc00707f
#define MATCH_VLUXSEG3EI32V  0x44006007
#define MASK_VLUXSEG3EI32V   0xfc00707f
#define MATCH_VSUXSEG3EI32V  0x44006027
#define MASK_VSUXSEG3EI32V   0xfc00707f
#define MATCH_VLUXSEG4EI32V  0x64006007
#define MASK_VLUXSEG4EI32V   0xfc00707f
#define MATCH_VSUXSEG4EI32V  0x64006027
#define MASK_VSUXSEG4EI32V   0xfc00707f
#define MATCH_VLUXSEG5EI32V  0x84006007
#define MASK_VLUXSEG5EI32V   0xfc00707f
#define MATCH_VSUXSEG5EI32V  0x84006027
#define MASK_VSUXSEG5EI32V   0xfc00707f
#define MATCH_VLUXSEG6EI32V  0xa4006007
#define MASK_VLUXSEG6EI32V   0xfc00707f
#define MATCH_VSUXSEG6EI32V  0xa4006027
#define MASK_VSUXSEG6EI32V   0xfc00707f
#define MATCH_VLUXSEG7EI32V  0xc4006007
#define MASK_VLUXSEG7EI32V   0xfc00707f
#define MATCH_VSUXSEG7EI32V  0xc4006027
#define MASK_VSUXSEG7EI32V   0xfc00707f
#define MATCH_VLUXSEG8EI32V  0xe4006007
#define MASK_VLUXSEG8EI32V   0xfc00707f
#define MATCH_VSUXSEG8EI32V  0xe4006027
#define MASK_VSUXSEG8EI32V   0xfc00707f

#define MATCH_VLOXSEG2EI64V  0x2c007007
#define MASK_VLOXSEG2EI64V   0xfc00707f
#define MATCH_VSOXSEG2EI64V  0x2c007027
#define MASK_VSOXSEG2EI64V   0xfc00707f
#define MATCH_VLOXSEG3EI64V  0x4c007007
#define MASK_VLOXSEG3EI64V   0xfc00707f
#define MATCH_VSOXSEG3EI64V  0x4c007027
#define MASK_VSOXSEG3EI64V   0xfc00707f
#define MATCH_VLOXSEG4EI64V  0x6c007007
#define MASK_VLOXSEG4EI64V   0xfc00707f
#define MATCH_VSOXSEG4EI64V  0x6c007027
#define MASK_VSOXSEG4EI64V   0xfc00707f
#define MATCH_VLOXSEG5EI64V  0x8c007007
#define MASK_VLOXSEG5EI64V   0xfc00707f
#define MATCH_VSOXSEG5EI64V  0x8c007027
#define MASK_VSOXSEG5EI64V   0xfc00707f
#define MATCH_VLOXSEG6EI64V  0xac007007
#define MASK_VLOXSEG6EI64V   0xfc00707f
#define MATCH_VSOXSEG6EI64V  0xac007027
#define MASK_VSOXSEG6EI64V   0xfc00707f
#define MATCH_VLOXSEG7EI64V  0xcc007007
#define MASK_VLOXSEG7EI64V   0xfc00707f
#define MATCH_VSOXSEG7EI64V  0xcc007027
#define MASK_VSOXSEG7EI64V   0xfc00707f
#define MATCH_VLOXSEG8EI64V  0xec007007
#define MASK_VLOXSEG8EI64V   0xfc00707f
#define MATCH_VSOXSEG8EI64V  0xec007027
#define MASK_VSOXSEG8EI64V   0xfc00707f

#define MATCH_VLUXSEG2EI64V  0x24007007
#define MASK_VLUXSEG2EI64V   0xfc00707f
#define MATCH_VSUXSEG2EI64V  0x24007027
#define MASK_VSUXSEG2EI64V   0xfc00707f
#define MATCH_VLUXSEG3EI64V  0x44007007
#define MASK_VLUXSEG3EI64V   0xfc00707f
#define MATCH_VSUXSEG3EI64V  0x44007027
#define MASK_VSUXSEG3EI64V   0xfc00707f
#define MATCH_VLUXSEG4EI64V  0x64007007
#define MASK_VLUXSEG4EI64V   0xfc00707f
#define MATCH_VSUXSEG4EI64V  0x64007027
#define MASK_VSUXSEG4EI64V   0xfc00707f
#define MATCH_VLUXSEG5EI64V  0x84007007
#define MASK_VLUXSEG5EI64V   0xfc00707f
#define MATCH_VSUXSEG5EI64V  0x84007027
#define MASK_VSUXSEG5EI64V   0xfc00707f
#define MATCH_VLUXSEG6EI64V  0xa4007007
#define MASK_VLUXSEG6EI64V   0xfc00707f
#define MATCH_VSUXSEG6EI64V  0xa4007027
#define MASK_VSUXSEG6EI64V   0xfc00707f
#define MATCH_VLUXSEG7EI64V  0xc4007007
#define MASK_VLUXSEG7EI64V   0xfc00707f
#define MATCH_VSUXSEG7EI64V  0xc4007027
#define MASK_VSUXSEG7EI64V   0xfc00707f
#define MATCH_VLUXSEG8EI64V  0xe4007007
#define MASK_VLUXSEG8EI64V   0xfc00707f
#define MATCH_VSUXSEG8EI64V  0xe4007027
#define MASK_VSUXSEG8EI64V   0xfc00707f

#define MATCH_VLSEG2E8FFV  0x21000007
#define MASK_VLSEG2E8FFV   0xfdf0707f
#define MATCH_VLSEG3E8FFV  0x41000007
#define MASK_VLSEG3E8FFV   0xfdf0707f
#define MATCH_VLSEG4E8FFV  0x61000007
#define MASK_VLSEG4E8FFV   0xfdf0707f
#define MATCH_VLSEG5E8FFV  0x81000007
#define MASK_VLSEG5E8FFV   0xfdf0707f
#define MATCH_VLSEG6E8FFV  0xa1000007
#define MASK_VLSEG6E8FFV   0xfdf0707f
#define MATCH_VLSEG7E8FFV  0xc1000007
#define MASK_VLSEG7E8FFV   0xfdf0707f
#define MATCH_VLSEG8E8FFV  0xe1000007
#define MASK_VLSEG8E8FFV   0xfdf0707f

#define MATCH_VLSEG2E16FFV  0x21005007
#define MASK_VLSEG2E16FFV   0xfdf0707f
#define MATCH_VLSEG3E16FFV  0x41005007
#define MASK_VLSEG3E16FFV   0xfdf0707f
#define MATCH_VLSEG4E16FFV  0x61005007
#define MASK_VLSEG4E16FFV   0xfdf0707f
#define MATCH_VLSEG5E16FFV  0x81005007
#define MASK_VLSEG5E16FFV   0xfdf0707f
#define MATCH_VLSEG6E16FFV  0xa1005007
#define MASK_VLSEG6E16FFV   0xfdf0707f
#define MATCH_VLSEG7E16FFV  0xc1005007
#define MASK_VLSEG7E16FFV   0xfdf0707f
#define MATCH_VLSEG8E16FFV  0xe1005007
#define MASK_VLSEG8E16FFV   0xfdf0707f

#define MATCH_VLSEG2E32FFV  0x21006007
#define MASK_VLSEG2E32FFV   0xfdf0707f
#define MATCH_VLSEG3E32FFV  0x41006007
#define MASK_VLSEG3E32FFV   0xfdf0707f
#define MATCH_VLSEG4E32FFV  0x61006007
#define MASK_VLSEG4E32FFV   0xfdf0707f
#define MATCH_VLSEG5E32FFV  0x81006007
#define MASK_VLSEG5E32FFV   0xfdf0707f
#define MATCH_VLSEG6E32FFV  0xa1006007
#define MASK_VLSEG6E32FFV   0xfdf0707f
#define MATCH_VLSEG7E32FFV  0xc1006007
#define MASK_VLSEG7E32FFV   0xfdf0707f
#define MATCH_VLSEG8E32FFV  0xe1006007
#define MASK_VLSEG8E32FFV   0xfdf0707f

#define MATCH_VLSEG2E64FFV  0x21007007
#define MASK_VLSEG2E64FFV   0xfdf0707f
#define MATCH_VLSEG3E64FFV  0x41007007
#define MASK_VLSEG3E64FFV   0xfdf0707f
#define MATCH_VLSEG4E64FFV  0x61007007
#define MASK_VLSEG4E64FFV   0xfdf0707f
#define MATCH_VLSEG5E64FFV  0x81007007
#define MASK_VLSEG5E64FFV   0xfdf0707f
#define MATCH_VLSEG6E64FFV  0xa1007007
#define MASK_VLSEG6E64FFV   0xfdf0707f
#define MATCH_VLSEG7E64FFV  0xc1007007
#define MASK_VLSEG7E64FFV   0xfdf0707f
#define MATCH_VLSEG8E64FFV  0xe1007007
#define MASK_VLSEG8E64FFV   0xfdf0707f

#define MATCH_VL1RE8V    0x02800007
#define MASK_VL1RE8V     0xfff0707f
#define MATCH_VL1RE16V   0x02805007
#define MASK_VL1RE16V    0xfff0707f
#define MATCH_VL1RE32V   0x02806007
#define MASK_VL1RE32V    0xfff0707f
#define MATCH_VL1RE64V   0x02807007
#define MASK_VL1RE64V    0xfff0707f

#define MATCH_VL2RE8V    0x22800007
#define MASK_VL2RE8V     0xfff0707f
#define MATCH_VL2RE16V   0x22805007
#define MASK_VL2RE16V    0xfff0707f
#define MATCH_VL2RE32V   0x22806007
#define MASK_VL2RE32V    0xfff0707f
#define MATCH_VL2RE64V   0x22807007
#define MASK_VL2RE64V    0xfff0707f

#define MATCH_VL4RE8V    0x62800007
#define MASK_VL4RE8V     0xfff0707f
#define MATCH_VL4RE16V   0x62805007
#define MASK_VL4RE16V    0xfff0707f
#define MATCH_VL4RE32V   0x62806007
#define MASK_VL4RE32V    0xfff0707f
#define MATCH_VL4RE64V   0x62807007
#define MASK_VL4RE64V    0xfff0707f

#define MATCH_VL8RE8V    0xe2800007
#define MASK_VL8RE8V     0xfff0707f
#define MATCH_VL8RE16V   0xe2805007
#define MASK_VL8RE16V    0xfff0707f
#define MATCH_VL8RE32V   0xe2806007
#define MASK_VL8RE32V    0xfff0707f
#define MATCH_VL8RE64V   0xe2807007
#define MASK_VL8RE64V    0xfff0707f

#define MATCH_VS1RV  0x02800027
#define MASK_VS1RV   0xfff0707f
#define MATCH_VS2RV  0x22800027
#define MASK_VS2RV   0xfff0707f
#define MATCH_VS4RV  0x62800027
#define MASK_VS4RV   0xfff0707f
#define MATCH_VS8RV  0xe2800027
#define MASK_VS8RV   0xfff0707f

/* Temporary AMO encoding info

width
010 AMO*.W
011 AMO*.D
100 AMO*.Q
000 VAMO*EI8.V
101 VAMO*EI16.V
110 VAMO*EI32.V
111 VAMO*EI64.V

amoop
00001 vamoswap
00000 vamoadd
00100 vamoxor
01100 vamoand
01000 vamoor
10000 vamomin
10100 vamomax
11000 vamominu
11100 vamomaxu

   31-27 26 25 24-20 19-15 14-12 11-7    6-0
   amoop wd vm  vs2   rs1  width vs3/vd  opcode
   00001 x 1 xxxxx xxxxx 110 xxxxx 0101111
   0000 1x1x xxxx xxxx x110 xxxx x010 1111
   1111 1010 0000 0000 0111 0000 0111 1111 */

#define MATCH_VAMOADDEI8V   0x0000002f
#define MASK_VAMOADDEI8V    0xf800707f
#define MATCH_VAMOSWAPEI8V  0x0800002f
#define MASK_VAMOSWAPEI8V   0xf800707f
#define MATCH_VAMOXOREI8V   0x2000002f
#define MASK_VAMOXOREI8V    0xf800707f
#define MATCH_VAMOANDEI8V   0x6000002f
#define MASK_VAMOANDEI8V    0xf800707f
#define MATCH_VAMOOREI8V    0x4000002f
#define MASK_VAMOOREI8V     0xf800707f
#define MATCH_VAMOMINEI8V   0x8000002f
#define MASK_VAMOMINEI8V    0xf800707f
#define MATCH_VAMOMAXEI8V   0xa000002f
#define MASK_VAMOMAXEI8V    0xf800707f
#define MATCH_VAMOMINUEI8V  0xc000002f
#define MASK_VAMOMINUEI8V   0xf800707f
#define MATCH_VAMOMAXUEI8V  0xe000002f
#define MASK_VAMOMAXUEI8V   0xf800707f

#define MATCH_VAMOADDEI16V   0x0000502f
#define MASK_VAMOADDEI16V    0xf800707f
#define MATCH_VAMOSWAPEI16V  0x0800502f
#define MASK_VAMOSWAPEI16V   0xf800707f
#define MATCH_VAMOXOREI16V   0x2000502f
#define MASK_VAMOXOREI16V    0xf800707f
#define MATCH_VAMOANDEI16V   0x6000502f
#define MASK_VAMOANDEI16V    0xf800707f
#define MATCH_VAMOOREI16V    0x4000502f
#define MASK_VAMOOREI16V     0xf800707f
#define MATCH_VAMOMINEI16V   0x8000502f
#define MASK_VAMOMINEI16V    0xf800707f
#define MATCH_VAMOMAXEI16V   0xa000502f
#define MASK_VAMOMAXEI16V    0xf800707f
#define MATCH_VAMOMINUEI16V  0xc000502f
#define MASK_VAMOMINUEI16V   0xf800707f
#define MATCH_VAMOMAXUEI16V  0xe000502f
#define MASK_VAMOMAXUEI16V   0xf800707f

#define MATCH_VAMOADDEI32V   0x0000602f
#define MASK_VAMOADDEI32V    0xf800707f
#define MATCH_VAMOSWAPEI32V  0x0800602f
#define MASK_VAMOSWAPEI32V   0xf800707f
#define MATCH_VAMOXOREI32V   0x2000602f
#define MASK_VAMOXOREI32V    0xf800707f
#define MATCH_VAMOANDEI32V   0x6000602f
#define MASK_VAMOANDEI32V    0xf800707f
#define MATCH_VAMOOREI32V    0x4000602f
#define MASK_VAMOOREI32V     0xf800707f
#define MATCH_VAMOMINEI32V   0x8000602f
#define MASK_VAMOMINEI32V    0xf800707f
#define MATCH_VAMOMAXEI32V   0xa000602f
#define MASK_VAMOMAXEI32V    0xf800707f
#define MATCH_VAMOMINUEI32V  0xc000602f
#define MASK_VAMOMINUEI32V   0xf800707f
#define MATCH_VAMOMAXUEI32V  0xe000602f
#define MASK_VAMOMAXUEI32V   0xf800707f

#define MATCH_VAMOADDEI64V   0x0000702f
#define MASK_VAMOADDEI64V    0xf800707f
#define MATCH_VAMOSWAPEI64V  0x0800702f
#define MASK_VAMOSWAPEI64V   0xf800707f
#define MATCH_VAMOXOREI64V   0x2000702f
#define MASK_VAMOXOREI64V    0xf800707f
#define MATCH_VAMOANDEI64V   0x6000702f
#define MASK_VAMOANDEI64V    0xf800707f
#define MATCH_VAMOOREI64V    0x4000702f
#define MASK_VAMOOREI64V     0xf800707f
#define MATCH_VAMOMINEI64V   0x8000702f
#define MASK_VAMOMINEI64V    0xf800707f
#define MATCH_VAMOMAXEI64V   0xa000702f
#define MASK_VAMOMAXEI64V    0xf800707f
#define MATCH_VAMOMINUEI64V  0xc000702f
#define MASK_VAMOMINUEI64V   0xf800707f
#define MATCH_VAMOMAXUEI64V  0xe000702f
#define MASK_VAMOMAXUEI64V   0xf800707f

/* Temporary ALU encoding info

funct3
000 OPIVV vv
001 OPFVV vv
010 OPMVV vv
011 OPIVI vi  simm[4:0]
100 OPIVX vx  GPR x-reg rs1
101 OPFVF vf  FP f-reg rs1
110 OPMVX vx  GPR x-reg rs1
111 OPCFG si  GPR x-reg rs1 & rs2/imm

INT OPI
funct6
000000 vadd
000001
000010 vsub
000011 vrsub
000100 vminu
000101 vmin
000110 vmaxu
000111 vmax
001000
001001 vand
001010 vor
001011 vxor
001100 vrgather
001101
001110 vslideup, vrgatherei16
001111 vslidedown
010000 vadc
010001 vmadc
010010 vsbc
010011 vmsbc
010100
010101
010110
010111 vmerge/vmv
011000 vmseq
011001 vmsne
011010 vmsltu
011011 vmslt
011100 vmsleu
011101 vmsle
011110 vmsgtu
011111 vmsgt
100000 vsaddu
100001 vsadd
100010 vssubu
100011 vssub
100100
100101 vsll
100110
100111 vmv<nf>r (nf = 1, 2, 4, 8)
101000 vsrl
101001 vsra
101010 vssrl
101011 vssra
101100 vnsrl
101101 vnsra
101110 vnclipu
101111 vnclip
110000 vwredsumu
110001 vwredsum
110010
110011
110100
110101
110110
110111
111000 vdotu **
111001 vdot **
111010
111011
111100 vqmaccu
111101 vqmacc
111110 vqmaccus
111111 vqmaccsu

INT OPM
funct6
000000 vredsum
000001 vredand
000010 vredor
000011 vredxor
000100 vredminu
000101 vredmin
000110 vredmaxu
000111 vredmax
001000 vaaddu
001001 vaadd
001010 vasubu
001011 vasub
001100
001101
001110 vslide1up
001111 vslide1down
010000 VRXUNARY0/VWXUNARY0
010001
010010 VXUNARY0
010011
010100 VMUNARY0
010101
010110
010111 vcompress
011000 vmandnot
011001 vmand
011010 vmor
011011 vmxor
011100 vmornot
011101 vmnand
011110 vmnor
011111 vmxnor
100000 vdivu
100001 vdiv
100010 vremu
100011 vrem
100100 vmulhu
100101 vmul
100110 vmulhsu
100111 vmulh
101000
101001 vmadd
101010
101011 vnmsub
101100
101101 vmacc
101110
101111 vnmsac
110000 vwaddu
110001 vwadd
110010 vwsubu
110011 vwsub
110100 vwaddu.w
110101 vwadd.w
110110 vwsubu.w
110111 vwsub.w
111000 vwmulu
111001
111010 vwmulsu
111011 vwmul
111100 vwmaccu
111101 vwmacc
111110 vwmaccus
111111 vwmaccsu

VRXUNARY0
vs2, funct3=X
00000 vmv.s.x

VWXUNARY0
vs1, funct3=V
00000 vmv.x.s
10000 vpopc
10001 vfirst

VXUNARY0
vs1, funct3=V
00010 vzext.vf8
00011 vsext.vf8
00100 vzext.vf4
00101 vsext.vf4
00110 vzext.vf2
00111 vsext.vf2

VMUNARY0
rs1
00001 vmsbf
00010 vmsof
00011 vmsif
10000 viota
10001 vid

VFLOAT
funct6
000000 vfadd
000001 vfredsum
000010 vfsub
000011 vfredosum
000100 vfmin
000101 vfredmin
000110 vfmax
000111 vfredmax
001000 vfsgnj
001001 vfsgnn
001010 vfsgnx
001011
001100
001101
001110 vfslide1up
001111 vfslide1down
010000 VRFUNARY0/VWFUNARY0
010001
010010 VFUNARY0
010011 VFUNARY1
010100
010101
010110
010111 vfmerge/vfmv
011000 vmfeq
011001 vmfle
011010
011011 vmflt
011100 vmfne
011101 vmfgt
011110
011111 vmfge
100000 vfdiv
100001 vfrdiv
100010
100011
100100 vfmul
100101
100110
100111 vfrsub
101000 vfmadd
101001 vfnmadd
101010 vfmsub
101011 vfnmsub
101100 vfmacc
101101 vfnmacc
101110 vfmsac
101111 vfnmsac
110000 vfwadd
110001 vfwredsum
110010 vfwsub
110011 vfwredosum
110100 vfwadd.w
110101
110110 vfwsub.w
110111
111000 vfwmul
111001 vfdot
111010
111011
111100 vfwmacc
111101 vfwnmacc
111110 vfwmsac
111111 vfwnmsac

VRFUNARY0
vs2, funct3=F
00000 vfmv.s.f

VWFUNARY0
vs1, funct3=V
00000 vfmv.f.s

VFUNARY0
vs1
00000 vfcvt.xu.f.v
00001 vfcvt.x.f.v
00010 vfcvt.f.xu.v
00011 vfcvt.f.x.v
00110 vfcvt.rtz.xu.f.v
00111 vfcvt.rtz.x.f.v

01000 vfwcvt.xu.f.v
01001 vfwcvt.x.f.v
01010 vfwcvt.f.xu.v
01011 vfwcvt.f.x.v
01100 vfwcvt.f.f.v
01110 vfwcvt.rtz.xu.f.v
01111 vfwcvt.rtz.x.f.v

10000 vfncvt.xu.f.w
10001 vfncvt.x.f.w
10010 vfncvt.f.xu.w
10011 vfncvt.f.x.w
10100 vfncvt.f.f.w
10101 vfncvt.rod.f.f.w
10110 vfncvt.rtz.xu.f.v
10111 vfncvt.rtz.x.f.v

VFUNARY1
vs1
00000 vfsqrt.v
00100 vfrsqrte7.v
00101 vfrece7.v
10000 vfclass.v

31-26 25 24-20   19-15     14-12 11-7 6-0
funct6 VM  VS2  VS1/RS1/IMM funct3 VD   opcode
010000 x xxxxx 00000 001 xxxxx 1010111
0100 00xx xxxx 0000 0001 xxxx x101 0111
*/

#define MATCH_VADDVV  0x00000057
#define MASK_VADDVV   0xfc00707f
#define MATCH_VADDVX  0x00004057
#define MASK_VADDVX   0xfc00707f
#define MATCH_VADDVI  0x00003057
#define MASK_VADDVI   0xfc00707f
#define MATCH_VSUBVV  0x08000057
#define MASK_VSUBVV   0xfc00707f
#define MATCH_VSUBVX  0x08004057
#define MASK_VSUBVX   0xfc00707f
#define MATCH_VRSUBVX 0x0c004057
#define MASK_VRSUBVX  0xfc00707f
#define MATCH_VRSUBVI 0x0c003057
#define MASK_VRSUBVI  0xfc00707f

#define MATCH_VWCVTXXV  0xc4006057
#define MASK_VWCVTXXV   0xfc0ff07f
#define MATCH_VWCVTUXXV 0xc0006057
#define MASK_VWCVTUXXV  0xfc0ff07f

#define MATCH_VWADDVV  0xc4002057
#define MASK_VWADDVV   0xfc00707f
#define MATCH_VWADDVX  0xc4006057
#define MASK_VWADDVX   0xfc00707f
#define MATCH_VWSUBVV  0xcc002057
#define MASK_VWSUBVV   0xfc00707f
#define MATCH_VWSUBVX  0xcc006057
#define MASK_VWSUBVX   0xfc00707f
#define MATCH_VWADDWV  0xd4002057
#define MASK_VWADDWV   0xfc00707f
#define MATCH_VWADDWX  0xd4006057
#define MASK_VWADDWX   0xfc00707f
#define MATCH_VWSUBWV  0xdc002057
#define MASK_VWSUBWV   0xfc00707f
#define MATCH_VWSUBWX  0xdc006057
#define MASK_VWSUBWX   0xfc00707f
#define MATCH_VWADDUVV  0xc0002057
#define MASK_VWADDUVV   0xfc00707f
#define MATCH_VWADDUVX  0xc0006057
#define MASK_VWADDUVX   0xfc00707f
#define MATCH_VWSUBUVV  0xc8002057
#define MASK_VWSUBUVV   0xfc00707f
#define MATCH_VWSUBUVX  0xc8006057
#define MASK_VWSUBUVX   0xfc00707f
#define MATCH_VWADDUWV  0xd0002057
#define MASK_VWADDUWV   0xfc00707f
#define MATCH_VWADDUWX  0xd0006057
#define MASK_VWADDUWX   0xfc00707f
#define MATCH_VWSUBUWV  0xd8002057
#define MASK_VWSUBUWV   0xfc00707f
#define MATCH_VWSUBUWX  0xd8006057
#define MASK_VWSUBUWX   0xfc00707f

#define MATCH_VZEXT_VF8 0x48012057
#define MASK_VZEXT_VF8  0xfc0ff07f
#define MATCH_VSEXT_VF8 0x4801a057
#define MASK_VSEXT_VF8  0xfc0ff07f
#define MATCH_VZEXT_VF4 0x48022057
#define MASK_VZEXT_VF4  0xfc0ff07f
#define MATCH_VSEXT_VF4 0x4802a057
#define MASK_VSEXT_VF4  0xfc0ff07f
#define MATCH_VZEXT_VF2 0x48032057
#define MASK_VZEXT_VF2  0xfc0ff07f
#define MATCH_VSEXT_VF2 0x4803a057
#define MASK_VSEXT_VF2  0xfc0ff07f

#define MATCH_VADCVVM  0x40000057
#define MASK_VADCVVM   0xfe00707f
#define MATCH_VADCVXM  0x40004057
#define MASK_VADCVXM   0xfe00707f
#define MATCH_VADCVIM  0x40003057
#define MASK_VADCVIM   0xfe00707f
#define MATCH_VMADCVVM 0x44000057
#define MASK_VMADCVVM  0xfe00707f
#define MATCH_VMADCVXM 0x44004057
#define MASK_VMADCVXM  0xfe00707f
#define MATCH_VMADCVIM 0x44003057
#define MASK_VMADCVIM  0xfe00707f
#define MATCH_VMADCVV  0x46000057
#define MASK_VMADCVV   0xfe00707f
#define MATCH_VMADCVX  0x46004057
#define MASK_VMADCVX   0xfe00707f
#define MATCH_VMADCVI  0x46003057
#define MASK_VMADCVI   0xfe00707f
#define MATCH_VSBCVVM  0x48000057
#define MASK_VSBCVVM   0xfe00707f
#define MATCH_VSBCVXM  0x48004057
#define MASK_VSBCVXM   0xfe00707f
#define MATCH_VMSBCVVM 0x4c000057
#define MASK_VMSBCVVM  0xfe00707f
#define MATCH_VMSBCVXM 0x4c004057
#define MASK_VMSBCVXM  0xfe00707f
#define MATCH_VMSBCVV  0x4e000057
#define MASK_VMSBCVV   0xfe00707f
#define MATCH_VMSBCVX  0x4e004057
#define MASK_VMSBCVX   0xfe00707f

#define MATCH_VNOTV   0x2c0fb057
#define MASK_VNOTV    0xfc0ff07f

#define MATCH_VANDVV  0x24000057
#define MASK_VANDVV   0xfc00707f
#define MATCH_VANDVX  0x24004057
#define MASK_VANDVX   0xfc00707f
#define MATCH_VANDVI  0x24003057
#define MASK_VANDVI   0xfc00707f
#define MATCH_VORVV   0x28000057
#define MASK_VORVV    0xfc00707f
#define MATCH_VORVX   0x28004057
#define MASK_VORVX    0xfc00707f
#define MATCH_VORVI   0x28003057
#define MASK_VORVI    0xfc00707f
#define MATCH_VXORVV  0x2c000057
#define MASK_VXORVV   0xfc00707f
#define MATCH_VXORVX  0x2c004057
#define MASK_VXORVX   0xfc00707f
#define MATCH_VXORVI  0x2c003057
#define MASK_VXORVI   0xfc00707f

#define MATCH_VSLLVV 0x94000057
#define MASK_VSLLVV  0xfc00707f
#define MATCH_VSLLVX 0x94004057
#define MASK_VSLLVX  0xfc00707f
#define MATCH_VSLLVI 0x94003057
#define MASK_VSLLVI  0xfc00707f
#define MATCH_VSRLVV 0xa0000057
#define MASK_VSRLVV  0xfc00707f
#define MATCH_VSRLVX 0xa0004057
#define MASK_VSRLVX  0xfc00707f
#define MATCH_VSRLVI 0xa0003057
#define MASK_VSRLVI  0xfc00707f
#define MATCH_VSRAVV 0xa4000057
#define MASK_VSRAVV  0xfc00707f
#define MATCH_VSRAVX 0xa4004057
#define MASK_VSRAVX  0xfc00707f
#define MATCH_VSRAVI 0xa4003057
#define MASK_VSRAVI  0xfc00707f

#define MATCH_VNCVTXXW 0xb0004057
#define MASK_VNCVTXXW  0xfc0ff07f

#define MATCH_VNSRLWV  0xb0000057
#define MASK_VNSRLWV   0xfc00707f
#define MATCH_VNSRLWX  0xb0004057
#define MASK_VNSRLWX   0xfc00707f
#define MATCH_VNSRLWI  0xb0003057
#define MASK_VNSRLWI   0xfc00707f
#define MATCH_VNSRAWV  0xb4000057
#define MASK_VNSRAWV   0xfc00707f
#define MATCH_VNSRAWX  0xb4004057
#define MASK_VNSRAWX   0xfc00707f
#define MATCH_VNSRAWI  0xb4003057
#define MASK_VNSRAWI   0xfc00707f

#define MATCH_VMSEQVV  0x60000057
#define MASK_VMSEQVV   0xfc00707f
#define MATCH_VMSEQVX  0x60004057
#define MASK_VMSEQVX   0xfc00707f
#define MATCH_VMSEQVI  0x60003057
#define MASK_VMSEQVI   0xfc00707f
#define MATCH_VMSNEVV  0x64000057
#define MASK_VMSNEVV   0xfc00707f
#define MATCH_VMSNEVX  0x64004057
#define MASK_VMSNEVX   0xfc00707f
#define MATCH_VMSNEVI  0x64003057
#define MASK_VMSNEVI   0xfc00707f
#define MATCH_VMSLTVV  0x6c000057
#define MASK_VMSLTVV   0xfc00707f
#define MATCH_VMSLTVX  0x6c004057
#define MASK_VMSLTVX   0xfc00707f
#define MATCH_VMSLTUVV 0x68000057
#define MASK_VMSLTUVV  0xfc00707f
#define MATCH_VMSLTUVX 0x68004057
#define MASK_VMSLTUVX  0xfc00707f
#define MATCH_VMSLEVV  0x74000057
#define MASK_VMSLEVV   0xfc00707f
#define MATCH_VMSLEVX  0x74004057
#define MASK_VMSLEVX   0xfc00707f
#define MATCH_VMSLEVI  0x74003057
#define MASK_VMSLEVI   0xfc00707f
#define MATCH_VMSLEUVV 0x70000057
#define MASK_VMSLEUVV  0xfc00707f
#define MATCH_VMSLEUVX 0x70004057
#define MASK_VMSLEUVX  0xfc00707f
#define MATCH_VMSLEUVI 0x70003057
#define MASK_VMSLEUVI  0xfc00707f
#define MATCH_VMSGTVX  0x7c004057
#define MASK_VMSGTVX   0xfc00707f
#define MATCH_VMSGTVI  0x7c003057
#define MASK_VMSGTVI   0xfc00707f
#define MATCH_VMSGTUVX 0x78004057
#define MASK_VMSGTUVX  0xfc00707f
#define MATCH_VMSGTUVI 0x78003057
#define MASK_VMSGTUVI  0xfc00707f

#define MATCH_VMINVV  0x14000057
#define MASK_VMINVV   0xfc00707f
#define MATCH_VMINVX  0x14004057
#define MASK_VMINVX   0xfc00707f
#define MATCH_VMAXVV  0x1c000057
#define MASK_VMAXVV   0xfc00707f
#define MATCH_VMAXVX  0x1c004057
#define MASK_VMAXVX   0xfc00707f
#define MATCH_VMINUVV 0x10000057
#define MASK_VMINUVV  0xfc00707f
#define MATCH_VMINUVX 0x10004057
#define MASK_VMINUVX  0xfc00707f
#define MATCH_VMAXUVV 0x18000057
#define MASK_VMAXUVV  0xfc00707f
#define MATCH_VMAXUVX 0x18004057
#define MASK_VMAXUVX  0xfc00707f

#define MATCH_VMULVV    0x94002057
#define MASK_VMULVV     0xfc00707f
#define MATCH_VMULVX    0x94006057
#define MASK_VMULVX     0xfc00707f
#define MATCH_VMULHVV   0x9c002057
#define MASK_VMULHVV    0xfc00707f
#define MATCH_VMULHVX   0x9c006057
#define MASK_VMULHVX    0xfc00707f
#define MATCH_VMULHUVV  0x90002057
#define MASK_VMULHUVV   0xfc00707f
#define MATCH_VMULHUVX  0x90006057
#define MASK_VMULHUVX   0xfc00707f
#define MATCH_VMULHSUVV 0x98002057
#define MASK_VMULHSUVV  0xfc00707f
#define MATCH_VMULHSUVX 0x98006057
#define MASK_VMULHSUVX  0xfc00707f

#define MATCH_VWMULVV   0xec002057
#define MASK_VWMULVV    0xfc00707f
#define MATCH_VWMULVX   0xec006057
#define MASK_VWMULVX    0xfc00707f
#define MATCH_VWMULUVV  0xe0002057
#define MASK_VWMULUVV   0xfc00707f
#define MATCH_VWMULUVX  0xe0006057
#define MASK_VWMULUVX   0xfc00707f
#define MATCH_VWMULSUVV 0xe8002057
#define MASK_VWMULSUVV  0xfc00707f
#define MATCH_VWMULSUVX 0xe8006057
#define MASK_VWMULSUVX  0xfc00707f

#define MATCH_VMACCVV  0xb4002057
#define MASK_VMACCVV   0xfc00707f
#define MATCH_VMACCVX  0xb4006057
#define MASK_VMACCVX   0xfc00707f
#define MATCH_VNMSACVV 0xbc002057
#define MASK_VNMSACVV  0xfc00707f
#define MATCH_VNMSACVX 0xbc006057
#define MASK_VNMSACVX  0xfc00707f
#define MATCH_VMADDVV  0xa4002057
#define MASK_VMADDVV   0xfc00707f
#define MATCH_VMADDVX  0xa4006057
#define MASK_VMADDVX   0xfc00707f
#define MATCH_VNMSUBVV 0xac002057
#define MASK_VNMSUBVV  0xfc00707f
#define MATCH_VNMSUBVX 0xac006057
#define MASK_VNMSUBVX  0xfc00707f

#define MATCH_VWMACCUVV  0xf0002057
#define MASK_VWMACCUVV   0xfc00707f
#define MATCH_VWMACCUVX  0xf0006057
#define MASK_VWMACCUVX   0xfc00707f
#define MATCH_VWMACCVV   0xf4002057
#define MASK_VWMACCVV    0xfc00707f
#define MATCH_VWMACCVX   0xf4006057
#define MASK_VWMACCVX    0xfc00707f
#define MATCH_VWMACCSUVV 0xfc002057
#define MASK_VWMACCSUVV  0xfc00707f
#define MATCH_VWMACCSUVX 0xfc006057
#define MASK_VWMACCSUVX  0xfc00707f
#define MATCH_VWMACCUSVX 0xf8006057
#define MASK_VWMACCUSVX  0xfc00707f

#define MATCH_VQMACCUVV  0xf0000057
#define MASK_VQMACCUVV   0xfc00707f
#define MATCH_VQMACCUVX  0xf0004057
#define MASK_VQMACCUVX   0xfc00707f
#define MATCH_VQMACCVV   0xf4000057
#define MASK_VQMACCVV    0xfc00707f
#define MATCH_VQMACCVX   0xf4004057
#define MASK_VQMACCVX    0xfc00707f
#define MATCH_VQMACCSUVV 0xfc000057
#define MASK_VQMACCSUVV  0xfc00707f
#define MATCH_VQMACCSUVX 0xfc004057
#define MASK_VQMACCSUVX  0xfc00707f
#define MATCH_VQMACCUSVX 0xf8004057
#define MASK_VQMACCUSVX  0xfc00707f

#define MATCH_VDIVVV  0x84002057
#define MASK_VDIVVV   0xfc00707f
#define MATCH_VDIVVX  0x84006057
#define MASK_VDIVVX   0xfc00707f
#define MATCH_VDIVUVV 0x80002057
#define MASK_VDIVUVV  0xfc00707f
#define MATCH_VDIVUVX 0x80006057
#define MASK_VDIVUVX  0xfc00707f
#define MATCH_VREMVV  0x8c002057
#define MASK_VREMVV   0xfc00707f
#define MATCH_VREMVX  0x8c006057
#define MASK_VREMVX   0xfc00707f
#define MATCH_VREMUVV 0x88002057
#define MASK_VREMUVV  0xfc00707f
#define MATCH_VREMUVX 0x88006057
#define MASK_VREMUVX  0xfc00707f

#define MATCH_VMERGEVVM 0x5c000057
#define MASK_VMERGEVVM  0xfe00707f
#define MATCH_VMERGEVXM 0x5c004057
#define MASK_VMERGEVXM  0xfe00707f
#define MATCH_VMERGEVIM 0x5c003057
#define MASK_VMERGEVIM  0xfe00707f

#define MATCH_VMVVV    0x5e000057
#define MASK_VMVVV     0xfff0707f
#define MATCH_VMVVX    0x5e004057
#define MASK_VMVVX     0xfff0707f
#define MATCH_VMVVI    0x5e003057
#define MASK_VMVVI     0xfff0707f

#define MATCH_VSADDUVV 0x80000057
#define MASK_VSADDUVV  0xfc00707f
#define MATCH_VSADDUVX 0x80004057
#define MASK_VSADDUVX  0xfc00707f
#define MATCH_VSADDUVI 0x80003057
#define MASK_VSADDUVI  0xfc00707f
#define MATCH_VSADDVV  0x84000057
#define MASK_VSADDVV   0xfc00707f
#define MATCH_VSADDVX  0x84004057
#define MASK_VSADDVX   0xfc00707f
#define MATCH_VSADDVI  0x84003057
#define MASK_VSADDVI   0xfc00707f
#define MATCH_VSSUBUVV 0x88000057
#define MASK_VSSUBUVV  0xfc00707f
#define MATCH_VSSUBUVX 0x88004057
#define MASK_VSSUBUVX  0xfc00707f
#define MATCH_VSSUBVV  0x8c000057
#define MASK_VSSUBVV   0xfc00707f
#define MATCH_VSSUBVX  0x8c004057
#define MASK_VSSUBVX   0xfc00707f

#define MATCH_VAADDUVV 0x20002057
#define MASK_VAADDUVV  0xfc00707f
#define MATCH_VAADDUVX 0x20006057
#define MASK_VAADDUVX  0xfc00707f
#define MATCH_VAADDVV  0x24002057
#define MASK_VAADDVV   0xfc00707f
#define MATCH_VAADDVX  0x24006057
#define MASK_VAADDVX   0xfc00707f
#define MATCH_VASUBUVV 0x28002057
#define MASK_VASUBUVV  0xfc00707f
#define MATCH_VASUBUVX 0x28006057
#define MASK_VASUBUVX  0xfc00707f
#define MATCH_VASUBVV  0x2c002057
#define MASK_VASUBVV   0xfc00707f
#define MATCH_VASUBVX  0x2c006057
#define MASK_VASUBVX   0xfc00707f

#define MATCH_VSMULVV  0x9c000057
#define MASK_VSMULVV   0xfc00707f
#define MATCH_VSMULVX  0x9c004057
#define MASK_VSMULVX   0xfc00707f

#define MATCH_VSSRLVV   0xa8000057
#define MASK_VSSRLVV    0xfc00707f
#define MATCH_VSSRLVX   0xa8004057
#define MASK_VSSRLVX    0xfc00707f
#define MATCH_VSSRLVI   0xa8003057
#define MASK_VSSRLVI    0xfc00707f
#define MATCH_VSSRAVV   0xac000057
#define MASK_VSSRAVV    0xfc00707f
#define MATCH_VSSRAVX   0xac004057
#define MASK_VSSRAVX    0xfc00707f
#define MATCH_VSSRAVI   0xac003057
#define MASK_VSSRAVI    0xfc00707f

#define MATCH_VNCLIPUWV 0xb8000057
#define MASK_VNCLIPUWV  0xfc00707f
#define MATCH_VNCLIPUWX 0xb8004057
#define MASK_VNCLIPUWX  0xfc00707f
#define MATCH_VNCLIPUWI 0xb8003057
#define MASK_VNCLIPUWI  0xfc00707f
#define MATCH_VNCLIPWV  0xbc000057
#define MASK_VNCLIPWV   0xfc00707f
#define MATCH_VNCLIPWX  0xbc004057
#define MASK_VNCLIPWX   0xfc00707f
#define MATCH_VNCLIPWI  0xbc003057
#define MASK_VNCLIPWI   0xfc00707f

#define MATCH_VFADDVV  0x00001057
#define MASK_VFADDVV   0xfc00707f
#define MATCH_VFADDVF  0x00005057
#define MASK_VFADDVF   0xfc00707f
#define MATCH_VFSUBVV  0x08001057
#define MASK_VFSUBVV   0xfc00707f
#define MATCH_VFSUBVF  0x08005057
#define MASK_VFSUBVF   0xfc00707f
#define MATCH_VFRSUBVF 0x9c005057
#define MASK_VFRSUBVF  0xfc00707f

#define MATCH_VFWADDVV  0xc0001057
#define MASK_VFWADDVV   0xfc00707f
#define MATCH_VFWADDVF  0xc0005057
#define MASK_VFWADDVF   0xfc00707f
#define MATCH_VFWSUBVV  0xc8001057
#define MASK_VFWSUBVV   0xfc00707f
#define MATCH_VFWSUBVF  0xc8005057
#define MASK_VFWSUBVF   0xfc00707f
#define MATCH_VFWADDWV  0xd0001057
#define MASK_VFWADDWV   0xfc00707f
#define MATCH_VFWADDWF  0xd0005057
#define MASK_VFWADDWF   0xfc00707f
#define MATCH_VFWSUBWV  0xd8001057
#define MASK_VFWSUBWV   0xfc00707f
#define MATCH_VFWSUBWF  0xd8005057
#define MASK_VFWSUBWF   0xfc00707f

#define MATCH_VFMULVV  0x90001057
#define MASK_VFMULVV   0xfc00707f
#define MATCH_VFMULVF  0x90005057
#define MASK_VFMULVF   0xfc00707f
#define MATCH_VFDIVVV  0x80001057
#define MASK_VFDIVVV   0xfc00707f
#define MATCH_VFDIVVF  0x80005057
#define MASK_VFDIVVF   0xfc00707f
#define MATCH_VFRDIVVF 0x84005057
#define MASK_VFRDIVVF  0xfc00707f

#define MATCH_VFWMULVV 0xe0001057
#define MASK_VFWMULVV  0xfc00707f
#define MATCH_VFWMULVF 0xe0005057
#define MASK_VFWMULVF  0xfc00707f

#define MATCH_VFMADDVV  0xa0001057
#define MASK_VFMADDVV   0xfc00707f
#define MATCH_VFMADDVF  0xa0005057
#define MASK_VFMADDVF   0xfc00707f
#define MATCH_VFNMADDVV 0xa4001057
#define MASK_VFNMADDVV  0xfc00707f
#define MATCH_VFNMADDVF 0xa4005057
#define MASK_VFNMADDVF  0xfc00707f
#define MATCH_VFMSUBVV  0xa8001057
#define MASK_VFMSUBVV   0xfc00707f
#define MATCH_VFMSUBVF  0xa8005057
#define MASK_VFMSUBVF   0xfc00707f
#define MATCH_VFNMSUBVV 0xac001057
#define MASK_VFNMSUBVV  0xfc00707f
#define MATCH_VFNMSUBVF 0xac005057
#define MASK_VFNMSUBVF  0xfc00707f
#define MATCH_VFMACCVV  0xb0001057
#define MASK_VFMACCVV   0xfc00707f
#define MATCH_VFMACCVF  0xb0005057
#define MASK_VFMACCVF   0xfc00707f
#define MATCH_VFNMACCVV 0xb4001057
#define MASK_VFNMACCVV  0xfc00707f
#define MATCH_VFNMACCVF 0xb4005057
#define MASK_VFNMACCVF  0xfc00707f
#define MATCH_VFMSACVV  0xb8001057
#define MASK_VFMSACVV   0xfc00707f
#define MATCH_VFMSACVF  0xb8005057
#define MASK_VFMSACVF   0xfc00707f
#define MATCH_VFNMSACVV 0xbc001057
#define MASK_VFNMSACVV  0xfc00707f
#define MATCH_VFNMSACVF 0xbc005057
#define MASK_VFNMSACVF  0xfc00707f

#define MATCH_VFWMACCVV  0xf0001057
#define MASK_VFWMACCVV   0xfc00707f
#define MATCH_VFWMACCVF  0xf0005057
#define MASK_VFWMACCVF   0xfc00707f
#define MATCH_VFWNMACCVV 0xf4001057
#define MASK_VFWNMACCVV  0xfc00707f
#define MATCH_VFWNMACCVF 0xf4005057
#define MASK_VFWNMACCVF  0xfc00707f
#define MATCH_VFWMSACVV  0xf8001057
#define MASK_VFWMSACVV   0xfc00707f
#define MATCH_VFWMSACVF  0xf8005057
#define MASK_VFWMSACVF   0xfc00707f
#define MATCH_VFWNMSACVV 0xfc001057
#define MASK_VFWNMSACVV  0xfc00707f
#define MATCH_VFWNMSACVF 0xfc005057
#define MASK_VFWNMSACVF  0xfc00707f

#define MATCH_VFSQRTV    0x4c001057
#define MASK_VFSQRTV     0xfc0ff07f
#define MATCH_VFRSQRT7V  0x4c021057
#define MASK_VFRSQRT7V   0xfc0ff07f
#define MATCH_VFREC7V    0x4c029057
#define MASK_VFREC7V     0xfc0ff07f
#define MATCH_VFCLASSV   0x4c081057
#define MASK_VFCLASSV    0xfc0ff07f

#define MATCH_VFMINVV  0x10001057
#define MASK_VFMINVV   0xfc00707f
#define MATCH_VFMINVF  0x10005057
#define MASK_VFMINVF   0xfc00707f
#define MATCH_VFMAXVV  0x18001057
#define MASK_VFMAXVV   0xfc00707f
#define MATCH_VFMAXVF  0x18005057
#define MASK_VFMAXVF   0xfc00707f

#define MATCH_VFSGNJVV  0x20001057
#define MASK_VFSGNJVV   0xfc00707f
#define MATCH_VFSGNJVF  0x20005057
#define MASK_VFSGNJVF   0xfc00707f
#define MATCH_VFSGNJNVV 0x24001057
#define MASK_VFSGNJNVV  0xfc00707f
#define MATCH_VFSGNJNVF 0x24005057
#define MASK_VFSGNJNVF  0xfc00707f
#define MATCH_VFSGNJXVV 0x28001057
#define MASK_VFSGNJXVV  0xfc00707f
#define MATCH_VFSGNJXVF 0x28005057
#define MASK_VFSGNJXVF  0xfc00707f

#define MATCH_VMFEQVV   0x60001057
#define MASK_VMFEQVV    0xfc00707f
#define MATCH_VMFEQVF   0x60005057
#define MASK_VMFEQVF    0xfc00707f
#define MATCH_VMFNEVV   0x70001057
#define MASK_VMFNEVV    0xfc00707f
#define MATCH_VMFNEVF   0x70005057
#define MASK_VMFNEVF    0xfc00707f
#define MATCH_VMFLTVV   0x6c001057
#define MASK_VMFLTVV    0xfc00707f
#define MATCH_VMFLTVF   0x6c005057
#define MASK_VMFLTVF    0xfc00707f
#define MATCH_VMFLEVV  0x64001057
#define MASK_VMFLEVV   0xfc00707f
#define MATCH_VMFLEVF  0x64005057
#define MASK_VMFLEVF   0xfc00707f
#define MATCH_VMFGTVF   0x74005057
#define MASK_VMFGTVF    0xfc00707f
#define MATCH_VMFGEVF  0x7c005057
#define MASK_VMFGEVF   0xfc00707f

#define MATCH_VFMERGEVFM 0x5c005057
#define MASK_VFMERGEVFM  0xfe00707f
#define MATCH_VFMVVF     0x5e005057
#define MASK_VFMVVF      0xfff0707f

#define MATCH_VFCVTXUFV 0x48001057
#define MASK_VFCVTXUFV  0xfc0ff07f
#define MATCH_VFCVTXFV 0x48009057
#define MASK_VFCVTXFV  0xfc0ff07f
#define MATCH_VFCVTFXUV 0x48011057
#define MASK_VFCVTFXUV  0xfc0ff07f
#define MATCH_VFCVTFXV 0x48019057
#define MASK_VFCVTFXV  0xfc0ff07f
#define MATCH_VFCVTRTZXUFV 0x48031057
#define MASK_VFCVTRTZXUFV  0xfc0ff07f
#define MATCH_VFCVTRTZXFV 0x48039057
#define MASK_VFCVTRTZXFV  0xfc0ff07f
#define MATCH_VFWCVTXUFV 0x48041057
#define MASK_VFWCVTXUFV  0xfc0ff07f
#define MATCH_VFWCVTXFV 0x48049057
#define MASK_VFWCVTXFV  0xfc0ff07f
#define MATCH_VFWCVTFXUV 0x48051057
#define MASK_VFWCVTFXUV  0xfc0ff07f
#define MATCH_VFWCVTFXV 0x48059057
#define MASK_VFWCVTFXV  0xfc0ff07f
#define MATCH_VFWCVTFFV 0x48061057
#define MASK_VFWCVTFFV  0xfc0ff07f
#define MATCH_VFWCVTRTZXUFV 0x48071057
#define MASK_VFWCVTRTZXUFV  0xfc0ff07f
#define MATCH_VFWCVTRTZXFV 0x48079057
#define MASK_VFWCVTRTZXFV  0xfc0ff07f
#define MATCH_VFNCVTXUFW 0x48081057
#define MASK_VFNCVTXUFW  0xfc0ff07f
#define MATCH_VFNCVTXFW 0x48089057
#define MASK_VFNCVTXFW  0xfc0ff07f
#define MATCH_VFNCVTFXUW 0x48091057
#define MASK_VFNCVTFXUW  0xfc0ff07f
#define MATCH_VFNCVTFXW 0x48099057
#define MASK_VFNCVTFXW  0xfc0ff07f
#define MATCH_VFNCVTFFW 0x480a1057
#define MASK_VFNCVTFFW  0xfc0ff07f
#define MATCH_VFNCVTRODFFW 0x480a9057
#define MASK_VFNCVTRODFFW  0xfc0ff07f
#define MATCH_VFNCVTRTZXUFW 0x480b1057
#define MASK_VFNCVTRTZXUFW  0xfc0ff07f
#define MATCH_VFNCVTRTZXFW 0x480b9057
#define MASK_VFNCVTRTZXFW  0xfc0ff07f

#define MATCH_VREDSUMVS  0x00002057
#define MASK_VREDSUMVS   0xfc00707f
#define MATCH_VREDMAXVS  0x1c002057
#define MASK_VREDMAXVS   0xfc00707f
#define MATCH_VREDMAXUVS 0x18002057
#define MASK_VREDMAXUVS  0xfc00707f
#define MATCH_VREDMINVS  0x14002057
#define MASK_VREDMINVS   0xfc00707f
#define MATCH_VREDMINUVS 0x10002057
#define MASK_VREDMINUVS  0xfc00707f
#define MATCH_VREDANDVS  0x04002057
#define MASK_VREDANDVS   0xfc00707f
#define MATCH_VREDORVS   0x08002057
#define MASK_VREDORVS    0xfc00707f
#define MATCH_VREDXORVS  0x0c002057
#define MASK_VREDXORVS   0xfc00707f

#define MATCH_VWREDSUMUVS 0xc0000057
#define MASK_VWREDSUMUVS  0xfc00707f
#define MATCH_VWREDSUMVS  0xc4000057
#define MASK_VWREDSUMVS   0xfc00707f

#define MATCH_VFREDOSUMVS 0x0c001057
#define MASK_VFREDOSUMVS  0xfc00707f
#define MATCH_VFREDSUMVS  0x04001057
#define MASK_VFREDSUMVS   0xfc00707f
#define MATCH_VFREDMAXVS  0x1c001057
#define MASK_VFREDMAXVS   0xfc00707f
#define MATCH_VFREDMINVS  0x14001057
#define MASK_VFREDMINVS   0xfc00707f

#define MATCH_VFWREDOSUMVS 0xcc001057
#define MASK_VFWREDOSUMVS  0xfc00707f
#define MATCH_VFWREDSUMVS  0xc4001057
#define MASK_VFWREDSUMVS   0xfc00707f

#define MATCH_VMANDMM    0x66002057
#define MASK_VMANDMM     0xfe00707f
#define MATCH_VMNANDMM   0x76002057
#define MASK_VMNANDMM    0xfe00707f
#define MATCH_VMANDNOTMM 0x62002057
#define MASK_VMANDNOTMM  0xfe00707f
#define MATCH_VMXORMM    0x6e002057
#define MASK_VMXORMM     0xfe00707f
#define MATCH_VMORMM     0x6a002057
#define MASK_VMORMM      0xfe00707f
#define MATCH_VMNORMM    0x7a002057
#define MASK_VMNORMM     0xfe00707f
#define MATCH_VMORNOTMM  0x72002057
#define MASK_VMORNOTMM   0xfe00707f
#define MATCH_VMXNORMM   0x7e002057
#define MASK_VMXNORMM    0xfe00707f

#define MATCH_VPOPCM   0x40082057
#define MASK_VPOPCM    0xfc0ff07f
#define MATCH_VFIRSTM  0x4008a057
#define MASK_VFIRSTM   0xfc0ff07f

#define MATCH_VMSBFM   0x5000a057
#define MASK_VMSBFM    0xfc0ff07f
#define MATCH_VMSIFM   0x5001a057
#define MASK_VMSIFM    0xfc0ff07f
#define MATCH_VMSOFM   0x50012057
#define MASK_VMSOFM    0xfc0ff07f
#define MATCH_VIOTAM   0x50082057
#define MASK_VIOTAM    0xfc0ff07f
#define MATCH_VIDV     0x5008a057
#define MASK_VIDV      0xfdfff07f

#define MATCH_VMVXS    0x42002057
#define MASK_VMVXS     0xfe0ff07f
#define MATCH_VMVSX    0x42006057
#define MASK_VMVSX     0xfff0707f

#define MATCH_VFMVFS   0x42001057
#define MASK_VFMVFS    0xfe0ff07f
#define MATCH_VFMVSF   0x42005057
#define MASK_VFMVSF    0xfff0707f

#define MATCH_VSLIDEUPVX   0x38004057
#define MASK_VSLIDEUPVX    0xfc00707f
#define MATCH_VSLIDEUPVI   0x38003057
#define MASK_VSLIDEUPVI    0xfc00707f
#define MATCH_VSLIDEDOWNVX 0x3c004057
#define MASK_VSLIDEDOWNVX  0xfc00707f
#define MATCH_VSLIDEDOWNVI 0x3c003057
#define MASK_VSLIDEDOWNVI  0xfc00707f

#define MATCH_VSLIDE1UPVX   0x38006057
#define MASK_VSLIDE1UPVX    0xfc00707f
#define MATCH_VSLIDE1DOWNVX 0x3c006057
#define MASK_VSLIDE1DOWNVX  0xfc00707f

#define MATCH_VFSLIDE1UPVF   0x38005057
#define MASK_VFSLIDE1UPVF    0xfc00707f
#define MATCH_VFSLIDE1DOWNVF 0x3c005057
#define MASK_VFSLIDE1DOWNVF  0xfc00707f

#define MATCH_VRGATHERVV      0x30000057
#define MASK_VRGATHERVV       0xfc00707f
#define MATCH_VRGATHERVX      0x30004057
#define MASK_VRGATHERVX       0xfc00707f
#define MATCH_VRGATHERVI      0x30003057
#define MASK_VRGATHERVI       0xfc00707f
#define MATCH_VRGATHEREI16VV  0x38000057
#define MASK_VRGATHEREI16VV   0xfc00707f

#define MATCH_VCOMPRESSVM   0x5e002057
#define MASK_VCOMPRESSVM    0xfe00707f

#define MATCH_VMV1RV 0x9e003057
#define MASK_VMV1RV  0xfe0ff07f
#define MATCH_VMV2RV 0x9e00b057
#define MASK_VMV2RV  0xfe0ff07f
#define MATCH_VMV4RV 0x9e01b057
#define MASK_VMV4RV  0xfe0ff07f
#define MATCH_VMV8RV 0x9e03b057
#define MASK_VMV8RV  0xfe0ff07f

#define MATCH_VDOTVV    0xe4000057
#define MASK_VDOTVV     0xfc00707f
#define MATCH_VDOTUVV   0xe0000057
#define MASK_VDOTUVV    0xfc00707f
#define MATCH_VFDOTVV   0xe4001057
#define MASK_VFDOTVV    0xfc00707f
/* END RVV */

#define MATCH_CUSTOM0 0xb
#define MASK_CUSTOM0  0x707f
#define MATCH_CUSTOM0_RS1 0x200b
#define MASK_CUSTOM0_RS1  0x707f
#define MATCH_CUSTOM0_RS1_RS2 0x300b
#define MASK_CUSTOM0_RS1_RS2  0x707f
#define MATCH_CUSTOM0_RD 0x400b
#define MASK_CUSTOM0_RD  0x707f
#define MATCH_CUSTOM0_RD_RS1 0x600b
#define MASK_CUSTOM0_RD_RS1  0x707f
#define MATCH_CUSTOM0_RD_RS1_RS2 0x700b
#define MASK_CUSTOM0_RD_RS1_RS2  0x707f
#define MATCH_CUSTOM1 0x2b
#define MASK_CUSTOM1  0x707f
#define MATCH_CUSTOM1_RS1 0x202b
#define MASK_CUSTOM1_RS1  0x707f
#define MATCH_CUSTOM1_RS1_RS2 0x302b
#define MASK_CUSTOM1_RS1_RS2  0x707f
#define MATCH_CUSTOM1_RD 0x402b
#define MASK_CUSTOM1_RD  0x707f
#define MATCH_CUSTOM1_RD_RS1 0x602b
#define MASK_CUSTOM1_RD_RS1  0x707f
#define MATCH_CUSTOM1_RD_RS1_RS2 0x702b
#define MASK_CUSTOM1_RD_RS1_RS2  0x707f
#define MATCH_CUSTOM2 0x5b
#define MASK_CUSTOM2  0x707f
#define MATCH_CUSTOM2_RS1 0x205b
#define MASK_CUSTOM2_RS1  0x707f
#define MATCH_CUSTOM2_RS1_RS2 0x305b
#define MASK_CUSTOM2_RS1_RS2  0x707f
#define MATCH_CUSTOM2_RD 0x405b
#define MASK_CUSTOM2_RD  0x707f
#define MATCH_CUSTOM2_RD_RS1 0x605b
#define MASK_CUSTOM2_RD_RS1  0x707f
#define MATCH_CUSTOM2_RD_RS1_RS2 0x705b
#define MASK_CUSTOM2_RD_RS1_RS2  0x707f
#define MATCH_CUSTOM3 0x7b
#define MASK_CUSTOM3  0x707f
#define MATCH_CUSTOM3_RS1 0x207b
#define MASK_CUSTOM3_RS1  0x707f
#define MATCH_CUSTOM3_RS1_RS2 0x307b
#define MASK_CUSTOM3_RS1_RS2  0x707f
#define MATCH_CUSTOM3_RD 0x407b
#define MASK_CUSTOM3_RD  0x707f
#define MATCH_CUSTOM3_RD_RS1 0x607b
#define MASK_CUSTOM3_RD_RS1  0x707f
#define MATCH_CUSTOM3_RD_RS1_RS2 0x707b
#define MASK_CUSTOM3_RD_RS1_RS2  0x707f

/* THEAD defines.  */
#define MATCH_DCACHE_CALL 0x0010000b
#define MASK_DCACHE_CALL 0xffffffff
#define MATCH_DCACHE_IALL 0x0020000b
#define MASK_DCACHE_IALL 0xffffffff
#define MATCH_DCACHE_CSW 0x0210000b
#define MASK_DCACHE_CSW 0xfff07fff
#define MATCH_DCACHE_ISW 0x0220000b
#define MASK_DCACHE_ISW 0xfff07fff
#define MATCH_DCACHE_CIALL 0x0030000b
#define MASK_DCACHE_CIALL 0xffffffff
#define MATCH_DCACHE_CISW 0x0230000b
#define MASK_DCACHE_CISW 0xfff07fff
#define MATCH_DCACHE_CVAL1 0x0240000b
#define MASK_DCACHE_CVAL1 0xfff07fff
#define MATCH_DCACHE_CVA 0x0250000b
#define MASK_DCACHE_CVA 0xfff07fff
#define MATCH_DCACHE_IVA 0x0260000b
#define MASK_DCACHE_IVA 0xfff07fff
#define MATCH_DCACHE_CIVA 0x0270000b
#define MASK_DCACHE_CIVA 0xfff07fff
#define MATCH_DCACHE_CPAL1 0x0280000b
#define MASK_DCACHE_CPAL1 0xfff07fff
#define MATCH_DCACHE_CPA 0x0290000b
#define MASK_DCACHE_CPA 0xfff07fff
#define MATCH_DCACHE_IPA 0x02a0000b
#define MASK_DCACHE_IPA 0xfff07fff
#define MATCH_DCACHE_CIPA 0x02b0000b
#define MASK_DCACHE_CIPA 0xfff07fff
#define MATCH_ICACHE_IALL 0x0100000b
#define MASK_ICACHE_IALL 0xffffffff
#define MATCH_ICACHE_IALLS 0x0110000b
#define MASK_ICACHE_IALLS 0xffffffff
#define MATCH_ICACHE_IVA 0x0300000b
#define MASK_ICACHE_IVA 0xfff07fff
#define MATCH_ICACHE_IPA 0x0380000b
#define MASK_ICACHE_IPA 0xfff07fff
#define MATCH_L2CACHE_CALL 0x0150000b
#define MASK_L2CACHE_CALL 0xffffffff
#define MATCH_L2CACHE_IALL 0x0160000b
#define MASK_L2CACHE_IALL 0xffffffff
#define MATCH_L2CACHE_CIALL 0x0170000b
#define MASK_L2CACHE_CIALL 0xffffffff
#define MATCH_SYNC 0x0180000b
#define MASK_SYNC 0xffffffff
#define MATCH_SYNC_S 0x0190000b
#define MASK_SYNC_S 0xffffffff
#define MATCH_SYNC_I 0x01a0000b
#define MASK_SYNC_I 0xffffffff
#define MATCH_SYNC_IS 0x01b0000b
#define MASK_SYNC_IS 0xffffffff
#define MATCH_SFENCE_VMAS 0x0400000b
#define MASK_SFENCE_VMAS 0xfe007fff
#define MATCH_TSTNBZ 0x8000100b
#define MASK_TSTNBZ 0xfff0707f
#define MATCH_MVEQZ 0x4000100b
#define MASK_MVEQZ 0xfe00707f
#define MATCH_MVNEZ 0x4200100b
#define MASK_MVNEZ 0xfe00707f
#define MATCH_MULA 0x2000100b
#define MASK_MULA 0xfe00707f
#define MATCH_MULS 0x2200100b
#define MASK_MULS 0xfe00707f
#define MATCH_MULAW 0x2400100b
#define MASK_MULAW 0xfe00707f
#define MATCH_MULSW 0x2600100b
#define MASK_MULSW 0xfe00707f
#define MATCH_MULAH 0x2800100b
#define MASK_MULAH 0xfe00707f
#define MATCH_MULSH 0x2a00100b
#define MASK_MULSH 0xfe00707f
#define MATCH_EXT 0x0000200b
#define MASK_EXT 0x0000707f
#define MATCH_EXTU 0x0000300b
#define MASK_EXTU 0x0000707f
#define MATCH_LRB 0x0000400b
#define MASK_LRB 0xf800707f
#define MATCH_LRH 0x2000400b
#define MASK_LRH 0xf800707f
#define MATCH_LRW 0x4000400b
#define MASK_LRW 0xf800707f
#define MATCH_LRD 0x6000400b
#define MASK_LRD 0xf800707f
#define MATCH_LRBU 0x8000400b
#define MASK_LRBU 0xf800707f
#define MATCH_LRHU 0xa000400b
#define MASK_LRHU 0xf800707f
#define MATCH_LRWU 0xc000400b
#define MASK_LRWU 0xf800707f
#define MATCH_LURB 0x1000400b
#define MASK_LURB 0xf800707f
#define MATCH_LURH 0x3000400b
#define MASK_LURH 0xf800707f
#define MATCH_LURW 0x5000400b
#define MASK_LURW 0xf800707f
#define MATCH_LURD 0x7000400b
#define MASK_LURD 0xf800707f
#define MATCH_LURBU 0x9000400b
#define MASK_LURBU 0xf800707f
#define MATCH_LURHU 0xb000400b
#define MASK_LURHU 0xf800707f
#define MATCH_LURWU 0xd000400b
#define MASK_LURWU 0xf800707f
#define MATCH_REV 0x8200100b
#define MASK_REV 0xfff0707f
#define MATCH_FF0 0x8400100b
#define MASK_FF0 0xfff0707f
#define MATCH_FF1 0x8600100b
#define MASK_FF1 0xfff0707f
#define MATCH_SRB 0x0000500b
#define MASK_SRB 0xf800707f
#define MATCH_SRH 0x2000500b
#define MASK_SRH 0xf800707f
#define MATCH_SRW 0x4000500b
#define MASK_SRW 0xf800707f
#define MATCH_SRD 0x6000500b
#define MASK_SRD 0xf800707f
#define MATCH_SURB 0x1000500b
#define MASK_SURB 0xf800707f
#define MATCH_SURH 0x3000500b
#define MASK_SURH 0xf800707f
#define MATCH_SURW 0x5000500b
#define MASK_SURW 0xf800707f
#define MATCH_SURD 0x7000500b
#define MASK_SURD 0xf800707f
#define MATCH_TST 0x8800100b
#define MASK_TST 0xfc00707f
#define MATCH_SRRIW 0x1400100b
#define MASK_SRRIW 0xfe00707f
#define MATCH_SRRI 0x1000100b
#define MASK_SRRI 0xfc00707f
#define MATCH_ADDSL 0x0000100b
#define MASK_ADDSL 0xf800707f
#define MATCH_SWD 0xe000500b
#define MASK_SWD 0xf800707f
#define MATCH_SDD 0xf800500b
#define MASK_SDD 0xf800707f
#define MATCH_SDIA 0x7800500b
#define MASK_SDIA 0xf800707f
#define MATCH_SDIB 0x6800500b
#define MASK_SDIB 0xf800707f
#define MATCH_SWIA 0x5800500b
#define MASK_SWIA 0xf800707f
#define MATCH_SWIB 0x4800500b
#define MASK_SWIB 0xf800707f
#define MATCH_SHIB 0x2800500b
#define MASK_SHIB 0xf800707f
#define MATCH_SHIA 0x3800500b
#define MASK_SHIA 0xf800707f
#define MATCH_SBIA 0x1800500b
#define MASK_SBIA 0xf800707f
#define MATCH_SBIB 0x0800500b
#define MASK_SBIB 0xf800707f
#define MATCH_LWUD 0xf000400b
#define MASK_LWUD 0xf800707f
#define MATCH_LWD 0xe000400b
#define MASK_LWD 0xf800707f
#define MATCH_LDD 0xf800400b
#define MASK_LDD 0xf800707f
#define MATCH_LWUIA 0xd800400b
#define MASK_LWUIA 0xf800707f
#define MATCH_LWUIB 0xc800400b
#define MASK_LWUIB 0xf800707f
#define MATCH_LHUIA 0xb800400b
#define MASK_LHUIA 0xf800707f
#define MATCH_LHUIB 0xa800400b
#define MASK_LHUIB 0xf800707f
#define MATCH_LBUIA 0x9800400b
#define MASK_LBUIA 0xf800707f
#define MATCH_LBUIB 0x8800400b
#define MASK_LBUIB 0xf800707f
#define MATCH_LDIA 0x7800400b
#define MASK_LDIA 0xf800707f
#define MATCH_LDIB 0x6800400b
#define MASK_LDIB 0xf800707f
#define MATCH_LWIA 0x5800400b
#define MASK_LWIA 0xf800707f
#define MATCH_LWIB 0x4800400b
#define MASK_LWIB 0xf800707f
#define MATCH_LHIA 0x3800400b
#define MASK_LHIA 0xf800707f
#define MATCH_LHIB 0x2800400b
#define MASK_LHIB 0xf800707f
#define MATCH_LBIA 0x1800400b
#define MASK_LBIA 0xf800707f
#define MATCH_LBIB 0x0800400b
#define MASK_LBIB 0xf800707f
#define MATCH_REVW 0x9000100b
#define MASK_REVW 0xfff0707f
#define MATCH_FSURD 0x7000700b
#define MASK_FSURD 0xf800707f
#define MATCH_FSURW 0x5000700b
#define MASK_FSURW 0xf800707f
#define MATCH_FSRD 0x6000700b
#define MASK_FSRD 0xf800707f
#define MATCH_FSRW 0x4000700b
#define MASK_FSRW 0xf800707f
#define MATCH_FLURD 0x7000600b
#define MASK_FLURD 0xf800707f
#define MATCH_FLURW 0x5000600b
#define MASK_FLURW 0xf800707f
#define MATCH_FLRD 0x6000600b
#define MASK_FLRD 0xf800707f
#define MATCH_FLRW 0x4000600b
#define MASK_FLRW 0xf800707f
#define MATCH_IPUSH 0x0040000b
#define MASK_IPUSH 0xffffffff
#define MATCH_IPOP 0x0050000b
#define MASK_IPOP 0xffffffff

/* H extensions.  */
#define MATCH_HFENCE_VVMA 0x22000073
#define MASK_HFENCE_VVMA 0xfe007fff
#define MATCH_HFENCE_GVMA 0x62000073
#define MASK_HFENCE_GVMA 0xfe007fff

/* T-Head crypto extensions.  */
#define MATCH_THEAD_ANDN 0x0800100b
#define MASK_THEAD_ANDN 0xfe00707f
#define MATCH_THEAD_ORN 0x0a00100b
#define MASK_THEAD_ORN 0xfe00707f
#define MATCH_XORN 0x0c00100b
#define MASK_XORN 0xfe00707f
#define MATCH_PACKL 0x1800100b
#define MASK_PACKL 0xfe00707f
#define MATCH_PACKH 0x1a00100b
#define MASK_PACKH 0xfe00707f
#define MATCH_PACKHL 0x1c00100b
#define MASK_PACKHL 0xfe00707f

/* T-HEAD security.  */
#define MATCH_WSC 0xcff01073
#define MASK_WSC  0xffffffff

/* T-HEAD Float for rv32.  */
#define MATCH_FMV_X_HW 0xc000100b
#define MASK_FMV_X_HW  0xfff0707f
#define MATCH_FMV_HW_X 0xa000100b
#define MASK_FMV_HW_X  0xfff0707f

/* T-HEAD FP16.  */
#define MATCH_FSH 0x1027
#define MASK_FSH  0x707f
#define MATCH_FLH 0x1007
#define MASK_FLH  0x707f
#define MATCH_FADD_H 0x04000053
#define MASK_FADD_H  0xfe00007f
#define MATCH_FSUB_H 0xc000053
#define MASK_FSUB_H  0xfe00007f
#define MATCH_FMUL_H 0x14000053
#define MASK_FMUL_H  0xfe00007f
#define MATCH_FDIV_H 0x1c000053
#define MASK_FDIV_H  0xfe00007f
#define MATCH_FSGNJ_H 0x24000053
#define MASK_FSGNJ_H  0xfe00707f
#define MATCH_FSGNJN_H 0x24001053
#define MASK_FSGNJN_H  0xfe00707f
#define MATCH_FSGNJX_H 0x24002053
#define MASK_FSGNJX_H  0xfe00707f
#define MATCH_FMIN_H 0x2c000053
#define MASK_FMIN_H  0xfe00707f
#define MATCH_FMAX_H 0x2c001053
#define MASK_FMAX_H  0xfe00707f
#define MATCH_FCVT_H_S 0x44000053
#define MASK_FCVT_H_S  0xfff0007f
#define MATCH_FCVT_S_H 0x40200053
#define MASK_FCVT_S_H  0xfff0007f
#define MATCH_FCVT_H_D 0x44100053
#define MASK_FCVT_H_D  0xfff0007f
#define MATCH_FCVT_D_H 0x42200053
#define MASK_FCVT_D_H  0xfff0007f
#define MATCH_FSQRT_H 0x5c000053
#define MASK_FSQRT_H  0xfff0007f
#define MATCH_FLE_H 0xa4000053
#define MASK_FLE_H  0xfe00707f
#define MATCH_FLT_H 0xa4001053
#define MASK_FLT_H  0xfe00707f
#define MATCH_FEQ_H 0xa4002053
#define MASK_FEQ_H  0xfe00707f
#define MATCH_FCVT_W_H 0xc4000053
#define MASK_FCVT_W_H  0xfff0007f
#define MATCH_FCVT_WU_H 0xc4100053
#define MASK_FCVT_WU_H  0xfff0007f
#define MATCH_FCVT_L_H 0xc4200053
#define MASK_FCVT_L_H  0xfff0007f
#define MATCH_FCVT_LU_H 0xc4300053
#define MASK_FCVT_LU_H  0xfff0007f
#define MATCH_FMV_X_H 0xe4000053
#define MASK_FMV_X_H  0xfff0707f
#define MATCH_FMV_X_HW 0xc000100b
#define MASK_FMV_X_HW  0xfff0707f
#define MATCH_FMV_HW_X 0xa000100b
#define MASK_FMV_HW_X  0xfff0707f
#define MATCH_FCLASS_H 0xe4001053
#define MASK_FCLASS_H  0xfff0707f
#define MATCH_FCVT_H_W 0xd4000053
#define MASK_FCVT_H_W  0xfff0007f
#define MATCH_FCVT_H_WU 0xd4100053
#define MASK_FCVT_H_WU  0xfff0007f
#define MATCH_FCVT_H_L 0xd4200053
#define MASK_FCVT_H_L  0xfff0007f
#define MATCH_FCVT_H_LU 0xd4300053
#define MASK_FCVT_H_LU  0xfff0007f
#define MATCH_FMV_H_X 0xf4000053
#define MASK_FMV_H_X  0xfff0707f
#define MATCH_FMADD_H 0x04000043
#define MASK_FMADD_H  0x600007f
#define MATCH_FMSUB_H 0x04000047
#define MASK_FMSUB_H  0x600007f
#define MATCH_FNMSUB_H 0x0400004b
#define MASK_FNMSUB_H  0x600007f
#define MATCH_FNMADD_H 0x0400004f
#define MASK_FNMADD_H  0x600007f

/* T-HEAD Vector extensions.  */
/*
   int8/int4
   https://yuque.antfin-inc.com/cpu/ywshmc/pdq7mt#W83zL

   name     func6[31:26]  vm[25]  rs2[24:20]  rs1[19:15]  func3[14:12]  rd[11:7] custom-0[6:0]
   vmaqa    10000v/x      vm         rs2          rs1        110          rd       0001011
   vmaqau   10001v/x      vm         rs2          rs1        110          rd       0001011
   vmaqasu  10010v/x      vm         rs2          rs1        110          rd       0001011
   vmaqaus  100111        vm         rs2          rs1        110          rd       0001011
   vpmaqa   10000v/x      vm         rs2          rs1        111          rd       0001011
   vpmaqau  10001v/x      vm         rs2          rs1        111          rd       0001011
   vpmaqasu 10010v/x      vm         rs2          rs1        111          rd       0001011
   vpmaqaus 100111        vm         rs2          rs1        111          rd       0001011
   vpnclip  10100v/x      vm         rs2          rs1        111          rd       0001011
   vpnclipu 10101v/x      vm         rs2          rs1        111          rd       0001011
   vpwadd   10110v/x      vm         rs2          rs1        111          rd       0001011
   vpwaddu  10111v/x      vm         rs2          rs1        111          rd       0001011
 */

#define MATCH_VMAQA_VV 0x8000600b
#define MASK_VMAQA_VV 0xfc00707f
#define MATCH_VMAQA_VX 0x8400600b
#define MASK_VMAQA_VX 0xfc00707f
#define MATCH_VMAQAU_VV 0x8800600b
#define MASK_VMAQAU_VV 0xfc00707f
#define MATCH_VMAQAU_VX 0x8c00600b
#define MASK_VMAQAU_VX 0xfc00707f
#define MATCH_VMAQASU_VV 0x9000600b
#define MASK_VMAQASU_VV 0xfc00707f
#define MATCH_VMAQASU_VX 0x9400600b
#define MASK_VMAQASU_VX 0xfc00707f
#define MATCH_VMAQAUS_VX 0x9c00600b
#define MASK_VMAQAUS_VX 0xfc00707f
#define MATCH_VPMAQA_VV 0x8000700b
#define MASK_VPMAQA_VV 0xfc00707f
#define MATCH_VPMAQA_VX 0x8400700b
#define MASK_VPMAQA_VX 0xfc00707f
#define MATCH_VPMAQAU_VV 0x8800700b
#define MASK_VPMAQAU_VV 0xfc00707f
#define MATCH_VPMAQAU_VX 0x8c00700b
#define MASK_VPMAQAU_VX 0xfc00707f
#define MATCH_VPMAQASU_VV 0x9000700b
#define MASK_VPMAQASU_VV 0xfc00707f
#define MATCH_VPMAQASU_VX 0x9400700b
#define MASK_VPMAQASU_VX 0xfc00707f
#define MATCH_VPMAQAUS_VX 0x9c00700b
#define MASK_VPMAQAUS_VX 0xfc00707f
#define MATCH_VPNCLIP_WV 0xa000700b
#define MASK_VPNCLIP_WV 0xfc00707f
#define MATCH_VPNCLIP_WX 0xa400700b
#define MASK_VPNCLIP_WX 0xfc00707f
#define MATCH_VPNCLIPU_WV 0xa800700b
#define MASK_VPNCLIPU_WV 0xfc00707f
#define MATCH_VPNCLIPU_WX 0xac00700b
#define MASK_VPNCLIPU_WX 0xfc00707f
#define MATCH_VPWADD_VV 0xb000700b
#define MASK_VPWADD_VV 0xfc00707f
#define MATCH_VPWADD_VX 0xb400700b
#define MASK_VPWADD_VX 0xfc00707f
#define MATCH_VPWADDU_VV 0xb800700b
#define MASK_VPWADDU_VV 0xfc00707f
#define MATCH_VPWADDU_VX 0xbc00700b
#define MASK_VPWADDU_VX 0xfc00707f

/* RV P  */
/* SMID 16-bit Add/Sub.  */
#define MATCH_ADD16 0x40000077
#define MASK_ADD16 0xfe00707f
#define MATCH_RADD16 0x00000077
#define MASK_RADD16 0xfe00707f
#define MATCH_URADD16 0x20000077
#define MASK_URADD16 0xfe00707f
#define MATCH_KADD16 0x10000077
#define MASK_KADD16 0xfe00707f
#define MATCH_UKADD16 0x30000077
#define MASK_UKADD16 0xfe00707f
#define MATCH_SUB16 0x42000077
#define MASK_SUB16 0xfe00707f
#define MATCH_RSUB16 0x02000077
#define MASK_RSUB16 0xfe00707f
#define MATCH_URSUB16 0x22000077
#define MASK_URSUB16 0xfe00707f
#define MATCH_KSUB16 0x12000077
#define MASK_KSUB16 0xfe00707f
#define MATCH_UKSUB16 0x32000077
#define MASK_UKSUB16 0xfe00707f
#define MATCH_CRAS16 0x44000077
#define MASK_CRAS16 0xfe00707f
#define MATCH_RCRAS16 0x04000077
#define MASK_RCRAS16 0xfe00707f
#define MATCH_URCRAS16 0x24000077
#define MASK_URCRAS16 0xfe00707f
#define MATCH_KCRAS16 0x14000077
#define MASK_KCRAS16 0xfe00707f
#define MATCH_UKCRAS16 0x34000077
#define MASK_UKCRAS16 0xfe00707f
#define MATCH_CRSA16 0x46000077
#define MASK_CRSA16 0xfe00707f
#define MATCH_RCRSA16 0x06000077
#define MASK_RCRSA16 0xfe00707f
#define MATCH_URCRSA16 0x26000077
#define MASK_URCRSA16 0xfe00707f
#define MATCH_KCRSA16 0x16000077
#define MASK_KCRSA16 0xfe00707f
#define MATCH_UKCRSA16 0x36000077
#define MASK_UKCRSA16 0xfe00707f
#define MATCH_STAS16 0xf4002077
#define MASK_STAS16 0xfe00707f
#define MATCH_RSTAS16 0xb4002077
#define MASK_RSTAS16 0xfe00707f
#define MATCH_URSTAS16 0xd4002077
#define MASK_URSTAS16 0xfe00707f
#define MATCH_KSTAS16 0xc4002077
#define MASK_KSTAS16 0xfe00707f
#define MATCH_UKSTAS16 0xe4002077
#define MASK_UKSTAS16 0xfe00707f
#define MATCH_STSA16 0xf6002077
#define MASK_STSA16 0xfe00707f
#define MATCH_RSTSA16 0xb6002077
#define MASK_RSTSA16 0xfe00707f
#define MATCH_URSTSA16 0xd6002077
#define MASK_URSTSA16 0xfe00707f
#define MATCH_KSTSA16 0xc6002077
#define MASK_KSTSA16 0xfe00707f
#define MATCH_UKSTSA16 0xe6002077
#define MASK_UKSTSA16 0xfe00707f
/* SMID 8-bit Add/Sub.  */
#define MATCH_ADD8 0x48000077
#define MASK_ADD8 0xfe00707f
#define MATCH_RADD8 0x08000077
#define MASK_RADD8 0xfe00707f
#define MATCH_URADD8 0x28000077
#define MASK_URADD8 0xfe00707f
#define MATCH_KADD8 0x18000077
#define MASK_KADD8 0xfe00707f
#define MATCH_UKADD8 0x38000077
#define MASK_UKADD8 0xfe00707f
#define MATCH_SUB8 0x4a000077
#define MASK_SUB8 0xfe00707f
#define MATCH_RSUB8 0x0a000077
#define MASK_RSUB8 0xfe00707f
#define MATCH_URSUB8 0x2a000077
#define MASK_URSUB8 0xfe00707f
#define MATCH_KSUB8 0x1a000077
#define MASK_KSUB8 0xfe00707f
#define MATCH_UKSUB8 0x3a000077
#define MASK_UKSUB8 0xfe00707f
/* SMID 16-bit shift.  */
#define MATCH_SRA16 0x50000077
#define MASK_SRA16 0xfe00707f
#define MATCH_SRAI16 0x70000077
#define MASK_SRAI16 0xff00707f
#define MATCH_SRA16_U 0x60000077
#define MASK_SRA16_U 0xfe00707f
#define MATCH_SRAI16_U 0x71000077
#define MASK_SRAI16_U 0xff00707f
#define MATCH_SRL16 0x52000077
#define MASK_SRL16 0xfe00707f
#define MATCH_SRLI16 0x72000077
#define MASK_SRLI16 0xff00707f
#define MATCH_SRL16_U 0x62000077
#define MASK_SRL16_U 0xfe00707f
#define MATCH_SRLI16_U 0x73000077
#define MASK_SRLI16_U 0xff00707f
#define MATCH_SLL16 0x54000077
#define MASK_SLL16 0xfe00707f
#define MATCH_SLLI16 0x74000077
#define MASK_SLLI16 0xff00707f
#define MATCH_KSLL16 0x64000077
#define MASK_KSLL16 0xfe00707f
#define MATCH_KSLLI16 0x75000077
#define MASK_KSLLI16 0xff00707f
#define MATCH_KSLRA16 0x56000077
#define MASK_KSLRA16 0xfe00707f
#define MATCH_KSLRA16_U 0x66000077
#define MASK_KSLRA16_U 0xfe00707f
/* SMID 8-bit shift.  */
#define MATCH_SRA8 0x58000077
#define MASK_SRA8 0xfe00707f
#define MATCH_SRAI8 0x78000077
#define MASK_SRAI8 0xff80707f
#define MATCH_SRA8_U 0x68000077
#define MASK_SRA8_U 0xfe00707f
#define MATCH_SRAI8_U 0x78800077
#define MASK_SRAI8_U 0xff80707f
#define MATCH_SRL8 0x5a000077
#define MASK_SRL8 0xfe00707f
#define MATCH_SRLI8 0x7a000077
#define MASK_SRLI8 0xff80707f
#define MATCH_SRL8_U 0x6a000077
#define MASK_SRL8_U 0xfe00707f
#define MATCH_SRLI8_U 0x7a800077
#define MASK_SRLI8_U 0xff80707f
#define MATCH_SLL8 0x5c000077
#define MASK_SLL8 0xfe00707f
#define MATCH_SLLI8 0x7c000077
#define MASK_SLLI8 0xff80707f
#define MATCH_KSLL8 0x6c000077
#define MASK_KSLL8 0xfe00707f
#define MATCH_KSLLI8 0x7c800077
#define MASK_KSLLI8 0xff80707f
#define MATCH_KSLRA8 0x5e000077
#define MASK_KSLRA8 0xfe00707f
#define MATCH_KSLRA8_U 0x6e000077
#define MASK_KSLRA8_U 0xfe00707f
/* SIMD 16-bit Compare.  */
#define MATCH_CMPEQ16 0x4c000077
#define MASK_CMPEQ16 0xfe00707f
#define MATCH_SCMPLT16 0x0c000077
#define MASK_SCMPLT16 0xfe00707f
#define MATCH_SCMPLE16 0x1c000077
#define MASK_SCMPLE16 0xfe00707f
#define MATCH_UCMPLT16 0x2c000077
#define MASK_UCMPLT16 0xfe00707f
#define MATCH_UCMPLE16 0x3c000077
#define MASK_UCMPLE16 0xfe00707f
/* SIMD 8-bit Compare.  */
#define MATCH_CMPEQ8 0x4e000077
#define MASK_CMPEQ8 0xfe00707f
#define MATCH_SCMPLT8 0x0e000077
#define MASK_SCMPLT8 0xfe00707f
#define MATCH_SCMPLE8 0x1e000077
#define MASK_SCMPLE8 0xfe00707f
#define MATCH_UCMPLT8 0x2e000077
#define MASK_UCMPLT8 0xfe00707f
#define MATCH_UCMPLE8 0x3e000077
#define MASK_UCMPLE8 0xfe00707f
/* SIMD 16-bit Multiply.  */
#define MATCH_SMUL16 0xa0000077
#define MASK_SMUL16 0xfe00707f
#define MATCH_SMULX16 0xa2000077
#define MASK_SMULX16 0xfe00707f
#define MATCH_UMUL16 0xb0000077
#define MASK_UMUL16 0xfe00707f
#define MATCH_UMULX16 0xb2000077
#define MASK_UMULX16 0xfe00707f
#define MATCH_KHM16 0x86000077
#define MASK_KHM16 0xfe00707f
#define MATCH_KHMX16 0x96000077
#define MASK_KHMX16 0xfe00707f
/* SIMD 8-bit Multiply.  */
#define MATCH_SMUL8 0xa8000077
#define MASK_SMUL8 0xfe00707f
#define MATCH_SMULX8 0xaa000077
#define MASK_SMULX8 0xfe00707f
#define MATCH_UMUL8 0xb8000077
#define MASK_UMUL8 0xfe00707f
#define MATCH_UMULX8 0xba000077
#define MASK_UMULX8 0xfe00707f
#define MATCH_KHM8 0x8e000077
#define MASK_KHM8 0xfe00707f
#define MATCH_KHMX8 0x9e000077
#define MASK_KHMX8 0xfe00707f
/* SIMD 16-bit Miscellaneous.  */
#define MATCH_SMIN16 0x80000077
#define MASK_SMIN16 0xfe00707f
#define MATCH_UMIN16 0x90000077
#define MASK_UMIN16 0xfe00707f
#define MATCH_SMAX16 0x82000077
#define MASK_SMAX16 0xfe00707f
#define MATCH_UMAX16 0x92000077
#define MASK_UMAX16 0xfe00707f
#define MATCH_SCLIP16 0x84000077
#define MASK_SCLIP16 0xff00707f
#define MATCH_UCLIP16 0x85000077
#define MASK_UCLIP16 0xff00707f
#define MATCH_KABS16 0xad100077
#define MASK_KABS16 0xfff0707f
#define MATCH_CLRS16 0xae800077
#define MASK_CLRS16 0xfff0707f
#define MATCH_CLZ16 0xae900077
#define MASK_CLO16 0xfff0707f
#define MATCH_CLO16 0xaeb00077
#define MASK_CLZ16 0xfff0707f
/* SIMD 8-bit Miscellaneous.  */
#define MATCH_SMIN8 0x88000077
#define MASK_SMIN8 0xfe00707f
#define MATCH_UMIN8 0x98000077
#define MASK_UMIN8 0xfe00707f
#define MATCH_SMAX8 0x8a000077
#define MASK_SMAX8 0xfe00707f
#define MATCH_UMAX8 0x9a000077
#define MASK_UMAX8 0xfe00707f
#define MATCH_SCLIP8 0x8c000077
#define MASK_SCLIP8 0xff80707f
#define MATCH_UCLIP8 0x8d000077
#define MASK_UCLIP8 0xff80707f
#define MATCH_KABS8 0xad000077
#define MASK_KABS8 0xfff0707f
#define MATCH_CLRS8 0xae000077
#define MASK_CLRS8 0xfff0707f
#define MATCH_CLZ8 0xae100077
#define MASK_CLO8 0xfff0707f
#define MATCH_CLO8 0xae300077
#define MASK_CLZ8 0xfff0707f
/* SIMD 8-bit unpacking.  */
#define MATCH_SUNPKD810 0xac800077
#define MASK_SUNPKD810 0xfff0707f
#define MATCH_SUNPKD820 0xac900077
#define MASK_SUNPKD820 0xfff0707f
#define MATCH_SUNPKD830 0xaca00077
#define MASK_SUNPKD830 0xfff0707f
#define MATCH_SUNPKD831 0xacb00077
#define MASK_SUNPKD831 0xfff0707f
#define MATCH_SUNPKD832 0xad300077
#define MASK_SUNPKD832 0xfff0707f
#define MATCH_ZUNPKD810 0xacc00077
#define MASK_ZUNPKD810 0xfff0707f
#define MATCH_ZUNPKD820 0xacd00077
#define MASK_ZUNPKD820 0xfff0707f
#define MATCH_ZUNPKD830 0xace00077
#define MASK_ZUNPKD830 0xfff0707f
#define MATCH_ZUNPKD831 0xacf00077
#define MASK_ZUNPKD831 0xfff0707f
#define MATCH_ZUNPKD832 0xad700077
#define MASK_ZUNPKD832 0xfff0707f
#define MATCH_PKBB16 0x0e001077
#define MASK_PKBB16 0xfe00707f
#define MATCH_PKBT16 0x1e001077
#define MASK_PKBT16 0xfe00707f
#define MATCH_PKTT16 0x2e001077
#define MASK_PKTT16 0xfe00707f
#define MATCH_PKTB16 0x3e001077
#define MASK_PKTB16 0xfe00707f
/* Signed MSW 32x32 Multiply and Add .  */
#define MATCH_SMMUL 0x40001077
#define MASK_SMMUL 0xfe00707f
#define MATCH_SMMUL_U 0x50001077
#define MASK_SMMUL_U 0xfe00707f
#define MATCH_KMMAC 0x60001077
#define MASK_KMMAC 0xfe00707f
#define MATCH_KMMAC_U 0x70001077
#define MASK_KMMAC_U 0xfe00707f
#define MATCH_KMMSB 0x42001077
#define MASK_KMMSB 0xfe00707f
#define MATCH_KMMSB_U 0x52001077
#define MASK_KMMSB_U 0xfe00707f
#define MATCH_KWMMUL 0x62001077
#define MASK_KWMMUL 0xfe00707f
#define MATCH_KWMMUL_U 0x72001077
#define MASK_KWMMUL_U 0xfe00707f
/* Signed MSW 32x12 Multiply and Add .  */
#define MATCH_SMMWB 0x44001077
#define MASK_SMMWB 0xfe00707f
#define MATCH_SMMWB_U 0x54001077
#define MASK_SMMWB_U 0xfe00707f
#define MATCH_SMMWT 0x64001077
#define MASK_SMMWT 0xfe00707f
#define MATCH_SMMWT_U 0x74001077
#define MASK_SMMWT_U 0xfe00707f
#define MATCH_KMMAWB 0x46001077
#define MASK_KMMAWB 0xfe00707f
#define MATCH_KMMAWB_U 0x56001077
#define MASK_KMMAWB_U 0xfe00707f
#define MATCH_KMMAWT 0x66001077
#define MASK_KMMAWT 0xfe00707f
#define MATCH_KMMAWT_U 0x76001077
#define MASK_KMMAWT_U 0xfe00707f
#define MATCH_KMMWB2 0x8e001077
#define MASK_KMMWB2 0xfe00707f
#define MATCH_KMMWB2_U 0x9e001077
#define MASK_KMMWB2_U 0xfe00707f
#define MATCH_KMMWT2 0xae001077
#define MASK_KMMWT2 0xfe00707f
#define MATCH_KMMWT2_U 0xbe001077
#define MASK_KMMWT2_U 0xfe00707f
#define MATCH_KMMAWB2 0xce001077
#define MASK_KMMAWB2 0xfe00707f
#define MATCH_KMMAWB2_U 0xde001077
#define MASK_KMMAWB2_U 0xfe00707f
#define MATCH_KMMAWT2 0xee001077
#define MASK_KMMAWT2 0xfe00707f
#define MATCH_KMMAWT2_U 0xfe001077
#define MASK_KMMAWT2_U 0xfe00707f
/* Singed 16bit multiply with 32bit add/sub.  */
#define MATCH_SMBB16 0x08001077
#define MASK_SMBB16 0xfe00707f
#define MATCH_SMBT16 0x18001077
#define MASK_SMBT16 0xfe00707f
#define MATCH_SMTT16 0x28001077
#define MASK_SMTT16 0xfe00707f
#define MATCH_KMDA 0x38001077
#define MASK_KMDA 0xfe00707f
#define MATCH_KMXDA 0x3a001077
#define MASK_KMXDA 0xfe00707f
#define MATCH_SMDS 0x58001077
#define MASK_SMDS 0xfe00707f
#define MATCH_SMDRS 0x68001077
#define MASK_SMDRS 0xfe00707f
#define MATCH_SMXDS 0x78001077
#define MASK_SMXDS 0xfe00707f
#define MATCH_KMABB 0x5a001077
#define MASK_KMABB 0xfe00707f
#define MATCH_KMABT 0x6a001077
#define MASK_KMABT 0xfe00707f
#define MATCH_KMATT 0x7a001077
#define MASK_KMATT 0xfe00707f
#define MATCH_KMADA 0x48001077
#define MASK_KMADA 0xfe00707f
#define MATCH_KMAXDA 0x4a001077
#define MASK_KMAXDA 0xfe00707f
#define MATCH_KMADS 0x5c001077
#define MASK_KMADS 0xfe00707f
#define MATCH_KMADRS 0x6c001077
#define MASK_KMADRS 0xfe00707f
#define MATCH_KMAXDS 0x7c001077
#define MASK_KMAXDS 0xfe00707f
#define MATCH_KMSDA 0x4c001077
#define MASK_KMSDA 0xfe00707f
#define MATCH_KMSXDA 0x4e001077
#define MASK_KMSXDA 0xfe00707f
/* Singed 16bit multiply with 64bit add/sub.  */
#define MATCH_SMAL 0x5e001077
#define MASK_SMAL 0xfe00707f
/* Partial-SIMD Miscellaneous */
#define MATCH_SCLIP32 0xe4000077
#define MASK_SCLIP32 0xfe00707f
#define MATCH_UCLIP32 0xf4000077
#define MASK_UCLIP32 0xfe00707f
#define MATCH_CLRS32 0xaf800077
#define MASK_CLRS32 0xfff0707f
#define MATCH_CLZ32 0xaf900077
#define MASK_CLZ32 0xfff0707f
#define MATCH_CLO32 0xafb00077
#define MASK_CLO32 0xfff0707f
#define MATCH_PBSAD 0xfc000077
#define MASK_PBSAD 0xfe00707f
#define MATCH_PBSADA 0xfe000077
#define MASK_PBSADA 0xfe00707f
/* 8bit mul32 add. */
#define MATCH_SMAQA 0xc8000077
#define MASK_SMAQA 0xfe00707f
#define MATCH_UMAQA 0xcc000077
#define MASK_UMAQA 0xfe00707f
#define MATCH_SMAQA_SU 0xca000077
#define MASK_SMAQA_SU 0xfe00707f
/* SMID 64-bit Add/Sub.  */
#define MATCH_ADD64 0xc0001077
#define MASK_ADD64 0xfe00707f
#define MATCH_RADD64 0x80001077
#define MASK_RADD64 0xfe00707f
#define MATCH_URADD64 0xa0001077
#define MASK_URADD64 0xfe00707f
#define MATCH_KADD64 0x90001077
#define MASK_KADD64 0xfe00707f
#define MATCH_UKADD64 0xb0001077
#define MASK_UKADD64 0xfe00707f
#define MATCH_SUB64 0xc2001077
#define MASK_SUB64 0xfe00707f
#define MATCH_RSUB64 0x82001077
#define MASK_RSUB64 0xfe00707f
#define MATCH_URSUB64 0xa2001077
#define MASK_URSUB64 0xfe00707f
#define MATCH_KSUB64 0x92001077
#define MASK_KSUB64 0xfe00707f
#define MATCH_UKSUB64 0xb2001077
#define MASK_UKSUB64 0xfe00707f
/* 32bit multiply 64bit add/sub.  */
#define MATCH_SMAR64 0x84001077
#define MASK_SMAR64 0xfe00707f
#define MATCH_SMSR64 0x86001077
#define MASK_SMSR64 0xfe00707f
#define MATCH_UMAR64 0xa4001077
#define MASK_UMAR64 0xfe00707f
#define MATCH_UMSR64 0xa6001077
#define MASK_UMSR64 0xfe00707f
#define MATCH_KMAR64 0x94001077
#define MASK_KMAR64 0xfe00707f
#define MATCH_KMSR64 0x96001077
#define MASK_KMSR64 0xfe00707f
#define MATCH_UKMAR64 0xb4001077
#define MASK_UKMAR64 0xfe00707f
#define MATCH_UKMSR64 0xb6001077
#define MASK_UKMSR64 0xfe00707f
/* singed 16bit multiply 64bit add/sub.  */
#define MATCH_SMALBB 0x88001077
#define MASK_SMALBB 0xfe00707f
#define MATCH_SMALBT 0x98001077
#define MASK_SMALBT 0xfe00707f
#define MATCH_SMALTT 0xa8001077
#define MASK_SMALTT 0xfe00707f
#define MATCH_SMALDA 0x8c001077
#define MASK_SMALDA 0xfe00707f
#define MATCH_SMALXDA 0x9c001077
#define MASK_SMALXDA 0xfe00707f
#define MATCH_SMALDS 0x8a001077
#define MASK_SMALDS 0xfe00707f
#define MATCH_SMALDRS 0x9a001077
#define MASK_SMALDRS 0xfe00707f
#define MATCH_SMALXDS 0xaa001077
#define MASK_SMALXDS 0xfe00707f
#define MATCH_SMSLDA 0xac001077
#define MASK_SMSLDA 0xfe00707f
#define MATCH_SMSLXDA 0xbc001077
#define MASK_SMSLXDA 0xfe00707f
/* Non SIMD Q15 saturation ALU.  */
#define MATCH_KADDH 0x04001077
#define MASK_KADDH 0xfe00707f
#define MATCH_KSUBH 0x06001077
#define MASK_KSUBH 0xfe00707f
#define MATCH_KHMBB 0x0c001077
#define MASK_KHMBB 0xfe00707f
#define MATCH_KHMBT 0x1c001077
#define MASK_KHMBT 0xfe00707f
#define MATCH_KHMTT 0x2c001077
#define MASK_KHMTT 0xfe00707f
#define MATCH_UKADDH 0x14001077
#define MASK_UKADDH 0xfe00707f
#define MATCH_UKSUBH 0x16001077
#define MASK_UKSUBH 0xfe00707f
/* Non SIMD Q31 saturation ALU.  */
#define MATCH_KADDW 0x00001077
#define MASK_KADDW 0xfe00707f
#define MATCH_UKADDW 0x10001077
#define MASK_UKADDW 0xfe00707f
#define MATCH_KSUBW 0x02001077
#define MASK_KSUBW 0xfe00707f
#define MATCH_UKSUBW 0x12001077
#define MASK_UKSUBW 0xfe00707f
#define MATCH_KDMBB 0x0a001077
#define MASK_KDMBB 0xfe00707f
#define MATCH_KDMBT 0x1a001077
#define MASK_KDMBT 0xfe00707f
#define MATCH_KDMTT 0x2a001077
#define MASK_KDMTT 0xfe00707f
#define MATCH_KSLRAW 0x6e001077
#define MASK_KSLRAW 0xfe00707f
#define MATCH_KSLRAW_U 0x7e001077
#define MASK_KSLRAW_U 0xfe00707f
#define MATCH_KSLLW 0x26001077
#define MASK_KSLLW 0xfe00707f
#define MATCH_KSLLIW 0x36001077
#define MASK_KSLLIW 0xfe00707f
#define MATCH_KDMABB 0xd2001077
#define MASK_KDMABB 0xfe00707f
#define MATCH_KDMABT 0xe2001077
#define MASK_KDMABT 0xfe00707f
#define MATCH_KDMATT 0xf2001077
#define MASK_KDMATT 0xfe00707f
#define MATCH_KABSW 0xad400077
#define MASK_KABSW 0xfff0707f
/* 32bit Computation.  */
#define MATCH_RADDW 0x20001077
#define MASK_RADDW 0xfe00707f
#define MATCH_URADDW 0x30001077
#define MASK_URADDW 0xfe00707f
#define MATCH_RSUBW 0x22001077
#define MASK_RSUBW 0xfe00707f
#define MATCH_URSUBW 0x32001077
#define MASK_URSUBW 0xfe00707f
#define MATCH_MAXW 0xf2000077
#define MASK_MAXW 0xfe00707f
#define MATCH_MINW 0xf0000077
#define MASK_MINW 0xfe00707f
#define MATCH_MULR64 0xf0001077
#define MASK_MULR64 0xfe00707f
#define MATCH_MULSR64 0xe0001077
#define MASK_MULSR64 0xfe00707f
/* Overflow flag set/clear.  */
#define MATCH_RDOV 0x00902073
#define MASK_RDOV 0xfffff07f
#define MATCH_CLROV 0x0090f073
#define MASK_CLROV 0xffffffff
#define MATCH_AVE 0xe0000077
#define MASK_AVE 0xfe00707f
#define MATCH_SRA_U 0x24001077
#define MASK_SRA_U 0xfe00707f
#define MATCH_SRAI_U 0xd4001077
#define MASK_SRAI_U 0xfc00707f
#define MATCH_BITREV 0xe6000077
#define MASK_BITREV 0xfe00707f
#define MATCH_BITREVI 0xe8000077
#define MASK_BITREVI 0xfc00707f
#define MATCH_WEXT 0xce000077
#define MASK_WEXT 0xfe00707f
#define MATCH_WEXTI 0xde000077
#define MASK_WEXTI 0xfe00707f
#define MATCH_BPICK 0x00003077
#define MASK_BPICK 0x0600707f
#define MATCH_INSB 0xac000077
#define MASK_INSB 0xff80707f
#define MATCH_MADDR32 0xc4001077
#define MASK_MADDR32 0xfe00707f
/* SMID 32-bit Add/Sub, Only in rv64.  */
#define MATCH_ADD32 0x40002077
#define MASK_ADD32 0xfe00707f
#define MATCH_RADD32 0x00002077
#define MASK_RADD32 0xfe00707f
#define MATCH_URADD32 0x20002077
#define MASK_URADD32 0xfe00707f
#define MATCH_KADD32 0x10002077
#define MASK_KADD32 0xfe00707f
#define MATCH_UKADD32 0x30002077
#define MASK_UKADD32 0xfe00707f
#define MATCH_SUB32 0x42002077
#define MASK_SUB32 0xfe00707f
#define MATCH_RSUB32 0x02002077
#define MASK_RSUB32 0xfe00707f
#define MATCH_URSUB32 0x22002077
#define MASK_URSUB32 0xfe00707f
#define MATCH_KSUB32 0x12002077
#define MASK_KSUB32 0xfe00707f
#define MATCH_UKSUB32 0x32002077
#define MASK_UKSUB32 0xfe00707f
#define MATCH_CRAS32 0x44002077
#define MASK_CRAS32 0xfe00707f
#define MATCH_RCRAS32 0x04002077
#define MASK_RCRAS32 0xfe00707f
#define MATCH_URCRAS32 0x24002077
#define MASK_URCRAS32 0xfe00707f
#define MATCH_KCRAS32 0x14002077
#define MASK_KCRAS32 0xfe00707f
#define MATCH_UKCRAS32 0x34002077
#define MASK_UKCRAS32 0xfe00707f
#define MATCH_CRSA32 0x46002077
#define MASK_CRSA32 0xfe00707f
#define MATCH_RCRSA32 0x06002077
#define MASK_RCRSA32 0xfe00707f
#define MATCH_URCRSA32 0x26002077
#define MASK_URCRSA32 0xfe00707f
#define MATCH_KCRSA32 0x16002077
#define MASK_KCRSA32 0xfe00707f
#define MATCH_UKCRSA32 0x36002077
#define MASK_UKCRSA32 0xfe00707f
#define MATCH_STAS32 0xf0002077
#define MASK_STAS32 0xfe00707f
#define MATCH_RSTAS32 0xb0002077
#define MASK_RSTAS32 0xfe00707f
#define MATCH_URSTAS32 0xd0002077
#define MASK_URSTAS32 0xfe00707f
#define MATCH_KSTAS32 0xc0002077
#define MASK_KSTAS32 0xfe00707f
#define MATCH_UKSTAS32 0xe0002077
#define MASK_UKSTAS32 0xfe00707f
#define MATCH_STSA32 0xf2002077
#define MASK_STSA32 0xfe00707f
#define MATCH_RSTSA32 0xb2002077
#define MASK_RSTSA32 0xfe00707f
#define MATCH_URSTSA32 0xd2002077
#define MASK_URSTSA32 0xfe00707f
#define MATCH_KSTSA32 0xc2002077
#define MASK_KSTSA32 0xfe00707f
#define MATCH_UKSTSA32 0xe2002077
#define MASK_UKSTSA32 0xfe00707f
/* (RV64 Only) SIMD 32-bit Shift Instructions. */
#define MATCH_SRA32 0x50002077
#define MASK_SRA32 0xfe00707f
#define MATCH_SRAI32 0x70002077
#define MASK_SRAI32 0xfe00707f
#define MATCH_SRA32_U 0x60002077
#define MASK_SRA32_U 0xfe00707f
#define MATCH_SRAI32_U 0x80002077
#define MASK_SRAI32_U 0xfe00707f
#define MATCH_SRL32 0x52002077
#define MASK_SRL32 0xfe00707f
#define MATCH_SRLI32 0x72002077
#define MASK_SRLI32 0xfe00707f
#define MATCH_SRL32_U 0x62002077
#define MASK_SRL32_U 0xfe00707f
#define MATCH_SRLI32_U 0x82002077
#define MASK_SRLI32_U 0xfe00707f
#define MATCH_SLL32 0x54002077
#define MASK_SLL32 0xfe00707f
#define MATCH_SLLI32 0x74002077
#define MASK_SLLI32 0xfe00707f
#define MATCH_KSLL32 0x64002077
#define MASK_KSLL32 0xfe00707f
#define MATCH_KSLLI32 0x84002077
#define MASK_KSLLI32 0xfe00707f
#define MATCH_KSLRA32 0x56002077
#define MASK_KSLRA32 0xfe00707f
#define MATCH_KSLRA32_U 0x66002077
#define MASK_KSLRA32_U 0xfe00707f
/* (RV64 Only) SIMD 32-bit Miscellaneous Instructions. */
#define MATCH_SMIN32 0x90002077
#define MASK_SMIN32 0xfe00707f
#define MATCH_UMIN32 0xa0002077
#define MASK_UMIN32 0xfe00707f
#define MATCH_SMAX32 0x92002077
#define MASK_SMAX32 0xfe00707f
#define MATCH_UMAX32 0xa2002077
#define MASK_UMAX32 0xfe00707f
#define MATCH_KABS32 0xad200077
#define MASK_KABS32 0xfff0707f
/* (RV64 Only) SIMD Q15 saturating Multiply Instructions.  */
#define MATCH_KHMBB16 0xdc001077
#define MASK_KHMBB16 0xfe00707f
#define MATCH_KHMBT16 0xec001077
#define MASK_KHMBT16 0xfe00707f
#define MATCH_KHMTT16 0xfc001077
#define MASK_KHMTT16 0xfe00707f
#define MATCH_KDMBB16 0xda001077
#define MASK_KDMBB16 0xfe00707f
#define MATCH_KDMBT16 0xea001077
#define MASK_KDMBT16 0xfe00707f
#define MATCH_KDMTT16 0xfa001077
#define MASK_KDMTT16 0xfe00707f
#define MATCH_KDMABB16 0xd8001077
#define MASK_KDMABB16 0xfe00707f
#define MATCH_KDMABT16 0xe8001077
#define MASK_KDMABT16 0xfe00707f
#define MATCH_KDMATT16 0xf8001077
#define MASK_KDMATT16 0xfe00707f
/* (RV64 Only) 32-bit Multiply Instructions.  */
//#define MATCH_SMBB32 0x08002077 Same as mulsr64
//#define MASK_SMBB32 0xfe00707f
#define MATCH_SMBT32 0x18002077
#define MASK_SMBT32 0xfe00707f
#define MATCH_SMTT32 0x28002077
#define MASK_SMTT32 0xfe00707f
/* (RV64 Only) 32-bit Multiply & Add Instructions.  */
#define MATCH_KMABB32 0x5a002077
#define MASK_KMABB32 0xfe00707f
#define MATCH_KMABT32 0x6a002077
#define MASK_KMABT32 0xfe00707f
#define MATCH_KMATT32 0x7a002077
#define MASK_KMATT32 0xfe00707f
/* (RV64 Only) 32-bit Parallel Multiply & Add Instructions.  */
#define MATCH_KMDA32 0x38002077
#define MASK_KMDA32 0xfe00707f
#define MATCH_KMXDA32 0x3a002077
#define MASK_KMXDA32 0xfe00707f
//#define MATCH_KMADA32 0x4a001077 no encoding, same as kmar64
//#define MASK_KMADA32 0xfe00707f
#define MATCH_KMAXDA32 0x4a002077
#define MASK_KMAXDA32 0xfe00707f
#define MATCH_KMADS32 0x5c002077
#define MASK_KMADS32 0xfe00707f
#define MATCH_KMADRS32 0x6c002077
#define MASK_KMADRS32 0xfe00707f
#define MATCH_KMAXDS32 0x7c002077
#define MASK_KMAXDS32 0xfe00707f
#define MATCH_KMSDA32 0x4c002077
#define MASK_KMSDA32 0xfe00707f
#define MATCH_KMSXDA32 0x4e002077
#define MASK_KMSXDA32 0xfe00707f
#define MATCH_SMDS32 0x58002077
#define MASK_SMDS32 0xfe00707f
#define MATCH_SMDRS32 0x68002077
#define MASK_SMDRS32 0xfe00707f
#define MATCH_SMXDS32 0x78002077
#define MASK_SMXDS32 0xfe00707f
/* (RV64 Only) Non-SIMD 32-bit Shift Instructions.  */
#define MATCH_SRAIW_U 0x34001077
#define MASK_SRAIW_U 0xfe00707f
/* 32-bit Packing Instructions. */
#define MATCH_PKBB32 0x0e002077
#define MASK_PKBB32 0xfe00707f
#define MATCH_PKBT32 0x1e002077
#define MASK_PKBT32 0xfe00707f
#define MATCH_PKTT32 0x2e002077
#define MASK_PKTT32 0xfe00707f
#define MATCH_PKTB32 0x3e002077
#define MASK_PKTB32 0xfe00707f

#define MATCH_MSUBR32 0xc6001077
#define MASK_MSUBR32 0xfe00707f
#define MATCH_SWAP8 0xad800077
#define MASK_SWAP8 0xfff0707f

/* Cache management instructions.  */
#define MATCH_CBO_CLEAN 0x10200f
#define MASK_CBO_CLEAN 0xfff07fff
#define MATCH_CBO_FLUSH 0x20200f
#define MASK_CBO_FLUSH 0xfff07fff
#define MATCH_CBO_INVAL 0x200f
#define MASK_CBO_INVAL 0xfff07fff
#define MATCH_CBO_ZERO 0x40200f
#define MASK_CBO_ZERO 0xfff07fff
/* Zicbop hint instructions. */
#define MATCH_PREFETCH_I 0x6013
#define MASK_PREFETCH_I 0x1f07fff
#define MATCH_PREFETCH_R 0x106013
#define MASK_PREFETCH_R 0x1f07fff
#define MATCH_PREFETCH_W 0x306013
#define MASK_PREFETCH_W 0x1f07fff

/* Svinval instruction.  */
#define MATCH_SINVAL_VMA 0x16000073
#define MASK_SINVAL_VMA 0xfe007fff
#define MATCH_SFENCE_W_INVAL 0x18000073
#define MASK_SFENCE_W_INVAL 0xffffffff
#define MATCH_SFENCE_INVAL_IR 0x18100073
#define MASK_SFENCE_INVAL_IR 0xffffffff
#define MATCH_HINVAL_VVMA 0x26000073
#define MASK_HINVAL_VVMA 0xfe007fff
#define MATCH_HINVAL_GVMA 0x66000073
#define MASK_HINVAL_GVMA 0xfe007fff

/* Privileged CSR addresses (v1.11).  */
#define CSR_USTATUS 0x0
#define CSR_UIE 0x4
#define CSR_UTVEC 0x5
#define CSR_USCRATCH 0x40
#define CSR_UEPC 0x41
#define CSR_UCAUSE 0x42
#define CSR_UTVAL 0x43
#define CSR_UIP 0x44
#define CSR_CYCLE 0xc00
#define CSR_TIME 0xc01
#define CSR_INSTRET 0xc02
#define CSR_HPMCOUNTER3 0xc03
#define CSR_HPMCOUNTER4 0xc04
#define CSR_HPMCOUNTER5 0xc05
#define CSR_HPMCOUNTER6 0xc06
#define CSR_HPMCOUNTER7 0xc07
#define CSR_HPMCOUNTER8 0xc08
#define CSR_HPMCOUNTER9 0xc09
#define CSR_HPMCOUNTER10 0xc0a
#define CSR_HPMCOUNTER11 0xc0b
#define CSR_HPMCOUNTER12 0xc0c
#define CSR_HPMCOUNTER13 0xc0d
#define CSR_HPMCOUNTER14 0xc0e
#define CSR_HPMCOUNTER15 0xc0f
#define CSR_HPMCOUNTER16 0xc10
#define CSR_HPMCOUNTER17 0xc11
#define CSR_HPMCOUNTER18 0xc12
#define CSR_HPMCOUNTER19 0xc13
#define CSR_HPMCOUNTER20 0xc14
#define CSR_HPMCOUNTER21 0xc15
#define CSR_HPMCOUNTER22 0xc16
#define CSR_HPMCOUNTER23 0xc17
#define CSR_HPMCOUNTER24 0xc18
#define CSR_HPMCOUNTER25 0xc19
#define CSR_HPMCOUNTER26 0xc1a
#define CSR_HPMCOUNTER27 0xc1b
#define CSR_HPMCOUNTER28 0xc1c
#define CSR_HPMCOUNTER29 0xc1d
#define CSR_HPMCOUNTER30 0xc1e
#define CSR_HPMCOUNTER31 0xc1f
#define CSR_CYCLEH 0xc80
#define CSR_TIMEH 0xc81
#define CSR_INSTRETH 0xc82
#define CSR_HPMCOUNTER3H 0xc83
#define CSR_HPMCOUNTER4H 0xc84
#define CSR_HPMCOUNTER5H 0xc85
#define CSR_HPMCOUNTER6H 0xc86
#define CSR_HPMCOUNTER7H 0xc87
#define CSR_HPMCOUNTER8H 0xc88
#define CSR_HPMCOUNTER9H 0xc89
#define CSR_HPMCOUNTER10H 0xc8a
#define CSR_HPMCOUNTER11H 0xc8b
#define CSR_HPMCOUNTER12H 0xc8c
#define CSR_HPMCOUNTER13H 0xc8d
#define CSR_HPMCOUNTER14H 0xc8e
#define CSR_HPMCOUNTER15H 0xc8f
#define CSR_HPMCOUNTER16H 0xc90
#define CSR_HPMCOUNTER17H 0xc91
#define CSR_HPMCOUNTER18H 0xc92
#define CSR_HPMCOUNTER19H 0xc93
#define CSR_HPMCOUNTER20H 0xc94
#define CSR_HPMCOUNTER21H 0xc95
#define CSR_HPMCOUNTER22H 0xc96
#define CSR_HPMCOUNTER23H 0xc97
#define CSR_HPMCOUNTER24H 0xc98
#define CSR_HPMCOUNTER25H 0xc99
#define CSR_HPMCOUNTER26H 0xc9a
#define CSR_HPMCOUNTER27H 0xc9b
#define CSR_HPMCOUNTER28H 0xc9c
#define CSR_HPMCOUNTER29H 0xc9d
#define CSR_HPMCOUNTER30H 0xc9e
#define CSR_HPMCOUNTER31H 0xc9f
#define CSR_SSTATUS 0x100
#define CSR_SEDELEG 0x102
#define CSR_SIDELEG 0x103
#define CSR_SIE 0x104
#define CSR_STVEC 0x105
#define CSR_SCOUNTEREN 0x106
#define CSR_SENVCFG 0x10a
#define CSR_SSCRATCH 0x140
#define CSR_SEPC 0x141
#define CSR_SCAUSE 0x142
#define CSR_STVAL 0x143
#define CSR_SIP 0x144
#define CSR_STIMECMP 0x14d
#define CSR_STIMECMPH 0x15d
#define CSR_SATP 0x180
#define CSR_MVENDORID 0xf11
#define CSR_MARCHID 0xf12
#define CSR_MIMPID 0xf13
#define CSR_MHARTID 0xf14
#define CSR_MCONFIGPTR 0xf15
#define CSR_MSTATUS 0x300
#define CSR_MISA 0x301
#define CSR_MEDELEG 0x302
#define CSR_MIDELEG 0x303
#define CSR_MIE 0x304
#define CSR_MTVEC 0x305
#define CSR_MCOUNTEREN 0x306
#define CSR_MTVT 0x307
#define CSR_MENVCFG 0x30a
#define CSR_MSCRATCH 0x340
#define CSR_MEPC 0x341
#define CSR_MCAUSE 0x342
#define CSR_MTVAL 0x343
#define CSR_MIP 0x344
#define CSR_MNXTI 0x345
#define CSR_MINTSTATUS 0x346
#define CSR_MSCRATCHCSW 0x348
#define CSR_MSCRATCHCSWL 0x349
#define CSR_MCLICBASE 0x350
#define CSR_PMPCFG0 0x3a0
#define CSR_PMPCFG1 0x3a1
#define CSR_PMPCFG2 0x3a2
#define CSR_PMPCFG3 0x3a3
#define CSR_PMPCFG4 0x3a4
#define CSR_PMPCFG5 0x3a5
#define CSR_PMPCFG6 0x3a6
#define CSR_PMPCFG7 0x3a7
#define CSR_PMPCFG8 0x3a8
#define CSR_PMPCFG9 0x3a9
#define CSR_PMPCFG10 0x3aa
#define CSR_PMPCFG11 0x3ab
#define CSR_PMPCFG12 0x3ac
#define CSR_PMPCFG13 0x3ad
#define CSR_PMPCFG14 0x3ae
#define CSR_PMPCFG15 0x3af
#define CSR_PMPADDR0 0x3b0
#define CSR_PMPADDR1 0x3b1
#define CSR_PMPADDR2 0x3b2
#define CSR_PMPADDR3 0x3b3
#define CSR_PMPADDR4 0x3b4
#define CSR_PMPADDR5 0x3b5
#define CSR_PMPADDR6 0x3b6
#define CSR_PMPADDR7 0x3b7
#define CSR_PMPADDR8 0x3b8
#define CSR_PMPADDR9 0x3b9
#define CSR_PMPADDR10 0x3ba
#define CSR_PMPADDR11 0x3bb
#define CSR_PMPADDR12 0x3bc
#define CSR_PMPADDR13 0x3bd
#define CSR_PMPADDR14 0x3be
#define CSR_PMPADDR15 0x3bf
#define CSR_PMPADDR16 0x3c0
#define CSR_PMPADDR17 0x3c1
#define CSR_PMPADDR18 0x3c2
#define CSR_PMPADDR19 0x3c3
#define CSR_PMPADDR20 0x3c4
#define CSR_PMPADDR21 0x3c5
#define CSR_PMPADDR22 0x3c6
#define CSR_PMPADDR23 0x3c7
#define CSR_PMPADDR24 0x3c8
#define CSR_PMPADDR25 0x3c9
#define CSR_PMPADDR26 0x3ca
#define CSR_PMPADDR27 0x3cb
#define CSR_PMPADDR28 0x3cc
#define CSR_PMPADDR29 0x3cd
#define CSR_PMPADDR30 0x3ce
#define CSR_PMPADDR31 0x3cf
#define CSR_PMPADDR32 0x3d0
#define CSR_PMPADDR33 0x3d1
#define CSR_PMPADDR34 0x3d2
#define CSR_PMPADDR35 0x3d3
#define CSR_PMPADDR36 0x3d4
#define CSR_PMPADDR37 0x3d5
#define CSR_PMPADDR38 0x3d6
#define CSR_PMPADDR39 0x3d7
#define CSR_PMPADDR40 0x3d8
#define CSR_PMPADDR41 0x3d9
#define CSR_PMPADDR42 0x3da
#define CSR_PMPADDR43 0x3db
#define CSR_PMPADDR44 0x3dc
#define CSR_PMPADDR45 0x3dd
#define CSR_PMPADDR46 0x3de
#define CSR_PMPADDR47 0x3df
#define CSR_PMPADDR48 0x3e0
#define CSR_PMPADDR49 0x3e1
#define CSR_PMPADDR50 0x3e2
#define CSR_PMPADDR51 0x3e3
#define CSR_PMPADDR52 0x3e4
#define CSR_PMPADDR53 0x3e5
#define CSR_PMPADDR54 0x3e6
#define CSR_PMPADDR55 0x3e7
#define CSR_PMPADDR56 0x3e8
#define CSR_PMPADDR57 0x3e9
#define CSR_PMPADDR58 0x3ea
#define CSR_PMPADDR59 0x3eb
#define CSR_PMPADDR60 0x3ec
#define CSR_PMPADDR61 0x3ed
#define CSR_PMPADDR62 0x3ee
#define CSR_PMPADDR63 0x3ef
#define CSR_MCYCLE 0xb00
#define CSR_MINSTRET 0xb02
#define CSR_MHPMCOUNTER3 0xb03
#define CSR_MHPMCOUNTER4 0xb04
#define CSR_MHPMCOUNTER5 0xb05
#define CSR_MHPMCOUNTER6 0xb06
#define CSR_MHPMCOUNTER7 0xb07
#define CSR_MHPMCOUNTER8 0xb08
#define CSR_MHPMCOUNTER9 0xb09
#define CSR_MHPMCOUNTER10 0xb0a
#define CSR_MHPMCOUNTER11 0xb0b
#define CSR_MHPMCOUNTER12 0xb0c
#define CSR_MHPMCOUNTER13 0xb0d
#define CSR_MHPMCOUNTER14 0xb0e
#define CSR_MHPMCOUNTER15 0xb0f
#define CSR_MHPMCOUNTER16 0xb10
#define CSR_MHPMCOUNTER17 0xb11
#define CSR_MHPMCOUNTER18 0xb12
#define CSR_MHPMCOUNTER19 0xb13
#define CSR_MHPMCOUNTER20 0xb14
#define CSR_MHPMCOUNTER21 0xb15
#define CSR_MHPMCOUNTER22 0xb16
#define CSR_MHPMCOUNTER23 0xb17
#define CSR_MHPMCOUNTER24 0xb18
#define CSR_MHPMCOUNTER25 0xb19
#define CSR_MHPMCOUNTER26 0xb1a
#define CSR_MHPMCOUNTER27 0xb1b
#define CSR_MHPMCOUNTER28 0xb1c
#define CSR_MHPMCOUNTER29 0xb1d
#define CSR_MHPMCOUNTER30 0xb1e
#define CSR_MHPMCOUNTER31 0xb1f
#define CSR_MCYCLEH 0xb80
#define CSR_MINSTRETH 0xb82
#define CSR_MHPMCOUNTER3H 0xb83
#define CSR_MHPMCOUNTER4H 0xb84
#define CSR_MHPMCOUNTER5H 0xb85
#define CSR_MHPMCOUNTER6H 0xb86
#define CSR_MHPMCOUNTER7H 0xb87
#define CSR_MHPMCOUNTER8H 0xb88
#define CSR_MHPMCOUNTER9H 0xb89
#define CSR_MHPMCOUNTER10H 0xb8a
#define CSR_MHPMCOUNTER11H 0xb8b
#define CSR_MHPMCOUNTER12H 0xb8c
#define CSR_MHPMCOUNTER13H 0xb8d
#define CSR_MHPMCOUNTER14H 0xb8e
#define CSR_MHPMCOUNTER15H 0xb8f
#define CSR_MHPMCOUNTER16H 0xb90
#define CSR_MHPMCOUNTER17H 0xb91
#define CSR_MHPMCOUNTER18H 0xb92
#define CSR_MHPMCOUNTER19H 0xb93
#define CSR_MHPMCOUNTER20H 0xb94
#define CSR_MHPMCOUNTER21H 0xb95
#define CSR_MHPMCOUNTER22H 0xb96
#define CSR_MHPMCOUNTER23H 0xb97
#define CSR_MHPMCOUNTER24H 0xb98
#define CSR_MHPMCOUNTER25H 0xb99
#define CSR_MHPMCOUNTER26H 0xb9a
#define CSR_MHPMCOUNTER27H 0xb9b
#define CSR_MHPMCOUNTER28H 0xb9c
#define CSR_MHPMCOUNTER29H 0xb9d
#define CSR_MHPMCOUNTER30H 0xb9e
#define CSR_MHPMCOUNTER31H 0xb9f
#define CSR_MCOUNTINHIBIT 0x320
#define CSR_MHPMEVENT3 0x323
#define CSR_MHPMEVENT4 0x324
#define CSR_MHPMEVENT5 0x325
#define CSR_MHPMEVENT6 0x326
#define CSR_MHPMEVENT7 0x327
#define CSR_MHPMEVENT8 0x328
#define CSR_MHPMEVENT9 0x329
#define CSR_MHPMEVENT10 0x32a
#define CSR_MHPMEVENT11 0x32b
#define CSR_MHPMEVENT12 0x32c
#define CSR_MHPMEVENT13 0x32d
#define CSR_MHPMEVENT14 0x32e
#define CSR_MHPMEVENT15 0x32f
#define CSR_MHPMEVENT16 0x330
#define CSR_MHPMEVENT17 0x331
#define CSR_MHPMEVENT18 0x332
#define CSR_MHPMEVENT19 0x333
#define CSR_MHPMEVENT20 0x334
#define CSR_MHPMEVENT21 0x335
#define CSR_MHPMEVENT22 0x336
#define CSR_MHPMEVENT23 0x337
#define CSR_MHPMEVENT24 0x338
#define CSR_MHPMEVENT25 0x339
#define CSR_MHPMEVENT26 0x33a
#define CSR_MHPMEVENT27 0x33b
#define CSR_MHPMEVENT28 0x33c
#define CSR_MHPMEVENT29 0x33d
#define CSR_MHPMEVENT30 0x33e
#define CSR_MHPMEVENT31 0x33f
#define CSR_HSTATUS 0x200
#define CSR_HEDELEG 0x202
#define CSR_HIDELEG 0x203
#define CSR_HIE 0x204
#define CSR_HTVEC 0x205
#define CSR_HSCRATCH 0x240
#define CSR_HEPC 0x241
#define CSR_HCAUSE 0x242
#define CSR_HBADADDR 0x243
#define CSR_HIP 0x244
#define CSR_MBASE 0x380
#define CSR_MBOUND 0x381
#define CSR_MIBASE 0x382
#define CSR_MIBOUND 0x383
#define CSR_MDBASE 0x384
#define CSR_MDBOUND 0x385
#define CSR_MSCOUNTEREN 0x321
#define CSR_MHCOUNTEREN 0x322


/* PRIV_SPEC_CLASS_1P12 for XUANTIE C908.  */
#define CSR_MSECCFG 0x747

/* New CSR for C908.  */
#define CSR_MCPER 0x7d9

/* T-HEAD M mode CSR.  */
#define CSR_MXSTATUS 0x7c0
#define CSR_MHCR 0x7c1
#define CSR_MCOR 0x7c2
#define CSR_MCCR2 0x7c3
#define CSR_MCER2 0x7c4
#define CSR_MHINT 0x7c5
#define CSR_MRMR 0x7c6
#define CSR_MRVBR 0x7c7
#define CSR_MCER 0x7c8
#define CSR_MCOUNTERWEN 0x7c9
#define CSR_MCOUNTERINTEN 0x7ca
#define CSR_MCOUNTEROF 0x7cb
#define CSR_MHINT2 0x7cc
#define CSR_MHINT3 0x7cd
#define CSR_MRADDR 0x7e0
#define CSR_MEXSTATUS 0x7e1
#define CSR_MNMICAUSE 0x7e2
#define CSR_MNMIPC 0x7e3
#define CSR_MHPMCR 0x7f0
#define CSR_MHPMSR 0x7f1
#define CSR_MHPMER 0x7f2
#define CSR_MSMPR 0x7f3
#define CSR_MTEECFG 0x7f4
#define CSR_MZONEID 0x7f5
#define CSR_MLLCPID 0x7f6
#define CSR_MLLWP 0x7f7
#define CSR_MDTCMCR 0x7f8
#define CSR_MITCMCR 0x7f9
#define CSR_MIESR 0x7fa
#define CSR_MSBEPA 0x7fb
#define CSR_MSBEPA2 0x7fc
#define CSR_USP 0x7d1
#define CSR_MCINS 0x7d2
#define CSR_MCINDEX 0x7d3
#define CSR_MCDATA0 0x7d4
#define CSR_MCDATA1 0x7d5
#define CSR_MEICR 0x7d6
#define CSR_MEICR2 0x7d7
#define CSR_MBEADDR 0x7d8
#define CSR_MCPUID 0xfc0
#define CSR_MAPBADDR 0xfc1
#define CSR_MWMSR 0xfc2
#define CSR_MHALTCAUSE 0xfe0
#define CSR_MDBGINFO 0xfe1
#define CSR_MPCFIFO 0xfe2

#define CSR_SCONTEXT 0x5a8

/* T-HEAD S mode CSR.  */
#define CSR_SXSTATUS 0x5c0
#define CSR_SHCR 0x5c1
#define CSR_SCER2 0x5c2
#define CSR_SCER 0x5c3
#define CSR_SCOUNTERINTEN 0x5c4
#define CSR_SCOUNTEROF 0x5c5
#define CSR_SHINT 0x5c6
#define CSR_SHINT2 0x5c7
#define CSR_SHPMINHIBIT 0x5c8
#define CSR_SHPMCR 0x5c9
#define CSR_SHPMSR 0x5ca
#define CSR_SHPMER 0x5cb
#define CSR_SLLCPID 0x5cc
#define CSR_SLLWP 0x5cd
#define CSR_SIESR 0x5ce
#define CSR_SBEADDR 0x5d0
#define CSR_SSBEPA 0x5d1
#define CSR_SSBEPA2 0x5d2
#define CSR_SCYCLE 0x5e0
#define CSR_SHPMCOUNTER1 0x5e1
#define CSR_SHPMCOUNTER2 0x5e2
#define CSR_SHPMCOUNTER3 0x5e3
#define CSR_SHPMCOUNTER4 0x5e4
#define CSR_SHPMCOUNTER5 0x5e5
#define CSR_SHPMCOUNTER6 0x5e6
#define CSR_SHPMCOUNTER7 0x5e7
#define CSR_SHPMCOUNTER8 0x5e8
#define CSR_SHPMCOUNTER9 0x5e9
#define CSR_SHPMCOUNTER10 0x5ea
#define CSR_SHPMCOUNTER11 0x5eb
#define CSR_SHPMCOUNTER12 0x5ec
#define CSR_SHPMCOUNTER13 0x5ed
#define CSR_SHPMCOUNTER14 0x5ee
#define CSR_SHPMCOUNTER15 0x5ef
#define CSR_SHPMCOUNTER16 0x5f0
#define CSR_SHPMCOUNTER17 0x5f1
#define CSR_SHPMCOUNTER18 0x5f2
#define CSR_SHPMCOUNTER19 0x5f3
#define CSR_SHPMCOUNTER20 0x5f4
#define CSR_SHPMCOUNTER21 0x5f5
#define CSR_SHPMCOUNTER22 0x5f6
#define CSR_SHPMCOUNTER23 0x5f7
#define CSR_SHPMCOUNTER24 0x5f8
#define CSR_SHPMCOUNTER25 0x5f9
#define CSR_SHPMCOUNTER26 0x5fa
#define CSR_SHPMCOUNTER27 0x5fb
#define CSR_SHPMCOUNTER28 0x5fc
#define CSR_SHPMCOUNTER29 0x5fd
#define CSR_SHPMCOUNTER30 0x5fe
#define CSR_SHPMCOUNTER31 0x5ff

/* T-HEAD U mode CSR.  */
#define CSR_FXCR 0x800

/* T-HEAD MMU extentions.  */
#define CSR_SMIR 0x9c0
#define CSR_SMEL 0x9c1
#define CSR_SMEH 0x9c2
#define CSR_SMCIR 0x9c3

/* T-HEAD Security CSR(May be droped).  */
#define CSR_MEBR 0xbe0
#define CSR_NT_MSTATUS 0xbe1
#define CSR_NT_MIE 0xbe2
#define CSR_NT_MTVEC 0xbe3
#define CSR_NT_MTVT 0xbe4
#define CSR_NT_MEPC 0xbe5
#define CSR_NT_MCAUSE 0xbe6
#define CSR_NT_MIP 0xbe7
#define CSR_NT_MINTSTATE 0xbe8
#define CSR_NT_MXSTATUS 0xbe9
#define CSR_NT_MEBR 0xbea
#define CSR_NT_MSP 0xbeb
#define CSR_T_USP 0xbec
#define CSR_T_MDCR 0xbed
#define CSR_T_MPCR 0xbee
#define CSR_PMPTEECFG 0xbef

/* Unprivileged CSR addresses.  */
#define CSR_FFLAGS 0x1
#define CSR_FRM 0x2
#define CSR_FCSR 0x3
#define CSR_VSTART 0x008
#define CSR_VXSAT 0x009
#define CSR_VXRM 0x00a
#define CSR_VCSR 0x00f
#define CSR_VL 0xc20
#define CSR_VTYPE 0xc21
#define CSR_VLENB 0xc22
#define CSR_DCSR 0x7b0
#define CSR_DPC 0x7b1
#define CSR_DSCRATCH0 0x7b2
#define CSR_DSCRATCH1 0x7b3
#define CSR_TSELECT 0x7a0
#define CSR_TDATA1 0x7a1
#define CSR_TDATA2 0x7a2
#define CSR_TDATA3 0x7a3
#define CSR_TINFO 0x7a4
#define CSR_TCONTROL 0x7a5
#define CSR_MCONTEXT 0x7a8

// NOTICE: scontext has defined in privileaged 1.12, 0x5a8.
// #define CSR_SCONTEXT 0x7aa

/* New for sstc.  */
#define CSR_STIMECMP 0x14d
#define CSR_STIMECMPH 0x15d

#endif /* RISCV_ENCODING_H.  */
#ifdef DECLARE_INSN
DECLARE_INSN(slli_rv32, MATCH_SLLI_RV32, MASK_SLLI_RV32)
DECLARE_INSN(srli_rv32, MATCH_SRLI_RV32, MASK_SRLI_RV32)
DECLARE_INSN(srai_rv32, MATCH_SRAI_RV32, MASK_SRAI_RV32)
DECLARE_INSN(frflags, MATCH_FRFLAGS, MASK_FRFLAGS)
DECLARE_INSN(fsflags, MATCH_FSFLAGS, MASK_FSFLAGS)
DECLARE_INSN(fsflagsi, MATCH_FSFLAGSI, MASK_FSFLAGSI)
DECLARE_INSN(frrm, MATCH_FRRM, MASK_FRRM)
DECLARE_INSN(fsrm, MATCH_FSRM, MASK_FSRM)
DECLARE_INSN(fsrmi, MATCH_FSRMI, MASK_FSRMI)
DECLARE_INSN(fscsr, MATCH_FSCSR, MASK_FSCSR)
DECLARE_INSN(frcsr, MATCH_FRCSR, MASK_FRCSR)
DECLARE_INSN(rdcycle, MATCH_RDCYCLE, MASK_RDCYCLE)
DECLARE_INSN(rdtime, MATCH_RDTIME, MASK_RDTIME)
DECLARE_INSN(rdinstret, MATCH_RDINSTRET, MASK_RDINSTRET)
DECLARE_INSN(rdcycleh, MATCH_RDCYCLEH, MASK_RDCYCLEH)
DECLARE_INSN(rdtimeh, MATCH_RDTIMEH, MASK_RDTIMEH)
DECLARE_INSN(rdinstreth, MATCH_RDINSTRETH, MASK_RDINSTRETH)
DECLARE_INSN(scall, MATCH_SCALL, MASK_SCALL)
DECLARE_INSN(sbreak, MATCH_SBREAK, MASK_SBREAK)
DECLARE_INSN(beq, MATCH_BEQ, MASK_BEQ)
DECLARE_INSN(bne, MATCH_BNE, MASK_BNE)
DECLARE_INSN(blt, MATCH_BLT, MASK_BLT)
DECLARE_INSN(bge, MATCH_BGE, MASK_BGE)
DECLARE_INSN(bltu, MATCH_BLTU, MASK_BLTU)
DECLARE_INSN(bgeu, MATCH_BGEU, MASK_BGEU)
DECLARE_INSN(jalr, MATCH_JALR, MASK_JALR)
DECLARE_INSN(jal, MATCH_JAL, MASK_JAL)
DECLARE_INSN(lui, MATCH_LUI, MASK_LUI)
DECLARE_INSN(auipc, MATCH_AUIPC, MASK_AUIPC)
DECLARE_INSN(addi, MATCH_ADDI, MASK_ADDI)
DECLARE_INSN(slli, MATCH_SLLI, MASK_SLLI)
DECLARE_INSN(slti, MATCH_SLTI, MASK_SLTI)
DECLARE_INSN(sltiu, MATCH_SLTIU, MASK_SLTIU)
DECLARE_INSN(xori, MATCH_XORI, MASK_XORI)
DECLARE_INSN(srli, MATCH_SRLI, MASK_SRLI)
DECLARE_INSN(srai, MATCH_SRAI, MASK_SRAI)
DECLARE_INSN(ori, MATCH_ORI, MASK_ORI)
DECLARE_INSN(andi, MATCH_ANDI, MASK_ANDI)
DECLARE_INSN(add, MATCH_ADD, MASK_ADD)
DECLARE_INSN(sub, MATCH_SUB, MASK_SUB)
DECLARE_INSN(sll, MATCH_SLL, MASK_SLL)
DECLARE_INSN(slt, MATCH_SLT, MASK_SLT)
DECLARE_INSN(sltu, MATCH_SLTU, MASK_SLTU)
DECLARE_INSN(xor, MATCH_XOR, MASK_XOR)
DECLARE_INSN(srl, MATCH_SRL, MASK_SRL)
DECLARE_INSN(sra, MATCH_SRA, MASK_SRA)
DECLARE_INSN(or, MATCH_OR, MASK_OR)
DECLARE_INSN(and, MATCH_AND, MASK_AND)
DECLARE_INSN(addiw, MATCH_ADDIW, MASK_ADDIW)
DECLARE_INSN(slliw, MATCH_SLLIW, MASK_SLLIW)
DECLARE_INSN(srliw, MATCH_SRLIW, MASK_SRLIW)
DECLARE_INSN(sraiw, MATCH_SRAIW, MASK_SRAIW)
DECLARE_INSN(addw, MATCH_ADDW, MASK_ADDW)
DECLARE_INSN(subw, MATCH_SUBW, MASK_SUBW)
DECLARE_INSN(sllw, MATCH_SLLW, MASK_SLLW)
DECLARE_INSN(srlw, MATCH_SRLW, MASK_SRLW)
DECLARE_INSN(sraw, MATCH_SRAW, MASK_SRAW)
DECLARE_INSN(lb, MATCH_LB, MASK_LB)
DECLARE_INSN(lh, MATCH_LH, MASK_LH)
DECLARE_INSN(lw, MATCH_LW, MASK_LW)
DECLARE_INSN(ld, MATCH_LD, MASK_LD)
DECLARE_INSN(lbu, MATCH_LBU, MASK_LBU)
DECLARE_INSN(lhu, MATCH_LHU, MASK_LHU)
DECLARE_INSN(lwu, MATCH_LWU, MASK_LWU)
DECLARE_INSN(sb, MATCH_SB, MASK_SB)
DECLARE_INSN(sh, MATCH_SH, MASK_SH)
DECLARE_INSN(sw, MATCH_SW, MASK_SW)
DECLARE_INSN(sd, MATCH_SD, MASK_SD)
DECLARE_INSN(pause, MATCH_PAUSE, MASK_PAUSE)
DECLARE_INSN(fence, MATCH_FENCE, MASK_FENCE)
DECLARE_INSN(fence_i, MATCH_FENCE_I, MASK_FENCE_I)
DECLARE_INSN(mul, MATCH_MUL, MASK_MUL)
DECLARE_INSN(mulh, MATCH_MULH, MASK_MULH)
DECLARE_INSN(mulhsu, MATCH_MULHSU, MASK_MULHSU)
DECLARE_INSN(mulhu, MATCH_MULHU, MASK_MULHU)
DECLARE_INSN(div, MATCH_DIV, MASK_DIV)
DECLARE_INSN(divu, MATCH_DIVU, MASK_DIVU)
DECLARE_INSN(rem, MATCH_REM, MASK_REM)
DECLARE_INSN(remu, MATCH_REMU, MASK_REMU)
DECLARE_INSN(mulw, MATCH_MULW, MASK_MULW)
DECLARE_INSN(divw, MATCH_DIVW, MASK_DIVW)
DECLARE_INSN(divuw, MATCH_DIVUW, MASK_DIVUW)
DECLARE_INSN(remw, MATCH_REMW, MASK_REMW)
DECLARE_INSN(remuw, MATCH_REMUW, MASK_REMUW)
DECLARE_INSN(amoadd_w, MATCH_AMOADD_W, MASK_AMOADD_W)
DECLARE_INSN(amoxor_w, MATCH_AMOXOR_W, MASK_AMOXOR_W)
DECLARE_INSN(amoor_w, MATCH_AMOOR_W, MASK_AMOOR_W)
DECLARE_INSN(amoand_w, MATCH_AMOAND_W, MASK_AMOAND_W)
DECLARE_INSN(amomin_w, MATCH_AMOMIN_W, MASK_AMOMIN_W)
DECLARE_INSN(amomax_w, MATCH_AMOMAX_W, MASK_AMOMAX_W)
DECLARE_INSN(amominu_w, MATCH_AMOMINU_W, MASK_AMOMINU_W)
DECLARE_INSN(amomaxu_w, MATCH_AMOMAXU_W, MASK_AMOMAXU_W)
DECLARE_INSN(amoswap_w, MATCH_AMOSWAP_W, MASK_AMOSWAP_W)
DECLARE_INSN(lr_w, MATCH_LR_W, MASK_LR_W)
DECLARE_INSN(sc_w, MATCH_SC_W, MASK_SC_W)
DECLARE_INSN(amoadd_d, MATCH_AMOADD_D, MASK_AMOADD_D)
DECLARE_INSN(amoxor_d, MATCH_AMOXOR_D, MASK_AMOXOR_D)
DECLARE_INSN(amoor_d, MATCH_AMOOR_D, MASK_AMOOR_D)
DECLARE_INSN(amoand_d, MATCH_AMOAND_D, MASK_AMOAND_D)
DECLARE_INSN(amomin_d, MATCH_AMOMIN_D, MASK_AMOMIN_D)
DECLARE_INSN(amomax_d, MATCH_AMOMAX_D, MASK_AMOMAX_D)
DECLARE_INSN(amominu_d, MATCH_AMOMINU_D, MASK_AMOMINU_D)
DECLARE_INSN(amomaxu_d, MATCH_AMOMAXU_D, MASK_AMOMAXU_D)
DECLARE_INSN(amoswap_d, MATCH_AMOSWAP_D, MASK_AMOSWAP_D)
DECLARE_INSN(lr_d, MATCH_LR_D, MASK_LR_D)
DECLARE_INSN(sc_d, MATCH_SC_D, MASK_SC_D)
DECLARE_INSN(ecall, MATCH_ECALL, MASK_ECALL)
DECLARE_INSN(ebreak, MATCH_EBREAK, MASK_EBREAK)
DECLARE_INSN(uret, MATCH_URET, MASK_URET)
DECLARE_INSN(sret, MATCH_SRET, MASK_SRET)
DECLARE_INSN(hret, MATCH_HRET, MASK_HRET)
DECLARE_INSN(mret, MATCH_MRET, MASK_MRET)
DECLARE_INSN(dret, MATCH_DRET, MASK_DRET)
DECLARE_INSN(sfence_vm, MATCH_SFENCE_VM, MASK_SFENCE_VM)
DECLARE_INSN(sfence_vma, MATCH_SFENCE_VMA, MASK_SFENCE_VMA)
DECLARE_INSN(wfi, MATCH_WFI, MASK_WFI)
DECLARE_INSN(csrrw, MATCH_CSRRW, MASK_CSRRW)
DECLARE_INSN(csrrs, MATCH_CSRRS, MASK_CSRRS)
DECLARE_INSN(csrrc, MATCH_CSRRC, MASK_CSRRC)
DECLARE_INSN(csrrwi, MATCH_CSRRWI, MASK_CSRRWI)
DECLARE_INSN(csrrsi, MATCH_CSRRSI, MASK_CSRRSI)
DECLARE_INSN(csrrci, MATCH_CSRRCI, MASK_CSRRCI)
DECLARE_INSN(fadd_s, MATCH_FADD_S, MASK_FADD_S)
DECLARE_INSN(fsub_s, MATCH_FSUB_S, MASK_FSUB_S)
DECLARE_INSN(fmul_s, MATCH_FMUL_S, MASK_FMUL_S)
DECLARE_INSN(fdiv_s, MATCH_FDIV_S, MASK_FDIV_S)
DECLARE_INSN(fsgnj_s, MATCH_FSGNJ_S, MASK_FSGNJ_S)
DECLARE_INSN(fsgnjn_s, MATCH_FSGNJN_S, MASK_FSGNJN_S)
DECLARE_INSN(fsgnjx_s, MATCH_FSGNJX_S, MASK_FSGNJX_S)
DECLARE_INSN(fmin_s, MATCH_FMIN_S, MASK_FMIN_S)
DECLARE_INSN(fmax_s, MATCH_FMAX_S, MASK_FMAX_S)
DECLARE_INSN(fsqrt_s, MATCH_FSQRT_S, MASK_FSQRT_S)
DECLARE_INSN(fadd_d, MATCH_FADD_D, MASK_FADD_D)
DECLARE_INSN(fsub_d, MATCH_FSUB_D, MASK_FSUB_D)
DECLARE_INSN(fmul_d, MATCH_FMUL_D, MASK_FMUL_D)
DECLARE_INSN(fdiv_d, MATCH_FDIV_D, MASK_FDIV_D)
DECLARE_INSN(fsgnj_d, MATCH_FSGNJ_D, MASK_FSGNJ_D)
DECLARE_INSN(fsgnjn_d, MATCH_FSGNJN_D, MASK_FSGNJN_D)
DECLARE_INSN(fsgnjx_d, MATCH_FSGNJX_D, MASK_FSGNJX_D)
DECLARE_INSN(fmin_d, MATCH_FMIN_D, MASK_FMIN_D)
DECLARE_INSN(fmax_d, MATCH_FMAX_D, MASK_FMAX_D)
DECLARE_INSN(fcvt_s_d, MATCH_FCVT_S_D, MASK_FCVT_S_D)
DECLARE_INSN(fcvt_d_s, MATCH_FCVT_D_S, MASK_FCVT_D_S)
DECLARE_INSN(fsqrt_d, MATCH_FSQRT_D, MASK_FSQRT_D)
DECLARE_INSN(fadd_q, MATCH_FADD_Q, MASK_FADD_Q)
DECLARE_INSN(fsub_q, MATCH_FSUB_Q, MASK_FSUB_Q)
DECLARE_INSN(fmul_q, MATCH_FMUL_Q, MASK_FMUL_Q)
DECLARE_INSN(fdiv_q, MATCH_FDIV_Q, MASK_FDIV_Q)
DECLARE_INSN(fsgnj_q, MATCH_FSGNJ_Q, MASK_FSGNJ_Q)
DECLARE_INSN(fsgnjn_q, MATCH_FSGNJN_Q, MASK_FSGNJN_Q)
DECLARE_INSN(fsgnjx_q, MATCH_FSGNJX_Q, MASK_FSGNJX_Q)
DECLARE_INSN(fmin_q, MATCH_FMIN_Q, MASK_FMIN_Q)
DECLARE_INSN(fmax_q, MATCH_FMAX_Q, MASK_FMAX_Q)
DECLARE_INSN(fcvt_s_q, MATCH_FCVT_S_Q, MASK_FCVT_S_Q)
DECLARE_INSN(fcvt_q_s, MATCH_FCVT_Q_S, MASK_FCVT_Q_S)
DECLARE_INSN(fcvt_d_q, MATCH_FCVT_D_Q, MASK_FCVT_D_Q)
DECLARE_INSN(fcvt_q_d, MATCH_FCVT_Q_D, MASK_FCVT_Q_D)
DECLARE_INSN(fsqrt_q, MATCH_FSQRT_Q, MASK_FSQRT_Q)
DECLARE_INSN(fle_s, MATCH_FLE_S, MASK_FLE_S)
DECLARE_INSN(flt_s, MATCH_FLT_S, MASK_FLT_S)
DECLARE_INSN(feq_s, MATCH_FEQ_S, MASK_FEQ_S)
DECLARE_INSN(fle_d, MATCH_FLE_D, MASK_FLE_D)
DECLARE_INSN(flt_d, MATCH_FLT_D, MASK_FLT_D)
DECLARE_INSN(feq_d, MATCH_FEQ_D, MASK_FEQ_D)
DECLARE_INSN(fle_q, MATCH_FLE_Q, MASK_FLE_Q)
DECLARE_INSN(flt_q, MATCH_FLT_Q, MASK_FLT_Q)
DECLARE_INSN(feq_q, MATCH_FEQ_Q, MASK_FEQ_Q)
DECLARE_INSN(fcvt_w_s, MATCH_FCVT_W_S, MASK_FCVT_W_S)
DECLARE_INSN(fcvt_wu_s, MATCH_FCVT_WU_S, MASK_FCVT_WU_S)
DECLARE_INSN(fcvt_l_s, MATCH_FCVT_L_S, MASK_FCVT_L_S)
DECLARE_INSN(fcvt_lu_s, MATCH_FCVT_LU_S, MASK_FCVT_LU_S)
DECLARE_INSN(fmv_x_s, MATCH_FMV_X_S, MASK_FMV_X_S)
DECLARE_INSN(fclass_s, MATCH_FCLASS_S, MASK_FCLASS_S)
DECLARE_INSN(fcvt_w_d, MATCH_FCVT_W_D, MASK_FCVT_W_D)
DECLARE_INSN(fcvt_wu_d, MATCH_FCVT_WU_D, MASK_FCVT_WU_D)
DECLARE_INSN(fcvt_l_d, MATCH_FCVT_L_D, MASK_FCVT_L_D)
DECLARE_INSN(fcvt_lu_d, MATCH_FCVT_LU_D, MASK_FCVT_LU_D)
DECLARE_INSN(fmv_x_d, MATCH_FMV_X_D, MASK_FMV_X_D)
DECLARE_INSN(fclass_d, MATCH_FCLASS_D, MASK_FCLASS_D)
DECLARE_INSN(fcvt_w_q, MATCH_FCVT_W_Q, MASK_FCVT_W_Q)
DECLARE_INSN(fcvt_wu_q, MATCH_FCVT_WU_Q, MASK_FCVT_WU_Q)
DECLARE_INSN(fcvt_l_q, MATCH_FCVT_L_Q, MASK_FCVT_L_Q)
DECLARE_INSN(fcvt_lu_q, MATCH_FCVT_LU_Q, MASK_FCVT_LU_Q)
DECLARE_INSN(fmv_x_q, MATCH_FMV_X_Q, MASK_FMV_X_Q)
DECLARE_INSN(fclass_q, MATCH_FCLASS_Q, MASK_FCLASS_Q)
DECLARE_INSN(fcvt_s_w, MATCH_FCVT_S_W, MASK_FCVT_S_W)
DECLARE_INSN(fcvt_s_wu, MATCH_FCVT_S_WU, MASK_FCVT_S_WU)
DECLARE_INSN(fcvt_s_l, MATCH_FCVT_S_L, MASK_FCVT_S_L)
DECLARE_INSN(fcvt_s_lu, MATCH_FCVT_S_LU, MASK_FCVT_S_LU)
DECLARE_INSN(fmv_s_x, MATCH_FMV_S_X, MASK_FMV_S_X)
DECLARE_INSN(fcvt_d_w, MATCH_FCVT_D_W, MASK_FCVT_D_W)
DECLARE_INSN(fcvt_d_wu, MATCH_FCVT_D_WU, MASK_FCVT_D_WU)
DECLARE_INSN(fcvt_d_l, MATCH_FCVT_D_L, MASK_FCVT_D_L)
DECLARE_INSN(fcvt_d_lu, MATCH_FCVT_D_LU, MASK_FCVT_D_LU)
DECLARE_INSN(fmv_d_x, MATCH_FMV_D_X, MASK_FMV_D_X)
DECLARE_INSN(fcvt_q_w, MATCH_FCVT_Q_W, MASK_FCVT_Q_W)
DECLARE_INSN(fcvt_q_wu, MATCH_FCVT_Q_WU, MASK_FCVT_Q_WU)
DECLARE_INSN(fcvt_q_l, MATCH_FCVT_Q_L, MASK_FCVT_Q_L)
DECLARE_INSN(fcvt_q_lu, MATCH_FCVT_Q_LU, MASK_FCVT_Q_LU)
DECLARE_INSN(fmv_q_x, MATCH_FMV_Q_X, MASK_FMV_Q_X)
DECLARE_INSN(clz, MATCH_CLZ, MASK_CLZ)
DECLARE_INSN(ctz, MATCH_CTZ, MASK_CTZ)
DECLARE_INSN(cpop, MATCH_CPOP, MASK_CPOP)
DECLARE_INSN(min, MATCH_MIN, MASK_MIN)
DECLARE_INSN(minu, MATCH_MINU, MASK_MINU)
DECLARE_INSN(max, MATCH_MAX, MASK_MAX)
DECLARE_INSN(maxu, MATCH_MAXU, MASK_MAXU)
DECLARE_INSN(sext_b, MATCH_SEXT_B, MASK_SEXT_B)
DECLARE_INSN(sext_h, MATCH_SEXT_H, MASK_SEXT_H)
DECLARE_INSN(andn, MATCH_ANDN, MASK_ANDN)
DECLARE_INSN(orn, MATCH_ORN, MASK_ORN)
DECLARE_INSN(xnor, MATCH_XNOR, MASK_XNOR)
DECLARE_INSN(rol, MATCH_ROL, MASK_ROL)
DECLARE_INSN(ror, MATCH_ROR, MASK_ROR)
DECLARE_INSN(rori, MATCH_RORI, MASK_RORI)
DECLARE_INSN(clzw, MATCH_CLZW, MASK_CLZW)
DECLARE_INSN(ctzw, MATCH_CTZW, MASK_CTZW)
DECLARE_INSN(cpopw, MATCH_CPOPW, MASK_CPOPW)
DECLARE_INSN(rolw, MATCH_ROLW, MASK_ROLW)
DECLARE_INSN(rorw, MATCH_RORW, MASK_RORW)
DECLARE_INSN(roriw, MATCH_RORIW, MASK_RORIW)
DECLARE_INSN(sh1add, MATCH_SH1ADD, MASK_SH1ADD)
DECLARE_INSN(sh2add, MATCH_SH2ADD, MASK_SH2ADD)
DECLARE_INSN(sh3add, MATCH_SH3ADD, MASK_SH3ADD)
DECLARE_INSN(sh1add_uw, MATCH_SH1ADD_UW, MASK_SH1ADD_UW)
DECLARE_INSN(sh2add_uw, MATCH_SH2ADD_UW, MASK_SH2ADD_UW)
DECLARE_INSN(sh3add_uw, MATCH_SH3ADD_UW, MASK_SH3ADD_UW)
DECLARE_INSN(add_uw, MATCH_ADD_UW, MASK_ADD_UW)
DECLARE_INSN(slli_uw, MATCH_SLLI_UW, MASK_SLLI_UW)
DECLARE_INSN(clmul, MATCH_CLMUL, MASK_CLMUL)
DECLARE_INSN(clmulh, MATCH_CLMULH, MASK_CLMULH)
DECLARE_INSN(clmulr, MATCH_CLMULR, MASK_CLMULR)
DECLARE_INSN(bclr, MATCH_BCLR, MASK_BCLR)
DECLARE_INSN(bclri, MATCH_BCLRI, MASK_BCLRI)
DECLARE_INSN(bext, MATCH_BEXT, MASK_BEXT)
DECLARE_INSN(bexti, MATCH_BEXTI, MASK_BEXTI)
DECLARE_INSN(binv, MATCH_BINV, MASK_BINV)
DECLARE_INSN(binvi, MATCH_BINVI, MASK_BINVI)
DECLARE_INSN(bset, MATCH_BSET, MASK_BSET)
DECLARE_INSN(bseti, MATCH_BSETI, MASK_BSETI)
DECLARE_INSN(flw, MATCH_FLW, MASK_FLW)
DECLARE_INSN(fld, MATCH_FLD, MASK_FLD)
DECLARE_INSN(flq, MATCH_FLQ, MASK_FLQ)
DECLARE_INSN(fsw, MATCH_FSW, MASK_FSW)
DECLARE_INSN(fsd, MATCH_FSD, MASK_FSD)
DECLARE_INSN(fsq, MATCH_FSQ, MASK_FSQ)
DECLARE_INSN(fmadd_s, MATCH_FMADD_S, MASK_FMADD_S)
DECLARE_INSN(fmsub_s, MATCH_FMSUB_S, MASK_FMSUB_S)
DECLARE_INSN(fnmsub_s, MATCH_FNMSUB_S, MASK_FNMSUB_S)
DECLARE_INSN(fnmadd_s, MATCH_FNMADD_S, MASK_FNMADD_S)
DECLARE_INSN(fmadd_d, MATCH_FMADD_D, MASK_FMADD_D)
DECLARE_INSN(fmsub_d, MATCH_FMSUB_D, MASK_FMSUB_D)
DECLARE_INSN(fnmsub_d, MATCH_FNMSUB_D, MASK_FNMSUB_D)
DECLARE_INSN(fnmadd_d, MATCH_FNMADD_D, MASK_FNMADD_D)
DECLARE_INSN(fmadd_q, MATCH_FMADD_Q, MASK_FMADD_Q)
DECLARE_INSN(fmsub_q, MATCH_FMSUB_Q, MASK_FMSUB_Q)
DECLARE_INSN(fnmsub_q, MATCH_FNMSUB_Q, MASK_FNMSUB_Q)
DECLARE_INSN(fnmadd_q, MATCH_FNMADD_Q, MASK_FNMADD_Q)
DECLARE_INSN(c_addi4spn, MATCH_C_ADDI4SPN, MASK_C_ADDI4SPN)
DECLARE_INSN(c_fld, MATCH_C_FLD, MASK_C_FLD)
DECLARE_INSN(c_lw, MATCH_C_LW, MASK_C_LW)
DECLARE_INSN(c_flw, MATCH_C_FLW, MASK_C_FLW)
DECLARE_INSN(c_fsd, MATCH_C_FSD, MASK_C_FSD)
DECLARE_INSN(c_sw, MATCH_C_SW, MASK_C_SW)
DECLARE_INSN(c_fsw, MATCH_C_FSW, MASK_C_FSW)
DECLARE_INSN(c_addi, MATCH_C_ADDI, MASK_C_ADDI)
DECLARE_INSN(c_jal, MATCH_C_JAL, MASK_C_JAL)
DECLARE_INSN(c_li, MATCH_C_LI, MASK_C_LI)
DECLARE_INSN(c_lui, MATCH_C_LUI, MASK_C_LUI)
DECLARE_INSN(c_srli, MATCH_C_SRLI, MASK_C_SRLI)
DECLARE_INSN(c_srai, MATCH_C_SRAI, MASK_C_SRAI)
DECLARE_INSN(c_andi, MATCH_C_ANDI, MASK_C_ANDI)
DECLARE_INSN(c_sub, MATCH_C_SUB, MASK_C_SUB)
DECLARE_INSN(c_xor, MATCH_C_XOR, MASK_C_XOR)
DECLARE_INSN(c_or, MATCH_C_OR, MASK_C_OR)
DECLARE_INSN(c_and, MATCH_C_AND, MASK_C_AND)
DECLARE_INSN(c_subw, MATCH_C_SUBW, MASK_C_SUBW)
DECLARE_INSN(c_addw, MATCH_C_ADDW, MASK_C_ADDW)
DECLARE_INSN(c_j, MATCH_C_J, MASK_C_J)
DECLARE_INSN(c_beqz, MATCH_C_BEQZ, MASK_C_BEQZ)
DECLARE_INSN(c_bnez, MATCH_C_BNEZ, MASK_C_BNEZ)
DECLARE_INSN(c_slli, MATCH_C_SLLI, MASK_C_SLLI)
DECLARE_INSN(c_fldsp, MATCH_C_FLDSP, MASK_C_FLDSP)
DECLARE_INSN(c_lwsp, MATCH_C_LWSP, MASK_C_LWSP)
DECLARE_INSN(c_flwsp, MATCH_C_FLWSP, MASK_C_FLWSP)
DECLARE_INSN(c_mv, MATCH_C_MV, MASK_C_MV)
DECLARE_INSN(c_add, MATCH_C_ADD, MASK_C_ADD)
DECLARE_INSN(c_fsdsp, MATCH_C_FSDSP, MASK_C_FSDSP)
DECLARE_INSN(c_swsp, MATCH_C_SWSP, MASK_C_SWSP)
DECLARE_INSN(c_fswsp, MATCH_C_FSWSP, MASK_C_FSWSP)
DECLARE_INSN(c_nop, MATCH_C_NOP, MASK_C_NOP)
DECLARE_INSN(c_addi16sp, MATCH_C_ADDI16SP, MASK_C_ADDI16SP)
DECLARE_INSN(c_jr, MATCH_C_JR, MASK_C_JR)
DECLARE_INSN(c_jalr, MATCH_C_JALR, MASK_C_JALR)
DECLARE_INSN(c_ebreak, MATCH_C_EBREAK, MASK_C_EBREAK)
DECLARE_INSN(c_ld, MATCH_C_LD, MASK_C_LD)
DECLARE_INSN(c_sd, MATCH_C_SD, MASK_C_SD)
DECLARE_INSN(c_addiw, MATCH_C_ADDIW, MASK_C_ADDIW)
DECLARE_INSN(c_ldsp, MATCH_C_LDSP, MASK_C_LDSP)
DECLARE_INSN(c_sdsp, MATCH_C_SDSP, MASK_C_SDSP)
DECLARE_INSN(sinval_vma, MATCH_SINVAL_VMA, MASK_SINVAL_VMA)
DECLARE_INSN(sfence_w_inval, MATCH_SFENCE_W_INVAL, MASK_SFENCE_W_INVAL)
DECLARE_INSN(sfence_inval_ir, MATCH_SFENCE_INVAL_IR, MASK_SFENCE_INVAL_IR)
DECLARE_INSN(hinval_vvma, MATCH_HINVAL_VVMA, MASK_HINVAL_VVMA)
DECLARE_INSN(hinval_gvma, MATCH_HINVAL_GVMA, MASK_HINVAL_GVMA)
DECLARE_INSN(custom0, MATCH_CUSTOM0, MASK_CUSTOM0)
DECLARE_INSN(custom0_rs1, MATCH_CUSTOM0_RS1, MASK_CUSTOM0_RS1)
DECLARE_INSN(custom0_rs1_rs2, MATCH_CUSTOM0_RS1_RS2, MASK_CUSTOM0_RS1_RS2)
DECLARE_INSN(custom0_rd, MATCH_CUSTOM0_RD, MASK_CUSTOM0_RD)
DECLARE_INSN(custom0_rd_rs1, MATCH_CUSTOM0_RD_RS1, MASK_CUSTOM0_RD_RS1)
DECLARE_INSN(custom0_rd_rs1_rs2, MATCH_CUSTOM0_RD_RS1_RS2, MASK_CUSTOM0_RD_RS1_RS2)
DECLARE_INSN(custom1, MATCH_CUSTOM1, MASK_CUSTOM1)
DECLARE_INSN(custom1_rs1, MATCH_CUSTOM1_RS1, MASK_CUSTOM1_RS1)
DECLARE_INSN(custom1_rs1_rs2, MATCH_CUSTOM1_RS1_RS2, MASK_CUSTOM1_RS1_RS2)
DECLARE_INSN(custom1_rd, MATCH_CUSTOM1_RD, MASK_CUSTOM1_RD)
DECLARE_INSN(custom1_rd_rs1, MATCH_CUSTOM1_RD_RS1, MASK_CUSTOM1_RD_RS1)
DECLARE_INSN(custom1_rd_rs1_rs2, MATCH_CUSTOM1_RD_RS1_RS2, MASK_CUSTOM1_RD_RS1_RS2)
DECLARE_INSN(custom2, MATCH_CUSTOM2, MASK_CUSTOM2)
DECLARE_INSN(custom2_rs1, MATCH_CUSTOM2_RS1, MASK_CUSTOM2_RS1)
DECLARE_INSN(custom2_rs1_rs2, MATCH_CUSTOM2_RS1_RS2, MASK_CUSTOM2_RS1_RS2)
DECLARE_INSN(custom2_rd, MATCH_CUSTOM2_RD, MASK_CUSTOM2_RD)
DECLARE_INSN(custom2_rd_rs1, MATCH_CUSTOM2_RD_RS1, MASK_CUSTOM2_RD_RS1)
DECLARE_INSN(custom2_rd_rs1_rs2, MATCH_CUSTOM2_RD_RS1_RS2, MASK_CUSTOM2_RD_RS1_RS2)
DECLARE_INSN(custom3, MATCH_CUSTOM3, MASK_CUSTOM3)
DECLARE_INSN(custom3_rs1, MATCH_CUSTOM3_RS1, MASK_CUSTOM3_RS1)
DECLARE_INSN(custom3_rs1_rs2, MATCH_CUSTOM3_RS1_RS2, MASK_CUSTOM3_RS1_RS2)
DECLARE_INSN(custom3_rd, MATCH_CUSTOM3_RD, MASK_CUSTOM3_RD)
DECLARE_INSN(custom3_rd_rs1, MATCH_CUSTOM3_RD_RS1, MASK_CUSTOM3_RD_RS1)
DECLARE_INSN(custom3_rd_rs1_rs2, MATCH_CUSTOM3_RD_RS1_RS2, MASK_CUSTOM3_RD_RS1_RS2)
#endif /* DECLARE_INSN.  */
#ifdef DECLARE_CSR
/* Privileged.  */
DECLARE_CSR(ustatus, CSR_USTATUS, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(uie, CSR_UIE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(utvec, CSR_UTVEC, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(uscratch, CSR_USCRATCH, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(uepc, CSR_UEPC, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(ucause, CSR_UCAUSE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(utval, CSR_UTVAL, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(uip, CSR_UIP, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(cycle, CSR_CYCLE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(time, CSR_TIME, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(instret, CSR_INSTRET, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter3, CSR_HPMCOUNTER3, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter4, CSR_HPMCOUNTER4, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter5, CSR_HPMCOUNTER5, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter6, CSR_HPMCOUNTER6, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter7, CSR_HPMCOUNTER7, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter8, CSR_HPMCOUNTER8, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter9, CSR_HPMCOUNTER9, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter10, CSR_HPMCOUNTER10, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter11, CSR_HPMCOUNTER11, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter12, CSR_HPMCOUNTER12, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter13, CSR_HPMCOUNTER13, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter14, CSR_HPMCOUNTER14, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter15, CSR_HPMCOUNTER15, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter16, CSR_HPMCOUNTER16, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter17, CSR_HPMCOUNTER17, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter18, CSR_HPMCOUNTER18, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter19, CSR_HPMCOUNTER19, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter20, CSR_HPMCOUNTER20, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter21, CSR_HPMCOUNTER21, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter22, CSR_HPMCOUNTER22, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter23, CSR_HPMCOUNTER23, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter24, CSR_HPMCOUNTER24, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter25, CSR_HPMCOUNTER25, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter26, CSR_HPMCOUNTER26, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter27, CSR_HPMCOUNTER27, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter28, CSR_HPMCOUNTER28, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter29, CSR_HPMCOUNTER29, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter30, CSR_HPMCOUNTER30, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter31, CSR_HPMCOUNTER31, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(cycleh, CSR_CYCLEH, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(timeh, CSR_TIMEH, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(instreth, CSR_INSTRETH, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter3h, CSR_HPMCOUNTER3H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter4h, CSR_HPMCOUNTER4H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter5h, CSR_HPMCOUNTER5H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter6h, CSR_HPMCOUNTER6H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter7h, CSR_HPMCOUNTER7H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter8h, CSR_HPMCOUNTER8H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter9h, CSR_HPMCOUNTER9H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter10h, CSR_HPMCOUNTER10H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter11h, CSR_HPMCOUNTER11H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter12h, CSR_HPMCOUNTER12H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter13h, CSR_HPMCOUNTER13H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter14h, CSR_HPMCOUNTER14H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter15h, CSR_HPMCOUNTER15H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter16h, CSR_HPMCOUNTER16H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter17h, CSR_HPMCOUNTER17H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter18h, CSR_HPMCOUNTER18H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter19h, CSR_HPMCOUNTER19H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter20h, CSR_HPMCOUNTER20H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter21h, CSR_HPMCOUNTER21H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter22h, CSR_HPMCOUNTER22H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter23h, CSR_HPMCOUNTER23H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter24h, CSR_HPMCOUNTER24H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter25h, CSR_HPMCOUNTER25H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter26h, CSR_HPMCOUNTER26H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter27h, CSR_HPMCOUNTER27H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter28h, CSR_HPMCOUNTER28H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter29h, CSR_HPMCOUNTER29H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter30h, CSR_HPMCOUNTER30H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(hpmcounter31h, CSR_HPMCOUNTER31H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(sstatus, CSR_SSTATUS, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(sedeleg, CSR_SEDELEG, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(sideleg, CSR_SIDELEG, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(sie, CSR_SIE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(stvec, CSR_STVEC, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(scounteren, CSR_SCOUNTEREN, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(sscratch, CSR_SSCRATCH, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(sepc, CSR_SEPC, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(scause, CSR_SCAUSE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(stval, CSR_STVAL, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(sip, CSR_SIP, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(stimecmp, CSR_STIMECMP, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(stimecmph, CSR_STIMECMPH, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(satp, CSR_SATP, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(senvcfg, CSR_SENVCFG, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mvendorid, CSR_MVENDORID, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(marchid, CSR_MARCHID, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mimpid, CSR_MIMPID, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhartid, CSR_MHARTID, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mconfigptr, CSR_MCONFIGPTR, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mstatus, CSR_MSTATUS, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(misa, CSR_MISA, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(medeleg, CSR_MEDELEG, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mideleg, CSR_MIDELEG, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mie, CSR_MIE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mtvec, CSR_MTVEC, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mcounteren, CSR_MCOUNTEREN, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mtvt, CSR_MTVT, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(menvcfg, CSR_MENVCFG, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mscratch, CSR_MSCRATCH, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mepc, CSR_MEPC, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mcause, CSR_MCAUSE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mtval, CSR_MTVAL, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mip, CSR_MIP, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mnxti, CSR_MNXTI, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mintstatus, CSR_MINTSTATUS, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mscratchcsw, CSR_MSCRATCHCSW, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mscratchcswl, CSR_MSCRATCHCSWL, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mclicbase, CSR_MCLICBASE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg0, CSR_PMPCFG0, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg1, CSR_PMPCFG1, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg2, CSR_PMPCFG2, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg3, CSR_PMPCFG3, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg4, CSR_PMPCFG4, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg5, CSR_PMPCFG5, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg6, CSR_PMPCFG6, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg7, CSR_PMPCFG7, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg8, CSR_PMPCFG8, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg9, CSR_PMPCFG9, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg10, CSR_PMPCFG10, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg11, CSR_PMPCFG11, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg12, CSR_PMPCFG12, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg13, CSR_PMPCFG13, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg14, CSR_PMPCFG14, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpcfg15, CSR_PMPCFG15, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr0, CSR_PMPADDR0, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr1, CSR_PMPADDR1, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr2, CSR_PMPADDR2, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr3, CSR_PMPADDR3, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr4, CSR_PMPADDR4, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr5, CSR_PMPADDR5, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr6, CSR_PMPADDR6, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr7, CSR_PMPADDR7, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr8, CSR_PMPADDR8, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr9, CSR_PMPADDR9, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr10, CSR_PMPADDR10, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr11, CSR_PMPADDR11, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr12, CSR_PMPADDR12, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr13, CSR_PMPADDR13, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr14, CSR_PMPADDR14, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr15, CSR_PMPADDR15, CSR_CLASS_I, PRIV_SPEC_CLASS_1P10, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr16, CSR_PMPADDR16, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr17, CSR_PMPADDR17, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr18, CSR_PMPADDR18, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr19, CSR_PMPADDR19, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr20, CSR_PMPADDR20, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr21, CSR_PMPADDR21, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr22, CSR_PMPADDR22, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr23, CSR_PMPADDR23, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr24, CSR_PMPADDR24, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr25, CSR_PMPADDR25, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr26, CSR_PMPADDR26, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr27, CSR_PMPADDR27, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr28, CSR_PMPADDR28, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr29, CSR_PMPADDR29, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr30, CSR_PMPADDR30, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr31, CSR_PMPADDR31, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr32, CSR_PMPADDR32, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr33, CSR_PMPADDR33, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr34, CSR_PMPADDR34, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr35, CSR_PMPADDR35, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr36, CSR_PMPADDR36, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr37, CSR_PMPADDR37, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr38, CSR_PMPADDR38, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr39, CSR_PMPADDR39, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr40, CSR_PMPADDR40, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr41, CSR_PMPADDR41, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr42, CSR_PMPADDR42, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr43, CSR_PMPADDR43, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr44, CSR_PMPADDR44, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr45, CSR_PMPADDR45, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr46, CSR_PMPADDR46, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr47, CSR_PMPADDR47, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr48, CSR_PMPADDR48, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr49, CSR_PMPADDR49, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr50, CSR_PMPADDR50, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr51, CSR_PMPADDR51, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr52, CSR_PMPADDR52, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr53, CSR_PMPADDR53, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr54, CSR_PMPADDR54, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr55, CSR_PMPADDR55, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr56, CSR_PMPADDR56, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr57, CSR_PMPADDR57, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr58, CSR_PMPADDR58, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr59, CSR_PMPADDR59, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr60, CSR_PMPADDR60, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr61, CSR_PMPADDR61, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr62, CSR_PMPADDR62, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(pmpaddr63, CSR_PMPADDR63, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mcycle, CSR_MCYCLE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(minstret, CSR_MINSTRET, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter3, CSR_MHPMCOUNTER3, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter4, CSR_MHPMCOUNTER4, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter5, CSR_MHPMCOUNTER5, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter6, CSR_MHPMCOUNTER6, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter7, CSR_MHPMCOUNTER7, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter8, CSR_MHPMCOUNTER8, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter9, CSR_MHPMCOUNTER9, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter10, CSR_MHPMCOUNTER10, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter11, CSR_MHPMCOUNTER11, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter12, CSR_MHPMCOUNTER12, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter13, CSR_MHPMCOUNTER13, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter14, CSR_MHPMCOUNTER14, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter15, CSR_MHPMCOUNTER15, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter16, CSR_MHPMCOUNTER16, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter17, CSR_MHPMCOUNTER17, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter18, CSR_MHPMCOUNTER18, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter19, CSR_MHPMCOUNTER19, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter20, CSR_MHPMCOUNTER20, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter21, CSR_MHPMCOUNTER21, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter22, CSR_MHPMCOUNTER22, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter23, CSR_MHPMCOUNTER23, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter24, CSR_MHPMCOUNTER24, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter25, CSR_MHPMCOUNTER25, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter26, CSR_MHPMCOUNTER26, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter27, CSR_MHPMCOUNTER27, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter28, CSR_MHPMCOUNTER28, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter29, CSR_MHPMCOUNTER29, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter30, CSR_MHPMCOUNTER30, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter31, CSR_MHPMCOUNTER31, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mcycleh, CSR_MCYCLEH, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(minstreth, CSR_MINSTRETH, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter3h, CSR_MHPMCOUNTER3H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter4h, CSR_MHPMCOUNTER4H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter5h, CSR_MHPMCOUNTER5H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter6h, CSR_MHPMCOUNTER6H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter7h, CSR_MHPMCOUNTER7H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter8h, CSR_MHPMCOUNTER8H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter9h, CSR_MHPMCOUNTER9H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter10h, CSR_MHPMCOUNTER10H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter11h, CSR_MHPMCOUNTER11H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter12h, CSR_MHPMCOUNTER12H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter13h, CSR_MHPMCOUNTER13H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter14h, CSR_MHPMCOUNTER14H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter15h, CSR_MHPMCOUNTER15H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter16h, CSR_MHPMCOUNTER16H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter17h, CSR_MHPMCOUNTER17H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter18h, CSR_MHPMCOUNTER18H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter19h, CSR_MHPMCOUNTER19H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter20h, CSR_MHPMCOUNTER20H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter21h, CSR_MHPMCOUNTER21H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter22h, CSR_MHPMCOUNTER22H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter23h, CSR_MHPMCOUNTER23H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter24h, CSR_MHPMCOUNTER24H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter25h, CSR_MHPMCOUNTER25H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter26h, CSR_MHPMCOUNTER26H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter27h, CSR_MHPMCOUNTER27H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter28h, CSR_MHPMCOUNTER28H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter29h, CSR_MHPMCOUNTER29H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter30h, CSR_MHPMCOUNTER30H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmcounter31h, CSR_MHPMCOUNTER31H, CSR_CLASS_I_32, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mcountinhibit, CSR_MCOUNTINHIBIT, CSR_CLASS_I, PRIV_SPEC_CLASS_1P11, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent3, CSR_MHPMEVENT3, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent4, CSR_MHPMEVENT4, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent5, CSR_MHPMEVENT5, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent6, CSR_MHPMEVENT6, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent7, CSR_MHPMEVENT7, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent8, CSR_MHPMEVENT8, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent9, CSR_MHPMEVENT9, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent10, CSR_MHPMEVENT10, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent11, CSR_MHPMEVENT11, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent12, CSR_MHPMEVENT12, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent13, CSR_MHPMEVENT13, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent14, CSR_MHPMEVENT14, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent15, CSR_MHPMEVENT15, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent16, CSR_MHPMEVENT16, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent17, CSR_MHPMEVENT17, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent18, CSR_MHPMEVENT18, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent19, CSR_MHPMEVENT19, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent20, CSR_MHPMEVENT20, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent21, CSR_MHPMEVENT21, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent22, CSR_MHPMEVENT22, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent23, CSR_MHPMEVENT23, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent24, CSR_MHPMEVENT24, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent25, CSR_MHPMEVENT25, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent26, CSR_MHPMEVENT26, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent27, CSR_MHPMEVENT27, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent28, CSR_MHPMEVENT28, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent29, CSR_MHPMEVENT29, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent30, CSR_MHPMEVENT30, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
DECLARE_CSR(mhpmevent31, CSR_MHPMEVENT31, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_DRAFT)
/* Dropped.  */
DECLARE_CSR(hstatus, CSR_HSTATUS, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(hedeleg, CSR_HEDELEG, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(hideleg, CSR_HIDELEG, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(hie, CSR_HIE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(htvec, CSR_HTVEC, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(hscratch, CSR_HSCRATCH, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(hepc, CSR_HEPC, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(hcause, CSR_HCAUSE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(hbadaddr, CSR_HBADADDR, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(hip, CSR_HIP, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(mbase, CSR_MBASE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(mbound, CSR_MBOUND, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(mibase, CSR_MIBASE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(mibound, CSR_MIBOUND, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(mdbase, CSR_MDBASE, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(mdbound, CSR_MDBOUND, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(mscounteren, CSR_MSCOUNTEREN, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR(mhcounteren, CSR_MHCOUNTEREN, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)

/* PRIV_SPEC_CLASS_1P12.  */
DECLARE_CSR(mseccfg, CSR_MSECCFG, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)

/* New CSR for C908.  */
DECLARE_CSR(mcper, CSR_MCPER, CSR_CLASS_I, PRIV_SPEC_CLASS_1P12, PRIV_SPEC_CLASS_DRAFT)

/* T-HEAD extentions.  */
DECLARE_CSR(mxstatus, CSR_MXSTATUS, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mhcr, CSR_MHCR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcor, CSR_MCOR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mccr2, CSR_MCCR2, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcer2, CSR_MCER2, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mhint, CSR_MHINT, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mrmr, CSR_MRMR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mrvbr, CSR_MRVBR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcer, CSR_MCER, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcounterwen, CSR_MCOUNTERWEN, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcounterinten, CSR_MCOUNTERINTEN, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcounterof, CSR_MCOUNTEROF, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mhint2, CSR_MHINT2, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mhint3, CSR_MHINT3, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mraddr, CSR_MRADDR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mexstatus, CSR_MEXSTATUS, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mnmicause, CSR_MNMICAUSE, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mnmipc, CSR_MNMIPC, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mhpmcr, CSR_MHPMCR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mhpmsr, CSR_MHPMSR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mhpmer, CSR_MHPMER, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(msmpr, CSR_MSMPR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mteecfg, CSR_MTEECFG, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mzoneid, CSR_MZONEID, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mllcpid, CSR_MLLCPID, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mllwp, CSR_MLLWP, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mdtcmcr, CSR_MDTCMCR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mitcmcr, CSR_MITCMCR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(miesr, CSR_MIESR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(msbepa, CSR_MSBEPA, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(msbepa2, CSR_MSBEPA2, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(usp, CSR_USP, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcins, CSR_MCINS, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcindex, CSR_MCINDEX, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcdata0, CSR_MCDATA0, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcdata1, CSR_MCDATA1, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(meicr, CSR_MEICR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(meicr2, CSR_MEICR2, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mbeaddr, CSR_MBEADDR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mebr, CSR_MEBR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mstatus, CSR_NT_MSTATUS, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mtvec, CSR_NT_MTVEC, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mie, CSR_NT_MIE, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mtvt, CSR_NT_MTVT, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mepc, CSR_NT_MEPC, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mcause, CSR_NT_MCAUSE, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mip, CSR_NT_MIP, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mintstate, CSR_NT_MINTSTATE, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mxstatus, CSR_NT_MXSTATUS, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_mebr, CSR_NT_MEBR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(nt_msp, CSR_NT_MSP, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(t_usp, CSR_T_USP, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(t_mdcr, CSR_T_MDCR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(t_mpcr, CSR_T_MPCR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(pmpteecfg, CSR_PMPTEECFG, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcpuid, CSR_MCPUID, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mapbaddr, CSR_MAPBADDR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mwmsr, CSR_MWMSR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mhaltcause, CSR_MHALTCAUSE, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mdbginfo, CSR_MDBGINFO, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mpcfifo, CSR_MPCFIFO, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(fxcr, CSR_FXCR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(smir, CSR_SMIR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(smel, CSR_SMEL, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(smeh, CSR_SMEH, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(smcir, CSR_SMCIR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(sxstatus, CSR_SXSTATUS, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shcr, CSR_SHCR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(scer2, CSR_SCER2, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(scer, CSR_SCER, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(scounterinten , CSR_SCOUNTERINTEN, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(scounterof, CSR_SCOUNTEROF, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shint, CSR_SHINT, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shint2, CSR_SHINT2, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpminhibit, CSR_SHPMINHIBIT, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcr, CSR_SHPMCR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmsr, CSR_SHPMSR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmer, CSR_SHPMER, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(sllcpid, CSR_SLLCPID, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(sllwp, CSR_SLLWP, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(siesr, CSR_SIESR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(sbeaddr, CSR_SBEADDR, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(ssbepa, CSR_SSBEPA, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(ssbepa2, CSR_SSBEPA2, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(scycle, CSR_SCYCLE, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter1, CSR_SHPMCOUNTER1, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter2, CSR_SHPMCOUNTER2, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter3, CSR_SHPMCOUNTER3, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter4, CSR_SHPMCOUNTER4, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter5, CSR_SHPMCOUNTER5, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter6, CSR_SHPMCOUNTER6, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter7, CSR_SHPMCOUNTER7, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter8, CSR_SHPMCOUNTER8, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter9, CSR_SHPMCOUNTER9, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter10, CSR_SHPMCOUNTER10, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter11, CSR_SHPMCOUNTER11, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter12, CSR_SHPMCOUNTER12, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter13, CSR_SHPMCOUNTER13, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter14, CSR_SHPMCOUNTER14, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter15, CSR_SHPMCOUNTER15, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter16, CSR_SHPMCOUNTER16, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter17, CSR_SHPMCOUNTER17, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter18, CSR_SHPMCOUNTER18, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter19, CSR_SHPMCOUNTER19, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter20, CSR_SHPMCOUNTER20, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter21, CSR_SHPMCOUNTER21, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter22, CSR_SHPMCOUNTER22, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter23, CSR_SHPMCOUNTER23, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter24, CSR_SHPMCOUNTER24, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter25, CSR_SHPMCOUNTER25, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter26, CSR_SHPMCOUNTER26, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter27, CSR_SHPMCOUNTER27, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter28, CSR_SHPMCOUNTER28, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter29, CSR_SHPMCOUNTER29, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter30, CSR_SHPMCOUNTER30, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(shpmcounter31, CSR_SHPMCOUNTER31, CSR_CLASS_THEAD, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)

/* Unprivileged.  */
DECLARE_CSR(fflags, CSR_FFLAGS, CSR_CLASS_F, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(frm, CSR_FRM, CSR_CLASS_F, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(fcsr, CSR_FCSR, CSR_CLASS_F, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(vstart, CSR_VSTART, CSR_CLASS_V, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(vxsat, CSR_VXSAT, CSR_CLASS_V, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(vxrm, CSR_VXRM, CSR_CLASS_V, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(vcsr, CSR_VCSR, CSR_CLASS_V, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(vl, CSR_VL, CSR_CLASS_V, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(vtype, CSR_VTYPE, CSR_CLASS_V, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(vlenb, CSR_VLENB, CSR_CLASS_V, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(dcsr, CSR_DCSR, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(dpc, CSR_DPC, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(dscratch0, CSR_DSCRATCH0, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(dscratch1, CSR_DSCRATCH1, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(tselect, CSR_TSELECT, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(tdata1, CSR_TDATA1, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(tdata2, CSR_TDATA2, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(tdata3, CSR_TDATA3, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(tinfo, CSR_TINFO, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(tcontrol, CSR_TCONTROL, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(mcontext, CSR_MCONTEXT, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR(scontext, CSR_SCONTEXT, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
#endif /* DECLARE_CSR.  */
#ifdef DECLARE_CSR_ALIAS
DECLARE_CSR_ALIAS(ubadaddr, CSR_UTVAL, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR_ALIAS(sbadaddr, CSR_STVAL, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR_ALIAS(sptbr, CSR_SATP, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR_ALIAS(mbadaddr, CSR_MTVAL, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR_ALIAS(mucounteren, CSR_MCOUNTINHIBIT, CSR_CLASS_I, PRIV_SPEC_CLASS_1P9P1, PRIV_SPEC_CLASS_1P10)
DECLARE_CSR_ALIAS(dscratch, CSR_DSCRATCH0, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR_ALIAS(mcontrol, CSR_TDATA1, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR_ALIAS(icount, CSR_TDATA1, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR_ALIAS(itrigger, CSR_TDATA1, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR_ALIAS(etrigger, CSR_TDATA1, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR_ALIAS(textra32, CSR_TDATA3, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
DECLARE_CSR_ALIAS(textra64, CSR_TDATA3, CSR_CLASS_DEBUG, PRIV_SPEC_CLASS_NONE, PRIV_SPEC_CLASS_NONE)
#endif /* DECLARE_CSR_ALIAS.  */
