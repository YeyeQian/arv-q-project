from m5.objects import *


class O3_RISCV_c910_Load(FUDesc):
    opList = [ OpDesc(opClass='MemRead',opLat=2),
               OpDesc(opClass='FloatMemRead',opLat=2) ]
    count = 1

class O3_RISCV_c910_Store(FUDesc):
    opList = [ OpDesc(opClass='MemWrite',opLat=2),
               OpDesc(opClass='FloatMemWrite',opLat=2) ]
    count = 1

class O3_RISCV_c910_IntDiv(FUDesc):
    opList = [ OpDesc(opClass='IntDiv', opLat=12, pipelined=False) ]

class O3_RISCV_c910_FP_ALU(FUDesc):
    opList = [ OpDesc(opClass='FloatAdd', opLat=3),
               OpDesc(opClass='FloatCmp', opLat=3),
               OpDesc(opClass='FloatCvt', opLat=4) ]
    count = 2

class O3_RISCV_c910_IntALU(FUDesc):
    opList = [ OpDesc(opClass='IntMult', opLat=3), OpDesc(opClass='IntAlu') ]
    count = 2

class O3_RISCV_c910_FP_MultDiv(FUDesc):
    opList = [ OpDesc(opClass='FloatMult', opLat=4),
               OpDesc(opClass='FloatMultAcc', opLat=5),
               OpDesc(opClass='FloatMisc', opLat=3),
               OpDesc(opClass='FloatDiv', opLat=10, pipelined=False),
               OpDesc(opClass='FloatSqrt', opLat=10, pipelined=False) ]
    count = 2


class O3_RISCV_c910_FUP(FUPool):
    FUList = [ IntALU(), IntMultDiv(), FP_ALU(), FP_MultDiv(), ReadPort(),
               SIMD_Unit(), PredALU(), WritePort(), RdWrPort(), IprPort() ]

class O3_RISCV_c910(RiscvO3CPU):
    LQEntries = 16
    SQEntries = 16
    LSQDepCheckShift = 0
    LFSTSize = 1024
    SSITSize = 1024
    decodeToFetchDelay = 1
    renameToFetchDelay = 1
    iewToFetchDelay = 1
    commitToFetchDelay = 1
    renameToDecodeDelay = 1
    iewToDecodeDelay = 1
    commitToDecodeDelay = 1
    iewToRenameDelay = 1
    commitToRenameDelay = 1
    commitToIEWDelay = 1
    fetchWidth = 3
    fetchBufferSize = 16
    fetchToDecodeDelay = 3
    decodeWidth = 3
    decodeToRenameDelay = 2
    renameWidth = 3
    renameToIEWDelay = 1
    issueToExecuteDelay = 1
    dispatchWidth = 6
    issueWidth = 8
    wbWidth = 8
    fuPool = O3_RISCV_c910_FUP()
    iewToCommitDelay = 1
    renameToROBDelay = 1
    commitWidth = 3
    squashWidth = 8
    trapLatency = 13
    backComSize = 5
    forwardComSize = 5
    numPhysIntRegs = 92
    numPhysFloatRegs = 64
    numPhysVecRegs = 48
    numIQEntries = 32
    numROBEntries = 40