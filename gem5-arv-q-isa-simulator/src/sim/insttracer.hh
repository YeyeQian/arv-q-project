/*
 * Copyright (c) 2014, 2017, 2020 ARM Limited
 * All rights reserved
 *
 * The license below extends only to copyright in the software and shall
 * not be construed as granting a license to any other intellectual
 * property including but not limited to intellectual property relating
 * to a hardware implementation of the functionality of the software
 * licensed hereunder.  You may use the software subject to the license
 * terms below provided that you ensure that this notice is replicated
 * unmodified and in its entirety in all distributions of the software,
 * modified or unmodified, in source code or in binary form.
 *
 * Copyright (c) 2001-2005 The Regents of The University of Michigan
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer;
 * redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution;
 * neither the name of the copyright holders nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef __INSTRECORD_HH__
#define __INSTRECORD_HH__

#include <memory>

#include "arch/generic/pcstate.hh"
#include "arch/vecregs.hh"
#include "base/types.hh"
#include "config/the_isa.hh"
#include "cpu/inst_seq.hh"
#include "cpu/static_inst.hh"
#include "sim/sim_object.hh"

#include "arch/riscv/regs/keccak.hh"


namespace gem5
{

class ThreadContext;

namespace Trace {

class InstRecord
{
  protected:
    Tick when;

    // The following fields are initialized by the constructor and
    // thus guaranteed to be valid.
    ThreadContext *thread;
    // need to make this ref-counted so it doesn't go away before we
    // dump the record
    StaticInstPtr staticInst;
    std::unique_ptr<PCStateBase> pc;
    StaticInstPtr macroStaticInst;

    // The remaining fields are only valid for particular instruction
    // types (e.g, addresses for memory ops) or when particular
    // options are enabled (e.g., tracing full register contents).
    // Each data field has an associated valid flag to indicate
    // whether the data field is valid.

    /*** @defgroup mem
     * @{
     * Memory request information in the instruction accessed memory.
     * @see mem_valid
     */
    Addr addr = 0; ///< The address that was accessed
    Addr size = 0; ///< The size of the memory request
    unsigned flags = 0; ///< The flags that were assigned to the request.

    /** @} */

    /** @defgroup data
     * If this instruction wrote any data values they're recorded here
     * WARNING: Instructions are quite loose with with what they write
     * since many instructions write multiple values (e.g. destintation
     * register, flags, status, ...) This only captures the last write.
     * @TODO fix this and record all destintations that an instruction writes
     * @see data_status
     */
    union
    {
        uint64_t as_int;
        double as_double;
        TheISA::VecRegContainer* as_vec;
        TheISA::VecPredRegContainer* as_pred;
        TheISA::KeccakRegContainer * as_keccak_reg;
    } data = {0};

    /** @defgroup fetch_seq
     * This records the serial number that the instruction was fetched in.
     * @see fetch_seq_valid
     */
    InstSeqNum fetch_seq = 0;

    /** @defgroup commit_seq
     * This records the instruction number that was committed in the pipeline
     * @see cp_seq_valid
     */
    InstSeqNum cp_seq = 0;

    /** @ingroup data
     * What size of data was written?
     */
    enum DataStatus
    {
        DataInvalid = 0,
        DataInt8 = 1,   // set to equal number of bytes
        DataInt16 = 2,
        DataInt32 = 4,
        DataInt64 = 8,
        DataDouble = 3,
        DataVec = 5,
        DataVecPred = 6,
        DataKeccak = 7
    } data_status = DataInvalid;

    /** @ingroup memory
     * Are the memory fields in the record valid?
     */
    bool mem_valid = false;

    /** @ingroup fetch_seq
     * Are the fetch sequence number fields valid?
     */
    bool fetch_seq_valid = false;
    /** @ingroup commit_seq
     * Are the commit sequence number fields valid?
     */
    bool cp_seq_valid = false;

    /** is the predicate for execution this inst true or false (not execed)?
     */
    bool predicate = true;

    /**
     * Did the execution of this instruction fault? (requires ExecFaulting
     * to be enabled)
     */
    bool faulting = false;

  public:
    InstRecord(Tick _when, ThreadContext *_thread,
               const StaticInstPtr _staticInst, const PCStateBase &_pc,
               const StaticInstPtr _macroStaticInst=nullptr)
        : when(_when), thread(_thread), staticInst(_staticInst),
        pc(_pc.clone()), macroStaticInst(_macroStaticInst)
    {}

    virtual ~InstRecord()
    {
        if (data_status == DataVec) {
            assert(data.as_vec);
            delete data.as_vec;
        } else if (data_status == DataVecPred) {
            assert(data.as_pred);
            delete data.as_pred;
        } else if (data_status == DataKeccak) {
            assert(data.as_keccak_reg);
            delete data.as_keccak_reg;
        }
    }

    void setWhen(Tick new_when) { when = new_when; }
    void
    setMem(Addr a, Addr s, unsigned f)
    {
        addr = a;
        size = s;
        flags = f;
        mem_valid = true;
    }

    template <typename T, size_t N>
    void
    setData(std::array<T, N> d)
    {
        data.as_int = d[0];
        data_status = (DataStatus)sizeof(T);
        static_assert(sizeof(T) == DataInt8 || sizeof(T) == DataInt16 ||
                      sizeof(T) == DataInt32 || sizeof(T) == DataInt64,
                      "Type T has an unrecognized size.");
    }

    void
    setData(uint64_t d)
    {
        data.as_int = d;
        data_status = DataInt64;
    }
    void
    setData(uint32_t d)
    {
        data.as_int = d;
        data_status = DataInt32;
    }
    void
    setData(uint16_t d)
    {
        data.as_int = d;
        data_status = DataInt16;
    }
    void
    setData(uint8_t d)
    {
        data.as_int = d;
        data_status = DataInt8;
    }

    void setData(int64_t d) { setData((uint64_t)d); }
    void setData(int32_t d) { setData((uint32_t)d); }
    void setData(int16_t d) { setData((uint16_t)d); }
    void setData(int8_t d)  { setData((uint8_t)d); }

    void
    setData(double d)
    {
        data.as_double = d;
        data_status = DataDouble;
    }

    void
    setData(TheISA::VecRegContainer& d)
    {
        data.as_vec = new TheISA::VecRegContainer(d);
        data_status = DataVec;
    }

    void
    setData(TheISA::VecPredRegContainer& d)
    {
        data.as_pred = new TheISA::VecPredRegContainer(d);
        data_status = DataVecPred;
    }

    void
    setData(TheISA::KeccakRegContainer& d)
    {
        data.as_keccak_reg = new TheISA::KeccakRegContainer(d);
        data_status = DataKeccak;
    }    

    void
    setFetchSeq(InstSeqNum seq)
    {
        fetch_seq = seq;
        fetch_seq_valid = true;
    }

    void
    setCPSeq(InstSeqNum seq)
    {
        cp_seq = seq;
        cp_seq_valid = true;
    }

    void setPredicate(bool val) { predicate = val; }

    void setFaulting(bool val) { faulting = val; }

    virtual void dump() = 0;

  public:
    Tick getWhen() const { return when; }
    ThreadContext *getThread() const { return thread; }
    StaticInstPtr getStaticInst() const { return staticInst; }
    const PCStateBase &getPCState() const { return *pc; }
    StaticInstPtr getMacroStaticInst() const { return macroStaticInst; }

    Addr getAddr() const { return addr; }
    Addr getSize() const { return size; }
    unsigned getFlags() const { return flags; }
    bool getMemValid() const { return mem_valid; }

    uint64_t getIntData() const { return data.as_int; }
    double getFloatData() const { return data.as_double; }
    int getDataStatus() const { return data_status; }

    InstSeqNum getFetchSeq() const { return fetch_seq; }
    bool getFetchSeqValid() const { return fetch_seq_valid; }

    InstSeqNum getCpSeq() const { return cp_seq; }
    bool getCpSeqValid() const { return cp_seq_valid; }

    bool getFaulting() const { return faulting; }
};

class InstTracer : public SimObject
{
  public:
    InstTracer(const Params &p) : SimObject(p) {}

    virtual ~InstTracer() {}

    virtual InstRecord *
        getInstRecord(Tick when, ThreadContext *tc,
                const StaticInstPtr staticInst, const PCStateBase &pc,
                const StaticInstPtr macroStaticInst=nullptr) = 0;
};

} // namespace Trace
} // namespace gem5

#endif // __INSTRECORD_HH__
