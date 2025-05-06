#ifndef __ARCH_RISCV_PQC_UTILITY_HH__
#define __ARCH_RISCV_PQC_UTILITY_HH__

#include <cmath>
#include <cstdint>
#include <sstream>
#include <string>
#include <utility>

#include "arch/riscv/regs/float.hh"
#include "arch/riscv/regs/int.hh"
#include "arch/riscv/regs/vector.hh"
#include "base/types.hh"
#include "cpu/reg_class.hh"
#include "cpu/static_inst.hh"
#include "cpu/thread_context.hh"
#include "rvk.hh"

namespace gem5
{

namespace RiscvISA
{
/*************************************************
* Name:        central_reduce
*
* Description: For finite field element a, compute r \quiv a (mod Q)
*              such that -Q < r < Q 
* Input Limitation: -2Q < a < 2Q
**************************************************/
template <typename ElemType>
ElemType central_reduce(ElemType a, const ElemType modulus_q) {
    ElemType r = a;
    if(r >= modulus_q) {
        r -= modulus_q;
    }
    else if(r <= -modulus_q) {
        r += modulus_q;
    }

    return r;
}

/*************************************************
* Name:        positive_reduce
*
* Description: For finite field element a, compute r \quiv a (+mod Q)
*              such that 0 <= r < Q 
* Input Limitation: -Q <= a < Q
**************************************************/
template <typename ElemType>
ElemType positive_reduce(ElemType a, const ElemType modulus_q) {
    ElemType r = a;
    if(r >= modulus_q) {
        r -= modulus_q;
    }
    else if(r < 0) {
        r += modulus_q;
    }

    panic_if(a >= (2*modulus_q), "input of positive_reduce should in [-q, 2q)");
    panic_if(a < -modulus_q, "input of positive_reduce should in [-q, 2q)");
    return r;
}

/*************************************************
* Name:        modular_add
*
* Description: For finite field element a, b, compute r \quiv a+b (+mod Q)
*              such that 0 <= r < Q 
* Input Limitation: 
* * condition1: 0 <= a, b < Q
* * condition2: a = q, and -Q <= b < 0
* * condition3: b = q, and -Q <= a < 0
**************************************************/
template <typename ElemType>
ElemType modular_add(ElemType a, ElemType b, const ElemType modulus_q) {

    panic_if((a < 0) && (b != modulus_q), "when a in [-q, 0), b must be q!");
    panic_if((b < 0) && (a != modulus_q), "when b in [-q, 0), a must be q!");
    ElemType r = positive_reduce<ElemType>(a+b, modulus_q);
    return r;
}

/*************************************************
* Name:        modular_sub
*
* Description: For finite field element a, b, compute r \quiv a-b (+mod Q)
*              such that 0 <= r < Q 
* Input Limitation: 0 <= a,b < Q
**************************************************/
template <typename ElemType>
ElemType modular_sub(ElemType a, ElemType b, const ElemType modulus_q) {
    ElemType r = positive_reduce<ElemType>(a - b, modulus_q);
    return r;
}

/*************************************************
* Name:        montgomery_reduce
*
* Description: For finite field element a, compute r \quiv a * R^{-1} (+mod Q)
*              such that 0 <= r <= Q-1
* Note:        q_inv should be -q^{-1} +mod 2^16
**************************************************/
template <typename ElemType, typename DoubleElemType>
ElemType montgomery_reduce(DoubleElemType a, const ElemType modulus_q, const ElemType q_inv) {
    DoubleElemType t;
    ElemType u;
    uint16_t right_shift_num = sizeof(ElemType) * 8;
    
    u = (ElemType)a * q_inv;
    t = (DoubleElemType)u * modulus_q;
    t = a + t;
    t >>= right_shift_num;

    if(t < 0) {
        t += modulus_q;
    }
    else if(t >= modulus_q) {
        t -= modulus_q;
    }

    panic_if(t >= modulus_q, "output of montgomery_reduce should in [0, q)");
    panic_if(t < 0, "output of montgomery_reduce should in [0, q)");
    return (ElemType)t;
}

/*************************************************
* Name:        montgomery_mul
*
* Description: modular multiplication with Montgomery reduction; R = 2^16 or 2^32
* Return: integer in {0. ..., q-1} congruent to (a*b)*R^{-1} modulo q
* Input Limitation: 0 <= a,b < Q
**************************************************/
template <typename ElemType, typename DoubleElemType>
ElemType montgomery_mul(ElemType a, ElemType b, const ElemType modulus_q, const ElemType q_inv) {
    panic_if(a >= modulus_q, "input a of montgomery_mul should in [0, q)");
    panic_if(a < 0, "input a of montgomery_mul should in [0, q)");
    panic_if(b >= modulus_q, "input b of montgomery_mul should in [0, q)");
    panic_if(b < 0, "input b of montgomery_mul should in [0, q)");

    DoubleElemType c = (DoubleElemType)a * b;
    ElemType r = montgomery_reduce<ElemType, DoubleElemType>(c, modulus_q, q_inv);
    return r;
}

/*************************************************
* Name:        modular_div2
*
* Description: modular divide 2
* Return: if(x&1 == 1) (x>>1 + (Q+1)/2) else x>>1
**************************************************/
template <typename ElemType>
ElemType modular_div2(ElemType a, const ElemType modulus_q) {
    ElemType r = (a >> 1) + (a & 1) * ((modulus_q + 1) / 2);
    return r;
}

template <typename ElemType, typename DoubleElemType>
void butterfly_ct(ElemType r[2], ElemType u, ElemType v, ElemType zeta, const ElemType modulus_q, const ElemType q_inv) {
    ElemType temp = montgomery_mul<ElemType, DoubleElemType>(v, zeta, modulus_q, q_inv);
    r[0] = modular_add<ElemType>(u, temp, modulus_q);
    r[1] = modular_sub<ElemType>(u, temp, modulus_q);
}

template <typename ElemType, typename DoubleElemType>
void butterfly_gs(ElemType r[2], ElemType u, ElemType v, ElemType zeta, const ElemType modulus_q, const ElemType q_inv) {
    ElemType t0 = modular_add<ElemType>(u, v, modulus_q);
    ElemType t1 = modular_sub<ElemType>(u, v, modulus_q);

    r[0] = modular_div2<ElemType>(t0, modulus_q);
    r[1] = modular_div2<ElemType>(t1, modulus_q);

    r[1] = montgomery_mul<ElemType, DoubleElemType>(r[1], zeta, modulus_q, q_inv);
}


/*************************************************
* Name:        basemul 2x2
*
* Description: Multiplication of polynomials in Zq[X]/(X^2-zeta)
*              used for multiplication of elements in Rq in NTT domain
*
* Arguments:   - ElemType r[2]:       pointer to the output polynomial
*              - const ElemType a[2]: pointer to the first factor
*              - const ElemType b[2]: pointer to the second factor
*              - ElemType zeta:       integer defining the reduction polynomial
**************************************************/
template <typename ElemType, typename DoubleElemType>
void basemul_2x2(ElemType r[2], ElemType a[2], ElemType b[2], ElemType zeta, const ElemType modulus_q, const ElemType q_inv)
{
    // use Karatsuba algorithm to optimize computation
    // 4 mod mul and 5 mod add/sub
    ElemType s0, s1, m0, m1;
    s0 = modular_add<ElemType>(a[0], a[1], modulus_q);
    s1 = modular_add<ElemType>(b[0], b[1], modulus_q);
    m0 = montgomery_mul<ElemType, DoubleElemType>(a[0], b[0], modulus_q, q_inv);
    m1 = montgomery_mul<ElemType, DoubleElemType>(a[1], b[1], modulus_q, q_inv);

    ElemType s2, m2, m3;
    s2 = modular_add<ElemType>(m0, m1, modulus_q);
    m2 = montgomery_mul<ElemType, DoubleElemType>(s0, s1, modulus_q, q_inv);
    m3 = montgomery_mul<ElemType, DoubleElemType>(m1, zeta, modulus_q, q_inv);

    r[0] = modular_add<ElemType>(m0, m3, modulus_q);
    r[1] = modular_sub<ElemType>(m2, s2, modulus_q);    
}


inline uint16_t getByteLowBound (uint16_t elemIdx, uint16_t Actual_Width)
{
    return (elemIdx * Actual_Width) >> 3;
}
inline uint16_t getByteUpBound (uint16_t elemIdx, uint16_t Actual_Width)
{
    return ((elemIdx+1) * Actual_Width) >> 3;
}
inline uint16_t getFirstShiftOffset (uint16_t elemIdx, uint16_t byteIdx, uint16_t Actual_Width)
{
    return (elemIdx * Actual_Width - byteIdx * 8);
}
inline uint16_t getOtherShiftOffset (uint16_t elemIdx, uint16_t byteIdx, uint16_t Actual_Width)
{
    return (byteIdx * 8 - elemIdx * Actual_Width);
}
inline bool isPackedFirstInByte (uint16_t elemIdx, uint16_t byteIdx, uint16_t Actual_Width)
{
    // mainly designed to cope with situation that when packing a vecreg into vd, vd may not be empty in advance
    return ((byteIdx * 8 / Actual_Width) == elemIdx);
}
inline bool byteUpBoundDivisible (uint16_t elemIdx, uint16_t Actual_Width)
{
    return ((((elemIdx + 1) * Actual_Width) & 0b111) == 0);
}

/*************************************************
* Name:        clmul
*
* Description: Carry-less multiplication (vector-scalar)
**************************************************/
template <typename ElemType>
void clmul(ElemType vd[], ElemType vs2, ElemType rs1) {
    ElemType temp[64][2]={0};
    uint8_t len=sizeof(ElemType)<<3;
    for (int i = 0; i < len; i++) {
        if ((rs1 & 1) == 1) {
            temp[i][0] = i == 0 ? vs2 : vs2 << i;
            temp[i][1] = i==0?0: vs2 >> (len - i);
        }
        rs1 >>= 1;
    }

    for (int j = (len>>1); j > 0; j >>= 1) {
        for (int i = 0; i < j; i++) {
            temp[i][0] ^= temp[i + j][0];
            temp[i][1] ^= temp[i + j][1];
        }
    }
    vd[0] = temp[0][0];
    vd[1] = temp[0][1];
}

uint16_t gf_mul(uint16_t a, uint16_t b, uint16_t prim);
uint16_t gf_inv(uint16_t a, uint16_t prim);

}

}

#endif // __ARCH_RISCV_PQC_UTILITY_HH__