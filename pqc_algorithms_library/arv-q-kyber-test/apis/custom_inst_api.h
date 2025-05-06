#ifndef CUSTOM_INST_API_H
#define CUSTOM_INST_API_H

#include <stddef.h>
#include <assert.h>
#include <riscv_vector.h>

#define VLEN 1024

#if (VLEN == 256)
#define ELEMENT_SEW8_PER_VECREG 32
#define ELEMENT_SEW16_PER_VECREG 16
#define ELEMENT_SEW32_PER_VECREG 8
#define ELEMENT_SEW64_PER_VECREG 4
#elif (VLEN == 512)
#define ELEMENT_SEW8_PER_VECREG 64
#define ELEMENT_SEW16_PER_VECREG 32
#define ELEMENT_SEW32_PER_VECREG 16
#define ELEMENT_SEW64_PER_VECREG 8
#elif (VLEN == 1024)
#define ELEMENT_SEW8_PER_VECREG 128
#define ELEMENT_SEW16_PER_VECREG 64
#define ELEMENT_SEW32_PER_VECREG 32
#define ELEMENT_SEW64_PER_VECREG 16
# else
#error "VLEN must be 256/512/1024"
#endif



#define SHAKE_128_MODE 0
#define SHAKE_256_MODE 1
#define SHA3_256_MODE  2
#define SHA3_512_MODE  3

#define GET_ZETA_MODE(SAME_TIMES, REPEAT_TIMES) ((SAME_TIMES&0xff) | (REPEAT_TIMES<<8))

static inline uint64_t read_cycle() {
    uint64_t cycle;
    __asm__ __volatile__ ( "rdcycle %[dest]" : [dest] "=r"(cycle): );

    return cycle;
}



/*
**
**  CSR Read-Write Instruction 
** 
*/
static inline int32_t csr_modulusq_rw(int32_t modulusq) {
    int32_t old_modulusq;
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]"
    : [rd] "=r"(old_modulusq)
    : [rt] "r"(modulusq), [csrnum] "i"(0)
    );
    return old_modulusq;
}

static inline int32_t csr_qinv_rw(int32_t qinv) {
    int32_t old_qinv;
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]"
    : [rd] "=r"(old_qinv)
    : [rt] "r"(qinv), [csrnum] "i"(1)
    );
    return old_qinv;    
}

static inline uint8_t csr_keccakmode_rw(uint8_t keccakmode) {
    uint8_t old_keccakmode;
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]"
    : [rd] "=r"(old_keccakmode)
    : [rt] "r"(keccakmode), [csrnum] "i"(2)
    );
    return old_keccakmode;     
}

static inline uint32_t csr_validnum_rw(void) {
    uint32_t validnum;
    uint32_t validnum_reset = 0;
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]"
    : [rd] "=r"(validnum)
    : [rt] "r"(validnum_reset), [csrnum] "i"(3)
    );
    return validnum;     
}

static inline uint16_t csr_zeta_selMode_rw(uint16_t zeta_selMode) {
    uint16_t old_zeta_selMode;
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]"
    : [rd] "=r"(old_zeta_selMode)
    : [rt] "r"(zeta_selMode), [csrnum] "i"(4)
    );
    return old_zeta_selMode;
}

// in FALCON512, set is to 9, while in FALCON1024, set it to 10
static inline uint16_t csr_selGaussMode_rw(uint16_t gauss_mode) {
    assert((gauss_mode == 9) || (gauss_mode == 10));
    uint16_t old_gauss_Mode;
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]"
    : [rd] "=r"(old_gauss_Mode)
    : [rt] "r"(gauss_mode), [csrnum] "i"(5)
    );
    return old_gauss_Mode;
}

static inline uint16_t csr_BSeqSftOffset_rw(uint16_t bSeqSftOffset) {
    assert((bSeqSftOffset >= 1) && (bSeqSftOffset <= 63));
    uint16_t old_bSeqSftOffset;
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]"
    : [rd] "=r"(old_bSeqSftOffset)
    : [rt] "r"(bSeqSftOffset), [csrnum] "i"(6)
    );
    return old_bSeqSftOffset;
}

/*
**
**  Gaussian sampling instruction
**
*/
// mu is fs, aka rs1, srcv1, srcf1;
// isigma is ft, aka rs2, srcv0, srcf0;
static inline double samplerz_double(double mu, double isigma) {
    double res;
    __asm__ __volatile__ ( "samplerz %[fd], %[fs], %[ft]"
    : [fd] "=f"(res)
    : [fs] "f"(mu), [ft] "f"(isigma)
    );
    return res;
}


/*
**
**  Vector Butterfly Instruction 
** 
*/

// vx is BFU input a(top half), vy is BFU input b(bottom half)
static inline vint16m2_t vbutterfly_ct_vvm_i16m2(vint16m1_t vx, vint16m1_t vy, vint16m1_t zetas) {
    vint16m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint16m4_t vbutterfly_ct_vvm_i16m4(vint16m2_t vx, vint16m2_t vy, vint16m1_t zetas) {
    vint16m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint16m8_t vbutterfly_ct_vvm_i16m8(vint16m4_t vx, vint16m4_t vy, vint16m1_t zetas) {
    vint16m8_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}


static inline vint32m2_t vbutterfly_ct_vvm_i32m2(vint32m1_t vx, vint32m1_t vy, vint16m1_t zetas) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint32m4_t vbutterfly_ct_vvm_i32m4(vint32m2_t vx, vint32m2_t vy, vint16m1_t zetas) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint32m8_t vbutterfly_ct_vvm_i32m8(vint32m4_t vx, vint32m4_t vy, vint16m1_t zetas) {
    vint32m8_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}


static inline vint16m2_t vbutterfly_gs_vvm_i16m2(vint16m1_t vx, vint16m1_t vy, vint16m1_t zetas) {
    vint16m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;
}

static inline vint16m4_t vbutterfly_gs_vvm_i16m4(vint16m2_t vx, vint16m2_t vy, vint16m1_t zetas) {
    vint16m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint16m8_t vbutterfly_gs_vvm_i16m8(vint16m4_t vx, vint16m4_t vy, vint16m1_t zetas) {
    vint16m8_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint32m2_t vbutterfly_gs_vvm_i32m2(vint32m1_t vx, vint32m1_t vy, vint16m1_t zetas) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint32m4_t vbutterfly_gs_vvm_i32m4(vint32m2_t vx, vint32m2_t vy, vint16m1_t zetas) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint32m8_t vbutterfly_gs_vvm_i32m8(vint32m4_t vx, vint32m4_t vy, vint16m1_t zetas) {
    vint32m8_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

// vx is BFU input a(top half), vy is BFU input b(bottom half)
static inline vuint16m2_t vbutterfly_ct_vvm_u16m2(vuint16m1_t vx, vuint16m1_t vy, vuint16m1_t zetas) {
    vuint16m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vuint16m4_t vbutterfly_ct_vvm_u16m4(vuint16m2_t vx, vuint16m2_t vy, vuint16m1_t zetas) {
    vuint16m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vuint16m8_t vbutterfly_ct_vvm_u16m8(vuint16m4_t vx, vuint16m4_t vy, vuint16m1_t zetas) {
    vuint16m8_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}


static inline vuint32m2_t vbutterfly_ct_vvm_u32m2(vuint32m1_t vx, vuint32m1_t vy, vuint16m1_t zetas) {
    vuint32m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vuint32m4_t vbutterfly_ct_vvm_u32m4(vuint32m2_t vx, vuint32m2_t vy, vuint16m1_t zetas) {
    vuint32m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vuint32m8_t vbutterfly_ct_vvm_u32m8(vuint32m4_t vx, vuint32m4_t vy, vuint16m1_t zetas) {
    vuint32m8_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}


static inline vuint16m2_t vbutterfly_gs_vvm_u16m2(vuint16m1_t vx, vuint16m1_t vy, vuint16m1_t zetas) {
    vuint16m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;
}

static inline vuint16m4_t vbutterfly_gs_vvm_u16m4(vuint16m2_t vx, vuint16m2_t vy, vuint16m1_t zetas) {
    vuint16m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vuint16m8_t vbutterfly_gs_vvm_u16m8(vuint16m4_t vx, vuint16m4_t vy, vuint16m1_t zetas) {
    vuint16m8_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vuint32m2_t vbutterfly_gs_vvm_u32m2(vuint32m1_t vx, vuint32m1_t vy, vuint16m1_t zetas) {
    vuint32m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vuint32m4_t vbutterfly_gs_vvm_u32m4(vuint32m2_t vx, vuint32m2_t vy, vuint16m1_t zetas) {
    vuint32m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vuint32m8_t vbutterfly_gs_vvm_u32m8(vuint32m4_t vx, vuint32m4_t vy, vuint16m1_t zetas) {
    vuint32m8_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}


/*
**
**  Vector-vector modular Instruction 
** 
*/
static inline vint16m1_t vmod_add_vv_i16m1(vint16m1_t vx, vint16m1_t vy) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "vmod.add.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m1_t vmod_add_vv_i16m1_m(vbool16_t mask, vint16m1_t vx, vint16m1_t vy) {
    vint16m1_t vz;
    __asm__ __volatile__ (
    "vmod.add.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint16m2_t vmod_add_vv_i16m2(vint16m2_t vx, vint16m2_t vy) {
    vint16m2_t vz;
    __asm__ __volatile__ ( "vmod.add.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m2_t vmod_add_vv_i16m2_m(vbool8_t mask, vint16m2_t vx, vint16m2_t vy) {
    vint16m2_t vz;
    __asm__ __volatile__ (
    "vmod.add.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint16m4_t vmod_add_vv_i16m4(vint16m4_t vx, vint16m4_t vy) {
    vint16m4_t vz;
    __asm__ __volatile__ ( "vmod.add.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m4_t vmod_add_vv_i16m4_m(vbool4_t mask, vint16m4_t vx, vint16m4_t vy) {
    vint16m4_t vz;
    __asm__ __volatile__ (
    "vmod.add.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint16m8_t vmod_add_vv_i16m8(vint16m8_t vx, vint16m8_t vy) {
    vint16m8_t vz;
    __asm__ __volatile__ ( "vmod.add.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m8_t vmod_add_vv_i16m8_m(vbool2_t mask, vint16m8_t vx, vint16m8_t vy) {
    vint16m8_t vz;
    __asm__ __volatile__ (
    "vmod.add.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m1_t vmod_add_vv_i32m1(vint32m1_t vx, vint32m1_t vy) {
    vint32m1_t vz;
    __asm__ __volatile__ ( "vmod.add.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m1_t vmod_add_vv_i32m1_m(vbool32_t mask, vint32m1_t vx, vint32m1_t vy) {
    vint32m1_t vz;
    __asm__ __volatile__ (
    "vmod.add.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint32m2_t vmod_add_vv_i32m2(vint32m2_t vx, vint32m2_t vy) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vmod.add.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m2_t vmod_add_vv_i32m2_m(vbool16_t mask, vint32m2_t vx, vint32m2_t vy) {
    vint32m2_t vz;
    __asm__ __volatile__ (
    "vmod.add.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint32m4_t vmod_add_vv_i32m4(vint32m4_t vx, vint32m4_t vy) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vmod.add.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m4_t vmod_add_vv_i32m4_m(vbool8_t mask, vint32m4_t vx, vint32m4_t vy) {
    vint32m4_t vz;
    __asm__ __volatile__ (
    "vmod.add.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m8_t vmod_add_vv_i32m8(vint32m8_t vx, vint32m8_t vy) {
    vint32m8_t vz;
    __asm__ __volatile__ ( "vmod.add.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m8_t vmod_add_vv_i32m8_m(vbool4_t mask, vint32m8_t vx, vint32m8_t vy) {
    vint32m8_t vz;
    __asm__ __volatile__ (
    "vmod.add.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vuint16m1_t vmod_add_vv_u16m1(vuint16m1_t vx, vuint16m1_t vy) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "vmod.add.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

//***
// vmod sub is vd = vs2 - vs1， vx is vs2, vy is vs1
//***
static inline vint16m1_t vmod_sub_vv_i16m1(vint16m1_t vx, vint16m1_t vy) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "vmod.sub.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m1_t vmod_sub_vv_i16m1_m(vbool16_t mask, vint16m1_t vx, vint16m1_t vy) {
    vint16m1_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint16m2_t vmod_sub_vv_i16m2(vint16m2_t vx, vint16m2_t vy) {
    vint16m2_t vz;
    __asm__ __volatile__ ( "vmod.sub.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m2_t vmod_sub_vv_i16m2_m(vbool8_t mask, vint16m2_t vx, vint16m2_t vy) {
    vint16m2_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint16m4_t vmod_sub_vv_i16m4(vint16m4_t vx, vint16m4_t vy) {
    vint16m4_t vz;
    __asm__ __volatile__ ( "vmod.sub.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m4_t vmod_sub_vv_i16m4_m(vbool4_t mask, vint16m4_t vx, vint16m4_t vy) {
    vint16m4_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint16m8_t vmod_sub_vv_i16m8(vint16m8_t vx, vint16m8_t vy) {
    vint16m8_t vz;
    __asm__ __volatile__ ( "vmod.sub.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m8_t vmod_sub_vv_i16m8_m(vbool2_t mask, vint16m8_t vx, vint16m8_t vy) {
    vint16m8_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m1_t vmod_sub_vv_i32m1(vint32m1_t vx, vint32m1_t vy) {
    vint32m1_t vz;
    __asm__ __volatile__ ( "vmod.sub.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m1_t vmod_sub_vv_i32m1_m(vbool32_t mask, vint32m1_t vx, vint32m1_t vy) {
    vint32m1_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint32m2_t vmod_sub_vv_i32m2(vint32m2_t vx, vint32m2_t vy) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vmod.sub.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m2_t vmod_sub_vv_i32m2_m(vbool16_t mask, vint32m2_t vx, vint32m2_t vy) {
    vint32m2_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint32m4_t vmod_sub_vv_i32m4(vint32m4_t vx, vint32m4_t vy) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vmod.sub.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m4_t vmod_sub_vv_i32m4_m(vbool8_t mask, vint32m4_t vx, vint32m4_t vy) {
    vint32m4_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m8_t vmod_sub_vv_i32m8(vint32m8_t vx, vint32m8_t vy) {
    vint32m8_t vz;
    __asm__ __volatile__ ( "vmod.sub.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m8_t vmod_sub_vv_i32m8_m(vbool4_t mask, vint32m8_t vx, vint32m8_t vy) {
    vint32m8_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vuint16m1_t vmod_sub_vv_u16m1(vuint16m1_t vx, vuint16m1_t vy) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "vmod.sub.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m1_t vmod_mul_vv_i16m1(vint16m1_t vx, vint16m1_t vy) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "vmod.mul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m1_t vmod_mul_vv_i16m1_m(vbool16_t mask, vint16m1_t vx, vint16m1_t vy) {
    vint16m1_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)        
    );
    return vz;    
}

static inline vint16m2_t vmod_mul_vv_i16m2(vint16m2_t vx, vint16m2_t vy) {
    vint16m2_t vz;
    __asm__ __volatile__ ( "vmod.mul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m2_t vmod_mul_vv_i16m2_m(vbool8_t mask, vint16m2_t vx, vint16m2_t vy) {
    vint16m2_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint16m4_t vmod_mul_vv_i16m4(vint16m4_t vx, vint16m4_t vy) {
    vint16m4_t vz;
    __asm__ __volatile__ ( "vmod.mul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m4_t vmod_mul_vv_i16m4_m(vbool4_t mask, vint16m4_t vx, vint16m4_t vy) {
    vint16m4_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint16m8_t vmod_mul_vv_i16m8(vint16m8_t vx, vint16m8_t vy) {
    vint16m8_t vz;
    __asm__ __volatile__ ( "vmod.mul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint16m8_t vmod_mul_vv_i16m8_m(vbool2_t mask, vint16m8_t vx, vint16m8_t vy) {
    vint16m8_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m1_t vmod_mul_vv_i32m1(vint32m1_t vx, vint32m1_t vy) {
    vint32m1_t vz;
    __asm__ __volatile__ ( "vmod.mul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m1_t vmod_mul_vv_i32m1_m(vbool32_t mask, vint32m1_t vx, vint32m1_t vy) {
    vint32m1_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint32m2_t vmod_mul_vv_i32m2(vint32m2_t vx, vint32m2_t vy) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vmod.mul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m2_t vmod_mul_vv_i32m2_m(vbool16_t mask, vint32m2_t vx, vint32m2_t vy) {
    vint32m2_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint32m4_t vmod_mul_vv_i32m4(vint32m4_t vx, vint32m4_t vy) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vmod.mul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m4_t vmod_mul_vv_i32m4_m(vbool8_t mask, vint32m4_t vx, vint32m4_t vy) {
    vint32m4_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m8_t vmod_mul_vv_i32m8(vint32m8_t vx, vint32m8_t vy) {
    vint32m8_t vz;
    __asm__ __volatile__ ( "vmod.mul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vint32m8_t vmod_mul_vv_i32m8_m(vbool4_t mask, vint32m8_t vx, vint32m8_t vy) {
    vint32m8_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;
}

static inline vuint16m1_t vmod_mul_vv_u16m1(vuint16m1_t vx, vuint16m1_t vy) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "vmod.mul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}


/*
**
**  Vector-scalar modular Instruction 
** 
*/
static inline vint16m1_t vmod_add_vx_i16m1(vint16m1_t vx, int16_t op2) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m1_t vmod_add_vx_i16m1_m(vbool16_t mask, vint16m1_t vx, int16_t op2) {
    vint16m1_t vz;
    __asm__ __volatile__ (
    "vmod.add.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint16m2_t vmod_add_vx_i16m2(vint16m2_t vx, int16_t op2) {
    vint16m2_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m2_t vmod_add_vx_i16m2_m(vbool8_t mask, vint16m2_t vx, int16_t op2) {
    vint16m2_t vz;
    __asm__ __volatile__ (
    "vmod.add.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint16m4_t vmod_add_vx_i16m4(vint16m4_t vx, int16_t op2) {
    vint16m4_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m4_t vmod_add_vx_i16m4_m(vbool4_t mask, vint16m4_t vx, int16_t op2) {
    vint16m4_t vz;
    __asm__ __volatile__ (
    "vmod.add.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint16m8_t vmod_add_vx_i16m8(vint16m8_t vx, int16_t op2) {
    vint16m8_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m8_t vmod_add_vx_i16m8_m(vbool2_t mask, vint16m8_t vx, int16_t op2) {
    vint16m8_t vz;
    __asm__ __volatile__ (
    "vmod.add.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m1_t vmod_add_vx_i32m1(vint32m1_t vx, int32_t op2) {
    vint32m1_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m1_t vmod_add_vx_i32m1_m(vbool32_t mask, vint32m1_t vx, int32_t op2) {
    vint32m1_t vz;
    __asm__ __volatile__ (
    "vmod.add.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint32m2_t vmod_add_vx_i32m2(vint32m2_t vx, int32_t op2) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m2_t vmod_add_vx_i32m2_m(vbool16_t mask, vint32m2_t vx, int32_t op2) {
    vint32m2_t vz;
    __asm__ __volatile__ (
    "vmod.add.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m4_t vmod_add_vx_i32m4(vint32m4_t vx, int32_t op2) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m4_t vmod_add_vx_i32m4_m(vbool8_t mask, vint32m4_t vx, int32_t op2) {
    vint32m4_t vz;
    __asm__ __volatile__ (
    "vmod.add.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m8_t vmod_add_vx_i32m8(vint32m8_t vx, int32_t op2) {
    vint32m8_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m8_t vmod_add_vx_i32m8_m(vbool4_t mask, vint32m8_t vx, int32_t op2) {
    vint32m8_t vz;
    __asm__ __volatile__ (
    "vmod.add.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vuint16m1_t vmod_add_vx_u16m1(vuint16m1_t vx, uint16_t op2) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint16m1_t vmod_add_vx_u16m1_i16m1(vint16m1_t vx, uint16_t op2) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "vmod.add.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}


static inline vint16m1_t vmod_sub_vx_i16m1(vint16m1_t vx, int16_t op2) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "vmod.sub.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m1_t vmod_sub_vx_i16m1_m(vbool16_t mask, vint16m1_t vx, int16_t op2) {
    vint16m1_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint16m2_t vmod_sub_vx_i16m2(vint16m2_t vx, int16_t op2) {
    vint16m2_t vz;
    __asm__ __volatile__ ( "vmod.sub.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m2_t vmod_sub_vx_i16m2_m(vbool8_t mask, vint16m2_t vx, int16_t op2) {
    vint16m2_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint16m4_t vmod_sub_vx_i16m4(vint16m4_t vx, int16_t op2) {
    vint16m4_t vz;
    __asm__ __volatile__ ( "vmod.sub.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m4_t vmod_sub_vx_i16m4_m(vbool4_t mask, vint16m4_t vx, int16_t op2) {
    vint16m4_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint16m8_t vmod_sub_vx_i16m8(vint16m8_t vx, int16_t op2) {
    vint16m8_t vz;
    __asm__ __volatile__ ( "vmod.sub.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m8_t vmod_sub_vx_i16m8_m(vbool2_t mask, vint16m8_t vx, int16_t op2) {
    vint16m8_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m1_t vmod_sub_vx_i32m1(vint32m1_t vx, int32_t op2) {
    vint32m1_t vz;
    __asm__ __volatile__ ( "vmod.sub.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m1_t vmod_sub_vx_i32m1_m(vbool32_t mask, vint32m1_t vx, int32_t op2) {
    vint32m1_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint32m2_t vmod_sub_vx_i32m2(vint32m2_t vx, int32_t op2) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vmod.sub.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m2_t vmod_sub_vx_i32m2_m(vbool16_t mask, vint32m2_t vx, int32_t op2) {
    vint32m2_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint32m4_t vmod_sub_vx_i32m4(vint32m4_t vx, int32_t op2) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vmod.sub.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m4_t vmod_sub_vx_i32m4_m(vbool8_t mask, vint32m4_t vx, int32_t op2) {
    vint32m4_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m8_t vmod_sub_vx_i32m8(vint32m8_t vx, int32_t op2) {
    vint32m8_t vz;
    __asm__ __volatile__ ( "vmod.sub.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m8_t vmod_sub_vx_i32m8_m(vbool4_t mask, vint32m8_t vx, int32_t op2) {
    vint32m8_t vz;
    __asm__ __volatile__ (
    "vmod.sub.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint16m1_t vmod_mul_vx_i16m1(vint16m1_t vx, int16_t op2) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m1_t vmod_mul_vx_i16m1_m(vbool16_t mask, vint16m1_t vx, int16_t op2) {
    vint16m1_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)        
    );
    return vz;    
}

static inline vint16m2_t vmod_mul_vx_i16m2(vint16m2_t vx, int16_t op2) {
    vint16m2_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m2_t vmod_mul_vx_i16m2_m(vbool8_t mask, vint16m2_t vx, int16_t op2) {
    vint16m2_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint16m4_t vmod_mul_vx_i16m4(vint16m4_t vx, int16_t op2) {
    vint16m4_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m4_t vmod_mul_vx_i16m4_m(vbool4_t mask, vint16m4_t vx, int16_t op2) {
    vint16m4_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint16m8_t vmod_mul_vx_i16m8(vint16m8_t vx, int16_t op2) {
    vint16m8_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint16m8_t vmod_mul_vx_i16m8_m(vbool2_t mask, vint16m8_t vx, int16_t op2) {
    vint16m8_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m1_t vmod_mul_vx_i32m1(vint32m1_t vx, int32_t op2) {
    vint32m1_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m1_t vmod_mul_vx_i32m1_m(vbool32_t mask, vint32m1_t vx, int32_t op2) {
    vint32m1_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vint32m2_t vmod_mul_vx_i32m2(vint32m2_t vx, int32_t op2) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m2_t vmod_mul_vx_i32m2_m(vbool16_t mask, vint32m2_t vx, int32_t op2) {
    vint32m2_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}


static inline vint32m4_t vmod_mul_vx_i32m4(vint32m4_t vx, int32_t op2) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m4_t vmod_mul_vx_i32m4_m(vbool8_t mask, vint32m4_t vx, int32_t op2) {
    vint32m4_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vint32m8_t vmod_mul_vx_i32m8(vint32m8_t vx, int32_t op2) {
    vint32m8_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vint32m8_t vmod_mul_vx_i32m8_m(vbool4_t mask, vint32m8_t vx, int32_t op2) {
    vint32m8_t vz;
    __asm__ __volatile__ (
    "vmod.mul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;
}

static inline vuint16m1_t vmod_mul_vx_u16m1(vuint16m1_t vx, uint16_t op2) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint32m2_t vmod_mul_vx_u32m2(vuint32m2_t vx, uint32_t op2) {
    vuint32m2_t vz;
    __asm__ __volatile__ ( "vmod.mul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

/*
**
**  Modular Basemul Instruction
** 
*/
static inline vint16m1_t vmod_basemul_vvm_i16m1(vint16m1_t vx, vint16m1_t vy, vint16m1_t vzeta) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "vmod.basemul.vvm %[vd], %[vt], %[vs], %[vzeta]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vzeta] "vm"(vzeta)
    );
    return vz;    
}


/*
**
**  Rejection Sample Instruction
** 
*/
static inline vint8m1_t sample_rej_vx_i8m1(vint8m1_t vx, int8_t bound) {
    vint8m1_t vz;
    __asm__ __volatile__ ( "sample.rej.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(bound)
    );
    return vz;    
}

static inline vint16m1_t sample_rej_vx_i16m1(vint16m1_t vx, int16_t bound) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "sample.rej.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(bound)
    );
    return vz;
}

static inline vint32m1_t sample_rej_vx_i32m1(vint32m1_t vx, int32_t bound) {
    vint32m1_t vz;
    __asm__ __volatile__ ( "sample.rej.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(bound)
    );
    return vz;    
}

static inline vint64m1_t sample_rej_vx_i64m1(vint64m1_t vx, int64_t bound) {
    vint64m1_t vz;
    __asm__ __volatile__ ( "sample.rej.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(bound)
    );
    return vz;    
}

static inline vuint8m1_t sample_rej_vx_u8m1(vuint8m1_t vx, uint8_t bound) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "sample.rej.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(bound)
    );
    return vz;    
}

static inline vuint16m1_t sample_rej_vx_u16m1(vuint16m1_t vx, uint16_t bound) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "sample.rej.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(bound)
    );
    return vz;
}

static inline vuint32m1_t sample_rej_vx_u32m1(vuint32m1_t vx, uint32_t bound) {
    vuint32m1_t vz;
    __asm__ __volatile__ ( "sample.rej.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(bound)
    );
    return vz;    
}

static inline vuint64m1_t sample_rej_vx_u64m1(vuint64m1_t vx, uint64_t bound) {
    vuint64m1_t vz;
    __asm__ __volatile__ ( "sample.rej.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(bound)
    );
    return vz;    
}

/*
**
**  Pack Instruction
** 
*/
static inline vuint8m1_t pack_vx_i8m1(vint8m1_t vx, uint32_t validnum) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "pack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint8m1_t pack_vx_i16m1(vint16m1_t vx, uint32_t validnum) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "pack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint8m1_t pack_vx_i32m1(vint32m1_t vx, uint32_t validnum) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "pack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint8m1_t pack_vx_i64m1(vint64m1_t vx, uint32_t validnum) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "pack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}


static inline vuint8m1_t pack_vx_u8m1(vuint8m1_t vx, uint32_t validnum) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "pack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint8m1_t pack_vx_u16m1(vuint16m1_t vx, uint32_t validnum) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "pack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint8m1_t pack_vx_u32m1(vuint32m1_t vx, uint32_t validnum) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "pack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint8m1_t pack_vx_u64m1(vuint64m1_t vx, uint32_t validnum) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "pack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}


/*
**
**  Unpack Instruction 
** 
*/
static inline vint8m1_t unpack_vx_i8m1(vuint8m1_t vx, uint8_t validnum) {
    vint8m1_t vz;
    __asm__ __volatile__ ( "unpack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vint16m1_t unpack_vx_i16m1(vuint8m1_t vx, uint16_t validnum) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "unpack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vint32m1_t unpack_vx_i32m1(vuint8m1_t vx, uint32_t validnum) {
    vint32m1_t vz;
    __asm__ __volatile__ ( "unpack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vint64m1_t unpack_vx_i64m1(vuint8m1_t vx, uint64_t validnum) {
    vint64m1_t vz;
    __asm__ __volatile__ ( "unpack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint8m1_t unpack_vx_u8m1(vuint8m1_t vx, uint8_t validnum) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "unpack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint16m1_t unpack_vx_u16m1(vuint8m1_t vx, uint16_t validnum) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "unpack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint32m1_t unpack_vx_u32m1(vuint8m1_t vx, uint32_t validnum) {
    vuint32m1_t vz;
    __asm__ __volatile__ ( "unpack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vuint64m1_t unpack_vx_u64m1(vuint8m1_t vx, uint64_t validnum) {
    vuint64m1_t vz;
    __asm__ __volatile__ ( "unpack.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}


static inline vint8m1_t unpacks_vx_i8m1(vuint8m1_t vx, uint8_t validnum) {
    vint8m1_t vz;
    __asm__ __volatile__ ( "unpacks.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vint16m1_t unpacks_vx_i16m1(vuint8m1_t vx, uint16_t validnum) {
    vint16m1_t vz;
    __asm__ __volatile__ ( "unpacks.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vint32m1_t unpacks_vx_i32m1(vuint8m1_t vx, uint32_t validnum) {
    vint32m1_t vz;
    __asm__ __volatile__ ( "unpacks.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

static inline vint64m1_t unpacks_vx_i64m1(vuint8m1_t vx, uint64_t validnum) {
    vint64m1_t vz;
    __asm__ __volatile__ ( "unpacks.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(validnum)
    );
    return vz;    
}

/*
**
**  Keccak Instruction 
** 
*/
static inline void keccak_loadm_v_u8m1(vuint8m1_t vx) {
    __asm__ __volatile__ ( "keccak.loadm.v %[vt]"
    :
    : [vt] "vr"(vx)
    );
}

static inline void keccak_loadm_v_u8m2(vuint8m2_t vx) {
    __asm__ __volatile__ ( "keccak.loadm.v %[vt]"
    :
    : [vt] "vr"(vx)
    );
}

static inline void keccak_loadm_v_u8m4(vuint8m4_t vx) {
    __asm__ __volatile__ ( "keccak.loadm.v %[vt]"
    :
    : [vt] "vr"(vx)
    );
}

static inline void keccak_loadm_v_u8m8(vuint8m8_t vx) {
    __asm__ __volatile__ ( "keccak.loadm.v %[vt]"
    :
    : [vt] "vr"(vx)
    );
}

static inline vuint8m1_t keccak_stores_v_u8m1(void) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "keccak.stores.v %[vd]"
    : [vd] "=vr" (vz)
    :
    );
    return vz;
}

static inline vuint8m2_t keccak_stores_v_u8m2(void) {
    vuint8m2_t vz;
    __asm__ __volatile__ ( "keccak.stores.v %[vd]"
    : [vd] "=vr" (vz)
    :
    );
    return vz;
}

static inline vuint8m4_t keccak_stores_v_u8m4(void) {
    vuint8m4_t vz;
    __asm__ __volatile__ ( "keccak.stores.v %[vd]"
    : [vd] "=vr" (vz)
    :
    );
    return vz;
}

static inline vuint8m8_t keccak_stores_v_u8m8(void) {
    vuint8m8_t vz;
    __asm__ __volatile__ ( "keccak.stores.v %[vd]"
    : [vd] "=vr" (vz)
    :
    );
    return vz;
}

static inline void keccak_init(void) {
    __asm__ __volatile__ ( "keccak.init %[imm]"
    : 
    : [imm] "i" (0)
    );
}

static inline size_t keccak_absorb_one_block(size_t length) {
    size_t left_length;
    __asm__ __volatile__ ( "keccak.absorb %[rd], %[rs]"
    : [rd] "=r" (left_length)
    : [rs] "r" (length)
    );
    return left_length;
}

static inline void keccak_squeeze(void) {
    __asm__ __volatile__ ( "keccak.squeeze %[imm]"
    : 
    : [imm] "i" (0)
    );
}

static inline vuint8m1_t vcpop_v_u8m1(vuint8m1_t vx) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "vcpop.v %[vd], %[vt]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx)
    );
    return vz;    
}

static inline vint8m1_t vcpop_v_i8m1(vint8m1_t vx) {
    vint8m1_t vz;
    __asm__ __volatile__ ( "vcpop.v %[vd], %[vt]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx)
    );
    return vz;    
}


/*
**
**  chacha20 and samplez related instructions
** 
*/
static inline void chacha20_init_v_u8m1(vuint8m1_t vx) {
    __asm__ __volatile__ ( "chacha20.init.v %[vt]"
    :
    : [vt] "vr"(vx)
    );
}

static inline void chacha20_init_v_u8m2(vuint8m2_t vx) {
    __asm__ __volatile__ ( "chacha20.init.v %[vt]"
    :
    : [vt] "vr"(vx)
    );
}

static inline void chacha20_init_v_u8m4(vuint8m4_t vx) {
    __asm__ __volatile__ ( "chacha20.init.v %[vt]"
    :
    : [vt] "vr"(vx)
    );
}

static inline void chacha20_init_v_u8m8(vuint8m8_t vx) {
    __asm__ __volatile__ ( "chacha20.init.v %[vt]"
    :
    : [vt] "vr"(vx)
    );
}


#endif  // CUSTOM_INST_API_H
