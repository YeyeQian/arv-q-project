#ifndef CUSTOM_INST_API_H
#define CUSTOM_INST_API_H

#include <stddef.h>
#include <riscv_vector.h>

#define VLEN 512                     

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
#define SHA3_384_MODE  4

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

static inline uint8_t csr_bseqoffset_rw(uint8_t offset) {
    uint8_t old_offset;
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]"
    : [rd] "=r"(old_offset)
    : [rt] "r"(offset), [csrnum] "i"(6)
    );
    return old_offset;
}

static inline uint16_t csr_primpoly_rw(uint16_t primpoly) {
    uint16_t old_primpoly;
    __asm__ __volatile__ ( "pqccsrrw %[rd], %[rt], %[csrnum]"
    : [rd] "=r"(old_primpoly)
    : [rt] "r"(primpoly), [csrnum] "i"(7)
    );
    return old_primpoly;
}

/*
**
**  Vector Butterfly Instruction 
** 
*/
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


static inline vint32m2_t vbutterfly_ct_vvm_i32m2(vint32m1_t vx, vint32m1_t vy, vint32m1_t zetas) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint32m4_t vbutterfly_ct_vvm_i32m4(vint32m2_t vx, vint32m2_t vy, vint32m1_t zetas) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.ct.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint32m8_t vbutterfly_ct_vvm_i32m8(vint32m4_t vx, vint32m4_t vy, vint32m1_t zetas) {
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

static inline vint32m2_t vbutterfly_gs_vvm_i32m2(vint32m1_t vx, vint32m1_t vy, vint32m1_t zetas) {
    vint32m2_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint32m4_t vbutterfly_gs_vvm_i32m4(vint32m2_t vx, vint32m2_t vy, vint32m1_t zetas) {
    vint32m4_t vz;
    __asm__ __volatile__ ( "vbutterfly.gs.vvm %[vd], %[vt], %[vs], %[zetas]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [zetas] "vm"(zetas)
    );
    return vz;    
}

static inline vint32m8_t vbutterfly_gs_vvm_i32m8(vint32m4_t vx, vint32m4_t vy, vint32m1_t zetas) {
    vint32m8_t vz;
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

// vmod sub is vd = vs2 - vs1， vx is vs2, vy is vs1
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
**  BIKE Related Instruction 
** 
*/
static inline vuint8m1_t bitseqsll_vx_u8m1(vuint8m1_t vx, uint8_t offset) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "bitseqsll.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(offset)
    );
    return vz;    
}

static inline vuint64m2_t clmul_vv_u64m2(vuint64m1_t vx, vuint64m1_t vy) {
    vuint64m2_t vz;
    __asm__ __volatile__ ( "clmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint64m4_t clmul_vv_u64m4(vuint64m2_t vx, vuint64m2_t vy) {
    vuint64m4_t vz;
    __asm__ __volatile__ ( "clmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint64m8_t clmul_vv_u64m8(vuint64m4_t vx, vuint64m4_t vy) {
    vuint64m8_t vz;
    __asm__ __volatile__ ( "clmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=&vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

/*
**
**  GF Operation Related Instruction 
** 
*/
static inline vuint8m1_t vgfmul_vv_u8m1(vuint8m1_t vx, vuint8m1_t vy) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "gfmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint8m1_t vgfmul_vv_u8m1_m(vbool8_t mask, vuint8m1_t vx, vuint8m1_t vy) {
    vuint8m1_t vz;
    __asm__ __volatile__ (
    "gfmul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint8m2_t vgfmul_vv_u8m2(vuint8m2_t vx, vuint8m2_t vy) {
    vuint8m2_t vz;
    __asm__ __volatile__ ( "gfmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint8m2_t vgfmul_vv_u8m2_m(vbool4_t mask, vuint8m2_t vx, vuint8m2_t vy) {
    vuint8m2_t vz;
    __asm__ __volatile__ (
    "gfmul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint8m4_t vgfmul_vv_u8m4(vuint8m4_t vx, vuint8m4_t vy) {
    vuint8m4_t vz;
    __asm__ __volatile__ ( "gfmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint8m4_t vgfmul_vv_u8m4_m(vbool2_t mask, vuint8m4_t vx, vuint8m4_t vy) {
    vuint8m4_t vz;
    __asm__ __volatile__ (
    "gfmul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint8m8_t vgfmul_vv_u8m8(vuint8m8_t vx, vuint8m8_t vy) {
    vuint8m8_t vz;
    __asm__ __volatile__ ( "gfmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint8m8_t vgfmul_vv_u8m8_m(vbool1_t mask, vuint8m8_t vx, vuint8m8_t vy) {
    vuint8m8_t vz;
    __asm__ __volatile__ (
    "gfmul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint16m1_t vgfmul_vv_u16m1(vuint16m1_t vx, vuint16m1_t vy) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "gfmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint16m1_t vgfmul_vv_u16m1_m(vbool16_t mask, vuint16m1_t vx, vuint16m1_t vy) {
    vuint16m1_t vz;
    __asm__ __volatile__ (
    "gfmul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint16m2_t vgfmul_vv_u16m2(vuint16m2_t vx, vuint16m2_t vy) {
    vuint16m2_t vz;
    __asm__ __volatile__ ( "gfmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint16m2_t vgfmul_vv_u16m2_m(vbool8_t mask, vuint16m2_t vx, vuint16m2_t vy) {
    vuint16m2_t vz;
    __asm__ __volatile__ (
    "gfmul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint16m4_t vgfmul_vv_u16m4(vuint16m4_t vx, vuint16m4_t vy) {
    vuint16m4_t vz;
    __asm__ __volatile__ ( "gfmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint16m4_t vgfmul_vv_u16m4_m(vbool4_t mask, vuint16m4_t vx, vuint16m4_t vy) {
    vuint16m4_t vz;
    __asm__ __volatile__ (
    "gfmul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint16m8_t vgfmul_vv_u16m8(vuint16m8_t vx, vuint16m8_t vy) {
    vuint16m8_t vz;
    __asm__ __volatile__ ( "gfmul.vv %[vd], %[vt], %[vs]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy)
    );
    return vz;    
}

static inline vuint16m8_t vgfmul_vv_u16m8_m(vbool2_t mask, vuint16m8_t vx, vuint16m8_t vy) {
    vuint16m8_t vz;
    __asm__ __volatile__ (
    "gfmul.vv %[vd], %[vt], %[vs], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [vs] "vr"(vy), [vm] "vm" (mask)
    );
    return vz;    
}


static inline vuint8m1_t vgfmul_vx_u8m1(vuint8m1_t vx, uint8_t op2) {
    vuint8m1_t vz;
    __asm__ __volatile__ ( "gfmul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint8m1_t vgfmul_vx_u8m1_m(vbool8_t mask, vuint8m1_t vx, uint8_t op2) {
    vuint8m1_t vz;
    __asm__ __volatile__ (
    "gfmul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint8m2_t vgfmul_vx_u8m2(vuint8m2_t vx, uint8_t op2) {
    vuint8m2_t vz;
    __asm__ __volatile__ ( "gfmul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint8m2_t vgfmul_vx_u8m2_m(vbool4_t mask, vuint8m2_t vx, uint8_t op2) {
    vuint8m2_t vz;
    __asm__ __volatile__ (
    "gfmul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint8m4_t vgfmul_vx_u8m4(vuint8m4_t vx, uint8_t op2) {
    vuint8m4_t vz;
    __asm__ __volatile__ ( "gfmul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint8m4_t vgfmul_vx_u8m4_m(vbool2_t mask, vuint8m4_t vx, uint8_t op2) {
    vuint8m4_t vz;
    __asm__ __volatile__ (
    "gfmul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint8m8_t vgfmul_vx_u8m8(vuint8m8_t vx, uint8_t op2) {
    vuint8m8_t vz;
    __asm__ __volatile__ ( "gfmul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint8m8_t vgfmul_vx_u8m8_m(vbool1_t mask, vuint8m8_t vx, uint8_t op2) {
    vuint8m8_t vz;
    __asm__ __volatile__ (
    "gfmul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}


static inline vuint16m1_t vgfmul_vx_u16m1(vuint16m1_t vx, uint16_t op2) {
    vuint16m1_t vz;
    __asm__ __volatile__ ( "gfmul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint16m1_t vgfmul_vx_u16m1_m(vbool16_t mask, vuint16m1_t vx, uint16_t op2) {
    vuint16m1_t vz;
    __asm__ __volatile__ (
    "gfmul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint16m2_t vgfmul_vx_u16m2(vuint16m2_t vx, uint16_t op2) {
    vuint16m2_t vz;
    __asm__ __volatile__ ( "gfmul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint16m2_t vgfmul_vx_u16m2_m(vbool8_t mask, vuint16m2_t vx, uint16_t op2) {
    vuint16m2_t vz;
    __asm__ __volatile__ (
    "gfmul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint16m4_t vgfmul_vx_u16m4(vuint16m4_t vx, uint16_t op2) {
    vuint16m4_t vz;
    __asm__ __volatile__ ( "gfmul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint16m4_t vgfmul_vx_u16m4_m(vbool4_t mask, vuint16m4_t vx, uint16_t op2) {
    vuint16m4_t vz;
    __asm__ __volatile__ (
    "gfmul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline vuint16m8_t vgfmul_vx_u16m8(vuint16m8_t vx, uint16_t op2) {
    vuint16m8_t vz;
    __asm__ __volatile__ ( "gfmul.vx %[vd], %[vt], %[op2]"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2)
    );
    return vz;    
}

static inline vuint16m8_t vgfmul_vx_u16m8_m(vbool2_t mask, vuint16m8_t vx, uint16_t op2) {
    vuint16m8_t vz;
    __asm__ __volatile__ (
    "gfmul.vx %[vd], %[vt], %[op2], %[vm].t"
    : [vd] "=vr"(vz)
    : [vt] "vr"(vx), [op2] "r"(op2), [vm] "vm" (mask)
    );
    return vz;    
}

static inline uint8_t gfinv_u8(uint8_t op1) {
    uint8_t res;
    __asm__ __volatile__ ( "gfinv %[rd], %[rs]"
    : [rd] "=r" (res)
    : [rs] "r" (op1)
    );
    return res;
}

static inline uint16_t gfinv_u16(uint16_t op1) {
    uint16_t res;
    __asm__ __volatile__ ( "gfinv %[rd], %[rs]"
    : [rd] "=r" (res)
    : [rs] "r" (op1)
    );
    return res;
}

static inline uint8_t gfmul_u8(uint8_t op1, uint8_t op2) {
    uint8_t res;
    __asm__ __volatile__ ( "gfmul %[rd], %[rs], %[rt]"
    : [rd] "=r" (res)
    : [rs] "r" (op1), [rt] "r" (op2)
    );
    return res;
}

static inline uint16_t gfmul_u16(uint16_t op1, uint16_t op2) {
    uint16_t res;
    __asm__ __volatile__ ( "gfmul %[rd], %[rs], %[rt]"
    : [rd] "=r" (res)
    : [rs] "r" (op1), [rt] "r" (op2)
    );
    return res;
}
#endif  // CUSTOM_INST_API_H