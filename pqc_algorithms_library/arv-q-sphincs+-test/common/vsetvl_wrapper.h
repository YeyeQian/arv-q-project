#ifndef VSETVL_WRAPPER_H
#define VSETVL_WRAPPER_H

#include <stddef.h>

static inline size_t vsetvl_e8m1_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e8, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e8m2_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e8, m2, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e8m4_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e8, m4, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e8m8_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e8, m8, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e16m1_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e16, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e16m2_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e16, m2, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e16m4_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e16, m4, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e16m8_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e16, m8, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e32m1_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e32, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e32m2_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e32, m2, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e32m4_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e32, m4, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e32m8_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e32, m8, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e64m1_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e64, m1, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e64m2_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e64, m2, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e64m4_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e64, m4, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}


static inline size_t vsetvl_e64m8_wrapper(size_t avl) {
    size_t vl;
    __asm__ __volatile__ ( "vsetvli %[vl], %[n], e64, m8, tu, mu" : [vl] "=r"(vl) : [n] "r"(avl) );
    return vl;
}
#endif
