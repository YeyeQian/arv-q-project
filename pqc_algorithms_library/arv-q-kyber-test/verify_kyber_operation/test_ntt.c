#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>

#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"

bool absolute_justify(int16_t a, int16_t b) {
    bool flag;

    int16_t diff = a - b;

    flag = ((diff % KYBER_Q) == 0);

    return flag;
}

#if (VLEN == 512)
bool test_ntt_compact()
{
    int16_t r[256];
    int16_t r_ref[256];
    int i;

    for(i = 0; i < 256; i++) {
        r[i] = i;
        r_ref[i] = i;
    }

    ntt_CG(r_ref);

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);
    ntt_cg_custom_zeta_compact(r);

    bool flag = true;
    for(i = 0; i < 256; i++) {
        if(r[i] != r_ref[i]) {
            if(!absolute_justify(r[i], r_ref[i])) {
                printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
                flag = false;
            }
        }
    }    

    if(flag) {
        printf("test_ntt_compact test pass!\n");
    }
    else {
        printf("test_ntt_compact test fail..\n");
    }

    return flag;
}

bool test_ntt_compact_asm()
{
    int16_t r[256];
    int16_t r_ref[256];
    int i;
    uint64_t start, end;

    for(i = 0; i < 256; i++) {
        r[i] = i;
        r_ref[i] = i;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    start = read_cycle();
    ntt_cg_custom_zeta_compact(r_ref);
    end = read_cycle();
    printf("ntt_cg_custom_zeta_compact takes %ld cycles\n", end-start);

    ntt_cg_custom_asm(r);

    bool flag = true;
    for(i = 0; i < 256; i++) {
        if(r[i] != r_ref[i]) {
            if(!absolute_justify(r[i], r_ref[i])) {
                printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
                flag = false;
            }
        }
    }    

    if(flag) {
        printf("test_ntt_compact_asm test pass!\n");
    }
    else {
        printf("test_ntt_compact_asm test fail..\n");
    }

    return flag;
}

bool test_intt_reorder_asm()
{
    int16_t r[256];
    int16_t r_ref[256];
    int i;
    bool flag;

    for(i = 0; i < 256; i++) {
        r[i] = i;
        r_ref[i] = i;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    ntt_cg_custom_reorder(r_ref);
    ntt_cg_custom_reorder_asm(r);

    flag = true;
    for(i = 0; i < 256; i++) {
        if(r[i] != r_ref[i]) {
            if(!absolute_justify(r[i], r_ref[i])) {
                printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
                flag = false;
            }

        }
    }    

    if(flag) {
        printf("test_ntt_reorder_asm test pass with flag=%d\n", flag);
    }
    else {
        printf("test_ntt_reorder_asm test fail with flag=%d\n", flag);
    }

    invntt_cg_custom_reorder(r_ref);

    invntt_cg_custom_reorder_asm(r);

    flag = true;
    for(i = 0; i < 256; i++) {
        if(r[i] != r_ref[i]) {
            if(!absolute_justify(r[i], r_ref[i])) {
                printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
                flag = false;
            }
        }
    }

    if(flag) {
        printf("test_intt_reorder_asm test pass with flag=%d\n", flag);
    }
    else {
        printf("test_intt_reorder_asm test fail with flag=%d\n", flag);
    }    

    return flag;
}
#endif

bool test_ntt_reorder()
{
    int16_t r[256];
    int16_t r_ref[256];
    int i;

    for(i = 0; i < 256; i++) {
        r[i] = i;
        r_ref[i] = i;
    }

    ntt_CG(r_ref);

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);
    ntt_cg_custom_reorder(r);

    bitreverse_cg_to_standard_all(r_ref);

    bool flag = true;
    for(i = 0; i < 256; i++) {
        if(r[i] != r_ref[i]) {
            if(!absolute_justify(r[i], r_ref[i])) {
                printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
                flag = false;
            }

        }
    }    

    if(flag) {
        printf("test_ntt_reorder test pass!\n");
    }
    else {
        printf("test_ntt_reorder test fail..\n");
    }

    return flag;    
}

bool test_ntt_reorder_asm()
{
    int16_t r[256];
    int16_t r_ref[256];
    int i;

    for(i = 0; i < 256; i++) {
        r[i] = i;
        r_ref[i] = i;
    }

    csr_modulusq_rw(KYBER_Q);
    csr_qinv_rw(QINV_HW);

    ntt_cg_custom_reorder(r_ref);

    ntt_cg_custom_reorder_asm(r);

    bool flag = true;
    for(i = 0; i < 256; i++) {
        if(r[i] != r_ref[i]) {
            if(!absolute_justify(r[i], r_ref[i])) {
                printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
                flag = false;
            }

        }
    }    

    if(flag) {
        printf("test_ntt_reorder_asm test pass!\n");
    }
    else {
        printf("test_ntt_reorder_asm test fail..\n");
    }

    return flag;    
}

// bool test_intt_compact()
// {
//     int16_t r[256];
//     int16_t r_ref[256];
//     int i;
//     bool flag;

//     for(i = 0; i < 256; i++) {
//         r[i] = i;
//         r_ref[i] = i;
//     }

//     ntt_CG(r_ref);
//     bitreverse_cg_to_standard_all(r_ref);

//     csr_modulusq_rw(KYBER_Q);
//     csr_qinv_rw(QINV_HW);
//     ntt_cg_custom_reorder(r);

//     flag = true;
//     for(i = 0; i < 256; i++) {
//         if(r[i] != r_ref[i]) {
//             if(!absolute_justify(r[i], r_ref[i])) {
//                 printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
//                 flag = false;
//             }

//         }
//     }    

//     if(flag) {
//         printf("test_ntt_reorder test pass!\n\n");
//     }
//     else {
//         printf("test_ntt_reorder test fail..\n\n");
//         return false;
//     }

//     invntt_CG(r_ref);
//     invntt_cg_zeta_compact(r);

//     flag = true;
//     for(i = 0; i < 256; i++) {
//         if(r[i] != r_ref[i]) {
//             if(!absolute_justify(r[i], r_ref[i])) {
//                 printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
//                 flag = false;
//             }

//         }
//     }

//     if(flag) {
//         printf("test_intt_compact test pass!\n");
//     }
//     else {
//         printf("test_intt_compact test fail..\n");
//     }    

//     return flag;
// }

// bool test_intt_reorder()
// {
//     int16_t r[256];
//     int16_t r_ref[256];
//     int i;
//     bool flag;

//     for(i = 0; i < 256; i++) {
//         r[i] = i;
//         r_ref[i] = i;
//     }

//     ntt_CG(r_ref);
//     bitreverse_cg_to_standard_all(r_ref);

//     csr_modulusq_rw(KYBER_Q);
//     csr_qinv_rw(QINV_HW);
//     ntt_cg_custom_reorder(r);

//     flag = true;
//     for(i = 0; i < 256; i++) {
//         if(r[i] != r_ref[i]) {
//             if(!absolute_justify(r[i], r_ref[i])) {
//                 printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
//                 flag = false;
//             }

//         }
//     }    

//     if(flag) {
//         printf("test_ntt_reorder test pass!\n\n");
//     }
//     else {
//         printf("test_ntt_reorder test fail..\n\n");
//         return false;
//     }

//     invntt_CG(r_ref);
//     bitreverse_cg_to_standard_all(r_ref);

//     for(i = 0; i < 256; i++) {
//         if(r[i] < 0) {
//             printf("r[%d] = %d\n", i, r[i]);
//         }
//     }

//     invntt_cg_custom_reorder(r);

//     flag = true;
//     for(i = 0; i < 256; i++) {
//         if(r[i] != r_ref[i]) {
//             if(!absolute_justify(r[i], r_ref[i])) {
//                 printf("r_ref[%d]=%d, while r[%d]=%d\n", i, r_ref[i], i, r[i]);
//                 flag = false;
//             }
//         }
//     }

//     if(flag) {
//         printf("test_intt_reorder test pass!\n");
//     }
//     else {
//         printf("test_intt_reorder test fail..\n");
//     }    

//     return flag;
// }

int main()
{   
    // test_ntt_compact_asm();
    // test_ntt_reorder();
    test_ntt_reorder_asm();
    // test_intt_reorder_asm();
    return 0;
}