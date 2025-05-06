#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <riscv_vector.h>
#include "defs.h"
#include "types.h"
#include "conversions.h"

#define BIT(len)       (1ULL << (len))
#define MASK(len)      (BIT(len) - 1)
#define LAST_R_BYTE_LEAD  (R_BITS & MASK(3))
#define LAST_R_BYTE_TRAIL (8 - LAST_R_BYTE_LEAD)
#define LAST_R_BYTE_MASK  MASK(LAST_R_BYTE_LEAD)

#define LAST_R_QWORD_LEAD  (R_BITS & MASK(6))
#define LAST_R_QWORD_TRAIL (64 - LAST_R_QWORD_LEAD)
#define LAST_R_QWORD_MASK  MASK(LAST_R_QWORD_LEAD)

#define K_SQR_THR (64)

#ifdef PARAM64
static_assert(R_BITS==12323,"BIKE LEVEL1 R_BITS SHOULD BE 12323");

// MAX_I = floor(log(r-2)) + 1
#  define MAX_I (14)
#  define EXP0_K_VALS \
     1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192
#  define EXP0_L_VALS \
     6162, 3081, 3851, 5632, 22, 484, 119, 1838, 1742, 3106, 10650, 1608, 10157, 8816
#  define EXP1_K_VALS 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 33, 4129
#  define EXP1_L_VALS 0, 0, 0, 0, 0, 6162, 0, 0, 0, 0, 0, 0, 242, 5717

#elif defined PARAM96
// The parameters below are hard-coded for R=24659
static_assert(R_BITS == 24659, "BIKE LEVEL3 R_BITS SHOULD BE 24659");

// MAX_I = floor(log(r-2)) + 1
#  define MAX_I (15)
#  define EXP0_K_VALS \
     1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384
#  define EXP0_L_VALS                                                           \
     12330, 6165, 7706, 3564, 2711, 1139, 15053, 1258, 4388, 20524, 9538, 6393, \
     10486, 1715, 6804
#  define EXP1_K_VALS 0, 0, 0, 0, 1, 0, 17, 0, 0, 0, 0, 0, 0, 81, 8273
#  define EXP1_L_VALS 0, 0, 0, 0, 12330, 0, 13685, 0, 0, 0, 0, 0, 0, 23678, 19056

#else
// The parameters below are hard-coded for R=40973
static_assert(R_BITS == 40973, "BIKE LEVEL5 R_BITS SHOULD BE 40973");

// MAX_I = floor(log(r-2)) + 1
#  define MAX_I (16)
#  define EXP0_K_VALS \
     1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768
#  define EXP0_L_VALS                                                         \
     20487, 30730, 28169, 9443, 13001, 12376, 8302, 6618, 38760, 21582, 1660, \
     10409, 14669, 30338, 17745, 7520
#  define EXP1_K_VALS 0, 1, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 11, 0, 8203
#  define EXP1_L_VALS 0, 20487, 0, 15365, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6302, 0, 10058

#endif

////////////////////////////////////
//      Portable Version
////////////////////////////////////
void gf2x_sqr_port(OUT dbl_pad_r_t* c, IN const pad_r_t* a);
void gf2x_red_port(OUT pad_r_t* c, IN const dbl_pad_r_t* a);

void gf2x_mul_base_port(OUT uint64_t* c,
    IN const uint64_t* a,
    IN const uint64_t* b);

void karatzuba_add1_port(OUT uint64_t* alah,
    OUT uint64_t* blbh,
    IN const uint64_t* a,
    IN const uint64_t* b,
    IN const size_t    qwords_len);

void karatzuba_add2_port(OUT uint64_t* z,
    IN const uint64_t* x,
    IN const uint64_t* y,
    IN const size_t    qwords_len);

void karatzuba_add3_port(OUT uint64_t* c,
    IN const uint64_t* mid,
    IN const size_t    qwords_len);

void karatzuba(OUT uint64_t* c,
    IN const uint64_t* a,
    IN const uint64_t* b,
    IN const size_t    qwords_len,
    IN const size_t    qwords_len_pad,
    uint64_t* sec_buf);

void gf2x_mod_mul_with_ctx(OUT pad_r_t* c,
    IN const pad_r_t* a,
    IN const pad_r_t* b);

void gf2x_mod_inv_wrapper(OUT uint8_t res_bin[R_SIZE], IN const uint8_t a_bin[R_SIZE]);

void mul2_512(OUT uint64_t* h, OUT uint64_t* l, IN const uint64_t* a, IN const uint64_t* b);

void gf2x_mul8_512_int(OUT uint64_t* h, OUT uint64_t* l, IN const uint64_t* a, IN const uint64_t* b);

void gf2x_mul_base_vpclmul(OUT uint64_t *c,IN const uint64_t *a,IN const uint64_t *b);

void k_sqr_port(OUT pad_r_t* c, IN const pad_r_t* a, IN const size_t l_param);

////////////////////////////////////
//      Customized Version
////////////////////////////////////
void gf2x_sqr_custom(OUT dbl_pad_r_t* c, IN const pad_r_t* a);
void gf2x_red_custom(OUT pad_r_t* c, IN const dbl_pad_r_t* a);

void gf2x_mod_mul_with_ctx_custom(OUT pad_r_t* c,
    IN const pad_r_t* a,
    IN const pad_r_t* b);

vuint64m2_t mul2_512_custom(IN vuint64m1_t va, IN vuint64m1_t vb);

vuint64m2_t gf2x_mul8_512_int_custom(IN vuint64m1_t va, IN vuint64m1_t vb);

void gf2x_mul_base_vpclmul_custom(OUT uint64_t *c,IN const uint64_t *a,IN const uint64_t *b);

void k_sqr_custom(OUT pad_r_t *c, IN const pad_r_t *a, IN const size_t l_param);

void gf2x_mod_inv_wrapper_custom(OUT uint8_t res_bin[R_SIZE], IN const uint8_t a_bin[R_SIZE]);