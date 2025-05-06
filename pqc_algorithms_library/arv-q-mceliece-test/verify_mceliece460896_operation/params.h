#ifndef PARAMS_H
#define PARAMS_H

#define CRYPTO_NAMESPACE(x) x
#define _CRYPTO_NAMESPACE(x) _##x
//#define KAT
#define KATNUM 10

#define GFBITS 13
#define SYS_N 4608
#define SYS_T 96

#define COND_BYTES ((1 << (GFBITS-4))*(2*GFBITS - 1))
#define IRR_BYTES (SYS_T * 2)

#define PK_NROWS (SYS_T*GFBITS) 
#define PK_NCOLS (SYS_N - PK_NROWS)
#define PK_ROW_BYTES ((PK_NCOLS + 7)/8)

#define SYND_BYTES ((PK_NROWS + 7)/8)

#define GFMASK ((1 << GFBITS) - 1)

#define CEIL_DIVIDE(a, b)  (((a)/(b)) + ((a) % (b) == 0 ? 0 : 1))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define MC_GF_POLY 0x201B

#endif

