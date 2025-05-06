#include <stdint.h>
#include <string.h>
#include <stdlib.h>

#include "address.h"
#include "utils.h"
#include "params.h"
#include "hash.h"
#include "fips202.h"

/*
 * Computes PRF(key, addr), given a secret key of SPX_N bytes and an address
 */
void prf_addr_custom(unsigned char *out, const unsigned char *key,
              const uint32_t addr[8])
{
    unsigned char buf[SPX_N + SPX_ADDR_BYTES];

    memcpy(buf, key, SPX_N);
    memcpy(buf + SPX_N, addr, SPX_ADDR_BYTES);

    shake256_custom(out, SPX_N, buf, SPX_N + SPX_ADDR_BYTES);
}

/**
 * Computes the message-dependent randomness R, using a secret seed and an
 * optional randomization value as well as the message.
 */
void gen_message_random_custom(unsigned char *R, const unsigned char *sk_prf,
                        const unsigned char *optrand,
                        const unsigned char *m, unsigned long long mlen)
{
    // uint64_t s_inc[26];
    // shake256_inc_init(s_inc);
    // shake256_inc_absorb(s_inc, sk_prf, SPX_N);
    // shake256_inc_absorb(s_inc, optrand, SPX_N);
    // shake256_inc_absorb(s_inc, m, mlen);
    // shake256_inc_finalize(s_inc);
    // shake256_inc_squeeze(R, SPX_N, s_inc);

    unsigned char *sk_prf_optrand_m;
    sk_prf_optrand_m = (unsigned char *) malloc(SPX_N + SPX_N + mlen);
    memcpy(sk_prf_optrand_m, sk_prf, SPX_N);
    memcpy(sk_prf_optrand_m+SPX_N, optrand, SPX_N);
    memcpy(sk_prf_optrand_m+SPX_N+SPX_N, m, mlen);
    shake256_custom(R, SPX_N, sk_prf_optrand_m, SPX_N + SPX_N + mlen);
}

/**
 * Computes the message hash using R, the public key, and the message.
 * Outputs the message digest and the index of the leaf. The index is split in
 * the tree index and the leaf index, for convenient copying to an address.
 */
void hash_message_custom(unsigned char *digest, uint64_t *tree, uint32_t *leaf_idx,
                  const unsigned char *R, const unsigned char *pk,
                  const unsigned char *m, unsigned long long mlen)
{
#define SPX_TREE_BITS (SPX_TREE_HEIGHT * (SPX_D - 1))
#define SPX_TREE_BYTES ((SPX_TREE_BITS + 7) / 8)
#define SPX_LEAF_BITS SPX_TREE_HEIGHT
#define SPX_LEAF_BYTES ((SPX_LEAF_BITS + 7) / 8)
#define SPX_DGST_BYTES (SPX_FORS_MSG_BYTES + SPX_TREE_BYTES + SPX_LEAF_BYTES)

    unsigned char buf[SPX_DGST_BYTES];
    unsigned char *bufp = buf;
    // uint64_t s_inc[26];

    // shake256_inc_init(s_inc);
    // shake256_inc_absorb(s_inc, R, SPX_N);
    // shake256_inc_absorb(s_inc, pk, SPX_PK_BYTES);
    // shake256_inc_absorb(s_inc, m, mlen);
    // shake256_inc_finalize(s_inc);
    // shake256_inc_squeeze(buf, SPX_DGST_BYTES, s_inc);

    unsigned char *R_pk_m;
    R_pk_m = (unsigned char *) malloc(SPX_N + SPX_PK_BYTES + mlen);
    memcpy(R_pk_m, R, SPX_N);
    memcpy(R_pk_m+SPX_N, pk, SPX_PK_BYTES);
    memcpy(R_pk_m+SPX_N+SPX_PK_BYTES, m, mlen);
    shake256_custom(buf, SPX_DGST_BYTES, R_pk_m, SPX_N + SPX_PK_BYTES + mlen);

    memcpy(digest, bufp, SPX_FORS_MSG_BYTES);
    bufp += SPX_FORS_MSG_BYTES;

#if SPX_TREE_BITS > 64
    #error For given height and depth, 64 bits cannot represent all subtrees
#endif

    *tree = bytes_to_ull(bufp, SPX_TREE_BYTES);
    *tree &= (~(uint64_t)0) >> (64 - SPX_TREE_BITS);
    bufp += SPX_TREE_BYTES;

    *leaf_idx = bytes_to_ull(bufp, SPX_LEAF_BYTES);
    *leaf_idx &= (~(uint32_t)0) >> (32 - SPX_LEAF_BITS);
}
