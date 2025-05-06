#include "gf2x_port.h"

// The k-squaring function computes c = a^(2^k) % (x^r - 1),
// By [1](Observation 1), if
//     a = sum_{j in supp(a)} x^j,
// then
//     a^(2^k) % (x^r - 1) = sum_{j in supp(a)} x^((j * 2^k) % r).
// Therefore, k-squaring can be computed as permutation of the bits of "a":
//     pi0 : j --> (j * 2^k) % r.
// For improved performance, we compute the result by inverted permutation pi1:
//     pi1 : (j * 2^-k) % r --> j.
// Input argument l_param is defined as the value (2^-k) % r.
void k_sqr_port(OUT pad_r_t* c, IN const pad_r_t* a, IN const size_t l_param)
{
    memset(c->val.raw, 0, sizeof(c->val));

    // Compute the result byte by byte
    size_t idx = 0;
    for (size_t i = 0; i < R_SIZE; i++) {
        for (size_t j = 0; j < 8; j++, idx++) {
            // Bit of "c" at position idx is set to the value of
            // the bit of "a" at position pi1(idx) = (l_param * idx) % R_BITS.
            size_t pos = (l_param * idx) % R_BITS;

            size_t  pos_byte = pos >> 3;
            size_t  pos_bit = pos & 7;
            uint8_t bit = (a->val.raw[pos_byte] >> pos_bit) & 1;

            c->val.raw[i] |= (bit << j);
        }
    }
    c->val.raw[R_SIZE - 1] &= LAST_R_BYTE_MASK;
}

#define LSB3(x) ((x)&7)

// 64x64 bit multiplication
// The algorithm is based on the windowing method, for example as in:
// Brent, R. P., Gaudry, P., Thomé, E., & Zimmermann, P. (2008, May), "Faster
// multiplication in GF (2)[x]". In: International Algorithmic Number Theory
// Symposium (pp. 153-166). Springer, Berlin, Heidelberg. In this implementation,
// the last three bits are multiplied using a schoolbook multiplication.
void gf2x_mul_base_port(OUT uint64_t* c,
    IN const uint64_t* a,
    IN const uint64_t* b)
{
    uint64_t       h = 0, l = 0, g1, g2, u[8];
    const uint64_t w = 64;
    const uint64_t s = 3;
    const uint64_t a0 = a[0];
    const uint64_t b0 = b[0];

    // Multiplying 64 bits by 7 can results in an overflow of 3 bits.
    // Therefore, these bits are masked out, and are treated in step 3.
    const uint64_t b0m = b0 & MASK(61);

    // Step 1: Calculate a multiplication table with 8 entries.
    u[0] = 0;
    u[1] = b0m;
    u[2] = u[1] << 1;
    u[3] = u[2] ^ b0m;
    u[4] = u[2] << 1;
    u[5] = u[4] ^ b0m;
    u[6] = u[3] << 1;
    u[7] = u[6] ^ b0m;

    // Step 2: Multiply two elements in parallel in positions i, i+s
    l = u[LSB3(a0)] ^ (u[LSB3(a0 >> 3)] << 3);
    h = (u[LSB3(a0 >> 3)] >> 61);

    for (size_t i = (2 * s); i < w; i += (2 * s)) {
        const size_t i2 = (i + s);

        g1 = u[LSB3(a0 >> i)];
        g2 = u[LSB3(a0 >> i2)];

        l ^= (g1 << i) ^ (g2 << i2);
        h ^= (g1 >> (w - i)) ^ (g2 >> (w - i2));
    }

    // Step 3: Multiply the last three bits.
    for (size_t i = 61; i < 64; i++) {
        uint64_t mask = (-((b0 >> i) & 1));
        l ^= ((a0 << i) & mask);
        h ^= ((a0 >> (w - i)) & mask);
    }

    c[0] = l;
    c[1] = h;
}

// c = a^2
void gf2x_sqr_port(OUT dbl_pad_r_t* c, IN const pad_r_t* a)
{
    const uint64_t* a64 = (const uint64_t*)a;
    uint64_t* c64 = (uint64_t*)c;

    for (size_t i = 0; i < R_QWORDS; i++) {
        gf2x_mul_base_port(&c64[2 * i], &a64[i], &a64[i]);
    }
}


#  define REG_T uint64_t
#define REG_QWORDS (sizeof(REG_T) / sizeof(uint64_t))
#  define LOAD(mem)       (mem)[0]
#  define STORE(mem, val) (mem)[0] = val

#  define SLLI_I64(a, imm) ((a) << (imm))
#  define SRLI_I64(a, imm) ((a) >> (imm))

// c = a mod (x^r - 1)
void gf2x_red_port(OUT pad_r_t* c, IN const dbl_pad_r_t* a)
{
    const uint64_t* a64 = (const uint64_t*)a;
    uint64_t* c64 = (uint64_t*)c;

    for (size_t i = 0; i < R_QWORDS; i += REG_QWORDS) {
        REG_T vt0 = LOAD(&a64[i]);
        REG_T vt1 = LOAD(&a64[i + R_QWORDS]);
        REG_T vt2 = LOAD(&a64[i + R_QWORDS - 1]);

        vt1 = SLLI_I64(vt1, LAST_R_QWORD_TRAIL);
        vt2 = SRLI_I64(vt2, LAST_R_QWORD_LEAD);

        vt0 ^= (vt1 | vt2);

        STORE(&c64[i], vt0);
    }

    c64[R_QWORDS - 1] &= LAST_R_QWORD_MASK;

    // Clean the secrets from the upper part of c
    memset((uint8_t*)&c64[R_QWORDS],0,(R_PADDED_QWORDS - R_QWORDS) * sizeof(uint64_t));
}

// a = a^2 mod (x^r - 1)
_INLINE_ void gf2x_mod_sqr_in_place(IN OUT pad_r_t* a,
    OUT dbl_pad_r_t* secure_buffer)
{
    gf2x_sqr_port(secure_buffer, a);
    gf2x_red_port(a, secure_buffer);
}

// c = a^2^2^num_sqrs
_INLINE_ void repeated_squaring(OUT pad_r_t* c,
    IN pad_r_t* a,
    IN const size_t num_sqrs,
    OUT dbl_pad_r_t* sec_buf)
{
    c->val = a->val;

    for (size_t i = 0; i < num_sqrs; i++) {
        gf2x_mod_sqr_in_place(c, sec_buf);
    }
}

// The secure buffer size required for Karatsuba is computed by:
//    size(n) = 3*n/2 + size(n/2) = 3*sum_{i}{n/2^i} < 3n
#define SECURE_BUFFER_QWORDS (3 * R_PADDED_QWORDS)

void karatzuba_add1_port(OUT uint64_t* alah,
    OUT uint64_t* blbh,
    IN const uint64_t* a,
    IN const uint64_t* b,
    IN const size_t    qwords_len)
{
    assert(qwords_len % REG_QWORDS == 0);

    REG_T va0, va1, vb0, vb1;

    for (size_t i = 0; i < qwords_len; i += REG_QWORDS) {
        va0 = LOAD(&a[i]);
        va1 = LOAD(&a[i + qwords_len]);
        vb0 = LOAD(&b[i]);
        vb1 = LOAD(&b[i + qwords_len]);

        STORE(&alah[i], va0 ^ va1);
        STORE(&blbh[i], vb0 ^ vb1);
    }
}

void karatzuba_add2_port(OUT uint64_t* z,
    IN const uint64_t* x,
    IN const uint64_t* y,
    IN const size_t    qwords_len)
{
    assert(qwords_len % REG_QWORDS == 0);

    REG_T vx, vy;

    for (size_t i = 0; i < qwords_len; i += REG_QWORDS) {
        vx = LOAD(&x[i]);
        vy = LOAD(&y[i]);

        STORE(&z[i], vx ^ vy);
    }
}

void karatzuba_add3_port(OUT uint64_t* c,
    IN const uint64_t* mid,
    IN const size_t    qwords_len)
{
    assert(qwords_len % REG_QWORDS == 0);

    REG_T vr0, vr1, vr2, vr3, vt;

    uint64_t* c0 = c;
    uint64_t* c1 = &c[qwords_len];
    uint64_t* c2 = &c[2 * qwords_len];
    uint64_t* c3 = &c[3 * qwords_len];

    for (size_t i = 0; i < qwords_len; i += REG_QWORDS) {
        vr0 = LOAD(&c0[i]);
        vr1 = LOAD(&c1[i]);
        vr2 = LOAD(&c2[i]);
        vr3 = LOAD(&c3[i]);
        vt = LOAD(&mid[i]);

        STORE(&c1[i], vt ^ vr0 ^ vr1);
        STORE(&c2[i], vt ^ vr2 ^ vr3);
    }
}

void karatzuba(OUT uint64_t* c,
    IN const uint64_t* a,
    IN const uint64_t* b,
    IN const size_t    qwords_len,
    IN const size_t    qwords_len_pad,
    uint64_t* sec_buf)
{
    if (qwords_len <= 1) {
        gf2x_mul_base_port(c, a, b);
        return;
    }

    const size_t half_qw_len = qwords_len_pad >> 1;

    // Split a and b into low and high parts of size n_padded/2
    const uint64_t* a_lo = a;
    const uint64_t* b_lo = b;
    const uint64_t* a_hi = &a[half_qw_len];
    const uint64_t* b_hi = &b[half_qw_len];

    // Split c into 4 parts of size n_padded/2 (the last ptr is not needed)
    uint64_t* c0 = c;
    uint64_t* c1 = &c[half_qw_len];
    uint64_t* c2 = &c[half_qw_len * 2];

    // Allocate 3 ptrs of size n_padded/2  on sec_buf
    uint64_t* alah = sec_buf;
    uint64_t* blbh = &sec_buf[half_qw_len];
    uint64_t* tmp = &sec_buf[half_qw_len * 2];

    // Move sec_buf ptr to the first free location for the next recursion call
    sec_buf = &sec_buf[half_qw_len * 3];

    // Compute a_lo*b_lo and store the result in (c1|c0)
    karatzuba(c0, a_lo, b_lo, half_qw_len, half_qw_len, sec_buf);

    // If the real number of digits n is less or equal to n_padded/2 then:
    //     a_hi = 0 and b_hi = 0
    // and
    //     (a_hi|a_lo)*(b_hi|b_lo) = a_lo*b_lo
    // so we can skip the remaining two multiplications
    if (qwords_len > half_qw_len) {
        // Compute a_hi*b_hi and store the result in (c3|c2)
        karatzuba(c2, a_hi, b_hi, qwords_len - half_qw_len, half_qw_len, sec_buf);

        // Compute alah = (a_lo + a_hi) and blbh = (b_lo + b_hi)
        karatzuba_add1_port(alah, blbh, a, b, half_qw_len);

        // Compute (c1 + c2) and store the result in tmp
        karatzuba_add2_port(tmp, c1, c2, half_qw_len);

        // Compute alah*blbh and store the result in (c2|c1)
        karatzuba(c1, alah, blbh, half_qw_len, half_qw_len, sec_buf);

        // Add (tmp|tmp) and (c3|c0) to (c2|c1)
        karatzuba_add3_port(c0, tmp, half_qw_len);
    }
}

void gf2x_mod_mul_with_ctx(OUT pad_r_t* c,
    IN const pad_r_t* a,
    IN const pad_r_t* b)
{
    static_assert((R_PADDED_BYTES % 2 == 0), "Karatsuba is odd");

    dbl_pad_r_t t = { 0 };
    uint64_t secure_buffer[SECURE_BUFFER_QWORDS];

    karatzuba((uint64_t*)&t, (const uint64_t*)a, (const uint64_t*)b, R_QWORDS,
        R_PADDED_QWORDS, secure_buffer);

    gf2x_red_port(c, &t);

    memset((uint8_t*)secure_buffer,0, sizeof(secure_buffer));
}

// Inversion in F_2[x]/(x^R - 1), [1](Algorithm 2).
// c = a^{-1} mod x^r-1
void gf2x_mod_inv(OUT pad_r_t* c, IN const pad_r_t* a)
{

    // Note that exp0/1_k/l are predefined constants that depend only on the value
    // of R. This value is public. Therefore, branches in this function, which
    // depends on R, are also "public". Code that releases these branches
    // (taken/not-taken) does not leak secret information.
    const size_t exp0_k[MAX_I] = { EXP0_K_VALS };
    const size_t exp0_l[MAX_I] = { EXP0_L_VALS };
    const size_t exp1_k[MAX_I] = { EXP1_K_VALS };
    const size_t exp1_l[MAX_I] = { EXP1_L_VALS };

    pad_r_t f = { 0 };
    pad_r_t g = { 0 };
    pad_r_t t = { 0 };
    dbl_pad_r_t sec_buf = { 0 };

    // Steps 2 and 3 in [1](Algorithm 2)
    f.val = a->val;
    t.val = a->val;

    for (size_t i = 1; i < MAX_I; i++) {
        // Step 5 in [1](Algorithm 2), exponentiation 0: g = f^2^2^(i-1)
        if (exp0_k[i - 1] <= K_SQR_THR) {
            repeated_squaring(&g, &f, exp0_k[i - 1], &sec_buf);
        }
        else {
            k_sqr_port(&g, &f, exp0_l[i - 1]);
        }

        // Step 6, [1](Algorithm 2): f = f*g
        gf2x_mod_mul_with_ctx(&f, &g, &f);

        if (exp1_k[i] != 0) {
            // Step 8, [1](Algorithm 2), exponentiation 1: g = f^2^((r-2) % 2^i)
            if (exp1_k[i] <= K_SQR_THR) {
                repeated_squaring(&g, &f, exp1_k[i], &sec_buf);
            }
            else {
                k_sqr_port(&g, &f, exp1_l[i]);
            }

            // Step 9, [1](Algorithm 2): t = t*g;
            gf2x_mod_mul_with_ctx(&t, &g, &t);
        }
    }

    // Step 10, [1](Algorithm 2): c = t^2
    gf2x_mod_sqr_in_place(&t, &sec_buf);
    c->val = t.val;
}

void gf2x_mod_inv_wrapper(OUT uint8_t res_bin[R_SIZE], IN const uint8_t a_bin[R_SIZE]) {
    pad_r_t a_bin_padded = { 0 };
    memcpy(a_bin_padded.val.raw,a_bin,R_SIZE);
    pad_r_t res_bin_padded = { 0 };
    gf2x_mod_inv(&res_bin_padded,&a_bin_padded);
    memcpy(res_bin,res_bin_padded.val.raw,R_SIZE);
}


/*****************************************
 * 
 *    Karatsuba Multiplication Related
 *    As Godlen Models in Custom Test
 * 
*****************************************/
void mul2_512(OUT uint64_t* h, OUT uint64_t* l, IN const uint64_t* a, IN const uint64_t* b){
    uint64_t s1[8]={0};
    uint64_t s2[8]={0};
    for(int i=0;i<4;i++){
        s1[i*2]=s1[i*2+1]=a[2*i]^a[2*i+1];
        s2[i*2]=s2[i*2+1]=b[2*i]^b[2*i+1];
    }

    uint64_t lq[8]={0};
    for(int i=0;i<4;i++){
        gf2x_mul_base_port(&lq[2*i],&a[2*i],&b[2*i]);
    }
    uint64_t hq[8]={0};
    for(int i=0;i<4;i++){
        gf2x_mul_base_port(&hq[2*i],&a[2*i+1],&b[2*i+1]);
    }

    uint64_t abq[8]={0};
    for(int i=0;i<4;i++){
        gf2x_mul_base_port(&abq[2*i],&s1[2*i],&s2[2*i]);
    }
    for(int i=0;i<8;i++){
        abq[i]^=lq[i]^hq[i];
    }

    for(int i=0;i<4;i++){
        uint64_t temp=abq[2*i];
        abq[2*i]=abq[2*i+1];
        abq[2*i+1]=temp;
    }

    for(int i=0;i<8;i++){
        if((i&1)==0){
            l[i]=lq[i];
        }else{
            l[i]=lq[i]^abq[i];
        }
    }

    for(int i=0;i<8;i++){
        if((i&1)==0){
            h[i]=hq[i]^abq[i];
        }else{
            h[i]=hq[i];
        }
    }
}

void gf2x_mul8_512_int(OUT uint64_t* h, OUT uint64_t* l, IN const uint64_t* a, IN const uint64_t* b){
    const uint64_t mask0[8]   = {0,1,8,9,4,5,12,13};
    const uint64_t mask1[8]   = {2,3,10,11,6,7,14,15};
    const uint64_t mask2[8]   = {4,5,6,7,0,1,2,3};
    const uint64_t mask3[8]   = {0,1,2,3,8,9,10,11};
    const uint64_t mask4[8]   = {4,5,6,7,12,13,14,15};
    const uint64_t mask_s1[8] = {2,3,0,1,4,5,6,7};
    const uint64_t mask_s2[8] = {0,1,4,5,6,7,2,3};

    uint64_t xl[8], xh[8], xabl[8], xabh[8], xab[8], xab1[8], xab2[8], oxh[8], oxl[8];
    uint64_t yl[8], yh[8], yabl[8], yabh[8], yab[8];
    uint64_t t[4][8];

    // t[0] = PERMXVAR_I64(mask_s1, a) ^ PERMXVAR_I64(mask_s2, a);
    // t[1] = PERMXVAR_I64(mask_s1, b) ^ PERMXVAR_I64(mask_s2, b);
    for(int i=0;i<8;i++){
        t[0][i]=a[mask_s1[i]]^a[mask_s2[i]];
        t[1][i]=b[mask_s1[i]]^b[mask_s2[i]];
    }

    // t[2] = t[0] ^ VALIGN(t[0], t[0], 4);
    // t[3] = t[1] ^ VALIGN(t[1], t[1], 4);
    uint64_t temp1=a[0]^a[2]^a[4]^a[6];
    uint64_t temp2=a[1]^a[3]^a[5]^a[7];
    for(int i=0;i<4;i++){
        t[2][2*i]=temp1;
        t[2][2*i+1]=temp2;
    }
    temp1=b[0]^b[2]^b[4]^b[6];
    temp2=b[1]^b[3]^b[5]^b[7];
    for(int i=0;i<4;i++){
        t[3][2*i]=temp1;
        t[3][2*i+1]=temp2;
    }

    // mul2_512(&xh, &xl, a, b);
    // mul2_512(&xabh, &xabl, t[0], t[1]);
    // mul2_512(&yabh, &yabl, t[2], t[3]);
    mul2_512(xh,xl,a,b);
    mul2_512(xabh,xabl,t[0],t[1]);
    mul2_512(yabh,yabl,t[2],t[3]);

    // xab  = xl ^ xh ^ PERMX2VAR_I64(xabl, mask0, xabh);
    for(int i=0;i<8;i++){
        xab[i]=xl[i]^xh[i]^(mask0[i]>=8?xabh[mask0[i]-8]:xabl[mask0[i]]);
    }

    // yl   = PERMX2VAR_I64(xl, mask3, xh);
    for(int i=0;i<8;i++){
        yl[i]=(mask3[i]>=8?xh[mask3[i]-8]:xl[mask3[i]]);
    }

    //yh   = PERMX2VAR_I64(xl, mask4, xh);
    for(int i=0;i<8;i++){
        yh[i]=(mask4[i]>=8?xh[mask4[i]-8]:xl[mask4[i]]);
    }

    //xab1 = VALIGN(xab, xab, 6);
    for(int i=0;i<8;i++){
        xab1[i]=xab[(i+6)&7];
    }

    //xab2 = VALIGN(xab, xab, 2);
    for(int i=0;i<8;i++){
        xab2[i]=xab[(i+2)&7];
    }

    uint8_t mask=0x3c;
    //yl   = MXOR_I64(yl, 0x3c, yl, xab1);
    for(int i=0;i<8;i++){
        if((mask&1)==1){
            yl[i]^=xab1[i];
        }
        mask>>=1;
    }

    mask=0x3c;
    //yh   = MXOR_I64(yh, 0x3c, yh, xab2);
    for(int i=0;i<8;i++){
        if((mask&1)==1){
            yh[i]^=xab2[i];
        }
        mask>>=1;
    }

    //__m512i oxh = PERMX2VAR_I64(xabl, mask1, xabh);
    for(int i=0;i<8;i++){
        oxh[i]=(mask1[i]>=8?xabh[mask1[i]-8]:xabl[mask1[i]]);
    }

    //__m512i oxl = VALIGN(oxh, oxh, 4);
    for(int i=0;i<8;i++){
        oxl[i]=oxh[(i+4)&7];
    }

    //yab         = oxl ^ oxh ^ PERMX2VAR_I64(yabl, mask0, yabh);
    for(int i=0;i<8;i++){
        yab[i]=oxl[i]^oxh[i]^(mask0[i]>=8?yabh[mask0[i]-8]:yabl[mask0[i]]);
    }

    //yab         = MXOR_I64(oxh, 0x3c, oxh, VALIGN(yab, yab, 2));
    uint64_t yab_temp[8];
    for(int i=0;i<8;i++){
        yab_temp[i]=yab[(i+2)&7];
    }
    mask=0x3c;
    for(int i=0;i<8;i++){
        if((mask&1)==1){
            yab[i]=oxh[i]^yab_temp[i];
        }else{
            yab[i]=oxh[i];
        }
        mask>>=1;
    }

    //yab ^= yl ^ yh;
    for(int i=0;i<8;i++){
        yab[i]^=yl[i]^yh[i];
    }

    //yab = PERMXVAR_I64(mask2, yab);
    for(int i=0;i<8;i++){
        yab_temp[i]=yab[mask2[i]];
    }

    //*zl = MXOR_I64(yl, 0xf0, yl, yab);
    mask=0xf0;
    for(int i=0;i<8;i++){
        if((mask&1)==1){
            l[i]=yl[i]^yab_temp[i];
        }else{
            l[i]=yl[i];
        }
        mask>>=1;
    }

    //*zh = MXOR_I64(yh, 0x0f, yh, yab);
    mask=0x0f;
    for(int i=0;i<8;i++){
        if((mask&1)==1){
            h[i]=yh[i]^yab_temp[i];
        }else{
            h[i]=yh[i];
        }
        mask>>=1;
    }
}


// 1024x1024 bit multiplication performed by Karatsuba algorithm.
// Here, a and b are considered as having 16 digits of size 64 bits.
void gf2x_mul_base_vpclmul(OUT uint64_t *c,IN const uint64_t *a,IN const uint64_t *b){
    uint64_t hi[2][8], lo[2][8], mi[2][8];
    uint64_t a_tmp[8],b_tmp[8];
    for(int i=0;i<8;i++){
        a_tmp[i]=a[i]^a[i+8];
        b_tmp[i]=b[i]^b[i+8];
    }

    gf2x_mul8_512_int(lo[1], lo[0], a, b);
    gf2x_mul8_512_int(hi[1], hi[0], &a[8], &b[8]);
    gf2x_mul8_512_int(mi[1], mi[0], a_tmp, b_tmp);

    uint64_t m[8];
    for(int i=0;i<8;i++){
        m[i]=lo[1][i]^hi[0][i];
    }

    //STORE(&c[0 * QWORDS_IN_ZMM], lo[0]);
    for(int i=0;i<8;i++){
        c[i]=lo[0][i];
    }

    //STORE(&c[1 * QWORDS_IN_ZMM], mi[0] ^ lo[0] ^ m);
    for(int i=0;i<8;i++){
        c[i+8]=mi[0][i]^lo[0][i]^m[i];
    }

    //STORE(&c[2 * QWORDS_IN_ZMM], mi[1] ^ hi[1] ^ m);
    for(int i=0;i<8;i++){
        c[i+16]=mi[1][i]^hi[1][i]^m[i];
    }

    //STORE(&c[3 * QWORDS_IN_ZMM], hi[1]);
    for(int i=0;i<8;i++){
        c[i+24]=hi[1][i];
    }
}