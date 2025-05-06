#include "operations.h"

#include "controlbits.h"
#include "crypto_hash.h"
#include "encrypt.h"
#include "decrypt.h"
#include "params.h"
#include "sk_gen.h"
#include "pk_gen.h"
#include "util.h"

#include <stdint.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>

int crypto_kem_enc_custom(
       unsigned char *c,
       unsigned char *key,
       const unsigned char *pk
)
{
	unsigned char e[ SYS_N/8 ]={0};
	unsigned char one_ec[ 1 + SYS_N/8 + SYND_BYTES ] = {1};

	//

	encrypt_custom(c, pk, e);

	memcpy(one_ec + 1, e, SYS_N/8);
	memcpy(one_ec + 1 + SYS_N/8, c, SYND_BYTES);

	crypto_hash_32b_custom(key, one_ec, sizeof(one_ec));

	return 0;
}

int crypto_kem_dec_custom(
       unsigned char *key,
       const unsigned char *c,
       const unsigned char *sk
)
{
	int i;
    size_t vl,avl;
    vuint8m8_t va,vb;

	unsigned char ret_decrypt = 0;

	uint16_t m;

	unsigned char e[ SYS_N/8 ]={0};
	unsigned char preimage[ 1 + SYS_N/8 + SYND_BYTES ];
	unsigned char *x = preimage;
	const unsigned char *s = sk + 40 + IRR_BYTES + COND_BYTES;

	//

	ret_decrypt = decrypt_custom(e, sk + 40, c);

	m = ret_decrypt;
	m -= 1;
	m >>= 8;

	*x++ = m & 1;
    avl=SYS_N>>3;
    const uint8_t* src1_addr=s;
    uint8_t* src2_addr=e;
	while(avl>0){
        vl=vsetvl_e8m8(avl);
        va=vle8_v_u8m8(src1_addr,vl);
        vb=vle8_v_u8m8(src2_addr,vl);
        va=vand_vx_u8m8(va,~m,vl);
        vb=vand_vx_u8m8(vb,m,vl);
        va=vor_vv_u8m8(va,vb,vl);
        vse8_v_u8m8(x,va,vl);
        x+=vl,src1_addr+=vl,src2_addr+=vl,avl-=vl;
    }

	for (i = 0; i < SYND_BYTES; i++) 
		*x++ = c[i];

	crypto_hash_32b_custom(key, preimage, sizeof(preimage)); 

	return 0;
}

int crypto_kem_keypair_custom
(
       unsigned char *pk,
       unsigned char *sk 
)
{
	int i;
	unsigned char seed[ 33 ] = {64};
	unsigned char r[ SYS_N/8 + (1 << GFBITS)*sizeof(uint32_t) + SYS_T*2 + 32 ];
	unsigned char *rp, *skp;

	gf f[ SYS_T ]; // element in GF(2^mt)
	gf irr[ SYS_T ]; // Goppa polynomial
	uint32_t perm[ 1 << GFBITS ]; // random permutation as 32-bit integers
	int16_t pi[ 1 << GFBITS ]; // random permutation

	srand(time(NULL));
	for (int i = 0; i < 32; i++){
		seed[i + 1] = rand();
		//seed[i + 1] = i+1;//for test
	}

	while (1)
	{
		rp = &r[ sizeof(r)-32 ];
		skp = sk;

		// expanding and updating the seed

		shake_custom(r, sizeof(r), seed, 33);
		memcpy(skp, seed+1, 32);
		skp += 32 + 8;
		memcpy(seed+1, &r[ sizeof(r)-32 ], 32);

		// generating irreducible polynomial

		rp -= sizeof(f); 

		for (i = 0; i < SYS_T; i++) 
			f[i] = load_gf(rp + i*2); 

		if (genpoly_gen_custom(irr, f)) 
			continue;

		for (i = 0; i < SYS_T; i++)
			store_gf(skp + i*2, irr[i]);

		skp += IRR_BYTES;

		// generating permutation

		rp -= sizeof(perm);

		for (i = 0; i < (1 << GFBITS); i++) 
			perm[i] = load4(rp + i*4); 

		if (pk_gen_custom(pk, skp - IRR_BYTES, perm, pi))
			continue;

		controlbitsfrompermutation(skp, pi, GFBITS, 1 << GFBITS);
		skp += COND_BYTES;

		// storing the random string s

		rp -= SYS_N/8;
		memcpy(skp, rp, SYS_N/8);

		// storing positions of the 32 pivots

		store8(sk + 32, 0xFFFFFFFF);

		break;
	}

	return 0;
}
