#include <stdio.h>
#include "decrypt.h"

#include "params.h"
#include "benes.h"
#include "util.h"
#include "synd.h"
#include "root.h"
#include "gf.h"
#include "bm.h"
#include "../apis/custom_inst_api.h"
#include "../apis/vsetvl_wrapper.h"
#include <riscv_vector.h>
#include <string.h>

int decrypt_custom(unsigned char *e, const unsigned char *sk, const unsigned char *c)   //VLEN Must be greater equal than 256bit
{
	int i, w = 0; 
	uint16_t check;	
    size_t vl;
    vuint16m8_t va,vb;
    vuint16m1_t vc;

	unsigned char r[ SYS_N/8 ]={0};

	gf g[ SYS_T+1 ]={0};
	gf L[ SYS_N ]={0};

	gf s[ SYS_T*2 ]={0};
	gf s_cmp[ SYS_T*2 ]={0};
	gf locator[ SYS_T+1 ]={0};
	gf images[ SYS_N ]={0};

	gf t;

	//

	memcpy(r,c,SYND_BYTES);

	for (i = 0; i < SYS_T; i++) { g[i] = load_gf(sk); sk += 2; } g[ SYS_T ] = 1;

	support_gen(L, sk);

	synd_custom(s, g, L, r);

	bm_custom(locator, s);

	root_custom(images, locator, L);

	for (i = 0; i < SYS_N; i++)
	{
		t = gf_iszero(images[i]) & 1;

		e[ i/8 ] |= t << (i%8);
		w += t;

	}

#ifdef KAT
  {
    int k;
    printf("decrypt e: positions");
    for (k = 0;k < SYS_N;++k)
      if (e[k/8] & (1 << (k&7)))
        printf(" %d",k);
    printf("\n");
  }
#endif
	
	synd_custom(s_cmp, g, L, e);

	//

	check = w;
	check ^= SYS_T;

	vl=vsetvl_e16m1(1);
    vc=vmv_s_x_u16m1(vc,check,vl);
    vl=vsetvl_e16m8(SYS_T<<1);
    va=vle16_v_u16m8(s,vl);
    vb=vle16_v_u16m8(s_cmp,vl);
    va=vxor_vv_u16m8(va,vb,vl);
    vc=vredor_vs_u16m8_u16m1(vc,va,vc,vl);
    check=vmv_x_s_u16m1_u16(vc);

	check -= 1;
	check >>= 15;

	return check ^ 1;
}