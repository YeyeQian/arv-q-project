#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "api.h"
#include "compiler.h"
#include "fips202.h"
#include "symmetric.h"
#include "param.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>
#include "../common/debug.h"

bool test_shake256(){
    size_t mu_sig_len=CRHBYTES+K*POLYW1_PACKEDBYTES;
    uint8_t* mu_sig=(uint8_t*)malloc(mu_sig_len*sizeof(uint8_t));
    for(int i=0;i<mu_sig_len;i++)mu_sig[i]=rand()%256;
    uint8_t sig[SEEDBYTES];
    shake256_custom(sig,SEEDBYTES,mu_sig,mu_sig_len);
    free(mu_sig);
}

int main(){
    test_shake256();
    return 0;
}