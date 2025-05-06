#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "api.h"
#include "compiler.h"
#include "cbd.h"
#include "params.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>


int test_cbd_eta1_custom()
{   
    uint8_t a[KYBER_ETA1*KYBER_N/4];
    poly r;
    bool flag;
    int i;

    srand((unsigned)time(NULL));
    for(i=0;i<(KYBER_ETA1*KYBER_N/4);i++){
        a[i]=rand()&255;
    }

    cbd_eta1_custom(&r,a);

    for(int i=0;i<KYBER_N;i++){
        printf("r[%d]=%d\n",i,r.coeffs[i]);
    }

    return 0;
}

int main(){
    test_cbd_eta1_custom();
    return 0;
}