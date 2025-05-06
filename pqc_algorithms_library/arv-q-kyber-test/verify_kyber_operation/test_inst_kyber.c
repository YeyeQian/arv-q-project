#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <stdlib.h>
#include <time.h>
#include "api.h"
#include "compiler.h"
#include "fips202.h"
#include "symmetric.h"
#include "poly.h"
#include "params.h"
#include "../apis/custom_inst_api.h"
#include <riscv_vector.h>


#define VMOD_ADD_VV_M1_AVL 32
int test_vmod_add_vv_m1()
{   
  size_t avl = VMOD_ADD_VV_M1_AVL;
  size_t vl;
  int16_t a_vec[VMOD_ADD_VV_M1_AVL];
  int16_t b_vec[VMOD_ADD_VV_M1_AVL];
  int16_t r_vec[VMOD_ADD_VV_M1_AVL];
  int16_t r_vec_ref[VMOD_ADD_VV_M1_AVL];

  csr_modulusq_rw(KYBER_Q);
  vint16m1_t va;
  vint16m1_t vb;
  vint16m1_t vr;
  while(avl > 0){
    vl = vsetvl_e16m1(avl);
    va = vle16_v_i16m1(a_ptr, vl);
    vb = vle16_v_i16m1(b_ptr, vl);
    vr = vmod_add_vv_i16m8(va, vb);
    vse16_v_i16m8(r_ptr, vr, vl);
    r_ptr += vl;
    a_ptr += vl;
    b_ptr += vl;
    avl -= vl;
  }
}


int main(){
  test_poly_tomont_custom();
  return 0;
}