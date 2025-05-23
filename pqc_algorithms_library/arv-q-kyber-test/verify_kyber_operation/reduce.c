#include <stdint.h>
#include "params.h"
#include "reduce.h"

/*************************************************
* Name:        montgomery_reduce
*
* Description: Montgomery reduction; given a 32-bit integer a, computes
*              16-bit integer congruent to a * R^-1 mod q,
*              where R=2^16
*
* Arguments:   - int32_t a: input integer to be reduced;
*                           has to be in {-q2^15,...,q2^15-1}
*
* Returns:     integer in {-q+1,...,q-1} congruent to a * R^-1 modulo q.
**************************************************/
int16_t montgomery_reduce(int32_t a)
{
  int32_t t;
  int16_t u;

  u = a*QINV;
  t = (int32_t)u*KYBER_Q;
  t = a - t;
  t >>= 16;
  return t;
}

/*************************************************
* Name:        barrett_reduce
*
* Description: Barrett reduction; given a 16-bit integer a, computes
*              16-bit integer congruent to a mod q in {0,...,q}
*
* Arguments:   - int16_t a: input integer to be reduced
*
* Returns:     integer in {0,...,q} congruent to a modulo q.
**************************************************/
int16_t barrett_reduce(int16_t a) {
  int16_t t;
  const int16_t v = ((1U << 26) + KYBER_Q/2)/KYBER_Q;

  t  = (int32_t)v*a >> 26;
  t *= KYBER_Q;
  return a - t;
}

/*************************************************
* Name:        csubq
*
* Description: Conditionallly subtract q
*
* Arguments:   - int16_t x: input integer
*
* Returns:     a - q if a >= q, else a
**************************************************/
int16_t csubq(int16_t a) {
  a -= KYBER_Q;
  a += (a >> 15) & KYBER_Q;
  return a;
}

/*************************************************
* Name:        mod_div2
*
* Description: modular divide 2
*
* Arguments:   - int16_t x: input integer
*
* Returns:     if(x & 1== 1) (x>>1 + (Q+1)/2) else (x>>1) 
**************************************************/
int16_t mod_div2(int16_t a) {
	int16_t t;

	t = (a >> 1) + (a & 1) * (KYBER_Q + 1) / 2;

	return t;
}

/*************************************************
* Name:        central_reduce_oneQ
*
* Description:  For finite field element a,
*              compute r \equiv a (mod Q) such that -Q < r < Q.
*
* Arguments:   - int16_t: finite field element a
*
* Returns r.
**************************************************/
int16_t central_reduce_oneQ(int16_t a) {
	int16_t r = a;
	if (a >= KYBER_Q) {
		r = a - KYBER_Q;
	}
	else if (a <= -KYBER_Q) {
		r = a + KYBER_Q;
	}

	return r;
}