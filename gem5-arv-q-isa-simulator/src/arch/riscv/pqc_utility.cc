#include "pqc_utility.hh"

namespace gem5
{

namespace RiscvISA
{

uint16_t degree(uint16_t a) {
	uint16_t r = 0;
	while (a >>= 1) {
		r++;
	}
	return r;
}

std::pair<uint16_t, uint16_t> divide(uint16_t a, uint16_t b) {

	uint16_t quot = 0; // quotient
	uint16_t rem = a; // remainder 

	if (b > a) { // then a is not divisible by b, and quot = 0
		return std::pair<uint16_t, uint16_t>(quot, rem);
	}

	while ((rem >= b) && (rem != 0)) {
		unsigned exp = degree(rem) - degree(b);
		quot |= (1 << exp);
		uint16_t btmp = (b << exp);
		rem ^= btmp;
	}
	return std::pair<uint16_t, uint16_t>(quot, rem);
};

uint16_t gf_add(uint16_t a, uint16_t b) {
	return a ^ b;
}

uint16_t gf_mul(uint16_t a, uint16_t b, uint16_t prim) {
	uint16_t res = 0;
	uint16_t m = degree(prim);
	while (a && b) {
		if (b & 1) res ^= a;
		if (a & (((uint16_t)1) << (m - 1))) // GF modulo: if a >= 2^m, it overflows when shifted left, so reduce 
			a = (a << 1) ^ prim; // XOR with the primitive polynomial
		else
			a <<= 1; // a*2 /* equivalent to a*2 
		b >>= 1; // b / 2	
	}
	return res;
}

uint16_t gf_inv(uint16_t a, uint16_t prim) {
	uint16_t rem1 = prim;
	uint16_t rem2 = a;
	uint16_t aux1 = 0;
	uint16_t aux2 = 1;

	while (rem2 != 0) {
		std::pair< uint16_t, uint16_t > res = divide(rem1, rem2);

		uint16_t aux_new = gf_add(gf_mul(res.first, aux2, prim),aux1);

		// prepare for the next step
		rem1 = rem2;
		rem2 = res.second;
		aux1 = aux2;
		aux2 = aux_new;

	}
	//aux1.multiply_monomial(pow(rem1.poly[0],-1), 0);
	return aux1;
}

}
}