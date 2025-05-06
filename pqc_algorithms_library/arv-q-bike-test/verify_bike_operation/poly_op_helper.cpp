#include <stdint.h>
#include "poly_op_helper.h"
//some helper functions
int hamming_count(uint32_t x) {
	int count = 0;
	while (x) {
		if (x & 0x1)++count;
		x = (x >> 1);
	}
	return count;
}
//quick bitreverse for 32bit number
void reverse(uint32_t* x)
{
	*x = (((*x & 0xaaaaaaaa) >> 1) | ((*x & 0x55555555) << 1));
	*x = (((*x & 0xcccccccc) >> 2) | ((*x & 0x33333333) << 2));
	*x = (((*x & 0xf0f0f0f0) >> 4) | ((*x & 0x0f0f0f0f) << 4));
	*x = (((*x & 0xff00ff00) >> 8) | ((*x & 0x00ff00ff) << 8));

	*x = ((*x >> 16) | (*x << 16));
}

void reverse(uint8_t* x) {
	*x = (((*x & 0xaa) >> 1) | ((*x & 0x55) << 1));
	*x = (((*x & 0xcc) >> 2) | ((*x & 0x33) << 2));

	*x = ((*x >> 4) | (*x << 4));
}

void reverse(uint16_t* x) {
	*x = (((*x & 0xaaaa) >> 1) | ((*x & 0x5555) << 1));
	*x = (((*x & 0xcccc) >> 2) | ((*x & 0x3333) << 2));
	*x = (((*x & 0xf0f0) >> 4) | ((*x & 0x0f0f) << 4));

	*x = ((*x >> 8) | (*x << 8));
}

//Check_parity, check whether two array with element of uint8_t and length of R_SIZE is the same 
bool Check_parity(IN const uint8_t res[R_SIZE], IN const uint8_t ref[R_SIZE]) {
	bool wrong_flag = false;
	for (int i = 0; i < R_SIZE; i++) {
		if (res[i] != ref[i]) {
			wrong_flag = true;
			break;
		}
	}
	return wrong_flag;
}