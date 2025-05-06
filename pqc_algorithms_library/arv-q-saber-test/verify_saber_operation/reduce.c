#include "reduce.h"


unsigned int montgomery_reduce(unsigned int a,unsigned int b) {
	unsigned long long t0 = (unsigned long long)a * (unsigned long long)b;
	unsigned int t1 = (unsigned int)(t0 & 0xffffffff);
	unsigned int t2 = ((unsigned long long)t1 *(unsigned long long) SABER_Q_EXT_INV) & 0xffffffff;
	unsigned long long t3 = t0 + (unsigned long long)t2 * (unsigned long long)SABER_Q_EXT;
	unsigned int res = (unsigned int)((unsigned long long)t3 >> 32);
	return res>=SABER_Q_EXT?res-SABER_Q_EXT:res;
}

unsigned int mod_add(unsigned int a, unsigned int b) {//a+b mod q, input>=0, output>=0
	unsigned int tmp = a + b;
	return tmp >= SABER_Q_EXT ? tmp - SABER_Q_EXT : tmp;
}

unsigned int mod_sub(unsigned int a, unsigned int b) {//a-b mod q, input>=0, output>=0
	int tmp = a - b;
	return tmp < 0 ? tmp + SABER_Q_EXT : tmp;
}