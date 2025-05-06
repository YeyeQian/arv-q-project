#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include "inner.h"
#include "params_custom.h"

void copy_array(uint16_t a[FALCON_N], uint16_t temp[FALCON_N]) {
	uint32_t i;
	for (i = 0; i < FALCON_N; i++) {
		a[i] = temp[i];
	}
}

uint16_t reverse_bit(uint16_t value, uint16_t REV) {
	uint16_t ret = 0;
	for (size_t i = 0; i < REV; i++) {
		ret <<= 1;
		ret |= ((value >> i) & 1);
	}
	return ret;
}

bool check_same_uint16(uint16_t* r_ref, uint16_t* r, size_t length)
{
	size_t i;
	bool same = true;
	for (i = 0; i < length; i++) {
		if (r_ref[i] != r[i]) {
			printf("pos=%ld, r_ref=%d, while r=%d\n", i, r_ref[i], r[i]);
			same = false;
		}
	}
	if (same == true) {
		printf("check pass, r_ref == r!\n");
	}
	else {
		printf("check failed, r_ref != r!\n");
	}

	return same;
}

bool check_same_int16(int16_t* r_ref, int16_t* r, size_t length)
{
	size_t i;
	bool same = true;
	for (i = 0; i < length; i++) {
		if (r_ref[i] != r[i]) {
			printf("pos=%ld, r_ref=%d, while r=%d\n", i, r_ref[i], r[i]);
			same = false;
		}
	}
	if (same == true) {
		printf("check pass, r_ref == r!\n");
	}
	else {
		printf("check failed, r_ref != r!\n");
	}

	return same;
}

bool check_same_uint8(uint8_t* r_ref, uint8_t* r, size_t length)
{
	size_t i;
	bool same = true;
	for (i = 0; i < length; i++) {
		if (r_ref[i] != r[i]) {
			printf("pos=%ld, r_ref=%d, while r=%d\n", i, r_ref[i], r[i]);
			same = false;
		}
	}
	if (same == true) {
		printf("check pass, r_ref == r!\n");
	}
	else {
		printf("check failed, r_ref != r!\n");
	}

	return same;
}

bool check_same_int8(int8_t* r_ref, int8_t* r, size_t length)
{
	size_t i;
	bool same = true;
	for (i = 0; i < length; i++) {
		if (r_ref[i] != r[i]) {
			printf("pos=%ld, r_ref=%d, while r=%d\n", i, r_ref[i], r[i]);
			same = false;
		}
	}
	if (same == true) {
		printf("check pass, r_ref == r!\n");
	}
	else {
		printf("check failed, r_ref != r!\n");
	}

	return same;
}

/*************************************************
* Name:        bitreverse_standard_pos_transfer
*
* Description: Transfer a efficient array in bitreverse order to standard order
*              or transfer a efficient array in standard order to bitreverse order 
*
* Arguments:   - uint16_t a[FALCON_N]: input/output coefficient array
**************************************************/
void bitreverse_standard_pos_transfer(uint16_t a[FALCON_N]) {
	uint32_t i;
	uint16_t temp[FALCON_N] = { 0 };

	for (i = 0; i < FALCON_N; i++) {
		temp[i] = a[reverse_bit(i, LOGN)];
	}
	copy_array(a, temp);
}

void generate_random_uint8(unsigned char *x, unsigned long long xlen, uint32_t bound) {
    size_t i;

    srand((unsigned int)time(NULL));
    for (i = 0; i < xlen; ++i) {
		if(bound == 0) {
			x[i] = (unsigned char)rand();
		}
        else {
			x[i] = (unsigned char)(rand() % bound);
		}
    }
}

void generate_inorder_uint8(unsigned char *x, unsigned long long xlen, uint32_t bound) {
    size_t i;

    for (i = 0; i < xlen; ++i) {
		if(bound == 0) {
			x[i] = (unsigned char)i;
		}
        else {
			x[i] = (unsigned char)(i % bound);
		}
    }
}

void generate_random_uint16(uint16_t *x, unsigned long long xlen, uint32_t bound) {
    size_t i;

    srand((unsigned int)time(NULL));
    for (i = 0; i < xlen; ++i) {
		if(bound == 0) {
        	x[i] = (uint16_t)rand();
		}
		else {
			x[i] = (uint16_t)(rand() % bound);
		}
    }
}

void generate_inorder_uint16(uint16_t *x, unsigned long long xlen, uint32_t bound) {
    size_t i;

    for (i = 0; i < xlen; ++i) {
		if(bound == 0) {
        	x[i] = (uint16_t)i;
		}
		else {
			x[i] = (uint16_t)(i % bound);
		}
    }
}


void generate_random_int16(int16_t *x, unsigned long long xlen, uint32_t bound) {
    size_t i;

    srand((unsigned int)time(NULL));
    for (i = 0; i < xlen; ++i) {
		if(bound == 0) {
        	int16_t a = (int16_t)rand();
			int16_t b = (int16_t)rand();
			x[i] = (int16_t) (a - b);
		}
		else {
        	int16_t a = (int16_t)rand();
			int16_t b = (int16_t)rand();
			x[i] = (int16_t) ((a - b) % bound);
		}
    }
}


void generate_inorder_int16(int16_t *x, unsigned long long xlen, uint32_t bound) {
    size_t i;

    for (i = 0; i < xlen; ++i) {
		if(bound == 0) {
        	x[i] = (int16_t)i;
		}
		else {
			x[i] = (int16_t)(i % bound);
		}
    }
}

void print_array_int8(int8_t* a, size_t length){
	int i;
	printf("start print\n");
	for(i = 0; i < length; i++) {
		if(i != (length-1)) {
			printf("%d,", a[i]);
		}
		else {
			printf("%d", a[i]);
		}
		if((i+1)%32 == 0) {
			printf("\n");
		}
	}
}

void print_array_uint16(uint16_t* a, size_t length){
	int i;
	printf("start print\n");
	for(i = 0; i < length; i++) {
		if(i != (length-1)) {
			printf("%d,", a[i]);
		}
		else {
			printf("%d", a[i]);
		}
		if((i+1)%32 == 0) {
			printf("\n");
		}
	}
}