#include "benes.h"

static inline uint32_t load4(const unsigned char *src)
{
	uint32_t a;

	a = src[3];
	a = src[2] | (a<<8);
	a = src[1] | (a<<8);
	a = src[0] | (a<<8);
	return a;
}

static inline void transpose_8x8_x4_in(uint32_t * out)
{
  uint32_t i0=out[0];
  uint32_t i1=out[1];
  uint32_t i2=out[2];
  uint32_t i3=out[3];
  uint32_t i4=out[4];
  uint32_t i5=out[5];
  uint32_t i6=out[6];
  uint32_t i7=out[7];
  uint32_t x,y;
  uint32_t mask;
  // 4x4
  mask = 0x0f0f0f0f;
  x = (i0&mask)|((i4&mask)<<4);
  y = ((i0&~mask)>>4)|(i4&~mask);
  i0 = x; i4 = y;
  x = (i1&mask)|((i5&mask)<<4);
  y = ((i1&~mask)>>4)|(i5&~mask);
  i1 = x; i5 = y;
  x = (i2&mask)|((i6&mask)<<4);
  y = ((i2&~mask)>>4)|(i6&~mask);
  i2 = x; i6 = y;
  x = (i3&mask)|((i7&mask)<<4);
  y = ((i3&~mask)>>4)|(i7&~mask);
  i3 = x; i7 = y;
  // 2x2
  mask = 0x33333333;
  x = (i0&mask)|((i2&mask)<<2);
  y = ((i0&~mask)>>2)|(i2&~mask);
  i0 = x; i2 = y;
  x = (i1&mask)|((i3&mask)<<2);
  y = ((i1&~mask)>>2)|(i3&~mask);
  i1 = x; i3 = y;
  x = (i4&mask)|((i6&mask)<<2);
  y = ((i4&~mask)>>2)|(i6&~mask);
  i4 = x; i6 = y;
  x = (i5&mask)|((i7&mask)<<2);
  y = ((i5&~mask)>>2)|(i7&~mask);
  i5 = x; i7 = y;
  // 1x1
  mask = 0x55555555;
  x = (i0&mask)|((i1&mask)<<1);
  y = ((i0&~mask)>>1)|(i1&~mask);
  i0 = x; i1 = y;
  x = (i2&mask)|((i3&mask)<<1);
  y = ((i2&~mask)>>1)|(i3&~mask);
  i2 = x; i3 = y;
  x = (i4&mask)|((i5&mask)<<1);
  y = ((i4&~mask)>>1)|(i5&~mask);
  i4 = x; i5 = y;
  x = (i6&mask)|((i7&mask)<<1);
  y = ((i6&~mask)>>1)|(i7&~mask);
  i6 = x; i7 = y;
  out[0] = i0;
  out[1] = i1;
  out[2] = i2;
  out[3] = i3;
  out[4] = i4;
  out[5] = i5;
  out[6] = i6;
  out[7] = i7;
}


/* input: out, a 32x32 matrix over GF(2) */
/* output: out, transpose of out */
void transpose_32x32_in(uint32_t * out)
{
  int i;
  uint32_t x, y;
  uint32_t mask;
  transpose_8x8_x4_in(out);
  transpose_8x8_x4_in(out+8);
  transpose_8x8_x4_in(out+16);
  transpose_8x8_x4_in(out+24);
  mask = 0x00ff00ff;
  for(i=0;i<8;i++) {
    x = (out[i]&mask)|((out[8+i]&mask)<<8);
    y = ((out[i]&~mask)>>8)|((out[8+i]&~mask));
    out[i] = x; out[8+i] = y;
  }
  for(i=0;i<8;i++) {
    x = (out[16+i]&mask)|((out[16+8+i]&mask)<<8);
    y = ((out[16+i]&~mask)>>8)|((out[16+8+i]&~mask));
    out[16+i] = x; out[16+8+i] = y;
  }
  mask = 0x0000ffff;
  for(i=0;i<16;i++) {
    x = (out[i]&mask)|((out[16+i]&mask)<<16);
    y = ((out[i]&~mask)>>16)|((out[16+i]&~mask));
    out[i] = x; out[16+i] = y;
  }
}

static void layer_u32(uint32_t data[], const uint32_t * bits, int lgs)
{
	int i, j, s;

	uint32_t d;

	s = 1 << lgs;

	for (i = 0; i < 128; i += s*2)
	for (j = i; j < i+s; j++)
	{
		d = (data[j+0] ^ data[j+s]);
		d &= (*bits++);
		data[j+0] ^= d;
		data[j+s] ^= d;
	}
}

#define _MULTI_LAYERS_

#if defined(_MULTI_LAYERS_)
static void layer_u32_x2_inc(uint32_t data[], const uint32_t * bits0, const uint32_t *bits1, int lgs0 )
{
	int i, j, s0, s1;

	s0 = 1 << lgs0;
	s1 = s0 << 1;

	for (i = 0; i < 128; i += s1*2) {
		const uint32_t *bits0_ = bits0 + s0;
		const uint32_t *bits1_ = bits1 + s0;
		for (j = i; j < i+s0; j++) {
			uint32_t d0 = data[j+0];
			uint32_t d1 = data[j+s0];
			uint32_t d2 = data[j+s1];
			uint32_t d3 = data[j+s1+s0];
			uint32_t tmp0, tmp1;

			tmp0 = d0^d1;
			tmp1 = d2^d3;
			tmp0 &= (*bits0++);
			tmp1 &= (*bits0_++);

			d0 ^= tmp0;
			d1 ^= tmp0;
			d2 ^= tmp1;
			d3 ^= tmp1;

			tmp0 = d0^d2;
			tmp1 = d1^d3;
			tmp0 &= (*bits1++);
			tmp1 &= (*bits1_++);

			d0 ^= tmp0;
			d1 ^= tmp1;
			d2 ^= tmp0;
			d3 ^= tmp1;

			data[j+0] = d0;
			data[j+s0] = d1;
			data[j+s1] = d2;
			data[j+s1+s0] = d3;
		}
		bits0 = bits0_;
		bits1 = bits1_;
	}
}

static void layer_u32_x2_dec(uint32_t data[], const uint32_t * bits0, const uint32_t *bits1, int lgs0 )
{
	int i, j, s0, s1;

	s0 = 1 << (lgs0-1);
	s1 = s0 << 1;

	for (i = 0; i < 128; i += s1*2) {
		const uint32_t *bits0_ = bits0 + s0;
		const uint32_t *bits1_ = bits1 + s0;
		for (j = i; j < i+s0; j++) {
			uint32_t d0 = data[j+0];
			uint32_t d1 = data[j+s0];
			uint32_t d2 = data[j+s1];
			uint32_t d3 = data[j+s1+s0];
			uint32_t tmp0,tmp1;

			tmp0 = d0^d2;
			tmp1 = d1^d3;
			tmp0 &= (*bits0++);
			tmp1 &= (*bits0_++);
			d0 ^= tmp0;
			d1 ^= tmp1;
			d2 ^= tmp0;
			d3 ^= tmp1;

			tmp0 = d0^d1;
			tmp1 = d2^d3;
			tmp0 &= (*bits1++);
			tmp1 &= (*bits1_++);
			d0 ^= tmp0;
			d1 ^= tmp0;
			d2 ^= tmp1;
			d3 ^= tmp1;

			data[j+0] = d0;
			data[j+s0] = d1;
			data[j+s1] = d2;
			data[j+s1+s0] = d3;
		}
		bits0 = bits0_;
		bits1 = bits1_;
	}
}
#else
#endif

void benes_4096(uint32_t * r32, const unsigned char * bits, int rev)
{
	int i, iter, inc;

	const unsigned char *bits_ptr;

	if (rev) { bits_ptr = bits + (2*12-2)*256; inc = -512; }
	else     { bits_ptr = bits;         inc = 0;    }

	uint32_t data[128];
	uint32_t conbit[64];

////////////////
	for (int j=0;j<128;j+=64) {
		for (i = 0; i < 32; i++)
		{
			data[j+ i]    = r32[j+ i*2 + 0];
			data[j+ 32+i] = r32[j+ i*2 + 1];
		}
		transpose_32x32_in( &data[j] );
		transpose_32x32_in( &data[j+32] );
	}

	for (iter = 0; iter <= 4; iter++)
	{
		for (i = 0; i < 64; i++) { conbit[i] = load4(bits_ptr); bits_ptr += 4; }
		bits_ptr += inc;

		transpose_32x32_in( conbit );
		transpose_32x32_in( &conbit[32] );

		layer_u32( data, conbit, iter);
	}
	for (int j=0;j<128;j+=64) {
		transpose_32x32_in( &data[j] );
		transpose_32x32_in( &data[j+32] );

		for (i = 0; i < 32; i++)
		{
			r32[j+ i*2 + 0] = data[j+ i];
			r32[j+ i*2 + 1] = data[j+ 32+i];
		}
	}

///////////////////////////////

#if defined(_MULTI_LAYERS_)
	for (iter = 0; iter < 6; iter+=2)
	{
		uint32_t conbit2[64];
		for (i = 0; i < 64; i++) { conbit[i] = load4(bits_ptr); bits_ptr += 4; }
		bits_ptr += inc;
		for (i = 0; i < 64; i++) { conbit2[i] = load4(bits_ptr); bits_ptr += 4; }
		bits_ptr += inc;

		layer_u32_x2_inc(r32, conbit, conbit2, iter);
	}
	for (iter = 6; iter <= 6; iter++)
#else
	for (iter = 0; iter <= 6; iter++)
#endif
	{
		for (i = 0; i < 64; i++) { conbit[i] = load4(bits_ptr); bits_ptr += 4; }
		bits_ptr += inc;

		layer_u32(r32, conbit, iter);
	}

#if defined(_MULTI_LAYERS_)
	for (iter = 5; iter >= 0; iter-=2)
	{
		uint32_t conbit2[64];
		for (i = 0; i < 64; i++) { conbit[i] = load4(bits_ptr); bits_ptr += 4; }
		bits_ptr += inc;
		for (i = 0; i < 64; i++) { conbit2[i] = load4(bits_ptr); bits_ptr += 4; }
		bits_ptr += inc;

		layer_u32_x2_dec(r32, conbit, conbit2, iter );
	}
#else
	for (iter = 5; iter >= 0; iter--)
	{
		for (i = 0; i < 64; i++) { conbit[i] = load4(bits_ptr); bits_ptr += 4; }
		bits_ptr += inc;

		layer_u32(r32, conbit, iter);
	}
#endif

////////////////////////////////

	for (int j=0;j<128;j+=64) {
		for (i = 0; i < 32; i++)
		{
			data[j+ i]    = r32[j+ i*2 + 0];
			data[j+ 32+i] = r32[j+ i*2 + 1];
		}
		transpose_32x32_in( &data[j] );
		transpose_32x32_in( &data[j+32] );
	}

	for (iter = 4; iter >=0; iter--)
	{
		for (i = 0; i < 64; i++) { conbit[i] = load4(bits_ptr); bits_ptr += 4; }

		bits_ptr += inc;

		transpose_32x32_in( conbit );
		transpose_32x32_in( &conbit[32] );

		layer_u32( data, conbit, iter);
	}
	for (int j=0;j<128;j+=64) {
		transpose_32x32_in( &data[j] );
		transpose_32x32_in( &data[j+32] );

		for (i = 0; i < 32; i++)
		{
			r32[j+ i*2 + 0] = data[j+ i];
			r32[j+ i*2 + 1] = data[j+ 32+i];
		}
	}
}