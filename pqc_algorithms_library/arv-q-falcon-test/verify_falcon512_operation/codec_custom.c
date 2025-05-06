#include "inner.h"
#include <riscv_vector.h>
#include "../apis/custom_inst_api.h"
#include "params_custom.h"

/* see inner.h */
size_t
Zf(modq_encode_custom)(
	void *out, size_t max_out_len,
	const uint16_t *x, unsigned logn)
{
	size_t n, out_len;
	uint8_t *buf;

	n = (size_t)1 << logn;
	out_len = ((n * 14) + 7) >> 3;
	if (out == NULL) {
		return out_len;
	}
	if (out_len > max_out_len) {
		return 0;
	}
	buf = out;
    uint8_t bits = 14;
    uint32_t vlmax = vsetvlmax_e16m1();
    uint32_t vlcp = (vlmax * bits) >> 3;
    for (size_t vl; n > 0; n -= vl, x += vl, buf += vlcp) {
        vl = vsetvl_e16m1(n);
        vuint16m1_t vec_value = vle16_v_u16m1(x, vl);
		vbool16_t mask = vmsgtu_vx_u16m1_b16(vec_value, 12288, vl);
		uint32_t r = vcpop_m_b16(mask, vl);
		if (r > 0) {
			return 0;
		}
        vuint8m1_t vec_out;
        vec_out = pack_vx_u16m1(vec_value, bits);
        vse8_v_u8m1(buf, vec_out, vlcp);
    }

	return out_len;
}

/* see inner.h */
size_t
Zf(modq_decode_custom)(
	uint16_t *x, unsigned logn,
	const void *in, size_t max_in_len)
{
	size_t n, in_len;
	const uint8_t *buf;

	n = (size_t)1 << logn;
	in_len = ((n * 14) + 7) >> 3;
	if (in_len > max_in_len) {
		return 0;
	}
	buf = in;
    uint8_t bits = 14;
    uint32_t vlmax = vsetvlmax_e16m1();
    uint32_t vlcp = (vlmax * bits) >> 3;
    for (size_t vl; n > 0; n -= vl, x += vl, buf += vlcp) {
        vuint8m1_t vec_value = vle8_v_u8m1(buf, vlcp);
        vl = vsetvl_e16m1(n);
        vuint16m1_t vec_out;
        vec_out = unpack_vx_u16m1(vec_value, bits);
		vbool16_t mask = vmsgtu_vx_u16m1_b16(vec_out, 12288, vl);
		uint32_t r = vcpop_m_b16(mask, vl);
		if (r > 0) {
			return 0;
		}
        vse16_v_u16m1(x, vec_out, vl);
    }

	return in_len;
}


/* see inner.h */
size_t
Zf(trim_i16_encode_custom)(
	void *out, size_t max_out_len,
	const int16_t *x, unsigned logn, unsigned bits)
{
	size_t n, out_len;
	int16_t minv, maxv;
	uint8_t *buf;

	n = (size_t)1 << logn;
	maxv = (1 << (bits - 1)) - 1;
	minv = -maxv;
	out_len = ((n * bits) + 7) >> 3;
	if (out == NULL) {
		return out_len;
	}
	if (out_len > max_out_len) {
		return 0;
	}
	buf = out;
    uint32_t vlmax = vsetvlmax_e16m1();
    uint32_t vlcp = (vlmax * bits) >> 3;
    for (size_t vl; n > 0; n -= vl, x += vl, buf += vlcp) {
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_value = vle16_v_i16m1(x, vl);
		vbool16_t mask = vmsgt_vx_i16m1_b16(vec_value, maxv, vl);
		uint32_t r = vcpop_m_b16(mask, vl);
		mask = vmslt_vx_i16m1_b16(vec_value, minv, vl);
		r += vcpop_m_b16(mask, vl);
		if (r > 0) {
			return 0;
		}
        vuint8m1_t vec_out;
        vec_out = pack_vx_i16m1(vec_value, bits);
        vse8_v_u8m1(buf, vec_out, vlcp);
    }

	return out_len;
}


/* see inner.h */
size_t
Zf(trim_i16_decode_custom)(
	int16_t *x, unsigned logn, unsigned bits,   
	const void *in, size_t max_in_len)
{
	size_t n, in_len;
	const uint8_t *buf;

	n = (size_t)1 << logn;
	in_len = ((n * bits) + 7) >> 3;
	if (in_len > max_in_len) {
		return 0;
	}
	buf = in;
    uint32_t vlmax = vsetvlmax_e16m1();
    uint32_t vlcp = (vlmax * bits) >> 3;
    for (size_t vl; n > 0; n -= vl, x += vl, buf += vlcp) {
        vuint8m1_t vec_value = vle8_v_u8m1(buf, vlcp);
        vl = vsetvl_e16m1(n);
        vint16m1_t vec_out;
		vec_out = unpacks_vx_i16m1(vec_value, bits);
        vse16_v_i16m1(x, vec_out, vl);
    }

	return in_len;
}


/* see inner.h */
size_t
Zf(trim_i8_encode_custom)(
	void *out, size_t max_out_len,
	const int8_t *x, unsigned logn, unsigned bits)
{
	size_t n, out_len;
	int8_t minv, maxv;
	uint8_t *buf;

	n = (size_t)1 << logn;
	maxv = (1 << (bits - 1)) - 1;
	minv = -maxv;
	out_len = ((n * bits) + 7) >> 3;
	if (out == NULL) {
		return out_len;
	}
	if (out_len > max_out_len) {
		return 0;
	}
	buf = out;
    uint32_t vlmax = vsetvlmax_e8m1();
    uint32_t vlcp = (vlmax * bits) >> 3;
    for (size_t vl; n > 0; n -= vl, x += vl, buf += vlcp) {
        vl = vsetvl_e8m1(n);
        vint8m1_t vec_value = vle8_v_i8m1(x, vl);
		vbool8_t mask = vmsgt_vx_i8m1_b8(vec_value, maxv, vl);
		uint32_t r = vcpop_m_b8(mask, vl);
		mask = vmslt_vx_i8m1_b8(vec_value, minv, vl);
		r += vcpop_m_b8(mask, vl);
		if (r > 0) {
			return 0;
		}
        vuint8m1_t vec_out;
        vec_out = pack_vx_i8m1(vec_value, bits);
        vse8_v_u8m1(buf, vec_out, vlcp);
    }

	return out_len;
}


/* see inner.h */
size_t
Zf(trim_i8_decode_custom)(
	int8_t *x, unsigned logn, unsigned bits,
	const void *in, size_t max_in_len)
{
	size_t n, in_len;
	const uint8_t *buf;

	n = (size_t)1 << logn;
	in_len = ((n * bits) + 7) >> 3;
	if (in_len > max_in_len) {
		return 0;
	}

	if(bits != 8) {
		buf = in;
    	uint32_t vlmax = vsetvlmax_e8m1();
    	uint32_t vlcp = (vlmax * bits) >> 3;
		// int8_t mask1 = (uint8_t)1 << (bits - 1);
		// int8_t mask2 = -mask1;
    	for (size_t vl; n > 0; n -= vl, x += vl, buf += vlcp) {
        	vuint8m1_t vec_value = vle8_v_u8m1(buf, vlcp);
        	vl = vsetvl_e8m1(n);
        	vint8m1_t vec_out;
			vec_out = unpacks_vx_i8m1(vec_value, bits);
        	vse8_v_i8m1(x, vec_out, vl);
    	}
	}
	else {
		memcpy(x, in, n);
	}
    
	return in_len;
}