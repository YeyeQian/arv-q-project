#ifndef CBD_H
#define CBD_H

#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "../apis/custom_inst_api.h"

#define cbd_even_index KYBER_NAMESPACE(_cbd_even_index)
extern const uint8_t cbd_even_index[ELEMENT_SEW8_PER_VECREG/2];

#define cbd_odd_index KYBER_NAMESPACE(_cbd_odd_index)
extern const uint8_t cbd_odd_index[ELEMENT_SEW8_PER_VECREG/2];

#define cbd_eta1 KYBER_NAMESPACE(_cbd_eta1)
void cbd_eta1(poly *r, const uint8_t buf[KYBER_ETA1*KYBER_N/4]);

#define cbd_eta2 KYBER_NAMESPACE(_cbd_eta2)
void cbd_eta2(poly *r, const uint8_t buf[KYBER_ETA2*KYBER_N/4]);

#define cbd_eta1_custom KYBER_NAMESPACE(_cbd_eta1_custom)
void cbd_eta1_custom(poly *r, const uint8_t buf[KYBER_ETA1*KYBER_N/4]);

#define cbd_eta2_custom KYBER_NAMESPACE(_cbd_eta2_custom)
void cbd_eta2_custom(poly *r, const uint8_t buf[KYBER_ETA1*KYBER_N/4]);

#endif
