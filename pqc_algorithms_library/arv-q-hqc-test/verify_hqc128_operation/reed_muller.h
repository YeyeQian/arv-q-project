#ifndef REED_MULLER_H
#define REED_MULLER_H

/**
 * @file reed_muller.h
 * @brief Header file of reed_muller.cpp
 */

#include "parameters.h"
#include <stddef.h>
#include <stdint.h>

void reed_muller_encode(uint64_t* cdw, const uint64_t* msg);
void reed_muller_decode(uint64_t* msg, const uint64_t* cdw);

//Customized Versions
void reed_muller_encode_custom(uint64_t *cdw, const uint64_t *msg);
void reed_muller_decode_custom(uint64_t *msg, const uint64_t *cdw);
#endif
