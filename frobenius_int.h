#ifndef FROBENIUS_INT_H
#define FROBENIUS_INT_H

#include "common.h"

Primality steps_3_4_5_int(const uint64_t n, const uint64_t b, const uint64_t c);
Primality QFT_int(const unsigned n, const unsigned b, const unsigned c);
Primality RQFT_int(const unsigned n, const unsigned k);

#endif
