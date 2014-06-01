#ifndef FROBENIUS_H
#define FROBENIUS_H

#include <gmp.h>

#include "common.h"

Primality QFT(const mpz_t n, const mpz_t b, const mpz_t c);
Primality RQFT(const mpz_t n, const unsigned k);

#endif

