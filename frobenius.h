#ifndef FROBENIUS_H
#define FROBENIUS_H

#include <gmp.h>

#include "common.h"

#define MODULUS n, b, c
#define MODULUS_ARGS const mpz_t n, const mpz_t b, const mpz_t c

#define POLY(name) name##_x, name##_1
#define CONST_POLY_ARGS(name) const mpz_t name##_x, const mpz_t name##_1
#define POLY_ARGS(name) mpz_t name##_x, mpz_t name##_1

Primality QFT(const mpz_t n, const mpz_t b, const mpz_t c);
Primality RQFT(const mpz_t n, const unsigned k);

#endif

