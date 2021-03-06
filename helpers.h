/* helpers.h -- helper functions for the gmp version of tests
 *
 * Copyright 2014 by Colin Benner <colin-software@yzhs.de>
 *
 * This file is part of frobenius-test.
 *
 * frobenius-test is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * frobenius-test is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with frobenius-test.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef HELPERS_H
#define HELPERS_H
#include <stdlib.h>
#include <gmp.h>

#include "common.h"

/*
 * A few macros for passing commonly used arguments and declaring the
 * corresponding parameters.
 */
#define MODULUS n, b, c
#define MODULUS_ARGS const mpz_t n, const mpz_t b, const mpz_t c

#define POLY(name) name##_x, name##_1
#define CONST_POLY_ARGS(name) const mpz_t name##_x, const mpz_t name##_1
#define POLY_ARGS(name) mpz_t name##_x, mpz_t name##_1

// Some macros to simplify handling polyonmial arithmetic.
#define right_shift mpz_fdiv_q_2exp

#define poly_mul(res, f, g) mult_mod(res##_x, res##_1, f##_x, f##_1, g##_x, g##_1, MODULUS)
#define poly_sqr(res, f) square_mod(res##_x, res##_1, f##_x, f##_1, MODULUS)
#define poly_powm(res, f, exponent) powm(res##_x, res##_1, f##_x, f##_1, (exponent), MODULUS)
#define poly_invert(res, f) invert(res##_x, res##_1, f##_x, f##_1, MODULUS)

#define poly_eq(f, g) (mpz_cmp(f##_x, g##_x) == 0 && mpz_cmp(f##_1, g##_1) == 0)

#define poly_mul_x(res, f) mult_x_mod(res##_x, res##_1, f##_x, f##_1, MODULUS)
#define poly_pow_x(res, exponent) power_of_x(res##_x, res##_1, (exponent), MODULUS)
#define poly_sigma(res, f) sigma(res##_x, res##_1, f##_x, f##_1, MODULUS)

// Some temporary variables.
extern mpz_t bb4c, base_x, base_1, tmp0, tmp1, tmp2;

// State of GMPs random number generator
extern gmp_randstate_t r_state;


// Seed random number generators.
int init(void);

// Whatever is to be done before program shutdown.
int cleanup(void);


// Compute s, d such that 2^s * d == num - 1.
void split(uint64_t *s, mpz_t d, const mpz_t num);

// Randomly select an integer`result` in the interval [2, n-2].
void get_random(mpz_t result, const mpz_t n);

// Returm a randomly selected integer from the interval [low, high].
unsigned randint(const unsigned low, const unsigned high);

#endif
