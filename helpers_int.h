/* helpers_int.h -- helper functions for the long long version of tests
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

#ifndef HELPERS_INT_H
#define HELPERS_INT_H
#include <stdbool.h>
#include <stdio.h>

#include "common.h"


/* Multiply two numbers modulo n. */
#define mul(x, y) mul_mod_n(x, y, n)
uint64_t mul_mod_n(uint64_t x, uint64_t y, uint64_t n);

/*
 * A few macros for passing commonly used arguments and declaring the
 * corresponding parameters.
 */
#define POLY_int(f) f##_x, f##_1
#define POLY_ARGS_int(f) uint64_t *f##_x, uint64_t *f##_1
#define CONST_POLY_ARGS_int(f) const uint64_t f##_x, const uint64_t f##_1
#define MODULUS_int n, b, c
#define MODULUS_ARGS_int const uint64_t n, const uint64_t b, const uint64_t c


/* Initialises the random number generator using a (static) seed. */
int init_int(void);

/* Replacements for some GMP functions. */
#define even(n) ((n) % 2lu == 0)
#define odd(n) ((n) % 2lu == 1)

/*
 * Calculate the greatest common divisor of a and b using the Euclidean
 * algorithm.
 */
uint64_t gcd(uint64_t a, uint64_t b);

/* Compute the jacobi symbol (x/y) of x over y. */
int jacobi(uint64_t x, uint64_t y);


/*
 * Given an integer 0 ≤ n < 2³², find the largest integer r such that r² ≤ n.
 */
uint64_t int_sqrt(const uint64_t n);

/*
 * Use the integer square root function to figure out whether a number is a
 * perfect square.
 */
bool is_square(const uint64_t n);


/* This function generates a random integer between low and high. Both low and
 * high are included. */
uint64_t get_random_int(uint64_t low, uint64_t high);


/* Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd. */
void split_int(uint64_t *s, uint64_t *d, const uint64_t n);

#endif
