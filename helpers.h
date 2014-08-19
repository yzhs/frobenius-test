#ifndef HELPERS_H
#define HELPERS_H
#include <stdlib.h>
#include <gmp.h>

#include "common.h"

#define MODULUS n, b, c
#define MODULUS_ARGS const mpz_t n, const mpz_t b, const mpz_t c

#define POLY(name) name##_x, name##_1
#define CONST_POLY_ARGS(name) const mpz_t name##_x, const mpz_t name##_1
#define POLY_ARGS(name) mpz_t name##_x, mpz_t name##_1

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
