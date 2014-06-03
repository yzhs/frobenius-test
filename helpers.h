#ifndef HELPERS_H
#define HELPERS_H
#include <stdlib.h>
#include <gmp.h>

#include "common.h"

// State of GMPs random number generator
extern gmp_randstate_t r_state;


// Seed random number generators.
int init(void);

// Whatever is to be done before program shutdown.
int cleanup(void);


// Compute s, d such that 2^s * d == num - 1.
void split(unsigned long *s, mpz_t d, const mpz_t num);

// Randomly select an integer`result` in the interval [2, n-2].
void get_random(mpz_t result, const mpz_t n);

// Returm a randomly selected integer from the interval [low, high].
unsigned randint(const unsigned low, const unsigned high);

#endif
