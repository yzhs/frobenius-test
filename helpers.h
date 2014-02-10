#ifndef HELPERS_H
#define HELPERS_H

#define error(...) fprintf(stderr, __VA_ARGS__)

#include <gmp.h>

unsigned long split(mpz_t d, mpz_t n);

unsigned randint(unsigned low, unsigned high);

gmp_randstate_t r_state;

void init();

void cleanup();

#endif
