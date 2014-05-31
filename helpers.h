#ifndef HELPERS_H
#define HELPERS_H
#include <stdlib.h>
#include <gmp.h>

#include "common.h"

extern gmp_randstate_t r_state;

void split(unsigned long *s, mpz_t d, const mpz_t n);
void get_random(mpz_t result, const mpz_t n);
unsigned randint(const unsigned low, const unsigned high);

int init(void);
int cleanup(void);

#endif
