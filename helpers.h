#ifndef HELPERS_H
#define HELPERS_H
#include <stdlib.h>
#include <gmp.h>

#include "common.h"

extern gmp_randstate_t r_state;

void split(unsigned long *s, mpz_t d, mpz_t n);
void get_random(mpz_t result, mpz_t n);
unsigned randint(unsigned low, unsigned high);

int init(void);
int cleanup(void);

#endif
