#ifndef HELPERS_H
#define HELPERS_H
#include <stdlib.h>
#include <gmp.h>

#define die(...) do { fprintf(stderr, __VA_ARGS__); exit(1); } while (0)

#define len(a) (sizeof(a)/sizeof(a[0]))

extern gmp_randstate_t r_state;

void split(unsigned long *s, mpz_t d, mpz_t n);
void get_random(mpz_t n, mpz_t result);
unsigned randint(unsigned low, unsigned high);

int init(void);
int cleanup(void);

#endif
