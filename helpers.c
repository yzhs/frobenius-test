#include <stdlib.h>
#include <sys/time.h>

#include "helpers.h"


/* 
 * Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd.
 */
unsigned long split(mpz_t d, mpz_t n)
{
	unsigned long s = 0;
	mpz_sub_ui(d, n, 1);

	while (mpz_odd_p(d)) {
		s++;
		mpz_fdiv_q_2exp(d, d, 1); // divide by 2^1
	}
	
	return s;
}

/*
 * Generate a random integer between [low] and [high].
 */
unsigned randint(unsigned low, unsigned high)
{
	// TODO make sure the distribution is uniform
	return low + (unsigned)rand() % (high - low);
}

gmp_randstate_t r_state;

void init()
{
	unsigned long int seed;
	struct timeval tv;

	seed = 123457;
	gettimeofday(&tv, NULL);
	seed = tv.tv_usec;

	gmp_randinit_default (r_state);
	gmp_randseed_ui(r_state, seed);

	srand(seed);
}

void cleanup()
{
	gmp_randclear(r_state);
}
