#include <stdlib.h>
#include <gmp.h>

/* 
 * Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd.
 */
unsigned long split(mpz_t n, mpz_t d)
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
