#include <stdlib.h>
#include <sys/time.h>

#include "helpers.h"

gmp_randstate_t r_state;

mpz_t bb4c, base_x, base_1, tmp0, tmp1, tmp2;

/*
 * Seed random number generators.
 */
int init(void)
{
	long int seed;
	struct timeval tv;

	//seed = 123457;
	gettimeofday(&tv, NULL);
	seed = tv.tv_usec;

	gmp_randinit_default(r_state);
	gmp_randseed_ui(r_state, seed);

	srand((unsigned int)seed);
	mpz_inits(bb4c, base_x, base_1, tmp0, tmp1, tmp2, NULL);

	return 0;
}

/*
 * Whatever is to be done before program shutdown.
 */
int cleanup(void)
{
	mpz_clears(bb4c, base_x, base_1, tmp0, tmp1, tmp2, NULL);
	gmp_randclear(r_state);
	return 0;
}


/*
 * Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd.
 */
void split(uint64_t *s, mpz_t d, const mpz_t n)
{
	// Find the least significant bit that is set, except for the 2^0 bit.
	*s = mpz_scan1(n, 1);

	mpz_fdiv_q_2exp(d, n, *s);
}

/**
 * This function generates a random integer between 2 and n-2.
 */
void get_random(mpz_t result, const mpz_t n)
{
	mpz_t tmp;

	mpz_init(tmp);
	mpz_sub_ui(tmp, n, 3);

	// Generate a random number between 0 and tmp-1.  The results of
	// mpz_urandomm are uniformly distributed.
	mpz_urandomm(result, r_state, tmp);

	mpz_add_ui(result, result, 2);
	mpz_clear(tmp);
}

/*
 * Generate a random integer between [low] and [high]+1, both included.
 */
unsigned randint(const unsigned low, const unsigned high)
{
	// TODO make sure the distribution is uniform
	return low + (unsigned)rand() % (high - low + 1);
}
