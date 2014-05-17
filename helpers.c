#include <stdlib.h>
#include <sys/time.h>

#include "helpers.h"

gmp_randstate_t r_state;


/*
 * Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd.
 */
void split(unsigned long *s, mpz_t d, mpz_t n)
{
	*s = 0;
	mpz_sub_ui(d, n, 1);

	while (mpz_even_p(d)) {
		(*s)++;
		mpz_fdiv_q_2exp(d, d, 1); // divide by 2^1
	}
}

/**
 * This function generates a random integer between 2 and n-2.
 */
void get_random(mpz_t result, mpz_t n)
{
	mpz_t tmp;

	mpz_init(tmp);
	mpz_sub_ui(tmp, n, 3);

	/* generate a random number between 0 and tmp-1 */
	mpz_urandomm(result, r_state, tmp);

	mpz_add_ui(result, result, 2);
	mpz_clear(tmp);
}

/*
 * Generate a random integer between [low] and [high]+1, both included.
 */
unsigned randint(unsigned low, unsigned high)
{
	// TODO make sure the distribution is uniform
	return low + (unsigned)rand() % (high - low + 1);
}

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

	return 0;
}

int cleanup(void)
{
	gmp_randclear(r_state);
	return 0;
}
