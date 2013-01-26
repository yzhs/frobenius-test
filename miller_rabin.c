#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <gmp.h>

#ifndef DEBUG
#define assert(x)
#endif

gmp_randstate_t r_state;

/**
 * This function has to be called before using any of the functions below.  It
 * initialises the random number generator using a (static) seed.
 */
void init()
{
	unsigned long int seed;

	seed = 123456;

	gmp_randinit_default (r_state);
	gmp_randseed_ui(r_state, seed);

//	gmp_randclear(r_state);
}

/**
 * This function generates a random integer between 2 and n-2.
 */
void get_random(mpz_t n, mpz_t result)
{
	mpz_t tmp;
	mpz_inits(tmp, result, NULL);
	mpz_sub_ui(tmp, n, 3);

	/* generate a random number between 0 and tmp-1 */
	mpz_urandomm(result, r_state, tmp);

	mpz_add_ui(result, result, 2);
	mpz_clear(tmp);
}

/** 
 * Find numbers s and d such that n-1 == 2^s*d.
 */
void split_n_minus_1(mpz_t n, unsigned long *s, mpz_t d)
{
	mpz_t nm1;
	mpz_init(nm1);

	mpz_sub_ui(nm1, n, 1);

	/* count the number of trailing zeros */
	*s = mpz_scan1(nm1, 0);
	/* divide n-1 by 2^s yealding d */
	mpz_tdiv_q_2exp(d, nm1, *s);
}

/**
 * This function checks whether a given number n is a prime or not, using the
 * Miller-Rabin primality test.  This is a probabilistic test which randomly
 * chooses an integer a as a base and checks whether n satisfies a certain
 * property (which depends on b).  If it does, n is a prime for at least three
 * out of four of the possible values of a, if it does not, it is certainly not
 * prime.
 * The implementation is taken from the following pseudo code found on
 * http://en.wikipedia.org/wiki/Miller-Rabin_primality_test:
 *
 *   Input: n > 3, an odd integer to be tested for primality;
 *   Input: k, a parameter that determines the accuracy of the test
 *   Output: composite if n is composite, otherwise probably prime
 *   write n − 1 as 2^s·d with d odd by factoring powers of 2 from n − 1
 *   LOOP: repeat k times:
 *      pick a random integer a in the range 2, n − 2
 *      x ← a^d mod n
 *      if x = 1 or x = n − 1 then do next LOOP
 *      for r = 1 .. s
 *         x ← x^2 mod n
 *         if x = 1 then return composite
 *         if x = n − 1 then do next LOOP
 *      return composite
 *   return probably prime
 *
 * The function returns true if it found no evidence, that n might be composite
 * and false if it found a counter example.
 */
bool miller_rabin(mpz_t n, unsigned long k)
{
	unsigned long s;
	mpz_t a, d, x, nm1;
	mpz_inits(a, d, x, nm1, NULL);
	mpz_sub_ui(nm1, n, 1);

	/* We need an odd integer */
	assert(mpz_odd_p(n));
        /* greater than 3 */
	assert(mpz_cmp_ui(n, 3) > 0);

	/* compute s and d s.t. n-1=2^s*d */
	split_n_minus_1(n, &s, d);

	/* Repeat the test itself k times to increase the accuracy */
	for (unsigned long i = 0; i < k; i++) {
		get_random(n, a);

		/* compute a^d mod n */
		mpz_powm(x, a, d, n);

		if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, nm1) == 0)
			continue;

		for (unsigned long r = 1; r <= s; r++) {
			mpz_powm_ui(x, x, 2, n);
			if (mpz_cmp_ui(x, 1) == 0) return false;
			if (mpz_cmp(x, nm1) == 0) break;
		}
		if (mpz_cmp(x, nm1) != 0)
			return false;

	}

	return true;
}


int main()
{
	unsigned long i, k = 1000000;
	mpz_t foo, tmp;
	mpz_inits(foo, tmp, NULL);
	init();

	mpz_set_ui(foo, 23456789);
	mpz_out_str(stdout, 10, foo);
	if (miller_rabin(foo, k))
		puts(" is probably prime");
	else
		puts(" is composite");

	mpz_set_ui(foo, 23456799);
	mpz_out_str(stdout, 10, foo);
	if (miller_rabin(foo, k))
		puts(" is probably prime");
	else
		puts(" is composite");
}
