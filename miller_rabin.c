#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "helpers.h"

#ifndef DEBUG
#undef assert
#define assert(x)
#endif

Primality miller_rabin(mpz_t n, int k);

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
Primality miller_rabin(mpz_t n, int k)
{
	Primality result = probably_prime;
	unsigned long s;
	int foo;
	mpz_t a, d, x, nm1;

	/* We need an odd integer */
	if (mpz_even_p(n))
		return mpz_cmp_ui(n, 2) == 0 ? prime : composite;

	/* greater than 3 */
	foo = mpz_cmp_ui(n, 3);
	if (foo == 0)
		return prime;
	else if (foo < 0)
		return composite;

	mpz_inits(a, d, x, nm1, NULL);
	mpz_sub_ui(nm1, n, 1);

	/* compute s and d s.t. n-1=2^s*d */
	split(&s, d, n);

	/* Repeat the test itself k times to increase the accuracy */
	for (int i = 0; i < k; i++) {
		get_random(a, n);

		/* compute a^d mod n */
		mpz_powm(x, a, d, n);

		if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, nm1) == 0)
			continue;

		for (unsigned long r = 1; r <= s; r++) {
			mpz_powm_ui(x, x, 2, n);
			if (mpz_cmp_ui(x, 1) == 0) {
				result = composite;
				goto exit;
			}
			if (mpz_cmp(x, nm1) == 0)
				break;
		}
		if (mpz_cmp(x, nm1) != 0) {
			result = composite;
			goto exit;
		}
	}

exit:
	mpz_clears(a, d, x, nm1, NULL);
	return result;
}

#ifndef TEST
int main()
{
	mpz_t foo, tmp;
	unsigned long upper_bound = (1lu << 32) - 1, dots_every = 1lu << 25;
	unsigned long counter = 0, i;

	//FILE *primes = fopen("primes_worker.txt", "a");

	mpz_inits(foo, tmp, NULL);
	init();

	for (i = 5; i < upper_bound; i += 2) {
		if (i % dots_every == 1) {  /* we need to check for == 1 since i will always be odd */
			printf(".");
			(void)fflush(stdout);
		}
		mpz_set_ui(foo, i);
		if (miller_rabin(foo, 1))
			counter++;
			//	fprintf(primes, "%lu\n", i);
	}

	//(void)fclose(primes);

	printf("\n%lu\n", counter);
	cleanup();
	mpz_clears(foo, tmp, NULL);
}
#endif
