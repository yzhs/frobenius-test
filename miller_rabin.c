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

/*
 * This function checks whether a given number n is a prime or not, using the
 * Miller-Rabin primality test.  This is a probabilistic test which randomly
 * chooses an integer a as a base and checks whether n satisfies a certain
 * property (which depends on b).  If it does, n is a prime for at least three
 * out of four of the possible values of a, if it does not, it is certainly not
 * prime.
 * The implementation is taken from the pseudo code found on
 * http://en.wikipedia.org/wiki/Miller-Rabin_primality_test.
 * The function returns `probably_prime` if it found no evidence, that n might
 * be composite and `composite` if it did find a counter example.
 */
Primality miller_rabin(const mpz_t n, const unsigned k)
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
	for (unsigned i = 0; i < k; i++) {
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
