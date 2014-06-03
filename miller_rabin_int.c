#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

#include "helpers_int.h"
#include "miller_rabin_int.h"

/*
 * Raise b to the e'th power modulo m.  This uses 64-bit registers to hold the
 * results of the multipliations.  Therefore, the results will be wrong if m is
 * greater than 2^32-1
 */
static unsigned long powm(unsigned long b, unsigned long e, unsigned m)
{
	unsigned long result = 1;

	while (LIKELY(e != 0)) {
		if (odd(e))
			result = (result * b) % m;
		b = (b * b) % m;
		e /= 2;
	}
	return result;
}

/*
 * This function checks whether a given number n is a prime or not, using the
 * Miller-Rabin primality test.  This is a probabilistic test which randomly
 * chooses an integer a as a base and checks whether n satisfies a certain
 * property (which depends on b).  If it does, n is a prime for at least three
 * out of four of the possible values of a, if it does not, it is certainly not
 * prime.
 *
 * The implementation is taken from the pseudo code found on
 * http://en.wikipedia.org/wiki/Miller-Rabin_primality_test.
 *
 * The function returns `probably_prime` if it found no evidence, that n might
 * be composite and `composite` if it did find a counter example.
 */
Primality miller_rabin_int(const unsigned n, const unsigned k)
{
	unsigned long s;
	unsigned long a, d, x, nm1;

	/* We need an odd integer greater than 3 */
	if (even(n))
		return n == 2 ? prime : composite;
	if (n == 3)
		return prime;
	else if (n < 3)
		return composite;

	nm1 = n - 1;

	/* compute s and d s.t. n-1=2^s*d */
	split_int(&s, &d, n);

	/* Repeat the test itself k times to increase the accuracy */
	for (unsigned i = 0; i < k; i++) {
		a = get_random_int(2, n - 2);

		/* compute a^d mod n */
		x = powm(a, d, n);

		if (x == 1 || x == nm1)
			continue;

		for (unsigned long r = 1; r <= s; r++) {
			//x = powm(x, 2, n);
			x = (x * x) % n;
			if (x == 1)
				return composite;
			if (x == nm1)
				break;
		}

		if (x != nm1)
			return composite;
	}

	return probably_prime;
}
