#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

#include "helpers_int.h"
#include "small_primes.h"

#define N 1

int enable_logging = 0;

#undef debug
#define debug(...) do {} while (0)


Primality QFT_int(unsigned long n, unsigned long b, unsigned long c);
Primality RQFT_int(unsigned long n, unsigned k);

static const unsigned B = 50000;


/*
 * Return f(x)*g(x) mod (n, x^2 - b*x - c) where f(x) = d*x + e and g(x) = f*x + g in the return arguments res0 and
 * res1, representing the polynomial res0*x + res1.
 */
static void mult_mod_int(unsigned long *res0, unsigned long *res1,
			 unsigned long d, unsigned long e,
			 unsigned long f, unsigned long g,
			 unsigned long n, unsigned long b, unsigned long c)
{
	unsigned long df, ef = (e * f) % n, eg = (e * g) % n;

	/*
	 * If deg f = 1, the whole thing amounts to multiplying the coefficients of g with a constant and reducing them
	 * modulo n.
	 */
	if (d == 0) {
		*res0 = ef;
		*res1 = eg;
		return;
	}
	df = (d * f) % n;

	// All these modulo operations are pretty expensive...
	*res0 = (((df * b) % n + (d * g) % n) % n + ef) % n;
	*res1 = ((df * c) % n + eg) % n;
}

/*
 * Calculate (dx + e)^2 mod (n, x^2-bx-c) returning the result as (*res0) * x + (*res1).
 */
static void square_mod_int(unsigned long *res0, unsigned long *res1, /* The resulting linear polynomial. */
			   unsigned long d, unsigned long e,
			   unsigned long n, unsigned long b, unsigned long c)
{
	unsigned long ee = (e * e) % n;
	unsigned long dd;

	if (d == 0) {
		*res0 = 0;
		*res1 = ee;
		return;
	}
	dd = (d * d) % n;
	// compute res0 = d^2*b+2*d*e
	*res0 = ((dd * b) % n + (2 * (d * e) % n) % n) % n;

	// and res1 = d^2*c+e^2
	*res1 = ((dd * c) % n + ee) % n;
}

/*
 * Calculate (base0 * x + base1)^exp mod (n, x^2-bx-c) returning the result as (*res0) * x + (*res1).
 */
static void powm_int(unsigned long *res0, unsigned long *res1,
		     unsigned long base0, unsigned long base1, unsigned long exp,
		     unsigned long n, unsigned long b, unsigned long c)
{
	*res0 = 0;
	*res1 = 1;

	while (exp != 0) {
		if (exp % 2 == 1)
			mult_mod_int(res0, res1, base0, base1, *res0, *res1, n, b, c);
		square_mod_int(&base0, &base1, base0, base1, n, b, c);
		exp /= 2;
	}
}

/*
 * Like QFT, return 'probably_prime' if n might be prime, 'prime' if n is
 * certainly prime and 'composite' if a proof for n's compositeness was found.
 */
static Primality steps_1_2_int(unsigned long n)
{
	unsigned long sqrt = int_sqrt(n);

	/*  (2) If n is a square, it can obviously not be prime. */
	if (sqrt * sqrt == n) {
		debug("found a square: %lu\n", n);
		return composite;
	}

	for (unsigned long i = 0; i < len(prime_list) && prime_list[i] <= sqrt; i++)
		if (n % prime_list[i] == 0) {
			debug("found a prime factor: %u divides %lu\n", prime_list[i], n);
			return composite;
		}

	/* If the given number is small enough, there cannot be a non-trivial
	 * divisor of n, whence n is prime. */
	return (sqrt < B) ? prime : probably_prime;
}

/*
 * Execute steps (3) through (5) of the Quadratic Frobenius Test.
 */
static Primality steps_3_4_5_int(unsigned long n, unsigned long b, unsigned long c)
{
	unsigned long x0, x1, s, tmp, foo0, foo1;
	unsigned long r, i;

	x0 = 1;
	x1 = s = tmp = foo0 = foo1 = 0;

	// Suppose n>1 is odd, (b^2+4c over n)=-1 and (-c over n)=1.

	/*
	 * (3) Compute x^((n+1)/2) mod (n, x^2-bx-c).  If x^((n+1)/2) not in ℤ/nℤ,
	 * declare n to be composite and stop.
	 */
	tmp = n + 1;    // tmp = n+1
	tmp = tmp / 2;  // tmp = (n+1)/2
	powm_int(&foo0, &foo1, x0, x1, tmp, n, b, c);
	if (foo0 != 0) {// check whether x^((n+1)/2) has degree 1
		debug("sqrt(-c) is not an integer: sqrt(-c) == %lu*x+%lu, n=%lu, b=%lu, c=%lu\n", foo0, foo1, n, b, c);
		return composite;
	}

	/*
	 * (4) Compute x^(n+1) mod (n, x^2-bx-c).  If x^(n+1) not congruent -c,
	 * declare n to be composite and stop.
	 */
	foo1 = foo1 * foo1;
	tmp = n - c;
	if (foo1 % n != tmp) {
		debug("x^(n+1) != c: x^(n+1) == %lu*x+%lu, n=%lu, b=%lu, c=%lu\n", foo0, foo1, n, b, c);
		return composite;
	}

	debug("composite = %i, probably_prime = %i, prime = %i\n", composite, probably_prime, prime);

	/*
	 * (5) Let n^2-1=2^r*s, where s is odd.  If x^s not congruent 1 mod (n,
	 * x^2-bx-c), and x^(2^j*s) not congruent -1 mod (n, x^2-bx-c) for all
	 * 0≤j≤r-2, declare n to be composite and stop.
	 * If n is not declared composite in Steps 1—5, declare n to be a
	 * probable prime.
	 */
	tmp = n * n;
	split_int(&r, &s, tmp); // calculate r,s such that 2^r*s + 1 == n^2
	if (tmp - 1 != (1lu << r) * s)
		die("split failed");
	powm_int(&foo0, &foo1, x0, x1, s, n, b, c);
	tmp = n - 1;
	if (foo0 == 0 && foo1 == 1) {
		debug("x^s == 1 mod (n, x^2-bx-c) where (b,c)=(%lu,%lu), so n=%lu is probably prime\n", b, c, n);
		return probably_prime;
	}

	for (i = 0; i < r - 1; i++) {
		if (foo0 == 0 && foo1 % n == tmp) {
			debug("x^(2^%lu*s) == n-1 mod (n, x^2-%lux-%lu), so n=%lu is probably prime\n", i, b, c, n);
			return probably_prime;
		}
		square_mod_int(&foo0, &foo1, foo0, foo1, n, b, c);
	}

	debug("n is certainly composite\n");
	return composite;
}

/*
 * Execute the Quadratic Frobenius Test with parameters (b,c).
 *
 * Returns 'prime' if n is certainly prime, 'probably_prime' if no evidence
 * could be found that n might be composite and 'composite' otherwise.
 */
Primality QFT_int(unsigned long n, unsigned long b, unsigned long c)
{
	Primality result = steps_1_2_int(n);

	if (result != probably_prime)
		return result;

	return steps_3_4_5_int(n, b, c);
}

/*
 * Run the randomized quadratic frobenius test to check whether [n] is a a
 * prime.  The Parameter [k] determines how many times the test will be run at
 * most.  If the test returns "composite", it will not be run again.
 */
Primality RQFT_int(unsigned long n, unsigned k)
{
	Primality result;
	unsigned long b = 0, c = 0;
	unsigned long bb4c, tmp;
	int j1 = 0, j2 = 0;

	if (n < 3)
		return n == 2 ? prime : composite;
	if (even(n))
		return composite;

	result = steps_1_2_int(n);
	/* If the number is found to be either composite or certainly prime, we can
	 * return that result immediately. */
	if (result != probably_prime)
		return result;

#ifdef check_non_trivial_divisor
#undef check_non_trivial_divisor
#endif
#define check_non_trivial_divisor(num) do { \
		tmp = gcd(num, n); \
		if (tmp != 1 && tmp != n) \
			return composite; \
} while (0)

	for (unsigned j = 0; j < k; j++) {
		for (unsigned i = 0; i < B; i++) {
			b = get_random_int(2, n - 2);
			c = get_random_int(2, n - 2);
			bb4c = ((b * b) % n + c * 4) % n;
			j1 = jacobi(bb4c, n);
			j2 = jacobi(n - c, n); /* Warning: n-c is not congruent to -c, since -c is interpreted as 2⁶⁴-c !!! */
			if (j1 == -1 && j2 == 1) {
				check_non_trivial_divisor(bb4c);
				check_non_trivial_divisor(b);
				check_non_trivial_divisor(c);
				break;
			}
		}
		if (j1 != -1 || j2 != 1) {
			printf("Found no suitable pair (b,c) modulo n=%lu. "
			       "This is highly unlikely unless the programme is wrong. "
			       "Assuming %lu is a prime...\n", n, n);
			return probably_prime;
		}

		result = steps_3_4_5_int(n, b, c);
		if (result == composite)
			return composite;
	}

	return probably_prime;
}

#ifndef TEST
static void *run(void *worker_id_cast_to_void_star)
{
	unsigned long upper_bound = (1lu << 32) - 1, dots_every = 1lu << 25;
	unsigned long counter = 0, i;
	unsigned long worker_id = (unsigned long)worker_id_cast_to_void_star;

#ifdef WRITE_PRIMES
	char str[32];
	FILE *primes;
	(void)snprintf(str, 32, "primes_worker_%lu.txt", worker_id);
	primes = fopen(str, "a");
#endif

	for (i = B * B + 1 + 2 * worker_id; i < upper_bound; i += 2 * N) {
		if (i % dots_every == 1) {  /* we need to check for == 1 since i will always be odd */
			printf(".");
			(void)fflush(stdout);
		}

		if (RQFT_int(i, 1)) {
			counter++;
#ifdef WRITE_PRIMES
			fprintf(primes, "%lu\n", i);
#endif
		}
	}

#ifdef WRITE_PRIMES
	(void)fclose(primes);
#endif
	return (void *)counter;
}

int main()
{
	unsigned i;
	unsigned counter = 0;
	pthread_t threads[N];

	init_int();

	for (i = 0; i < N; i++)
		if (0 != pthread_create(&threads[i], NULL, run, (void *)(unsigned long)i))
			die("failed to create thread %u, exiting\n", i);

	for (i = 0; i < N; i++) {
		unsigned long tmp;
		(void)pthread_join(threads[i], (void **)&tmp);
		counter += tmp;
	}

	printf("\n%u\n", counter);

	return 0;
}
#endif
