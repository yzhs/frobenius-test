#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <gmp.h>

#include "helpers.h"
#include "small_primes.h"

Primality QFT(mpz_t n, mpz_t b, mpz_t c);
Primality RQFT(mpz_t n, unsigned B);

static mpz_t tmp0, tmp1, tmp2;

/*
 * Return f(x)*g(x) mod (n, x^2 - b*x - c) where f(x) = d*x + e and g(x) = f*x + g in the return arguments res0 and
 * res1, representing the polynomial res0*x + res1.
 */
static void mult_mod(mpz_t res0, mpz_t res1, mpz_t d, mpz_t e, mpz_t f, mpz_t g, mpz_t n, mpz_t b, mpz_t c)
{
	/*
	 * If deg f = 1, the whole thing amounts to multiplying the coefficients of g with a constant and reducing them
	 * modulo n.
	 */
	if (mpz_sgn(d) == 0) {
		mpz_mul(res0, e, f);
		mpz_mod(res0, res0, n);
		return;
	}

	// res0 = (d*f*b + d*g + e*f) % n
	mpz_mul(tmp2, d, f);
	mpz_mul(tmp0, tmp2, b);
	mpz_addmul(tmp0, d, g);
	mpz_addmul(tmp0, e, f);

	// res1 = (d*f*c + e*g) % n
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, e, g);

	mpz_mod(res0, tmp0, n);
	mpz_mod(res1, tmp1, n);
}

static void square_mod(mpz_t res0, mpz_t res1, mpz_t d, mpz_t e, mpz_t n, mpz_t b, mpz_t c)
{
	if (mpz_sgn(d) == 0) {
		mpz_set_ui(res0, 0);
		mpz_mul(res1, e, e);
		mpz_mod(res1, res1, n);
		return;
	}

	// compute res0 = d^2*b + 2*d*e
	mpz_mul(tmp2, d, d);
	mpz_mul(tmp0, tmp2, b);
	mpz_mul(tmp1, d, e);
	mpz_addmul_ui(tmp0, tmp1, 2);

	// and res1 = d^2*c + e^2
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, e, e);

	mpz_mod(res0, tmp0, n);
	mpz_mod(res1, tmp1, n);
}

static void pomw(mpz_t res0, mpz_t res1, mpz_t base0, mpz_t base1, mpz_t exp, mpz_t n, mpz_t b, mpz_t c)
{
	mpz_set_ui(res0, 0);
	mpz_set_ui(res1, 1);

	while (mpz_sgn(exp) != 0) {
		if (mpz_odd_p(exp))
			mult_mod(res0, res1, base0, base1, res0, res1, n, b, c);
		square_mod(base0, base1, base0, base1, n, b, c);
		mpz_fdiv_q_ui(exp, exp, 2);
	}
}

static Primality steps_1_2(mpz_t n)
{
#define tmp tmp0
	/*  (2) If n is a square, it can obviously not be prime. */
	if (mpz_perfect_square_p(n))
		return composite;

	/* Every number larger than 2^31 is certainly larger then B = 50000, whence
	 * the full list of small primes has to be used in trial division. */
	if (mpz_fits_sint_p(n)) {
		unsigned long sqrt;
		mpz_sqrt(tmp, n);
		sqrt = mpz_get_ui(tmp);

		for (unsigned long i = 0; i < len(prime_list) && sqrt < prime_list[i]; i++)
			if (mpz_divisible_ui_p(n, prime_list[i]))
				return composite;

		/* If n < B^2, we have either found a prime factor already or n itself
		 * is prime. */
		if (sqrt < 50000)
			return prime;
	} else {
		for (unsigned long i = 0; i < len(prime_list); i++)
			if (mpz_divisible_ui_p(n, prime_list[i]))
				return composite;
	}

	return probably_prime;
#undef tmp
}

#define ret(x) do { result = (x); goto exit; } while (0)
static Primality steps_3_4_5(mpz_t n, mpz_t b, mpz_t c)
{
	mpz_t x0, x1, s, tmp, foo0, foo1;
	unsigned long r;
	Primality result = composite;

	mpz_inits(x0, x1, s, tmp, foo0, foo1, NULL);

	mpz_set_ui(x0, 1);
	// x1 is initialized as 0 by mpz_inits

	/*
	 * (3) Compute x^((n+1)/2) mod (n, x^2-bx-c).  If x^((n+1)/2) not in ℤ/nℤ,
	 * declare n to be composite and stop.
	 */
	mpz_add_ui(tmp, n, 1);          // tmp = n+1
	mpz_fdiv_q_2exp(tmp, tmp, 1);   // tmp = (n+1)/2
	pomw(foo0, foo1, x0, x1, tmp, n, b, c);

	/* check whether x^((n+1)/2) has degree 1 */
	if (mpz_sgn(foo0) != 0)
		ret(composite);

	/*
	 * (4) Compute x^(n+1) mod (n, x^2-bx-c).  If x^(n+1) not congruent -c,
	 * declare n to be composite and stop.
	 */
	mpz_mul(foo1, foo1, foo1);
	mpz_sub(tmp, n, c);
	if (!mpz_congruent_p(foo1, tmp, n))
		ret(composite);

	/*
	 * (5) Let n^2-1=2^r*s, where s is odd.  If x^s not congruent 1 mod (n,
	 * x^2-bx-c), and x^(2^j*s) not congruent -1 mod (n, x^2-bx-c) for all
	 * 0≤j≤r-2, declare n to be composite and stop.
	 * If n is not declared composite in Steps 1—5, declare n to be a probable
	 * prime.
	 */
	mpz_mul(tmp, n, n);
	/* calculate r,s such that 2^r*s + 1 == n^2 */
	split(&r, s, tmp);
	pomw(foo0, foo1, x0, x1, s, n, b, c);
	mpz_sub_ui(tmp, n, 1);

	if (mpz_sgn(foo0) == 0 && mpz_cmp_ui(foo1, 1) == 0)
		ret(probably_prime);

	for (unsigned long i = 0; i < r - 1; i++) {
		if (mpz_sgn(foo0) == 0 && mpz_congruent_p(foo1, tmp, n))
			ret(probably_prime);
		square_mod(foo0, foo1, foo0, foo1, n, b, c);
	}

exit:
	// cleanup
	mpz_clears(x0, x1, s, tmp, foo0, foo1, NULL);
	return result;
}

/*
 * The Quadratic Frobenius Test (QFT) with parameters (b,c) consists of the
 * following.
 */
Primality QFT(mpz_t n, mpz_t b, mpz_t c)
{
	Primality result = steps_1_2(n);

	if (result != probably_prime)
		return result;

	return steps_3_4_5(n, b, c);
}

Primality RQFT(mpz_t n, unsigned B)
{
	mpz_t b, c, nm1;
	mpz_t bb4c, neg_c, tmp;
	int j1 = 0, j2 = 0;
	Primality result;

	if (mpz_even_p(n)) {
		/*  2 is the only odd prime... */
		if (mpz_cmp_ui(n, 2) == 0)
			return prime;
		else
			return composite;
	}

	assert(mpz_cmp_ui(n, 1) > 0);

	mpz_inits(b, c, nm1, bb4c, neg_c, tmp, NULL);
	mpz_sub_ui(nm1, n, 1);

	result = steps_1_2(n);
	/* If the number is found to be either composite or certainly prime, we can
	 * return that result immediately. */
	if (result != probably_prime)
		ret(result);

#ifdef check_non_trivial_divisor
#undef check_non_trivial_divisor
#endif
#define check_non_trivial_divisor(num) do { \
		mpz_gcd(tmp, num, n); \
		if (mpz_cmp_ui(tmp, 1) != 0 && mpz_cmp(tmp, n) != 0) \
			ret(composite); \
} while (0)

	for (unsigned i = 0; i < B; i++) {
		get_random(b, n);
		get_random(c, n);

		mpz_mul(bb4c, b, b);
		mpz_addmul_ui(bb4c, c, 4);
		j1 = mpz_jacobi(bb4c, n);
		mpz_neg(neg_c, c);
		j2 = mpz_jacobi(neg_c, n);
		if (j1 == -1 && j2 == 1) {
			check_non_trivial_divisor(bb4c);
			check_non_trivial_divisor(b);
			check_non_trivial_divisor(c);
			break;
		}
	}
	if (j1 != -1 || j2 != 1) {
		gmp_printf("Found no suitable pair (b,c) modulo n=%Zd.  This is highly " \
			   "unlikely unless the programme is wrong.  Assuming n is a prime...\n", n);
		ret(probably_prime);
	}

	result = steps_3_4_5(n, b, c);
exit:
	mpz_clears(b, c, nm1, bb4c, neg_c, tmp, NULL);
	return result;
}

#ifndef TEST
int main()
{
	mpz_t tmp;

	mpz_inits(tmp, tmp0, tmp1, tmp2, NULL);
	init();

	mpz_set_ui(tmp, 1215239);
	assert(RQFT(tmp, 1) != composite);
	mpz_set_ui(tmp, 1215237);
	assert(RQFT(tmp, 1) == composite);
	mpz_set_str(tmp, "2147483659", 10);
	assert(RQFT(tmp, 1) == probably_prime);
	mpz_set_str(tmp, "32317006071311007300714876688669951960444102669715484032130"\
			"34542752465513886789089319720141152291346368871796092189801949411955"\
			"91504909210950881523864482831206308773673009960917501977503896521067"\
			"96057638384067568276792218642619756161838094338476170470581645852036"\
			"30504288757589154106580860755239912393038552191433338966834242068497"\
			"47865645694948561760353263220580778056593310261927084603141502585928"\
			"64177116725943603718461857357598351152301645904403697613233287231227"\
			"12568471082020972515710172693132346967854258065669793504599726835299"\
			"86382155251663894373355436021354332296046453184786049521481935558536"\
			"11059596231637", 10);
	assert(RQFT(tmp, 1) == probably_prime);

	cleanup();
	mpz_clears(tmp, tmp0, tmp1, tmp2, NULL);
}
#endif
