#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <gmp.h>

#include "helpers.h"
#include "small_primes.h"
#include "frobenius.h"

#define MODULUS n, b, c
#define poly(name) name##_x, name##_1

static mpz_t tmp0, tmp1, tmp2;

static mpz_t poly(base), exponent;

unsigned long multiplications;

/*
 * Return f(x)*g(x) mod (n, x^2 - b*x - c) where f(x) = d*x + e and g(x) = f*x + g in the return arguments res_x and
 * res_1, representing the polynomial res_x*x + res_1.
 */
static void mult_mod(mpz_t res_x, mpz_t res_1,
                     const mpz_t d, const mpz_t e, const mpz_t f, const mpz_t g,
                     const mpz_t n, const mpz_t b, const mpz_t c)
{
	/*
	 * If deg f = 1, the whole thing amounts to multiplying the coefficients of g with a constant and reducing them
	 * modulo n.
	 */
	if (mpz_sgn(d) == 0) {
		mpz_mul(res_x, e, f);
		mpz_mul(res_1, e, g);
		mpz_mod(res_x, res_x, n);
		mpz_mod(res_1, res_1, n);

		multiplications += 2;

		return;
	}

	// res_x = (d*f*b + d*g + e*f) % n
	mpz_mul(tmp2, d, f);
	mpz_mul(tmp0, tmp2, b);
	mpz_addmul(tmp0, d, g);
	mpz_addmul(tmp0, e, f);

	// res_1 = (d*f*c + e*g) % n
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, e, g);

	mpz_mod(res_x, tmp0, n);
	mpz_mod(res_1, tmp1, n);

	multiplications += 6;
}

static void square_mod(mpz_t res_x, mpz_t res_1,
                       const mpz_t d, const mpz_t e,
                       const mpz_t n, const mpz_t b, const mpz_t c)
{
	if (mpz_sgn(d) == 0) {
		mpz_set_ui(res_x, 0);
		mpz_mul(res_1, e, e);
		mpz_mod(res_1, res_1, n);

		multiplications += 1;

		return;
	}

	// compute res_x = d^2*b + 2*d*e
	mpz_mul(tmp2, d, d);
	mpz_mul(tmp0, tmp2, b);
	mpz_mul(tmp1, d, e);
	mpz_add(tmp1, tmp1, tmp1);
	mpz_add(tmp0, tmp0, tmp1);

	// and res_1 = d^2*c + e^2
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, e, e);

	mpz_mod(res_x, tmp0, n);
	mpz_mod(res_1, tmp1, n);

	multiplications += 5;
}

static void powm(mpz_t res_x, mpz_t res_1,
                 const mpz_t b0, const mpz_t b1, const mpz_t e,
                 const mpz_t n, const mpz_t b, const mpz_t c)
{

	// Copy all input parameters that will be changed in this function.
	mpz_set(base_x, b0);
	mpz_set(base_1, b1);
	mpz_set(exponent, e);

	// Initialize the return value.
	mpz_set_ui(res_x, 0);
	mpz_set_ui(res_1, 1);

	while (mpz_sgn(exponent) != 0) {
		if (mpz_odd_p(exponent))
			mult_mod(poly(res), poly(base), poly(res), MODULUS);
		square_mod(poly(base), poly(base), MODULUS);
		mpz_fdiv_q_2exp(exponent, exponent, 1);
	}
}

static Primality steps_1_2(const mpz_t n)
{
#define tmp tmp0
	/*  (2) If n is a square, it can obviously not be prime. */
	if (mpz_perfect_square_p(n))
		return composite;

	/* Every number larger than 2^31 is certainly larger than B, whence the
	 * full list of small primes has to be used in trial division. */
	if (mpz_fits_sint_p(n)) {
		unsigned long sqrt;
		mpz_sqrt(tmp, n);
		sqrt = mpz_get_ui(tmp);

		// Start from prime_list[1] == 3 stead of prime_list[0] == 2.
		for (unsigned long i = 1; i < len(prime_list) && prime_list[i] <= sqrt; i++)
			if (mpz_divisible_ui_p(n, prime_list[i]))
				return composite;

		/* If n < B^2, we have either found a prime factor already or n itself
		 * is prime. */
		if (sqrt < B)
			return prime;
	} else {
		// Start from prime_list[1] == 3 stead of prime_list[0] == 2.
		for (unsigned long i = 1; i < len(prime_list); i++)
			if (mpz_divisible_ui_p(n, prime_list[i]))
				return composite;
	}

	return probably_prime;
#undef tmp
}

static void mult_x_mod(mpz_t res_x, mpz_t res_1, const mpz_t d, const mpz_t e, const mpz_t n, const mpz_t b, const mpz_t c)
{
	// In case res_1 and d point to the same memory, we have to make a copy.
	mpz_set(tmp0, d);

	mpz_mul(res_1, c, d);
	mpz_mod(res_1, res_1, n);

	mpz_mul(res_x, b, tmp0);
	mpz_add(res_x, res_x, e);
	mpz_mod(res_x, res_x, n);
}

static void invert(mpz_t res_x, mpz_t res_1, const mpz_t d, const mpz_t e, const mpz_t n, const mpz_t b, const mpz_t c)
{
	// (dx+e)^(-1) = (bde-cd^2+e)^(-1)(-dx + bde)
	mpz_t foo;
	mpz_init(foo);

	mpz_neg(res_x, d);

	mpz_mul(res_1, b, d);
	mpz_mul(res_1, res_1, e);
	mpz_mod(res_1, res_1, n);

	mpz_mul(foo, d, d);
	mpz_mul(foo, foo, c);
	mpz_sub(foo, e, foo);
	mpz_add(foo, foo, res_1);
	mpz_invert(foo, foo, n);

	mpz_mul(res_x, res_x, foo);
	mpz_mul(res_1, res_1, foo);

	mpz_mod(res_x, res_x, n);
	mpz_mod(res_1, res_1, n);

	mpz_clear(foo);
}

#define ret(x) do { result = (x); goto exit; } while (0)

static Primality steps_3_4_5(const mpz_t n, const mpz_t b, const mpz_t c)
{
	mpz_t poly(x), poly(x_t), s, t, tmp, poly(foo);
	unsigned long r;
	Primality result = composite;

	mpz_inits(poly(x), poly(x_t), s, t, tmp, poly(foo), NULL);

	mpz_set_ui(x_x, 1);
	// x_1 is initialized as 0 by mpz_inits

	/*
	 * (3) Compute x^((n+1)/2) mod (n, x^2-bx-c).  If x^((n+1)/2) not in ℤ/nℤ,
	 * declare n to be composite and stop.
	 */
	mpz_add_ui(tmp, n, 1);		// tmp = n+1
	mpz_fdiv_q_2exp(tmp, tmp, 1);   // tmp = (n+1)/2
	powm(foo0, foo1, x0, x1, tmp, n, b, c);
	powm(poly(foo), poly(x), tmp, MODULUS);

	/* check whether x^((n+1)/2) has degree 1 */
	if (mpz_sgn(foo_x) != 0)
		ret(composite);

	/*
	 * (4) Compute x^(n+1) mod (n, x^2-bx-c).  If x^(n+1) not congruent -c,
	 * declare n to be composite and stop.
	 */
	mpz_mul(foo_1, foo_1, foo_1);
	mpz_sub(tmp, n, c);
	if (!mpz_congruent_p(foo_1, tmp, n))
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
	powm(poly(foo), poly(x), s, MODULUS);
	mpz_sub_ui(tmp, n, 1);

	if (mpz_sgn(foo_x) == 0 && mpz_cmp_ui(foo_1, 1) == 0)
		ret(probably_prime);

	for (unsigned long i = 0; i < r - 1; i++) {
		if (mpz_sgn(foo_x) == 0 && mpz_congruent_p(foo_1, tmp, n))
			ret(probably_prime);
		square_mod(poly(foo), poly(foo), MODULUS);
	}

exit:
	// cleanup
	mpz_clears(poly(x), poly(x_t), s, tmp, poly(foo), NULL);
	return result;
}

/*
 * The Quadratic Frobenius Test (QFT) with parameters (b,c) consists of the
 * following.
 */
Primality QFT(const mpz_t n, const mpz_t b, const mpz_t c)
{
	Primality result = steps_1_2(n);

	if (result != probably_prime)
		return result;

	return steps_3_4_5(MODULUS);
}

#define check_non_trivial_divisor(num) do { \
		mpz_gcd(tmp, num, n); \
		if (mpz_cmp_ui(tmp, 1) != 0 && mpz_cmp(tmp, n) != 0) \
			ret(composite); \
} while (0)

Primality RQFT(const mpz_t n, const unsigned k)
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

	for (unsigned j = 0; j < k; j++) {
		do {
			get_random(c, n);
			mpz_mul(c, c, c);
			mpz_mod(c, c, n);
			mpz_sub(c, n, c);
		} while (mpz_cmp_ui(c, 3) < 0);
		check_non_trivial_divisor(c);
		j2 = 1;

		for (unsigned i = 0; i < B; i++) {
			get_random(b, n);

			mpz_mul(bb4c, b, b);
			mpz_addmul_ui(bb4c, c, 4);
			j1 = mpz_jacobi(bb4c, n);
			//mpz_sub(neg_c, n, c);
			//j2 = mpz_jacobi(neg_c, n);
			if (j1 == -1 && j2 == 1) {
				check_non_trivial_divisor(bb4c);
				check_non_trivial_divisor(b);
				break;
			}
		}
		if (j1 != -1 || j2 != 1) {
			gmp_printf("Found no suitable pair (b,c) modulo n=%Zd.  This is highly " \
					"unlikely unless the programme is wrong.  Assuming n is a prime...\n", n);
		} else {
			result = steps_3_4_5(MODULUS);
			if (result != probably_prime)
				ret(result);
		}
	}

exit:
	mpz_clears(b, c, nm1, bb4c, neg_c, tmp, NULL);
	return result;
}

#undef check_non_trivial_divisor
