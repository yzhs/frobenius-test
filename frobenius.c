#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <gmp.h>

#include "helpers.h"
#include "small_primes.h"
#include "frobenius.h"

static mpz_t tmp0, tmp1, tmp2;

static mpz_t POLY(base), exponent;

unsigned long multiplications;

/*
 * Return g_x(x)*g_1(x) mod (n, x^2 - b*x - c) where g_x(x) = d*x + e and g_1(x) = g_x*x + g_1 in the return arguments res_x and
 * res_1, representing the polynomial res_x*x + res_1.
 */
static void mult_mod(POLY_ARGS(res), CONST_POLY_ARGS(f), CONST_POLY_ARGS(g), MODULUS_ARGS)
{
	/*
	 * If deg g_x = 1, the whole thing amounts to multiplying the coefficients of g_1 with a constant and reducing them
	 * modulo n.
	 */
	if (mpz_sgn(f_x) == 0) {
		mpz_mul(res_x, f_1, g_x);
		mpz_mul(res_1, f_1, g_1);
		mpz_mod(res_x, res_x, n);
		mpz_mod(res_1, res_1, n);

		multiplications += 2;

		return;
	}

	// res_x = (f_x*g_x*b + f_x*g_1 + f_1*g_x) % n
	mpz_mul(tmp2, f_x, g_x);
	mpz_mul(tmp0, tmp2, b);
	mpz_addmul(tmp0, f_x, g_1);
	mpz_addmul(tmp0, f_1, g_x);

	// res_1 = (f_x*g_x*c + f_1*g_1) % n
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, f_1, g_1);

	mpz_mod(res_x, tmp0, n);
	mpz_mod(res_1, tmp1, n);

	multiplications += 6;
}

static void square_mod(POLY_ARGS(res), CONST_POLY_ARGS(f), MODULUS_ARGS)
{
	if (mpz_sgn(f_x) == 0) {
		mpz_set_ui(res_x, 0);
		mpz_mul(res_1, f_1, f_1);
		mpz_mod(res_1, res_1, n);

		multiplications += 1;

		return;
	}

	// compute res_x = f_x^2*b + 2*f_x*f_1
	mpz_mul(tmp2, f_x, f_x);
	mpz_mul(tmp0, tmp2, b);
	mpz_mul(tmp1, f_x, f_1);
	mpz_add(tmp1, tmp1, tmp1);
	mpz_add(tmp0, tmp0, tmp1);

	// and res_1 = f_x^2*c + f_1^2
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, f_1, f_1);

	mpz_mod(res_x, tmp0, n);
	mpz_mod(res_1, tmp1, n);

	multiplications += 5;
}

static void powm(POLY_ARGS(res), CONST_POLY_ARGS(b), const mpz_t e, MODULUS_ARGS)
{

	// Copy all input parameters that will be changed in this function.
	mpz_set(base_x, b_x);
	mpz_set(base_1, b_1);
	mpz_set(exponent, e);

	// Initialize the return value.
	mpz_set_ui(res_x, 0);
	mpz_set_ui(res_1, 1);

	while (mpz_sgn(exponent) != 0) {
		if (mpz_odd_p(exponent))
			mult_mod(POLY(res), POLY(base), POLY(res), MODULUS);
		square_mod(POLY(base), POLY(base), MODULUS);
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

static void mult_x_mod(POLY_ARGS(res), CONST_POLY_ARGS(f), MODULUS_ARGS)
{
	// In case res_1 and f_x point to the same memory, we have to make a copy.
	mpz_set(tmp0, f_x);
	mpz_set(tmp1, f_1);

	mpz_mul(res_1, c, f_x);
	mpz_mod(res_1, res_1, n);

	mpz_mul(res_x, b, tmp0);
	mpz_add(res_x, res_x, tmp1);
	mpz_mod(res_x, res_x, n);
}

static void sigma(POLY_ARGS(res), CONST_POLY_ARGS(f), MODULUS_ARGS)
{
	mpz_sub(res_x, n, f_x);
	mpz_mul(res_1, f_x, b);
	mpz_add(res_1, res_1, f_1);
	mpz_mod(res_1, res_1, n);
}

#define ret(x) do { result = (x); goto exit; } while (0)

static Primality steps_3_4_5(const mpz_t n, const mpz_t b, const mpz_t c)
{
	mpz_t POLY(x), POLY(x_t), s, t, tmp, POLY(foo);
	bool n_is_1_mod_4;
	unsigned long r;
	Primality result = composite;

	mpz_inits(POLY(x), POLY(x_t), s, t, tmp, POLY(foo), NULL);

	mpz_set_ui(x_x, 1);
	// x_1 is initialized as 0 by mpz_inits

	/*
	 * (3) Compute x^((n+1)/2) mod (n, x^2-bx-c).  If x^((n+1)/2) not in ℤ/nℤ,
	 * declare n to be composite and stop.
	 */

	// According to Grantham, Theorem 3.4, we have to differentiate the
	// cases where n=1 mod 4 and n=3 mod 4
	// So we check whether n = 1 mod 2^2 (x_x is 1 at this point anyway).
	n_is_1_mod_4 = mpz_congruent_2exp_p(n, x_x, 2);
	if (n_is_1_mod_4)
		mpz_sub_ui(tmp, n, 1);
	else
		mpz_add_ui(tmp, n, 1);

	split(&r, s, tmp);
	mpz_fdiv_q_2exp(t, s, 1);  // t = (s-1)/2

	// Calculate x_t_x and x_t_1, such that (x_t_x*x+x_t_1) = x^t mod (n, x^2-bx-c).
	powm(POLY(x_t), POLY(x), t, MODULUS);

	// Calculate (x^t)^2 = x^(s-1)
	square_mod(POLY(foo), POLY(x_t), MODULUS);

	// Now compute x * x^(s-1) = x^s
	mult_x_mod(POLY(foo), POLY(foo), MODULUS);

	// We now have foo_x * x + foo_1 = x^s.  All we have to do, to
	// calculate x^(n-1)/2 or x^(n+1)/2, is to square this polynomial r-1
	// times.
	for (unsigned long i = 0; i < r-1; i++)
		square_mod(POLY(foo), POLY(foo), MODULUS);

	if (n_is_1_mod_4)
		// At this point, foo_x * x + foo_1 = x^(n-1)/2.  We need to
		// calculate x^(n+1)/2, so we multiply the result by x again.
		mult_x_mod(POLY(foo), POLY(foo), MODULUS);

	mpz_add_ui(tmp, n, 1);
	mpz_fdiv_q_2exp(tmp, tmp, 1);
	powm(POLY(foo), POLY(x), tmp, MODULUS);

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
	powm(POLY(foo), POLY(x), s, MODULUS);
	mpz_sub_ui(tmp, n, 1);

	if (mpz_sgn(foo_x) == 0 && mpz_cmp_ui(foo_1, 1) == 0)
		ret(probably_prime);

	for (unsigned long i = 0; i < r - 1; i++) {
		if (mpz_sgn(foo_x) == 0 && mpz_congruent_p(foo_1, tmp, n))
			ret(probably_prime);
		square_mod(POLY(foo), POLY(foo), MODULUS);
	}

exit:
	// cleanup
	mpz_clears(POLY(x), POLY(x_t), s, tmp, POLY(foo), NULL);
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
