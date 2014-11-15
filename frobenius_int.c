/* frobenius_int.c -- the quadratic frobenius test (long long version)
 *
 * Copyright 2014 by Colin Benner <colin-software@yzhs.de>
 *
 * This file is part of frobenius-test.
 *
 * frobenius-test is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * frobenius-test is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with frobenius-test.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

#include "helpers_int.h"
#include "small_primes.h"
#include "frobenius_int.h"

// Compute b^2+4c once for every n.
static uint64_t bb4c_int;

// Count the number of multiplications mod n performed.
static uint64_t multiplications_int;

#ifdef TEST
/*
 * Return f(x)*g(x) mod (n, x^2 - b*x - c) where f(x) = f_x*x + f_1 and
 * g(x) = * g_x*x + g_1 in the return arguments res_x and res_1, representing
 * the polynomial res_x*x + res_1.
 */
static void mult_mod_int(POLY_ARGS_int(res), CONST_POLY_ARGS_int(f), CONST_POLY_ARGS_int(g), MODULUS_ARGS_int)
{
	uint64_t fxgx, f1gx = mul(f_1, g_x), f1g1 = mul(f_1, g_1);
	multiplications_int += 2;

	// If deg f = 1, the whole thing amounts to multiplying the
	// coefficients of f with a constant and reducing them modulo n.
	if (f_x == 0) {
		*res_x = f1gx;
		*res_1 = f1g1;

		return;
	}
	fxgx = mul(f_x, g_x);

	*res_x = (mul(fxgx, b) + mul(f_x, g_1) + f1gx) % n;
	*res_1 = mul(fxgx, c) + f1g1;
	if (*res_1 > n)
		*res_1 -= n;

	multiplications_int += 4;
}
#endif

/*
 * Calculate (dx + e)^2 mod (n, x^2-bx-c) returning the result as (*res_x) * x + (*res_1).
 */
static void square_mod_int(POLY_ARGS_int(res), CONST_POLY_ARGS_int(f), MODULUS_ARGS_int)
{
	uint64_t f1f1 = mul(f_1, f_1);
	uint64_t fxfx;

	multiplications_int++;

	if (f_x == 0) {
		*res_x = 0;
		*res_1 = f1f1;
		return;
	}
	fxfx = mul(f_x, f_x);

	// compute res_x = f_x^2*b+2*f_x*f_1
	*res_x = (mul(fxfx, b) + 2 * mul(f_x, f_1)) % n;

	// and res_1 = f_x^2*c+f_1^2
	*res_1 = (mul(fxfx, c) + f1f1) % n;

	multiplications_int += 4;
}

#ifdef TEST
/*
 * Calculate (base_x * x + base_1)^exp mod (n, x^2-bx-c) returning the result as
 * (*res_x) * x + (*res_1).  The computation is done using exponentiation by
 * squaring.
 *
 * This function is only used to insure that power_of_x_int does the right thing.
 */
static void powm_int(POLY_ARGS_int(res), CONST_POLY_ARGS_int(b), uint64_t exp, MODULUS_ARGS_int)
{
	uint64_t POLY_int(base);
	base_x = b_x;
	base_1 = b_1;

	*res_x = 0;
	*res_1 = 1;

	while (exp != 0) {
		if (odd(exp))
			mult_mod_int(POLY_int(res), POLY_int(base), *res_x, *res_1, MODULUS_int);
		square_mod_int(&base_x, &base_1, POLY_int(base), MODULUS_int);
		exp /= 2;
	}
}
#endif

/*
 * Compute x^exponent mod (n, x² - bx + c) (basically using Lucas sequences) as
 * described in section 3.4 of the paper.
 */
static void power_of_x_int(POLY_ARGS_int(res), const uint64_t exponent, MODULUS_ARGS_int)
{
	int j_even = false; // We only need j to compute (-1)^j, so all we care about is whether j is odd or even.
	uint64_t A_j = b, B_j = 1, C_j = c, tmp0;

	// The following thre values are needed for the chain addition steps.
	// Values that are already stored in variables available locally, are
	// #defined, the remaining values is stored in a global temporary
	// variable.
#define A_1 b
	uint64_t B_1 = 1;
#define C_1 c

	int k;

	// Skip the leading 1 bit and convert to 0 based indexing
	for (k = 8*sizeof(uint64_t) - 1; k >= 0 && ((1lu << k) & exponent) == 0; k--);

	for (; k >= 0; k--) {
		/*
		 * Doubling
		 */

		B_j = mul(A_j, B_j);

		tmp0 = mul(2, C_j);

		A_j = mul(A_j, A_j);
		C_j = mul(C_j, C_j);

		if (j_even) // TODO A conditional branch in a tight inner loop is a bad idea. Rewrite this!
			tmp0 = -tmp0;

		A_j = (A_j + tmp0) % n;

		j_even = true;
		multiplications_int += 3;

		if (exponent & (1lu << k)) {
			/*
			 * Chain addition
			 */

			// Compute A_{j+1}
			tmp0 = mul(mul(bb4c_int, B_1), B_j) + mul(A_1, A_j);
			// Divide by 2
			if (tmp0 % 2 == 1)
				tmp0 += n;
			tmp0 >>= 1; // At this point tmp0 is certainly below n.

			// Compute B_{j+1}
			B_j = (mul(A_1, B_j) + mul(A_j, B_1)) % n;
			// Divide by 2
			if (B_j % 2 == 1)
				B_j += n;
			B_j >>= 1;  // Now B_j < n.

			// Set the new A_j
			A_j = tmp0;

			// Compute C_{j+1}
			C_j = mul(C_j, C_1);

			j_even = false;
			multiplications_int += 8;
		}
	}

	// Compute the polynomial x^j = res_x * x + res_1 from A_j and B_j.
	*res_x = B_j;

	*res_1 = mul(b, B_j);
	if (*res_1 < A_j)
		*res_1 = A_j - *res_1;
	else
		*res_1 = n + A_j - *res_1;
	// Divide *res_1 by 2
	if (*res_1 % 2 == 1)
		*res_1 += n;
	*res_1 >>= 1;
#undef B_1
}



/*
 * Like QFT, return 'probably_prime' if n might be prime, 'prime' if n is
 * certainly prime and 'composite' if a proof for n's compositeness was found.
 */
static Primality steps_1_2_int(const uint64_t n)
{
	uint64_t sqrt = int_sqrt(n);

	/*
	 * Step (2) If n is a square, it can obviously not be prime.
	 */
	if (sqrt * sqrt == n)
		return composite;

	/*
	 * Step (1) Trial division.
	 */
	// Start from prime_list[1] = 3 stead of prime_list[0] = 2.
	for (uint64_t i = 1; i < len(prime_list) && prime_list[i] <= sqrt; i++)
		if (n % prime_list[i] == 0)
			return composite;

	// If the given number is small enough, there cannot be a non-trivial
	// divisor of n, whence n is certainly prime.
	return (sqrt < B) ? prime : probably_prime;
}

/*
 * Execute steps (3) through (5) of the Quadratic Frobenius Test.
 */
Primality steps_3_4_5_int(MODULUS_ARGS_int)
{
	uint64_t POLY_int(foo), s, tmp;
	uint64_t r, i;

	s = tmp = foo_x = foo_1 = 0;

	/*
	 * Step (3) Check, whether -c is a square modulo n.
	 */
	split_int(&r, &s, n+1);
	tmp = n + 1;
	tmp = tmp / 2;
	power_of_x_int(&foo_x, &foo_1, tmp, MODULUS_int);
	if (foo_x != 0)  // Check whether x^((n+1)/2) has degree 1.
		return composite;

	/*
	 * Step (4) Check, whether x^(n+1) = -c mod (n, x^2-bx-c).
	 */
	foo_1 = mul(foo_1, foo_1);
	tmp = n - c;
	if (foo_1 != tmp)
		return composite;

	/*
	 * Step (5)
	 */
	split_int(&r, &s, n * n); // Calculate r,s such that 2^r*s + 1 == n^2
	power_of_x_int(&foo_x, &foo_1, s, MODULUS_int); // TODO make this more efficient using σ : x ↦ b-x.

	if (foo_x == 0 && foo_1 == 1)
		return probably_prime;

	tmp = n - 1;
	for (i = 0; i < r - 1; i++) {
		if (foo_x == 0 && foo_1 % n == tmp)
			return probably_prime;

		square_mod_int(&foo_x, &foo_1, foo_x, foo_1, MODULUS_int);
	}

	return composite;
}

/*
 * Execute the Quadratic Frobenius Test with parameters (b,c).
 *
 * Returns 'prime' if n is certainly prime, 'probably_prime' if no evidence
 * could be found that n might be composite and 'composite' otherwise.
 */
Primality QFT_int(const unsigned n_, const unsigned b_, const unsigned c_)
{
	const uint64_t n = n_, b = b_, c = c_;
	// Running step (2) is not really necessary, as (b,c) is a valid
	// parameter pair, so n cannot possibly be a perfect square.
	Primality result = steps_1_2_int(n);

	if (result != probably_prime)
		return result;

	return steps_3_4_5_int(MODULUS_int);
}

/*
 * Check whether gcd(n, num) is either n or 1.  Otherwise return composite from
 * RQFT.
 */
#define check_non_trivial_divisor(num) do { \
		tmp = gcd(num, n); \
		if (tmp != 1 && tmp != n) \
			return composite; \
} while (0)

/*
 * Run the randomized quadratic frobenius test to check whether [n] is a a
 * prime.  The Parameter [k] determines how many times the test will be run at
 * most.  If the test returns "composite", it will not be run again.
 */
Primality RQFT_int(const unsigned n_, const unsigned k)
{
	uint64_t n = n_;
	Primality result;
	uint64_t b = 0, c = 0;
	uint64_t tmp; // This is used to store the greatest commond divisors we calculate.

	// These variables store the values of the Jacobi symbols (b²+4c/n) and
	// (-c/n) respectivelyg
	int j_bb4c = 0, j_c = 0;

	if (n < 3)
		return n == 2 ? prime : composite;
	if (even(n))
		return composite;

	result = steps_1_2_int(n);
	// If the number is found to be either composite or certainly prime, we
	// can return that result immediately.
	if (result != probably_prime)
		return result;

	for (unsigned j = 0; j < k; j++) {
		for (unsigned i = 0; i < B; i++) {
			b = get_random_int(2, n - 2);
			c = get_random_int(2, n - 2);

			bb4c_int = (mul(b, b) + c * 4) % n;

			j_bb4c = jacobi(bb4c_int, n);
			j_c = jacobi(n - c, n); // NOTE: n-c is not congruent to -c, since -c is understood as 2⁶⁴-c !!!

			if (j_bb4c == -1 && j_c == 1) {
				check_non_trivial_divisor(bb4c_int);
				check_non_trivial_divisor(b);
				check_non_trivial_divisor(c);
				break;
			}
		}

		if (j_bb4c != -1 || j_c != 1) {
			printf("Found no suitable pair (b,c) modulo n=%lu. "
			       "This is highly unlikely unless the programme is wrong. "
			       "Assuming %lu is a prime...\n", n, n);
			return probably_prime;
		}

		result = steps_3_4_5_int(MODULUS_int);
		if (result == composite)
			return composite;
	}

	return probably_prime;
}

#undef check_non_trivial_divisor
