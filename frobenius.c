#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <gmp.h>

#include "helpers.h"
#include "small_primes.h"
#include "frobenius.h"

uint64_t multiplications;

/*
 * Return f(x) * g(x) mod (n, x² - b*x - c) where f(x) = f_x*x + f_1 and
 * g(x) = g_x*x + g_1 in the return arguments res_x and res_1, representing the
 * polynomial res_x*x + res_1.
 */
static void mult_mod(POLY_ARGS(res), CONST_POLY_ARGS(f), CONST_POLY_ARGS(g), MODULUS_ARGS)
{
	// If deg g_x = 1, the whole thing amounts to multiplying the
	// coefficients of g_1 with a constant and reducing them modulo n.
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
	mpz_mod(res_x, tmp0, n);

	// res_1 = (f_x*g_x*c + f_1*g_1) % n
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, f_1, g_1);
	mpz_mod(res_1, tmp1, n);

	multiplications += 6;
}

/*
 * Compute the square of f, that is (f_x*x + f_1)² mod (n, x² - bx - c).
 */
static void square_mod(POLY_ARGS(res), CONST_POLY_ARGS(f), MODULUS_ARGS)
{
	if (mpz_sgn(f_x) == 0) {
		mpz_set_ui(res_x, 0);
		mpz_mul(res_1, f_1, f_1);
		mpz_mod(res_1, res_1, n);

		multiplications += 1;

		return;
	}

	// Compute res_x = f_x²*b + 2*f_x*f_1
	mpz_mul(tmp2, f_x, f_x);
	mpz_mul(tmp0, tmp2, b);
	mpz_mul(tmp1, f_x, f_1);
	mpz_add(tmp1, tmp1, tmp1);
	mpz_add(tmp0, tmp0, tmp1);
	mpz_mod(res_x, tmp0, n);

	// and res_1 = f_x²*c + f_1²
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, f_1, f_1);
	mpz_mod(res_1, tmp1, n);

	multiplications += 5;
}

/*
 * Compute b^exponent mod (n, x² - bx - c) where b is the polynomial b_x*x + b_1.
 */
static void powm(POLY_ARGS(res), CONST_POLY_ARGS(b), const mpz_t exponent, MODULUS_ARGS)
{

	// Copy all input parameters that will be changed in this function.
	mpz_set(base_x, b_x);
	mpz_set(base_1, b_1);

	// Initialize the return value.
	mpz_set_ui(res_x, 0);
	mpz_set_ui(res_1, 1);

	for (uint64_t k = mpz_sizeinbase(exponent, 2) - 1; k < (1lu << 63); k--) {
		square_mod(POLY(res), POLY(res), MODULUS);
		if (mpz_tstbit(exponent, k))
			mult_mod(POLY(res), POLY(base), POLY(res), MODULUS);
	}
}

/*
 * Compute x^exponent mod (n, x² - bx + c) using Lucas sequences.
 */
static void power_of_x(POLY_ARGS(res), const mpz_t exponent, MODULUS_ARGS)
{
	bool j_even = false; // We only need j to compute (-1)^j, so all we care about is whether j is odd or even.
	mpz_t A_j, B_j, C_j;
	mpz_inits(A_j, B_j, C_j, NULL);

	// Start with A_1 = x^1 + (b-x)^1 = b
	mpz_set(A_j, b);
	// and B_1 = (x^1 - (b - x)^1)/(2x-b) = 1.
	mpz_set_ui(B_j, 1);
	// Obviously C_1 = c.
	mpz_set(C_j, c);

	// Skip the leading 1 bit and convert to 0 based indexing
	for (uint64_t k = mpz_sizeinbase(exponent, 2) - 1 - 1; k < (1lu << 63); k--) {
		/*
		 * Doubling
		 */

		// Compute B_{2j}
		mpz_mul(B_j, B_j, A_j);
		mpz_mod(B_j, B_j, n);

		// Compute A_{2j}
		mpz_mul(A_j, A_j, A_j);
		mpz_add(tmp1, C_j, C_j);
		// TODO A conditional branch in a tight inner loop is a bad idea. Rewrite this!
		if (j_even)
			mpz_sub(A_j, A_j, tmp1);
		else
			mpz_add(A_j, A_j, tmp1);
		mpz_mod(A_j, A_j, n);

		// Compute C_{2j}
		mpz_mul(C_j, C_j, C_j);
		mpz_mod(C_j, C_j, n);

		j_even = true;
		multiplications += 3;

		if (mpz_tstbit(exponent, k)) {
			/*
			 * Chain addition
			 */

			// Compute A_{j+1}
			mpz_mul(tmp1, bb4c, B_j); // B_1 is 1, so we can just ignore it.
			mpz_addmul(tmp1, b, A_j); // A_1 = b
			mpz_mod(tmp1, tmp1, n);
			// If tmp1 is odd, tmp1+n is even.
			if (mpz_odd_p(tmp1))
				mpz_add(tmp1, tmp1, n);
			// Now tmp1 is even, so we can just do a right shift to
			// divide by 2.
			mpz_fdiv_q_2exp(tmp1, tmp1, 1);

			// Compute B_{j+1}
			mpz_mul(B_j, b, B_j);    // A_1 = b
			mpz_add(B_j, B_j, A_j);  // Use the old A_j, not A_{j+1}; B_1 = 1
			mpz_mod(B_j, B_j, n);
			if (mpz_odd_p(B_j))
				mpz_add(B_j, B_j, n);
			mpz_fdiv_q_2exp(B_j, B_j, 1);

			// Set the new A_j
			mpz_set(A_j, tmp1);

			// Compute C_{j+1}
			mpz_mul(C_j, C_j, c);
			mpz_mod(C_j, C_j, n);

			j_even = false;
			multiplications += 4;
		}
	}

	// Compute the polynomial x^j = res_x * x + res_1 from A_j and B_j.
	mpz_set(res_x, B_j);

	mpz_set(res_1, A_j);
	mpz_submul(res_1, b, B_j);
	if (mpz_odd_p(res_1))
		mpz_add(res_1, res_1, n);
	mpz_fdiv_q_2exp(res_1, res_1, 1);
	mpz_mod(res_1, res_1, n);

	multiplications++;

	mpz_clears(A_j, B_j, C_j, NULL);
}

/*
 * Perform the deterministic steps of the QFT, that is trial division and the
 * test whether n is a perfect square.
 */
Primality steps_1_2(const mpz_t n)
{
	/**********************************************************************\
	* Step (2)                                                             *
	\**********************************************************************/

	if (mpz_perfect_square_p(n))
		return composite;

	/**********************************************************************\
	* Step (1)                                                             *
	\**********************************************************************/

	// Every number larger than 2^31 is certainly larger than B, whence the
	// full list of small primes has to be used in trial division.
	if (mpz_fits_sint_p(n)) {
		uint64_t sqrt;
		mpz_sqrt(tmp0, n);
		sqrt = mpz_get_ui(tmp0);

		// Start from prime_list[1] == 3 stead of prime_list[0] == 2.
		for (uint64_t i = 1; i < len(prime_list) && prime_list[i] <= sqrt; i++)
			if (mpz_divisible_ui_p(n, prime_list[i]))
				return composite;

		// If n < B², we have either found a prime factor already or n
		// itself is prime.
		if (sqrt < B)
			return prime;
	} else {
		// Start from prime_list[1] == 3 instead of prime_list[0] == 2.
		for (uint64_t i = 1; i < len(prime_list); i++)
			if (mpz_divisible_ui_p(n, prime_list[i]))
				return composite;
	}

	// No factors found and n is large enougth that it might still factor
	// into multiple primes larger B.
	return probably_prime;
}

/*
 * Compute f*x = f_x * x² + f_1 * x = (b * f_x + f_1) * x + c * f_x for a
 * given polynomial f. Thus res_x = b * f_x + f_1 and res_1 = c * f_x.
 */
static void mult_x_mod(POLY_ARGS(res), CONST_POLY_ARGS(f), MODULUS_ARGS)
{
	// In case res_x and f_x or res_1 and f_1  point to the same memory, we
	// have to make a copy.
	mpz_set(tmp0, f_x);
	mpz_set(tmp1, f_1);

	mpz_mul(res_1, c, f_x);
	mpz_mod(res_1, res_1, n);

	mpz_mul(res_x, b, tmp0);
	mpz_add(res_x, res_x, tmp1);
	mpz_mod(res_x, res_x, n);

	multiplications += 2;
}

/*
 * Apply the homomorphism given by x |--> b - x.
 */
static void sigma(POLY_ARGS(res), CONST_POLY_ARGS(f), MODULUS_ARGS)
{
	mpz_set(res_1, f_1);
	mpz_addmul(res_1, f_x, b);
	mpz_mod(res_1, res_1, n);

	mpz_sub(res_x, n, f_x);

	multiplications += 1;
}

/*
 * Perform the non-deterministic steps of the Quadratic Frobenius Test.
 */
static Primality steps_3_4_5(MODULUS_ARGS)
{
	Primality result = composite;

	bool n_is_1_mod_4;

	// At first, 2^r*s = n ± 1 (sign depending on whether n is 1 mod 4).
	// At that time, r and s correspond to the variables $r'$ and $s'$ as
	// mentioned in the thesis.
	// Later, 2^r*s = n² - 1, at which point the r and s correspond to the
	// variables $r$ and $s$.
	uint64_t r;
	mpz_t s;

	// If 2^r*s = n ± 1, t = (s-1)/2.
	mpz_t t;

	mpz_t POLY(x);       // The polynomial x
	mpz_t POLY(x_t);     // The polynomial x^t reduced modulo (n, x² - bx - c)
	mpz_t POLY(x_n_1_2); // x^((n+1)/2) reduced mod (n, x² - bx -c)
	mpz_t POLY(foo);     // Temporary polynomial.  Used to store x^(n+1), x^s, ...

	// Allocate memory for long integers
	mpz_inits(POLY(x), POLY(x_t), POLY(x_n_1_2), POLY(foo), s, t, NULL);

	// x_1 is initialized as 0 by mpz_inits
	mpz_set_ui(x_x, 1);

	/**********************************************************************\
	* Step (3). Check whether -c is a square mod (n, x² - bx - c).         *
	*                                                                      *
	* The following calculations could be replaced by                      *
	*                                                                      *
	*       mpz_cdiv_q_2exp(tmp0, n, 1);  // tmp0 = ceil(n/2) = (n+1)/2    *
	*       power_of_x(POLY(foo), tmp0, MODULUS);                          *
	*                                                                      *
	* That would, however, mean that we will have to do more work for      *
	* step (5).  We can avoid this by precomputing x^t (and x^((n+1)/2)    *
	* if n is 1 mod 4).                                                    *
	\**********************************************************************/

	// According to Grantham, Theorem 3.4, we have to differentiate the
	// cases where n=1 mod 4 and n=3 mod 4
	// So we check whether n = 1 mod 2² (x_x is 1 at this point anyway).
	n_is_1_mod_4 = mpz_fdiv_ui(n, 4) == 1;
	if (n_is_1_mod_4)
		mpz_sub_ui(tmp0, n, 1);
	else
		mpz_add_ui(tmp0, n, 1);

	split(&r, s, tmp0);
	mpz_fdiv_q_2exp(t, s, 1);  // t = (s-1)/2

	// Calculate x_t_x and x_t_1, such that (x_t_x*x+x_t_1) = x^t mod (n, x²-bx-c).
	power_of_x(POLY(x_t), t, MODULUS);

	// Calculate (x^t)² = x^(s-1)
	square_mod(POLY(foo), POLY(x_t), MODULUS);

	// Now compute x * x^(s-1) = x^s
	mult_x_mod(POLY(foo), POLY(foo), MODULUS);

	// We now have foo_x * x + foo_1 = x^s.  All we have to do, to
	// calculate x^(n-1)/2 or x^(n+1)/2, is to square this polynomial r-1
	// times.
	for (uint64_t i = 0; i < r - 1; i++)
		square_mod(POLY(foo), POLY(foo), MODULUS);

	if (n_is_1_mod_4) {
		// At this point, foo_x * x + foo_1 = x^(n-1)/2.  We need to
		// calculate x^((n+1)/2), so we multiply the result by x again.
		mult_x_mod(POLY(foo), POLY(foo), MODULUS);

		// Store a copy for later use.
		mpz_set(x_n_1_2_x, foo_x);
		mpz_set(x_n_1_2_1, foo_1);
	}

	// Check whether x^((n+1)/2) has degree 1
	if (mpz_sgn(foo_x) != 0)
		ret(composite);

	/**********************************************************************\
	* Step (4). Check, whether x^n is -c mod (n, x² - bx - c).             *
	\**********************************************************************/
	mpz_mul(foo_1, foo_1, foo_1);
	multiplications += 1;
	mpz_sub(tmp0, n, c);
	if (!mpz_congruent_p(foo_1, tmp0, n))
		ret(composite);

	/**********************************************************************\
	* Step (5).  Use the precomputed values for computing x^s to reduce    *
	* run time.                                                            *
	\**********************************************************************/
	mpz_mul(tmp0, n, n);
	multiplications += 1;
	mpz_sub_ui(tmp0, tmp0, 1);

	// Calculate r,s such that 2^r*s == n² - 1.
	split(&r, s, tmp0);
	if (n_is_1_mod_4) {
		sigma(POLY(foo), POLY(x_t), MODULUS);
		mult_mod(POLY(foo), POLY(foo), POLY(x_t), MODULUS);
		mult_mod(POLY(foo), POLY(foo), POLY(x_n_1_2), MODULUS);
	} else {
		//sigma(POLY(foo), POLY(x_t), MODULUS);
		//mult_mod(POLY(foo), POLY(foo), POLY(x_t), MODULUS);
		powm(POLY(foo), POLY(x), s, MODULUS);
	}

	mpz_sub_ui(tmp0, n, 1);

	if (mpz_sgn(foo_x) == 0 && mpz_cmp_ui(foo_1, 1) == 0)
		ret(probably_prime);

	for (uint64_t i = 0; i < r - 1; i++) {
		if (mpz_sgn(foo_x) == 0 && mpz_cmp(foo_1, tmp0) == 0)
			ret(probably_prime);
		square_mod(POLY(foo), POLY(foo), MODULUS);
	}

exit:
	// Deallocate the local variables.
	mpz_clears(POLY(x), POLY(x_t), POLY(x_n_1_2), POLY(foo), s, t, NULL);
	return result;
}

/*
 * Execute the Quadratic Frobenius Test with parameters (b,c).
 *
 * Returns 'prime' if n is certainly prime, 'probably_prime' if no evidence
 * could be found that n might be composite and 'composite' otherwise.
 */
Primality QFT(MODULUS_ARGS)
{
	Primality result = steps_1_2(n);

	if (result != probably_prime)
		return result;

	return steps_3_4_5(MODULUS);
}

/*
 * Check whether gcd(n, num) is either n or 1.  Otherwise return composite from
 * RQFT.
 */
#define check_non_trivial_divisor(num) do { \
		mpz_gcd(tmp0, num, n); \
		if (mpz_cmp_ui(tmp0, 1) != 0 && mpz_cmp(tmp0, n) != 0) \
			ret(composite); \
} while (0)

/*
 * The randomized quadratic Frobenius Test (RQFT).
 */
Primality RQFT(const mpz_t n, const unsigned k)
{
	// Temporary storage for results from helper functions
	Primality result;

	// The Jacobi symbol (b²+4c/n)
	int j_bb4c = 0;

	// The pair of parameters
	mpz_t b, c;

	if (mpz_even_p(n)) {
		// 2 is the only odd prime...
		if (mpz_cmp_ui(n, 2) == 0)
			return prime;
		else
			return composite;
	}

	assert(mpz_cmp_ui(n, 1) > 0);

	mpz_inits(b, c, NULL);

	result = steps_1_2(n);
	// If the number is found to be either composite or certainly prime, we
	// can return that result immediately.
	if (result != probably_prime)
		ret(result);

	for (unsigned j = 0; j < k; j++) {
		do {
			get_random(c, n);
			mpz_mul(c, c, c);
			multiplications += 1;
			mpz_mod(c, c, n);
			mpz_sub(c, n, c);
		} while (mpz_cmp_ui(c, 3) < 0);
		// At this point, -c is a square by construction, so we don't
		// have to compute the Jacobi symbol (-c/n).
		check_non_trivial_divisor(c);

		for (unsigned i = 0; i < B; i++) {
			// Try choosing a `b` such that (b,c) is a valid pair
			// of parameters.
			get_random(b, n);

			mpz_mul(bb4c, b, b);
			mpz_addmul_ui(bb4c, c, 4);
			multiplications += 2;
			j_bb4c = mpz_jacobi(bb4c, n);
			if (j_bb4c == -1) {
				check_non_trivial_divisor(bb4c);
				check_non_trivial_divisor(b);
				break;
			}
		}
		if (j_bb4c != -1) {
			// This case is incredibly unlikely to ever get
			// executed.  The probability is provably less than
			// (3/4)^B < 10^(-5616).
			gmp_printf("Found no suitable pair (b,c) modulo n=%Zd.  This is highly " \
			           "unlikely unless the programme is wrong.  Assuming n is a prime...\n", n);
		} else {
			// If we did find a valid pair (b,c), execute the
			// non-deterministic steps (3), (4) and (5) with these
			// parameters (b,c).
			result = steps_3_4_5(MODULUS);

			// If the result is `composite`, stop immediately,
			// otherwise we might have to do another iteration.
			if (result != probably_prime)
				ret(result);
		}
	}

exit:
	mpz_clears(b, c, NULL);
	return result;
}

#undef check_non_trivial_divisor
