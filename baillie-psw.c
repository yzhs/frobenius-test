#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include <gmp.h>

#include "helpers.h"
#include "baillie-psw.h"
#include "frobenius.h"
#include "miller_rabin.h"

/*
 * Compute (U_{2k}, V_{2k}) from (U_k, V_k).
 */
static void doubling(mpz_t U, mpz_t V, mpz_t Q,
                     const mpz_t U_k, const mpz_t V_k, const mpz_t Q_k,
		     const mpz_t n)
{

	mpz_mul(tmp1, V_k, V_k);
	mpz_mod(tmp1, tmp1, n);
	mpz_add(tmp0, Q_k, Q_k);
	mpz_sub(tmp1, tmp1, tmp0);

	mpz_mul(tmp0, U_k, V_k);

	mpz_mod(U, tmp0, n);
	mpz_mod(V, tmp1, n);

	mpz_mul(Q, Q_k, Q_k);
	mpz_mod(Q, Q, n);
}

/*
 * Compute (U_{k+1}, V_{k+1}) from (U_k, V_k).
 */
static void chain_addition(mpz_t U, mpz_t V, mpz_t Q_new,
                           const mpz_t U_k, const mpz_t V_k, const mpz_t Q_k,
			   const mpz_t n, const mpz_t D, const mpz_t Q)
{
	mpz_add(tmp0, U_k, V_k);

	mpz_mul(tmp1, D, U_k);
	mpz_add(tmp1, tmp1, V_k);

	// Since n is odd, either tmp0 or tmp0+n is even, so we can divide by
	// two using a right shift.
	if (mpz_odd_p(tmp0))
		mpz_add(tmp0, tmp0, n);
	mpz_fdiv_q_2exp(tmp0, tmp0, 1);
	mpz_mod(U, tmp0, n);

	// The same goes for tmp1.
	if (mpz_odd_p(tmp1))
		mpz_add(tmp1, tmp1, n);
	mpz_fdiv_q_2exp(tmp1, tmp1, 1);
	mpz_mod(V, tmp1, n);

	mpz_mul(tmp0, Q_k, Q);
	mpz_mod(Q_new, tmp0, n);
}

static Primality lucas_test(const mpz_t n, const mpz_t D, const mpz_t Q)
{
	Primality result = probably_prime;
#define VARIABLES delta, U_1, V_1, U_k, V_k, Q_k
	mpz_t VARIABLES;

	mpz_inits(VARIABLES, NULL);

	// P=1, so P and n are relatively prime.  Since the Jacobi symbol
	// (D/n)=-1, D and n cannot share a prime factor; but Q and n might.
	mpz_gcd(tmp0, n, Q);
	if (mpz_cmp(tmp0, n) != 0 && mpz_cmp_ui(tmp0, 1) != 0) 
		ret(composite);

	// We force the Jacobi symbol (D/n) to be -1, so we might as well write n+1 instead of n-(D/n).
	mpz_add_ui(delta, n, 1);

	mpz_set_ui(U_1, 1);

	mpz_set_ui(V_1, 1);

	// Left-to-right binary
	for (uint64_t k = mpz_sizeinbase(delta, 2) - 1 - 1; k < (1lu << 63); k--) {
		doubling(U_k, V_k, Q_k, U_k, V_k, Q_k, n);
		if (mpz_tstbit(delta, k))
			chain_addition(U_k, V_k, Q_k, U_k, V_k, Q_k, n, D, Q);
	}
	if (mpz_sgn(U_k) != 0)
		result = composite;

exit:
	mpz_clears(VARIABLES, NULL);
	return result;
#undef VARIABLES
}

/*
 * Perform the Baillie-PSW primality test.
 *
 * Based on the algorithm described at
 * http://en.wikipedia.org/wiki/Baillie–PSW_primality_test
 */
Primality Baillie_PSW(const mpz_t n)
{
#define VARIABLES Q, D
	mpz_t VARIABLES;
	long i, sign;
	// Perform trial division up to B and check whether n is a square.
	Primality result = steps_1_2(n);
	if (result != probably_prime)
		return result;

	mpz_inits(VARIABLES, tmp0, tmp1, NULL);

	mpz_set_ui(tmp0, 2);
	result = miller_rabin_base(n, tmp0);
	if (result != probably_prime)
		ret(result);

	/*
	 * Perform the Lucas test with the parameters chosen as described in
	 * the Wikipedia article.
	 */
	// The upper bound is used to make sure that the test will, at some
	// point, terminate.
	sign = 1;
	for (i = 5; i < 0xffffff; i+=2) {
		mpz_set_si(D, sign*i);
		sign = -sign;
		if (mpz_jacobi(D, n) == -1)
			break;
	}
	if (i == 0xffffff)
		return probably_prime;

	mpz_sub_ui(tmp0, n, 1);
	mpz_sub(tmp0, D, tmp0);
	mpz_fdiv_q_2exp(Q, tmp0, 2);

	result = lucas_test(n, D, Q);

exit:
	mpz_clears(VARIABLES, tmp0, tmp1, NULL);
	return result;
#undef VARIABLES
}
