#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include "helpers.h"

#ifdef DEBUG
#define debug(...) printf(__VA_ARGS__)
#else
#define debug(...)
#endif

#define len(a) (sizeof(a)/sizeof(a[0]))

#include "small_primes.c"

mpz_t tmp0, tmp1, tmp2;

/*
 * Return f(x)*g(x) mod (n, x^2 - b*x - c) where f(x) = d*x + e and g(x) = f*x + g in the return arguments res0 and
 * res1, representing the polynomial res0*x + res1.
 */
void mult_mod(mpz_t res0, mpz_t res1, mpz_t d, mpz_t e, mpz_t f, mpz_t g, mpz_t n, mpz_t b, mpz_t c)
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

	// res0 = (d * f * b + d * g + e * f) % n
	mpz_mul(tmp2, d, f);
	mpz_mul(tmp0, tmp2, b);
	mpz_addmul(tmp0, d, g);
	mpz_addmul(tmp0, e, f);

	// res1 = (d * f * c + e * g) % n
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, e, g);

	mpz_mod(res0, tmp0, n);
	mpz_mod(res1, tmp1, n);
}

void square_mod(mpz_t res0, mpz_t res1, mpz_t d, mpz_t e, mpz_t n, mpz_t b, mpz_t c)
{
	if (mpz_sgn(d) == 0) {
		mpz_set_ui(res0, 0);
		mpz_mul(res1, e, e);
		mpz_mod(res1, res1, n);
		return;
	}

	// compute res0 = d^2*b+2*d*e
	mpz_mul(tmp2, d, d);
	mpz_mul(tmp0, tmp2, b);
	mpz_mul(tmp1, d, e);
	mpz_addmul_ui(tmp0, tmp1, 2);

	// and res1 = d^2*c+e^2
	mpz_mul(tmp1, tmp2, c);
	mpz_addmul(tmp1, e, e);

	mpz_mod(res0, tmp0, n);
	mpz_mod(res1, tmp1, n);
}

void pow_mod(mpz_t res0, mpz_t res1, mpz_t base0, mpz_t base1, mpz_t exp, mpz_t n, mpz_t b, mpz_t c)
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

int no_nontrivial_small_prime_divisor(mpz_t n, unsigned B)
{
#define tmp tmp0
#define sqrt tmp1

	if (mpz_cmp_ui(n, B) < 0) {
		unsigned m = mpz_get_ui(n);
		int low = 0, high = len(prime_list);

		while (low + 1 < high) {
			int i = (low + high) / 2;
			if (prime_list[i] == m)
				return -1;
			else if (prime_list[i] > m)
				high = i;
			else
				low = i;
		}
	}

	mpz_sqrt(sqrt, n);

	for (int i = 0; i < len(prime_list) && mpz_cmp_ui(sqrt, prime_list[i]) < 0; i++)
		if (mpz_divisible_ui_p(n, prime_list[i]))
			return 0;

	return 1;
#undef tmp
#undef sqrt
}

/*
 * The Quadratic Frobenius Test (QFT) with parameters (b,c) consists of the
 * following.
 */
int QFT(mpz_t n, mpz_t b, mpz_t c, unsigned B, int use_rqft)
{
#define ret(x) do { result = (x); goto exit; } while (0)
	mpz_t x0, x1, s, tmp, foo0, foo1;
	unsigned long r;
	int result = 0;
	mpz_inits(x0, x1, s, tmp, foo0, foo1, NULL);

	// Suppose n>1 is odd, (b^2+4c over n)=-1 and (-c over n)=1.

	if (!use_rqft) {
		/*
		 * (1) Test n for divisibility by primes less than or equal to
		 * min{B, sqrt(n)}.  If it is divisible by one of these primes,
		 * declare n to be composite and stop.
		 */
		int foo = no_nontrivial_small_prime_divisor(n, B);
		if (!foo)
			ret(0); // composite
		else if (foo == -1)
			ret(1); // n occurs in small_primes and is therefore certainly prime

		/*
		 * (2) Test whether sqrt(n) in ℤ.  If it is, declare n to be composite and
		 * stop.
		 */
		if (mpz_perfect_square_p(n)) {
			ret(0); // composite
		}
	}

	mpz_set_ui(x0, 1);
	// x1 is 0 automatically

	/*
	 * (3) Compute x^((n+1)/2) mod (n, x^2-bx-c).  If x^((n+1)/2) not in ℤ/nℤ,
	 * declare n to be composite and stop.
	 */
	mpz_add_ui(tmp, n, 1); // tmp = n+1
	mpz_fdiv_q_2exp(tmp, tmp, 1); // tmp = (n+1)/2
	pow_mod(foo0, foo1, x0, x1, tmp, n, b, c);
	if (mpz_sgn(foo0) != 0) { // check whether x^((n+1)/2) has degree 1
		ret(0); // composite
	}

	/*
	 * (4) Compute x^(n+1) mod (n, x^2-bx-c).  If x^(n+1) not congruent -c,
	 * declare n to be composite and stop.
	 */
	mpz_mul(foo1, foo1, foo1);
	mpz_sub(tmp, n, c);
	if (!mpz_congruent_p(foo1, tmp, n)) {
		ret(0); // composite
	}

	/*
	 * (5) Let n^2-1=2^r*s, where s is odd.  If x^s not congruent 1 mod (n,
	 * x^2-bx-c), and x^(2^j*s) not congruent -1 mod (n, x^2-bx-c) for all
	 * 0≤j≤r-2, declare n to be composite and stop.
	 * If n is not declared composite in Steps 1—5, declare n to be a probable
	 * prime.
	 */
	mpz_mul(tmp, n, n);
	r = split(s, tmp); // calculate r,s such that 2^r*s + 1 == n^2
	pow_mod(foo0, foo1, x0, x1, s, n, b, c);
	mpz_sub_ui(tmp, n, 1);
	if (mpz_sgn(foo0) == 0 && mpz_cmp_ui(foo1, 1) == 0)
		ret(1); // probably prime
	for (unsigned long i = 0; i < r-1; i++) {
		if (mpz_sgn(foo0) == 0 && mpz_congruent_p(foo1, tmp, n))
			ret(1); // probably prime
		square_mod(foo0, foo1, foo0, foo1, n, b, c);
	}
exit:
	mpz_clears(x0, x1, s, tmp, foo0, foo1, NULL);
	return result; // composite
}

int RQFT(mpz_t n, unsigned B)
{
	mpz_t b, c, n1;
	mpz_t bb4c, neg_c, tmp;
	int j1 = 0, j2 = 0, foo, result;

	assert(mpz_cmp_ui(n, 1) > 0);
	assert(mpz_odd_p(n));

	mpz_inits(b, c, n1, bb4c, neg_c, tmp, NULL);
	mpz_sub_ui(n1, n, 1);

	foo = no_nontrivial_small_prime_divisor(n, B);
	if (!foo)
		ret(0); // composite
	else if (foo == -1)
		ret(1); // n occurs in small_primes and is therefore certainly prime

	if (mpz_perfect_square_p(n))
		ret(0);

	for (unsigned i = 0; i < B; i++) {
		mpz_urandomm(b, r_state, n1);
		mpz_urandomm(c, r_state, n1);
		mpz_mul(bb4c, b, b);
		mpz_addmul_ui(bb4c, c, 4);
		j1 = mpz_jacobi(bb4c, n);
		mpz_neg(neg_c, c);
		j2 = mpz_jacobi(neg_c, n);
		if (j1 == -1 && j2 == 1) {
			mpz_gcd(tmp, bb4c, n);
			if (mpz_cmp_ui(tmp, 1) != 0 && mpz_cmp(tmp, n) != 0)
				ret(0);
			mpz_gcd(tmp, b, n);
			if (mpz_cmp_ui(tmp, 1) != 0 && mpz_cmp(tmp, n) != 0)
				ret(0);
			mpz_gcd(tmp, c, n);
			if (mpz_cmp_ui(tmp, 1) != 0 && mpz_cmp(tmp, n) != 0)
				ret(0);
			break;
		}
	}
	if (j1 != -1 || j2 != 1) {
		gmp_printf("Found no suitable pair (b,c) modulo n=%Zd.  This is highly unlikely unless the programme is wrong.  Assuming n is a prime...\n", n);
		ret(1);
	}
	
	result = QFT(n, b, c, B, 1);
exit:
	mpz_clears(b, c, n1, bb4c, neg_c, tmp, NULL);
	return result;
}

void frob_init(unsigned *large_primes, mpz_t *composites)
{
	FILE *file;

	init();

	file = fopen("primelist.txt", "r");
	while (EOF != fscanf(file, "%d\n", large_primes++));
	fclose(file);

	file = fopen("composites.txt", "r");
	do {
		mpz_init(*composites);
	} while (EOF != gmp_fscanf(file, "%Zd\n", *(composites++)));
	fclose(file);

	mpz_inits(tmp0, tmp1, tmp2, NULL);
}

int main()
{
	unsigned B = 50000, k = 1000;
	unsigned long n;
	static unsigned large_primes[3069262];
	static mpz_t composites[1013];
	mpz_t tmp;

	mpz_init(tmp);

	frob_init(large_primes, composites);

	for (int counter = 0; counter < 10; counter++) {
		// test primes
		for (int i = 0; i < k; i++) {
			n = randint(0, len(large_primes)-1);
			mpz_set_ui(tmp, large_primes[n]);
			//assert(RQFT(tmp, B));
		}

		// test composites
		for (int i = 0; i < k; i++) {
			int j = randint(0, len(large_primes)-1);
			if (large_primes[j] + 2 == large_primes[j+1])
				j++;
			do
				n = randint(large_primes[j]+1, large_primes[j+1]-1);
			while (n % 2 == 0);
			mpz_set_ui(tmp, n);
			assert(! RQFT(tmp, B));
		}

		// test some stuff
		mpz_set_ui(tmp, 7);
		assert(RQFT(tmp, B));

		mpz_set_ui(tmp, 611879L * 611957 * 611957);
		mpz_mul(tmp, tmp, tmp);
		assert(! RQFT(tmp, B));

		mpz_set_ui(tmp, 1);
		mpz_mul_2exp(tmp, tmp, 1<<15);
		mpz_add_ui(tmp, tmp, 1);
		//mpz_set_str(tmp, "1235790412356789098765432827498274392743929834792843282734279348239482349", 10);
		//mpz_pow_ui(tmp, tmp, 9);
		assert(! RQFT(tmp, B));

		// test squarefree composites
		for (int i = 0; i < len(composites); i++) {
			mpz_set(tmp, composites[i]);
			assert(! RQFT(tmp, B));
		}
	}

	mpz_clear(tmp);
	cleanup();
	mpz_clears(tmp0, tmp1, tmp2, NULL);
}
