#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "test_frobenius_long.h"

#include "../frobenius.c"
#include "../frobenius_int.c"

static int num_iterations = 10000;

void test_frobenius_power(void)
{
	unsigned long r, n_;
	mpz_t MODULUS, s, t, poly(foo), poly(bar), poly(x), tmp1, tmp2;
	mpz_inits(MODULUS, s, t, poly(foo), poly(bar), poly(x), tmp1, tmp2, NULL);

	// Initialize poly(x)
	mpz_set_ui(x_x, 1);
	n_ = 2500000001;
	// n
	mpz_set_ui(n, n_);

	// Make sure the basics are working
	mult_x_mod(poly(foo), poly(x), MODULUS);
	CU_ASSERT_FATAL(mpz_cmp(foo_x, b) == 0);
	CU_ASSERT_FATAL(mpz_cmp(foo_1, c) == 0);

	for (unsigned long c_ = 1; c_ < 100; c_++) {
		if (jacobi(n_ - c_, n_) != 1)
			continue;
		mpz_set_ui(c, c_);

		for (unsigned long b_ = 1; b_ < 100; b_++) {
			if (jacobi(b_*b_+4*c_, n_) != -1)
				continue;
			mpz_set_ui(b, b_);

			// Compute r' and s' such that 2^r' s' = n-1 (as n is 1 mod 4)
			mpz_sub_ui(tmp1, n, 1);
			split(&r, s, tmp1);

			// Directly compute poly(foo) = x^s'
			powm(poly(foo), poly(x), s, MODULUS);

			// And compute it via x^t
			mpz_fdiv_q_2exp(t, s, 1);
			powm(poly(bar), poly(x), t, MODULUS);
			square_mod(poly(bar), poly(bar), MODULUS);
			mult_x_mod(poly(bar), poly(bar), MODULUS);

			// Make sure that the two methods agree...
			CU_ASSERT(mpz_cmp(foo_x, bar_x) == 0);
			CU_ASSERT(mpz_cmp(foo_1, bar_1) == 0);

			// Calculate x^(n+1) directly
			mpz_add_ui(tmp1, n, 1);
			powm(poly(foo), poly(x), tmp1, MODULUS);

			// Given x^s, comput x^((n-1)/2) = (x^s')^(2^r'-1)
			for (unsigned long i = 0; i < r-1; i++)
				square_mod(poly(bar), poly(bar), MODULUS);
			// Now x^((n+1)/2) = x^((n-1)/2) * x
			mult_x_mod(poly(bar), poly(bar), MODULUS);
			// and x^(n+1) = (x^((n+1)/2))^2
			square_mod(poly(bar), poly(bar), MODULUS);

			CU_ASSERT(mpz_cmp(foo_x, bar_x) == 0);
			CU_ASSERT(mpz_cmp(foo_1, bar_1) == 0);

			CU_ASSERT(mpz_sgn(foo_x) == 0);
			mpz_sub(tmp1, n, c);
			CU_ASSERT(mpz_cmp(foo_1, tmp1) == 0);
		}
	}

	mpz_clears(MODULUS, s, t, poly(foo), poly(bar), poly(x), tmp1, tmp2, NULL);
}

void test_frobenius_split(void)
{
	unsigned long max = 1000*1000*1000;

	unsigned long n, s, d;

	unsigned long s_gmp;
	mpz_t n_gmp, d_gmp;
	mpz_inits(n_gmp, d_gmp, NULL);

	for (int i = 0; i < num_iterations; i++) {
		n = 2 * ((unsigned long)rand() % max) + 1;
		mpz_set_ui(n_gmp, n);

		split_int(&s, &d, n);
		split(&s_gmp, d_gmp, n_gmp);

		CU_ASSERT_EQUAL(s_gmp, s);
		CU_ASSERT_TRUE(mpz_cmp_ui(d_gmp, d) == 0);
	}
}

void test_frobenius_mult_mod(void)
{
	unsigned long b, c, n;
	unsigned long d, e, f, g;
	unsigned long res0, res1;

	mpz_t b_, c_, n_;
	mpz_t d_, e_, f_, g_;
	mpz_t res0_, res1_;
	mpz_inits(b_, c_, n_, d_, e_, f_, g_, res0_, res1_, NULL);

	n = 131071;

	mpz_set_ui(n_, n);

	b = c = 0;
	for (int i = 0; i < num_iterations; i++) {
		while (jacobi(n-c, n) != 1)
			c = (unsigned long)rand() % n;
		while (jacobi(b*b+4*c, n) != -1)
			b = (unsigned long)rand() % n;
		d = (unsigned long)rand() % n;
		e = (unsigned long)rand() % n;
		f = (unsigned long)rand() % n;
		g = (unsigned long)rand() % n;

		mpz_set_ui(b_, b);
		mpz_set_ui(c_, c);
		mpz_set_ui(d_, d);
		mpz_set_ui(e_, e);
		mpz_set_ui(f_, f);
		mpz_set_ui(g_, g);

		mult_mod_int(&res0, &res1, d, e, f, g, n, b, c);
		mult_mod(res0_, res1_, d_, e_, f_, g_, n_, b_, c_);

		CU_ASSERT_FATAL(mpz_cmp_ui(res0_, res0) == 0);
		CU_ASSERT_FATAL(mpz_cmp_ui(res1_, res1) == 0);

		square_mod_int(&res0, &res1, d, e, n, b, c);
		square_mod(res0_, res1_, d_, e_, n_, b_, c_);

		CU_ASSERT_FATAL(mpz_cmp_ui(res0_, res0) == 0);
		CU_ASSERT_FATAL(mpz_cmp_ui(res1_, res1) == 0);
	}
}


void test_frobenius_powm_mod(void)
{
	unsigned long b, c, n;
	unsigned long d, e, k;
	unsigned long res0, res1;

	mpz_t b_, c_, n_;
	mpz_t d_, e_, k_;
	mpz_t res0_, res1_;
	mpz_inits(b_, c_, n_, d_, e_, k_, res0_, res1_, NULL);

	n = 131071;

	mpz_set_ui(n_, n);

	for (int i = 0; i < num_iterations; i++) {
		do
			c = (unsigned long)rand() % n;
		while (jacobi(n-c, n) != 1);
		do
			b = (unsigned long)rand() % n;
		while (jacobi(b*b+4*c, n) != -1);
		d = (unsigned long)rand() % n;
		e = (unsigned long)rand() % n;
		k = (unsigned long)rand() % n;

		mpz_set_ui(b_, b);
		mpz_set_ui(c_, c);
		mpz_set_ui(d_, d);
		mpz_set_ui(e_, e);
		mpz_set_ui(k_, k);

		powm_int(&res0, &res1, d, e, k, n, b, c);
		powm(res0_, res1_, d_, e_, k_, n_, b_, c_);

		CU_ASSERT(mpz_cmp_ui(res0_, res0) == 0);
		CU_ASSERT(mpz_cmp_ui(res1_, res1) == 0);
	}
}


void test_frobenius_squares(void)
{
	for (int i = 0; i < num_iterations; i++) {
		unsigned long n = 2 + (unsigned long)rand() % ((1<<15) - 2);
		mpz_t n_;
		n = n * n;
		mpz_init_set_ui(n_, n);

		CU_ASSERT_FALSE(steps_1_2(n_));
		CU_ASSERT_FALSE(RQFT(n_, 1));
	}
}


void test_frobenius_trial_division(void)
{
	for (int i = 0; i < num_iterations; i++) {
		unsigned long n = 2 + (unsigned long)rand() % ((1lu<<31) - 2);
		mpz_t n_;
		mpz_init_set_ui(n_, n);

		CU_ASSERT_EQUAL_FATAL(steps_1_2(n_), steps_1_2_int(n));
	}
}


void test_frobenius_rqft_small_primes(void)
{
	for (int i = 0; i < num_iterations; i++) {
		unsigned long n = 2 + (unsigned long)rand() % ((1lu<<30) - 2);
		mpz_t n_;
		mpz_init_set_ui(n_, n);

		CU_ASSERT_EQUAL_FATAL(RQFT(n_, 1), RQFT_int(n, 1));
	}
}



void test_frobenius_primelist(void)
{
	FILE *fp = fopen(TEST_DATA_PATH "primelist.txt", "r");
	unsigned p;
	unsigned long i = 0;
	static unsigned large_primes[3069262];
	mpz_t n;

	if (NULL == fp)
		die(TEST_DATA_PATH "primelist.txt: %s\n", strerror(errno));

	while (EOF != fscanf(fp, "%u\n", &p))
		large_primes[i++] = p;

	fclose(fp);

	mpz_init(n);
	for (i = 0; i < len(large_primes); i++) {
		Primality foo;
		mpz_set_ui(n, large_primes[i]);
		foo = RQFT(n, 1);
		if (foo == composite && large_primes[i] % 16 != 15)
			printf("%x\n", large_primes[i]);
		CU_ASSERT_NOT_EQUAL(foo, composite);
	}
}


void test_frobenius_larger_primes(void)
{
	mpz_t n;
	mpz_init_set_ui(n, 2500000001);
	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);

	mpz_set_ui(n, 2500000033);
	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);

	mpz_set_ui(n, 2500000039);
	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);

	mpz_set_ui(n, 2500000043);
	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);

	mpz_set_ui(n, 2500000057);
	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);
}


