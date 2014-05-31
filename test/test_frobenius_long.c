#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "test_frobenius_long.h"

#define TEST
#include "../frobenius.c"
#include "../frobenius_int.c"

static int num_iterations = 10000;

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

	b = c = 0;
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
		unsigned long n = 1 + (unsigned long)rand() % ((1<<15) - 1);
		mpz_t n_;
		mpz_init(n_);

		n = n * n;
		mpz_set_ui(n_, n);

		CU_ASSERT_FALSE(steps_1_2(n_));
		CU_ASSERT_FALSE(RQFT(n_, 1));
	}
}


void test_frobenius_trial_division(void)
{
	for (int i = 0; i < num_iterations; i++) {
	}
}


void test_frobenius_rqft_small_primes(void)
{
	for (int i = 0; i < num_iterations; i++) {
	}
}


void test_frobenius_rqft_small_composites(void)
{
	for (int i = 0; i < num_iterations; i++) {
	}
}


void test_frobenius_problematic_primes(void)
{
	for (int i = 0; i < num_iterations; i++) {
	}
}


void test_frobenius_primelist(void)
{
	for (int i = 0; i < num_iterations; i++) {
	}
}


void test_frobenius_larger_primes(void)
{
	for (int i = 0; i < num_iterations; i++) {
	}
}


