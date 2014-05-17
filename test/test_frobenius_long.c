#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "test_frobenius_long.h"

#define TEST
#include "../frobenius.c"
#include "../frobenius_int.c"

static int num_iterations = 100;

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


void test_frobenius_get_random(void)
{
	for (int i = 0; i < num_iterations; i++) {
	}
}


void test_frobenius_mult_mod(void)
{
	for (int i = 0; i < num_iterations; i++) {
	}
}


void test_frobenius_powm_mod(void)
{
	for (int i = 0; i < num_iterations; i++) {
	}
}


void test_frobenius_squares(void)
{
	for (int i = 0; i < num_iterations; i++) {
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


