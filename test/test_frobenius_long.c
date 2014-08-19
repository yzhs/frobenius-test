#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "test_frobenius_long.h"

#include "../frobenius.c"
#include "../frobenius_int.c"

static int num_iterations = 10000;

void test_frobenius_sigma(void)
{
	unsigned long n_, b_, c_, p;
	mpz_t MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz;
	mpz_inits(MODULUS, POLY(foo), POLY(f), POLY(bar), POLY(x), baz, NULL);

	// Initialize POLY(x)
	mpz_set_ui(x_x, 1);

	n_ = 0x7ffffff;
	mpz_set_ui(n, n_);
	mpz_urandomm(n, r_state, n);
	mpz_add(n, n, n);
	mpz_add_ui(n, n, 0x70000001);
	n_ = mpz_get_ui(n);

	mpz_nextprime(n, n);
	CU_ASSERT_FATAL(mpz_probab_prime_p(n, 50));
	n_ = mpz_get_ui(n);

	for (c_ = 1; c_ < 50; c_++) {
		if (jacobi(n_ - c_, n_) != 1)
			continue;
		mpz_set_ui(c, c_);

		for (b_ = 1; b_ < 20; b_++) {
			if (jacobi(b_*b_+4*c_, n_) != -1)
				continue;
			mpz_set_ui(b, b_);
		}

		for (p = 0; p < 20; p++) {
			mpz_set_ui(f_x, 0);
			mpz_set_ui(f_1, p);
			powm(POLY(foo), POLY(f), n, MODULUS);
			CU_ASSERT(mpz_cmp_ui(foo_x, 0) == 0);
			CU_ASSERT(mpz_cmp_ui(foo_1, p) == 0);
		}
	}

	CU_ASSERT_FATAL(mpz_probab_prime_p(n, 50));

	for (unsigned long c_ = 1; c_ < 1000; c_++) {
		if (jacobi(n_ - c_, n_) != 1)
			continue;
		mpz_set_ui(c, c_);

		for (unsigned long b_ = 1; b_ < 100; b_++) {
			b_ %= n_;
			if (jacobi(b_*b_+4*c_, n_) != -1)
				continue;
			mpz_set_ui(b, b_);

			mpz_urandomm(f_x, r_state, n);
			mpz_urandomm(f_1, r_state, n);

			// Make sure σ is indeed an involution
			sigma(POLY(foo), POLY(f), MODULUS);
			sigma(POLY(foo), POLY(foo), MODULUS);
			CU_ASSERT_FATAL(mpz_cmp(foo_x, f_x) == 0);
			CU_ASSERT_FATAL(mpz_cmp(foo_1, f_1) == 0);

			sigma(POLY(foo), POLY(x), MODULUS);
			mpz_sub_ui(baz, n, 1);
			CU_ASSERT_FATAL(mpz_cmp(foo_x, baz) == 0);
			CU_ASSERT_FATAL(mpz_cmp(foo_1, b) == 0);

			powm(POLY(foo), POLY(x), n, MODULUS);
			CU_ASSERT_FATAL(mpz_cmp(foo_x, baz) == 0);
			CU_ASSERT_FATAL(mpz_cmp(foo_1, b) == 0);

			power_of_x(POLY(foo), n, MODULUS);
			CU_ASSERT_FATAL(mpz_cmp(foo_x, baz) == 0);
			CU_ASSERT_FATAL(mpz_cmp(foo_1, b) == 0);

			// Check whether x^n is the same as σ(x)
			sigma(POLY(foo), POLY(x), MODULUS);
			powm(POLY(bar), POLY(x), n, MODULUS);
			if (mpz_cmp(foo_x, bar_x) != 0 || mpz_cmp(foo_1, bar_1) != 0) {
				gmp_printf("\nError: σ(f) = %Zd x + %Zd, but\n"
						"        f^n = %Zd x + %Zd\n",\
						foo_x, foo_1, bar_x, bar_1);
				gmp_printf("n = %Zd = %d mod 4, b = %Zd, c = %Zd\n", n, mpz_fdiv_ui(n, 4), b, c);
				CU_ASSERT_FATAL(false);
			} else CU_ASSERT(true);
		}
	}

	mpz_clears(MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz, NULL);
}

void test_frobenius_power_x_lucas(void)
{
	unsigned long n_;
	mpz_t MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz;
	mpz_inits(MODULUS, POLY(foo), POLY(f), POLY(bar), POLY(x), baz, NULL);

	// Initialize POLY(x)
	n_ = 0x7ffffffff;
	mpz_set_ui(n, n_);
	mpz_urandomm(n, r_state, n);
	mpz_add(n, n, n);
	mpz_add_ui(n, n, 0x70000001);
	n_ = mpz_get_ui(n);
	mpz_set_ui(x_x, 1);

	mpz_nextprime(n, n);
	CU_ASSERT_FATAL(mpz_probab_prime_p(n, 50));

	for (unsigned long c_ = 1; c_ < 50; c_++) {
		if (jacobi(n_ - c_, n_) != 1)
			continue;
		mpz_set_ui(c, c_);

		for (unsigned long b_ = 1; b_ < 20; b_++) {
			if (jacobi(b_*b_+4*c_, n_) != -1)
				continue;
			mpz_set_ui(b, b_);

			mpz_urandomm(f_x, r_state, n);
			mpz_urandomm(f_1, r_state, n);

			for (unsigned long k = 1; k < 1000; k++) {
				mpz_set_ui(baz, k);
				power_of_x(POLY(bar), baz, MODULUS);
				powm(POLY(foo), POLY(x), baz, MODULUS);
				CU_ASSERT(mpz_cmp(foo_x, bar_x) == 0);
				CU_ASSERT(mpz_cmp(foo_1, bar_1) == 0);
			}


		}
	}

	mpz_clears(MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz, NULL);
}

void test_frobenius_power(void)
{
	unsigned long r, n_;
	mpz_t MODULUS, s, t, POLY(foo), POLY(bar), POLY(x), tmp1, tmp2;
	mpz_inits(MODULUS, s, t, POLY(foo), POLY(bar), POLY(x), tmp1, tmp2, NULL);

	// Initialize POLY(x)
	mpz_set_ui(x_x, 1);
	n_ = 2500000001;
	// n
	mpz_set_ui(n, n_);

	// Make sure the basics are working
	mult_x_mod(POLY(foo), POLY(x), MODULUS);
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

			// Directly compute POLY(foo) = x^s'
			powm(POLY(foo), POLY(x), s, MODULUS);

			// And compute it via x^t
			mpz_fdiv_q_2exp(t, s, 1);
			powm(POLY(bar), POLY(x), t, MODULUS);
			square_mod(POLY(bar), POLY(bar), MODULUS);
			mult_x_mod(POLY(bar), POLY(bar), MODULUS);

			// Make sure that the two methods agree...
			CU_ASSERT(mpz_cmp(foo_x, bar_x) == 0);
			CU_ASSERT(mpz_cmp(foo_1, bar_1) == 0);

			// Calculate x^(n+1) directly
			mpz_add_ui(tmp1, n, 1);
			powm(POLY(foo), POLY(x), tmp1, MODULUS);

			// Given x^s, comput x^((n-1)/2) = (x^s')^(2^r'-1)
			for (unsigned long i = 0; i < r-1; i++)
				square_mod(POLY(bar), POLY(bar), MODULUS);
			// Now x^((n+1)/2) = x^((n-1)/2) * x
			mult_x_mod(POLY(bar), POLY(bar), MODULUS);
			// and x^(n+1) = (x^((n+1)/2))^2
			square_mod(POLY(bar), POLY(bar), MODULUS);

			CU_ASSERT(mpz_cmp(foo_x, bar_x) == 0);
			CU_ASSERT(mpz_cmp(foo_1, bar_1) == 0);

			CU_ASSERT(mpz_sgn(foo_x) == 0);
			mpz_sub(tmp1, n, c);
			CU_ASSERT(mpz_cmp(foo_1, tmp1) == 0);
		}
	}

	mpz_clears(MODULUS, s, t, POLY(foo), POLY(bar), POLY(x), tmp1, tmp2, NULL);
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

		mult_mod_int(&res0, &res1, d, e, f, g, MODULUS);
		mult_mod(res0_, res1_, d_, e_, f_, g_, n_, b_, c_);

		CU_ASSERT_FATAL(mpz_cmp_ui(res0_, res0) == 0);
		CU_ASSERT_FATAL(mpz_cmp_ui(res1_, res1) == 0);

		square_mod_int(&res0, &res1, d, e, MODULUS);
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

		powm_int(&res0, &res1, d, e, k, MODULUS);
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
	static unsigned large_primes[23006166];
	mpz_t n;

	if (NULL == fp)
		die(TEST_DATA_PATH "primelist.txt: %s\n", strerror(errno));

	while (EOF != fscanf(fp, "%u\n", &p))
		large_primes[i++] = p;

	fclose(fp);

	mpz_init(n);
	for (i = 0; i < len(large_primes); i++) {
		Primality foo;
		if (large_primes[i] < B*B)
			continue;
		mpz_set_ui(n, large_primes[i]);
		foo = RQFT(n, 1);
		if (foo == composite)
			printf("%u\n", large_primes[i]);
		CU_ASSERT_NOT_EQUAL_FATAL(foo, composite);
	}
}


void test_frobenius_larger_primes(void)
{
	mpz_t n;
//	mpz_init_set_ui(n, 2500000001);
//	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);

//	mpz_set_ui(n, 2500000033);
//	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);

	mpz_set_ui(n, 2500000039);
	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);

	mpz_set_ui(n, 2500000043);
	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);

//	mpz_set_ui(n, 2500000057);
//	CU_ASSERT_EQUAL(RQFT(n, 1), probably_prime);
}


void test_frobenius_composites(void)
{
	static mpz_t composites[1013];
	FILE *fp = fopen(TEST_DATA_PATH "composites.txt", "r");
	unsigned long i;

	for (i = 0; i < len(composites); i++)
		gmp_fscanf(fp, "%Zd\n", composites[i]);

	fclose(fp);

	for (i = 0; i < len(composites); i++)
		CU_ASSERT_EQUAL(RQFT(composites[i], 10), composite);
}

