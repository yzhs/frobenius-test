#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "test_frobenius_long.h"

#include "../frobenius.c"
#include "../frobenius_int.c"

/*
 * Compute b^exponent mod (n, x² - bx - c) where b is the polynomial b_x*x + b_1.
 */
static void powm(POLY_ARGS(res), CONST_POLY_ARGS(b), const mpz_t exponent, MODULUS_ARGS)
{
	mpz_t POLY(base);
	mpz_inits(POLY(base), NULL);

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

	mpz_clears(POLY(base), NULL);
}

static int num_iterations = 10000;

void frob_mult_x(void)
{
	uint64_t n_;
#define VARIABLES MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

	// Initialize POLY(x)
	mpz_set_ui(x_x, 1);

	n_ = 0x7ffffff;
	mpz_set_ui(n, n_);
	mpz_urandomm(n, r_state, n);
	mpz_add(n, n, n);
	mpz_add_ui(n, n, 0x70000001);
	n_ = mpz_get_ui(n);

	for (uint64_t c_ = 1; c_ < 1000; c_++) {
		if (jacobi(n_ - c_, n_) != 1)
			continue;
		mpz_set_ui(c, c_);

		for (uint64_t b_ = 1; b_ < 100; b_++) {
			b_ %= n_;
			if (jacobi(b_*b_+4*c_, n_) != -1)
				continue;
			mpz_set_ui(b, b_);

			for (int i = 0; i < 100; i++) {
				mpz_urandomm(f_1, r_state, n);
				mpz_set_ui(f_x, 0);

				mpz_mul(bb4c, b, b);
				mpz_addmul_ui(bb4c, c, 4);

				// Ensure that (0x+a)*x = a*x for all a in ℤ_n
				mult_x_mod(POLY(foo), POLY(f), MODULUS);
				CU_ASSERT(mpz_cmp(foo_x, f_1) == 0);
				CU_ASSERT(mpz_sgn(foo_1) == 0);

				mpz_urandomm(f_x, r_state, n);

				mult_x_mod(POLY(foo), POLY(f), MODULUS);
				mult_mod(POLY(bar), POLY(f), POLY(x), MODULUS);
				CU_ASSERT(mpz_cmp(foo_x, bar_x) == 0);
				CU_ASSERT(mpz_cmp(foo_1, bar_1) == 0);
			}
		}
	}

	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}

void frob_sigma_basics(void)
{
	uint64_t n_;
#define VARIABLES MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

	// Initialize POLY(x)
	mpz_set_ui(x_x, 1);

	n_ = 0x7ffffff;
	mpz_set_ui(n, n_);
	mpz_urandomm(n, r_state, n);
	mpz_add(n, n, n);
	mpz_add_ui(n, n, 0x70000001);
	n_ = mpz_get_ui(n);

	for (uint64_t c_ = 1; c_ < 1000; c_++) {
		if (jacobi(n_ - c_, n_) != 1)
			continue;
		mpz_set_ui(c, c_);

		for (uint64_t b_ = 1; b_ < 100; b_++) {
			b_ %= n_;
			if (jacobi(b_*b_+4*c_, n_) != -1)
				continue;
			mpz_set_ui(b, b_);

			mpz_mul(bb4c, b, b);
			mpz_addmul_ui(bb4c, c, 4);

			mpz_urandomm(f_x, r_state, n);
			mpz_urandomm(f_1, r_state, n);

			// Test whether sigma(x) = b - x.
			sigma(POLY(foo), POLY(x), MODULUS);
			mpz_sub_ui(baz, n, 1);
			CU_ASSERT(mpz_cmp(foo_x, baz) == 0);
			CU_ASSERT(mpz_cmp(foo_1, b) == 0);

			// Make sure σ is indeed an involution
			sigma(POLY(foo), POLY(f), MODULUS);
			sigma(POLY(bar), POLY(foo), MODULUS);
			CU_ASSERT(mpz_cmp(bar_x, f_x) == 0);
			CU_ASSERT(mpz_cmp(bar_1, f_1) == 0);
			sigma(POLY(bar), POLY(bar), MODULUS);
			CU_ASSERT(mpz_cmp(bar_x, foo_x) == 0);
			CU_ASSERT(mpz_cmp(bar_1, foo_1) == 0);
		}
	}

	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}

void frob_sigma_short_integer(void)
{
	uint64_t n_, b_, c_, p;
#define VARIABLES MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x)
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

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

		mpz_mul(bb4c, b, b);
		mpz_addmul_ui(bb4c, c, 4);

		for (p = 0; p < 20; p++) {
			mpz_set_ui(f_x, 0);
			mpz_set_ui(f_1, p);
			poly_powm(foo, f, n);
			CU_ASSERT(mpz_cmp_ui(foo_x, 0) == 0);
			CU_ASSERT(mpz_cmp_ui(foo_1, p%n_) == 0);
		}
	}

	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}

void frob_sigma_power(void)
{
	uint64_t n_;
#define VARIABLES MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

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

	for (uint64_t c_ = 1; c_ < 1000; c_++) {
		if (jacobi(n_ - c_, n_) != 1)
			continue;
		mpz_set_ui(c, c_);

		for (uint64_t b_ = 1; b_ < 100; b_++) {
			b_ %= n_;
			if (jacobi(b_*b_+4*c_, n_) != -1)
				continue;
			mpz_set_ui(b, b_);

			mpz_mul(bb4c, b, b);
			mpz_addmul_ui(bb4c, c, 4);

			mpz_urandomm(f_x, r_state, n);
			mpz_urandomm(f_1, r_state, n);

			// Check whether x^n is the same as σ(x)
			poly_sigma(foo, x);
			poly_powm(bar, x, n);
			CU_ASSERT(mpz_cmp(foo_x, bar_x) == 0);
			CU_ASSERT(mpz_cmp(foo_1, bar_1) == 0);
		}
	}

	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}

void frob_power_x_lucas(void)
{
	uint64_t n_;
#define VARIABLES MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

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

	for (uint64_t c_ = 1; c_ < 50; c_++) {
		if (jacobi(n_ - c_, n_) != 1)
			continue;
		mpz_set_ui(c, c_);

		for (uint64_t b_ = 1; b_ < 20; b_++) {
			if (jacobi(b_*b_+4*c_, n_) != -1)
				continue;
			mpz_set_ui(b, b_);

			mpz_mul(bb4c, b, b);
			mpz_addmul_ui(bb4c, c, 4);

			mpz_urandomm(f_x, r_state, n);
			mpz_urandomm(f_1, r_state, n);

			for (uint64_t k = 1; k < 1000; k++) {
				// Compare x*x^k and x^(k+1), both computed
				// using power_of_x.
				mpz_set_ui(baz, k);
				poly_pow_x(foo, baz);
				poly_mul_x(foo, foo);

				mpz_set_ui(baz, k+1);
				poly_pow_x(bar, baz);

				CU_ASSERT(poly_eq(foo, bar));
			}

			for (uint64_t k = 1; k < 1000; k++) {
				mpz_set_ui(baz, k);
				poly_pow_x(bar, baz);
				poly_powm(foo, x, baz);

				CU_ASSERT(poly_eq(foo, bar));
			}


		}
	}

	mpz_clears(VARIABLES, NULL);
}

void frob_power_basics(void)
{
	uint64_t n_;
#define VARIABLES MODULUS, POLY(f), POLY(foo), POLY(bar), POLY(x), baz
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

	// Initialize POLY(x)
	n_ = 0x7fffffff;
	mpz_set_ui(n, n_);
	mpz_urandomm(n, r_state, n);
	mpz_add(n, n, n);
	mpz_add_ui(n, n, 0x70000001);
	n_ = mpz_get_ui(n);
	mpz_set_ui(x_x, 1);

	// Make sure the basics are working
	mult_x_mod(POLY(foo), POLY(x), MODULUS);
	CU_ASSERT_FATAL(mpz_cmp(foo_x, b) == 0);
	CU_ASSERT_FATAL(mpz_cmp(foo_1, c) == 0);

	mpz_nextprime(n, n);
	CU_ASSERT_FATAL(mpz_probab_prime_p(n, 50));
	n_ = mpz_get_ui(n);

	// Make sure the basics are working
	mult_x_mod(POLY(foo), POLY(x), MODULUS);
	CU_ASSERT_FATAL(mpz_cmp(foo_x, b) == 0);
	CU_ASSERT_FATAL(mpz_cmp(foo_1, c) == 0);

	for (uint64_t c_ = 1; c_ < 50; c_++) {
		if (jacobi(n_ - c_, n_) != 1)
			continue;
		mpz_set_ui(c, c_);

		for (uint64_t b_ = 1; b_ < 20; b_++) {
			if (jacobi(b_*b_+4*c_, n_) != -1)
				continue;
			mpz_set_ui(b, b_);

			mpz_mul(bb4c, b, b);
			mpz_addmul_ui(bb4c, c, 4);

			mpz_urandomm(f_x, r_state, n);
			mpz_urandomm(f_1, r_state, n);

			// Make sure that x^n = b - x
			mpz_sub_ui(baz, n, 1);
			poly_powm(foo, x, n);
			CU_ASSERT(mpz_cmp(foo_x, baz) == 0);
			CU_ASSERT(mpz_cmp(foo_1, b) == 0);
		}
	}

	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}

void frob_inverse(void)
{
	uint64_t n_;
#define VARIABLES MODULUS, POLY(foo), POLY(bar), POLY(baz), POLY(x), tmp1, tmp2
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

	// Initialize POLY(x)
	mpz_set_ui(x_x, 1);
	for (n_ = 2500000001; n_ < 2500000001 + 1000; n_+=4) {
		// n
		mpz_set_ui(n, n_);

		for (uint64_t c_ = 1; c_ < 100; c_++) {
			if (jacobi(n_ - c_, n_) != 1)
				continue;
			mpz_set_ui(c, c_);

			for (uint64_t b_ = 1; b_ < 100; b_++) {
				if (jacobi(b_*b_+4*c_, n_) != -1)
					continue;
				mpz_set_ui(b, b_);

				mpz_mul(bb4c, b, b);
				mpz_addmul_ui(bb4c, c, 4);

				mpz_urandomm(foo_x, r_state, n);
				mpz_urandomm(foo_1, r_state, n);

				if (!poly_invert(bar, foo))
					continue;
				if (!poly_invert(baz, bar))
					continue;

				CU_ASSERT(poly_eq(foo, baz));

				poly_mul(baz, foo, bar);

				CU_ASSERT(mpz_sgn(baz_x) == 0);
				CU_ASSERT(mpz_cmp_ui(baz_1, 1) == 0);
			}
		}
	}

	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}

void frob_fast_algorithm1(void)
{
	uint64_t r, n_;
#define VARIABLES MODULUS, s, t, POLY(foo), POLY(bar), POLY(baz), POLY(x), tmp1, tmp2
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

	// Initialize POLY(x)
	mpz_set_ui(x_x, 1);
	for (n_ = 2500000001; n_ < 2500000040; n_+=38) {
		// n
		mpz_set_ui(n, n_);

		for (uint64_t c_ = 1; c_ < 100; c_++) {
			if (jacobi(n_ - c_, n_) != 1)
				continue;
			mpz_set_ui(c, c_);

			for (uint64_t b_ = 1; b_ < 100; b_++) {
				if (jacobi(b_*b_+4*c_, n_) != -1)
					continue;
				mpz_set_ui(b, b_);

				mpz_mul(bb4c, b, b);
				mpz_addmul_ui(bb4c, c, 4);

				// Compute r' and s' such that 2^r' s' = n-1 (as n is 1 mod 4)
				mpz_sub_ui(tmp1, n, 1);
				split(&r, s, tmp1);

				// Directly compute POLY(foo) = x^s'
				poly_powm(foo, x, s);

				// And compute it via x^t
				mpz_fdiv_q_2exp(t, s, 1);
				poly_powm(bar, x, t);
				poly_pow_x(baz, t);
				CU_ASSERT(poly_eq(bar, baz));

				poly_sqr(bar, bar);
				poly_mul_x(bar, bar);

				// Make sure that the two methods agree...
				CU_ASSERT(poly_eq(foo, bar));

				// Calculate x^(n+1) directly
				mpz_add_ui(tmp1, n, 1);
				poly_powm(foo, x, tmp1);
				poly_pow_x(baz, tmp1);
				CU_ASSERT(poly_eq(foo, baz));

				// Given x^s, comput x^((n-1)/2) = (x^s')^(2^r'-1)
				for (uint64_t i = 0; i < r-1; i++)
					poly_sqr(bar, bar);
				// Now x^((n+1)/2) = x^((n-1)/2) * x
				poly_mul_x(bar, bar);
				// and x^(n+1) = (x^((n+1)/2))^2
				poly_sqr(bar, bar);

				CU_ASSERT(poly_eq(foo, bar));

				CU_ASSERT(mpz_sgn(foo_x) == 0);
				mpz_sub(tmp1, n, c);
				CU_ASSERT(mpz_cmp(foo_1, tmp1) == 0);
			}
		}
	}

	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}

void frob_fast_algorithm2(void)
{
	uint64_t r_prime, r, n_;
#define VARIABLES n, b, c, s_prime, s, t, POLY(foo), POLY(bar), POLY(baz), POLY(x), tmp1, tmp2, tmp3
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

	// Initialize POLY(x)
	mpz_set_ui(x_x, 1);
	for (n_ = 7; n_ < 20000; n_+=4) {
		// n
		mpz_set_ui(n, n_);
		if (!mpz_probab_prime_p(n, 100))
			continue;

		// Make sure that the exponents behave as expected.
		mpz_add_ui(tmp1, n, 1); // n+1
		split(&r_prime, s_prime, tmp1);
		mpz_fdiv_q_2exp(t, s_prime, 1); // (s'-1)/2

		mpz_mul(tmp2, n, n);
		mpz_sub_ui(tmp2, tmp2, 1); // n^2-1
		split(&r, s, tmp2);

		mpz_fdiv_q_2exp(tmp1, tmp1, 1); // (n+1)/2
		mpz_addmul(tmp1, n, t); // + nt
		mpz_sub(tmp1, tmp1, t); // - t
		mpz_sub_ui(tmp1, tmp1, 1); // - 1

		CU_ASSERT_EQUAL_FATAL(r, r_prime+1);
		CU_ASSERT_FATAL(mpz_cmp(tmp1, s) == 0);
	}
	mpz_clears(VARIABLES, NULL);
}

void frob_fast_algorithm3(void)
{
	int counter = 0;
	uint64_t r, n_;
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);
	mpz_set_ui(x_x, 1);

	for (n_ = 7; n_ < 200; n_+=4) {
		// n
		mpz_set_ui(n, n_);
		if (!mpz_probab_prime_p(n, 100))
			continue;

		mpz_mul(tmp2, n, n);
		mpz_sub_ui(tmp2, tmp2, 1); // n^2-1
		split(&r, s, tmp2);
		// TODO Fix this test

		for (uint64_t c_ = 1; c_ < 100 && c_ < n_; c_++) {
			if (jacobi(n_ - c_, n_) != 1)
				continue;
			mpz_set_ui(c, c_);

			for (uint64_t b_ = 1; b_ < 100 && b_ < n_; b_++) {
				if (jacobi(b_*b_+4*c_, n_) != -1)
					continue;
				mpz_set_ui(b, b_);

				mpz_mul(bb4c, b, b);
				mpz_addmul_ui(bb4c, c, 4);

				poly_pow_x(foo, s); // foo = x^s

				poly_pow_x(bar, t); // bar = x^t
				poly_powm(bar, bar, n); // bar = x^(nt)
				mpz_cdiv_q_2exp(tmp1, n, 1);
				poly_pow_x(baz, tmp1); // baz = x^((n+1)/2)
				poly_mul(bar, bar, baz); // bar = x^(nt + (n+1)/2)

				poly_pow_x(baz, t); // baz = x^t
				poly_mul_x(baz, baz); // baz = x^(t+1)
				poly_invert(baz, baz); // baz = x^(-(t+1))

				poly_mul(bar, bar, baz); // bar = x^(nt+(n+1)/2-(t+1))

				if (!poly_eq(foo, bar)) {
					if (counter == 0)
						puts("");
					gmp_printf("%2Zd x + %2Zd ≠ %2Zd x + %2Zd, n = %2Zd, b = %2Zd, c = %2Zd\n",
					           foo_x, foo_1, bar_x, bar_1, MODULUS);
					counter++;
				}
				CU_ASSERT(poly_eq(foo, bar));

				if (counter == 20)
					goto exit;
			}
		}
	}

exit:
	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}

void frob_split(void)
{
	uint64_t max = 1000*1000*1000;

	uint64_t n, s, d;

	uint64_t s_gmp;
#define VARIABLES n_gmp, d_gmp
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

	for (int i = 0; i < num_iterations; i++) {
		n = 2 * ((uint64_t)rand() % max) + 1;
		mpz_set_ui(n_gmp, n);

		split_int(&s, &d, n);
		split(&s_gmp, d_gmp, n_gmp);

		CU_ASSERT_EQUAL(s_gmp, s);
		CU_ASSERT_TRUE(mpz_cmp_ui(d_gmp, d) == 0);
	}
	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}

void frob_mult_mod(void)
{
	uint64_t b, c, n;
	uint64_t d, e, f, g;
	uint64_t res0, res1;

#define VARIABLES b_, c_, n_, d_, e_, f_, g_, res0_, res1_
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

	n = 131071;

	mpz_set_ui(n_, n);

	b = c = 0;
	for (int i = 0; i < num_iterations; i++) {
		while (jacobi(n-c, n) != 1)
			c = (uint64_t)rand() % n;
		while (jacobi(b*b+4*c, n) != -1)
			b = (uint64_t)rand() % n;
		d = (uint64_t)rand() % n;
		e = (uint64_t)rand() % n;
		f = (uint64_t)rand() % n;
		g = (uint64_t)rand() % n;

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
	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}


void frob_powm_mod(void)
{
	uint64_t b, c, n;
	uint64_t d, e, k;
	uint64_t res0, res1;

#define VARIABLES b_, c_, n_, d_, e_, k_, res0_, res1_
	mpz_t VARIABLES;
	mpz_inits(VARIABLES, NULL);

	n = 131071;

	mpz_set_ui(n_, n);

	for (int i = 0; i < num_iterations; i++) {
		do
			c = (uint64_t)rand() % n;
		while (jacobi(n-c, n) != 1);
		do
			b = (uint64_t)rand() % n;
		while (jacobi(b*b+4*c, n) != -1);
		d = (uint64_t)rand() % n;
		e = (uint64_t)rand() % n;
		k = (uint64_t)rand() % n;

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

	mpz_clears(VARIABLES, NULL);
#undef VARIABLES
}


void frob_squares(void)
{
	for (int i = 0; i < num_iterations; i++) {
		uint64_t n = 2 + (uint64_t)rand() % ((1<<15) - 2);
		mpz_t n_;
		n = n * n;
		mpz_init_set_ui(n_, n);

		CU_ASSERT_FALSE(steps_1_2(n_));
		CU_ASSERT_FALSE(RQFT(n_, 1));
		mpz_clear(n_);
	}
}


void frob_trial_division(void)
{
	for (int i = 0; i < num_iterations; i++) {
		uint64_t n = 2 + (uint64_t)rand() % ((1lu<<31) - 2);
		mpz_t n_;
		mpz_init_set_ui(n_, n);

		CU_ASSERT_EQUAL_FATAL(steps_1_2(n_), steps_1_2_int(n));
		mpz_clear(n_);
	}
}


void frob_rqft_small_primes(void)
{
	for (int i = 0; i < num_iterations; i++) {
		uint64_t n = 2 + (uint64_t)rand() % ((1lu<<30) - 2);
		mpz_t n_;
		mpz_init_set_ui(n_, n);

		CU_ASSERT_EQUAL_FATAL(RQFT(n_, 1), RQFT_int(n, 1));
		mpz_clear(n_);
	}
}



void frob_primelist(void)
{
	FILE *fp = fopen(TEST_DATA_PATH "primelist.txt", "r");
	unsigned p;
	uint64_t i = 0;
	static unsigned large_primes[23006167];
	mpz_t n;

	if (NULL == fp)
		die(TEST_DATA_PATH "primelist.txt: %s\n", strerror(errno));

	while (EOF != fscanf(fp, "%u\n", &p) && i < len(large_primes))
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
	mpz_clear(n);
}


void frob_larger_primes(void)
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

	mpz_clear(n);
}


void frob_composites(void)
{
	static mpz_t composites[1013];
	FILE *fp = fopen(TEST_DATA_PATH "composites.txt", "r");
	uint64_t i;

	for (i = 0; i < len(composites); i++)
		gmp_fscanf(fp, "%Zd\n", composites[i]);

	fclose(fp);

	for (i = 0; i < len(composites); i++)
		CU_ASSERT_EQUAL(RQFT(composites[i], 10), composite);
}

