#include <CUnit/CUnit.h>

#include "test_frobenius_int.h"

#define TEST
#include "../frobenius_int.c"

/*
 * Helper functions
 */
void test_frobenius_int_int_sqrt(void)
{
	CU_ASSERT_EQUAL(int_sqrt(0), 0);
	CU_ASSERT_EQUAL(int_sqrt(1), 1);
	CU_ASSERT_EQUAL(int_sqrt(15), 3);
	CU_ASSERT_EQUAL(int_sqrt(16), 4);
	CU_ASSERT_EQUAL(int_sqrt(27), 5);
	CU_ASSERT_EQUAL(int_sqrt(1l << 32), 1 << 16);
	CU_ASSERT_EQUAL(int_sqrt(1lu << 62), 1lu << 31);

	for (int i = 0; i < 1000000; i++) {
		unsigned long n = get_random_int(2, (1lu << 32) - 2);
		unsigned long sqrt = int_sqrt(n);
		CU_ASSERT_TRUE(sqrt * sqrt <= n);
		CU_ASSERT_TRUE((sqrt + 1) * (sqrt + 1) > n);
	}
}

void test_frobenius_int_gcd(void)
{
	CU_ASSERT_EQUAL(gcd(0, 1), 1);
	CU_ASSERT_EQUAL(gcd(1, 4), 1);
	CU_ASSERT_EQUAL(gcd(2, 4), 2);
	CU_ASSERT_EQUAL(gcd(3, 2), 1);
	CU_ASSERT_EQUAL(gcd(4, 6), 2);
	CU_ASSERT_EQUAL(gcd(12, 18), 6);

	for (unsigned long i = 1; i < 100; i++)
		for (unsigned long j = 1; j < 100; j++) {
			CU_ASSERT_EQUAL(i % gcd(i, j), 0);
			CU_ASSERT_EQUAL(j % gcd(i, j), 0);
		}
}

void test_frobenius_int_is_square(void)
{
	CU_ASSERT_TRUE(is_square(0));
	CU_ASSERT_TRUE(is_square(1));
	CU_ASSERT_TRUE(is_square(16));
	CU_ASSERT_FALSE(is_square(15));
	CU_ASSERT_FALSE(is_square(27));
}

void test_frobenius_int_jacobi(void)
{
	CU_ASSERT_EQUAL(jacobi(1, 2), 1);
	CU_ASSERT_EQUAL(jacobi(1, 3), 1);
	CU_ASSERT_EQUAL(jacobi(1, 5), 1);
	CU_ASSERT_EQUAL(jacobi(3 - 1, 3), -1);
	CU_ASSERT_EQUAL(jacobi(5 - 1, 5), 1);
	CU_ASSERT_EQUAL(jacobi(3, 5), -1);
	CU_ASSERT_EQUAL(jacobi(21, 7), 0);
}

void test_frobenius_int_split(void)
{
	unsigned long s, d, n;

	for (n = 5; n < 1000; n += 2) {
		split_int(&s, &d, n * n);
		CU_ASSERT_EQUAL_FATAL((1lu << s) * d, n * n - 1);
	}
}

void test_frobenius_int_get_random_int(void)
{
	for (int i = 0; i < 10000; i++) {
		unsigned long rand = get_random_int(2, 1234 - 2);
		CU_ASSERT_TRUE(2 <= rand);
		CU_ASSERT_TRUE(rand <= 1234 - 2);
	}
	for (int i = 0; i < 10000; i++) {
		unsigned long rand = get_random_int(2, 7 - 2);
		CU_ASSERT_TRUE(2 <= rand);
		CU_ASSERT_TRUE(rand <= 7 - 2);
	}
}

/*
 * The main arithmetic functions
 */
void test_frobenius_mult_mod_int(void)
{
	unsigned long b, c, n;
	//unsigned long d, e, f, g;
	unsigned long res0, res1;

	b = 55516;
	c = 108625;
	n = 131071;

	//mult_mod_int(&res0, &res1, d, e, f, g, n, b, c);
	powm_int(&res0, &res1, 1, 0, (n + 1) / 2, n, b, c);

	CU_ASSERT_EQUAL(res0, 1);
	CU_ASSERT_EQUAL(res1, 1);
}


void test_frobenius_squares_int(void)
{
	for (unsigned long i = 1; i < len(prime_list); i++) {
		unsigned long n = prime_list[i];
		n = n * n;
		CU_ASSERT_FALSE(steps_1_2_int(n));
		CU_ASSERT_FALSE(RQFT_int(n, 1));
	}
	CU_ASSERT_FALSE(steps_1_2_int(50021 * 50021lu));
	CU_ASSERT_FALSE(RQFT_int(50021 * 50021lu, 1));
}

/*
 * The various parts of the frobenius test
 */
void test_frobenius_trial_division_int(void)
{
	for (int i = 1; i < 1000; i++)
		for (unsigned long j = prime_list[i] * prime_list[i]; j < 50000; j += prime_list[i])
			CU_ASSERT_EQUAL(steps_1_2_int(j), composite);
}

/*
 * Full RQFT tests
 */
void test_frobenius_int_rqft_small_primes(void)
{
	for (unsigned long i = 2; i < len(prime_list); i++)
		CU_ASSERT_TRUE(RQFT_int(prime_list[i], 1));
}

void test_frobenius_int_rqft_small_composites(void)
{
	unsigned long i = 2, n = prime_list[i]; /* n = 5 */

	for (; n < B; n += 2) {
		Primality foo;
		if (prime_list[i] == n) {
			i++;
			continue;
		}
		foo = RQFT_int(n, 100);
		if (foo)
			printf("\n%lu", n);
		CU_ASSERT_EQUAL(foo, composite);
	}
}

void test_frobenius_problematic_primes_int(void)
{
	static unsigned primes[] = { 94207, 106367, 131071, 195071, 342191, 524287, 917503, 6561791 };

	CU_ASSERT_EQUAL(jacobi(55516 * 55516lu + 4 * 108625, 131071), -1);
	CU_ASSERT_EQUAL(jacobi(131071 - 108625, 131071), 1);
	CU_ASSERT_EQUAL(steps_1_2_int(131071), prime);
	CU_ASSERT_NOT_EQUAL(QFT_int(131071, 55516, 108625), composite);

	for (unsigned long i = 0; i < len(primes); i++) {
		Primality foo = RQFT_int(primes[i], 100);
		CU_ASSERT_EQUAL(steps_1_2_int(primes[i]), prime);
		CU_ASSERT_NOT_EQUAL(foo, composite);
	}
}

void test_frobenius_primelist_int(void)
{
	FILE *fp = fopen("data/primelist.txt", "r");
	unsigned p;
	unsigned long i = 0;
	static unsigned large_primes[3069262];

	if (NULL == fp)
		die("data/primelist.txt: %s\n", strerror(errno));

	while (EOF != fscanf(fp, "%u\n", &p))
		large_primes[i++] = p;

	fclose(fp);

	for (i = 0; i < len(large_primes); i++) {
		Primality foo = RQFT_int(large_primes[i], 1);
		if (foo == composite && large_primes[i] % 16 != 15)
			printf("%x\n", large_primes[i]);
		CU_ASSERT_NOT_EQUAL(foo, composite);
	}
}

// Test the first few primes larger than B^2.
// These are the smallest primes that
void test_frobenius_larger_primes_int(void)
{
	CU_ASSERT_EQUAL(RQFT_int(2500000001, 1), probably_prime);
	CU_ASSERT_EQUAL(RQFT_int(2500000033, 1), probably_prime);
	CU_ASSERT_EQUAL(RQFT_int(2500000039, 1), probably_prime);
	CU_ASSERT_EQUAL(RQFT_int(2500000043, 1), probably_prime);
	CU_ASSERT_EQUAL(RQFT_int(2500000057, 1), probably_prime);
}
