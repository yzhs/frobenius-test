#include <CUnit/CUnit.h>
#include <gmp.h>
#include "test_miller_rabin_int.h"

#define TEST
#include "../miller_rabin_int.c"

static unsigned large_primes[3069262];

void test_miller_rabin_powm_int(void)
{
	mpz_t base, mod, result;

	for (unsigned i = 0; i < 30; i++)
		CU_ASSERT_EQUAL(powm(2, i, 1u << 31), 1lu << i);

	CU_ASSERT_EQUAL(powm(2, 16, 1 << 30), 0x10000);
	CU_ASSERT_EQUAL(powm(2, 16, 0xffff), 1);

	mpz_inits(base, mod, result, NULL);
	for (unsigned i = 0; i < 100; i++) {
		unsigned b = (unsigned)get_random_int(0, (1lu << 32) - 1);
		unsigned m = (unsigned)get_random_int(0, (1lu << 32) - 1);
		mpz_set_ui(base, b);
		mpz_set_ui(mod, m);
		for (unsigned e = 0; e < 10000; e++) {
			mpz_powm_ui(result, base, e, mod);
			unsigned long r = mpz_get_ui(result);
			CU_ASSERT_EQUAL(powm(b, e, m), r);
		}
	}
	mpz_clears(base, mod, result, NULL);
}

void test_miller_rabin_primes_int(void)
{
	FILE *fp = fopen("data/primelist.txt", "r");
	unsigned p;
	unsigned long i = 0;

	if (NULL == fp)
		die("data/primelist.txt: %s\n", strerror(errno));

	while (EOF != fscanf(fp, "%u\n", &p))
		large_primes[i++] = p;

	fclose(fp);

	for (i = 0; i < len(large_primes); i++)
		CU_ASSERT_NOT_EQUAL(miller_rabin_int(large_primes[i], 1), composite);
}
