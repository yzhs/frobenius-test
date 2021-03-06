/* test/test_miller_rabin_long.c -- tests for the Miller-Rabin test (GMP)
 *
 * Copyright 2014 by Colin Benner <colin-software@yzhs.de>
 *
 * This file is part of frobenius-test.
 *
 * frobenius-test is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * frobenius-test is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with frobenius-test.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <CUnit/CUnit.h>
#include "test_miller_rabin_long.h"

#include "../miller_rabin.c"

static unsigned large_primes[3069262];
static mpz_t composites[1013];

static const uint64_t n = 1000;

void mr_some_numbers(void)
{
	mpz_t tmp;

	mpz_init(tmp);

	uint64_t primes[] = { 5, 7, 11, 13, 17, 19, 23, 29, 31 };

	for (uint64_t i = 0; i < len(primes); i++) {
		mpz_set_ui(tmp, primes[i]);
		CU_ASSERT_NOT_EQUAL(miller_rabin(tmp, 1), composite);
	}

	// 611879² * 611957⁴ = 52506700005424014690584902271185441
	mpz_set_str(tmp, "52506700005424014690584902271185441", 10);
	CU_ASSERT_EQUAL(miller_rabin(tmp, 1), composite);

	mpz_pow_ui(tmp, tmp, 9);
	CU_ASSERT_EQUAL(miller_rabin(tmp, 1), composite);

	mpz_set_str(tmp, "1235790412356789098765432827498274392743929834792843282734279348239482349", 10);
	CU_ASSERT_EQUAL(miller_rabin(tmp, 1), composite);
}

void mr_primes(void)
{
	FILE *fp = fopen(TEST_DATA_PATH "primelist.txt", "r");
	unsigned p;
	uint64_t i = 0;
	mpz_t prime;

	if (NULL == fp)
		die(TEST_DATA_PATH "primelist.txt: %s\n", strerror(errno));

	while (EOF != fscanf(fp, "%u\n", &p))
		large_primes[i++] = p;

	fclose(fp);

	mpz_init(prime);
	for (int k = 0; k < 10000; k++) {
		i = (uint64_t)rand() % len(large_primes);
		mpz_set_ui(prime, large_primes[i]);
		CU_ASSERT_NOT_EQUAL_FATAL(miller_rabin(prime, 1), composite);
		if (i % 1000 == 999)
			fprintf(stderr, ".");
	}
	mpz_clear(prime);
}

void mr_composites(void)
{
	FILE *fp = fopen(TEST_DATA_PATH "composites.txt", "r");
	uint64_t i;

	for (i = 0; i < len(composites); i++)
		gmp_fscanf(fp, "%Zd\n", composites[i]);

	fclose(fp);

	for (i = 0; i < len(composites); i++)
		CU_ASSERT_EQUAL(miller_rabin(composites[i], 10), composite);
}

void mr_composites2(void)
{
	unsigned i, c, k;
	mpz_t tmp;

	mpz_init(tmp);

	for (i = 0; i < n; i++) {
		for (k = 0; k < ((large_primes[i + 1] - large_primes[i]) / 10 | 1); k++) {
			Primality res;
			c = randint(large_primes[i] + 1, large_primes[i + 1] - 1) | 1;
			if (c == large_primes[i + 1])
				continue;
			mpz_set_ui(tmp, c);
			res = miller_rabin(tmp, 10);
			if (res != composite)
				fprintf(stderr, "%i was not correctly recognized as composite\n", c);
			CU_ASSERT_EQUAL_FATAL(res, composite);
		}
	}
	mpz_clear(tmp);
}

extern Primality miller_rabin_int(unsigned n, int k);

void mr_both(void)
{
	mpz_t tmp;

	mpz_init(tmp);

	for (unsigned k = 5; k < 10000; k++) {
		mpz_set_ui(tmp, k);
		CU_ASSERT_EQUAL(miller_rabin(tmp, 100), miller_rabin_int(k, 100));
	}
}
