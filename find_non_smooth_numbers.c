#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include "small_primes.h"

#define len(a) (sizeof(a)/sizeof(a[0]))

static int has_small_prime_factor(mpz_t n)
{
	for (unsigned i = 0; i < len(prime_list); i++)
		if (mpz_divisible_ui_p(n, prime_list[i]))
			return 1;
	return 0;
}

int main()
{
	unsigned i = 0, p, length;
	unsigned ints[27];
	FILE *fp = fopen("composites.txt", "r");
	mpz_t n;

	if (NULL == fp)
		return 1;
	while (EOF != fscanf(fp, "%u\t\n", &p)) {
		ints[i++] = p;
	}
	fclose(fp);
	length = i;

	mpz_init(n);

	for (i = 0; i < length; i++) {
		printf("%u", ints[i]);
		mpz_ui_pow_ui(n, 2, ints[i]);
		mpz_add_ui(n, n, 1);

		while (mpz_cmp_ui(n, 50000*50000lu) < 0 || has_small_prime_factor(n) || mpz_probab_prime_p(n, 50)) {
			mpz_add_ui(n, n, 2);
		}
		gmp_printf("\t%Zd\n", n);
	}

	mpz_clear(n);
	return 0;
}
