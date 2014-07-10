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
	unsigned i = 0, p, length, counter = 0;
	unsigned ints[50];
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
		fflush(stdout);
		mpz_ui_pow_ui(n, 2, ints[i]);
		mpz_add_ui(n, n, 1);

		while (mpz_cmp_ui(n, 50000*50000lu) < 0 || has_small_prime_factor(n) || mpz_probab_prime_p(n, 20)) {
			mpz_add_ui(n, n, 2);
			counter++;
			if (counter >= 10) {
				fprintf(stderr, ".");
				fflush(stderr);
				counter = 0;
			}
		}
		gmp_printf("\t%Zd\n", n);
	}

	mpz_clear(n);
	return 0;
}
