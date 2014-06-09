#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "common.h"

int main(int argc, char *argv[])
{
	unsigned bits;
	mpz_t n;

	if (argc < 2)
		die("Usage: %s <bits>\n", argv[0]);

	bits = (unsigned)atoi(argv[1]);

	mpz_init(n);
	mpz_ui_pow_ui(n, 2, bits);
	mpz_add_ui(n, n, 1);
	while (!mpz_probab_prime_p(n, 50))
		mpz_add_ui(n, n, 2);

	gmp_printf("%Zd\n", n);
	return 0;
}

