#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

int main(int argc, char *argv[])
{
	unsigned bits;
	int counter = 0;
	mpz_t n;

	if (argc < 2) {
		fprintf(stderr, "Usage: %s <bits>\n", argv[0]);
		return 1;
	}

	bits = (unsigned)atoi(argv[1]);

	mpz_init(n);
	mpz_ui_pow_ui(n, 2, bits);
	mpz_add_ui(n, n, 1);

	while (!mpz_probab_prime_p(n, 25)) {
		mpz_add_ui(n, n, 2);
		counter++;
		if (counter % 200 == 0)
			fprintf(stderr, "tried %i numbers so far\n", counter);
	}

	gmp_printf("%Zd\n", n);
	return 0;
}

