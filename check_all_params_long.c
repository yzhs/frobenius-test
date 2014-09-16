#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include "frobenius.h"
#include "helpers.h"

/*
 * Run the randomized quadratic frobenius test to check whether [n] is a a
 * prime.  The Parameter [k] determines how many times the test will be run at
 * most.  If the test returns "composite", it will not be run again.
 */
int main()
{
	Primality result;

	uint64_t i, j, N = 251*257;
	mpz_t n, b, c, bb4c, tmp;
	mpz_init_set_ui(n, N); // The largest two prime factors such that their product is at most 2¹⁶.
	mpz_inits(b, c, bb4c, tmp, NULL);

	uint64_t valid_pairs = 0, false_positives = 0;

	for (i = N-1; i < N; i++) {
		mpz_set_ui(c, i);
		mpz_sub(tmp, n, c);
		if (mpz_jacobi(tmp, n) != 1)
			continue;
		for (j = 0; j < N; j++) {
			mpz_set_ui(b, j);
			mpz_set_ui(bb4c, (j*j+4*i)%N);
			if (mpz_jacobi(bb4c, n) != -1)
				continue;
			valid_pairs++;

			result = steps_3_4_5(n, b, c);
			if (result != composite) {
				false_positives++;
				printf("#");
				fflush(stdout);
				gmp_fprintf(stderr, "(GMP) Found a false positive: n = %Zd, b = %Zd, c = %Zd\n", n, b, c);
			}
		}
		if (i % 100 == 1) {
			printf(".");
			fflush(stdout);
		}
	}

	printf("Checked a total of %lu parameter pairs.  %lu of these were valid parameters.\n", N*N, valid_pairs);
	printf("A total number of %lu false positives were found\n", false_positives);

	return 0;
}
