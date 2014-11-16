/* check_all_small_numbers_long.c -- test the RQFT on a range of small numbers (GMP version)
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

#include <stdio.h>
#include <stdlib.h>

#include <gmp.h>

#include "frobenius.h"
#include "helpers.h"
#include "helpers_int.h"

/*
 * Uncomment the following to always select the smallest possible parameter c
 * instead of using c = n-4.
 */
//#define SMALL_C

/*
 * Run the randomized quadratic frobenius test to check whether [n] is a a
 * prime.  The Parameter [k] determines how many times the test will be run at
 * most.  If the test returns "composite", it will not be run again.
 */
int main(int argc, char *argv[])
{
	Primality result;
#define VARS n, b, c, bb4c
	mpz_t VARS;
	uint64_t n_, b_ = 0, c_;
	uint64_t false_positives = 0;

	uint64_t lower_bound = 5, upper_bound = 1000000;

	mpz_inits(VARS, NULL);

	if (argc > 2)
		lower_bound = strtoul(argv[1], NULL, 10);
	if (argc > 1)
		upper_bound = strtoul(argv[argc > 2 ? 2 : 1], NULL, 10);
	printf("%lu to %lu\n", lower_bound, upper_bound);

	for (n_ = lower_bound | 1; n_ <= upper_bound; n_+=2) {
		mpz_set_ui(n, n_);
		if (mpz_perfect_square_p(n))
			continue;
#ifdef SMALL_C
		for (c_ = 2; c_ < n_; c_++) {
			if (jacobi(n_ - c_, n_) != 1)
				continue;
#else
			c_ = n_ - 4;
			mpz_sub_ui(c, n, 4);
#endif
			for (b_ = 1; b_ < n_; b_++) {
				mpz_mul(bb4c, b, b);
				mpz_addmul_ui(bb4c, c, 4);
				mpz_mod(bb4c, bb4c, n);
				if (mpz_jacobi(bb4c, n) == -1)
					break;
				if (b_ % 1000 == 999) {
				       fprintf(stderr, "Could not find a valid parameter pair (b, c)\n");
				       goto next_n;
				}
			}
#ifdef SMALL_C
			break;
		}
#endif

		result = steps_3_4_5(n, b, c);
		if (result != composite) {
			false_positives++;
			fflush(stdout);
			printf("Found a false positive: n = %lu, b = %lu, c = %lu\n", n_, b_, c_);
		}
		if (n_ % 10000000 == 1) {
			fprintf(stderr, ".");
			fflush(stderr);
		}
next_n:
		;
	}


	printf("\nA total number of %lu false positives were found among the numbers %lu,...,%lu\n",
	       false_positives, lower_bound, upper_bound);
	mpz_clears(VARS, NULL);

	return 0;
}
