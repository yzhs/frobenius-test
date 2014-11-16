/* check_all_params_long.c -- test all params for the QFT on small numbers (GMP version)
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

/*
 * Run the randomized quadratic frobenius test to check whether [n] is a a
 * prime.  The Parameter [k] determines how many times the test will be run at
 * most.  If the test returns "composite", it will not be run again.
 */
int main()
{
	Primality result;
	uint64_t n_, b_ = 0, c_ = 0, bb4c_int;
	uint64_t valid_pairs = 0, false_positives = 0;
	mpz_t n, b, c, bb4c, tmp;
	mpz_inits(n, b, c, bb4c, tmp, NULL);

	for (n_ = 3; n_ < 10000; n_+=2) {
		// Trial division is enough, so it does not matter, whether we call RQFT_int(n, 0) or RQFT_int(n, 10000).
		mpz_set_ui(n, n_);
		if (RQFT(n, 0) == prime || mpz_perfect_square_p(n))
			continue;
		valid_pairs = false_positives = 0;

		for (c_ = 1; c_ < n_; c_++) {
			mpz_set_ui(c, c_);
			mpz_sub(tmp, n, c);
			if (mpz_jacobi(tmp, n) != 1)
				continue;
			for (b_ = 0; b_ < n_; b_++) {
				mpz_set_ui(b, b_);
				bb4c_int = ((b_ * b_) % n_ + c_ * 4) % n_;
				mpz_set_ui(bb4c, bb4c_int);
				if (mpz_jacobi(bb4c, n) != -1)
					continue;
				valid_pairs++;

				result = steps_3_4_5(n, b, c);
				if (result != composite) {
					false_positives++;
					//printf("#");
					fflush(stdout);
					fprintf(stdout, "n=%5lu Found a false positive: n = %5lu, b = %5lu, c = %5lu\n",
					        n_, n_, b_, c_);
				}
			}
			//if (c_ % 100 == 1) {
			//	printf(".");
			//	fflush(stdout);
			//}
		}

		printf("n=%5lu Checked a total of %lu parameter pairs modulo %lu.  %lu of these were valid parameters.\n",
		       n_, n_*n_, n_, valid_pairs);
		printf("n=%5lu A total number of %lu false positives were found\n", n_, false_positives);
	}

	mpz_clears(n, b, c, bb4c, tmp, NULL);
	return 0;
}
