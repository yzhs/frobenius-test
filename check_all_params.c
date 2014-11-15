/* check_all_params.c -- test all params for the QFT on small numbers
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

#include "frobenius_int.h"
#include "helpers_int.h"

/*
 * Run the randomized quadratic frobenius test to check whether [n] is a a
 * prime.  The Parameter [k] determines how many times the test will be run at
 * most.  If the test returns "composite", it will not be run again.
 */
int main()
{
	Primality result;

	uint64_t n = 251*257; // The largest two prime factors such that their product is at most 2¹⁶.
	uint64_t b = 0, c = 0;
	uint64_t bb4c_int;

	uint64_t valid_pairs = 0, false_positives = 0;

	for (c = n-1; c < n; c++) {
		if (jacobi(n - c, n) != 1)
			continue;
		for (b = 0; b < n; b++) {
			bb4c_int = ((b * b) % n + c * 4) % n;
			if (jacobi(bb4c_int, n) != -1)
				continue;
			valid_pairs++;

			result = steps_3_4_5_int(n, b, c);
			if (result != composite) {
				false_positives++;
				printf("#");
				fflush(stdout);
				fprintf(stderr, "Found a false positive: n = %lu, b = %lu, c = %lu\n", n, b, c);
			}
		}
		if (c % 100 == 1) {
			printf(".");
			fflush(stdout);
		}
	}

	printf("Checked a total of %lu parameter pairs.  %lu of these were valid parameters.\n", n*n, valid_pairs);
	printf("A total number of %lu false positives were found\n", false_positives);

	return 0;
}
