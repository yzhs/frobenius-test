/* check_all_small_numbers.c -- test the RQFT on a range of small numbers.
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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "frobenius_int.h"
#include "miller_rabin_int.h"
#include "helpers_int.h"

/*
 * Uncomment the following line to use the smallest possible c instead of
 * always choosing c=-4.
 */
//#define SMALL_C

/*
 * Given an integer 0 ≤ n < 2³², find the largest integer r such that r² ≤ n.
 */
uint64_t square_p(const uint64_t n)
{
	uint64_t root = (uint64_t)(sqrt((double)n));

	if (n == 0)
		return 0;
	if (root * root > n)
		root = n / root;
	while ((root + 1) * (root + 1) <= n)
		root++;
	return root * root == n;
}


/*
 * Run the randomized quadratic frobenius test to check whether [n] is a a
 * prime.  The Parameter [k] determines how many times the test will be run at
 * most.  If the test returns "composite", it will not be run again.
 */
int main(int argc, char *argv[])
{
	Primality result;
	uint64_t n, b = 0, c, bb4c_int;
	uint64_t false_positives = 0;

	uint64_t lower_bound = 5, upper_bound = 1000000;

	// Parse command line arguments giving an upper bound or a lower bound
	// and an upper bound.
	if (argc > 2)
		lower_bound = strtoul(argv[1], NULL, 10);
	if (argc > 1)
		upper_bound = strtoul(argv[argc > 2 ? 2 : 1], NULL, 10);
	printf("%lu to %lu\n", lower_bound, upper_bound);

	for (n = lower_bound | 1; n <= upper_bound; n+=2) {
		// There is no valid pair (b,c) if n is a perfect square.
		// Hence, we skip those n.  Since we want to see how many false
		// positives there are, we don't want to apply the test to
		// primes n.
		if (square_p(n) || miller_rabin_int(n, 100) != composite)
			continue;
#ifdef SMALL_C
		// Iterate over all possible 2≤c<n until we find one such that
		// -c is a square mod n.
		for (c = 2; c < n; c++) {
			if (jacobi(n - c, n) != 1)
				continue;
#else
			// -c is always a square if we set
			c = n - 4;
#endif
			// Find the smallest suitable parameter 1≤b<n.
			for (b = 1; b < n; b++) {
				bb4c_int = ((b * b) % n + 4 * c) % n;
				if (jacobi(bb4c_int, n) == -1)
					break;
			}
#ifdef SMALL_C
			break;
		}
#endif

		// Perform the QFT with parameters (b,c) without checking
		// whether n has a small primefactor.  That would obviously
		// always be true for 5≤n≤1000000<B^2.  Even if we go up to
		// n<2^32 (the largest value this implementation can handle),
		// we would otherwise only test a few numbers.
		result = steps_3_4_5_int(n, b, c);
		if (result != composite) {
			false_positives++;
			fflush(stdout);
			printf("Found a false positive: n = %lu, b = %lu, c = %lu\n", n, b, c);
		}
		if (n % 10000000 == 1) {
			fprintf(stderr, "%lu\n", n);
			fflush(stderr);
		}
	}


	printf("\nA total number of %lu false positives were found among the numbers %lu,...,%lu\n",
	       false_positives, lower_bound, upper_bound);

	return 0;
}
