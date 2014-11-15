/* nextprime.c -- find the smallest prime above 2^k for some large k
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

