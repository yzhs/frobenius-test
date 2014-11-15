/* add_epsilon_mod_4.c -- add a column with n-2^k to a cvs file
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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "common.h"

#define NUM_PRIMES 40
#define NUM_COMPOSITES 58
#define NUM_MERSENNE_NUMBERS 292
#define NUM_MERSENNE_PRIMES 25

// Test inputs
static unsigned bits_primes[NUM_PRIMES];
static unsigned bits_composites[NUM_COMPOSITES];
static unsigned bits_mersenne_numbers[NUM_MERSENNE_NUMBERS];
static unsigned bits_mersenne_primes[NUM_MERSENNE_PRIMES];

static mpz_t primes[NUM_PRIMES];
static mpz_t composites[NUM_COMPOSITES];
static mpz_t mersenne_numbers[NUM_MERSENNE_NUMBERS];
static mpz_t mersenne_primes[NUM_MERSENNE_PRIMES];

/*
 * Read a set of numbers for testing
 */
static unsigned load_numbers(unsigned num_bits[], mpz_t nums[], const char file[], size_t length)
{
	unsigned p, i = 0;
	FILE *fp = fopen(file, "r");
	mpz_t tmp;

	if (NULL == fp)
		return 0;

	mpz_init(tmp);

	while (EOF != gmp_fscanf(fp, "%u\t%Zd\n", &p, tmp)) {
		mpz_init(nums[i]);
		mpz_set(nums[i], tmp);
		if (NULL != num_bits)
			num_bits[i] = p;
		i++;
		if (i >= length)
			break;
	}

	fclose(fp);
	mpz_clear(tmp);

	return i;
}

static int count_commas(char *s)
{
	int counter = 0;

	for (char *c = s; *c != '\0'; c++) {
		counter += (*c == ',');
	}

	return counter;
}

/*
 * Find the number used in a test given which set the number belongs two and
 * the logarithm of the closest power of 2.
 */
static void get_num(mpz_t num, mpz_t numbers[], unsigned bits_numbers[], uint64_t bits)
{
	int i;
	for (i = 0; bits_numbers[i] < bits; i++);
	assert(bits_numbers[i] == bits);
	mpz_set(num, numbers[i]);
}

int main()
{
	FILE *input, *output;
	mpz_t n, tmp;

	load_numbers(bits_primes, primes, "primes.txt", NUM_PRIMES)
	    || die("failed to load primes\n");
	load_numbers(bits_composites, composites, "composites.txt", NUM_COMPOSITES)
	    || die("failed to load composites\n");
	load_numbers(bits_mersenne_numbers, mersenne_numbers, "mersenne_numbers.txt", NUM_MERSENNE_NUMBERS)
	    || die("failed to load Mersenne numbers\n");
	load_numbers(bits_mersenne_primes, mersenne_primes, "mersenne_primes.txt", NUM_MERSENNE_PRIMES)
	    || die("failed to load Mersenne primes\n");

	input = fopen("timings_20140708.csv", "r");
	if (NULL == input)
		die("failed to open input file");

	output = fopen("timings_new.csv", "w");
	if (NULL == output)
		die("failed to open output file");

	mpz_inits(n, tmp, NULL);

	/*
	 * Now we have to read the file with the measurements.  For every of
	 * its lines, we add the missing measurements, if there are any.  The
	 * first line containing the column headers is skipped.
	 */
	long read;
	char *line = NULL;
	size_t len;
	while ((read = getline(&line, &len, input)) != -1) {
		// Get rid of the final newline character.
		line[strlen(line)-1] = '\0';
		uint64_t bits;
		if (count_commas(line) == 10) {
			fprintf(output, "%s\n", line);
			continue;
		}
		sscanf(line, "%lu,", &bits);
		mpz_ui_pow_ui(tmp, 2, bits);
		if (NULL != strstr(line, ",primes,"))
			get_num(n, primes, bits_primes, bits);
		else if (NULL != strstr(line, ",composites,"))
			get_num(n, composites, bits_composites, bits);
		else if (NULL != strstr(line, ",Mersenne numbers,"))
			get_num(n, mersenne_numbers, bits_mersenne_numbers, bits);
		else if (NULL != strstr(line, ",composites,"))
			get_num(n, mersenne_primes, bits_mersenne_primes, bits);
		// We do not have to do anything otherwise, because it is the
		// line containing the column headers.
		mpz_sub(tmp, n, tmp);
		fprintf(output, "%lu,%ld%s,%lu\n", bits, mpz_get_si(tmp), strchr(line, ','), mpz_fdiv_ui(n, 4));
	}
	free(line);

	mpz_clears(n, tmp, NULL);
	fclose(input);
	fclose(output);

	return 0;
}
