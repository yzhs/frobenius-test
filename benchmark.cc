/* benchmark.cc -- measure runtime of the various algorithms and inputs
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

#include <cstdio>
#include <cstdlib>
#include <ctime>

#include <gmp.h>

extern "C" {
#include "common.h"
#include "helpers.h"
#include "helpers_int.h"

#include "frobenius.h"
#include "frobenius_int.h"
#include "miller_rabin.h"
#include "miller_rabin_int.h"
}

#define NUM_MEASUREMENTS 10

#define NUM_PRIMES 40
#define NUM_COMPOSITES 58
#define NUM_MERSENNE_NUMBERS 292
#define NUM_MERSENNE_PRIMES 25

#define log(...) fprintf(stderr, __VA_ARGS__)

#define TIME_IT(alg, set) \
	if (measure_full && measure_##set && measure_##alg) \
		time_it<alg,set>()

#define TIME_IT_PRECOMP(alg, set) \
	if (measure_prep && measure_##set && measure_##alg) \
		time_it<alg##_precomputation,set>()

extern uint64_t multiplications;

// Where to write the measurements
static FILE *output;

/*
 * Dummy variable used to make sure calls to the different primality tests are
 * not removed by the optimizer.
 */
static int phantom;

struct Primes {
	static const int first;
	static const int last;
	static const char name[];
	static unsigned bits[NUM_PRIMES];
	static mpz_t numbers[NUM_PRIMES];
};
const char Primes::name[] = "primes";
const int Primes::first = 0;
const int Primes::last = NUM_PRIMES-1;
unsigned Primes::bits[NUM_PRIMES];
mpz_t Primes::numbers[NUM_PRIMES];

struct Composites {
	static const int first;
	static const int last;
	static const char name[];
	static unsigned bits[NUM_COMPOSITES];
	static mpz_t numbers[NUM_COMPOSITES];
};
const char Composites::name[] = "composites";
const int Composites::first = 0;
const int Composites::last = 44;
unsigned Composites::bits[NUM_COMPOSITES];
mpz_t Composites::numbers[NUM_COMPOSITES];

struct MersenneNumbers {
	static const int first;
	static const int last;
	static const char name[];
	static unsigned bits[NUM_MERSENNE_NUMBERS];
	static mpz_t numbers[NUM_MERSENNE_NUMBERS];
};
const char MersenneNumbers::name[] = "Mersenne numbers";
const int MersenneNumbers::first = 0;
const int MersenneNumbers::last = NUM_MERSENNE_NUMBERS-1;
unsigned MersenneNumbers::bits[NUM_MERSENNE_NUMBERS];
mpz_t MersenneNumbers::numbers[NUM_MERSENNE_NUMBERS];

struct MersennePrimes {
	static const int first;
	static const int last;
	static const char name[];
	static unsigned bits[NUM_MERSENNE_PRIMES];
	static mpz_t numbers[NUM_MERSENNE_PRIMES];
};
const char MersennePrimes::name[] = "Mersenne primes";
const int MersennePrimes::first = 0;
const int MersennePrimes::last = 20;
unsigned MersennePrimes::bits[NUM_MERSENNE_PRIMES];
mpz_t MersennePrimes::numbers[NUM_MERSENNE_PRIMES];

/*
 * The part that is different for the different tests.
 */
struct GMP {
	static const char name[];
	static const char mode[];
	static Primality check(const mpz_t n) {
		return (Primality) mpz_probab_prime_p(n, 1);
	}
};
const char GMP::name[] = "mpz_probab_prime_p";
const char GMP::mode[] = "full";

struct MillerRabin {
	static const char name[];
	static const char mode[];
	static Primality check(const mpz_t n) {
		return miller_rabin(n, 1);
	}
};
const char MillerRabin::name[] = "Miller-Rabin";
const char MillerRabin::mode[] = "full";

struct Frobenius {
	static const char name[];
	static const char mode[];
	static Primality check(const mpz_t n) {
		return RQFT(n, 1);
	}
};
const char Frobenius::name[] = "Frobenius";
const char Frobenius::mode[] = "full";

struct GMP_precomputation {
	static const char name[];
	static const char mode[];
	static Primality check(const mpz_t n) {
		return (Primality) mpz_probab_prime_p(n, 0);
	}
};
const char GMP_precomputation::name[] = "mpz_probab_prime_p";
const char GMP_precomputation::mode[] = "prep";

struct MillerRabin_precomputation {
	static const char name[];
	static const char mode[];
	static Primality check(const mpz_t n) {
		return miller_rabin(n, 0);
	}
};
const char MillerRabin_precomputation::name[] = "Miller-Rabin";
const char MillerRabin_precomputation::mode[] = "prep";

struct Frobenius_precomputation {
	static const char name[];
	static const char mode[];
	static Primality check(const mpz_t n) {
		return RQFT(n, 0);
	}
};
const char Frobenius_precomputation::name[] = "Frobenius";
const char Frobenius_precomputation::mode[] = "prep";

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

/*
 * Calculate how many seconds a single iteration takes, if it takes from
 * [start] to [stop] to perform [its] iterations.
 */
static double get_duration(struct timespec start, struct timespec stop, unsigned its)
{
	return (stop.tv_sec - start.tv_sec +
		(stop.tv_nsec - start.tv_nsec) * 1e-9) / its;
}

/*
 * Figure out how often to repeat each test so that the overall runtime for
 * that particular test is about one second.  This minimizes error in
 * measurement.  For tests where each iteration takes longer than one second, a
 * single iteration will be measured.
 */
template<class T>
static unsigned get_number_of_iterations(const mpz_t num)
{
	unsigned its = 1;
	double duration = 0;
	struct timespec start, stop;

	// Figure out how many iterations to perform
	for (; duration < 1e-1 /* seconds */; its *= 2) {
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
		for (unsigned j = 0; j < its; j++) {
			phantom += T::check(num);
		}
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
		duration = get_duration(start, stop, 1);

		if (its > (1 << 30))
			die("Too many iterations needed, causing integer overflow\n");
	}
	return 1 + (unsigned)(its / duration);
}

/*
 * Compute how many iterations to perform, then measure how long each of these
 * iterations took on average.  Each measurement is repeated NUM_MEASUREMENTS
 * times.  The results are written to output.
 */
template<class T, class U>
static void time_it()
{
	unsigned its = 10000;
	int64_t epsilon;
	double duration;
	struct timespec start, stop;

	for (unsigned i = U::first; i <= U::last; i++) {
		its = get_number_of_iterations<T>(U::numbers[i]);
		log("will perform %d iterations to reach a runtime of about 1 second\n", its);

		for (unsigned k = 0; k < NUM_MEASUREMENTS; k++) {
			unsigned l = 0;
			multiplications = 0;
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
			for (unsigned j = 0; j < its; j++) {
				l += (T::check(U::numbers[i]) != composite);
			}
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
			duration = get_duration(start, stop, its);

			mpz_ui_pow_ui(tmp0, 2, U::bits[i]);
			mpz_sub(tmp0, U::numbers[i], tmp0);
			epsilon = mpz_get_si(tmp0);

			fprintf(output, "%d,%li,%lu,%lu,%E,%s,%s,%s,%u,%u,%lu\n",
				U::bits[i], epsilon, mpz_popcount(U::numbers[i]),
				multiplications, duration, T::name, U::name,
				T::mode, its, l, mpz_fdiv_ui(U::numbers[i], 4));
			fflush(output);
		}
	}
}

int main(int argc, char *argv[])
{
	init();
	//init_int();

	if (argc < 2)
		output = stdout;
	else
		output = fopen(argv[1], "a");

	if (output == stdout || ftell(output) == 0)
		fprintf(output, "Bits,Epsilon,HammingWeight,Multiplications,Time,Algorithm,Set,Mode,Iterations,IsPrime,nMod4\n");
	fflush(output);

	if (NULL == output)
		die("failed to open output file");

	load_numbers(Primes::bits, Primes::numbers, "primes.txt", NUM_PRIMES)
	    || die("failed to load primes\n");
	load_numbers(Composites::bits, Composites::numbers, "composites.txt", NUM_COMPOSITES)
	    || die("failed to load composites\n");
	load_numbers(MersenneNumbers::bits, MersenneNumbers::numbers, "mersenne_numbers.txt", NUM_MERSENNE_NUMBERS)
	    || die("failed to load Mersenne numbers\n");
	load_numbers(MersennePrimes::bits, MersennePrimes::numbers, "mersenne_primes.txt", NUM_MERSENNE_PRIMES)
	    || die("failed to load Mersenne primes\n");

	// Which tests to run
	static bool measure_full, measure_prep;
	static bool measure_Primes, measure_Composites, measure_MersenneNumbers, measure_MersennePrimes;
	static bool measure_GMP, measure_MillerRabin, measure_Frobenius;

	measure_full = true;
	measure_prep = true;

	measure_Primes = true;
	measure_Composites = true;
	measure_MersenneNumbers = true;
	measure_MersennePrimes = true;

	measure_GMP = true;
	measure_MillerRabin = true;
	measure_Frobenius = true;

	/*
	 * Precomputation
	 */

	// Apply the different tests to primes
	TIME_IT_PRECOMP(GMP, Primes);
	TIME_IT_PRECOMP(MillerRabin, Primes);
	TIME_IT_PRECOMP(Frobenius, Primes);

	// composites
	TIME_IT_PRECOMP(GMP, Composites);
	TIME_IT_PRECOMP(MillerRabin, Composites);
	TIME_IT_PRECOMP(Frobenius, Composites);

	// some Mersenne numbers
	TIME_IT_PRECOMP(GMP, MersenneNumbers);
	// Finished 2^317-1
	TIME_IT_PRECOMP(MillerRabin, MersenneNumbers);
	TIME_IT_PRECOMP(Frobenius, MersenneNumbers);

	// and 25 Mersenne primes
	TIME_IT_PRECOMP(GMP, MersennePrimes);
	TIME_IT_PRECOMP(MillerRabin, MersennePrimes);
	TIME_IT_PRECOMP(Frobenius, MersennePrimes);

	/*
	 * Full test
	 */

	// Apply the different tests to primes
	TIME_IT(GMP, Primes);
	TIME_IT(MillerRabin, Primes);
	TIME_IT(Frobenius, Primes);

	// composites
	TIME_IT(GMP, Composites);
	TIME_IT(MillerRabin, Composites);
	TIME_IT(Frobenius, Composites);

	// some Mersenne numbers
	TIME_IT(GMP, MersenneNumbers);
	TIME_IT(MillerRabin, MersenneNumbers);
	TIME_IT(Frobenius, MersenneNumbers);

	// and 25 Mersenne primes
	TIME_IT(GMP, MersennePrimes);
	TIME_IT(MillerRabin, MersennePrimes);
	TIME_IT(Frobenius, MersennePrimes);

	cleanup();
	// Make sure phantom is not removed by the optimizer.  This might yield
	// a bogus return value which can, however, easily be ignored.
	return phantom % 314159265 != 0;
}
