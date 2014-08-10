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

#define max(x,y) (((x) < (y)) ? (y) : (x))
#define log(...) fprintf(stderr, __VA_ARGS__)

extern unsigned long multiplications;

/*
 * How long to run the tests to figure out how many iterations to run to get to
 * at least about one second of runtime.
 */
static const double epsilon = 1e-1;

// Test inputs
static unsigned bits_primes[NUM_PRIMES];
static unsigned bits_composites[NUM_COMPOSITES];
static unsigned bits_mersenne_numbers[NUM_MERSENNE_NUMBERS];
static unsigned bits_mersenne_primes[NUM_MERSENNE_PRIMES];

static mpz_t primes[NUM_PRIMES];
static mpz_t composites[NUM_COMPOSITES];
static mpz_t mersenne_numbers[NUM_MERSENNE_NUMBERS];
static mpz_t mersenne_primes[NUM_MERSENNE_PRIMES];

// Where to write the measurements
static FILE *output;

/*
 * Dummy variable used to make sure calls to the different primality tests are
 * not removed by the optimizer.
 */
static int phantom;

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
	return (stop.tv_sec - start.tv_sec + (stop.tv_nsec - start.tv_nsec) * 1e-9) / its;
}

/*
 * The part that is different for the different tests.
 */
struct GMP {
	static const char name[];
	static const char mode[];
	static Primality check(const mpz_t n) {
		return (Primality)mpz_probab_prime_p(n, 1);
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
		return (Primality)mpz_probab_prime_p(n, 0);
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
 * Figure out how often to repeat each test so that the overall runtime for
 * that particular test is about one second.  This minimizes error in
 * measurement.  For tests where each iteration takes longer than one second, a
 * single iteration will be measured.
 */
template <class T>
static unsigned get_number_of_iterations(const mpz_t num)
{
	unsigned its = 1;
	double duration = 0;
	struct timespec start, stop;

	// Figure out how many iterations to perform
	for (; duration < epsilon; its *= 2) {
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
		for (unsigned j = 0; j < its; j++) {
			phantom += T::check(num);
		}
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
		duration = get_duration(start, stop, 1);

		if (its > (1<<30))
			die("Too many iterations needed, causing integer overflow\n");
	}
	return 1 + (unsigned)(its / duration);
}

/*
 * Compute how many iterations to perform, then measure how long each of these
 * iterations took on average.  Each measurement is repeated NUM_MEASUREMENTS
 * times.  The results are written to output.
 */
template <class T>
static void time_it(const unsigned *bits, const mpz_t *numbers,
	       	const unsigned long first, const unsigned long last,
	       	const char *num_name)
{
	unsigned its = 10000;
	double duration;
	struct timespec start, stop;

	for (unsigned i = first; i <= last; i++) {
		its = get_number_of_iterations<T>(numbers[i]);
		log("will perform %d iterations to reach a runtime of about 1 second\n", its);

		for (unsigned k = 0; k < NUM_MEASUREMENTS; k++) {
			unsigned l = 0;
			multiplications = 0;
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
			for (unsigned j = 0; j < its; j++) {
				l += (T::check(numbers[i]) != composite);
			}
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
			duration = get_duration(start, stop, its);

			fprintf(output, "%d,%lu,%lu,%E,%s,%s,%s,%u,%u\n", bits[i], mpz_popcount(numbers[i]),
					multiplications, duration, T::name, num_name, T::mode, its, l);
			fflush(output);
		}
	}
}

int main(int argc, char *argv[])
{
	unsigned long first_prime, last_prime, first_composite, last_composite,
		      first_mersenne_number, last_mersenne_number, first_mersenne_prime, last_mersenne_prime;

	static const char *primes_name = "primes";
	static const char *composites_name = "composites";
	static const char *mersenne_numbers_name = "Mersenne numbers";
	static const char *mersenne_primes_name = "Mersenne primes";

	init();
	init_int();

	if (argc < 2)
		output = stdout;
	else
		output = fopen(argv[1], "a");

	if (output != stdout && ftell(output) == 0)
		fprintf(output, "Bits,HammingWeight,Multiplications,Time,Algorithm,Set,Mode,Iterations,IsPrime\n");
	fflush(output);

	if (NULL == output)
		die("failed to open output file");

	load_numbers(bits_primes, primes, "primes.txt", NUM_PRIMES) || die("failed to load primes\n");
	load_numbers(bits_composites, composites, "composites.txt", NUM_COMPOSITES) || die("failed to load composites\n");
	load_numbers(bits_mersenne_numbers, mersenne_numbers, "mersenne_numbers.txt", NUM_MERSENNE_NUMBERS) || die("failed to load Mersenne numbers\n");
	load_numbers(bits_mersenne_primes, mersenne_primes, "mersenne_primes.txt", NUM_MERSENNE_PRIMES) || die("failed to load Mersenne primes\n");

	// Figure out how long the precomputation takes, so we can compute how long a single iteration really takes.

	first_prime = 0;
	last_prime = NUM_PRIMES - 1;

	first_composite = 0;
	last_composite = 44;

	first_mersenne_number = 0;
	last_mersenne_number = NUM_MERSENNE_NUMBERS - 1;

	first_mersenne_prime = 0;
	last_mersenne_prime = 20;

#define TIME_IT(alg, set) \
	time_it<alg>(bits_##set##s, set##s, first_##set, last_##set, set##s_name)

	// Apply the different tests to primes
	TIME_IT(GMP_precomputation, prime);
	TIME_IT(MillerRabin_precomputation, prime);
	TIME_IT(Frobenius_precomputation, prime);

	// composites
	TIME_IT(GMP_precomputation, composite);
	TIME_IT(MillerRabin_precomputation, composite);
	TIME_IT(Frobenius_precomputation, composite);

	// some Mersenne numbers
	TIME_IT(GMP_precomputation, mersenne_number);
	TIME_IT(MillerRabin_precomputation, mersenne_number);
	TIME_IT(Frobenius_precomputation, mersenne_number);

	// and 25 Mersenne primes
	TIME_IT(GMP_precomputation, mersenne_prime);
	TIME_IT(MillerRabin_precomputation, mersenne_prime);
	TIME_IT(Frobenius_precomputation, mersenne_prime);

	// Apply the different tests to primes
	TIME_IT(GMP, prime);
	TIME_IT(MillerRabin, prime);
	TIME_IT(Frobenius, prime);

	// composites
	TIME_IT(GMP, composite);
	TIME_IT(MillerRabin, composite);
	TIME_IT(Frobenius, composite);

	// some Mersenne numbers
	TIME_IT(GMP, mersenne_number);
	TIME_IT(MillerRabin, mersenne_number);
	TIME_IT(Frobenius, mersenne_number);

	// and 25 Mersenne primes
	TIME_IT(GMP, mersenne_prime);
	TIME_IT(MillerRabin, mersenne_prime);
	TIME_IT(Frobenius, mersenne_prime);

	// Make sure phantom is not removed by the optimizer.  This might yield
	// a bogus return value which can, however, easily be ignored.
	return phantom % 314159265 != 0;
}