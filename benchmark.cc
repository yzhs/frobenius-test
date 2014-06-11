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

#define NUM_PRIMES 52
#define NUM_COMPOSITES 52

#define max(x,y) (((x) < (y)) ? (y) : (x))
#define log(...) fprintf(stderr, __VA_ARGS__)

/*
 * How long to run the tests to figure out how many iterations to run to get to
 * at least about one second of runtime.
 */
static const double epsilon = 1e-1;

// Test inputs
static unsigned bits[max(NUM_PRIMES, NUM_COMPOSITES)];
static mpz_t primes[NUM_PRIMES];
static mpz_t composites[NUM_COMPOSITES];

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
static unsigned load_numbers(unsigned num_bits[], mpz_t nums[], const char file[])
{
	unsigned p, i = 0;
	FILE *fp = fopen(file, "r");
	mpz_t tmp;

	if (NULL == fp)
		die("Could not open file %s for reading\n", file);

	mpz_init(tmp);

	while (EOF != gmp_fscanf(fp, "%u\t%Zd\n", &p, tmp)) {
		mpz_init(nums[i]);
		mpz_set(nums[i], tmp);
		if (NULL != num_bits)
			num_bits[i] = p;
		i++;
	}

	fclose(fp);
	mpz_clear(tmp);

	return i;
}

/*
 * Given the times when the computation started and stopped, as well as the
 * number of iterations performed, compute the average time per iteration.
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
	static const unsigned iterations = 368; // 4^368 ~= 7710^57
	static Primality check(const mpz_t n) {
		return (Primality)mpz_probab_prime_p(n, iterations);
	}
};
const char GMP::name[] = "mpz_probab_prime_p";

struct MillerRabin {
	static const char name[];
	static const unsigned iterations = 368; // 4^368 ~= 7710^57
	static Primality check(const mpz_t n) {
		return miller_rabin(n, iterations);
	}
};
const char MillerRabin::name[] = "Miller-Rabin";

struct Frobenius {
	static const char name[];
	static const unsigned iterations = 57; // 4^368 ~= 7710^57
	static Primality check(const mpz_t n) {
		return RQFT(n, iterations);
	}
};
const char Frobenius::name[] = "Frobenius";

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
	for (; duration < epsilon; its *= 3) {
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
		for (unsigned j = 0; j < its; j++) {
			phantom += T::check(num);
		}
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
		duration = get_duration(start, stop, 1);

		if (its > (1<<30))
			die("Too many iterations needed causing integer overflow\n");
	}
	return 1 + (unsigned)(its / duration);
}

/*
 * Compute how many iterations to perform, then measure how long each of these
 * iterations took on average.  Each measurement is repeated NUM_MEASUREMENTS
 * times.  The results are written to output.
 */
template <class T>
static void time_it(const mpz_t *numbers, const unsigned long len, const char *num_name)
{
	unsigned its;
	double duration;
	struct timespec start, stop;

	for (unsigned i = 0; i < len; i++) {
		its = get_number_of_iterations<T>(numbers[i]);
		log("will perform %d iterations to reach a runtime of about 1 second\n", its);

		for (unsigned k = 0; k < NUM_MEASUREMENTS; k++) {

			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
			for (unsigned j = 0; j < its; j++) {
				phantom += T::check(numbers[i]);
			}
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
			duration = ((stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec) * 1e-9) / its;

			fprintf(output, "%d\t%E\t%s (%s)\n", bits[i], duration, T::name, num_name);
			fflush(output);
		}
	}
}

int main(int argc, char *argv[])
{
	static const char *primes_name = "primes";
	static const char *composites_name = "composites";

	init();
	init_int();

	if (argc < 2)
		output = stdout;
	else
		output = fopen(argv[1], "w");
	if (NULL == output)
		die("failed to open output file");

	load_numbers(bits, primes, "primes.txt") || die("failed to load primes\n");
	load_numbers(NULL, composites, "composites.txt") || die("failed to load composites\n");;

	fprintf(output, "Bits,Time,Algorithm\n");
	fflush(output);
	fflush(stdout);

	time_it<GMP>        (primes, len(primes), primes_name);
	time_it<MillerRabin>(primes, len(primes), primes_name);
	time_it<Frobenius>  (primes, len(primes), primes_name);

	time_it<GMP>        (composites, len(composites), composites_name);
	time_it<MillerRabin>(composites, len(composites), composites_name);
	time_it<Frobenius>  (composites, len(composites), composites_name);

	// Make sure phantom is not removed by the optimizer.  This might yield
	// a bogus return value which can, however, easily be ignored.
	return phantom % 1312313 != 0;
}
