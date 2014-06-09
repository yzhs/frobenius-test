#define _POSIX_C_SOURCE 201406L
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <gmp.h>

#include "common.h"
#include "helpers.h"
#include "helpers_int.h"

#include "frobenius.h"
#include "frobenius_int.h"
#include "miller_rabin.h"
#include "miller_rabin_int.h"

#define NUM_PRIMES 52
#define NUM_COMPOSITES 52

#define max(x,y) (((x) < (y)) ? (y) : (x))

static unsigned iterations = 100000;


static unsigned load_numbers(unsigned num_bits[], mpz_t nums[], const char file[])
{
	unsigned p, i = 0;
	FILE *fp = fopen(file, "r");
	mpz_t tmp;

	if (NULL == fp)
		exit(1);

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


int main(int argc, char *argv[])
{
	unsigned bits[max(NUM_PRIMES, NUM_COMPOSITES)];
	unsigned its;
	int phantom = 0;

	mpz_t primes[NUM_PRIMES];
	mpz_t composites[NUM_COMPOSITES];

	FILE *output;

	double factor;
	double duration = 0;

	double epsilon = 1;

	mpz_t num;

	mpz_init(num);
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

	// mpz_propab_prime_p (GMP) applied to primes
	factor = 1;
	its = 1;
	miller_rabin(primes[0], iterations);
	printf("%d iterations, each repeated %d times took %E seconds\n", iterations, its, duration);
	// Figure out how many iterations to perform
	for (; duration < epsilon; its *= 2) {
		struct timespec start, stop;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
		for (unsigned j = 0; j < its; j++) {
			phantom += mpz_probab_prime_p(primes[0], iterations) + (int)j;
		}
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
		duration = (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec) * 1e-9;
		printf("%u iterations, each repeated %u times took %E seconds\n", iterations, its, duration);
	}
	printf("performing %d iterations took %f seconds\n", its, duration);
	its = (unsigned)(its / duration);
	printf("will perform %d iterations to reach a duration of about 1 second\n", its);
	fflush(stdout);
	for (unsigned i = 0; i < len(primes); i++) {
		struct timespec start, stop;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &start);
		for (unsigned j = 0; j < its; j++) {
			phantom += mpz_probab_prime_p(primes[i], iterations) + (int)j;
		}
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &stop);
		duration = (stop.tv_sec - start.tv_sec) + (stop.tv_nsec - start.tv_nsec) * 1e-9;
		fprintf(output, "%d\t%E\tprobab_prime_p (prime)\n", bits[i], factor*duration);
		fflush(output);
		if (duration > 60 && iterations >= 10) {
			fprintf(stderr, "decreasing number of iterations by a factor of 10\n");
				iterations /= 10;
				factor *= 10;
		}
	}

	printf("%d\n", phantom);

#if 0
	// mpz_propab_prime_p (GMP) applied to composites
	factor = 1;
	miller_rabin(composites[15], iterations);
	for (unsigned i = 0; i < len(composites); i++) {
		fprintf(output, "%d	", bits[i]);
		time_it(mpz_probab_prime_p(composites[i], iterations), duration);
		fprintf(output, "%E	probab_prime_p (composite)\n", factor*duration);
		fflush(output);
		adjust_iterations();
	}


	// Miller-Rabin (GMP) applied to primes
	factor = 1;
	miller_rabin(primes[15], iterations);
	for (unsigned i = 0; i < len(primes); i++) {
		fprintf(output, "%d	", bits[i]);
		time_it(miller_rabin(primes[i], iterations), duration);
		fprintf(output, "%E	Miller-Rabin (prime)\n", factor*duration);
		fflush(output);
		adjust_iterations();
	}

	// Miller-Rabin (GMP) applied to composites
	factor = 1;
	miller_rabin(composites[15], iterations);
	for (unsigned i = 0; i < len(composites); i++) {
		fprintf(output, "%d	", bits[i]);
		time_it(miller_rabin(composites[i], iterations), duration);
		fprintf(output, "%E	Miller-Rabin (composite)\n", factor*duration);
		fflush(output);
		adjust_iterations();
	}

	// Frobenius (GMP) applied to primes
	factor = 1;
	RQFT(primes[15], iterations);
	for (unsigned i = 0; i < len(primes); i++) {
		fprintf(output, "%d	", bits[i]);
		time_it(RQFT(primes[i], iterations), duration);
		fprintf(output, "%E	Frobenius (prime)\n", factor*duration);
		fflush(output);
		adjust_iterations();
	}

	// Frobenius (GMP) applied to composites
	factor = 1;
	RQFT(composites[15], iterations);
	for (unsigned i = 0; i < len(composites); i++) {
		fprintf(output, "%d	", bits[i]);
		time_it(RQFT(composites[i], iterations), duration);
		fprintf(output, "%E	Frobenius (composite)\n", factor*duration);
		fflush(output);
		adjust_iterations();
	}
#endif

	mpz_clear(num);
	return 0;
}
