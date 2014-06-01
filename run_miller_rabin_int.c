#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "helpers_int.h"
#include "miller_rabin_int.h"

#define N 4

static void *run(void *worker_id_cast_to_void_star)
{
	unsigned long i;
	unsigned long worker_id = (unsigned long)worker_id_cast_to_void_star;
	unsigned long counter = 0;

#ifdef WRITE_PRIMES
	char str[32];
	(void)snprintf(str, 32, "primes_worker_%lu.txt", worker_id);
	FILE *primes = fopen(str, "w");
#endif

	for (i = 5 + 2 * worker_id; i < (1lu << 32) - 1; i += 2 * N) {
		if (i % (1 << 25) == 1) {
			printf(".");
			(void)fflush(stdout);
		}
		if (miller_rabin_int((unsigned)i, 1)) {
			counter++;
#ifdef WRITE_PRIMES
			fprintf(primes, "%lu\n", i);
#endif
		}
	}

#ifdef WRITE_PRIMES
	(void)fclose(primes);
#endif
	return (void *)counter;
}

int main()
{
	unsigned long i;
	unsigned long counter = 0;
	pthread_t threads[N];

	init_int();

	for (i = 0; i < N; i++)
		if (0 != pthread_create(&threads[i], NULL, run, (void *)i))
			die("failed to create thread %lu, exiting\n", i);

	for (i = 0; i < N; i++) {
		unsigned long tmp;
		(void)pthread_join(threads[i], (void **)&tmp);
		counter += tmp;
	}

	printf("\n%lu\n", counter);

	return 0;
}
