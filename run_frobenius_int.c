#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "helpers_int.h"
#include "frobenius_int.h"

#define N 1

static void *run(void *worker_id_cast_to_void_star)
{
	uint64_t upper_bound = /*1lu << 32*/50000lu*50000+(1lu<<27), dots_every = 1lu << 25;
	uint64_t counter = 0;
	unsigned worker_id = (unsigned)(long)worker_id_cast_to_void_star;

#ifdef WRITE_PRIMES
	char str[32];
	FILE *primes;
	(void)snprintf(str, 32, "primes_worker_%lu.txt", worker_id);
	primes = fopen(str, "a");
#endif

	for (unsigned i = B * B + 1 + 2 * worker_id; i < upper_bound; i += 2 * N) {
		if ((i & (dots_every - 1)) == 1) {  /* we need to check for == 1 since i will always be odd */
			printf(".");
			(void)fflush(stdout);
		}

		if (RQFT_int(i, 1)) {
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
	unsigned i;
	unsigned counter = 0;
	pthread_t threads[N];

	init_int();

	for (i = 0; i < N; i++)
		if (0 != pthread_create(&threads[i], NULL, run, (void *)(uint64_t)i))
			die("failed to create thread %u, exiting\n", i);

	for (i = 0; i < N; i++) {
		uint64_t tmp;
		(void)pthread_join(threads[i], (void **)&tmp);
		counter += tmp;
	}

	printf("\n%u\n", counter);

	return 0;
}
