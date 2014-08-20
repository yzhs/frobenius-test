#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

#include "helpers_int.h"
#include "frobenius_int.h"

#define N 1

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
