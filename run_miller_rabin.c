#include <stdio.h>
#include <stdlib.h>

#include "miller_rabin.h"
#include "helpers.h"

int main()
{
	mpz_t foo, tmp;
	unsigned long upper_bound = (1lu << 32) - 1, dots_every = 1lu << 25;
	unsigned long counter = 0, i;

	mpz_inits(foo, tmp, NULL);
	init();

	for (i = 5; i < upper_bound; i += 2) {
		if (i % dots_every == 1) {  /* we need to check for == 1 since i will always be odd */
			printf(".");
			(void)fflush(stdout);
		}
		mpz_set_ui(foo, i);
		if (miller_rabin(foo, 1))
			counter++;
	}

	printf("\n%lu\n", counter);
	cleanup();
	mpz_clears(foo, tmp, NULL);
}
