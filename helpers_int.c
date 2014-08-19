#include "helpers_int.h"
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

/*
 * Given an integer 0 ≤ n < 2³², find the largest integer r such that r² ≤ n.
 */
uint64_t int_sqrt(const uint64_t n)
{
	uint64_t root = (uint64_t)(sqrt((double)n));

	if (n == 0)
		return 0;
	if (root * root > n)
		root = n / root;
	while ((root + 1) * (root + 1) <= n)
		root++;
	return root;
}

/*
 * Calculate the greatest common divisor of a and b using the Euclidean
 * algorithm.
 */
uint64_t gcd(uint64_t a, uint64_t b)
{
	uint64_t t;

	while (b != 0) {
		t = b;
		b = a % b;
		a = t;
	}
	return a;
}

/*
 * Use the integer square root function to figure out whether a number is a
 * perfect square.
 */
bool is_square(const uint64_t n)
{
	uint64_t sqrt = int_sqrt(n);

	return sqrt * sqrt == n;
}

/*
 * Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd.
 */
void split_int(uint64_t *s, uint64_t *d, const uint64_t n)
{
	*s = 0;
	*d = n - 1;

	while (even(*d)) {
		(*s)++;
		*d /= 2;
	}
}

/*
 * This function generates a random integer between 2 and n-2.  This function
 * will fail if n is equal to 3 and produce weird results for smaller ns.
 */
uint64_t get_random_int(const uint64_t low, const uint64_t high)
{
	return (uint64_t)rand() % (high - low + 1) + low;
}

/*
 * Compute the jacobi symbol (x/y) of x over y.  This is a direct translation
 * of the algorithm given in Otto Forster: Algorithmische Zahlentheorie.
 */
int jacobi(uint64_t x, uint64_t y)
{
	uint64_t t;
	int res = 1, m8;

	for (;; ) {
		x %= y;
		if (x == 0)
			return 0;
		while (even(x)) {
			x /= 2;
			m8 = y % 8;
			if (m8 == 3 || m8 == 5)
				res = -res;
		}
		if (x == 1)
			return res;
		t = x;
		x = y;
		y = t;
		if (x % 4 == 3 && y % 4 == 3)
			res = -res;
	}
}

/*
 * Initialises the random number generator using a (static) seed.
 */
void init_int(void)
{
	unsigned seed;

	//struct timeval tv;

	//gettimeofday(&tv, NULL);
	//seed = tv.tv_usec;
	seed = 123456;

	srand(seed);
}
