#include "helpers_int.h"
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

/*
 * Given an integer 0 ≤ n < 2³², find the largest integer r such that r² ≤ n.
 */
unsigned long int_sqrt(const unsigned long n)
{
	unsigned long root = (unsigned long)(sqrt((double)n));

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
unsigned long gcd(unsigned long a, unsigned long b)
{
	unsigned long t;

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
bool is_square(const unsigned long n)
{
	unsigned long sqrt = int_sqrt(n);

	return sqrt * sqrt == n;
}

/*
 * Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd.
 */
void split_int(unsigned long *s, unsigned long *d, const unsigned long n)
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
unsigned long get_random_int(const unsigned long low, const unsigned long high)
{
	return (unsigned long)rand() % (high - low + 1) + low;
}

/*
 * Compute the jacobi symbol (x/y) of x over y.  This is a direct translation
 * of the algorithm given in Otto Forster: Algorithmische Zahlentheorie.
 */
int jacobi(unsigned long x, unsigned long y)
{
	unsigned long t;
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
