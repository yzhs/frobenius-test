#include "helpers_int.h"
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

/*
 * Given an integer 0 ≤ n < 2³², find the largest integer r such that r² ≤ n.
 */
unsigned long int_sqrt(unsigned long n)
{
	unsigned long root = (unsigned long)sqrt((double)n), new_root;
	return root;
	do {
		new_root = (root + n / root) / 2;
	} while (labs((long)new_root - root) < 1);
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
bool is_square(unsigned long n)
{
	unsigned long sqrt = int_sqrt(n);
	return sqrt*sqrt == n;
}

/*
 * Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd.
 */
unsigned long split(unsigned long *d, unsigned long n)
{
	unsigned long s = 0;
	*d = n - 1;

	while (*d % 2 == 1) {
		s++;
		*d /= 2;
	}

	return s;
}

/*
 * This function generates a random integer between 2 and n-2.  This function
 * will fail if n is equal to 3 and produce weird results for smaller ns.
 */
unsigned long get_random(unsigned long n)
{
	return (unsigned long)rand() % (n - 3) + 2;
}

/*
 * Compute the jacobi symbol (x/y) of x over y.  This algorithm is taken from
 * Otto Forster: Algorithmische Zahlentheorie.
 */
int jacobi(unsigned long x, unsigned long y)
{
	unsigned long t;
	int res = 1, m8;

	for (;;) {
		x %= y;
		if (x == 0)
			return 0;
		while (x % 2 == 0) {
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
void init()
{
	unsigned seed;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	seed = tv.tv_usec;
	seed = 123456;

	srand(seed);
}
