#include "helpers_int.h"
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>

/*
 * Given an integer 0 ≤ n < 2³², find the largest integer r such that r² ≤ n.
 */
unsigned long int_sqrt(unsigned long n)
{
	unsigned long root = (unsigned long) (sqrt((double)n));
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
bool is_square(unsigned long n)
{
	unsigned long sqrt = int_sqrt(n);
	return sqrt*sqrt == n;
}

/*
 * Calculate [s], [d] such that [n-1=2^s*d] where [d] is odd.
 */
void split_int(unsigned long *s, unsigned long *d, unsigned long n)
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
unsigned long get_random_int(unsigned long low, unsigned long high)
{
	return (unsigned long)rand() % (high - low + 1) + low;
}

#if 1
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

#else

/*
 * The following implementation of the jacobi function was compied straight
 * from the soucre code of PARI.
 */

/*  Compute the 2-adic valuation of z. */
long vals(unsigned long z)
{
  static char tab[64]={-1,0,1,12,2,6,-1,13,3,-1,7,-1,-1,-1,-1,14,10,4,-1,-1,8,-1,-1,25,-1,-1,-1,-1,-1,21,27,15,31,11,5,-1,-1,-1,-1,-1,9,-1,-1,24,-1,-1,20,26,30,-1,-1,-1,-1,23,-1,19,29,-1,22,18,28,17,16,-1};
  long s;

  if (!z) return -1;
  if (! (z&0xffffffff)) { s = 32; z >>=32; } else s = 0;
  z |= ~z + 1;
  z += z << 4;
  z += z << 6;
  z ^= z << 16; /* or  z -= z<<16 */
  return s + tab[(z&0xffffffff)>>26];
}

/* Figure out whether t is 3 or 5 modulo 8 */
#define  ome(t) (labs(((t)&7)-4) == 1)

int jacobi(unsigned long x, unsigned long y)
{
	unsigned long z;
	int s = 1;

	while (x) {
		long r = vals(x);
		if (r) {
			if (odd(r) && ome(y))
				s = -s;
			x >>= r;
		}
		if (x & y & 2)
			s = -s;
		z = y % x;
		y = x;
		x = z;
	}

	return (y == 1) ? s : 0;
}
#endif

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
